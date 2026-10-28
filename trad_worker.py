# -*- coding: utf-8 -*-
"""
trad_worker.py

PDAL-delen av DetekteraTrad.pyt: läser ett block ur Lantmäteriets COPC-filer
eller ur lokala LAS/LAZ-filer, räknar höjd över mark, bygger en ythöjdmodell
över mark (CHM) och hittar trädtoppar i den.

Körs som egen process, en per samtidigt block, och importerar aldrig arcpy.
Flera PDAL-pipelines i trådar i samma process som arcpy kraschade processen
(access violation i arcpy, uppmätt i Pro 3.6; båda laddar Esris gdal_e.dll).

Protokoll: en uppgift per rad som JSON på stdin, ett svar per rad på stdout
med prefixet RESULT_PREFIX (annan utskrift från PDAL ignoreras av anroparen).
Processen avslutas när stdin stängs.

Uppgifter:
    {"op": "header", "paths": [...]}
        Huvudet i lokala LAS/LAZ-filer: utbredning, punktantal, koordinatsystem.
    {"op": "block", "block": {...}, "token": "...", "params": {...}, "workdir": "..."}
        Trädtoppar i ett block, sparade som .npz i workdir.
Svar: {"ok": true, "result": {...}} eller {"ok": false, "error": "traceback"}

Metod (lokala maxima i en utjämnad CHM med sökfönster som växer med höjden,
Popescu & Wynne 2004; samma princip som lmf() i R-paketet lidR):
1. Höjd över mark: markpunkterna (klass 2) trianguleras och rastreras till
   DTM_RES, som interpoleras bilinjärt i varje punkts läge (som nrd_worker.py).
2. CHM: förstareturerna trianguleras med höjd över mark som Z och rastreras
   till CHM_RES. En triangulering ger en sammanhängande yta; ett raster med
   högsta punkten per cell hade 47 % tomma celler vid 0,5 m (uppmätt, ca 2
   förstareturer per m²).
3. Gaussisk utjämning med standardavvikelsen params["sigma"] i meter.
4. En cell är en topp om den är högst inom en cirkel med diametern
   a * h + b (h = cellens utjämnade höjd), avrundad till hela halvmetrar i
   radie, och minst hmin. Sammanhängande celler med samma toppvärde blir en topp.
5. Toppens läge och höjd tas från den högsta förstareturen inom APEX_RADIUS
   från cellen, eftersom utjämningen flyttar maximum något.
6. Toppen behålls om läget ligger i blocket (halvöppet, x0 <= X < x1), så
   att varje topp hamnar i exakt ett block trots att blocken läses med marginal.
"""

import json
import os
import sys
import traceback

import numpy as np

RESULT_PREFIX = "@@TRAD_RESULT "

NODATA = -9999.0
CLASS_GROUND = 2
NOISE_CLASSES = (7, 18)

# Upplösning (m) för det interna markrutnätet som punkthöjderna räknas från.
DTM_RES = 0.5
# Upplösning (m) för CHM:en som topparna söks i.
CHM_RES = 0.5
# Toppens läge: högsta förstareturen inom så här många meter från toppcellen.
APEX_RADIUS = 1.0
# Sökfönstrets radie avrundas till hela halvmetrar inom detta intervall.
MIN_RADIUS = 1.0
MAX_RADIUS = 6.0
RADIUS_STEP = 0.5

# Läsning av lokala filer: punkter per strömmad bit.
STREAM_CHUNK = 1000000


def import_pdal():
    # Pro:s PDAL-bygge har curl utan CA-certifikat. Utan detta misslyckas varje
    # HTTPS-anslutning och arbiter försöker om i all oändlighet.
    try:
        import certifi
        os.environ.setdefault("ARBITER_CA_INFO", certifi.where())
    except ImportError:
        pass
    import pdal
    return pdal


def peak_memory_bytes():
    """Processens högsta arbetsminne hittills, eller None."""
    try:
        import psutil
        return int(psutil.Process().memory_info().peak_wset)
    except Exception:
        return None


# ── Läsning ───────────────────────────────────────────────────────────────────

def _keep(a):
    """X, Y, Z, klass och returnummer som float64, utan brus."""
    keep = ~np.isin(a["Classification"], NOISE_CLASSES)
    return np.stack([a["X"][keep], a["Y"][keep], a["Z"][keep],
                     a["Classification"][keep].astype(np.float64),
                     a["ReturnNumber"][keep].astype(np.float64)], axis=1)


def read_copc(pdal, href, bounds, token):
    """Punkterna inom bounds ur en fjärr-COPC-fil. Bara de oktree-noder som behövs hämtas."""
    stages = [{"type": "readers.copc", "bounds": bounds,
               "filename": {"path": href, "headers": {"Authorization": "Bearer " + token}}}]
    pipe = pdal.Pipeline(json.dumps(stages))
    pipe.execute()
    out = [_keep(a) for a in pipe.arrays]
    del pipe
    return out


def read_las(pdal, path, bounds):
    """
    Punkterna inom bounds ur en lokal LAS/LAZ-fil. Strömmas i bitar, så att
    hela filen aldrig ligger i minnet; readers.las har ingen bounds-option och
    läser annars in alla punkter innan filters.crop.
    """
    stages = [{"type": "readers.las", "filename": path},
              {"type": "filters.crop", "bounds": bounds}]
    pipe = pdal.Pipeline(json.dumps(stages))
    out = []
    for a in pipe.iterator(chunk_size=STREAM_CHUNK):
        if len(a):
            out.append(_keep(a))
    del pipe
    return out


def las_header(pdal, path):
    """Utbredning, punktantal och horisontellt koordinatsystem (WKT) ur filhuvudet."""
    pipe = pdal.Pipeline(json.dumps([{"type": "readers.las", "filename": path, "count": 0}]))
    pipe.execute()
    md = pipe.metadata["metadata"]["readers.las"]
    srs = md.get("srs") or {}
    return {"path": path, "bbox": [md["minx"], md["miny"], md["maxx"], md["maxy"]],
            "count": int(md["count"]), "wkt": srs.get("horizontal") or "",
            "year": md.get("creation_year"), "doy": md.get("creation_doy")}


# ── Rutnät ────────────────────────────────────────────────────────────────────

def tin_grid(pdal, pts, x0, y0, x1, y1, res, path):
    """
    Triangulera pts (strukturerad array med X, Y, Z) och rastrera till res.
    Returnerar (array, x0, y_top, res) med NaN utanför trianguleringen, rad 0 i
    norr. En enda array: filters.delaunay triangulerar varje array för sig.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    width = int(np.ceil((x1 - x0) / res))
    height = int(np.ceil((y1 - y0) / res))
    stages = [
        {"type": "filters.delaunay"},
        {"type": "filters.faceraster", "resolution": res, "origin_x": x0, "origin_y": y0,
         "width": width, "height": height},
        {"type": "writers.raster", "filename": path, "data_type": "float32", "nodata": NODATA},
    ]
    pdal.Pipeline(json.dumps(stages), arrays=[pts]).execute()
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    arr = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    ds = None
    gdal.GetDriverByName("GTiff").Delete(path)
    arr[arr == NODATA] = np.nan
    # writers.raster skriver norr uppåt: gt[3] är övre kanten, gt[5] negativ.
    return arr, gt[0], gt[3], gt[1]


def sample_bilinear(grid, gx0, gytop, res, x, y):
    """Bilinjär interpolation i cellmittpunkterna. NaN där någon granne saknas."""
    fx = (x - gx0) / res - 0.5
    fy = (gytop - y) / res - 0.5
    h, w = grid.shape
    c0 = np.clip(np.floor(fx).astype(np.int64), 0, w - 2)
    r0 = np.clip(np.floor(fy).astype(np.int64), 0, h - 2)
    tx = np.clip(fx - c0, 0.0, 1.0)
    ty = np.clip(fy - r0, 0.0, 1.0)
    z00 = grid[r0, c0]
    z01 = grid[r0, c0 + 1]
    z10 = grid[r0 + 1, c0]
    z11 = grid[r0 + 1, c0 + 1]
    return (z00 * (1 - tx) * (1 - ty) + z01 * tx * (1 - ty)
            + z10 * (1 - tx) * ty + z11 * tx * ty)


def _structured(x, y, z):
    a = np.zeros(len(x), dtype=[("X", "f8"), ("Y", "f8"), ("Z", "f8")])
    a["X"], a["Y"], a["Z"] = x, y, z
    return a


# ── Toppar ────────────────────────────────────────────────────────────────────

def window_radius(h, a, b):
    """Sökfönstrets radie (m) för höjden h: diameter a*h + b, i hela halvmetrar."""
    r = np.round((a * h + b) / 2.0 / RADIUS_STEP) * RADIUS_STEP
    return np.clip(r, MIN_RADIUS, MAX_RADIUS)


def local_maxima(chm, res, sigma, a, b, hmin):
    """
    Toppceller i chm (rad 0 i norr, NaN = ingen yta). Returnerar (rader,
    kolumner) som flyttal (tyngdpunkten av en platå) och den utjämnade CHM:en.
    """
    from scipy import ndimage
    valid = np.isfinite(chm)
    surf = np.where(valid, chm, 0.0)
    sm = ndimage.gaussian_filter(surf, sigma / res) if sigma > 0 else surf
    want = window_radius(sm, a, b)
    tops = np.zeros(sm.shape, bool)
    for rad in np.unique(want[valid & (sm >= hmin)]):
        k = int(np.ceil(rad / res))
        yy, xx = np.mgrid[-k:k + 1, -k:k + 1]
        footprint = (xx * xx + yy * yy) * res * res <= rad * rad
        mx = ndimage.maximum_filter(sm, footprint=footprint, mode="nearest")
        tops |= (want == rad) & (sm >= mx)
    tops &= valid & (sm >= hmin)
    lab, n = ndimage.label(tops, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros(0), np.zeros(0), sm
    rc = np.array(ndimage.center_of_mass(tops, lab, range(1, n + 1)), ndmin=2)
    return rc[:, 0], rc[:, 1], sm


def find_tops(fx, fy, fh, gx0, gytop, res, rows, cols, sm):
    """
    Toppens läge och höjd: den högsta förstareturen inom APEX_RADIUS från
    toppcellens mitt. Saknas en sådan används cellens mitt och den utjämnade
    höjden. Två toppar som får samma retur blir en.
    """
    from scipy.spatial import cKDTree
    cx = gx0 + (cols + 0.5) * res
    cy = gytop - (rows + 0.5) * res
    ch = sm[np.round(rows).astype(int), np.round(cols).astype(int)]
    if len(cx) == 0:
        return cx, cy, ch
    tree = cKDTree(np.column_stack([fx, fy]))
    out_x, out_y, out_h, used = [], [], [], set()
    for x, y, h, idx in zip(cx, cy, ch, tree.query_ball_point(np.column_stack([cx, cy]),
                                                            APEX_RADIUS)):
        if idx:
            j = idx[int(np.argmax(fh[idx]))]
            if j in used:
                continue
            used.add(j)
            out_x.append(fx[j])
            out_y.append(fy[j])
            out_h.append(fh[j])
        else:
            out_x.append(x)
            out_y.append(y)
            out_h.append(h)
    return np.array(out_x), np.array(out_y), np.array(out_h)


# ── Block ─────────────────────────────────────────────────────────────────────

def process_block(pdal, task):
    """Läs ett block, räkna höjd över mark, bygg CHM och hitta trädtoppar."""
    blk = task["block"]
    prm = task["params"]
    rx0, ry0, rx1, ry1 = blk["read_bbox"]
    tag = "{:03d}_{:03d}".format(blk["row"], blk["col"])
    workdir = task["workdir"]

    arrays = []
    for rd in blk["reads"]:
        if rd.get("path"):
            arrays.extend(read_las(pdal, rd["path"], rd["bounds"]))
        else:
            arrays.extend(read_copc(pdal, rd["href"], rd["bounds"], task["token"]))
    pts = np.concatenate(arrays) if arrays else np.zeros((0, 5))
    del arrays

    x, y = pts[:, 0], pts[:, 1]
    inside = (x >= blk["x0"]) & (x < blk["x1"]) & (y >= blk["y0"]) & (y < blk["y1"])
    is_ground = pts[:, 3] == CLASS_GROUND
    is_first = pts[:, 4] == 1
    result = {"points": int(inside.sum()), "ground": int((inside & is_ground).sum()),
              "first": int((inside & is_first).sum()), "no_height": 0, "tops": 0,
              "npz": None, "peak_bytes": None}
    if not (inside & is_first).any() or is_ground.sum() < 3:
        result["peak_bytes"] = peak_memory_bytes()
        return result

    ground = _structured(pts[is_ground, 0], pts[is_ground, 1], pts[is_ground, 2])
    first = pts[is_first]
    del pts, x, y, inside, is_ground, is_first
    dtm, dx0, dytop, dres = tin_grid(pdal, ground, rx0, ry0, rx1, ry1, DTM_RES,
                                     os.path.join(workdir, "dtm_{}.tif".format(tag)).replace("\\", "/"))
    del ground
    fx, fy = first[:, 0], first[:, 1]
    fh = first[:, 2] - sample_bilinear(dtm, dx0, dytop, dres, fx, fy)
    del first
    ok = np.isfinite(fh)
    fin = (fx >= blk["x0"]) & (fx < blk["x1"]) & (fy >= blk["y0"]) & (fy < blk["y1"])
    result["no_height"] = int((fin & ~ok).sum())
    fx, fy, fh = fx[ok], fy[ok], fh[ok]
    if len(fx) < 3:
        result["peak_bytes"] = peak_memory_bytes()
        return result

    chm, cx0, cytop, cres = tin_grid(pdal, _structured(fx, fy, np.maximum(fh, 0.0)),
                                     rx0, ry0, rx1, ry1, CHM_RES,
                                     os.path.join(workdir, "chm_{}.tif".format(tag)).replace("\\", "/"))
    rows, cols, sm = local_maxima(chm, cres, prm["sigma"], prm["a"], prm["b"], prm["hmin"])
    del chm
    tx, ty, th = find_tops(fx, fy, fh, cx0, cytop, cres, rows, cols, sm)
    del sm
    keep = ((tx >= blk["x0"]) & (tx < blk["x1"]) & (ty >= blk["y0"]) & (ty < blk["y1"])
            & (th >= prm["hmin"]))
    tx, ty, th = tx[keep], ty[keep], th[keep]
    tz = sample_bilinear(dtm, dx0, dytop, dres, tx, ty)
    result["tops"] = int(len(tx))
    if len(tx):
        npz = os.path.join(workdir, "tops_{}.npz".format(tag))
        np.savez(npz, x=tx, y=ty, h=th, ground=tz)
        result["npz"] = npz
    result["peak_bytes"] = peak_memory_bytes()
    return result


def handle(pdal, task):
    if task.get("op") == "header":
        return {"headers": [las_header(pdal, p) for p in task["paths"]]}
    return process_block(pdal, task)


def main():
    # Sökvägar kan innehålla å, ä och ö; läs och skriv alltid UTF-8.
    sys.stdin.reconfigure(encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    pdal = import_pdal()
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            out = {"ok": True, "result": handle(pdal, json.loads(line))}
        except Exception:
            out = {"ok": False, "error": traceback.format_exc()}
        sys.stdout.write(RESULT_PREFIX + json.dumps(out) + "\n")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
