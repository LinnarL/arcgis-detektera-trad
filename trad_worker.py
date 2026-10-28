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
7. Varje topp märks om den står på en hård, plan yta (hard_surface_mask: tak,
   cisterner) eller är ett tunt föremål (låg täckning och få returer:
   ledningar). Anroparen tar bort dem om användaren valt det. Broar (klass
   17) ingår aldrig i CHM:en.
8. Kronor (om params["crowns"]): markörstyrd vattendelare på den utjämnade
   CHM:en med topparna som markörer (crown_labels), polygoniserade med GDAL
   och förenklade (crown_polygons). Skickas som WKB i npz-filen.
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
# Bro (17) ger plana, höga ytor och hålls utanför CHM:en.
CLASS_BRIDGE = 17
# Ytmått runt toppen (top_features): returer minst så här högt över mark,
# inom sökfönstrets radie men minst FEATURE_MIN_RADIUS.
FEATURE_MIN_H = 2.0
FEATURE_MIN_RADIUS = 1.5
# Hårda, plana ytor (hard_surface_mask).
HARD_WINDOW_M = 2.5
HARD_MIN_H = 2.0
HARD_MAX_ROUGH = 0.2
HARD_MIN_SINGLE = 0.7
HARD_MIN_POINTS = 4
HARD_MIN_AREA_M2 = 10.0
HARD_GROW_M = 1.5
# Ledningar och andra tunna föremål: en topp med låg täckning och få returer.
# Uppmätt på 2026 års skanning (ca 3 punkter per m²): ledningstoppar hade
# täckning 0,05-0,15 (median) och 1-3 returer, kronor i tät skog minst 0,21 och
# i häckar och gårdar 0,38 (1:a percentilen). Regeln tog 89-95 % av topparna
# längs två ledningar och 2 av 1640 toppar i tät skog.
THIN_MAX_COVER = 0.3
THIN_MAX_POINTS = 8
# Kronor (crown_labels): en cell hör till en krona om den är minst
# CROWN_MIN_REL av trädets höjd och högst CROWN_MAX_R_FACTOR x höjden från
# toppen (silva2016 i lidR: exclusion 0,3, max_cr_factor 0,6 som diameter).
# Ingen cell under CROWN_MIN_H över mark.
CROWN_MIN_REL = 0.5
CROWN_MAX_R_FACTOR = 0.3
CROWN_MIN_H = 2.0
CROWN_LEVEL_M = 0.25
# Polygonernas trappsteg (cellkanter) förenklas med denna tolerans i meter.
CROWN_SIMPLIFY_M = 0.3
MIN_RING_AREA_M2 = 0.05

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
    """X, Y, Z, klass, returnummer och antal returer som float64, utan brus."""
    keep = ~np.isin(a["Classification"], NOISE_CLASSES)
    return np.stack([a["X"][keep], a["Y"][keep], a["Z"][keep],
                     a["Classification"][keep].astype(np.float64),
                     a["ReturnNumber"][keep].astype(np.float64),
                     a["NumberOfReturns"][keep].astype(np.float64)], axis=1)


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


def hard_surface_mask(chm, fx, fy, fh, fsingle, gx0, gytop, res):
    """
    Celler (rad 0 i norr, samma rutnät som chm) som ligger på byggnader och
    andra hårda, plana ytor ovan mark, t.ex. tak, cisterner och bropelare.

    En cell räknas som hård yta om, i ett fönster på HARD_WINDOW_M:
      - CHM:en är minst HARD_MIN_H hög,
      - CHM:en avviker högst HARD_MAX_ROUGH m (RMS) från ett anpassat plan,
        vilket ett tak gör men en krona inte gör, och
      - minst HARD_MIN_SINGLE av förstareturerna där kommer från pulser med
        bara en retur: ett tak delar inte pulsen, grenar gör det.
    Sammanhängande ytor under HARD_MIN_AREA_M2 tas bort, och masken växer med
    HARD_GROW_M så att takkanter, nockar och utrustning på taket kommer med;
    de ger toppar men inte plana fönster.
    """
    from scipy import ndimage
    valid = np.isfinite(chm)
    z = np.where(valid, chm, 0.0)
    k = max(3, int(round(HARD_WINDOW_M / res)) | 1)
    half = k // 2
    u = (np.arange(k) - half) * res
    ku = np.tile(u, (k, 1))           # x-led
    kv = ku.T                          # y-led
    n = float(k * k)
    mean = ndimage.uniform_filter(z, k, mode="nearest")
    mean2 = ndimage.uniform_filter(z * z, k, mode="nearest")
    su = ndimage.correlate(z, ku, mode="nearest")
    sv = ndimage.correlate(z, kv, mode="nearest")
    suu = float((ku * ku).sum())
    var = mean2 - mean * mean - (su * su) / (n * suu) - (sv * sv) / (n * suu)
    rough = np.sqrt(np.maximum(var, 0.0))

    # Andel enkelreturer bland förstareturerna minst HARD_MIN_H över mark.
    h, w = chm.shape
    high = fh >= HARD_MIN_H
    col = np.floor((fx[high] - gx0) / res).astype(np.int64)
    row = np.floor((gytop - fy[high]) / res).astype(np.int64)
    ok = (col >= 0) & (col < w) & (row >= 0) & (row < h)
    idx = row[ok] * w + col[ok]
    n_all = np.bincount(idx, minlength=h * w).reshape(h, w).astype(np.float64)
    n_one = np.bincount(idx, weights=fsingle[high][ok], minlength=h * w).reshape(h, w)
    s_all = ndimage.uniform_filter(n_all, k, mode="constant") * n
    s_one = ndimage.uniform_filter(n_one, k, mode="constant") * n
    share = np.where(s_all >= HARD_MIN_POINTS, s_one / np.maximum(s_all, 1), 0.0)

    hard = valid & (z >= HARD_MIN_H) & (rough <= HARD_MAX_ROUGH) & (share >= HARD_MIN_SINGLE)
    lab, nlab = ndimage.label(hard)
    if nlab:
        sizes = ndimage.sum(hard, lab, range(1, nlab + 1)) * res * res
        hard = np.isin(lab, np.nonzero(sizes >= HARD_MIN_AREA_M2)[0] + 1)
    g = int(round(HARD_GROW_M / res))
    if g > 0 and hard.any():
        yy, xx = np.mgrid[-g:g + 1, -g:g + 1]
        hard = ndimage.binary_dilation(hard, structure=(xx * xx + yy * yy) <= g * g)
    return hard


def crown_labels(sm, hard, tx, ty, th, gx0, gytop, res):
    """
    Kronor som etikettraster (0 = ingen krona, j + 1 = toppen j): en
    vattendelare som fyller den utjämnade CHM:en uppifrån, CROWN_LEVEL_M i
    taget, från en markör per topp. På varje nivå växer kronorna cell för cell
    in i celler som nått nivån, så att gränsen mellan två kronor hamnar i
    sadeln mellan dem. En cell tas bara med om den är minst CROWN_MIN_H över
    mark, inte hård yta, minst CROWN_MIN_REL av toppens höjd i den utjämnade
    CHM:en och inom CROWN_MAX_R_FACTOR x trädhöjden från toppen. Villkoren
    prövas när cellen läggs till, så varje krona hänger ihop med sin topp.

    scipy.ndimage.watershed_ift är ingen vanlig vattendelare: med två toppar
    fick den högre bara sin egen cell och resten gick till den lägre
    (uppmätt), och med bakgrundsmarkörer gick en fjärdedel av krontaket till
    bakgrunden. Därför denna egna fyllning.
    """
    h, w = sm.shape
    labels = np.zeros(h * w, np.int32)
    if len(tx) == 0:
        return labels.reshape(h, w)
    col = np.clip(np.floor((tx - gx0) / res).astype(np.int64), 0, w - 1)
    row = np.clip(np.floor((gytop - ty) / res).astype(np.int64), 0, h - 1)
    labels[row * w + col] = np.arange(1, len(tx) + 1, dtype=np.int32)
    flat = sm.ravel()
    h_top = flat[row * w + col]
    min_h = np.concatenate([[0.0], np.maximum(CROWN_MIN_REL * h_top, CROWN_MIN_H)])
    max_r = np.concatenate([[0.0], CROWN_MAX_R_FACTOR * th + res])
    ax = np.concatenate([[0.0], tx])
    ay = np.concatenate([[0.0], ty])
    eligible = ((sm >= CROWN_MIN_H) & ~hard).ravel()
    offsets = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    top = float(flat[eligible].max()) if eligible.any() else CROWN_MIN_H
    for level in np.arange(top, CROWN_MIN_H - 1e-9, -CROWN_LEVEL_M):
        # Bara cellerna som väntar på en etikett på den här nivån gås igenom,
        # inte hela rastret: det gjorde ett block på 1 km² 4 gånger långsammare.
        free = np.flatnonzero(eligible & (flat >= level) & (labels == 0))
        while free.size:
            r, c = np.divmod(free, w)
            best = np.zeros(free.size, np.int32)
            for dr, dc in offsets:
                rr, cc = r + dr, c + dc
                inb = (rr >= 0) & (rr < h) & (cc >= 0) & (cc < w)
                nb = np.zeros(free.size, np.int32)
                nb[inb] = labels[rr[inb] * w + cc[inb]]
                np.maximum(best, nb, out=best)
            new = best > 0
            if not new.any():
                break
            idx, lab = free[new], best[new]
            r, c = r[new], c[new]
            cx = gx0 + (c + 0.5) * res
            cy = gytop - (r + 0.5) * res
            ok = (flat[idx] >= min_h[lab]) & (np.hypot(cx - ax[lab], cy - ay[lab]) <= max_r[lab])
            labels[idx[ok]] = lab[ok]
            # En cell som inte klarar villkoren för sin granne prövas inte igen.
            eligible[idx[~ok]] = False
            free = free[~new]
    return labels.reshape(h, w)


def _ring_area(pts):
    """Ringens yta med tecken (positiv moturs), skosnöresformeln."""
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def _ring_wkb(pts, ccw):
    """En ring som WKB-del, orienterad moturs (ccw) eller medurs."""
    import struct
    if (_ring_area(pts) > 0) != ccw:
        pts = pts[::-1]
    return struct.pack("<I", len(pts)) + np.ascontiguousarray(pts, "<f8").tobytes()


def crown_polygons(labels, keep_idx, gx0, gytop, res):
    """
    Polygoner för kronorna till topparna keep_idx (index i etikettrastret).
    Returnerar (wkb, offs): kronan k är wkb[offs[k]:offs[k + 1]] som WKB
    MultiPolygon (tom om toppen saknar krona). Ringarna orienteras enligt OGC
    (yttre moturs, hål medurs); fel orientering lagras annars tyst som negativ
    yta (se CLAUDE.md).
    """
    import struct
    from osgeo import gdal, ogr
    gdal.UseExceptions()
    n = len(keep_idx)
    lut = np.zeros(int(labels.max()) + 1 if labels.size else 1, np.int32)
    lut[keep_idx + 1] = np.arange(1, n + 1, dtype=np.int32)
    out = lut[labels]
    parts = [[] for _ in range(n)]
    if out.any():
        h, w = out.shape
        ds = gdal.GetDriverByName("MEM").Create("", w, h, 1, gdal.GDT_Int32)
        ds.SetGeoTransform([gx0, res, 0.0, gytop, 0.0, -res])
        band = ds.GetRasterBand(1)
        band.WriteArray(out)
        vds = gdal.GetDriverByName("MEM").Create("", 0, 0, 0, gdal.GDT_Unknown)
        lyr = vds.CreateLayer("kronor", geom_type=ogr.wkbPolygon)
        lyr.CreateField(ogr.FieldDefn("k", ogr.OFTInteger))
        gdal.Polygonize(band, band, lyr, 0, ["8CONNECTED=8"])
        for feat in lyr:
            k = feat.GetField(0) - 1
            geom = feat.GetGeometryRef().SimplifyPreserveTopology(CROWN_SIMPLIFY_M)
            if geom is None or geom.IsEmpty():
                continue
            polys = [geom] if geom.GetGeometryType() == ogr.wkbPolygon else \
                [geom.GetGeometryRef(i) for i in range(geom.GetGeometryCount())]
            for poly in polys:
                # Förenklingen kan lämna ringar utan yta: en sådan yttre ring
                # stryker polygonen, ett sådant hål stryks.
                rings = [np.array(poly.GetGeometryRef(i).GetPoints())[:, :2]
                         for i in range(poly.GetGeometryCount())]
                if not rings or len(rings[0]) < 4 or \
                        abs(_ring_area(rings[0])) < MIN_RING_AREA_M2:
                    continue
                holes = [r for r in rings[1:]
                         if len(r) >= 4 and abs(_ring_area(r)) >= MIN_RING_AREA_M2]
                parts[k].append([rings[0]] + holes)
        vds = ds = None
    chunks, offs = [], [0]
    for polys in parts:
        if not polys:
            offs.append(offs[-1])
            continue
        b = struct.pack("<BII", 1, 6, len(polys))
        for rings in polys:
            b += struct.pack("<BII", 1, 3, len(rings))
            b += b"".join(_ring_wkb(r, i == 0) for i, r in enumerate(rings))
        chunks.append(b)
        offs.append(offs[-1] + len(b))
    wkb = np.frombuffer(b"".join(chunks), np.uint8) if chunks else np.zeros(0, np.uint8)
    return wkb, np.array(offs, np.int64)


def top_features(fx, fy, fh, fsingle, tx, ty, th, a, b):
    """
    Mått på ytan runt varje topp, ur förstareturerna inom sökfönstrets radie
    och minst FEATURE_MIN_H över mark:
        n       antal sådana returer
        single  andel av dem från pulser med bara en retur (tak och mark nära 1,
                kronor lägre eftersom pulsen delas i grenarna)
        rough   RMS-avvikelse (m) från ett plan anpassat till dem (tak lågt)
        cover   andel av alla förstareturer inom radien som är minst
                FEATURE_MIN_H över mark: en krona fyller fönstret, en
                ledning eller en stolpe nästan inte alls
    """
    from scipy.spatial import cKDTree
    n_t = len(tx)
    n = np.zeros(n_t, np.int32)
    single = np.full(n_t, np.nan)
    rough = np.full(n_t, np.nan)
    cover = np.full(n_t, np.nan)
    if n_t == 0:
        return n, single, rough, cover
    high = fh >= FEATURE_MIN_H
    hx, hy, hh, hs = fx[high], fy[high], fh[high], fsingle[high]
    if len(hx) == 0:
        return n, single, rough, cover
    tree = cKDTree(np.column_stack([hx, hy]))
    tree_all = cKDTree(np.column_stack([fx, fy]))
    radius = np.maximum(window_radius(th, a, b), FEATURE_MIN_RADIUS)
    for k in range(n_t):
        idx = tree.query_ball_point((tx[k], ty[k]), radius[k])
        n[k] = len(idx)
        n_all = len(tree_all.query_ball_point((tx[k], ty[k]), radius[k]))
        cover[k] = n[k] / n_all if n_all else np.nan
        if not idx:
            continue
        single[k] = float(hs[idx].mean())
        if len(idx) >= 6:
            px, py, pz = hx[idx] - tx[k], hy[idx] - ty[k], hh[idx]
            A = np.column_stack([px, py, np.ones(len(idx))])
            coef, *_ = np.linalg.lstsq(A, pz, rcond=None)
            rough[k] = float(np.sqrt(np.mean((pz - A @ coef) ** 2)))
    return n, single, rough, cover


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
    pts = np.concatenate(arrays) if arrays else np.zeros((0, 6))
    del arrays

    x, y = pts[:, 0], pts[:, 1]
    inside = (x >= blk["x0"]) & (x < blk["x1"]) & (y >= blk["y0"]) & (y < blk["y1"])
    is_ground = pts[:, 3] == CLASS_GROUND
    is_first = (pts[:, 4] == 1) & (pts[:, 3] != CLASS_BRIDGE)
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
    fsingle = (first[:, 5] == 1).astype(np.float64)
    del first
    ok = np.isfinite(fh)
    fin = (fx >= blk["x0"]) & (fx < blk["x1"]) & (fy >= blk["y0"]) & (fy < blk["y1"])
    result["no_height"] = int((fin & ~ok).sum())
    fx, fy, fh, fsingle = fx[ok], fy[ok], fh[ok], fsingle[ok]
    if len(fx) < 3:
        result["peak_bytes"] = peak_memory_bytes()
        return result

    chm, cx0, cytop, cres = tin_grid(pdal, _structured(fx, fy, np.maximum(fh, 0.0)),
                                     rx0, ry0, rx1, ry1, CHM_RES,
                                     os.path.join(workdir, "chm_{}.tif".format(tag)).replace("\\", "/"))
    rows, cols, sm = local_maxima(chm, cres, prm["sigma"], prm["a"], prm["b"], prm["hmin"])
    hard = hard_surface_mask(chm, fx, fy, fh, fsingle, cx0, cytop, cres)
    del chm
    tx, ty, th = find_tops(fx, fy, fh, cx0, cytop, cres, rows, cols, sm)
    keep = ((tx >= blk["x0"]) & (tx < blk["x1"]) & (ty >= blk["y0"]) & (ty < blk["y1"])
            & (th >= prm["hmin"]))
    crowns = None
    if prm.get("crowns"):
        # Alla toppar i läsområdet är markörer, även de i marginalen, så att en
        # krona i blocket inte breder ut sig över grannträdets.
        labels = crown_labels(sm, hard, tx, ty, th, cx0, cytop, cres)
        crowns = crown_polygons(labels, np.nonzero(keep)[0], cx0, cytop, cres)
        del labels
    del sm
    tx, ty, th = tx[keep], ty[keep], th[keep]
    tz = sample_bilinear(dtm, dx0, dytop, dres, tx, ty)
    fn, fs, fr, fc = top_features(fx, fy, fh, fsingle, tx, ty, th, prm["a"], prm["b"])
    hr = np.clip(np.floor((cytop - ty) / cres).astype(np.int64), 0, hard.shape[0] - 1)
    hc = np.clip(np.floor((tx - cx0) / cres).astype(np.int64), 0, hard.shape[1] - 1)
    on_hard = hard[hr, hc] if len(tx) else np.zeros(0, bool)
    thin = (fc < THIN_MAX_COVER) & (fn <= THIN_MAX_POINTS) & ~on_hard
    result["tops"] = int(len(tx))
    result["tops_hard"] = int(on_hard.sum())
    result["tops_thin"] = int(thin.sum())
    if len(tx):
        npz = os.path.join(workdir, "tops_{}.npz".format(tag))
        extra = {"crown_wkb": crowns[0], "crown_offs": crowns[1]} if crowns else {}
        np.savez(npz, x=tx, y=ty, h=th, ground=tz, n=fn, single=fs, rough=fr, cover=fc,
                 hard=on_hard, thin=thin, **extra)
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
