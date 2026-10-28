# -*- coding: utf-8 -*-
"""
Test av DetekteraTrad.pyt genom geobehandlingsramverket (arcpy.ImportToolbox).

Kör med ArcGIS Pro:s Python, utanför Pro:
    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" test_detektera_trad.py <fil.las>
eller i Pro:s Python-fönster:
    import sys; sys.argv = ["", r"<fil.las>"]; exec(open(r"<sökväg>\\test_detektera_trad.py", encoding="utf-8").read())

Använder en lokal LAS/LAZ-fil i SWEREF 99 TM (Lantmäteriets laserdata) och
kräver ingen inloggning. Utdata skrivs till en temporär geodatabas som tas
bort efteråt.

1. Detektera träd två gånger i samma process, med utbredning som inskrivna
   koordinater och med de valfria parametrarna tomma respektive ifyllda.
2. Utvärdera noggrannhet mot syntetiska referensträd: de detekterade träden
   i en provyta, förskjutna slumpvis (0,5 m), där 20 % tagits bort och 10 %
   påhittade lagts till. Förväntat: precision ca 0,8, träffsäkerhet ca 0,8/0,9.
"""

import os
import random
import shutil
import sys
import tempfile

import arcpy

HERE = os.path.dirname(os.path.abspath(__file__))
PYT = os.path.join(HERE, "DetekteraTrad.pyt")


def main(las):
    work = tempfile.mkdtemp(prefix="trad_test_")
    gdb = arcpy.management.CreateFileGDB(work, "test.gdb")[0]
    arcpy.env.overwriteOutput = True
    arcpy.ImportToolbox(PYT)
    tb = arcpy.detektera_trad
    print("Verktygstips lästa:", "Lägsta trädhöjd" in (tb.DetekteraTrad.__doc__ or "") or
          "trädhöjd" in (tb.DetekteraTrad.__doc__ or ""))

    # Mitt i filen, 400 x 400 m.
    import json
    import subprocess
    hdr = subprocess.run([os.path.join(sys.exec_prefix, "python.exe"), "-c",
                          "import json,sys,pdal;p=pdal.Pipeline(json.dumps([{'type':'readers.las',"
                          "'filename':sys.argv[1],'count':0}]));p.execute();m=p.metadata['metadata']"
                          "['readers.las'];print(json.dumps([m['minx'],m['miny'],m['maxx'],m['maxy']]))",
                          las], capture_output=True, text=True, check=True)
    x0, y0, x1, y1 = json.loads(hdr.stdout.strip().splitlines()[-1])
    cx, cy = round((x0 + x1) / 2), round((y0 + y1) / 2)
    ext = "{} {} {} {}".format(cx - 200, cy - 200, cx + 200, cy + 200)

    out1 = os.path.join(gdb, "toppar1")
    r = tb.DetekteraTrad("Utbredning (kartvy, lager eller koordinater)", None, ext,
                         "Lokala LAS/LAZ-filer", las, out1)
    print(r.getMessages())
    n1 = int(arcpy.management.GetCount(out1)[0])

    out2 = os.path.join(gdb, "toppar2")
    tb.DetekteraTrad("Utbredning (kartvy, lager eller koordinater)", None, ext,
                     "Lokala LAS/LAZ-filer", las, out2, 4, 1.0, 0.1, 3, None, None, 2, 2)
    n2 = int(arcpy.management.GetCount(out2)[0])
    print("Toppar, standard:", n1, " utjämning 1 m:", n2)
    assert n1 > 0 and 0 < n2 < n1, "fler toppar med mer utjämning"
    sr = arcpy.Describe(out1).spatialReference
    print("Utdata:", sr.name, [f.name for f in arcpy.ListFields(out1)])
    assert sr.factoryCode == 3006

    # Syntetiska referensträd i en provyta (100 x 100 m mitt i området).
    plots = arcpy.management.CreateFeatureclass(gdb, "provytor", "POLYGON",
                                                spatial_reference=sr)[0]
    arcpy.management.AddField(plots, "NAMN", "TEXT", field_length=20)
    with arcpy.da.InsertCursor(plots, ["SHAPE@", "NAMN"]) as cur:
        for name, (px, py) in (("A", (cx - 120, cy - 120)), ("B", (cx + 20, cy + 20))):
            cur.insertRow((arcpy.Polygon(arcpy.Array([arcpy.Point(px, py), arcpy.Point(px, py + 100),
                                                      arcpy.Point(px + 100, py + 100),
                                                      arcpy.Point(px + 100, py)]), sr), name))
    ref = arcpy.management.CreateFeatureclass(gdb, "referens", "POINT", spatial_reference=sr)[0]
    arcpy.management.AddField(ref, "HOJD", "DOUBLE")
    rnd = random.Random(7)
    kept = dropped = added = 0
    with arcpy.da.SearchCursor(plots, ["SHAPE@"]) as pc:
        polys = [p for (p,) in pc]
    with arcpy.da.SearchCursor(out1, ["SHAPE@XY", "TRADHOJD"]) as cur, \
            arcpy.da.InsertCursor(ref, ["SHAPE@XY", "HOJD"]) as ins:
        for (x, y), h in cur:
            if not any(p.contains(arcpy.PointGeometry(arcpy.Point(x, y), sr)) for p in polys):
                continue
            if rnd.random() < 0.2:
                dropped += 1
                continue
            ins.insertRow(((x + rnd.gauss(0, 0.5), y + rnd.gauss(0, 0.5)), h + rnd.gauss(0, 1)))
            kept += 1
        for p in polys:
            e = p.extent
            for _ in range(int(0.1 * kept / len(polys))):
                ins.insertRow(((rnd.uniform(e.XMin, e.XMax), rnd.uniform(e.YMin, e.YMax)), 15.0))
                added += 1
    print("Referens: {} behållna, {} borttagna, {} påhittade".format(kept, dropped, added))

    out_eval = os.path.join(gdb, "utvardering")
    out_tab = os.path.join(gdb, "matt")
    r = tb.UtvarderaNoggrannhet(out1, ref, plots, out_eval, out_tab, 2.0, "TRADHOJD", "HOJD",
                                30, "NAMN")
    print(r.getMessages())
    with arcpy.da.SearchCursor(out_tab, ["PROVYTA", "TRAFFSAKERHET", "PRECISION", "F1"]) as cur:
        rows = list(cur)
    print(rows)
    alla = [x for x in rows if x[0] == "Alla"][0]
    exp_p = kept / (kept + dropped)
    print("Förväntad precision ca {:.2f}, fick {:.2f}".format(exp_p, alla[2]))
    assert abs(alla[2] - exp_p) < 0.05
    # Utan valfria parametrar.
    tb.UtvarderaNoggrannhet(out1, ref, plots, os.path.join(gdb, "utv2"), None, 2.0)
    del r
    arcpy.management.Delete(gdb)
    shutil.rmtree(work, ignore_errors=True)
    print("OK")


if __name__ == "__main__":
    main(sys.argv[1])
