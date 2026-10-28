# Detektera träd

ArcGIS Pro Python toolbox with two tools:

- **Detektera träd** finds tree tops (crown apices) in airborne laser data and writes them as points
  with tree height, and optionally delineates each tree's crown as a polygon. Points come from
  Lantmäteriet's Laserdata Nedladdning, skog for the area, or from local LAS/LAZ files.
- **Utvärdera noggrannhet** compares detected trees with reference trees in field plots and reports
  detection rate, precision, F1, stems per hectare and position error.

**The position is the crown top, not a measured stem base.** At Lantmäteriet's density (2-4 points
per m2 counting all returns) almost no pulses hit the stem, and a top can lie a metre or two from
the stem base on leaning trees or asymmetric crowns. Mainly dominant and co-dominant trees are
found; suppressed trees under other crowns rarely are. Published benchmarks at similar density
find roughly half of all trees.

**The default parameters are not calibrated against Swedish reference data.** They follow the
`lmf()` example in the R package lidR. Use Utvärdera noggrannhet with your own plots to tune them
before relying on counts.

## Requirements

- ArcGIS Pro 3.x, Basic licence is enough. Developed on 3.6 with Python 3.13.
- No extra packages. `pdal`, `scipy`, `numpy`, `certifi` and GDAL's Python bindings ship with the
  default `arcgispro-py3` environment.
- For Lantmäteriet data: a consumer key and secret from Lantmäteriet's API portal with both the API
  `STAC-hojd` and an order of
  [Laserdata Nedladdning, skog](https://geotorget.lantmateriet.se/geodataprodukter/laserdata-nedladdning-skog-api).
  The same key works for the Laserdata höjdmodeller and NRD toolboxes. Local files need no key.

## Install

1. Clone or download this repo. Keep `DetekteraTrad.pyt`, `trad_worker.py` and the two `.lyrx`
   files in the same folder; the toolbox runs the worker file in separate processes, and the layer
   files give the outputs their symbols.
2. In ArcGIS Pro: Catalog, Toolboxes, Add Toolbox, select `DetekteraTrad.pyt`.

## Detektera träd

| Parameter | Default | Notes |
|---|---|---|
| Avgränsa området med | Polygoner | Polygons (layer or drawn in the map) or extent (display, layer, drawn or typed; typed coordinates use the active map's coordinate system) |
| Intresseområde (polygoner) / Utbredning | - | The layer's selection is honoured. Only tops inside the polygons are kept |
| Punktmoln | Lantmäteriets Laserdata Skog | Or local LAS/LAZ files in SWEREF 99 TM with ground classified (class 2) |
| LAS/LAZ-filer | - | Only for local files. Files must not overlap with the same points |
| Utdata, trädtoppar | `tradtoppar` in the project geodatabase | New point feature class, drawn as a small red cross. Honours Overwrite output and Output Coordinate System |
| Utdata, trädkronor | `tradkronor` in the project geodatabase | Optional polygon feature class, drawn with a dashed outline and no fill. Clear it to skip crowns (saves about a fifth of the run time) |
| Lägsta trädhöjd (m) | 4 | Under Detektering. 0-50 m |
| Utjämning (m) | 0.5 | Under Detektering. Gaussian sigma in metres, 0-3. More smoothing gives fewer tops |
| Sökfönster, A and B | 0.1 and 3 | Under Detektering. A cell is a top if it is highest within a circle of diameter A x height + B metres |
| Ta bort toppar på byggnader och ledningar | on | Under Detektering. See "Buildings and power lines" below |
| Consumer key / secret | - | Under Inloggning. Only for Lantmäteriet data |
| Minnesbudget (GB) / Block samtidigt | automatic | Under Avancerat. Same scheme as the Laserdata toolboxes |

Tree top fields: `TRAD_ID`, `TRADHOJD` (tree height, m above ground), `MARKHOJD` (ground elevation
at the top, m RH 2000), `KALLA` (tile or file), `SKANNING` (scanning area), `INSAMLAD` (capture
period, empty for local files since a LAS header holds only the processing date).

Crown fields: `TRAD_ID` (the same as the tree top's, for joining), `TRADHOJD`, `KRONAREA` (m2) and
`KRONDIAM` (m, the diameter of a circle with the same area).

The symbols come from `tradtoppar.lyrx` and `tradkronor.lyrx`, which Pro applies when it adds the
outputs to the map. `skapa_lagerfiler.py` regenerates them with arcpy if the colours or sizes are to
be changed.

### Crowns

Crowns are grown from the tops through the smoothed canopy height model, from the top down in
0.25 m steps, so the border between two touching crowns falls in the saddle between them (a
watershed with one marker per top). A cell joins a crown only if it is at least half the top's
height, within 0.3 x the tree height of the top (the defaults of `silva2016` in lidR), at least 2 m
above ground and not part of the building mask. Each crown is therefore connected to its own top.
Edges are simplified with a 0.3 m tolerance to remove the 0.5 m cell steps, so neighbouring crowns
can overlap or leave a gap of a few decimetres.

- In 6 ha of dense forest every top got a crown, median 19.6 m2 (5-95 %: 7.5-47.5 m2), and every
  top lies inside its own crown.
- Crowns never spread over roofs. In the town test the trees kept after filtering had a median
  crown of 18 m2.
- The whole crown is kept also where it reaches outside the area of interest.
- At 2-4 points per m2 crown edges are approximate. Suppressed trees without a detected top are
  absorbed into a neighbour's crown or left out.
- Adds about 4 s per km2 to a run (16 s without crowns, 20 s with, for a 1 km2 forest block).

### How it works

1. With Lantmäteriet data, the open STAC catalogue is searched and the newest scan per place is
   used (the 2.5 km and older large tiles of one scan are never both read). Only points inside the
   area's bounding box plus a 20 m margin are read, directly from the cloud-optimised files.
2. The area is split into blocks that fit the memory budget, processed in parallel worker
   processes. Noise (classes 7 and 18) is dropped.
3. Ground points are triangulated and rasterised to 0.5 m; each first return gets its height above
   ground by bilinear interpolation.
4. Canopy height model: first returns triangulated with height above ground as Z, rasterised to
   0.5 m. A per-cell maximum grid left 47 % of cells empty at about 2 first returns per m2, the
   triangulation does not.
5. Gaussian smoothing, then a cell is a top if it is the highest within a circle of diameter
   A x height + B (radius rounded to half metres, 1-6 m) and at least the minimum height.
6. Each top is moved to the highest first return within 1 m, since smoothing shifts the maximum
   slightly. That return gives position and height.
7. A top is kept by the block it lies in (half-open bounds), so blocks never produce duplicates.
   On a test area, one 400 m block and four 200 m blocks gave 4 356 and 4 354 tops; the handful of
   differences were scattered, not on the seams, and come from Delaunay triangulation not being
   unique where points are co-circular.

### Buildings and power lines

Lantmäteriet's point cloud has no building or vegetation class, so roofs, roof ridges and power
lines all produce "tops". With the option on (default) they are removed using the laser data
alone; no building map is needed, and buildings newer than any map are caught too.

- **Buildings and other hard surfaces.** A roof is flat and returns one echo per pulse; a crown is
  rough and splits the pulse among branches. Cells at least 2 m high where the canopy model is
  within 0.2 m (RMS) of a plane over 2.5 m, and at least 70 % of first returns come from
  single-return pulses, form a mask. Patches under 10 m2 are dropped and the mask is grown by
  1.5 m so roof edges, ridges and rooftop equipment are covered. Tops inside it are removed.
- **Power lines and other thin objects.** A top is removed if at most 8 returns above 2 m support
  it and they are under 30 % of all first returns within its search window. A crown fills its
  window; a wire hardly does.
- Bridges (class 17) are kept out of the canopy model.

Tested on Lantmäteriet's 2026 scan (26A019, April) of about 100 km2 of farmland, forest, a small
town and coast in south-east Småland (643 000 tops), checked visually against Esri World Imagery
(0.5 m, April 2024):

| Setting | Tops before | Removed as hard surface | Removed as thin object |
|---|---|---|---|
| Industrial and residential area, 16 ha | 769 | 413 | 7 |
| Villa area with terraced houses, 12 ha | 336 | 148 | 0 |
| Power line across fields and forest | 1 585 | 0 | 278 (all tops on the line) |
| Dense forest, 6 ha | 1 640 | 0 | 2 |
| Trees in hedgerows and farmyards (roof tops not counted) | 516 | - | 1 |

With every second point dropped (about 1 first return per m2, like older scans) the result was
nearly the same: 5 of 1 542 forest tops removed, 242 of the wire tops, and the same roofs.

Known limits: a tree standing within about 1.5 m of a building may be removed with the roof; a
few tops remain on roof equipment and lattice pylons; and the thresholds were set on one scan in
one region.

## Utvärdera noggrannhet

| Parameter | Default | Notes |
|---|---|---|
| Detekterade träd | - | Point layer, e.g. the output above. Selection honoured |
| Referensträd | - | Point layer of reference trees, ideally stems mapped in the field. Selection honoured |
| Provytor | - | Polygons where all reference trees were mapped |
| Utdata, parade träd | `utvardering_trad` | Points of both sets with status (Träff, Missad, Felaktig, Utanför provytan), paired OID, distance and dx/dy |
| Utdata, mått per provyta | none | Optional table, one row per plot plus `Alla` |
| Sökradie (m) | 2 | 0.1-20 m |
| Höjdfält | `TRADHOJD` if present / none | Under Höjd |
| Största höjdskillnad (%) | 30 | Under Höjd. Applied only when both height fields are set |
| Provytans ID-fält | OID | Under Provytor |

Matching is one-to-one and minimises the total distance (Hungarian method,
`scipy.optimize.linear_sum_assignment`, solved per connected group of candidate pairs). It was
checked against brute force on 300 random cases: same number of pairs and same total distance in
every case.

- Only reference trees inside a plot count. Detected trees up to the search radius outside a plot
  may be paired with a reference tree inside; unpaired ones outside are not counted as false.
- Träffsäkerhet (recall) = paired reference trees / reference trees. Precision = paired detected
  trees inside the plots / detected trees inside the plots. Both are reported with a 95 % Wilson
  interval.
- Also reported: F1, stems per hectare, mean dx/dy and RMSE of paired positions, height bias and
  RMSE, recall per reference height class (under 10, 10-20, 20 m and over), counts of possible
  split crowns and merged trees, and results at 1, 1.5, 2 and 3 m search radius so the choice of
  radius is visible.
- Position error includes the error of the reference. GNSS under canopy is often off by 1-3 m;
  map stems from control points with a total station or tape and compass where possible.
- With a detection rate near 0.6, about 370 reference trees give a 95 % interval of about plus or
  minus 5 percentage points.

## Notes on the data

- Lantmäteriet's nominal density is 1-2 points per m2; counting all returns it is 2-4. The 2025-26
  scans fly lower and are somewhat denser.
- Many areas are scanned twice. The tool uses the newest scan; its date and leaf state are in the
  messages and the output.
- Southern Sweden is mostly scanned without leaves, so deciduous crowns return fewer first returns
  than conifers.
- Coverage is about 75 % of Sweden. The mountains are not included.

## Tests

`test_detektera_trad.py` runs both tools through `arcpy.ImportToolbox` on a local LAS file: two
detection runs in one process (typed extent, optional parameters empty and set), a check that every
crown matches its tree top by `TRAD_ID` and height and contains it, optionally a built-up area with
the building filter on and off, then an evaluation against synthetic reference trees with a known
expected precision. Set `TRAD_TEST_KEEP=1` to keep the outputs.
