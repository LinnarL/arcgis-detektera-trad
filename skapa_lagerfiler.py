# -*- coding: utf-8 -*-
"""
Skapar lagerfilerna som DetekteraTrad.pyt ger sina utdata:

    tradtoppar.lyrx  trädtopparna som ett litet, skarpt rött kors
    tradkronor.lyrx  kronorna med streckad kontur och ingen fyllning

Pro skriver filerna själv: ett lager sparas med SaveToLayerFile och bara
symbolen byts via CIM, så att formatet aldrig skrivs för hand. Körs med
ArcGIS Pro:s Python när symbolerna ska ändras:

    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" skapa_lagerfiler.py
"""

import os
import shutil
import tempfile

import arcpy

HERE = os.path.dirname(os.path.abspath(__file__))

CROSS_COLOR = [230, 0, 0, 100]      # RGB + alfa i procent
CROSS_SIZE_PT = 7
CROSS_WIDTH_PT = 1.0
CROWN_COLOR = [255, 176, 0, 100]
CROWN_WIDTH_PT = 1.0
CROWN_DASH_PT = [4, 2.5]


def cim(name):
    return arcpy.cim.CreateCIMObjectFromClassName(name, "V3")


def rgb(values):
    c = cim("CIMRGBColor")
    c.values = list(values)
    return c


def stroke(color, width, dashes=None):
    s = cim("CIMSolidStroke")
    s.color = rgb(color)
    s.width = width
    s.enable = True
    s.capStyle = "Butt"
    s.joinStyle = "Miter"
    s.miterLimit = 10
    if dashes:
        d = cim("CIMGeometricEffectDashes")
        d.dashTemplate = list(dashes)
        d.lineDashEnding = "NoConstraint"
        s.effects = [d]
    return s


def cross_symbol():
    """Ett kors av två linjer med raka ändar, utan kontur eller fyllning."""
    line = cim("CIMLineSymbol")
    line.symbolLayers = [stroke(CROSS_COLOR, CROSS_WIDTH_PT)]
    graphic = cim("CIMMarkerGraphic")
    graphic.geometry = {"paths": [[[-5, 0], [5, 0]], [[0, -5], [0, 5]]]}
    graphic.symbol = line
    marker = cim("CIMVectorMarker")
    marker.enable = True
    marker.size = CROSS_SIZE_PT
    marker.frame = {"xmin": -5, "ymin": -5, "xmax": 5, "ymax": 5}
    marker.markerGraphics = [graphic]
    marker.scaleSymbolsProportionally = True
    marker.respectFrame = True
    marker.anchorPointUnits = "Relative"
    sym = cim("CIMPointSymbol")
    sym.symbolLayers = [marker]
    sym.haloSize = 1
    sym.scaleX = 1
    sym.angleAlignment = "Display"
    return sym


def crown_symbol():
    """Streckad kontur, ingen fyllning."""
    sym = cim("CIMPolygonSymbol")
    sym.symbolLayers = [stroke(CROWN_COLOR, CROWN_WIDTH_PT, CROWN_DASH_PT)]
    return sym


def make(work, name, geom, symbol, layer_name):
    fc = arcpy.management.CreateFeatureclass(work, name, geom,
                                             spatial_reference=arcpy.SpatialReference(3006))[0]
    lyr = arcpy.management.MakeFeatureLayer(fc, layer_name)[0]
    out = os.path.join(HERE, name + ".lyrx")
    arcpy.management.SaveToLayerFile(lyr, out, "RELATIVE")
    lf = arcpy.mp.LayerFile(out)
    layer = lf.listLayers()[0]
    d = layer.getDefinition("V3")
    d.renderer.symbol.symbol = symbol
    layer.setDefinition(d)
    lf.save()
    # Läs tillbaka, så att ett fel i symbolen syns här och inte först i Pro.
    back = arcpy.mp.LayerFile(out).listLayers()[0].getDefinition("V3").renderer.symbol.symbol
    print(out, type(back).__name__, [type(s).__name__ for s in back.symbolLayers])


def main():
    work = tempfile.mkdtemp(prefix="lagerfiler_")
    gdb = arcpy.management.CreateFileGDB(work, "mall.gdb")[0]
    try:
        make(gdb, "tradtoppar", "POINT", cross_symbol(), "Trädtoppar")
        make(gdb, "tradkronor", "POLYGON", crown_symbol(), "Trädkronor")
    finally:
        arcpy.management.ClearWorkspaceCache(gdb)
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
