# -*- coding: utf-8 -*-
"""
DetekteraTrad.pyt

Två verktyg:

Detektera träd
    Hittar trädtoppar (krontoppar) i Lantmäteriets Laserdata Nedladdning, skog
    eller i lokala LAS/LAZ-filer och skriver dem som en punktfeatureklass med
    trädhöjd och markhöjd. Läget är krontoppens läge, inte en uppmätt
    stambas: vid 2-4 punkter per m² träffar nästan inga pulser stammen, och
    toppen kan ligga en bit från stambasen hos lutande träd och osymmetriska
    kronor.

Utvärdera noggrannhet
    Jämför detekterade träd med referensträd (t.ex. inmätta stammar) i
    provytor: parar ihop dem ett-till-ett (Ungerska metoden) inom en sökradie
    och redovisar träffsäkerhet (recall), precision, F1, stamantal per ha och
    lägesfel.

Datakälla (Detektera träd)
--------------------------
STAC-katalog (öppen, ingen inloggning):
    https://api.lantmateriet.se/stac-hojd/v1
    Samling: dsm-skoglig-copc (Laserdata Skog)
    Sök:     POST /search  {"collections": [...], "bbox": [...], "limit": 1000}
             Nästa sida via länken rel="next" (method + body). Sökningen
             svarar 429 efter ett tiotal snabba anrop i följd (inget
             Retry-After), så stora sidor och få anrop.

Varje item är en ruta i SWEREF 99 TM + RH 2000 (EPSG:5845), med en asset
"data" som pekar på en COPC-fil (.copc.laz). Sedan oktober 2026 finns varje
skanning i rutor på 2,5 x 2,5 km och dessutom i de äldre större rutorna med
samma punkter; _drop_retiled ser till att bara en av indelningarna läses.
properties: proj:bbox, pc:count, datetime, start_datetime, end_datetime,
skanningsomrade, punkttathet, flyghojd. Många områden är skannade två gånger;
verktyget tar den senaste skanningen per plats.

Nedladdning kräver OAuth2 (client credentials):
    POST https://apimanager.lantmateriet.se/oauth2/token
    Basic-auth med consumer key/secret, grant_type=client_credentials.
    Token gäller 3600 s. Nyckeln måste ha både API:et STAC-hojd och en
    beställning av Laserdata Nedladdning, skog, annars svarar dl1 med 403.

Hela rutor laddas aldrig ned. PDAL läser COPC-filerna direkt över HTTP och
hämtar bara punkterna inom områdets utbredning. Lokala LAS/LAZ-filer (t.ex.
utdrag skrivna av verktyget Laserdata höjdmodeller) läses i stället med
readers.las och klipps till varje block. Området delas i block som vart och
ett ryms i minnesbudgeten; blocken bearbetas i egna processer
(trad_worker.py), som returnerar topparna. Metoden beskrivs i trad_worker.py.

Verktygstips skrivs till DetekteraTrad.<Verktyg>.pyt.xml från TOOLTIPS nedan
när verktygslådan laddas.

Krav: ArcGIS Pro 3.x. arcpy, numpy, scipy, pdal och certifi ingår i
arcgispro-py3. Ingen licensnivå utöver Basic behövs.
"""

import base64
import datetime
import importlib
import importlib.util
import json
import math
import os
import queue
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from xml.sax.saxutils import escape

import arcpy
import numpy as np

# PDAL-delen ligger i trad_worker.py bredvid verktygslådan och körs i egna
# processer. reload: Pro håller modulen i minnet mellan körningar.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import trad_worker as _tw  # noqa: E402
importlib.reload(_tw)

# ── Konstanter ────────────────────────────────────────────────────────────────

TOOL_VERSION = "1.0"
REPO_NAME = "arcgis-detektera-trad"

STAC_SEARCH_URL = "https://api.lantmateriet.se/stac-hojd/v1/search"
TOKEN_URL = "https://apimanager.lantmateriet.se/oauth2/token"
COLLECTION = "dsm-skoglig-copc"
GEOTORGET_URL = "https://geotorget.lantmateriet.se/geodataprodukter/laserdata-nedladdning-skog-api"
TERMS_URL = ("https://www.lantmateriet.se/globalassets/geodata/geodataprodukter/"
             "anvandningsvillkor-for-laserdata-nedladdning-skog.pdf")

USER_AGENT = "arcgis-detektera-trad/" + TOOL_VERSION
HTTP_TIMEOUT = 60
HTTP_RETRIES = 3

# Token gäller 3600 s. Hämta en ny innan nästa block om den är äldre än så här.
TOKEN_MAX_AGE_S = 50 * 60

SWEREF99TM_WKID = 3006

# Varje block läses med denna marginal: största sökfönstret (radie 6 m), den
# gaussiska utjämningen och trianguleringen vid blockkanten ryms med god råge.
READ_MARGIN_M = 20.0
BLOCK_GRID_M = 1.0

# Minnestopp per läst punkt i en arbetsprocess och fast minne per process.
# Uppmätt: 269 MB topp för ett block med ca 0,8 miljoner lästa punkter, alltså
# ca 185 byte per punkt utöver ca 120 MB fast. Planeras med marginal.
PEAK_BYTES_PER_POINT = 260
WORKER_BASE_BYTES = 150e6

MIN_BLOCK_M = 250
MAX_BLOCK_M = 5000

# Automatiskt läge (tomma fält): hälften av det lediga arbetsminnet, men minst
# 2 GB kvar åt Pro och Windows, och en process per processortråd utom en.
AUTO_MEMORY_SHARE = 0.5
AUTO_MEMORY_RESERVE_GB = 2.0
MIN_MEMORY_GB = 0.5
MIN_GB_PER_WORKER = 0.75
MAX_WORKERS = 16

# Varna över denna yta: körningen blir lång.
LARGE_AREA_KM2 = 500

# Detektering. Fönstret är lmf-exemplet i lidR (diameter 0,1 * h + 3 m).
# Ingen av standardvärdena är kalibrerade mot svenska referensdata.
DEFAULT_MIN_HEIGHT = 4.0
DEFAULT_SIGMA = 0.5
DEFAULT_WIN_A = 0.1
DEFAULT_WIN_B = 3.0

AOI_POLYGONS = "Polygoner (lager eller ritade i kartan)"
AOI_EXTENT = "Utbredning (kartvy, lager eller koordinater)"
AOI_EMPTY = "Rita minst en polygon i kartan eller välj ett polygonlager."

LYRX_TOPS = "tradtoppar.lyrx"
LYRX_CROWNS = "tradkronor.lyrx"

SRC_STAC = "Lantmäteriets Laserdata Skog (hämtas över nätet)"
SRC_LOCAL = "Lokala LAS/LAZ-filer"

# Utvärdering
DEFAULT_SEARCH_RADIUS = 2.0
DEFAULT_MAX_DH_PCT = 30.0
SWEEP_RADII = (1.0, 1.5, 2.0, 3.0)
HEIGHT_CLASSES = ((0.0, 10.0), (10.0, 20.0), (20.0, 1e9))

STATUS_HIT = "Träff"
STATUS_MISS = "Missad"
STATUS_FALSE = "Felaktig"
STATUS_OUTSIDE = "Utanför provytan"

MONTHS = ["januari", "februari", "mars", "april", "maj", "juni", "juli", "augusti",
          "september", "oktober", "november", "december"]

DETECT_SUMMARY = (
    "Hittar trädtoppar (krontoppar) i laserdata och skriver dem som punkter med trädhöjd. "
    "Punkterna hämtas ur Lantmäteriets Laserdata Nedladdning, skog för området, eller läses "
    "ur lokala LAS/LAZ-filer. Läget är krontoppens, inte en uppmätt stambas. Vid "
    "Lantmäteriets punkttäthet hittas främst de dominerande träden; undertryckta träd under "
    "andra kronor syns sällan."
)
EVAL_SUMMARY = (
    "Jämför detekterade träd med referensträd i provytor. Varje referensträd paras med högst "
    "ett detekterat träd inom sökradien och tvärtom, så att summan av avstånden blir så liten "
    "som möjligt. Redovisar träffsäkerhet (andel referensträd som hittades), precision (andel "
    "detekterade träd som är riktiga), F1, stamantal per ha och lägesfel."
)

# Nyckeln kan sparas i Windows Autentiseringshanteraren. Samma post används av alla
# verktyg med Lantmäteriets nyckel, så namnen ska vara lika i verktygslådorna.
CREDENTIAL_SERVICE = "Lantmateriet API (ArcGIS-verktyg)"
CREDENTIAL_USER = "consumer"
LOGIN_NO_KEY = "Ange en consumer key från Lantmäteriets API-portal."
LOGIN_NO_SECRET = "Ange consumer secret."

# Verktygstips per parameter, visas i verktygsdialogen. Se _write_tool_metadata.
TOOLTIPS_DETECT = {
    "aoi_mode": (
        "Hur området avgränsas. 'Polygoner' använder polygoner ur ett lager eller polygoner "
        "som du ritar i kartan, och bara toppar inom själva polygonerna sparas. 'Utbredning' "
        "ger en rektangel: aktuell kartvy, utbredningen av ett lager, en ritad rektangel "
        "eller inskrivna koordinater."
    ),
    "aoi": (
        "Polygoner som avgränsar området. Välj ett polygonlager i listan, eller rita en eller "
        "flera polygoner i kartan med ritverktyget bredvid fältet (dubbelklicka för att "
        "avsluta en polygon). Från ett lager används alla objekt, eller bara de valda om "
        "lagret har ett urval. Alla polygoner slås ihop till ett område. Punkter läses inom "
        "polygonernas utbredning (bounding box); bara toppar inom polygonerna sparas."
    ),
    "aoi_extent": (
        "Rektangel som avgränsar området. I listan kan du välja kartvyns aktuella "
        "utbredning, eller ett lager för att använda dess utbredning. Du kan också rita en "
        "rektangel i kartan eller skriva in koordinater. Inskrivna koordinater tolkas i den "
        "aktiva kartans koordinatsystem; vilket som användes står i meddelandena."
    ),
    "source": (
        "Var punktmolnet kommer ifrån. 'Lantmäteriets Laserdata Skog' läser bara punkterna "
        "inom området direkt ur Lantmäteriets molnoptimerade filer, med den senaste "
        "skanningen per plats; det kräver en consumer key och secret under Inloggning. "
        "'Lokala LAS/LAZ-filer' läser filer du redan har, t.ex. utdrag från verktyget "
        "Laserdata höjdmodeller. De måste vara i SWEREF 99 TM och ha markpunkter klassade "
        "(klass 2)."
    ),
    "las_files": (
        "LAS- eller LAZ-filer som täcker området. Bara filer som överlappar området läses, och "
        "bara punkterna inom området används. Filerna får inte överlappa varandra med samma "
        "punkter, eftersom punkter då räknas två gånger. Krävs bara när punktmolnet kommer "
        "från lokala filer."
    ),
    "out_features": (
        "Punktfeatureklassen som skapas, en punkt per trädtopp, visad som ett litet rött kors. "
        "Fält: träd-ID (samma som i kronorna), trädhöjd (m över mark), "
        "markhöjd vid toppen (m ö.h., RH 2000), källfil eller ruta, skanningsområde och "
        "insamlingsperiod. Följer miljöinställningen Skriv över utdata. Koordinatsystemet är "
        "SWEREF 99 TM om inte miljön Utdatakoordinatsystem anger ett annat."
    ),
    "out_crowns": (
        "Valfri polygonfeatureklass med en krona per trädtopp, kopplad till punkten med fältet "
        "TRAD_ID. Fält: trädhöjd, kronans yta (m²) och diameter (m, cirkel med samma yta). "
        "Kronorna växer uppifrån i den utjämnade ythöjdmodellen från varje topp tills de möter "
        "grannens krona, sjunker under halva trädets höjd eller når 0,3 x trädhöjden från "
        "toppen, och breder aldrig ut sig över byggnader. Hela kronan behålls även där den når "
        "utanför området. Kronor runt toppar som togs bort som byggnad eller ledning tas också "
        "bort. Visas med streckad kontur utan fyllning. Töm fältet för att bara få topparna; "
        "det sparar ungefär en femtedel av tiden."
    ),
    "min_height": (
        "Lägsta trädhöjd i meter över mark. Lägre toppar sparas inte. Under ungefär 2-4 m blir "
        "buskar, block och brus i markmodellen toppar. Standard 4 m."
    ),
    "sigma": (
        "Hur mycket ythöjdmodellen jämnas ut innan topparna söks, som standardavvikelse i meter "
        "för ett gaussiskt filter. Större värde ger färre toppar: kronor med flera spetsar "
        "(tall, lövträd) delas mer sällan, men små träd tätt intill stora försvinner. 0 = "
        "ingen utjämning. Standard 0,5 m. Värdet är inte kalibrerat mot svenska referensdata; "
        "använd Utvärdera noggrannhet för att välja."
    ),
    "win_a": (
        "Sökfönstrets tillväxt med trädhöjden. En punkt i ythöjdmodellen är en topp om den är "
        "högst inom en cirkel med diametern A x höjden + B meter. Standard A = 0,1 och B = 3, "
        "alltså 4,5 m diameter för ett 15 m högt träd och 5,5 m för ett 25 m högt. Radien "
        "avrundas till hela halvmetrar, 1-6 m. Större fönster ger färre toppar."
    ),
    "win_b": (
        "Sökfönstrets diameter i meter vid höjden 0, se A. Standard 3 m."
    ),
    "filter_objects": (
        "Ta bort toppar som inte är träd, med hjälp av laserdata. Standard på.\n"
        "Byggnader och andra hårda ytor: tak är plana och ger en retur per puls, kronor är "
        "ojämna och delar pulsen i grenarna. Ytor minst 2 m höga som är plana inom 0,2 m och "
        "där minst 70 % av pulserna gav en enda retur blir en mask, som utvidgas 1,5 m så att "
        "takkanter och nockar kommer med. Träd som står tätt intill en byggnad kan försvinna.\n"
        "Ledningar och andra tunna föremål: en topp där högst 8 returer ligger över 2 m och "
        "de är under 30 % av alla returer inom sökfönstret. Stolpar och fackverksmaster ger "
        "fortfarande ofta en topp.\n"
        "Antalet borttagna toppar står i meddelandena."
    ),
    "consumer_key": (
        "Consumer key från Lantmäteriets API-portal. Nyckeln behöver både API:et STAC-hojd "
        "och en beställning av Laserdata Nedladdning, skog på Geotorget. Samma nyckel som till "
        "verktygen Laserdata höjdmodeller och Laserdata NRD. Behövs inte för lokala filer."
    ),
    "consumer_secret": (
        "Consumer secret som hör till nyckeln. Visas dold i dialogen och skrivs aldrig till "
        "meddelandena. Kan lämnas tomt om nyckeln är sparad."
    ),
    "save_login": (
        "Spara consumer key och secret i Windows Autentiseringshanteraren när körningen startar, "
        "krypterat för ditt Windows-konto och bara på den här datorn. Nyckeln sparas bara om "
        "Lantmäteriet godkänner den. Sedan fylls consumer key i av sig själv och secret kan "
        "lämnas tomt, i alla verktyg som använder Lantmäteriets nyckel. Fyll i och markera igen "
        "för att byta nyckel. Ta bort den under Autentiseringshanteraren, Windows-autentiseringsuppgifter."
    ),
    "memory_gb": (
        "Ungefär hur mycket arbetsminne verktyget får använda utöver ArcGIS Pro självt. "
        "Tomt = automatiskt: hälften av det lediga arbetsminnet när körningen startar, men "
        "minst 2 GB lämnas kvar åt Pro och Windows. Området delas i block som var för sig ryms "
        "i minnet, så ytan begränsas inte av minnet. Varje process behöver ungefär 260 byte "
        "per läst punkt plus 150 MB. Minst 0,5 GB."
    ),
    "workers": (
        "Antal block som bearbetas samtidigt, vart och ett i en egen process. Tomt = "
        "automatiskt: en process per processortråd utom en, högst 16, och högst en per "
        "0,75 GB i minnesbudgeten. Aldrig fler än det finns block. 1 till 32."
    ),
}

TOOLTIPS_EVAL = {
    "detected": (
        "Punktlager med de detekterade träden, t.ex. resultatet av Detektera träd. Har lagret "
        "ett urval används bara de valda punkterna."
    ),
    "reference": (
        "Punktlager med referensträden, helst stammar inmätta i fält. Positionsfelet i "
        "referensen ingår i lägesfelet som redovisas, så GNSS under krontak (ofta 1-3 m fel) "
        "ger ett sämre resultat än det verkliga. Har lagret ett urval används bara de valda."
    ),
    "plots": (
        "Polygoner som avgränsar de ytor där alla referensträd är inmätta. Bara referensträd "
        "inom en provyta räknas. Ett detekterat träd utanför provytan men inom sökradien från "
        "kanten får paras med ett referensträd innanför; blir det inte parat räknas det inte "
        "som felaktigt. Provytorna får inte ligga närmare varandra än sökradien."
    ),
    "out_features": (
        "Punktfeatureklass med både referensträd och detekterade träd, med status (Träff, "
        "Missad, Felaktig, Utanför provytan), provyta, det parade trädets objekt-ID, avstånd "
        "och förskjutning i x och y. Följer miljöinställningen Skriv över utdata."
    ),
    "out_table": (
        "Valfri tabell med måtten per provyta och för alla provytor tillsammans. Måtten står "
        "också i meddelandena."
    ),
    "search_radius": (
        "Största avstånd i meter mellan ett referensträd och ett detekterat träd för att de ska "
        "kunna paras. Standard 2 m. Krontoppen kan ligga 1-2 m från stambasen hos lutande "
        "träd, så en mindre radie mäter delvis den skillnaden. Meddelandena visar resultatet "
        "även för 1, 1,5, 2 och 3 m, så att valet syns. 0,1 till 20 m."
    ),
    "det_height_field": (
        "Fält med trädhöjd i meter i de detekterade träden. Tomt = TRADHOJD om fältet finns. "
        "Används för att inte para ihop träd med mycket olika höjd och för höjdfelet."
    ),
    "ref_height_field": (
        "Fält med trädhöjd i meter i referensträden. Tomt = ingen höjdkontroll och inga "
        "måttindelningar per höjdklass."
    ),
    "max_height_diff": (
        "Största höjdskillnad för att två träd ska kunna paras, i procent av referensträdets "
        "höjd. Används bara när båda höjdfälten är kända. Standard 30 %. Tomt = ingen gräns."
    ),
    "plot_id_field": (
        "Fält med provytans namn eller nummer. Tomt = provytans objekt-ID."
    ),
}


# =============================================================================
# Förlopp och format
# =============================================================================

def _fmt_duration(seconds):
    seconds = max(0, int(round(seconds)))
    if seconds < 60:
        return "{} s".format(seconds)
    if seconds < 3600:
        return "{} min".format(int(round(seconds / 60.0)))
    return "{} h {} min".format(seconds // 3600, (seconds % 3600) // 60)


def _fmt_count(n):
    return "{:,}".format(int(n)).replace(",", " ")


def _fmt_m(v):
    """Tal med decimalkomma och utan onödiga nollor: 0,15, 2, 1,85."""
    return "{:g}".format(round(float(v), 3)).replace(".", ",")


def _fmt_pct(v):
    return "{} %".format(_fmt_m(round(100.0 * v, 1))) if v is not None else "-"


class _Steps:
    """Numrerade steg i förloppsindikatorn och i meddelandena."""

    def __init__(self, total, messages):
        self.total = total
        self.messages = messages
        self.k = 0
        self.t0 = time.time()
        self.t_step = self.t0

    def next(self, text):
        self.k += 1
        self.t_step = time.time()
        label = "Steg {} av {}: {}".format(self.k, self.total, text)
        arcpy.SetProgressor("default", label)
        self.messages.addMessage(label)

    def label(self, text):
        arcpy.SetProgressorLabel("Steg {} av {}: {}".format(self.k, self.total, text))

    def done(self, text=None):
        msg = "    klart på {}".format(_fmt_duration(time.time() - self.t_step))
        self.messages.addMessage(msg + (". " + text if text else "."))


# =============================================================================
# HTTP
# =============================================================================

def _http(req):
    """urlopen med omförsök vid tillfälliga fel. Returnerar (status, body)."""
    delay = 2
    for attempt in range(HTTP_RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as exc:
            # HTTPError först: den är också URLError och OSError.
            if exc.code in (408, 429) or exc.code >= 500:
                if attempt < HTTP_RETRIES - 1:
                    time.sleep(delay)
                    delay *= 2
                    continue
            raise
        except (urllib.error.URLError, OSError):
            if attempt < HTTP_RETRIES - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise


def _get_token(key, secret):
    auth = base64.b64encode("{}:{}".format(key, secret).encode("utf-8")).decode("ascii")
    req = urllib.request.Request(
        TOKEN_URL,
        data=urllib.parse.urlencode({"grant_type": "client_credentials"}).encode("ascii"),
        headers={"Authorization": "Basic " + auth, "User-Agent": USER_AGENT,
                 "Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    try:
        _status, body = _http(req)
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 401):
            raise ValueError(
                "Lantmäteriet godkände inte consumer key/secret (HTTP {}). Kontrollera "
                "nyckeln i Lantmäteriets API-portal.".format(exc.code)
            )
        raise
    return json.loads(body.decode("utf-8"))["access_token"]


def _saved_login():
    """(key, secret) sparade i Windows Autentiseringshanteraren, annars tomma."""
    try:
        import keyring
        raw = keyring.get_password(CREDENTIAL_SERVICE, CREDENTIAL_USER)
        saved = json.loads(raw) if raw else {}
    except Exception:
        saved = {}
    return saved.get("key", ""), saved.get("secret", "")


def _login(p):
    """
    Key och secret: det som står i fälten, annars det sparade. Varje värde för
    sig, så att en ny secret fungerar även om Pro har tömt fältet för key.
    """
    saved_key, saved_secret = _saved_login()
    key = (p["consumer_key"].valueAsText or "").strip() or saved_key
    secret = (p["consumer_secret"].valueAsText or "").strip() or saved_secret
    if not key:
        raise ValueError(LOGIN_NO_KEY)
    if not secret:
        raise ValueError(LOGIN_NO_SECRET)
    return key, secret


def _check_login_fields(p):
    """Felmeddelanden i dialogen för det som varken är ifyllt eller sparat."""
    saved_key, saved_secret = _saved_login()
    if not (p["consumer_key"].valueAsText or "").strip() and not saved_key:
        p["consumer_key"].setErrorMessage(LOGIN_NO_KEY)
    if not (p["consumer_secret"].valueAsText or "").strip() and not saved_secret:
        p["consumer_secret"].setErrorMessage(LOGIN_NO_SECRET)


def _remember_login(key, secret, messages):
    """Spara nyckeln, men bara om Lantmäteriet godkänner den (annars ValueError)."""
    _get_token(key, secret)
    try:
        import keyring
        keyring.set_password(CREDENTIAL_SERVICE, CREDENTIAL_USER,
                             json.dumps({"key": key, "secret": secret}))
    except Exception as exc:
        messages.addWarningMessage("Nyckeln kunde inte sparas: {}".format(exc))
    else:
        messages.addMessage("Nyckeln är sparad i Windows Autentiseringshanteraren.")


def _check_access(url, token):
    """Hämta en byte av första filen, för ett begripligt fel i stället för PDAL:s."""
    req = urllib.request.Request(
        url, headers={"Authorization": "Bearer " + token, "Range": "bytes=0-0",
                      "User-Agent": USER_AGENT},
    )
    try:
        _http(req)
    except urllib.error.HTTPError as exc:
        if exc.code in (401, 403):
            raise ValueError(
                "Nyckeln saknar behörighet till punktmolnen (HTTP {}). Den behöver både "
                "API:et STAC-hojd i Lantmäteriets API-portal och en beställning av "
                "Laserdata Nedladdning, skog på Geotorget:\n{}".format(exc.code, GEOTORGET_URL)
            )
        raise


def _stac_search(bbox_wgs84):
    """Alla items i samlingen som skär bbox, över alla sidor."""
    body = {"collections": [COLLECTION], "bbox": list(bbox_wgs84), "limit": 1000}
    url, method = STAC_SEARCH_URL, "POST"
    items = []
    seen_pages = set()
    while url:
        data = json.dumps(body).encode("utf-8") if method == "POST" else None
        req = urllib.request.Request(
            url, data=data, method=method,
            headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
        )
        _status, raw = _http(req)
        page = json.loads(raw.decode("utf-8"))
        feats = page.get("features", [])
        # Skydd mot en pager som returnerar samma sida om och om igen.
        page_key = tuple(f.get("id") for f in feats)
        if page_key in seen_pages:
            break
        seen_pages.add(page_key)
        items.extend(feats)

        nxt = next((l for l in page.get("links", []) if l.get("rel") == "next"), None)
        if not nxt or not feats:
            break
        url = nxt["href"]
        method = nxt.get("method", "GET").upper()
        body = nxt.get("body", body)
    return items


class _Token:
    """Token som förnyas innan den går ut."""

    def __init__(self, key, secret):
        self._key, self._secret = key, secret
        self._token = _get_token(key, secret)
        self._time = time.time()

    def get(self):
        if time.time() - self._time > TOKEN_MAX_AGE_S:
            self._token = _get_token(self._key, self._secret)
            self._time = time.time()
        return self._token


# =============================================================================
# Geometri
# =============================================================================

def _tm():
    return arcpy.SpatialReference(SWEREF99TM_WKID)


def _rect_polygon(x0, y0, x1, y1):
    return arcpy.Polygon(arcpy.Array([arcpy.Point(x0, y0), arcpy.Point(x0, y1),
                                      arcpy.Point(x1, y1), arcpy.Point(x1, y0)]), _tm())


def _polygon_schema():
    """
    Tom polygonfeatureklass i SWEREF 99 TM, standardvärde för Feature Set-
    parametern så att dialogens ritverktyg ritar polygoner. Unikt namn i memory:
    memory-arbetsytan delas av hela Pro-sessionen, och Exists/Delete på ett fast
    namn har gett SQL-fel där.
    """
    name = "aoi_schema_{}".format(uuid.uuid4().hex[:12])
    arcpy.management.CreateFeatureclass("memory", name, "POLYGON", spatial_reference=_tm())
    return "memory/" + name


def _is_layer(value):
    """
    Sant för ett Layer-objekt (ett lager valt i listan). arcpy.mp har ingen
    Layer-klass att jämföra med (AttributeError), klassen ligger i arcpy._mp.
    """
    return type(value).__name__ == "Layer"


def _has_features(value):
    """
    True om värdet har minst ett objekt (ett lagers urval räknas), False om det
    är tomt, None om det inte går att läsa. Läser bara första raden.
    """
    try:
        with arcpy.da.SearchCursor(value, ["OID@"]) as cur:
            return next(iter(cur), None) is not None
    except Exception:
        return None


def _aoi_geometry(value, messages=None):
    """
    Alla polygoner i parametervärdet, sammanslagna, i SWEREF 99 TM.

    value är det som Feature Set-parametern (GPFeatureRecordSetLayer) ger i
    execute(): ett Layer-objekt när ett lager valts (då läses bara urvalet om
    det finns ett), ett record set med polygoner ritade i kartan, eller en
    featureklass som sökväg. SearchCursor och Describe fungerar på alla tre.
    """
    if _has_features(value) is False:
        raise ValueError(AOI_EMPTY)
    # Kontrollera koordinatsystemet innan något omprojiceras: utan det skulle
    # polygonerna tyst läsas som om de redan vore i SWEREF 99 TM.
    sr_in = arcpy.Describe(value).spatialReference
    if sr_in is None or not (sr_in.factoryCode or sr_in.exportToString()):
        raise ValueError("Intresseområdet saknar koordinatsystem.")
    geom = None
    n = 0
    with arcpy.da.SearchCursor(value, ["SHAPE@"], spatial_reference=_tm()) as cur:
        for (shape,) in cur:
            if shape is None or shape.area <= 0:
                continue
            n += 1
            geom = shape if geom is None else geom.union(shape)
    if geom is None:
        raise ValueError("Intresseområdet innehåller inga polygoner med yta.")
    if messages is not None:
        messages.addMessage("Intresseområde: {} polygon(er) i {}, sammanlagt {} ha.".format(
            n, sr_in.name, _fmt_m(round(geom.area / 1e4, 1))))
    return geom


def _active_map_sr():
    try:
        m = arcpy.mp.ArcGISProject("CURRENT").activeMap
        if m is not None and m.spatialReference is not None:
            return m.spatialReference, m.name
    except Exception:
        pass
    return None, None


def _densified_rect(xmin, ymin, xmax, ymax, sr):
    """Rektangel med förtätade kanter, omprojicerad till SWEREF 99 TM."""
    n = 16
    pts = ([(xmin + (xmax - xmin) * i / n, ymin) for i in range(n)]
           + [(xmax, ymin + (ymax - ymin) * i / n) for i in range(n)]
           + [(xmax - (xmax - xmin) * i / n, ymax) for i in range(n)]
           + [(xmin, ymax - (ymax - ymin) * i / n) for i in range(n)])
    poly = arcpy.Polygon(arcpy.Array([arcpy.Point(x, y) for x, y in pts]), sr)
    if sr.factoryCode != SWEREF99TM_WKID:
        poly = poly.projectAs(_tm())
    return poly


def _extent_polygon(value, text):
    """
    En GPExtent-parameter som polygon i SWEREF 99 TM. Returnerar (polygon, text om
    vilket koordinatsystem som användes).

    .value är ett geoprocessing-extentobjekt, inte arcpy.Extent. Hörnen är
    vanliga tal. Koordinatsystemet följer bara med när utbredningen kommer från
    ett lager eller en datakälla, och då bara som WKT2 i valueAsText efter de
    fyra talen. Inskrivna koordinater har inget; de tolkas i den aktiva kartans
    koordinatsystem.
    """
    xmin, ymin, xmax, ymax = (float(value.XMin), float(value.YMin),
                              float(value.XMax), float(value.YMax))
    if not (xmax > xmin and ymax > ymin):
        raise ValueError("Utbredningen har ingen yta.")

    sr = None
    parts = (text or "").split(" ", 4)
    if len(parts) == 5 and parts[4].strip():
        sr = arcpy.SpatialReference()
        try:
            sr.loadFromString(parts[4].strip())
        except Exception:
            sr = None
    if sr is not None and (sr.factoryCode or sr.exportToString()):
        source = "utbredningens eget koordinatsystem"
    else:
        sr, map_name = _active_map_sr()
        if sr is not None and (sr.factoryCode or sr.exportToString()):
            source = "den aktiva kartans koordinatsystem ({})".format(map_name)
        else:
            sr = _tm()
            source = "ingen aktiv karta, så SWEREF 99 TM antas"
    poly = _densified_rect(xmin, ymin, xmax, ymax, sr)
    if poly is None or poly.area <= 0:
        raise ValueError("Utbredningen kunde inte omvandlas till SWEREF 99 TM.")
    return poly, "Utbredningen tolkas i {}: {}.".format(sr.name, source)


def _grid(extent, cell):
    """Rutnät justerat till jämna multiplar av cellstorleken."""
    x0 = math.floor(extent.XMin / cell) * cell
    y0 = math.floor(extent.YMin / cell) * cell
    width = int(math.ceil((extent.XMax - x0) / cell))
    height = int(math.ceil((extent.YMax - y0) / cell))
    return {"resolution": cell, "origin_x": x0, "origin_y": y0,
            "width": width, "height": height}


# =============================================================================
# Rutor i katalogen
# =============================================================================

def _actual_density(tile):
    """Punkter per m² med alla returer, ur punktantalet och rutans yta."""
    x0, y0, x1, y1 = tile["bbox"]
    area = (x1 - x0) * (y1 - y0)
    return tile["count"] / area if tile["count"] and area > 0 else None


def _tile_from_item(it):
    props = it.get("properties", {})
    asset = it.get("assets", {}).get("data")
    pb = props.get("proj:bbox") or (asset or {}).get("proj:bbox")
    if not asset or not pb:
        return None
    dt = props.get("datetime") or ""
    return {"id": it["id"], "href": asset["href"], "path": None, "datetime": dt,
            "bbox": tuple(pb[:4]),
            "count": props.get("pc:count") or 0,
            "start": props.get("start_datetime") or dt,
            "end": props.get("end_datetime") or dt,
            "area": props.get("skanningsomrade") or ""}


def _last_day(tile):
    """Sista flygdagen. end_datetime är midnatt efter sista dagen, så backa en dag då."""
    start, last = tile["start"][:10], tile["end"][:10]
    if last > start and tile["end"][11:19] == "00:00:00":
        try:
            last = (datetime.date.fromisoformat(last) - datetime.timedelta(days=1)).isoformat()
        except ValueError:
            pass
    return last or start


def _capture_period(tile):
    """'2021-03-07 - 2021-04-01', ett enda datum, eller 'okänt'."""
    start = tile["start"][:10]
    if not start:
        return "okänt"
    last = _last_day(tile)
    return start if last in ("", start) else "{} - {}".format(start, last)


def _months(tiles):
    """Alla månader (år, månad) som rutornas insamling spänner över."""
    out = set()
    for t in tiles:
        try:
            d0 = datetime.date.fromisoformat(t["start"][:10])
            d1 = datetime.date.fromisoformat(_last_day(t))
        except ValueError:
            continue
        y, m = d0.year, d0.month
        while (y, m) <= (d1.year, d1.month):
            out.add((y, m))
            y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return sorted(out)


def _season(months):
    """Lövläge ur insamlingsmånaderna: 'med löv', 'utan löv' eller 'lövläge osäkert'."""
    if not months:
        return "lövläge okänt"
    kinds = set()
    for _y, m in months:
        kinds.add("med" if 6 <= m <= 8 else "utan" if m in (11, 12, 1, 2, 3, 4) else "osäkert")
    if kinds == {"med"}:
        return "med löv"
    if kinds == {"utan"}:
        return "utan löv"
    return "lövläge osäkert"


def _overlaps(aoi, bbox):
    """
    Sant om rutan och aoi har gemensam yta. disjoint() räcker inte: en ruta som
    bara nuddar områdets kant skulle då räknas med.
    """
    rect = _rect_polygon(*bbox)
    return not aoi.disjoint(rect) and aoi.intersect(rect, 4).area > 1.0


def _scan_of(tile):
    return (tile["area"] or tile["id"].split("-")[0]).upper()


def _rect_area(tile):
    x0, y0, x1, y1 = tile["bbox"]
    return (x1 - x0) * (y1 - y0)


def _contains(outer, inner):
    """Sant om bbox inner ligger helt inom bbox outer, på en halv meter när."""
    return (outer[0] <= inner[0] + 0.5 and outer[1] <= inner[1] + 0.5
            and outer[2] >= inner[2] - 0.5 and outer[3] >= inner[3] - 0.5)


def _union(tiles):
    union = None
    for t in tiles:
        r = _rect_polygon(*t["bbox"])
        union = r if union is None else union.union(r)
    return union


def _covered(geom, union):
    """Sant om geom ligger inom union, på en kvadratmeter när."""
    return union is not None and geom.difference(union).area < 1.0


def _children_complete(parent):
    """
    Sant om de mindre rutorna inom parent, ur samma skanning, tillsammans har
    lika många punkter som den. Söker över hela rutan, eftersom sökningen för
    området bara ger de rutor som skär området.
    """
    x0, y0, x1, y1 = parent["bbox"]
    wgs = _rect_polygon(x0, y0, x1, y1).densify("DISTANCE", 500, 1).projectAs(
        arcpy.SpatialReference(4326)).extent
    counts = {}
    for it in _stac_search((wgs.XMin, wgs.YMin, wgs.XMax, wgs.YMax)):
        t = _tile_from_item(it)
        if (t and t["id"] != parent["id"] and _scan_of(t) == _scan_of(parent)
                and _rect_area(t) < _rect_area(parent) and _contains(parent["bbox"], t["bbox"])):
            counts[t["id"]] = t["count"]
    return bool(counts) and sum(counts.values()) == parent["count"]


def _drop_retiled(tiles, aoi):
    """
    Samma punkter får inte läsas två gånger. I oktober 2026 delade Lantmäteriet
    skanningarna i rutor på 2,5 x 2,5 km men behöll de större rutorna i
    katalogen, med samma punkter. En ruta som rymmer mindre rutor ur samma
    skanning ersätts av dem, om de täcker rutans del av området eller har alla
    dess punkter. Annars behålls den stora rutan och de mindre inom den tas bort.
    """
    by_scan = {}
    for t in tiles:
        by_scan.setdefault(_scan_of(t), []).append(t)
    drop = set()
    for scan_tiles in by_scan.values():
        # Största först, så att en ruta som behålls tar bort allt inom sig.
        for big in sorted(scan_tiles, key=_rect_area, reverse=True):
            if big["id"] in drop:
                continue
            kids = [t for t in scan_tiles
                    if t["id"] not in drop and _rect_area(t) < _rect_area(big)
                    and _contains(big["bbox"], t["bbox"])]
            if not kids:
                continue
            part = aoi.intersect(_rect_polygon(*big["bbox"]), 4)
            if part.area <= 0 or _covered(part, _union(kids)) or _children_complete(big):
                drop.add(big["id"])
            else:
                drop.update(k["id"] for k in kids)
    return [t for t in tiles if t["id"] not in drop]


def _newest_per_place(tiles, aoi):
    """
    Den senaste insamlingen per plats: en ruta tas bort om dess del av området
    redan täcks av rutor ur nyare skanningar.
    """
    kept, union = [], None
    for t in sorted(tiles, key=lambda t: (t["datetime"], t["id"]), reverse=True):
        rect = _rect_polygon(*t["bbox"])
        part = aoi.intersect(rect, 4)
        if _covered(part if part.area > 0 else rect, union):
            continue
        kept.append(t)
        union = rect if union is None else union.union(rect)
    return sorted(kept, key=lambda t: t["id"])


def _pick_tiles(items, aoi):
    """Rutor som skär själva polygonen, den senaste insamlingen per plats."""
    tiles = []
    for it in items:
        t = _tile_from_item(it)
        if t is not None and _overlaps(aoi, t["bbox"]):
            tiles.append(t)
    return _newest_per_place(_drop_retiled(tiles, aoi), aoi)


def _local_tiles(headers, aoi, messages):
    """
    Lokala filer som rutor: de som överlappar området, med utbredning och
    punktantal ur filhuvudet. Kontrollerar att filerna är i SWEREF 99 TM.
    """
    tiles = []
    for hd in headers:
        name = os.path.basename(hd["path"])
        wkt = hd.get("wkt") or ""
        if wkt:
            sr = arcpy.SpatialReference()
            try:
                sr.loadFromString(wkt)
            except Exception:
                sr = None
            if sr is None or (sr.factoryCode != SWEREF99TM_WKID
                              and "SWEREF99_TM" not in (sr.name or "").upper()
                              and "SWEREF99 TM" not in (sr.name or "").upper()):
                raise ValueError("{} är inte i SWEREF 99 TM ({}). Verktyget kräver SWEREF 99 TM."
                                 .format(name, sr.name if sr is not None else "okänt"))
        else:
            messages.addWarningMessage(
                "{} saknar koordinatsystem i filhuvudet; SWEREF 99 TM antas.".format(name))
        x0, y0, x1, y1 = hd["bbox"]
        if not (x1 > x0 and y1 > y0) or not _overlaps(aoi, (x0, y0, x1, y1)):
            continue
        tiles.append({"id": name, "href": None, "path": hd["path"], "datetime": "",
                      "bbox": (x0, y0, x1, y1), "count": hd["count"], "start": "", "end": "",
                      "area": ""})
    return tiles


# =============================================================================
# Block
# =============================================================================

def _available_memory_gb():
    """Ledigt arbetsminne i GB, eller None. psutil ingår i arcgispro-py3."""
    try:
        import psutil
        return psutil.virtual_memory().available / 1e9
    except Exception:
        return None


def _cpu_threads():
    """Processortrådar som processen får använda (logiska kärnor)."""
    n = getattr(os, "process_cpu_count", os.cpu_count)()
    return max(1, n or 1)


def _resolve_resources(memory_gb, workers):
    """
    Minnesbudget och antal processer, där None betyder automatiskt. Returnerar
    (minne i GB, processer, text om vad som valdes automatiskt eller "").
    """
    notes = []
    if memory_gb is None:
        avail = _available_memory_gb()
        if avail is None:
            memory_gb = 4.0
            notes.append("minnesbudget 4 GB (ledigt minne kunde inte läsas)")
        else:
            memory_gb = max(MIN_MEMORY_GB, min(AUTO_MEMORY_SHARE * avail,
                                               avail - AUTO_MEMORY_RESERVE_GB))
            notes.append("minnesbudget {} GB av {} GB ledigt".format(
                _fmt_m(round(memory_gb, 1)), _fmt_m(round(avail, 1))))
    memory_gb = max(MIN_MEMORY_GB, float(memory_gb))
    if workers is None:
        threads = _cpu_threads()
        workers = max(1, min(threads - 1, MAX_WORKERS, int(memory_gb / MIN_GB_PER_WORKER)))
        notes.append("{} processer ({} processortrådar)".format(workers, threads))
    workers = max(1, min(32, int(workers)))
    note = "Automatiskt: {}.".format(", ".join(notes)) if notes else ""
    fit = max(1, int(memory_gb / MIN_GB_PER_WORKER))
    if workers > fit:
        note = (note + " " if note else "") + (
            "{} processer ryms inte i {} GB; använder {} (minst {} GB per process).".format(
                workers, _fmt_m(round(memory_gb, 1)), fit, _fmt_m(MIN_GB_PER_WORKER)))
        workers = fit
    return memory_gb, workers, note


def _plan_blocks(aoi, grid, tiles, memory_gb, workers):
    """
    Dela områdets rutnät i block som vart och ett ryms i minnesbudgeten delat
    med antalet samtidiga block. Blockstorleken räknas från den tätaste rutans
    punkttäthet. Blocken görs också så små att de blir minst lika många som
    processerna (men inte under MIN_BLOCK_M). Block som inte skär polygonen
    hoppas över.

    Returnerar (block, blocksida i m). Varje block har sina hörn, läsområdet
    med marginal och de rutor det behöver läsa.
    """
    cell = grid["resolution"]
    margin = READ_MARGIN_M
    density = max(_actual_density(t) or 0.0 for t in tiles) or 1.0
    per_worker = memory_gb * 1e9 / max(1, workers) - WORKER_BASE_BYTES
    points_per_block = max(1e5, per_worker / PEAK_BYTES_PER_POINT)
    side = math.sqrt(points_per_block / density) - 2 * margin
    area = grid["width"] * grid["height"] * cell * cell
    side = min(side, math.sqrt(area / max(1, workers)))
    side = min(MAX_BLOCK_M, max(MIN_BLOCK_M, side))
    side_cells = max(1, int(side // cell))

    blocks = []
    n_cols = int(math.ceil(grid["width"] / side_cells))
    n_rows = int(math.ceil(grid["height"] / side_cells))
    for ry in range(n_rows):
        for cx in range(n_cols):
            c0, r0 = cx * side_cells, ry * side_cells
            w = min(side_cells, grid["width"] - c0)
            h = min(side_cells, grid["height"] - r0)
            x0 = grid["origin_x"] + c0 * cell
            y0 = grid["origin_y"] + r0 * cell
            x1, y1 = x0 + w * cell, y0 + h * cell
            if aoi.disjoint(_rect_polygon(x0, y0, x1, y1)):
                continue
            mx0, my0 = x0 - margin, y0 - margin
            mx1, my1 = x1 + margin, y1 + margin
            reads = []
            expected = 0.0
            for t in tiles:
                tx0, ty0, tx1, ty1 = t["bbox"]
                bx0, bx1 = max(tx0, mx0), min(tx1, mx1)
                by0, by1 = max(ty0, my0), min(ty1, my1)
                if bx1 <= bx0 or by1 <= by0:
                    continue
                bounds = "([{:.2f},{:.2f}],[{:.2f},{:.2f}])".format(bx0, bx1, by0, by1)
                reads.append((t, bounds))
                expected += t["count"] * (bx1 - bx0) * (by1 - by0) / ((tx1 - tx0) * (ty1 - ty0))
            if not reads:
                continue
            blocks.append({
                # Rad räknas från norr, som man läser en karta.
                "row": n_rows - ry, "col": cx + 1,
                "x0": x0, "y0": y0, "x1": x1, "y1": y1,
                "read_bbox": [mx0, my0, mx1, my1],
                "reads": reads, "expected": expected,
            })
    return blocks, side_cells * cell


# =============================================================================
# Arbetsprocesser
# =============================================================================

def _python_exe():
    """python.exe i den aktiva miljön. I Pro är sys.executable ArcGISPro.exe."""
    for cand in (os.path.join(sys.exec_prefix, "python.exe"), sys.executable):
        if cand and os.path.basename(cand).lower() == "python.exe" and os.path.isfile(cand):
            return cand
    raise ValueError("Hittar inte python.exe i den aktiva miljön ({}).".format(sys.exec_prefix))


class _WorkerPool:
    """
    Ett fast antal trad_worker-processer som tar ett block i taget.

    Varför processer och inte trådar: flera PDAL-pipelines i trådar i samma
    process som arcpy kraschade processen (access violation, uppmätt i Pro 3.6).
    Processerna startas en gång och återanvänds. En läsartråd per process
    väntar på svar; trådarna gör bara I/O mot rören, aldrig PDAL eller arcpy.
    """

    def __init__(self, n, workdir):
        self._results = queue.Queue()
        self._procs = []
        self._logs = []
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        exe = _python_exe()
        script = os.path.join(_HERE, "trad_worker.py")
        env = dict(os.environ, PYTHONIOENCODING="utf-8")
        for i in range(n):
            log_path = os.path.join(workdir, "worker_{}.log".format(i))
            log = open(log_path, "w", encoding="utf-8", errors="replace")
            proc = subprocess.Popen(
                [exe, "-u", script], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=log, encoding="utf-8", errors="replace", creationflags=flags,
                env=env, cwd=workdir)
            self._procs.append(proc)
            self._logs.append((log, log_path))
            threading.Thread(target=self._reader, args=(i, proc), daemon=True).start()

    def _reader(self, i, proc):
        for line in proc.stdout:
            if line.startswith(_tw.RESULT_PREFIX):
                self._results.put((i, json.loads(line[len(_tw.RESULT_PREFIX):])))
        self._results.put((i, None))  # processen har avslutats

    def submit(self, i, task):
        self._procs[i].stdin.write(json.dumps(task) + "\n")
        self._procs[i].stdin.flush()

    def get(self, timeout):
        """(processindex, svar) eller None vid timeout. Svar None = processen dog."""
        try:
            return self._results.get(timeout=timeout)
        except queue.Empty:
            return None

    def call(self, task):
        """Kör en uppgift i process 0 och vänta på svaret (bara innan blocken startar)."""
        self.submit(0, task)
        while True:
            got = self.get(timeout=1.0)
            if getattr(arcpy.env, "isCancelled", False):
                raise ValueError("Avbrutet av användaren.")
            if got is None:
                continue
            _i, reply = got
            if reply is None:
                raise RuntimeError("Arbetsprocessen avslutades oväntat:\n" + self.log_tail(0))
            if not reply["ok"]:
                raise RuntimeError(reply["error"])
            return reply["result"]

    def log_tail(self, i, n=15):
        log, path = self._logs[i]
        log.flush()
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                return "".join(fh.readlines()[-n:])
        except OSError:
            return ""

    def close(self, kill=False):
        for proc in self._procs:
            try:
                if kill:
                    proc.kill()
                else:
                    proc.stdin.close()
            except OSError:
                pass
        for proc in self._procs:
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
        for log, _path in self._logs:
            log.close()


def _block_task(blk, token, params, workdir):
    """Blocket som JSON-uppgift till en arbetsprocess (bara det den behöver)."""
    return {
        "op": "block",
        "block": {k: blk[k] for k in ("row", "col", "x0", "y0", "x1", "y1", "read_bbox")}
        | {"reads": [{"href": t["href"], "path": t["path"], "id": t["id"], "bounds": bounds}
                     for t, bounds in blk["reads"]]},
        "token": token, "params": params, "workdir": workdir,
    }


# =============================================================================
# Utdata och metadata
# =============================================================================

def _is_network_path(path):
    """Sant om sökvägen ligger på en nätverksresurs (UNC eller mappad nätverksenhet)."""
    norm = (path or "").replace("/", "\\")
    if norm.startswith("\\\\"):
        return True
    drive = os.path.splitdrive(norm)[0]
    if len(drive) == 2 and drive[1] == ":":
        try:
            import ctypes
            return ctypes.windll.kernel32.GetDriveTypeW(drive + "\\") == 4   # DRIVE_REMOTE
        except Exception:
            return False
    return False


def _local_scratch_folder(messages):
    """
    Mappen som arbetsmappen läggs i: arcpy.env.scratchFolder om den är lokal.
    Pekar den på en nätverksresurs eller en synkad mapp används systemets
    temp-mapp. Utdata byggs här och skrivs till målet med en enda
    ExportFeatures: på en filserver kostar varje geobehandlingsanrop 16-27 s.
    """
    folder = arcpy.env.scratchFolder
    synced = ("onedrive", "sharepoint", "dropbox", "google drive")
    if folder and not _is_network_path(folder) and \
            not any(h in folder.lower() for h in synced):
        return folder
    fallback = tempfile.gettempdir()
    messages.addMessage(
        "    Scratch Workspace ligger på en nätverksresurs eller i en synkad mapp ({}). "
        "Arbetsfilerna läggs i {} i stället.".format(folder, fallback))
    return fallback


def _default_gdb():
    try:
        gdb = arcpy.mp.ArcGISProject("CURRENT").defaultGeodatabase
        if gdb:
            return gdb
    except Exception:
        pass
    return None


def _norm_path(path):
    return os.path.normcase(os.path.normpath(str(path)))


def _release_from_maps(paths, messages):
    """
    Ta bort lager i projektets kartor som visar data som ska skrivas över. Ett
    lager håller ett schemalås, och att skriva över misslyckas då med ERROR 000464.
    """
    targets = {_norm_path(p) for p in paths}
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
    except Exception:
        return
    removed = 0
    for m in aprx.listMaps():
        for lyr in m.listLayers():
            try:
                if (lyr.supports("DATASOURCE") and not lyr.isGroupLayer
                        and _norm_path(lyr.dataSource) in targets):
                    m.removeLayer(lyr)
                    removed += 1
            except Exception:
                continue
    if removed:
        messages.addMessage("    Tog bort {} lager ur kartan som visade data som skrivs "
                            "över.".format(removed))


def _out_sr():
    """Utdatakoordinatsystem: miljön om den är satt, annars SWEREF 99 TM."""
    sr = arcpy.env.outputCoordinateSystem
    if sr is not None and (sr.factoryCode or sr.exportToString()):
        return sr
    return _tm()


def _memory_name(stem):
    """
    Unikt namn i memory för mellanresultat. En filgeodatabas i arbetsmappen
    gick inte att ta bort: en tabell förblev låst tills processen avslutades,
    även efter ClearWorkspaceCache (uppmätt), så varje körning lämnade en mapp
    kvar. memory delas av hela Pro-sessionen, därav uuid-namnet.
    """
    return "memory/{}_{}".format(stem, uuid.uuid4().hex[:12])


def _delete_quietly(*paths):
    """Delete på mellanresultat i memory; kan misslyckas sporadiskt i Pro (se CLAUDE.md)."""
    for path in paths:
        if not path:
            continue
        try:
            arcpy.management.Delete(path)
        except Exception:
            pass


def _remove_workdir(workdir):
    """Ta bort arbetsmappen (arbetsprocessernas loggar och rastren de skrev)."""
    shutil.rmtree(workdir, ignore_errors=True)


def _export(src, dst, out_sr, messages, where=None):
    """
    Skriv det lokalt byggda resultatet till målet med en enda ExportFeatures.
    EnvManager hindrar att miljöerna klipper eller projicerar om något annat
    än det som är avsett.
    """
    if arcpy.Exists(dst):
        _release_from_maps([dst], messages)
    with arcpy.EnvManager(outputCoordinateSystem=out_sr, extent="MAXOF"):
        arcpy.conversion.ExportFeatures(src, dst, where_clause=where or "")


def _write_tool_metadata(tool_cls, toolbox_alias, tooltips, summary):
    """
    Skriv verktygets metadatafil med parameterförklaringar.

    Pro läser verktygstipsen i dialogen från <verktygslåda>.<verktyg>.pyt.xml
    (elementet dialogReference per parameter). Det finns inget attribut på
    arcpy.Parameter för detta. Filen skrivs bara om innehållet har ändrats.
    """
    toolbox = os.path.splitext(os.path.basename(__file__))[0]
    path = os.path.join(_HERE, "{}.{}.pyt.xml".format(toolbox, tool_cls.__name__))

    def html(text):
        body = escape(text).replace("\n", "</SPAN></P><P><SPAN>")
        return escape('<DIV STYLE="text-align:Left;"><P><SPAN>{}</SPAN></P></DIV>'.format(body))

    tool = tool_cls()
    params = []
    for p in tool.getParameterInfo():
        tip = tooltips.get(p.name)
        if not tip:
            continue
        params.append(
            '<param name="{n}" displayname="{d}" type="{t}" direction="{r}">'
            "<dialogReference>{h}</dialogReference>"
            "<pythonReference>{h}</pythonReference></param>".format(
                n=p.name, d=escape(p.displayName, {'"': "&quot;"}),
                t=p.parameterType, r=p.direction, h=html(tip))
        )
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<metadata xml:lang="sv"><Esri><ArcGISFormat>1.0</ArcGISFormat></Esri>'
        '<tool name="{name}" displayname="{label}" toolboxalias="{alias}" xmlns="">'
        "<parameters>{params}</parameters><summary>{summary}</summary></tool>"
        "<dataIdInfo><idCitation><resTitle>{label}</resTitle></idCitation>"
        "<idAbs>{summary}</idAbs></dataIdInfo></metadata>\n"
    ).format(name=tool_cls.__name__, label=escape(tool.label), alias=toolbox_alias,
             params="".join(params), summary=html(summary))

    try:
        with open(path, encoding="utf-8") as fh:
            if fh.read() == xml:
                return
    except OSError:
        pass
    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(xml)
    except OSError:
        # Skrivskyddad plats: verktyget fungerar ändå, bara utan verktygstips.
        pass


def _set_symbology(parameter, lyrx):
    """
    Lagerfil som Pro använder när utdata läggs till i kartan. Saknas filen
    (t.ex. om bara .pyt kopierats) får lagret Pro:s standardsymbol.
    """
    path = os.path.join(_HERE, lyrx)
    if os.path.isfile(path):
        parameter.symbology = path


def _unique_output(name):
    """Ett ledigt namn i projektets geodatabas, som standardvärde för utdata."""
    gdb = _default_gdb()
    if not gdb:
        return None
    try:
        return arcpy.CreateUniqueName(name, gdb)
    except Exception:
        return os.path.join(gdb, name)


# =============================================================================
# Utvärdering: parning och mått (numpy och scipy, inget arcpy)
# =============================================================================

def _match(ref_xy, det_xy, radius, ref_h=None, det_h=None, max_dh_pct=None):
    """
    Ett-till-ett-parning mellan referensträd och detekterade träd inom radius,
    som minimerar summan av avstånden (Ungerska metoden,
    scipy.optimize.linear_sum_assignment). Med höjder och max_dh_pct får två
    träd bara paras om |h_det - h_ref| <= max_dh_pct % av h_ref.

    Löses per sammanhängande komponent i grafen av möjliga par, så att stora
    provytor inte kräver en full avståndsmatris.
    Returnerar en lista (i_ref, j_det, avstånd).
    """
    from scipy.optimize import linear_sum_assignment
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    n_r, n_d = len(ref_xy), len(det_xy)
    if n_r == 0 or n_d == 0:
        return []
    cand = cKDTree(ref_xy).query_ball_tree(cKDTree(det_xy), radius)
    ii, jj = [], []
    for i, js in enumerate(cand):
        for j in js:
            if max_dh_pct is not None and ref_h is not None and det_h is not None:
                hr, hd = ref_h[i], det_h[j]
                if np.isfinite(hr) and np.isfinite(hd) and abs(hd - hr) > max_dh_pct / 100.0 * hr:
                    continue
            ii.append(i)
            jj.append(j)
    if not ii:
        return []
    ii, jj = np.array(ii), np.array(jj)
    dist = np.hypot(ref_xy[ii, 0] - det_xy[jj, 0], ref_xy[ii, 1] - det_xy[jj, 1])
    # Graf med referensträden som nod 0..n_r-1 och de detekterade som n_r..
    g = coo_matrix((np.ones(len(ii)), (ii, jj + n_r)), shape=(n_r + n_d, n_r + n_d))
    _n, labels = connected_components(g, directed=False)
    pairs = []
    big = radius * 1e3 + 1.0
    comp = labels[ii]
    for c in np.unique(comp):
        sel = comp == c
        ri, dj, dd = ii[sel], jj[sel], dist[sel]
        ur, ir = np.unique(ri, return_inverse=True)
        ud, jd = np.unique(dj, return_inverse=True)
        cost = np.full((len(ur), len(ud)), big)
        cost[ir, jd] = dd
        rows, cols = linear_sum_assignment(cost)
        for r, k in zip(rows, cols):
            if cost[r, k] < big:
                pairs.append((int(ur[r]), int(ud[k]), float(cost[r, k])))
    return pairs


def _wilson(k, n, z=1.96):
    """95 % konfidensintervall (Wilson) för andelen k/n, eller (None, None)."""
    if n <= 0:
        return None, None
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, mid - half), min(1.0, mid + half)


def _scores(n_ref, n_det_in, hits_ref, hits_det_in):
    """Träffsäkerhet, precision och F1. hits_det_in: parade detekterade inom provytan."""
    recall = hits_ref / n_ref if n_ref else None
    precision = hits_det_in / n_det_in if n_det_in else None
    if recall is None or precision is None:
        return recall, precision, None
    f1 = 2 * recall * precision / (recall + precision) if recall + precision > 0 else 0.0
    return recall, precision, f1


def _split_merge(ref_xy, det_xy, radius):
    """
    Antal referensträd med två eller fler detekterade träd inom radien (möjlig
    uppdelning av en krona) och detekterade träd med två eller fler
    referensträd inom radien (möjlig sammanslagning).
    """
    from scipy.spatial import cKDTree
    if len(ref_xy) == 0 or len(det_xy) == 0:
        return 0, 0
    a = cKDTree(ref_xy).query_ball_tree(cKDTree(det_xy), radius)
    b = cKDTree(det_xy).query_ball_tree(cKDTree(ref_xy), radius)
    return sum(1 for x in a if len(x) >= 2), sum(1 for x in b if len(x) >= 2)


# =============================================================================
# Toolbox
# =============================================================================

class Toolbox:
    def __init__(self):
        self.label = "Detektera träd"
        self.alias = "detektera_trad"
        self.tools = [DetekteraTrad, UtvarderaNoggrannhet]
        _write_tool_metadata(DetekteraTrad, self.alias, TOOLTIPS_DETECT, DETECT_SUMMARY)
        _write_tool_metadata(UtvarderaNoggrannhet, self.alias, TOOLTIPS_EVAL, EVAL_SUMMARY)


class DetekteraTrad:
    def __init__(self):
        self.label = "Detektera träd"
        self.description = DETECT_SUMMARY
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_mode = arcpy.Parameter(
            displayName="Avgränsa området med", name="aoi_mode", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_mode.filter.type = "ValueList"
        p_mode.filter.list = [AOI_POLYGONS, AOI_EXTENT]
        p_mode.value = AOI_POLYGONS

        # Båda är Optional i ramverket; updateMessages kräver den som valts.
        p_aoi = arcpy.Parameter(
            displayName="Intresseområde (polygoner)", name="aoi",
            datatype="GPFeatureRecordSetLayer", parameterType="Optional", direction="Input",
        )
        p_aoi.filter.list = ["Polygon"]
        p_aoi.value = _polygon_schema()

        p_extent = arcpy.Parameter(
            displayName="Utbredning", name="aoi_extent", datatype="GPExtent",
            parameterType="Optional", direction="Input",
        )
        p_extent.enabled = False

        p_source = arcpy.Parameter(
            displayName="Punktmoln", name="source", datatype="GPString",
            parameterType="Required", direction="Input",
        )
        p_source.filter.type = "ValueList"
        p_source.filter.list = [SRC_STAC, SRC_LOCAL]
        p_source.value = SRC_STAC

        p_files = arcpy.Parameter(
            displayName="LAS/LAZ-filer", name="las_files", datatype="DEFile",
            parameterType="Optional", direction="Input", multiValue=True,
        )
        p_files.filter.list = ["las", "laz"]
        p_files.enabled = False

        p_out = arcpy.Parameter(
            displayName="Utdata, trädtoppar", name="out_features", datatype="DEFeatureClass",
            parameterType="Required", direction="Output",
        )
        default_out = _unique_output("tradtoppar")
        if default_out:
            p_out.value = default_out
        _set_symbology(p_out, LYRX_TOPS)

        p_crowns = arcpy.Parameter(
            displayName="Utdata, trädkronor", name="out_crowns", datatype="DEFeatureClass",
            parameterType="Optional", direction="Output",
        )
        default_crowns = _unique_output("tradkronor")
        if default_crowns:
            p_crowns.value = default_crowns
        _set_symbology(p_crowns, LYRX_CROWNS)

        p_min = arcpy.Parameter(
            displayName="Lägsta trädhöjd (m)", name="min_height", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Detektering",
        )
        p_min.value = DEFAULT_MIN_HEIGHT

        p_sigma = arcpy.Parameter(
            displayName="Utjämning (m)", name="sigma", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Detektering",
        )
        p_sigma.value = DEFAULT_SIGMA

        p_a = arcpy.Parameter(
            displayName="Sökfönster, A (m per m trädhöjd)", name="win_a", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Detektering",
        )
        p_a.value = DEFAULT_WIN_A

        p_b = arcpy.Parameter(
            displayName="Sökfönster, B (m)", name="win_b", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Detektering",
        )
        p_b.value = DEFAULT_WIN_B

        p_filter = arcpy.Parameter(
            displayName="Ta bort toppar på byggnader och ledningar", name="filter_objects",
            datatype="GPBoolean", parameterType="Optional", direction="Input",
            category="Detektering",
        )
        p_filter.value = True

        p_key = arcpy.Parameter(
            displayName="Consumer key", name="consumer_key", datatype="GPString",
            parameterType="Optional", direction="Input", category="Inloggning",
        )
        p_key.value = _saved_login()[0] or None
        p_secret = arcpy.Parameter(
            displayName="Consumer secret", name="consumer_secret", datatype="GPStringHidden",
            parameterType="Optional", direction="Input", category="Inloggning",
        )

        p_save = arcpy.Parameter(
            displayName="Spara nyckeln på den här datorn", name="save_login", datatype="GPBoolean",
            parameterType="Optional", direction="Input", category="Inloggning",
        )
        p_save.value = False

        p_mem = arcpy.Parameter(
            displayName="Minnesbudget (GB)", name="memory_gb", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Avancerat",
        )
        p_workers = arcpy.Parameter(
            displayName="Block samtidigt", name="workers", datatype="GPLong",
            parameterType="Optional", direction="Input", category="Avancerat",
        )
        p_workers.filter.type = "Range"
        p_workers.filter.list = [1, 32]

        return [p_mode, p_aoi, p_extent, p_source, p_files, p_out, p_crowns, p_min, p_sigma,
                p_a, p_b, p_filter, p_key, p_secret, p_mem, p_workers, p_save]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        p = {q.name: q for q in parameters}
        by_extent = p["aoi_mode"].valueAsText == AOI_EXTENT
        p["aoi"].enabled = not by_extent
        p["aoi_extent"].enabled = by_extent
        local = p["source"].valueAsText == SRC_LOCAL
        p["las_files"].enabled = local
        p["consumer_key"].enabled = not local
        p["consumer_secret"].enabled = not local
        p["save_login"].enabled = not local

    def updateMessages(self, parameters):
        p = {q.name: q for q in parameters}

        if p["aoi_mode"].valueAsText == AOI_EXTENT:
            if not p["aoi_extent"].valueAsText:
                p["aoi_extent"].setErrorMessage("Ange en utbredning.")
        elif not p["aoi"].valueAsText:
            p["aoi"].setErrorMessage(AOI_EMPTY)
        elif not _is_layer(p["aoi"].value) and _has_features(p["aoi"].value) is False:
            # Tom polygonmall eller tom featureklass. Ett valt lager kontrolleras
            # först vid körningen.
            p["aoi"].setErrorMessage(AOI_EMPTY)

        if p["source"].valueAsText == SRC_LOCAL:
            if not p["las_files"].valueAsText:
                p["las_files"].setErrorMessage("Välj minst en LAS- eller LAZ-fil.")
        else:
            _check_login_fields(p)

        mh = p["min_height"].value
        if mh is not None and not (0 <= mh <= 50):
            p["min_height"].setErrorMessage("Ange 0 till 50 m.")
        elif mh is not None and mh < 2:
            p["min_height"].setWarningMessage(
                "Under 2 m blir buskar, block och små fel i markmodellen ofta toppar.")
        sg = p["sigma"].value
        if sg is not None and not (0 <= sg <= 3):
            p["sigma"].setErrorMessage("Ange 0 till 3 m.")
        a = p["win_a"].value
        if a is not None and not (0 <= a <= 1):
            p["win_a"].setErrorMessage("Ange 0 till 1.")
        b = p["win_b"].value
        if b is not None and not (0 <= b <= 12):
            p["win_b"].setErrorMessage("Ange 0 till 12 m.")

        mem = p["memory_gb"].value
        if mem is not None and mem < MIN_MEMORY_GB:
            p["memory_gb"].setErrorMessage("Ange minst 0,5 GB, eller lämna tomt för automatiskt.")
        elif mem is not None:
            avail = _available_memory_gb()
            if avail is not None and mem > avail:
                p["memory_gb"].setWarningMessage(
                    "Mer än det lediga arbetsminnet ({} GB). Windows börjar då skriva till disk "
                    "och körningen blir mycket långsam.".format(_fmt_m(round(avail, 1))))

    def execute(self, parameters, messages):
        p = {q.name: q for q in parameters}
        try:
            if p["aoi_mode"].valueAsText == AOI_EXTENT:
                if not p["aoi_extent"].valueAsText:
                    raise ValueError("Ange en utbredning.")
                aoi, note = _extent_polygon(p["aoi_extent"].value, p["aoi_extent"].valueAsText)
                messages.addMessage(note)
            else:
                if not p["aoi"].valueAsText:
                    raise ValueError(AOI_EMPTY)
                aoi = _aoi_geometry(p["aoi"].value, messages)
            local = p["source"].valueAsText == SRC_LOCAL
            files = None
            if local:
                files = [v for v in (p["las_files"].valueAsText or "").split(";") if v.strip()]
                files = [f.strip().strip("'\"") for f in files]
                if not files:
                    raise ValueError("Välj minst en LAS- eller LAZ-fil.")
            params = {
                "hmin": DEFAULT_MIN_HEIGHT if p["min_height"].value is None
                else float(p["min_height"].value),
                "sigma": DEFAULT_SIGMA if p["sigma"].value is None else float(p["sigma"].value),
                "a": DEFAULT_WIN_A if p["win_a"].value is None else float(p["win_a"].value),
                "b": DEFAULT_WIN_B if p["win_b"].value is None else float(p["win_b"].value),
                "filter": p["filter_objects"].value is not False,
            }
            key, secret = "", ""
            if p["source"].valueAsText != SRC_LOCAL:
                key, secret = _login(p)
                if p["save_login"].value:
                    _remember_login(key, secret, messages)
            _run_detect(
                aoi, files,
                key,
                secret,
                p["out_features"].valueAsText, params,
                p["memory_gb"].value or None, p["workers"].value or None, messages,
                out_crowns=p["out_crowns"].valueAsText or None,
            )
        except ValueError as exc:
            error = str(exc)
        else:
            error = None
        # Utanför except-blocket, så att Pro inte skriver ut hela kedjan av
        # undantag under det läsbara felmeddelandet.
        if error:
            messages.addErrorMessage(error)
            raise arcpy.ExecuteError

    def postExecute(self, parameters):
        return


class UtvarderaNoggrannhet:
    def __init__(self):
        self.label = "Utvärdera noggrannhet"
        self.description = EVAL_SUMMARY
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_det = arcpy.Parameter(
            displayName="Detekterade träd", name="detected", datatype="GPFeatureLayer",
            parameterType="Required", direction="Input",
        )
        p_det.filter.list = ["Point"]
        p_ref = arcpy.Parameter(
            displayName="Referensträd", name="reference", datatype="GPFeatureLayer",
            parameterType="Required", direction="Input",
        )
        p_ref.filter.list = ["Point"]
        p_plots = arcpy.Parameter(
            displayName="Provytor", name="plots", datatype="GPFeatureLayer",
            parameterType="Required", direction="Input",
        )
        p_plots.filter.list = ["Polygon"]

        p_out = arcpy.Parameter(
            displayName="Utdata, parade träd", name="out_features", datatype="DEFeatureClass",
            parameterType="Required", direction="Output",
        )
        default_out = _unique_output("utvardering_trad")
        if default_out:
            p_out.value = default_out
        p_table = arcpy.Parameter(
            displayName="Utdata, mått per provyta", name="out_table", datatype="DETable",
            parameterType="Optional", direction="Output",
        )

        p_radius = arcpy.Parameter(
            displayName="Sökradie (m)", name="search_radius", datatype="GPDouble",
            parameterType="Required", direction="Input",
        )
        p_radius.value = DEFAULT_SEARCH_RADIUS

        p_dh = arcpy.Parameter(
            displayName="Höjdfält, detekterade träd", name="det_height_field", datatype="Field",
            parameterType="Optional", direction="Input", category="Höjd",
        )
        p_dh.parameterDependencies = [p_det.name]
        p_dh.filter.list = ["Short", "Long", "Float", "Double"]
        p_rh = arcpy.Parameter(
            displayName="Höjdfält, referensträd", name="ref_height_field", datatype="Field",
            parameterType="Optional", direction="Input", category="Höjd",
        )
        p_rh.parameterDependencies = [p_ref.name]
        p_rh.filter.list = ["Short", "Long", "Float", "Double"]
        p_maxdh = arcpy.Parameter(
            displayName="Största höjdskillnad (%)", name="max_height_diff", datatype="GPDouble",
            parameterType="Optional", direction="Input", category="Höjd",
        )
        p_maxdh.value = DEFAULT_MAX_DH_PCT

        p_pid = arcpy.Parameter(
            displayName="Provytans ID-fält", name="plot_id_field", datatype="Field",
            parameterType="Optional", direction="Input", category="Provytor",
        )
        p_pid.parameterDependencies = [p_plots.name]

        return [p_det, p_ref, p_plots, p_out, p_table, p_radius, p_dh, p_rh, p_maxdh, p_pid]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        p = {q.name: q for q in parameters}
        # Föreslå TRADHOJD (utdata från Detektera träd) om fältet finns och inget är valt.
        if p["detected"].altered and not p["det_height_field"].valueAsText and \
                not p["det_height_field"].altered:
            try:
                names = [f.name.upper() for f in arcpy.ListFields(p["detected"].value)]
                if "TRADHOJD" in names:
                    p["det_height_field"].value = "TRADHOJD"
            except Exception:
                pass

    def updateMessages(self, parameters):
        p = {q.name: q for q in parameters}
        r = p["search_radius"].value
        if r is not None and not (0.1 <= r <= 20):
            p["search_radius"].setErrorMessage("Ange 0,1 till 20 m.")
        dh = p["max_height_diff"].value
        if dh is not None and dh <= 0:
            p["max_height_diff"].setErrorMessage("Ange ett positivt tal, eller lämna tomt.")
        if dh is not None and not (p["det_height_field"].valueAsText
                                   and p["ref_height_field"].valueAsText):
            p["max_height_diff"].setWarningMessage(
                "Används bara när både de detekterade trädens och referensträdens höjdfält är "
                "angivna.")

    def execute(self, parameters, messages):
        p = {q.name: q for q in parameters}
        try:
            _run_evaluate(
                p["detected"].value, p["reference"].value, p["plots"].value,
                p["out_features"].valueAsText, p["out_table"].valueAsText or None,
                float(p["search_radius"].value),
                p["det_height_field"].valueAsText or None,
                p["ref_height_field"].valueAsText or None,
                None if p["max_height_diff"].value is None else float(p["max_height_diff"].value),
                p["plot_id_field"].valueAsText or None,
                messages,
            )
        except ValueError as exc:
            error = str(exc)
        else:
            error = None
        if error:
            messages.addErrorMessage(error)
            raise arcpy.ExecuteError

    def postExecute(self, parameters):
        return


# =============================================================================
# Detektera träd (separat funktion - går att testa utanför Pro)
# =============================================================================

TOP_FIELDS = [
    ["TRAD_ID", "LONG", "Träd-ID"],
    ["TRADHOJD", "DOUBLE", "Trädhöjd (m)"],
    ["MARKHOJD", "DOUBLE", "Markhöjd (m ö.h., RH 2000)"],
    ["KALLA", "TEXT", "Källa", 120],
    ["SKANNING", "TEXT", "Skanningsområde", 16],
    ["INSAMLAD", "TEXT", "Insamlad", 32],
]


CROWN_FIELDS = [
    ["TRAD_ID", "LONG", "Träd-ID"],
    ["TRADHOJD", "DOUBLE", "Trädhöjd (m)"],
    ["KRONAREA", "DOUBLE", "Kronans yta (m²)"],
    ["KRONDIAM", "DOUBLE", "Kronans diameter (m, cirkel med samma yta)"],
]


class _TopSink:
    """
    Tar emot topparna (och kronorna) block för block och skriver dem till
    featureklasser i memory: bara toppar inom området, och utan byggnader och
    ledningar om filtret är på. TRAD_ID numreras i den ordning de skrivs och
    är detsamma i punkterna och kronorna.
    """

    POINT_NAMES = ["SHAPE@XY", "TRAD_ID", "TRADHOJD", "MARKHOJD", "KALLA", "SKANNING",
                   "INSAMLAD"]
    CROWN_NAMES = ["SHAPE@", "TRAD_ID", "TRADHOJD", "KRONAREA", "KRONDIAM"]

    def __init__(self, tiles, aoi, params, temps):
        self.tiles, self.aoi, self.params = tiles, aoi, params
        self.tile_text = [(t["id"][:120], (t["area"] or "")[:16], _capture_period(t)[:32])
                          for t in tiles]
        self.n = self.n_hard = self.n_thin = 0
        self.heights, self.crown_areas = [], []
        self.points = _memory_name("toppar")
        temps.append(self.points)
        arcpy.management.CreateFeatureclass(os.path.dirname(self.points),
                                            os.path.basename(self.points), "POINT",
                                            spatial_reference=_tm())
        arcpy.management.AddFields(self.points, TOP_FIELDS)
        self._pcur = arcpy.da.InsertCursor(self.points, self.POINT_NAMES)
        self.crowns, self._ccur = None, None
        if params.get("crowns"):
            self.crowns = _memory_name("kronor")
            temps.append(self.crowns)
            arcpy.management.CreateFeatureclass(os.path.dirname(self.crowns),
                                                os.path.basename(self.crowns), "POLYGON",
                                                spatial_reference=_tm())
            arcpy.management.AddFields(self.crowns, CROWN_FIELDS)
            self._ccur = arcpy.da.InsertCursor(self.crowns, self.CROWN_NAMES)

    def add(self, z):
        x, y, h, g = z["x"], z["y"], z["h"], z["ground"]
        inside = _in_polygon(self.aoi, x, y)
        if self.params["filter"]:
            hard, thin = z["hard"], z["thin"]
            self.n_hard += int((hard & inside).sum())
            self.n_thin += int((thin & inside).sum())
            keep = inside & ~(hard | thin)
        else:
            keep = inside
        wkb = offs = None
        if self._ccur is not None and "crown_offs" in z:
            wkb, offs = z["crown_wkb"].tobytes(), z["crown_offs"]
        tidx = _tile_index(self.tiles, x, y)
        sr = _tm()
        for k in np.flatnonzero(keep):
            self.n += 1
            hk = round(float(h[k]), 2)
            kal, skn, ins = self.tile_text[tidx[k]] if tidx[k] >= 0 else ("", "", "")
            self._pcur.insertRow(((float(x[k]), float(y[k])), self.n, hk,
                                  None if not np.isfinite(g[k]) else round(float(g[k]), 2),
                                  kal, skn, ins))
            self.heights.append(hk)
            if wkb is not None and offs[k + 1] > offs[k]:
                geom = arcpy.FromWKB(bytearray(wkb[offs[k]:offs[k + 1]]), sr)
                area = geom.area
                if area > 0:
                    self._ccur.insertRow((geom, self.n, hk, round(area, 2),
                                          round(2 * math.sqrt(area / math.pi), 2)))
                    self.crown_areas.append(area)

    def close(self):
        # Cursorerna måste släppas innan featureklasserna exporteras eller tas bort.
        for name in ("_pcur", "_ccur"):
            cur = getattr(self, name)
            if cur is not None:
                del cur
                setattr(self, name, None)


def _in_polygon(poly, x, y):
    """
    Sant för punkterna inom polygonen (jämn-udda-regeln över alla ringar, så
    hål och flera delar hanteras). Ringarna läses med __geo_interface__, som
    ger varje ring för sig (se CLAUDE.md).
    """
    gi = poly.__geo_interface__
    parts = gi["coordinates"] if gi["type"] == "MultiPolygon" else [gi["coordinates"]]
    inside = np.zeros(len(x), bool)
    for part in parts:
        for ring in part:
            r = np.asarray(ring, float)
            x0, y0, x1, y1 = r[:-1, 0], r[:-1, 1], r[1:, 0], r[1:, 1]
            for a, b, c, d in zip(x0, y0, x1, y1):
                cross = (b > y) != (d > y)
                with np.errstate(divide="ignore", invalid="ignore"):
                    xi = a + (y - b) * (c - a) / (d - b)
                inside ^= cross & (x < xi)
    return inside


def _tile_index(tiles, x, y):
    """Index för den första rutan som innehåller varje punkt, -1 om ingen."""
    idx = np.full(len(x), -1, np.int64)
    for k, t in enumerate(tiles):
        x0, y0, x1, y1 = t["bbox"]
        sel = (idx < 0) & (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        idx[sel] = k
    return idx


def _run_detect(aoi, files, key, secret, out_fc, params, memory_gb, workers, messages,
                out_crowns=None):
    """
    aoi: polygon i SWEREF 99 TM. files: lokala LAS/LAZ-filer, eller None för
    Lantmäteriets STAC. params: {"hmin", "sigma", "a", "b", "filter"}.
    out_crowns: featureklass för kronorna, eller None för inga kronor.
    Returnerar sökvägen till punktfeatureklassen.
    """
    if not out_fc:
        raise ValueError("Ange utdata.")
    params = dict(params, crowns=bool(out_crowns))
    if out_crowns and os.path.normcase(out_crowns) == os.path.normcase(out_fc):
        raise ValueError("Trädtopparna och kronorna kan inte ha samma utdatanamn.")
    if files is None and not (key and secret):
        raise ValueError("Ange consumer key och secret under Inloggning, eller välj lokala "
                         "LAS/LAZ-filer.")
    if files:
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            raise ValueError("Filen finns inte: {}".format(", ".join(missing)))
    for path in (out_fc, out_crowns):
        out_dir = os.path.dirname(path) if path else ""
        if out_dir and not arcpy.Exists(out_dir):
            raise ValueError("Utdata-arbetsytan finns inte: {}".format(out_dir))
    if importlib.util.find_spec("pdal") is None or importlib.util.find_spec("scipy") is None:
        raise ValueError(
            "Python-paketen pdal och scipy behövs och saknas i den aktiva miljön ({}). De ingår "
            "i standardmiljön arcgispro-py3.".format(sys.prefix))
    memory_gb, workers, auto_note = _resolve_resources(memory_gb, workers)
    if auto_note:
        messages.addMessage(auto_note)

    ext = aoi.extent
    grid = _grid(ext, BLOCK_GRID_M)
    area_km2 = (ext.XMax - ext.XMin) * (ext.YMax - ext.YMin) / 1e6
    messages.addMessage("Intresseområdets utbredning: {} km² (SWEREF 99 TM).".format(
        _fmt_m(round(area_km2, 2))))
    messages.addMessage(
        "Detektering: lägsta höjd {} m, utjämning {} m, sökfönstrets diameter {} x höjd + {} m."
        .format(_fmt_m(params["hmin"]), _fmt_m(params["sigma"]), _fmt_m(params["a"]),
                _fmt_m(params["b"])))
    if area_km2 > LARGE_AREA_KM2:
        messages.addWarningMessage(
            "Stort område ({:.0f} km²). Alla punkter inom utbredningen läses, ungefär 3 "
            "miljoner per km², så körningen kan ta lång tid.".format(area_km2))

    steps = _Steps(5, messages)
    workdir = os.path.join(_local_scratch_folder(messages), "trad_" + uuid.uuid4().hex[:8])
    os.makedirs(workdir)
    pool = None
    sink = None
    temps = []
    try:
        # ── 1. Rutor ─────────────────────────────────────────────────────────
        if files is None:
            steps.next("söker rutor i Lantmäteriets STAC-katalog")
            wgs = aoi.projectAs(arcpy.SpatialReference(4326)).extent
            tiles = _pick_tiles(_stac_search((wgs.XMin, wgs.YMin, wgs.XMax, wgs.YMax)), aoi)
            if not tiles:
                raise ValueError(
                    "Inga laserdata för området. Laserdata Skog täcker ungefär 75 % av Sverige, "
                    "men inte fjällen.")
        else:
            steps.next("läser huvudet i {} LAS/LAZ-fil(er)".format(len(files)))
            pool = _WorkerPool(1, workdir)
            try:
                headers = pool.call({"op": "header", "paths": files})["headers"]
            finally:
                pool.close()
                pool = None
            tiles = _local_tiles(headers, aoi, messages)
            if not tiles:
                raise ValueError("Ingen av filerna överlappar området.")
        blocks, side = _plan_blocks(aoi, grid, tiles, memory_gb, workers)
        expected = sum(b["expected"] for b in blocks)
        steps.done("{} ruta/rutor. Ungefär {} miljoner punkter väntas.".format(
            len(tiles), _fmt_m(round(expected / 1e6, 1))))
        for t in tiles:
            d = _actual_density(t)
            if files is None:
                messages.addMessage(
                    "    {}: skanningsområde {}, insamlad {} ({}), {} punkter per m² med alla "
                    "returer.".format(t["id"], t["area"] or "okänt", _capture_period(t),
                                      _season(_months([t])),
                                      _fmt_m(round(d, 1)) if d else "?"))
            else:
                messages.addMessage("    {}: {} punkter per m² med alla returer.".format(
                    t["id"], _fmt_m(round(d, 1)) if d else "?"))
        if files is None and len({_capture_period(t) for t in tiles}) > 1:
            messages.addWarningMessage(
                "Rutorna är skannade vid olika tillfällen. Antalet toppar per ha kan skilja sig "
                "vid rutgränsen, särskilt om årstid eller år skiljer sig.")
        workers = max(1, min(workers, len(blocks)))
        messages.addMessage(
            "    {} block på upp till {:.0f} x {:.0f} m, {} åt gången, inom en minnesbudget "
            "på {} GB.".format(len(blocks), side, side, workers, _fmt_m(round(memory_gb, 1))))

        # ── 2. Behörighet ────────────────────────────────────────────────────
        token = None
        steps.next("hämtar token och kontrollerar behörighet" if files is None
                   else "kontrollerar filerna")
        if files is None:
            token = _Token(key, secret)
            _check_access(tiles[0]["href"], token.get())
        steps.done()

        # ── 3. Blocken ───────────────────────────────────────────────────────
        steps.next("läser punkter och söker trädtoppar, {} block, {} åt gången".format(
            len(blocks), workers))
        arcpy.SetProgressor("step", "", 0, 100, 1)
        # Toppar och kronor skrivs block för block till featureklasser i
        # memory, så att de aldrig ligger som Python-objekt för hela området.
        sink = _TopSink(tiles, aoi, params, temps)
        n_points = n_ground = n_first = n_nohgt = 0
        peak = 0
        done_expected = 0.0
        t_blocks = time.time()
        pool = _WorkerPool(workers, workdir)
        failed = True
        try:
            pending = list(blocks)
            running = {}
            for i in range(workers):
                b = pending.pop(0)
                running[i] = b
                pool.submit(i, _block_task(b, token.get() if token else "", params, workdir))
            k = 0
            while running:
                got = pool.get(timeout=1.0)
                if getattr(arcpy.env, "isCancelled", False):
                    raise ValueError("Avbrutet av användaren.")
                if got is None:
                    continue
                i, reply = got
                b = running.pop(i)
                if reply is None:
                    raise RuntimeError("Arbetsprocessen avslutades oväntat under block rad {}, "
                                       "kolumn {}:\n{}".format(b["row"], b["col"], pool.log_tail(i)))
                if not reply["ok"]:
                    raise RuntimeError("Fel i block rad {}, kolumn {}:\n{}".format(
                        b["row"], b["col"], reply["error"]))
                if pending:
                    nb = pending.pop(0)
                    running[i] = nb
                    pool.submit(i, _block_task(nb, token.get() if token else "", params, workdir))

                res = reply["result"]
                k += 1
                n_points += res["points"]
                n_ground += res["ground"]
                n_first += res["first"]
                n_nohgt += res["no_height"]
                peak = max(peak, res.get("peak_bytes") or 0)
                if res["npz"]:
                    with np.load(res["npz"]) as z:
                        sink.add(z)
                    os.remove(res["npz"])

                done_expected += b["expected"]
                elapsed = time.time() - t_blocks
                left = elapsed / done_expected * (expected - done_expected) if done_expected else 0
                messages.addMessage(
                    "    Block {} av {} (rad {}, kolumn {}): {} punkter, {} toppar.".format(
                        k, len(blocks), b["row"], b["col"], _fmt_count(res["points"]),
                        _fmt_count(res["tops"])))
                steps.label("block {} av {} klara{}".format(
                    k, len(blocks), ", ca {} kvar".format(_fmt_duration(left)) if running else ""))
                arcpy.SetProgressorPosition(min(100, int(100 * done_expected / max(expected, 1))))
            failed = False
        finally:
            pool.close(kill=failed)
            pool = None

        if n_points == 0:
            raise ValueError("Inga punkter inom området.")
        note = "{} punkter, varav {} markpunkter och {} förstareturer.".format(
            _fmt_count(n_points), _fmt_count(n_ground), _fmt_count(n_first))
        if n_nohgt:
            note += " {} förstareturer saknade markpunkter runt omkring och användes inte.".format(
                _fmt_count(n_nohgt))
        if peak:
            note += " Största minnesåtgång i en process: {} MB.".format(int(peak / 1e6))
        steps.done(note)

        # ── 4. Utdata ────────────────────────────────────────────────────────
        steps.next("skriver toppar{} till utdata".format(" och kronor" if out_crowns else ""))
        sink.close()
        n_in = sink.n
        _export(sink.points, out_fc, _out_sr(), messages)
        if out_crowns:
            _export(sink.crowns, out_crowns, _out_sr(), messages)
        area_ha = aoi.area / 1e4
        heights = np.array(sink.heights)
        note = "{} toppar inom området, {} per ha{}.".format(
            _fmt_count(n_in), _fmt_m(round(n_in / area_ha, 1)) if area_ha > 0 else "?",
            ", trädhöjd median {} m (5-95 %: {}-{} m)".format(
                _fmt_m(round(float(np.median(heights)), 1)),
                _fmt_m(round(float(np.percentile(heights, 5)), 1)),
                _fmt_m(round(float(np.percentile(heights, 95)), 1))) if heights.size else "")
        note += (" Borttagna: {} på byggnader och andra hårda ytor, {} på ledningar och andra "
                 "tunna föremål.".format(_fmt_count(sink.n_hard), _fmt_count(sink.n_thin))
                 if params["filter"] else " Toppar på byggnader och ledningar togs inte bort.")
        if out_crowns:
            areas = np.array(sink.crown_areas)
            note += " {} kronor, yta median {} m² (5-95 %: {}-{} m²).".format(
                _fmt_count(areas.size),
                *[_fmt_m(round(float(v), 1)) for v in (np.median(areas),
                                                        *np.percentile(areas, [5, 95]))]
            ) if areas.size else " Inga kronor."
        steps.done(note)

        # ── 5. Metadata ──────────────────────────────────────────────────────
        steps.next("skriver metadata")
        for path, kind in ((out_fc, "toppar"), (out_crowns, "kronor")):
            if not path:
                continue
            try:
                _write_detect_metadata(path, kind, tiles, params, files is not None, n_in,
                                       area_ha)
            except Exception as exc:
                messages.addWarningMessage("Kunde inte skriva metadata för {}: {}".format(
                    path, exc))
        steps.done()
        messages.addMessage("Klart på {}.".format(_fmt_duration(time.time() - steps.t0)))
        return out_fc
    finally:
        if pool is not None:
            pool.close(kill=True)
        if sink is not None:
            sink.close()
        arcpy.ResetProgressor()
        _delete_quietly(*temps)
        _remove_workdir(workdir)


def _write_detect_metadata(path, kind, tiles, params, local, n, area_ha):
    """kind: "toppar" eller "kronor"."""
    h = escape
    settings = [
        ("Metod", "Lokala maxima i en ythöjdmodell över mark (CHM) med sökfönster som växer med "
                  "höjden (Popescu & Wynne 2004, som lmf() i R-paketet lidR)"),
        ("CHM", "Förstareturernas höjd över mark, triangulerad och rastrerad till {} m".format(
            _fmt_m(_tw.CHM_RES))),
        ("Höjd över mark", "Markpunkterna (klass 2) triangulerade och rastrerade till {} m, "
                           "bilinjärt interpolerade".format(_fmt_m(_tw.DTM_RES))),
        ("Utjämning", "Gaussiskt filter, standardavvikelse {} m".format(_fmt_m(params["sigma"]))),
        ("Sökfönster", "Diameter {} x höjd + {} m, radie 1-6 m i halvmetersteg".format(
            _fmt_m(params["a"]), _fmt_m(params["b"]))),
        ("Lägsta trädhöjd", "{} m".format(_fmt_m(params["hmin"]))),
        ("Byggnader och ledningar", (
            "Borttagna: toppar på plana ytor med enkelreturer (minst {} m höga, högst {} m "
            "ojämnhet, minst {} % enkelreturer, mask utvidgad {} m) och toppar med högst {} "
            "returer över 2 m som är under {} % av returerna i sökfönstret".format(
                _fmt_m(_tw.HARD_MIN_H), _fmt_m(_tw.HARD_MAX_ROUGH),
                int(_tw.HARD_MIN_SINGLE * 100), _fmt_m(_tw.HARD_GROW_M),
                _tw.THIN_MAX_POINTS, int(_tw.THIN_MAX_COVER * 100))
            if params["filter"] else "Inte borttagna")),
        ("Toppens läge", "Högsta förstareturen inom {} m från toppen i den utjämnade "
                         "modellen".format(_fmt_m(_tw.APEX_RADIUS))),
        ("Antal toppar", "{} ({} per ha)".format(_fmt_count(n), _fmt_m(round(n / area_ha, 1))
                                                 if area_ha > 0 else "?")),
        ("Kronor", (
            "Fyllning uppifrån i den utjämnade CHM:en från en markör per topp (vattendelare), "
            "{} m i taget. En cell hör till kronan om den är minst {} % av toppens höjd, högst "
            "{} x trädhöjden från toppen, minst {} m över mark och inte hård yta. Kanterna "
            "förenklade med {} m tolerans. TRAD_ID kopplar kronan till sin topp".format(
                _fmt_m(_tw.CROWN_LEVEL_M), int(_tw.CROWN_MIN_REL * 100),
                _fmt_m(_tw.CROWN_MAX_R_FACTOR), _fmt_m(_tw.CROWN_MIN_H),
                _fmt_m(_tw.CROWN_SIMPLIFY_M))
            if params.get("crowns") else "Inte avgränsade")),
        ("Skapad", "{} med verktyget Detektera träd {} ({})".format(
            datetime.date.today().isoformat(), TOOL_VERSION, REPO_NAME)),
    ]
    rows = "".join("<tr><td>{}</td><td>{}</td></tr>".format(h(k), h(str(v))) for k, v in settings)
    tile_rows = "".join("<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(
        h(t["id"]), h(t["area"]), h(_capture_period(t))) for t in tiles)
    src = ("lokala LAS/LAZ-filer" if local
           else "Lantmäteriet, Laserdata Nedladdning, skog (STAC-samling {})".format(COLLECTION))
    desc = (
        "<p>Trädtoppar (krontoppar) hittade i laserdata. Läget är krontoppens, inte en uppmätt "
        "stambas. Vid 2-4 punkter per m² hittas främst dominerande träd. Standardvärdena är inte "
        "kalibrerade mot svenska referensdata.</p>"
        "<p><b>Inställningar</b></p><table border='1' cellpadding='3'>{rows}</table>"
        "<p><b>Källa:</b> {src}.</p>"
        "<table border='1' cellpadding='3'><tr><th>Ruta/fil</th><th>Skanningsområde</th>"
        "<th>Insamlad</th></tr>{tile_rows}</table>"
    ).format(rows=rows, src=h(src), tile_rows=tile_rows)
    md = arcpy.metadata.Metadata(path)
    if kind == "kronor":
        md.title = "Trädkronor från laserdata"
        md.summary = ("Trädkronor avgränsade i laserdata kring trädtopparna, med kronans yta "
                      "och diameter ({}).".format(src))
        md.tags = "trädkronor, kronavgränsning, laserskanning, Laserdata Skog, CHM"
    else:
        md.title = "Trädtoppar från laserdata"
        md.summary = "Trädtoppar (krontoppar) med trädhöjd, hittade i laserdata ({}).".format(src)
        md.tags = "trädtoppar, träddetektering, laserskanning, Laserdata Skog, CHM"
    md.description = desc
    if not local:
        md.credits = "© Lantmäteriet, Laserdata Nedladdning, skog."
        md.accessConstraints = (
            "Användningsvillkor för Laserdata Nedladdning, skog: {}".format(TERMS_URL))
    md.save()


# =============================================================================
# Utvärdera noggrannhet (separat funktion - går att testa utanför Pro)
# =============================================================================

def _read_points(layer, height_field, where_geom):
    """
    (oid, x, y, h) för punkterna i lagret som skär where_geom, i SWEREF 99 TM.
    Ett lagers urval och definitionsfråga följer med automatiskt.
    """
    fields = ["OID@", "SHAPE@XY"] + ([height_field] if height_field else [])
    out = []
    with arcpy.da.SearchCursor(layer, fields, spatial_reference=_tm(),
                               spatial_filter=where_geom,
                               spatial_relationship="INTERSECTS") as cur:
        for row in cur:
            if row[1] is None or row[1][0] is None:
                continue
            hv = row[2] if height_field else None
            out.append((row[0], row[1][0], row[1][1],
                        float(hv) if hv is not None else float("nan")))
    if not out:
        return np.zeros(0, np.int64), np.zeros((0, 2)), np.zeros(0)
    a = np.array(out, dtype=float)
    return a[:, 0].astype(np.int64), a[:, 1:3], a[:, 3]


def _check_sr(layer, label):
    sr = arcpy.Describe(layer).spatialReference
    if sr is None or not (sr.factoryCode or sr.exportToString()) or sr.name == "Unknown":
        raise ValueError("{} saknar koordinatsystem.".format(label))


def _run_evaluate(detected, reference, plots, out_fc, out_table, radius, det_hf, ref_hf,
                  max_dh_pct, plot_id_field, messages):
    """
    Para referensträd och detekterade träd per provyta och redovisa måtten.
    Returnerar en dict med måtten för alla provytor.
    """
    if not out_fc:
        raise ValueError("Ange utdata.")
    if radius <= 0:
        raise ValueError("Sökradien måste vara större än 0.")
    for lyr, label in ((detected, "Detekterade träd"), (reference, "Referensträd"),
                       (plots, "Provytor")):
        _check_sr(lyr, label)
    use_dh = bool(det_hf and ref_hf and max_dh_pct)
    steps = _Steps(3, messages)
    temps = []
    try:
        # ── 1. Provytor och träd ─────────────────────────────────────────────
        steps.next("läser provytor och träd")
        plot_list = []
        fields = ["OID@", "SHAPE@"] + ([plot_id_field] if plot_id_field else [])
        with arcpy.da.SearchCursor(plots, fields, spatial_reference=_tm()) as cur:
            for row in cur:
                if row[1] is None or row[1].area <= 0:
                    continue
                pid = str(row[2]) if plot_id_field and row[2] is not None else str(row[0])
                plot_list.append((pid, row[1]))
        if not plot_list:
            raise ValueError("Provytorna innehåller inga polygoner med yta.")
        data = []
        for pid, geom in plot_list:
            r_oid, r_xy, r_h = _read_points(reference, ref_hf, geom)
            near = geom.buffer(radius)
            d_oid, d_xy, d_h = _read_points(detected, det_hf, near)
            d_in = np.array([geom.contains(arcpy.PointGeometry(arcpy.Point(x, y), _tm()))
                             for x, y in d_xy], bool) if len(d_xy) else np.zeros(0, bool)
            data.append({"pid": pid, "geom": geom, "area_ha": geom.area / 1e4,
                         "r_oid": r_oid, "r_xy": r_xy, "r_h": r_h,
                         "d_oid": d_oid, "d_xy": d_xy, "d_h": d_h, "d_in": d_in})
        n_ref = sum(len(d["r_oid"]) for d in data)
        n_det = sum(int(d["d_in"].sum()) for d in data)
        steps.done("{} provytor, {} ha, {} referensträd och {} detekterade träd inom "
                   "provytorna.".format(len(data), _fmt_m(round(sum(d["area_ha"] for d in data), 2)),
                                        _fmt_count(n_ref), _fmt_count(n_det)))
        if n_ref == 0:
            raise ValueError("Inga referensträd ligger inom provytorna.")
        if len(plot_list) > 1:
            close = [(a[0], b[0]) for i, a in enumerate(plot_list) for b in plot_list[i + 1:]
                     if a[1].distanceTo(b[1]) < radius]
            if close:
                messages.addWarningMessage(
                    "Provytorna {} ligger närmare varandra än sökradien. Ett detekterat träd "
                    "mellan dem kan paras i båda.".format(
                        ", ".join("{} och {}".format(*c) for c in close[:5])))

        # ── 2. Parning ───────────────────────────────────────────────────────
        steps.next("parar ihop träd inom {} m".format(_fmt_m(radius)))

        def evaluate(rad):
            rows = []
            for d in data:
                pairs = _match(d["r_xy"], d["d_xy"], rad, d["r_h"], d["d_h"],
                               max_dh_pct if use_dh else None)
                rows.append(pairs)
            return rows

        all_pairs = evaluate(radius)
        per_plot = []
        tot = {"n_ref": 0, "n_det": 0, "hit_ref": 0, "hit_det": 0, "area": 0.0,
               "dx": [], "dy": [], "dh": [], "split": 0, "merge": 0}
        for d, pairs in zip(data, all_pairs):
            hit_ref = len(pairs)
            hit_det = sum(1 for _i, j, _ in pairs if d["d_in"][j])
            n_r, n_d = len(d["r_oid"]), int(d["d_in"].sum())
            dx = [d["d_xy"][j, 0] - d["r_xy"][i, 0] for i, j, _ in pairs]
            dy = [d["d_xy"][j, 1] - d["r_xy"][i, 1] for i, j, _ in pairs]
            dh = [d["d_h"][j] - d["r_h"][i] for i, j, _ in pairs
                  if np.isfinite(d["d_h"][j]) and np.isfinite(d["r_h"][i])]
            sp, mg = _split_merge(d["r_xy"], d["d_xy"][d["d_in"]], radius)
            rec, prec, f1 = _scores(n_r, n_d, hit_ref, hit_det)
            per_plot.append({"pid": d["pid"], "area_ha": d["area_ha"], "n_ref": n_r, "n_det": n_d,
                             "hits": hit_ref, "recall": rec, "precision": prec, "f1": f1,
                             "dx": dx, "dy": dy, "dh": dh})
            tot["n_ref"] += n_r
            tot["n_det"] += n_d
            tot["hit_ref"] += hit_ref
            tot["hit_det"] += hit_det
            tot["area"] += d["area_ha"]
            tot["dx"] += dx
            tot["dy"] += dy
            tot["dh"] += dh
            tot["split"] += sp
            tot["merge"] += mg
        rec, prec, f1 = _scores(tot["n_ref"], tot["n_det"], tot["hit_ref"], tot["hit_det"])
        rl, rh = _wilson(tot["hit_ref"], tot["n_ref"])
        pl, ph = _wilson(tot["hit_det"], tot["n_det"])
        dx, dy, dh = np.array(tot["dx"]), np.array(tot["dy"]), np.array(tot["dh"])
        summary = {
            "n_ref": tot["n_ref"], "n_det": tot["n_det"], "hits": tot["hit_ref"],
            "recall": rec, "precision": prec, "f1": f1,
            "recall_ci": (rl, rh), "precision_ci": (pl, ph),
            "ref_per_ha": tot["n_ref"] / tot["area"] if tot["area"] else None,
            "det_per_ha": tot["n_det"] / tot["area"] if tot["area"] else None,
            "bias_x": float(dx.mean()) if dx.size else None,
            "bias_y": float(dy.mean()) if dy.size else None,
            "rmse_xy": float(np.sqrt(np.mean(dx * dx + dy * dy))) if dx.size else None,
            "bias_h": float(dh.mean()) if dh.size else None,
            "rmse_h": float(np.sqrt(np.mean(dh * dh))) if dh.size else None,
            "split": tot["split"], "merge": tot["merge"], "per_plot": per_plot,
        }
        steps.done()

        m = messages.addMessage
        m("Resultat för alla provytor, sökradie {} m{}:".format(
            _fmt_m(radius), ", högst {} % höjdskillnad".format(_fmt_m(max_dh_pct)) if use_dh else ""))
        m("    Referensträd {}, detekterade inom provytorna {}, par {}.".format(
            _fmt_count(summary["n_ref"]), _fmt_count(summary["n_det"]), _fmt_count(summary["hits"])))
        m("    Träffsäkerhet (andel referensträd som hittades): {} (95 % KI {}-{}).".format(
            _fmt_pct(rec), _fmt_pct(rl), _fmt_pct(rh)))
        m("    Precision (andel detekterade som är riktiga): {} (95 % KI {}-{}).".format(
            _fmt_pct(prec), _fmt_pct(pl), _fmt_pct(ph)))
        m("    F1: {}.".format(_fmt_m(round(f1, 3)) if f1 is not None else "-"))
        if summary["ref_per_ha"] is not None:
            m("    Stamantal per ha: referens {}, detekterade {}.".format(
                _fmt_m(round(summary["ref_per_ha"], 1)), _fmt_m(round(summary["det_per_ha"], 1))))
        if summary["rmse_xy"] is not None:
            m("    Lägesfel för paren: medel dx {} m, dy {} m, RMSE {} m.".format(
                _fmt_m(round(summary["bias_x"], 2)), _fmt_m(round(summary["bias_y"], 2)),
                _fmt_m(round(summary["rmse_xy"], 2))))
        if summary["rmse_h"] is not None:
            m("    Höjdfel för paren (detekterad - referens): medel {} m, RMSE {} m.".format(
                _fmt_m(round(summary["bias_h"], 2)), _fmt_m(round(summary["rmse_h"], 2))))
        m("    Referensträd med flera detekterade inom radien (möjlig uppdelad krona): {}. "
          "Detekterade med flera referensträd inom radien (möjlig sammanslagning): {}.".format(
              summary["split"], summary["merge"]))
        if ref_hf:
            for lo, hi in HEIGHT_CLASSES:
                n_c = h_c = 0
                for d, pairs in zip(data, all_pairs):
                    inc = (d["r_h"] >= lo) & (d["r_h"] < hi)
                    n_c += int(inc.sum())
                    h_c += sum(1 for i, _j, _ in pairs if inc[i])
                if n_c:
                    cl, ch = _wilson(h_c, n_c)
                    m("    Referensträd {}: {} av {} hittade, {} (95 % KI {}-{}).".format(
                        "under {} m".format(_fmt_m(hi)) if lo == 0 else
                        "{} m och högre".format(_fmt_m(lo)) if hi > 1e8 else
                        "{}-{} m".format(_fmt_m(lo), _fmt_m(hi)),
                        h_c, n_c, _fmt_pct(h_c / n_c), _fmt_pct(cl), _fmt_pct(ch)))
        m("Känslighet för sökradien:")
        for rad in sorted(set(SWEEP_RADII) | {radius}):
            pr = evaluate(rad)
            hr = sum(len(x) for x in pr)
            hd = sum(1 for d, x in zip(data, pr) for _i, j, _ in x if d["d_in"][j])
            r2, p2, f2 = _scores(tot["n_ref"], tot["n_det"], hr, hd)
            m("    {} m: träffsäkerhet {}, precision {}, F1 {}.".format(
                _fmt_m(rad), _fmt_pct(r2), _fmt_pct(p2),
                _fmt_m(round(f2, 3)) if f2 is not None else "-"))

        # ── 3. Utdata ────────────────────────────────────────────────────────
        steps.next("skriver utdata")
        tmp_fc = _memory_name("parade")
        temps.append(tmp_fc)
        arcpy.management.CreateFeatureclass(os.path.dirname(tmp_fc), os.path.basename(tmp_fc),
                                            "POINT", spatial_reference=_tm())
        arcpy.management.AddFields(tmp_fc, [
            ["TYP", "TEXT", "Typ", 12],
            ["STATUS", "TEXT", "Status", 20],
            ["PROVYTA", "TEXT", "Provyta", 64],
            ["KALL_OID", "LONG", "Objekt-ID i indata"],
            ["PAR_OID", "LONG", "Parat objekt-ID"],
            ["AVSTAND", "DOUBLE", "Avstånd till paret (m)"],
            ["DX", "DOUBLE", "dx till paret (m)"],
            ["DY", "DOUBLE", "dy till paret (m)"],
            ["HOJD", "DOUBLE", "Höjd (m)"],
            ["HOJD_PAR", "DOUBLE", "Parets höjd (m)"],
        ])

        def num(v):
            return None if v is None or not np.isfinite(v) else round(float(v), 2)

        names = ["SHAPE@XY", "TYP", "STATUS", "PROVYTA", "KALL_OID", "PAR_OID", "AVSTAND",
                 "DX", "DY", "HOJD", "HOJD_PAR"]
        with arcpy.da.InsertCursor(tmp_fc, names) as cur:
            for d, pairs in zip(data, all_pairs):
                by_ref = {i: (j, dist) for i, j, dist in pairs}
                by_det = {j: (i, dist) for i, j, dist in pairs}
                for i in range(len(d["r_oid"])):
                    x, y = d["r_xy"][i]
                    if i in by_ref:
                        j, dist = by_ref[i]
                        cur.insertRow(((x, y), "Referens", STATUS_HIT, d["pid"], int(d["r_oid"][i]),
                                       int(d["d_oid"][j]), round(dist, 2),
                                       round(float(d["d_xy"][j, 0] - x), 2),
                                       round(float(d["d_xy"][j, 1] - y), 2),
                                       num(d["r_h"][i]), num(d["d_h"][j])))
                    else:
                        cur.insertRow(((x, y), "Referens", STATUS_MISS, d["pid"],
                                       int(d["r_oid"][i]), None, None, None, None,
                                       num(d["r_h"][i]), None))
                for j in range(len(d["d_oid"])):
                    x, y = d["d_xy"][j]
                    if j in by_det:
                        i, dist = by_det[j]
                        cur.insertRow(((x, y), "Detekterad", STATUS_HIT, d["pid"],
                                       int(d["d_oid"][j]), int(d["r_oid"][i]), round(dist, 2),
                                       round(float(d["r_xy"][i, 0] - x), 2),
                                       round(float(d["r_xy"][i, 1] - y), 2),
                                       num(d["d_h"][j]), num(d["r_h"][i])))
                    else:
                        status = STATUS_FALSE if d["d_in"][j] else STATUS_OUTSIDE
                        cur.insertRow(((x, y), "Detekterad", status, d["pid"],
                                       int(d["d_oid"][j]), None, None, None, None,
                                       num(d["d_h"][j]), None))
        _export(tmp_fc, out_fc, _out_sr(), messages)

        if out_table:
            tmp_tab = _memory_name("matt")
            temps.append(tmp_tab)
            arcpy.management.CreateTable(os.path.dirname(tmp_tab), os.path.basename(tmp_tab))
            arcpy.management.AddFields(tmp_tab, [
                ["PROVYTA", "TEXT", "Provyta", 64],
                ["AREAL_HA", "DOUBLE", "Areal (ha)"],
                ["N_REF", "LONG", "Referensträd"],
                ["N_DET", "LONG", "Detekterade inom provytan"],
                ["PAR", "LONG", "Par"],
                ["TRAFFSAKERHET", "DOUBLE", "Träffsäkerhet"],
                ["PRECISION", "DOUBLE", "Precision"],
                ["F1", "DOUBLE", "F1"],
                ["REF_PER_HA", "DOUBLE", "Referensträd per ha"],
                ["DET_PER_HA", "DOUBLE", "Detekterade per ha"],
                ["RMSE_XY", "DOUBLE", "Lägesfel RMSE (m)"],
                ["BIAS_X", "DOUBLE", "Medel dx (m)"],
                ["BIAS_Y", "DOUBLE", "Medel dy (m)"],
                ["HOJD_BIAS", "DOUBLE", "Medel höjdfel (m)"],
                ["HOJD_RMSE", "DOUBLE", "Höjdfel RMSE (m)"],
                ["SOKRADIE", "DOUBLE", "Sökradie (m)"],
            ])

            def stats(dx, dy, dh):
                dx, dy, dh = np.array(dx), np.array(dy), np.array(dh)
                return (float(np.sqrt(np.mean(dx * dx + dy * dy))) if dx.size else None,
                        float(dx.mean()) if dx.size else None,
                        float(dy.mean()) if dy.size else None,
                        float(dh.mean()) if dh.size else None,
                        float(np.sqrt(np.mean(dh * dh))) if dh.size else None)

            def r3(v):
                return None if v is None else round(v, 3)

            tab_fields = ["PROVYTA", "AREAL_HA", "N_REF", "N_DET", "PAR", "TRAFFSAKERHET",
                          "PRECISION", "F1", "REF_PER_HA", "DET_PER_HA", "RMSE_XY", "BIAS_X",
                          "BIAS_Y", "HOJD_BIAS", "HOJD_RMSE", "SOKRADIE"]
            with arcpy.da.InsertCursor(tmp_tab, tab_fields) as cur:
                for pp in per_plot:
                    a = pp["area_ha"]
                    cur.insertRow((pp["pid"], r3(a), pp["n_ref"], pp["n_det"], pp["hits"],
                                   r3(pp["recall"]), r3(pp["precision"]), r3(pp["f1"]),
                                   r3(pp["n_ref"] / a if a else None),
                                   r3(pp["n_det"] / a if a else None),
                                   *[r3(v) for v in stats(pp["dx"], pp["dy"], pp["dh"])], radius))
                cur.insertRow(("Alla", r3(tot["area"]), summary["n_ref"], summary["n_det"],
                               summary["hits"], r3(rec), r3(prec), r3(f1),
                               r3(summary["ref_per_ha"]), r3(summary["det_per_ha"]),
                               r3(summary["rmse_xy"]), r3(summary["bias_x"]),
                               r3(summary["bias_y"]), r3(summary["bias_h"]),
                               r3(summary["rmse_h"]), radius))
            if arcpy.Exists(out_table):
                _release_from_maps([out_table], messages)
            arcpy.conversion.ExportTable(tmp_tab, out_table)
        steps.done()
        return summary
    finally:
        arcpy.ResetProgressor()
        _delete_quietly(*temps)
