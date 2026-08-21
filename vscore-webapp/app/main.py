"""main.py -- FastAPI application for the vibrational mode scorer.

Deliberately stateless. Scoring naphthalene (18 atoms, 54 modes) takes ~11 ms,
so every request is handled synchronously: no job registry, no progress
polling, no results directory, nothing written to disk. That is what lets the
same code run unchanged on a normal server and on a read-only serverless
filesystem.

Downloads (CSV, .vsc) are generated at score time and embedded in the result
page, so the browser can produce them from a Blob without a second round trip
and without the server holding any run state.
"""

from __future__ import annotations

import csv
import io
import json

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .config import (IS_LOCAL, MAX_UPLOAD_BYTES, SCORE_KEYS, SITE_NOTE,
                     STATIC_DIR, TEMPLATES_DIR, VERSION, asset_version)
from .core.parsers import (ParseError, parse_connectivity, parse_gaussian_log,
                           parse_vsc, write_vsc)
from .core.pipeline import analyse, to_csv_rows

app = FastAPI(
    title="Vibrational Mode Scorer",
    description=("Scores each normal mode of a molecule against seven reference "
                 "motions (Tx, Ty, Tz, Rx, Ry, Rz, V_S) and classifies it. "
                 + SITE_NOTE),
    version=VERSION,
)
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))
templates.env.globals["version"] = VERSION
templates.env.globals["asset_v"] = asset_version
templates.env.globals["is_local"] = IS_LOCAL
templates.env.globals["site_note"] = SITE_NOTE
templates.env.globals["score_keys"] = SCORE_KEYS


# ----------------------------------------------------------------------
# Upload handling
# ----------------------------------------------------------------------
async def _read(upload: UploadFile | None, what: str) -> str | None:
    """Decode an upload to text, refusing oversized or binary files.

    The size cap matters on serverless hosts, which reject large request
    bodies at the platform layer with an opaque error; failing here gives the
    user a message that says what to do about it.
    """
    if upload is None or not upload.filename:
        return None
    data = await upload.read()
    if len(data) > MAX_UPLOAD_BYTES:
        raise ParseError(
            f"{what} is {len(data)/1e6:.1f} MB, over the {MAX_UPLOAD_BYTES/1e6:.1f} MB "
            "limit. Gaussian logs are mostly SCF iterations the scorer never reads -- "
            "convert to a .vsc file locally with make_vsc.py and upload that instead "
            "(typically a few kB).")
    try:
        return data.decode("utf-8", errors="replace")
    except Exception:
        raise ParseError(f"{what} could not be read as text.")


def _build(mode, log_text, com_text, vsc_text, log_name, com_name, vsc_name,
           pasted=None, arb_text=None, arb_name=None):
    """Route the three input paths onto one payload.

    Returns ``(payload, raw)``. ``raw`` keeps the geometry and displacements
    exactly as parsed -- BEFORE principal-axis alignment -- because that is
    what the .vsc writer must serialise. ``MIT()`` is not idempotent: its
    sign-fix heuristic negates the rotation when the heaviest atom's projected
    coordinates sum negative, so re-aligning an already-aligned geometry can
    flip axes and invert the T/R scores. Writing the original frame means a
    downloaded .vsc re-scores to the same numbers as the files it came from.
    """
    if mode == "paste":
        if not IS_LOCAL:
            raise ParseError(
                "This page is available on a local instance only. Use the "
                ".vsc upload above instead.")
        # Two ways in, one meaning: whatever arrives here is DECLARED to be a
        # selection of vibrations. mode_set is forced rather than detected, so
        # a count that happens to equal 3N or 3N-6 is still treated as a
        # selection and the ideal T/R references are still constructed to fill
        # the external slots.
        text = arb_text if (arb_text or "").strip() else pasted
        if not (text or "").strip():
            raise ParseError("No modes were entered or uploaded.")
        if len(text) > MAX_UPLOAD_BYTES:
            raise ParseError("That is too large; upload it as a file.")
        v = parse_vsc(text)
        raw = {"atoms": v["atoms"], "coords": v["coords"],
               "bonds": v["bonds"], "modes": v["modes"]}
        return analyse(v["atoms"], v["coords"], v["bonds"], v["modes"],
                       title=v["title"] or "arbitrary modes",
                       source=(arb_name or "typed in"),
                       mode_set="arbitrary"), raw

    if mode == "vsc":
        if not vsc_text:
            raise ParseError("No .vsc file was uploaded.")
        v = parse_vsc(vsc_text)
        raw = {"atoms": v["atoms"], "coords": v["coords"],
               "bonds": v["bonds"], "modes": v["modes"]}
        payload = analyse(v["atoms"], v["coords"], v["bonds"], v["modes"],
                          title=v["title"] or _stem(vsc_name),
                          source=vsc_name or "uploaded .vsc")
        return payload, raw

    if not log_text:
        raise ParseError("No Gaussian output file was uploaded.")
    if not com_text:
        raise ParseError(
            "No .com/.gjf input deck was uploaded. A Gaussian .log carries no "
            "connectivity -- the bond list lives in the input deck -- and the "
            "V-score is a sum over bonds, so both files are required.")
    g = parse_gaussian_log(log_text)
    bonds = parse_connectivity(com_text, len(g["atoms"]))
    raw = {"atoms": g["atoms"], "coords": g["coords"],
           "bonds": bonds, "modes": g["modes"]}
    # A Gaussian frequency job prints 3N-6 (or 3N-5) modes; T/R are projected
    # out at a stationary point. Detection confirms that rather than assuming it.
    payload = analyse(g["atoms"], g["coords"], bonds, g["modes"],
                      title=_stem(log_name),
                      source=f"{log_name} + {com_name}")
    return payload, raw


def _stem(name):
    if not name:
        return "molecule"
    return name.rsplit("/", 1)[-1].rsplit(".", 1)[0]


def _csv_text(payload, include_references=False):
    rows = to_csv_rows(payload, include_references)
    if not rows:
        return ""
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def _vsc_text(payload, raw):
    """Serialise the parsed input -- original frame -- to .vsc text.

    Deliberately built from ``raw``, not from the scored payload: the payload
    holds principal-axis-aligned coordinates, and re-aligning those on upload
    can flip axis signs (see ``_build``). Serialising the frame the data
    arrived in makes the round trip exact.
    """
    return write_vsc(raw["atoms"], raw["coords"],
                     [tuple(b) for b in raw["bonds"]], raw["modes"],
                     title=payload["title"],
                     source=payload.get("source", ""))


# ----------------------------------------------------------------------
# Pages
# ----------------------------------------------------------------------
@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(request, "index.html", {})


@app.post("/score", response_class=HTMLResponse)
async def score(
    request: Request,
    mode: str = Form("log"),
    logfile: UploadFile | None = File(None),
    comfile: UploadFile | None = File(None),
    vscfile: UploadFile | None = File(None),
    pasted: str = Form(""),
    arbfile: UploadFile | None = File(None),
):
    try:
        payload, raw = _build(
            mode,
            await _read(logfile, "The Gaussian output file"),
            await _read(comfile, "The input deck"),
            await _read(vscfile, "The .vsc file"),
            logfile.filename if logfile else None,
            comfile.filename if comfile else None,
            vscfile.filename if vscfile else None,
            pasted=pasted,
            arb_text=await _read(arbfile, "The .vsc file"),
            arb_name=arbfile.filename if arbfile else None,
        )
    except ParseError as exc:
        return templates.TemplateResponse(
            request, "index.html", {"error": str(exc)}, status_code=400)
    except Exception as exc:                       # unexpected: still readable
        return templates.TemplateResponse(
            request, "index.html",
            {"error": f"{type(exc).__name__}: {exc}"}, status_code=500)

    return templates.TemplateResponse(request, "result.html", {
        "p": payload,
        "payload_json": _embed(payload),
        "downloads_json": _embed({
            "csv": _csv_text(payload),
            "csvall": _csv_text(payload, include_references=True),
            "vsc": _vsc_text(payload, raw),
        }),
    })


def _embed(obj):
    """JSON for inlining in a <script> block.

    Downloadable text goes through JSON rather than straight into the element
    body: Jinja would otherwise HTML-escape it, turning an irrep like A' into
    A&#39; inside the file the user saves. Escaping '</' additionally stops any
    payload content from closing the script element early.
    """
    return json.dumps(obj).replace("</", "<\\/")


@app.get("/about", response_class=HTMLResponse)
def about(request: Request):
    return templates.TemplateResponse(request, "about.html", {})


@app.get("/format", response_class=HTMLResponse)
def format_page(request: Request):
    return templates.TemplateResponse(request, "format.html", {})


# ----------------------------------------------------------------------
# JSON API
# ----------------------------------------------------------------------
@app.get("/api/health")
def api_health():
    return {"status": "ok", "version": VERSION, "scores": list(SCORE_KEYS)}


@app.post("/api/score")
async def api_score(
    mode: str = Form("log"),
    logfile: UploadFile | None = File(None),
    comfile: UploadFile | None = File(None),
    vscfile: UploadFile | None = File(None),
    pasted: str = Form(""),
    arbfile: UploadFile | None = File(None),
):
    try:
        payload, _raw = _build(
            mode,
            await _read(logfile, "The Gaussian output file"),
            await _read(comfile, "The input deck"),
            await _read(vscfile, "The .vsc file"),
            logfile.filename if logfile else None,
            comfile.filename if comfile else None,
            vscfile.filename if vscfile else None,
            pasted=pasted,
            arb_text=await _read(arbfile, "The .vsc file"),
            arb_name=arbfile.filename if arbfile else None,
        )
    except ParseError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    return JSONResponse(payload)
