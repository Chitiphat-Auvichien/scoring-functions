# vscore — vibrational mode scoring, on the web

Upload a frequency calculation, get seven scores and a label for every normal
mode, and watch any mode animate in 3D.

Seven scores per mode: `s[Tx] s[Ty] s[Tz]` (translation), `s[Rx] s[Ry] s[Rz]`
(rotation), `s[V_S]` (stretch character). Then a label — **S** stretching,
**B** bending, **SB** mixed — with the largest |score| highlighted in each row.

The scoring core is a verbatim copy of the reference implementation in
`../scoring-functions`. All 24 molecules with published results reproduce them
with a worst-case deviation of **0.00e+00** and no label changes.

## Run locally

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m uvicorn app.main:app --reload --port 8021
```

Then open <http://127.0.0.1:8021>. API docs at `/docs`.

## Input

Three paths, all converging on the same four values.

| Path | Files |
|---|---|
| Gaussian | `.log` **+** `.com`/`.gjf` |
| All-in-one | `.vsc` |
| Convert | `.log` + `.com` → `.vsc`, from the results page or `make_vsc.py` |

**Both Gaussian files are required.** A `.log` contains no connectivity — the
bond list lives in the input deck, which must have `geom=connectivity` in its
route line. `s[V_S]` is a sum over bonds, so this app asks for them rather than
guessing. There is no distance cutoff anywhere in the code.

That refusal is not fussiness. With an empty bond list `s[V_S]` evaluates to
exactly `0.000` for every mode, and `0.000 ≤ τ_B`, so the entire table comes
back labelled **B** — a complete, plausible, silently wrong answer.

### The `.vsc` format

```
#VSCORE 1.0
#TITLE   HOCl  MP2/3-21G

[GEOMETRY] Angstrom
    1  O        0.038069       1.197522       0.000000
    2  H       -0.951724       1.378354       0.000000
    3  Cl       0.038069      -0.644619       0.000000

[CONNECTIVITY]
    1  2  3
    2
    3

[MODES] 3
  mode 1   freq=667.6406   mu=17.2105   k=4.5199   irrep=A'
    1      0.00088      0.85015     -0.00000
    2     -0.16190      0.30478     -0.00000
    3      0.00426     -0.39765      0.00000
```

Element symbols or atomic numbers; leading indices optional; `#` comments and
blank lines ignored; bond orders accepted and ignored. Full spec at `/format`.
Samples in `examples/`.

Typically **23× smaller** than the logs it came from (9.4 MB of Gaussian output
→ 408 kB of `.vsc` across the reference set), because the parts that dominate a
log — SCF iterations, the Hessian, integrals — are never read.

```bash
python make_vsc.py mol.log mol.com -o mol.vsc
python make_vsc.py --batch data/logs data/gjf -o vsc/
```

## What the scorer reads

| Data | Feeds | Needed for |
|---|---|---|
| element per atom | atomic mass | **all 7** |
| coordinates | COM, principal axes, bond vectors | **all 7** |
| displacements | the mode | **all 7** |
| bond list | per-bond terms | **`s[V_S]` only** |

Never read: frequency (a row label — zeroing it changes no score), reduced
mass, force constant, irrep, IR/Raman intensities, charge, multiplicity, basis
set, and the Cartesian Hessian — the largest block in a VEDA `.fmt`.

Element identity is load-bearing for **all seven** scores, not just `V_S`:
masses set the centre of mass and the principal axes, and the R-scores are
measured in that frame. Falsifying HOCl to all-carbon moves `Tx` from −0.300 to
+0.328.

## Two things worth knowing

**Displacement precision.** Write `freq=hpmodes` (5 dp). Standard 2-dp output
leaves the labels intact (`V_S` moves ≤0.03 against a stretch/bend gap of 0.73)
but shifts individual **T-scores by up to 0.31** on a [−1, 1] scale, because
`Tscore()` unit-normalises each atom's displacement before summing — so an atom
printed as `0.00 0.00 0.00` is either masked to zero or promoted to a full unit
vector on rounding noise. The app detects low precision and warns.

**Alignment is not idempotent.** `MIT()`'s sign-fix heuristic negates the
rotation when the heaviest atom's projected coordinates sum negative, so
re-aligning an already-aligned geometry can flip axes and invert T/R scores.
This is why `.vsc` stores the **original** frame, not the scored one. Guarded by
`test_mit_is_not_idempotent` and `test_webapp_vsc_download_roundtrips`.

## Visualization

3Dmol.js, in the browser. `vibrate(numFrames, amplitude, bothWays, arrowSpec)`
reads `dx/dy/dz` per atom — which *is* the parsed displacement vector, used
directly with no conversion — and the viewer is handed **our** bond list rather
than its distance-based guess, so the picture and the V-score are drawn from the
same connectivity.

**Not PyMOL, deliberately.** Tested: the PyPI macOS wheel is broken (hardcoded
rpath to the packager's own machine, `/Users/Martin/.local/share/mamba/...`),
and **no Linux wheel exists at all**, so there is no pip route to PyMOL on a
Linux host. Vercel functions additionally have no GL context, a 250 MB cap, and
are stateless — a working install would still only yield server-rendered stills,
one round trip per frame. 3Dmol.js is interactive, needs no server, and behaves
identically in both deployments.

## Architecture

Fully synchronous. Scoring naphthalene (18 atoms, 54 modes) takes **11 ms**, so
there is no job queue, no progress polling, no server-side run state and nothing
written to disk. CSV and `.vsc` downloads are generated at score time, embedded
in the result page as JSON, and turned into Blobs by the browser.

```
app/
  config.py            paths, constants, upload cap
  main.py              routes: pages, JSON API
  core/
    scoring.py         ) verbatim from scoring-functions,
    classifier.py      ) import paths adjusted only
    utils.py           )
    thresholds.json    ) calibrated τ_TR / τ_S / τ_B
    parsers.py         .log / .com / .vsc readers + .vsc writer
    pipeline.py        parse -> align -> score -> classify
  templates/           Jinja pages
  static/              app.css, app.js, 3Dmol-min.js (vendored, no CDN)
make_vsc.py            batch converter
tests/                 regression against the published results
```

Dependencies: FastAPI, Jinja2, python-multipart, numpy, scipy. **No pandas**
(stdlib `csv` writes the table), **no matplotlib** (the browser draws), no
PyMOL. ≈136 MB against Vercel's 250 MB cap. `scipy` earns its 92 MB in exactly
one place — `linear_sum_assignment` in `classifier.py` — and hand-rolling a
Hungarian solver to save it would risk silently changing labels.

## The vendored scoring core

`app/core/{scoring,classifier,utils}.py` and `app/core/thresholds.json` are
copies of `../src/{scoring,classifier,utils}.py` and
`../data/results/thresholds.json`. They are duplicated rather than imported so
that `api/index.py` is self-contained and deploys to a serverless function
without bundling the whole repository.

**That duplication is a drift hazard.** As of this commit the copies differ from
the originals in exactly two lines, both in `classifier.py`:

```
src/classifier.py                         app/core/classifier.py
  from src.scoring import ...          ->   from .scoring import ...
  DEFAULT_CALIBRATION_PATH =               DEFAULT_CALIBRATION_PATH =
    os.path.join("data", "results",          os.path.join(os.path.dirname(
                 "thresholds.json")            __file__), "thresholds.json")
```

`scoring.py`, `utils.py` and `thresholds.json` are byte-identical.

To re-sync after changing `src/`:

```bash
cp src/{scoring,utils}.py vscore-webapp/app/core/
cp data/results/thresholds.json vscore-webapp/app/core/
# classifier.py needs the two import lines re-applied by hand
cd vscore-webapp && .venv/bin/python -m pytest tests/ -q
```

The test suite is what actually guards this: it compares every score and label
against `../data/results/*_normal.csv`, so a stale copy fails loudly rather
than silently returning different numbers.

## Tests

```bash
.venv/bin/python -m pytest tests/ -q     # 79 passed
```

Compares every score and label against `scoring-functions/data/results/*.csv`,
round-trips every `.vsc` (including the one the results page hands the user),
and asserts the failure modes stay failures: missing connectivity, atom-count
mismatch, a log with no frequency job, an unknown element, low precision.

## Deploy

**Vercel.** `vercel.json` builds `api/index.py` with `@vercel/python` and routes
everything to it. The whole scoring core is vendored under `app/`, so the
function is self-contained.

```bash
npx vercel login
npx vercel --prod
```

Every feature works there — no statefulness is required anywhere. The only
platform limit that bites is the ~4.5 MB request body cap; the app refuses
oversized uploads with a message pointing at `.vsc` conversion.

**Anything long-running.** `uvicorn app.main:app --host 0.0.0.0 --port $PORT`.
