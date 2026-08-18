"""config.py -- paths and constants.

No sys.path wiring: the scoring core is vendored under ``app/core/`` rather
than imported from a sibling tree, so the app is self-contained and deploys
as one unit. ``app/core/{scoring,classifier,utils}.py`` and ``thresholds.json``
are verbatim copies from the reference implementation (scoring-functions), with
only import paths adjusted -- the numbers are checked against its published
results by tests/test_regression.py.
"""

from __future__ import annotations

from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
STATIC_DIR = APP_DIR / "static"
TEMPLATES_DIR = APP_DIR / "templates"

VERSION = "1.0"


def asset_version():
    """Cache-buster for app.css / app.js, from their modification times.

    Browsers hold on to these aggressively: an edited stylesheet kept rendering
    the previous layout until a hard refresh, which reads as "the CSS is
    broken" rather than "the CSS is cached". Appending ?v=<mtime> makes the URL
    change whenever the file does, so a normal reload picks it up.
    """
    stamp = 0
    for name in ("app.css", "app.js"):
        f = STATIC_DIR / name
        if f.exists():
            stamp = max(stamp, int(f.stat().st_mtime))
    return str(stamp)

# The seven scores, in display order.
SCORE_KEYS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_S")

# Serverless platforms cap request bodies (Vercel at ~4.5 MB). Refusing here
# gives a message that names the fix; the platform's own rejection does not.
MAX_UPLOAD_BYTES = 4_000_000

SITE_NOTE = (
    "Scores are computed by the reference implementation verbatim: elements, "
    "geometry, displacements and an explicit bond list in, seven scores and a "
    "label out. Connectivity is never guessed -- it comes from your input deck "
    "or your .vsc file, and is shown back to you on every result."
)
