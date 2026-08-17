"""Vercel serverless entrypoint.

@vercel/python serves an ASGI app exported as ``app``. The whole scoring core
is vendored under ``app/core/``, so nothing outside this directory is needed
at runtime and the function is self-contained.

Unlike a long-running reservoir job, scoring is ~11 ms and fully synchronous
with no server-side state and no disk writes, so every feature of this app
works identically here and on localhost.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from app.main import app  # noqa: E402

__all__ = ["app"]
