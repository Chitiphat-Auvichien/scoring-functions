#!/usr/bin/env python3
"""make_vsc.py -- convert Gaussian .log + .com pairs into .vsc files.

The webapp offers a .vsc download on every result page; this is the same
conversion for batch work, e.g. on the compute server.

    python make_vsc.py HOCl.log HOCl.com -o HOCl.vsc
    python make_vsc.py --batch data/logs data/gjf -o vsc/

In batch mode every ``<stem>.log`` is paired with ``<stem>.com`` or
``<stem>.gjf`` in the second directory. A log with no matching deck is
reported and skipped -- never scored with guessed bonds.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from app.core.parsers import (ParseError, parse_connectivity,  # noqa: E402
                              parse_gaussian_log, write_vsc)


def convert(log_path: Path, com_path: Path) -> str:
    g = parse_gaussian_log(log_path.read_text(errors="replace"))
    bonds = parse_connectivity(com_path.read_text(errors="replace"), len(g["atoms"]))
    return write_vsc(g["atoms"], g["coords"], bonds, g["modes"],
                     title=log_path.stem,
                     source=f"{log_path.name} + {com_path.name}")


def find_deck(gjf_dir: Path, stem: str) -> Path | None:
    for ext in (".com", ".gjf", ".inp"):
        p = gjf_dir / f"{stem}{ext}"
        if p.exists():
            return p
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", nargs="?", help="Gaussian output file")
    ap.add_argument("com", nargs="?", help="Gaussian input deck (connectivity)")
    ap.add_argument("-o", "--out", help="output .vsc file, or directory in --batch mode")
    ap.add_argument("--batch", nargs=2, metavar=("LOGDIR", "GJFDIR"),
                    help="convert every .log in LOGDIR paired with its deck in GJFDIR")
    args = ap.parse_args(argv)

    if args.batch:
        log_dir, gjf_dir = Path(args.batch[0]), Path(args.batch[1])
        out_dir = Path(args.out or "vsc")
        out_dir.mkdir(parents=True, exist_ok=True)
        ok = skipped = failed = 0
        for log in sorted(log_dir.glob("*.log")):
            deck = find_deck(gjf_dir, log.stem)
            if deck is None:
                print(f"  skip  {log.name}: no .com/.gjf deck -- connectivity is "
                      f"required and is never inferred")
                skipped += 1
                continue
            try:
                text = convert(log, deck)
            except ParseError as exc:
                print(f"  FAIL  {log.name}: {exc}")
                failed += 1
                continue
            dest = out_dir / f"{log.stem}.vsc"
            dest.write_text(text)
            print(f"  ok    {log.stem:20s} -> {dest}  ({len(text)/1000:.1f} kB)")
            ok += 1
        print(f"\n{ok} converted, {skipped} skipped (no deck), {failed} failed")
        return 1 if failed else 0

    if not args.log or not args.com:
        ap.error("give both a .log and a .com, or use --batch LOGDIR GJFDIR")

    try:
        text = convert(Path(args.log), Path(args.com))
    except ParseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {args.out} ({len(text)/1000:.1f} kB)")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
