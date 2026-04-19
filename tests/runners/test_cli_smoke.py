"""CLI smoke test — launch the standalone executable on the infra test."""

from __future__ import annotations

import argparse
import subprocess
import sys


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--binary", required=True, help="Path to the nextdftb CLI")
    args = p.parse_args()

    r = subprocess.run(
        [args.binary, "--log-level", "WARN"],
        check=False,
        capture_output=True,
        text=True,
    )
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)

    if r.returncode != 0:
        print(f"CLI exited with {r.returncode}", file=sys.stderr)
        return r.returncode

    if "test =" not in r.stdout:
        print("CLI output does not contain the expected 'test =' line", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
