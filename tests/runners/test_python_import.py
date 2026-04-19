"""Python import + infrastructure smoke test round-trip."""

from __future__ import annotations

import math
import sys

import nextdftb


def main() -> int:
    print(f"nextdftb version: {nextdftb.version()}")
    nextdftb.init()

    try:
        out = nextdftb.test()
        # Closed-form: sum_{i=1..1000} i = 1000 * 1001 / 2 = 500500
        expected = 500500.0
        if not math.isclose(out, expected, rel_tol=1e-12):
            print(f"test mismatch: got {out}, expected {expected}", file=sys.stderr)
            return 1
    finally:
        nextdftb.finalize()

    return 0


if __name__ == "__main__":
    sys.exit(main())
