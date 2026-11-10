
#!/usr/bin/env python3
"""Flag bedpe files with both high (>=98) and low (90-96) ANI bins as recombination events."""
import glob, os

HIGH, LOW_MIN, LOW_MAX = 98.0, 90.0, 96.0
import sys
MIN_HIGH_FRAC, MIN_LOW_FRAC = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (0.2, 0.1)
OUT = sys.argv[3] if len(sys.argv) > 3 else "recombination.txt"

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, OUT), "w") as out:
    for d in sorted(os.listdir(os.path.join(here, "dotplots2"))):
        files = glob.glob(os.path.join(here, "dotplots2", d, "*", "*.bedpe"))
        n = hi = lo = 0
        for fn in files:
            with open(fn) as fh:
                next(fh, None)
                for line in fh:
                    f = line.split("\t")
                    # Skip self-comparisons
                    if f[1] == f[4] and f[2] == f[5]:
                        continue
                    v = float(f[6])
                    n += 1
                    if v >= HIGH:
                        hi += 1
                    elif LOW_MIN <= v <= LOW_MAX:
                        lo += 1
        flag = int(n > 0 and hi / n >= MIN_HIGH_FRAC and lo / n >= MIN_LOW_FRAC)
        out.write(f"{d}\t{flag}\n")