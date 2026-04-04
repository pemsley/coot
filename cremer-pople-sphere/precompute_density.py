#!/usr/bin/env python3
"""
Precompute a conformation-density grid from a Cremer-Pople parameters CSV and
write it to a compact binary file that can be loaded in microseconds.

Usage:
    python precompute_density.py <input.csv> <output.bin> [lat=128] [lon=128]

Binary format (all little-endian):
    8 bytes       magic  "CPDENS\\x00\\x00"
    4 bytes       uint32 version = 1
    4 bytes       int32  lat
    4 bytes       int32  lon
    lat*lon*4     float32 values, log1p-normalised density, row-major
                  index = i*lon + j  where i=latitude bin (0=north), j=longitude bin
"""

import array
import csv
import math
import struct
import sys


def precompute(csv_path: str, out_path: str, lat: int = 128, lon: int = 128) -> None:
    counts = [0.0] * (lat * lon)
    n_loaded = n_skipped = 0

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            theta_deg = float(row["theta"])
            phi_deg   = float(row["phi"])

            # Sentinel value used for missing / invalid observations
            if theta_deg < 0.0 or phi_deg < 0.0:
                n_skipped += 1
                continue

            i = int(theta_deg / 180.0 * lat)
            j = int(phi_deg   / 360.0 * lon)
            i = max(0, min(i, lat - 1))
            j = max(0, min(j, lon - 1))
            counts[i * lon + j] += 1.0
            n_loaded += 1

    print(f"Loaded {n_loaded:,} conformations, skipped {n_skipped:,} invalid rows.")

    # Log1p-normalise so sparse and dense regions are both visible
    max_count = max(counts)
    if max_count > 0.0:
        log_max = math.log1p(max_count)
        counts = [math.log1p(c) / log_max for c in counts]

    # Write binary file
    MAGIC = b"CPDENS\x00\x00"
    with open(out_path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(struct.pack("<I", 1))           # version
        fh.write(struct.pack("<ii", lat, lon))   # grid dimensions
        fh.write(array.array("f", counts))       # float32 payload

    size_kb = (len(MAGIC) + 4 + 8 + lat * lon * 4) / 1024
    print(f"Written {lat}×{lon} grid ({size_kb:.1f} KB) → {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    _lat = int(sys.argv[3]) if len(sys.argv) > 3 else 128
    _lon = int(sys.argv[4]) if len(sys.argv) > 4 else 128
    precompute(sys.argv[1], sys.argv[2], _lat, _lon)
