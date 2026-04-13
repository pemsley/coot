#!/usr/bin/env python3
"""
3D Hilbert curve: convert between (x, y, z) and 1D index
for a 2^order x 2^order x 2^order grid.

Based directly on:
  T. Skilling, "Programming the Hilbert curve",
  AIP Conf. Proc. 707, 381 (2004)
"""


def _transpose_to_axes(X, order, ndim):
    """TransposetoAxes: convert transposed Hilbert index to spatial coordinates."""
    N = 1 << order

    # Gray decode: undo the Gray code
    t = X[ndim - 1] >> 1
    for i in range(ndim - 1, 0, -1):
        X[i] ^= X[i - 1]
    X[0] ^= t

    # Undo excess work
    Q = 2
    while Q != N:
        P = Q - 1
        for i in range(ndim - 1, -1, -1):
            if X[i] & Q:
                X[0] ^= P
            else:
                t = (X[0] ^ X[i]) & P
                X[0] ^= t
                X[i] ^= t
        Q <<= 1


def _axes_to_transpose(X, order, ndim):
    """AxestoTranspose: convert spatial coordinates to transposed Hilbert index."""
    M = 1 << (order - 1)

    # Inverse undo
    Q = M
    while Q > 1:
        P = Q - 1
        for i in range(ndim):
            if X[i] & Q:
                X[0] ^= P
            else:
                t = (X[0] ^ X[i]) & P
                X[0] ^= t
                X[i] ^= t
        Q >>= 1

    # Gray encode
    for i in range(1, ndim):
        X[i] ^= X[i - 1]
    t = 0
    Q = M
    while Q > 1:
        if X[ndim - 1] & Q:
            t ^= Q - 1
        Q >>= 1
    for i in range(ndim):
        X[i] ^= t


def _transpose_to_index(X, order, ndim):
    """Interleave bits from transposed form into a single integer."""
    d = 0
    for i in range(order - 1, -1, -1):
        for dim in range(ndim):
            d = (d << 1) | ((X[dim] >> i) & 1)
    return d


def _index_to_transpose(d, order, ndim):
    """De-interleave a single integer into transposed form."""
    X = [0] * ndim
    for i in range(order):
        for dim in range(ndim - 1, -1, -1):
            X[dim] |= (d & 1) << i
            d >>= 1
    return X


def xyz_to_hilbert_3d(x, y, z, order):
    """Convert (x, y, z) to 3D Hilbert curve index.

    Grid is 2^order on each side.  Returns an integer in [0, 8^order).
    """
    X = [x, y, z]
    _axes_to_transpose(X, order, 3)
    return _transpose_to_index(X, order, 3)


def hilbert_to_xyz_3d(d, order):
    """Convert 3D Hilbert curve index to (x, y, z).

    Grid is 2^order on each side.  d must be in [0, 8^order).
    """
    X = _index_to_transpose(d, order, 3)
    _transpose_to_axes(X, order, 3)
    return X[0], X[1], X[2]


# ------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------

def test_round_trip():
    """Verify xyz->d->xyz and d->xyz->d for all points at several orders."""
    for order in range(1, 5):
        n = 1 << order
        total = n * n * n
        seen = set()

        for x in range(n):
            for y in range(n):
                for z in range(n):
                    d = xyz_to_hilbert_3d(x, y, z, order)
                    assert 0 <= d < total, f"d={d} out of range for order={order}"
                    assert d not in seen, (
                        f"duplicate d={d} at ({x},{y},{z}), order={order}"
                    )
                    seen.add(d)

                    x2, y2, z2 = hilbert_to_xyz_3d(d, order)
                    assert (x2, y2, z2) == (x, y, z), (
                        f"round-trip failed: ({x},{y},{z}) -> d={d} -> ({x2},{y2},{z2})"
                    )

        for d in range(total):
            x, y, z = hilbert_to_xyz_3d(d, order)
            d2 = xyz_to_hilbert_3d(x, y, z, order)
            assert d2 == d, (
                f"inverse round-trip failed: d={d} -> ({x},{y},{z}) -> d2={d2}"
            )

        print(f"order {order}: all {total} points passed round-trip")


def test_adjacency():
    """Check consecutive Hilbert indices are spatial neighbours (Manhattan distance 1)."""
    for order in range(1, 5):
        n = 1 << order
        total = n * n * n
        prev = hilbert_to_xyz_3d(0, order)
        for d in range(1, total):
            cur = hilbert_to_xyz_3d(d, order)
            dist = abs(cur[0]-prev[0]) + abs(cur[1]-prev[1]) + abs(cur[2]-prev[2])
            assert dist == 1, (
                f"order={order}, d={d}: distance from {prev} to {cur} is {dist}, expected 1"
            )
            prev = cur
        print(f"order {order}: adjacency OK for {total} points")


if __name__ == "__main__":
    test_round_trip()
    test_adjacency()
    print("All tests passed.")
