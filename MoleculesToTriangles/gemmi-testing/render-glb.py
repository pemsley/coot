#!/usr/bin/env python3
"""Minimal headless .glb viewer: parse a binary glTF (as written by tinygltf),
pull POSITION / COLOR_0 / indices out of the first mesh primitive and render the
triangle mesh to a PNG with matplotlib. No external glTF library needed.

Usage: render-glb.py <in.glb> <out.png>
"""
import sys, json, struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

COMP_TYPE = {5120: ("b", 1), 5121: ("B", 1), 5122: ("h", 2), 5123: ("H", 2),
             5125: ("I", 4), 5126: ("f", 4)}
NCOMP = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4}


def load_glb(path):
    with open(path, "rb") as f:
        data = f.read()
    magic, version, length = struct.unpack_from("<III", data, 0)
    assert magic == 0x46546C67, "not a glb"
    off, jchunk, bchunk = 12, None, None
    while off < length:
        clen, ctype = struct.unpack_from("<II", data, off)
        off += 8
        chunk = data[off:off + clen]
        off += clen
        if ctype == 0x4E4F534A:   # JSON
            jchunk = json.loads(chunk)
        elif ctype == 0x004E4942:  # BIN
            bchunk = chunk
    return jchunk, bchunk


def accessor(gltf, bin_blob, idx):
    acc = gltf["accessors"][idx]
    bv = gltf["bufferViews"][acc["bufferView"]]
    fmt, size = COMP_TYPE[acc["componentType"]]
    ncomp = NCOMP[acc["type"]]
    base = bv.get("byteOffset", 0) + acc.get("byteOffset", 0)
    count = acc["count"]
    arr = np.frombuffer(bin_blob, dtype=np.dtype("<" + fmt),
                        count=count * ncomp, offset=base)
    return arr.reshape(count, ncomp) if ncomp > 1 else arr


def main():
    in_glb, out_png = sys.argv[1], sys.argv[2]
    gltf, blob = load_glb(in_glb)

    verts, cols, tris = [], [], []
    voff = 0
    for mesh in gltf.get("meshes", []):
        for prim in mesh["primitives"]:
            attr = prim["attributes"]
            p = accessor(gltf, blob, attr["POSITION"])
            idx = accessor(gltf, blob, prim["indices"]).astype(np.int64)
            if "COLOR_0" in attr:
                c = accessor(gltf, blob, attr["COLOR_0"]).astype(np.float64)
                if c.max() > 1.5:
                    c = c / 255.0
            else:
                c = np.tile([0.6, 0.6, 0.9, 1.0], (len(p), 1))
            verts.append(p)
            cols.append(c[:, :3])
            tris.append(idx.reshape(-1, 3) + voff)
            voff += len(p)

    V = np.vstack(verts)
    C = np.vstack(cols)
    T = np.vstack(tris)
    print(f"{in_glb}: {len(V)} verts, {len(T)} triangles")

    # one face colour = mean of its 3 vertex colours
    face_cols = C[T].mean(axis=1)

    ctr = V.mean(axis=0)
    rad = np.linalg.norm(V - ctr, axis=1).max()

    views = [(20, -60), (20, 30), (90, -90)]
    fig = plt.figure(figsize=(15, 5))
    for i, (elev, azim) in enumerate(views):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        coll = Poly3DCollection(V[T], facecolors=face_cols, edgecolors="none",
                                linewidths=0, shade=False)
        ax.add_collection3d(coll)
        ax.set_xlim(ctr[0] - rad, ctr[0] + rad)
        ax.set_ylim(ctr[1] - rad, ctr[1] + rad)
        ax.set_zlim(ctr[2] - rad, ctr[2] + rad)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(f"elev={elev} azim={azim}", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110, bbox_inches="tight")
    print("wrote", out_png)


if __name__ == "__main__":
    main()
