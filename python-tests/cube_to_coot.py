#!/usr/bin/env python3
"""Read a Gaussian/ORCA cube file and write out a model (PDB) and a map (CCP4)
so they can be opened in GUI Coot to check the result.

Usage:
    python3 cube_to_coot.py file.cube [output_prefix]

Writes <output_prefix>.pdb and <output_prefix>.map
(default output_prefix is the cube file name without its extension).

Then in Coot:  File -> Open Coordinates  (the .pdb)
               File -> Open Map          (the .map)
An orbital cube is written as a difference map, so it shows both lobes.
"""

import os
import sys

# Make the headless-api module importable. Set COOT_CHAPI_BUILD_DIR to override,
# otherwise try the build dir used in this project.
_build_dir = os.environ.get(
    "COOT_CHAPI_BUILD_DIR",
    os.path.expanduser("~/Projects/coot/git/coot-main/build-chapi-for-claude"))
if _build_dir and _build_dir not in sys.path:
    sys.path.insert(0, _build_dir)

import coot_headless_api as chapi


def cube_to_coot(cube_file_name, output_prefix):
    pdb_file_name = output_prefix + ".pdb"
    map_file_name = output_prefix + ".map"

    mc = chapi.molecules_container_t(False)
    imol_model, imol_map = mc.read_cube(cube_file_name)

    if imol_model == -1 and imol_map == -1:
        print("ERROR:: failed to read cube file %s" % cube_file_name)
        return False

    if mc.is_valid_model_molecule(imol_model):
        mc.write_coordinates(imol_model, pdb_file_name)
        print("Wrote model: %s" % pdb_file_name)
    else:
        print("WARNING:: no valid model molecule from %s" % cube_file_name)

    if mc.is_valid_map_molecule(imol_map):
        is_diff = mc.is_a_difference_map(imol_map)
        mc.write_map(imol_map, map_file_name)
        print("Wrote map:   %s   (difference map: %s)" % (map_file_name, is_diff))
    else:
        print("WARNING:: no valid map molecule from %s "
              "(non-orthogonal grid?)" % cube_file_name)

    return True


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    cube_file_name = sys.argv[1]
    if len(sys.argv) > 2:
        output_prefix = sys.argv[2]
    else:
        output_prefix = os.path.splitext(cube_file_name)[0]
    ok = cube_to_coot(cube_file_name, output_prefix)
    sys.exit(0 if ok else 1)
