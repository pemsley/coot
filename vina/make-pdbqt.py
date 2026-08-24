
"""Prepare receptor and ligand PDBQT files for AutoDock Vina using Coot.

  - reads a receptor PDB file and writes a rigid receptor .pdbqt
  - reads a ligand restraints dictionary (.cif), builds the ligand from it,
    and writes a flexible-ligand .pdbqt (with a torsion tree)

Usage:
    python make-pdbqt.py <receptor.pdb> <ligand.cif> [comp_id]

    comp_id is the ligand three-letter code (default "LIG").

Outputs (next to the inputs):
    <receptor>.pdbqt
    <ligand>.pdbqt

The coot_headless_api module must be importable: either install it, put it on
PYTHONPATH, or set COOT_HEADLESS_API_DIR to the build directory that contains
coot_headless_api*.so.
"""

import os
import sys


def import_chapi():
    # Prefer an explicit build directory over any installed copy.
    d = os.environ.get("COOT_HEADLESS_API_DIR")
    if d and os.path.isdir(d):
        sys.path.insert(0, d)
    try:
        import coot_headless_api
        return coot_headless_api
    except ImportError as e:
        sys.exit("cannot import coot_headless_api (%s); set COOT_HEADLESS_API_DIR "
                 "or PYTHONPATH" % e)


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: make-pdbqt.py <receptor.pdb> <ligand.cif> [comp_id]")

    receptor_pdb = sys.argv[1]
    ligand_cif   = sys.argv[2]
    comp_id      = sys.argv[3] if len(sys.argv) > 3 else "LIG"

    receptor_out = os.path.splitext(receptor_pdb)[0] + ".pdbqt"
    ligand_out   = os.path.splitext(ligand_cif)[0]   + ".pdbqt"

    chapi = import_chapi()
    mc = chapi.molecules_container_t(True)
    mc.set_use_gemmi(False)
    imol_enc_any = mc.get_imol_enc_any()

    # ---------------------------------------------------------------- receptor
    imol_rec = mc.read_pdb(receptor_pdb)
    if imol_rec < 0:
        sys.exit("failed to read receptor: " + receptor_pdb)

    # The receptor typing (polar HD) and charges need hydrogens; add them with
    # Coot's own reduce only if the model doesn't already have them.
    if mc.get_number_of_hydrogen_atoms(imol_rec) == 0:
        mc.add_hydrogen_atoms(imol_rec)

    n_rec = mc.export_molecule_as_pdbqt(imol_rec, receptor_out)
    if n_rec > 0:
        print("receptor: wrote %d atoms to %s" % (n_rec, receptor_out))
    else:
        sys.exit("failed to write receptor PDBQT")

    # ------------------------------------------------------------------ ligand
    mc.import_cif_dictionary(ligand_cif, imol_enc_any)
    imol_lig = mc.get_monomer(comp_id)
    if imol_lig < 0:
        sys.exit("failed to build monomer '%s' from %s "
                 "(is that the right three-letter code?)" % (comp_id, ligand_cif))

    n_lig = mc.export_ligand_as_pdbqt(imol_lig, "/*/*/*", ligand_out)
    if n_lig > 0:
        print("ligand:   wrote %d atoms to %s" % (n_lig, ligand_out))
    else:
        sys.exit("failed to write ligand PDBQT")


if __name__ == "__main__":
    main()
