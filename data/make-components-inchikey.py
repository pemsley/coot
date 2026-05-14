
"""
Build Components-inchikey.ich from a CCD components.cif file.

Output format (tab-separated, lowercase name):

    <InChIKey>\t<comp_id>\t<name>

Usage:
    make-components-inchikey.py components.cif > Components-inchikey.ich

Requires gemmi. The Coot chapi build has it:
    /home/paule/autobuild/build-for-chapi-arch/bin/python3
"""
import re
import sys
import gemmi


def normalise_name(s: str) -> str:
    """Collapse CIF text-block line wraps into a single-line name.

    Line wraps adjacent to a hyphen are treated as in-token continuations
    (e.g. "1H-\\ninden" → "1H-inden"); other wraps become a single space.
    """
    s = s.strip()
    s = re.sub(r"[ \t]*\n[ \t]*", " ", s)
    s = re.sub(r"[ \t]+", " ", s)
    # Treat a space adjacent to a hyphen as a wrap artifact; chemical names
    # don't have real spaces around hyphens.
    s = s.replace("- ", "-").replace(" -", "-")
    return s.lower()


def main(path: str) -> int:
    doc = gemmi.cif.read(path)
    for block in doc:
        comp_id = block.find_value("_chem_comp.id")
        name = block.find_value("_chem_comp.name")
        if comp_id is None or name is None:
            continue
        comp_id = gemmi.cif.as_string(comp_id).strip()
        name = normalise_name(gemmi.cif.as_string(name))

        inchikey = None
        table = block.find("_pdbx_chem_comp_descriptor.",
                           ["type", "descriptor"])
        for row in table:
            if gemmi.cif.as_string(row[0]) == "InChIKey":
                inchikey = gemmi.cif.as_string(row[1]).strip()
                break
        if not inchikey:
            continue

        sys.stdout.write(f"{inchikey}\t{comp_id}\t{name}\n")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write(f"usage: {sys.argv[0]} components.cif\n")
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
