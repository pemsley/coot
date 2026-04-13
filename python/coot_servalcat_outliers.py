import json
import re
import os

# this script reads the JSON output from Servalcat refinement (prefix_mod_refined_stats.json) and generates a new JSON file formatted for Coot, 
# containing outlier information for bonds, angles, torsions, chirality, planes, and van der Waals clashes.
# It reports outliers from the last refinement cycle.

def convert_outliers_json_inner(json_file_name, output_json_file_name):

    with open(json_file_name, 'r') as file:
        data = json.load(file)

    output = {
        "title": "Servalcat Validation",
        "sections": []
    }

    def parse_atom_string(atom_str):
        # Example: "B/SER 282I/CB" or "A/ZN 360/ZN.A"
        chain, rest = atom_str.split("/", 1)
        res_part, atom_part = rest.split("/")
        resname, resseq_raw = res_part.split()

        # Split resseq and insertion code (e.g., "282I" → 282, "I")
        match = re.match(r"(\d+)([A-Za-z]?)", resseq_raw)
        if match:
            resseq = int(match.group(1))
            icode = match.group(2)
        else:
            raise ValueError(f"Invalid residue format: {resseq_raw}")

        # Handle atom name + altloc
        if "." in atom_part:
            name, altloc = atom_part.split(".", 1)
        else:
            name = atom_part
            altloc = ""

        return {
            "chain_id": chain,
            "resname": resname,
            "resseq": resseq,
            "icode": icode,
            "atom_name": name,
            "altloc": altloc
        }


    def get_atom_spec_as_string(obj):
        parts = [
            obj["chain_id"],
            obj["resname"],
            str(obj["resseq"]),
            obj["icode"] if obj["icode"] else "",
            obj["atom_name"],
            obj["altloc"] if obj["altloc"] else ""
        ]

        return " ".join(parts)


    def get_atom_spec(obj):
        return [
            obj["chain_id"].strip(),
            int(obj["resseq"]),
            obj["icode"].strip() if obj["icode"] else "",
            obj["atom_name"].strip(),
            obj["altloc"].strip() if obj["altloc"] else ""
        ]


    #---------Van der Waals Clashes------------------------
    vdw_items = []

    try:
        vdw_outliers = data[-1]["geom"]["outliers"]["vdw"]    
        for vdw in vdw_outliers:      
            atom1 = parse_atom_string(vdw["atom1"])
            atom2 = parse_atom_string(vdw["atom2"])
            vdw_items.append({
                "type": "Van der Waals Clash",
                "label": (
                    f"Van der Waals Clash: {get_atom_spec_as_string(atom1)} - "
                    f"{get_atom_spec_as_string(atom2)}, "
                    f"value={round(vdw.get('value', 0), 3)}, "
                    f"ideal={round(vdw.get('ideal', 0), 3)}, "
                    f"sigma={round(vdw.get('sigma', 0), 3)}, "
                    f"z={round(vdw.get('z', 0), 1)}, "
                    f"clash_type={vdw.get('type', '')}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(atom1),
                "atom-2-spec": get_atom_spec(atom2),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Van der Waals Clashes",
            "items": vdw_items
        })

    except KeyError:
        # safely skip if 'vdw' outliers are missing
        print("No van der Waals clashes")


    #--------- Bonds -----------------------
    bonds_items = []

    try:
        bonds = data[-1]["geom"]["outliers"]["bond"]
        for bond in bonds:
            atom1 = parse_atom_string(bond["atom1"])
            atom2 = parse_atom_string(bond["atom2"])
            # base label with always-present fields
            label = (
                f"Bond Outlier: {get_atom_spec_as_string(atom1)} - "
                f"{get_atom_spec_as_string(atom2)}, "
                f"value={round(bond.get('value', 0), 3)}, "
                f"ideal={round(bond.get('ideal', 0), 3)}, "
                f"sigma={round(bond.get('sigma', 0), 3)}, "
                f"z={round(bond.get('z', 0), 1)}"
            )

            # append type/alpha only if type >= 2
            bond_type = bond.get("type", 0)
            if bond_type >= 2:
                label += f", type={bond_type}, alpha={round(bond.get('alpha', 0), 1)}"

            bonds_items.append({
                "type": "Bond Outlier",
                "label": label,
                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(atom1),
                "atom-2-spec": get_atom_spec(atom2),
                "action": ["sphere-refinement-action"]
            })


        output["sections"].append({
            "title": "Bonds outliers",
            "items": bonds_items
        })
    except KeyError:
        print("No bond outliers")


#---------Angle------------------------

    angle_items = []

    try: 
        angles = data[-1]["geom"]["outliers"]["angle"]
        for angle in angles:
            atom1 = parse_atom_string(angle["atom1"])
            atom2 = parse_atom_string(angle["atom2"])
            atom3 = parse_atom_string(angle["atom3"])
            angle_items.append({
                "type": "Angle Outlier",
                "label": (
                    f"Angle Outlier: {get_atom_spec_as_string(atom1)} - "
                    f"{get_atom_spec_as_string(atom2)} - "
                    f"{get_atom_spec_as_string(atom3)}, "
                    f"value={round(angle.get('value', 0), 3)}, "
                    f"ideal={round(angle.get('ideal', 0), 3)}, "
                    f"sigma={round(angle.get('sigma', 0), 3)}, "
                    f"z={round(angle.get('z', 0), 1)}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(atom1),
                "atom-2-spec": get_atom_spec(atom2),
                "atom-3-spec": get_atom_spec(atom3),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Angles outliers",
            "items": angle_items
        })

    except KeyError:
        # safely skip if 'angle' outliers are missing
        print("No angle outliers")


#----------Torsion------------------------

    torsion_items = []

    try:
        torsions = data[-1]["geom"]["outliers"]["torsion"]
        for torsion in torsions:
            atom1 = parse_atom_string(torsion["atom1"])
            atom2 = parse_atom_string(torsion["atom2"])
            atom3 = parse_atom_string(torsion["atom3"])
            atom4 = parse_atom_string(torsion["atom4"])
            torsion_items.append({
                "type": "Torsion Angle Outlier",
                "label": (
                    f"Torsion Angle Outlier: {torsion.get('label', '')}, "
                    f"{get_atom_spec_as_string(atom1)} - "
                    f"{get_atom_spec_as_string(atom2)} - "
                    f"{get_atom_spec_as_string(atom3)} - "
                    f"{get_atom_spec_as_string(atom4)}, "
                    f"value={round(torsion.get('value', 0), 3)}, "
                    f"ideal={round(torsion.get('ideal', 0), 3)}, "
                    f"sigma={round(torsion.get('sigma', 0), 3)}, "
                    f"periodicity={round(torsion.get('per', 0), 3)}, "
                    f"z={round(torsion.get('z', 0), 1)}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(atom1),
                "atom-2-spec": get_atom_spec(atom2),
                "atom-3-spec": get_atom_spec(atom3),
                "atom-4-spec": get_atom_spec(atom4),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Torsion Angle outliers",
            "items": torsion_items
        })

    except KeyError:
        # safely skip if 'torsion' outliers are missing
        print("No torsion angle outliers")


    #---------Chirality------------------------
    chiral_items = []

    try:
        chirs = data[-1]["geom"]["outliers"]["chir"]
        for chir in chirs:
            atomc = parse_atom_string(chir["atomc"])
            atom1 = parse_atom_string(chir["atom1"])
            atom2 = parse_atom_string(chir["atom2"])
            atom3 = parse_atom_string(chir["atom3"])
            chiral_items.append({
                "type": "Chiral Outlier",
                "label": (
                    f"Chiral Outlier: {get_atom_spec_as_string(atomc)} - "
                    f"{get_atom_spec_as_string(atom1)} - "
                    f"{get_atom_spec_as_string(atom2)} - "
                    f"{get_atom_spec_as_string(atom3)}, "
                    f"value={round(chir.get('value', 0), 3)}, "
                    f"ideal={round(chir.get('ideal', 0), 3)}, "
                    f"sigma={round(chir.get('sigma', 0), 3)}, "
                    f"both={chir.get('both')}, "
                    f"z={round(chir.get('z', 0), 1)}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(atomc),
                "atom-2-spec": get_atom_spec(atom1),
                "atom-3-spec": get_atom_spec(atom2),
                "atom-4-spec": get_atom_spec(atom3),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Chiral Outliers",
            "items": chiral_items
        })

    except KeyError:
        # safely skip if 'chir' outliers are missing
        print("No chiral outliers")


    #---------Plane------------------------
    plane_items = []

    try:
        planes = data[-1]["geom"]["outliers"]["plane"]
        for plane in planes:
            atom1 = parse_atom_string(plane["atom"])
            plane_items.append({
                "type": "Plane Outlier",
                "label": (
                    f"Plane Outlier: {plane.get('label', '')}, "
                    f"{get_atom_spec_as_string(atom1)}, "
                    f"deviation={round(plane.get('dev', 0), 3)}, "
                    f"sigma={round(plane.get('sigma', 0), 3)}, "
                    f"z={round(plane.get('z', 0), 1)}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(atom1),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Plane Outliers",
            "items": plane_items
        })

    except KeyError:
        # safely skip if 'plane' outliers are missing
        print("No plane outliers")


    #---------Staca------------------------
    staca_items = []

    try:
        stacas = data[-1]["geom"]["outliers"]["staca"]
        for staca in stacas:
            plane1 = parse_atom_string(staca["plane1"])
            plane2 = parse_atom_string(staca["plane2"])
            staca_items.append({
                "type": "Staca Outlier",
                "label": (
                    f"Staca Outlier: {get_atom_spec_as_string(plane1)} - "
                    f"{get_atom_spec_as_string(plane2)}, "
                    f"value={round(staca.get('value', 0), 3)}, "
                    f"ideal={round(staca.get('ideal', 0), 3)}, "
                    f"sigma={round(staca.get('sigma', 0), 3)}, "
                    f"z={round(staca.get('z', 0), 1)}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(plane1),
                "atom-2-spec": get_atom_spec(plane2),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Staca Outliers",
            "items": staca_items
        })

    except KeyError:
        # safely skip if 'plane' outliers are missing
        print("No staca outliers")


    #---------Stacd------------------------
    stacd_items = []

    try:
        stacds = data[-1]["geom"]["outliers"]["stacd"]
        for stacd in stacds:
            plane1 = parse_atom_string(stacd["plane1"])
            plane2 = parse_atom_string(stacd["plane2"])
            stacd_items.append({
                "type": "Stacd Outlier",
                "label": (
                    f"Stacd Outlier: {get_atom_spec_as_string(plane1)} - "
                    f"{get_atom_spec_as_string(plane2)}, "
                    f"value={round(stacd.get('value', 0), 3)}, "
                    f"ideal={round(stacd.get('ideal', 0), 3)}, "
                    f"sigma={round(stacd.get('sigma', 0), 3)}, "
                    f"z={round(stacd.get('z', 0), 1)}"
                ),

                "position-type": "by-atom-spec-pair",
                "atom-1-spec": get_atom_spec(plane1),
                "atom-2-spec": get_atom_spec(plane2),
                "action": ["sphere-refinement-action"]
            })

        output["sections"].append({
            "title": "Stacd Outliers",
            "items": stacd_items
        })

    except KeyError:
        # safely skip if 'plane' outliers are missing
        print("No stacd outliers")


    with open(output_json_file_name, "w") as f:
        json.dump(output, f, indent=4)


# json_file_name = "3LNN_mod_refined5_stats.json"

def convert_outliers_json(json_file_name):
    ff,ext = os.path.splitext(json_file_name)
    output_json_file_name = ff + "_for_coot.json"
    convert_outliers_json_inner(json_file_name, output_json_file_name)
    import coot
    coot.read_interesting_places_json_file(output_json_file_name)



