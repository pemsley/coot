import sys
import os

def parse_screenfile(file_name):
    r = {}
    if os.path.exists(file_name):
        with open(file_name) as f:
            lines = f.readlines()
            front_page           = []
            bond_length_outliers = []
            bond_angle_outliers  = []
            torsion_outliers     = []
            plane_outliers       = []
            unhappy_atoms        = []
            ideal_contact_outliers = []
            aniso_sph_outliers     = []
            aniso_sph_nonb_outliers = []
            front_page_is_live             = False
            bond_length_outliers_is_live   = False
            bond_angle_outliers_is_live    = False
            torsion_outliers_is_live       = False
            plane_outliers_is_live         = False
            ideal_contact_outliers_is_live = False
            unhappy_atom_is_live           = False
            aniso_sph_is_live              = False
            aniso_sph_nonb_is_live          = False
            for line in lines:

                # ---------- Front page ----------------
                if front_page_is_live:
                    if " Bond lenth outliers," in line:
                        front_page_is_live = False
                    if " Weighted rms " in line:
                        front_page_is_live = False
                if front_page_is_live:
                    front_page.append(line.strip())
                if " Screen output " in line:
                    front_page_is_live = True
                    front_page.append(line.strip())

                # ---------- Bonds ----------------
                if bond_length_outliers_is_live:
                    if line == " \n":
                        bond_length_outliers_is_live = False
                if bond_length_outliers_is_live:
                    if "Actual       Ideal       Delta       Sigma" in line:
                        pass
                    else:
                        parts = line.split()
                        if len(parts) == 7:
                            actual      = parts[0]
                            ideal       = parts[1]
                            delta       = parts[2]
                            sigma       = parts[3]
                            delta_sigma = parts[4]
                            atom_ids    = parts[5]
                            res_type    = parts[6]
                            blo = {"actual": actual, "ideal": ideal, "sigma": sigma, "delta_sigma": delta_sigma,
                                   "atom_ids": atom_ids, "res_type": res_type}
                            bond_length_outliers.append(blo)
                        if len(parts) == 0:
                            bond_length_outliers_is_live = False
                if " Bond length outliers, distance" in line:
                    bond_length_outliers_is_live = True

                # ---------- Angles ----------------
                if bond_angle_outliers_is_live:
                    if line == " \n":
                        bond_angle_outliers_is_live = False
                if bond_angle_outliers_is_live:
                    if "Torsion angles in degrees." in line:
                        bond_angle_outliers_is_live = False

                if bond_angle_outliers_is_live:
                    parts = line.split()
                    if len(parts) == 0:
                        bond_angle_outliers_is_live = False
                    else:
                        if len(parts) == 7:
                            actual = parts[0]
                            ideal  = parts[1]
                            delta  = parts[2]
                            sigma  = parts[3]
                            delta_sigma = parts[4]
                            atom_ids = parts[5]
                            res_type = parts[6]
                            bao = {"actual": actual, "ideal": ideal, "sigma": sigma, "delta_sigma": delta_sigma,
                                   "atom_ids": atom_ids, "res_type": res_type}
                            bond_angle_outliers.append(bao)

                if "Bond angle outliers, angles in degrees" in line:
                    bond_angle_outliers_is_live = True

                # ---------- Torsions ----------------
                if torsion_outliers_is_live:
                    if line == " \n":
                        torsion_outliers_is_live = False
                if torsion_outliers_is_live:
                    if "PLAN plane outliers," in line:
                        torsion_outliers_is_live = False
                if torsion_outliers_is_live:
                    if "Ideal-contact dist" in line:
                        torsion_outliers_is_live = False
                if torsion_outliers_is_live:
                    parts = line.split()
                    if len(parts) == 0:
                        torsion_outliers_is_live = False
                    else:
                        if len(parts) == 7:
                            actual = parts[0]
                            ideal  = parts[1]
                            delta  = parts[2]
                            sigma  = parts[3]
                            delta_sigma = parts[4]
                            atom_ids    = parts[5]
                            res_type    =  parts[6]
                            to = {"actual": actual, "ideal": ideal, "sigma": sigma, "delta_sigma": delta_sigma,
                                   "atom_ids": atom_ids, "res_type": res_type}
                            torsion_outliers.append(to)
                if "Torsion angles in degrees." in line:
                    torsion_outliers_is_live = True;

                # ---------- Planes ----------------
                if plane_outliers_is_live:
                    if line == " \n":
                        plane_outliers_is_live = False
                if plane_outliers_is_live:
                    if "Ideal-contact dist" in line:
                        plane_outliers_is_live = False
                if plane_outliers_is_live:
                    parts = line.split()
                    if len(parts) == 0:
                        plane_outliers_is_live = False
                    else:
                        if len(parts) >= 4:
                            actual       = parts[0]
                            sigma        = parts[1]
                            actual_sigma = parts[2]
                            atom_ids     = parts[3]
                            po = {"actual": actual, "sigma": sigma, "actual_sigma": actual_sigma,
                                   "atom_ids": atom_ids}
                            if len(parts) == 5:
                                res_type = parts[4]
                                po["res_type"] = res_type
                            plane_outliers.append(po)
                if " PLAN plane outliers, r.m.s. devia from the plane in Angs" in line:
                    plane_outliers_is_live = True;

                # ---------- Ideal_Contacts ----------------
                if ideal_contact_outliers_is_live:
                    if line == " \n":
                        ideal_contact_outliers_is_live = False
                if ideal_contact_outliers_is_live:
                    if "Ideal-contact distance outliers, total number" in line:
                        ideal_contact_outliers_is_live = False
                if ideal_contact_outliers_is_live:
                    parts = line.split()
                    if len(parts) == 0:
                        ideal_contact_outliers_is_live = False
                    else:
                        if len(parts) > 6:
                            actual      = parts[0]
                            contactD    = parts[1]
                            delta       = parts[2]
                            sigma       = parts[3]
                            delta_sigma = parts[4]
                            atom_ids    = line[60:].strip()
                            co = {"actual": actual, "contactD": contactD, "delta": delta,
                                  "sigma": sigma, "delta_sigma": delta_sigma, "atom_ids": atom_ids}
                            ideal_contact_outliers.append(co)
                if "Ideal-contact distance outliers, distances in " in line:
                    ideal_contact_outliers_is_live = True;

                # ---------- Unhappy Atoms ----------------
                if unhappy_atom_is_live:
                    if line == " \n":
                        unhappy_atom_is_live = False
                if unhappy_atom_is_live:
                    if "maximum number of reports in this screen output" in line:
                        unhappy_atom_is_live = False
                if unhappy_atom_is_live:
                    parts = line.split()
                    if len(parts) == 2:
                        f_over_median = parts[0]
                        atom_id = parts[1]
                        ua = {"f_over_median": f_over_median, "atom_id": atom_id}
                        unhappy_atoms.append(ua)
                if "    Report the (function contribution)/(median contribution). Median=" in line:
                    unhappy_atom_is_live = True

                # ------------- anisotropic atoms ------------------------

                if aniso_sph_is_live:
                    if line == " \n":
                        aniso_sph_is_live = False

                if "Atom" in line: continue

                if aniso_sph_is_live:
                    parts = line.split()
                    if len(parts) == 3:
                        parts = line.split()
                        value    = parts[0]
                        function = parts[1]
                        atom     = parts[2]
                        asph = {"value": value, "function": function, "atom": atom}
                        aniso_sph_outliers.append(asph)

                if "Anisotropic restraint outliers for SPH " in line:
                    aniso_sph_is_live = True

                # ------------- anisotropic atoms non-bonded ------------------------

                if aniso_sph_nonb_is_live:
                    if line == " \n":
                        aniso_sph_nonb_is_live = False

                if aniso_sph_nonb_is_live:
                    parts = line.split()
                    if len(parts) == 3:
                        parts = line.split()
                        value    = parts[0]
                        function = parts[1]
                        atom     = parts[2]
                        asph = {"value": value, "function": function, "atom": atom}
                        aniso_sph_nonb_outliers.append(asph)

                if "Anisotropic restraint outliers for SPH[nonb] " in line:
                    aniso_sph_nonb_is_live = True

            # consolidate:
            if front_page:              r["Front page"]     = front_page
            if bond_length_outliers:    r["Bond length"]    = bond_length_outliers
            if bond_angle_outliers:     r["Bond angle"]     = bond_angle_outliers
            if torsion_outliers:        r["Torsion"]        = torsion_outliers
            if plane_outliers:          r["Plane"]          = plane_outliers
            if ideal_contact_outliers:  r["Ideal Contact"]  = ideal_contact_outliers
            if unhappy_atoms:           r["Unhappy Atom"]   = unhappy_atoms
            if aniso_sph_outliers:      r["Aniso SPH"]      = aniso_sph_outliers
            if aniso_sph_nonb_outliers: r["Aniso SPH nonb"] = aniso_sph_nonb_outliers

    else:
        print("file does not exist", file_name)
    return r



if __name__ == "__main__":

    if "coot-1" in sys.argv[0]:

        pass
        # file_name = "../../../../../../buster-tutorial/MapOnly/01-BUSTER/Cycle-2/shell.01/screen_final.txt"
        # r = parse_screenfile(file_name)
        # coot.global_phasing_screen(0, r)

    else:

        file_name = sys.argv[1]
        r = parse_screenfile(file_name)
        for key in r:
            print(key, ":::")
            data = r[key]
            for item in data:
                print(item)
