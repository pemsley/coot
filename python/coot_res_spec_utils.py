

# rename coot_utils.py functions

def residue_spec_to_chain_id(rs):
    if not isinstance(rs, list):
        return False
    else:
        if (len(rs) == 3):
            return rs[0]
        else:
            if (len(rs) == 4):
                return rs[1]
            else:
                return False

def residue_spec_to_res_no(rs):
    if not isinstance(rs, list):
        return False
    else:
        if (len(rs) == 3):
            return rs[1]
        else:
            if (len(rs) == 4):
                return rs[2]
        return False

def residue_spec_to_ins_code(rs):
    if not isinstance(rs, list):
        return False
    else:
        if (len(rs) == 3):
            return rs[2]
        else:
            if (len(rs) == 4):
                return rs[3]
        return False


def residue_atom_to_position(ra):
    if not isinstance(ra, list):
        return False
    else:
        return ra[2]

