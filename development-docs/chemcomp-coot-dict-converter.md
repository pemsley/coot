

On the converter: bidirectional, but the two directions aren't symmetric

Yes, build it bidirectional — but for different reasons and to different completeness bars:

Forward (Coot dictionary_residue_restraints_t → gemmi ChemComp): required now, and it's the hard/correctness-critical direction. classify_atoms keys almost entirely on bond order + aromaticity + formal charge, so those must be faithful or the AceDRG lookup silently picks wrong rows. Must carry: atoms (id, element, formal charge, type_energy↔chem_type), bonds (id1/id2, bond type, aromaticity flag), angles (id1/id2/id3), and include hydrogens (H-connectivity affects typing).

Reverse (gemmi ChemComp → Coot dict): split into two things.
- Value merge — copy filled bond/angle value+esd back into the existing Coot dict. Needed now. Do it name-keyed (atom-name tuples, order-insensitive), not by array index the way the existing patch does — gemmi/classify may reorder or drop entries, and name-keying matches how Coot's own conservatively_replace_with already works.
- Full rebuild — construct a whole Coot dict from an arbitrary ChemComp (torsions, chirals, planes, coords, metadata). Broadly useful as Coot goes gemmi-native, and it's what makes a round-trip test (dict → ChemComp → dict) possible — which is the best way to validate forward fidelity. But it's more surface than the AceDRG task needs, so build it as the reusable component, not a blocker.

So: full forward + name-keyed value-merge reverse gets AceDRG working; the full reverse is worth doing for the round-trip test and Coot's broader gemmi direction, just decoupled from this task.

Field-mapping gotchas the converter must nail (these are where a naive converter goes subtly wrong):
- Bond type/aromaticity: Coot has type() (string) and a separate aromaticity enum; gemmi has BondType enum + aromatic bool. Aromatic can be encoded in either — map both consistently.
- Formal vs partial charge: classify wants formal charge for hybridization. The existing patch feeds partial_charge into gemmi's charge — likely wrong; keep them distinct.
- Nuclear vs electron distances: gemmi Bond carries value/esd and value_nucleus/esd_nucleus; Coot has nuclear fields + a use_nuclear flag. Decide which the merge writes (pyrogen sets set_use_nuclear_distances(true)).
- Atom-name form: Coot atom_id vs atom_id_4c (padded) vs gemmi id — normalize so name-keyed matching is reliable.


