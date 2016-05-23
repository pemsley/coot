
// header-here

#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coot-utils/coot-coord-utils.hh"
#include "chemical-feature-clusters.hh"

int main(int argc, char **argv) {

   typedef std::pair<std::string, coot::residue_spec_t> file_lig;

   std::vector<file_lig> l;
   l.push_back(file_lig("../src/coot-download/pdb4zzn.ent", coot::residue_spec_t("A", 1355, "")));
   l.push_back(file_lig("../src/coot-download/pdb4zzm.ent", coot::residue_spec_t("A", 1357, "")));
   l.push_back(file_lig("../src/coot-download/pdb4zzo.ent", coot::residue_spec_t("A", 1357, "")));

   std::vector<coot::chem_feat_solvated_ligand_spec> ligands;
   coot::protein_geometry pg;
   std::vector<coot::residue_spec_t> neighbs_residues;

   for (unsigned int il=0; il<l.size(); il++) {

      const coot::residue_spec_t &ligand_spec = l[il].second;

      atom_selection_container_t atom_sel =
	 get_atom_selection(l[il].first, true, true);

      if (atom_sel.read_success) {

	 float radius = 6;
	 std::vector<coot::residue_spec_t> neighbs_raw =
	    coot::residues_near_residue(l[il].second, atom_sel.mol, radius);

	 std::vector<coot::residue_spec_t> neighbs_waters;

	 for (unsigned int i=0; i<neighbs_raw.size(); i++) { 
	    mmdb::Residue *res = coot::util::get_residue(neighbs_raw[i],
							 atom_sel.mol);
	    if (res) {
	       std::string res_name = res->GetResName();
	       if (res_name == "HOH")
		  neighbs_waters.push_back(neighbs_raw[i]);
	       else
		  if (il == 0)
		     neighbs_residues.push_back(neighbs_raw[i]);
	    }
	 }
	 
	 std::cout << "INFO:: found " << neighbs_residues.size()
		   << " residue neighbour of the ligand and "
		   << neighbs_waters.size() << " water neighbours for ligand "
		   << ligand_spec << std::endl;

	 coot::chem_feat_solvated_ligand_spec lig(ligand_spec, neighbs_waters,
						  atom_sel.mol);
	 ligands.push_back(lig);

      }
   }

   coot::chem_feat_clust cl(neighbs_residues, ligands, &pg);
   cl.cluster_waters();

   return 0;
}
