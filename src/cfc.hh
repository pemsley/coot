
// header here

#ifndef CFC_HH
#define CFC_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#ifdef USE_PYTHON
#include "Python.h"

#include <map>
#include <algorithm>

#include <clipper/core/coords.h>
#include "coot-utils/residue-and-atom-specs.hh"

// return a Python object that contains (with indexing)
//   water positions around ligand
//   chemical features of ligands
//
// environment_residues_py is a list of residue specs
// solvated_ligand_info_py is a list of
//   list mol_no ligand_spec
// which, with radius will be used to find the waters
// 
PyObject *chemical_feature_clusters_py(PyObject *environment_residues_py,
				       PyObject *solvated_ligand_info_py,
				       double radius_1, double radius_2);
// scipy has done some clustering
// we get the cluster info here
// 
void chemical_feature_clusters_accept_info_py(PyObject *env_residue_py,
					      PyObject *mol_ligand_specs_py,
					      PyObject *cluster_info_py);

namespace cfc {

   // the water cluster spec
   //
   class water_cluster_info_from_python {
   public:
      water_cluster_info_from_python(const clipper::Coord_orth &pos_in,
				     double w_in,
				     double r_in) {
	 pos = pos_in;
	 weight = w_in;
	 radius = r_in;
      }
      water_cluster_info_from_python() {weight = 0;}
      clipper::Coord_orth pos;
      double weight;
      double radius;
   };

   // the container for the water residues, annotated by which cluster the water is in
   // 
   class clustered_feature_info_from_python {
   public:
      clustered_feature_info_from_python(int imol_in, coot::residue_spec_t &spec_in, unsigned int cluster_number_in) {
	 imol = imol_in;
	 residue_spec = spec_in;
	 cluster_number = cluster_number_in;
      }
      clustered_feature_info_from_python() { imol = -1; cluster_number = -1;}
      int imol;
      unsigned int cluster_number; // unsigned because it indexes into the cluster means vector
      coot::residue_spec_t residue_spec;
   };

   class extracted_cluster_info_from_python {
      
      void extract_water_info(PyObject *cluster_info_py);
      void extract_chemical_feature_info(PyObject *cf_py);

      std::vector<clipper::Coord_orth> extract_cluster_means(PyObject *means_py);

   public:
      std::vector<water_cluster_info_from_python>   wc;
      std::vector<clustered_feature_info_from_python> cw;
      extracted_cluster_info_from_python(PyObject *cluster_info_py);
      unsigned int n_water_structures() const;
      unsigned int n_pharmacophore_structures() const; 
      std::vector<int> water_structures_vec() const;  // for structure buttons (water cluster)
      std::vector<int> pharmacophore_structures_vec() const; // for structure buttons (pharma)
      std::vector<std::pair<int, coot::residue_spec_t> >
      pharmacophore_structures_and_specs_vec() const; // for structure buttons (pharma)
      
      std::vector<std::pair<int, coot::residue_spec_t> > water_cluster_imol_residue_spec_vec() const;
      std::vector<std::pair<int, coot::residue_spec_t> > pharmacophore_cluster_imol_residue_spec_vec(const std::string &type, unsigned int cluster_idx);
      unsigned int water_cluster_idx_max() const;

      std::map<std::string, std::vector<clustered_feature_info_from_python> > pharmacophore;
      std::map<std::string, std::vector<clipper::Coord_orth> > pharmacophore_model_cluster_means;
      
      // sort by the number of molecules (structures) are present in this cluster
      static bool cluster_vector_sorter(const std::pair<std::vector<int>, water_cluster_info_from_python> &v1,
					const std::pair<std::vector<int>, water_cluster_info_from_python> &v2) {
	 return (v2.first.size() < v1.first.size());
      }
      // return the generic display object index
      int show_water_balls() const;
   };

}
   



#endif // USE_PYTHON

#endif // MAKE_ENHANCED_LIGAND_TOOLS

#endif // CFC_HH
