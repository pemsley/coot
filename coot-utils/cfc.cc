

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <fstream>
#include <memory>

#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"
#include "coot-coord-utils.hh"
#include "dirichlet-process.hh"
#include "gmm.hh"
#include "json.hpp"
#include "cfc.hh"
#include "lidia-core/rdkit-interface.hh"

#include "utils/logging.hh"
extern logging logger;


std::pair<std::vector<cfc::typed_cluster_t>, std::vector<std::vector<cfc::water_info_t> > >
cfc::chemical_feature_clustering(const std::vector<cfc::input_info_t> &mol_infos,
                                 const coot::protein_geometry &geom,
                                 const cfc::clustering_options_t &options) {

   // the ligand is the residue called res_name. If there is more than one of them,
   // use the one closest to pt (if use_pt), otherwise the first one.
   auto find_ligand_residue = [] (mmdb::Manager *mol, const std::string &res_name,
                                  bool use_pt, const clipper::Coord_orth &pt) {

      double best_dd = 1.0e30;
      mmdb::Residue *r = nullptr;
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string rn = residue_p->GetResName();
                     if (rn == res_name) {
                        if (! use_pt) return residue_p;
                        auto rc = coot::util::get_residue_centre(residue_p);
                        if (rc.first) {
                           double dd = (rc.second - pt).lengthsq();
                           if (dd < best_dd) {
                              best_dd = dd;
                              r = residue_p;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      return r;
   };

   // main-chain LSQ matches for the protein residues of mol_ref near pt - the same
   // residue numbers are expected in the moving molecules
   auto make_lsq_matchers = [] (mmdb::Manager *mol_ref, const clipper::Coord_orth &pt, double radius) {
      std::vector<coot::lsq_range_match_info_t> lsq_matchers;
      std::vector<mmdb::Residue *> residues = coot::residues_near_position(pt, mol_ref, radius);
      for (mmdb::Residue *residue_p : residues) {
         std::string rn = residue_p->GetResName();
         if (coot::util::is_standard_amino_acid_name(rn)) {
            std::string chain_id = residue_p->GetChainID();
            int res_no = residue_p->GetSeqNum();
            coot::lsq_range_match_info_t m(res_no, res_no, chain_id, res_no, res_no, chain_id, coot::lsq_t::MAIN);
            lsq_matchers.push_back(m);
         }
      }
      return lsq_matchers;
   };

   auto transform_point = [] (const clipper::RTop_orth &rtop, const RDGeom::Point3D &p) {
      clipper::Coord_orth pt(p.x, p.y, p.z);
      clipper::Coord_orth t = rtop * pt;
      return RDGeom::Point3D(t.x(), t.y(), t.z());
   };

   auto file_to_string = [] (const std::string &file_name) {

      std::string s;
      std::string line;
      std::ifstream f(file_name.c_str());
      if (!f) {
         std::cout << "WARNING:: file_to_string(): Failed to open " << file_name << std::endl;
      } else {
         while (std::getline(f, line)) {
            s += line;
            s += "\n";
         }
      }
      return s;
   };

   class feature_info_t {
   public:
      std::string family;
      std::string type;
      int imol;
      coot::residue_spec_t residue_spec;
      std::string residue_name;
      RDGeom::Point3D pos;
      feature_info_t(const std::string &f, const std::string &t, int imol_in,
                     const coot::residue_spec_t &rs, const RDGeom::Point3D &p) :
         family(f), type(t), imol(imol_in), residue_spec(rs), pos(p) {}
   };

   auto waters_near_ligand_site = [] (int imol, mmdb::Manager *mol, const clipper::Coord_orth &pt_ref) {
      double dist_crit = 6.0;
      double dist_crit_sq = dist_crit * dist_crit;
      std::vector<water_info_t> waters;
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string rn = residue_p->GetResName();
                     if (rn == "HOH") {
                        mmdb::Atom **residue_atoms = 0;
                        int n_residue_atoms = 0;
                        residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                        for (int iat=0; iat<n_residue_atoms; iat++) {
                           mmdb::Atom *at = residue_atoms[iat];
                           if (! at->isTer()) {
                              clipper::Coord_orth posc(at->x, at->y, at->z);
                              double dd = (posc - pt_ref).lengthsq();
                              if (dd < dist_crit_sq) {
                                 RDGeom::Point3D pos(at->x, at->y, at->z);
                                 coot::residue_spec_t res_spec(residue_p);
                                 water_info_t wi(imol, res_spec, pos);
                                 waters.push_back(wi);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      return waters;
   };

   // Look for BaseFeatures.fdef in the current directory, then the file named by
   // RDKIT_BASEFEATURES, then $prefix/share/RDKit/Data/BaseFeatures.fdef.
   // Return the empty string if it was not found.
   auto find_feature_definitions_file = [] () {
      std::string fn = "BaseFeatures.fdef";
      if (coot::file_exists(fn)) return fn;
      const char *e = getenv("RDKIT_BASEFEATURES");
      if (e) {
         if (coot::file_exists(e)) return std::string(e);
         logger.log(log_t::WARNING, {"RDKIT_BASEFEATURES file", std::string(e), "not found"});
      }
      std::string data_dir = coot::util::append_dir_dir(coot::rdkit_package_data_dir(), "Data");
      std::string p = coot::util::append_dir_file(data_dir, fn);
      if (coot::file_exists(p)) return p;
      return std::string();
   };

   auto chemical_features = [] (mmdb::Residue *residue_p,
                                int imol,
                                int imol_enc,
                                const coot::protein_geometry &geom,
                                RDKit::MolChemicalFeatureFactory *factory) {

      std::vector<feature_info_t> feature_infos; // for this residue
      try {

         std::cout << "   Chemical features for " << coot::residue_spec_t(residue_p) << std::endl;

         RDKit::RWMol m = coot::rdkit_mol_sanitized(residue_p, imol_enc, geom);
         int n_atoms = m.getNumHeavyAtoms();
         // std::cout << "   rdkit mol m has " << n_atoms << " (non-hydrogen) atoms" << std::endl;
         if (n_atoms > 0) {

            int n_conf = m.getNumConformers();
            if (n_conf > 0) {
               int iconf = 0;
               RDKit::Conformer &conf = m.getConformer(iconf);
               // now interrogate the conformer
               if (false) {
                  std::cout << "   Conformer " << iconf << " has " << conf.getNumAtoms()
                            << " atoms" << std::endl;
                  for (unsigned int i=0; i<conf.getNumAtoms(); i++) {
                     RDKit::Atom *at = m.getAtomWithIdx(i);
                     RDGeom::Point3D pos = conf.getAtomPos(i);
                     std::cout << "   atom " << i << " " << at->getIdx() << " "
                               << at->getSymbol() << " " << at->getAtomicNum() << " "
                               << at->getFormalCharge() << " " << at->getIsotope() << " "
                               << at->getHybridization() << " " << at->getDegree() << " "
                               << at->getTotalNumHs() << " " << at->getIsAromatic() << " "
                               << "pos " << pos << std::endl;
                  }
               }
            }

            RDKit::FeatSPtrList features = factory->getFeaturesForMol(m);
            std::cout << "DEBUG:: Found " << features.size() << " features" << std::endl;
            if (features.empty()) {
               std::string residue_name = residue_p->GetResName();
               // std::cout << "INFO:: -- boo! no features found for a " << residue_name
               //           << " (surprising)" << std::endl;
               logger.log(log_t::INFO, "-- boo! no features found for a", residue_name, "(surprising)");
            } else {
               coot::residue_spec_t res_spec(residue_p);
               for (const auto &feature : features) {
                  std::cout << "Feature: " << feature->getFamily() << " " << feature->getType()
                            << " at pos " << feature->getPos() << std::endl;
                  feature_info_t fi(feature->getFamily(), feature->getType(), imol, res_spec, feature->getPos());
                  fi.residue_name = residue_p->GetResName();
                  feature_infos.push_back(fi);
               }
            }

         } else {
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (! at->isTer()) {
                  std::cout << "   " << iat << " " << coot::atom_spec_t(at) << std::endl;
               }
            }
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
      return feature_infos;
   };

   auto make_test_points = [] (const std::vector<feature_info_t> &v) {
      std::vector<glm::vec3> points;
      for (const auto &fi : v) {
         glm::vec3 p(fi.pos.x, fi.pos.y, fi.pos.z);
         points.push_back(p);
      }
      return points;
   };

   auto get_n_clusters = [] (std::vector<unsigned int> &clustered_points) {
      std::set<int> cluster_indices;
      for (const auto &c : clustered_points) {
         cluster_indices.insert(c);
      }
      return cluster_indices.size();
   };

   auto split_string = [] (const std::string &s, const std::string &delim) {
      std::pair<std::string, std::string> p;
      size_t pos = s.find(delim);
      if (pos != std::string::npos) {
         p.first = s.substr(0, pos);
         p.second = s.substr(pos + delim.length());
      }
      return p;
   };

   // manually set the num_clusters
   auto cluster_waters_gmm_old = [] (const std::vector<water_info_t> &water_infos) {

      std::map<int, water_info_t> iwat_to_water_info_index;
      std::vector<glm::vec3> water_positions;
      for (const auto &water : water_infos) {
         glm::vec3 p(water.pos.x, water.pos.y, water.pos.z);
         int iwat = water_positions.size();
         water_positions.push_back(p);
         iwat_to_water_info_index[iwat] = water;
      }

      unsigned int num_clusters = 3;
      std::vector<glm::vec3> data;
      for (const auto &water : water_infos) {
         glm::vec3 pos(water.pos.x, water.pos.y, water.pos.z);
         data.push_back(pos);
      }
      // Create and fit the GMM.
      GMM gmm(num_clusters);
      gmm.fit(data, 200);
      gmm.printParameters("Waters Final");

      std::vector<std::vector<water_info_t> > clusters;
      return clusters;

   };

   auto cluster_waters_gmm = [] (const std::vector<water_info_t> &water_infos) {

      std::map<int, water_info_t> iwat_to_water_info_index;
      std::vector<glm::vec3> water_positions;
      for (const auto &water : water_infos) {
         glm::vec3 p(water.pos.x, water.pos.y, water.pos.z);
         int iwat = water_positions.size();
         water_positions.push_back(p);
         iwat_to_water_info_index[iwat] = water;
      }

      // Use BIC to automatically select optimal number of clusters
      // Test k from 1 to min(30, number of waters)
      int k_min = 1;
      int k_max = std::min(30, static_cast<int>(water_infos.size()));
      int max_iterations = 200;

      std::cout << "Testing GMM with BIC for optimal number of clusters..." << std::endl;
      GMMFitResult gmm_result = fitGMMWithBIC(water_positions, k_min, k_max, max_iterations, true);

      std::cout << "Selected " << gmm_result.optimal_k << " clusters for waters (BIC="
                << gmm_result.best_bic << ")" << std::endl;
      if (gmm_result.optimal_k == k_max && k_max < static_cast<int>(water_infos.size()))
         logger.log(log_t::WARNING, logging::function_name_t("cluster_waters_gmm"),
                    {"the BIC chose k = k_max =", k_max, "- the optimal number of water clusters may be larger"});

      // Convert cluster assignments to the expected output format
      std::vector<std::vector<water_info_t>> clusters(gmm_result.optimal_k);

      for (size_t iwat = 0; iwat < water_positions.size(); ++iwat) {
         int cluster_index = gmm_result.cluster_assignments[iwat];
         auto wi = water_infos[iwat];
         clusters[cluster_index].push_back(wi);
      }

      // Print cluster statistics
      if (true) {
         std::cout << "DEBUG:: cluster_waters_gmm(): Water cluster statistics:" << std::endl;
         for (unsigned int ic = 0; ic < clusters.size(); ic++) {
            std::cout << "  Cluster " << ic << " has " << clusters[ic].size()
                      << " waters" << std::endl;
         }
      }

      return clusters;
   };


   auto cluster_waters_dirichlet = [get_n_clusters] (const std::vector<water_info_t> &water_infos) {

      std::map<int, water_info_t> iwat_to_water_info_index;
      std::vector<glm::vec3> water_positions;
      for (const auto &water : water_infos) {
         glm::vec3 p(water.pos.x, water.pos.y, water.pos.z);
         int iwat = water_positions.size();
         water_positions.push_back(p);
         iwat_to_water_info_index[iwat] = water;
      }
      double alpha = 2.0; // was 9.1;
      double beta  = 0.0085; // was 0.01;
      DirichletProcessClustering dpc(alpha, beta);
      std::vector<unsigned int> clustered_points = dpc.fit(water_positions);
      if (clustered_points.size() == water_positions.size()) {
         unsigned int n_clusters = get_n_clusters(clustered_points);
         std::cout << ":::::::::: dirichletprocess n_clusters for waters " << n_clusters << std::endl;
         std::vector<std::vector<water_info_t> > clusters(n_clusters);
         for (unsigned int iwat=0; iwat<water_positions.size(); iwat++) {
            unsigned int cluster_index = clustered_points[iwat];
            std::cout << "debug:: iwat " << iwat << " cluster_index " << cluster_index
                      << std::endl;
            auto wi = water_infos[iwat];
            std::cout << "adding water info with imol " << wi.imol << " to cluster number "
                      << cluster_index << std::endl;
            clusters[cluster_index].push_back(wi);
         }
         if (true) {
            for (unsigned int ic=0; ic<clusters.size(); ic++) {
               std::cout << "post: water cluster " << ic << " has " << clusters[ic].size()
                         << " waters" << std::endl;
            }
         }

         if (false)
            std::cout << "------- in chemical_feature_clustering() cluster_waters() end --- "
                      << std::endl;

         for (unsigned int i=0; i<clusters.size(); i++) {
            const auto &wi = clusters[i];
            std::cout << "cluster " << i << " has " << wi.size() << " contributons"
                      << std::endl;
            std::set<int> imols_in_cluster;
            for (unsigned int jj=0; jj<wi.size(); jj++) {
               int imol = wi[jj].imol;
               imols_in_cluster.insert(imol);
               if (false)
                  std::cout << "cluster " << i <<  ": adding imol " << imol << std::endl;
            }
            if (false)
               std::cout << "debug:: water_clusters [" << i << "] has "
                         << imols_in_cluster.size() << " waters"
                         << std::endl;
         }

         if (false)
            std::cout << "------------" << std::endl;


         return clusters;
      } else {
         logger.log(log_t::INFO, logging::function_name_t(__FUNCTION__), " bad cluster length");
         std::vector<std::vector<water_info_t> > dummy;
         return dummy;
      }
   };

   auto cluster_features = [make_test_points, get_n_clusters, split_string]
      (const std::vector<feature_info_t> &feature_infos) {

      std::vector<typed_cluster_t> typed_clusters;
      if (true) {
         std::cout << "--------- features infos at start of cluster_features" << std::endl;
         for (const auto &fi : feature_infos) {
            std::cout << "debug:: fi " <<fi.family << " " << fi.type << " " << fi.imol
                      << " " << fi.residue_spec << " " << fi.pos << std::endl;
         }
      }

      std::map<std::string, std::vector<feature_info_t> > feature_info_map;

      for (const auto &fi : feature_infos) {
         std::string key = fi.family + "_" + fi.type;
         feature_info_map[key].push_back(fi);
      }

      if (true) {
         for (const auto &pair : feature_info_map) {
            std::string key = pair.first;
            std::cout << "debug:: for key: " << key << ": imols: ";
            const std::vector<feature_info_t> &v = pair.second;
            for (const auto &fi : v)
               std::cout << fi.imol << " ";
            std::cout << std::endl;
            for (const auto &fi : v)
               std::cout << "    fi pos " << fi.pos << std::endl;
         }
      }

      for (const auto &pair : feature_info_map) {
         std::string key = pair.first;
         const std::vector<feature_info_t> &vv = pair.second;
         if (true)
            std::cout << "DEBUG:: feature_info_map key " << key
                      << " has " << vv.size() << " features " << std::endl;
         double alpha = 1.0;
         double beta  = 0.01;
         std::vector<glm::vec3> v = make_test_points(vv);
         DirichletProcessClustering dpc(alpha, beta);
         std::vector<unsigned int> clustered_points = dpc.fit(v);
         std::cout << "debug:: dirichletprocess v in: " << v.size()
                   << " clustered_points out " << clustered_points.size() << std::endl;
         if (true) {
            std::cout << "   in cluster ";
            for (const auto &item : clustered_points)
               std::cout << item << " ";
            std::cout << std::endl;
         }
         unsigned int n_clusters = get_n_clusters(clustered_points);
         logger.log(log_t::DEBUG, logging::function_name_t("cluster_features"),
                    {"Cluster", key, "had", n_clusters, "clusters"});
         std::pair<std::string, std::string> p = split_string(key, "_");
         std::string family = p.first;
         std::string type   = p.second;
         for (unsigned int iclust=0; iclust<n_clusters; iclust++) {
            typed_cluster_t tc(family, type, iclust);
            RDGeom::Point3D pos_sum = RDGeom::Point3D(0,0,0);
            std::vector<RDGeom::Point3D> contributing_points;
            unsigned int n_contributors = 0;
            for (unsigned int i_feat_info=0; i_feat_info<vv.size(); i_feat_info++) {
               if (clustered_points[i_feat_info] == iclust) {
                  const feature_info_t &fi = vv[i_feat_info];
                  coot::residue_spec_t residue_spec = fi.residue_spec;
                  int imol_feat = fi.imol;
                  n_contributors++;
                  pos_sum += fi.pos;
                  contributing_points.push_back(fi.pos);
                  bool found = false;
                  for (unsigned int ii=0; ii<tc.imols_with_specs.size(); ii++) {
                     if (imol_feat == tc.imols_with_specs[ii].first) {
                        if (true) { // check the spec too?
                           found = true;
                        }
                     }
                  }
                  if (! found)
                     tc.add_imol(imol_feat, residue_spec);
               }
            }
            if (n_contributors > 0) {
               tc.pos = pos_sum / static_cast<double>(n_contributors);
               tc.contributing_points = contributing_points;
               std::cout << "debug:: n_contributors " << n_contributors
                         << " average pos " << pos_sum << std::endl;
            }
            typed_clusters.push_back(tc);
         }
      }

      if (true) {
         for (const auto &tc : typed_clusters) {
            std::cout << "DEBUG:: typed_cluster: " << tc.family << " " << tc.type
                      << " cluster-index: " << tc.idx << " imols:\n";
            for (const auto &imol_with_spec : tc.imols_with_specs) {
               int imol                          = imol_with_spec.first;
               coot::residue_spec_t residue_spec = imol_with_spec.second;
               std::cout << "      " << imol << " " << residue_spec << std::endl;
            }
            std::cout << " at " << tc.pos << std::endl;
         }
      }

      return typed_clusters;
   };

   auto get_mol = [] (int imol, const std::vector<cfc::input_info_t> &mol_infos) {
      mmdb::Manager *mol = nullptr;
      for (const auto &m : mol_infos) {
         if (m.imol == imol) {
            mol = m.mol;
         }
      }
      return mol;
   };

   auto output_clusters = [get_mol] (const std::vector<typed_cluster_t> &typed_clusters,
                                     const std::vector<cfc::input_info_t> &mol_infos) {

      // create a json string using nlohmann::json

      nlohmann::json j;
      for (const auto &tc : typed_clusters) {
         nlohmann::json j_cluster;
         j_cluster["family"] = tc.family;
         j_cluster["type"]   = tc.type;
         j_cluster["idx"]    = tc.idx;
         nlohmann::json j_imols = nlohmann::json::array();
         for (const auto &imol_with_specs : tc.imols_with_specs) {
            int imol = imol_with_specs.first;
            const auto &spec = imol_with_specs.second;
            nlohmann::json j_imol     = imol_with_specs.first;
            nlohmann::json j_chain_id = imol_with_specs.second.chain_id;
            nlohmann::json j_res_no   = imol_with_specs.second.res_no;
            nlohmann::json j_imol_and_res_spec;
            j_imol_and_res_spec["imol"] = j_imol;
            j_imol_and_res_spec["chain_id"] = j_chain_id;
            j_imol_and_res_spec["res_no"] = j_res_no;
            mmdb::Manager *mol = get_mol(imol, mol_infos);
            if (mol) {
               mmdb::Residue *residue_p = coot::util::get_residue(spec, mol);
               if (residue_p) {
                  std::string rn = residue_p->GetResName();
                  j_imol_and_res_spec["res_name"] = rn;
               }
            }
            j_imols.push_back(j_imol_and_res_spec);
         }
         j_cluster["imols"] = j_imols;
         j_cluster["pos"] = {tc.pos.x, tc.pos.y, tc.pos.z};
         j["clusters"].push_back(j_cluster);
      }
      std::string sjson = j.dump(4);
      // now write sjson to a file
      std::string fn = "clusters.json";
      std::ofstream f(fn.c_str());
      if (!f) {
         // std::cout << "ERROR:: output_clusters(): Failed to open " << fn << std::endl;
         logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                    "output_clusters() failed to open", fn);
      } else {
         f << sjson;
         f.close();
         logger.log(log_t::DEBUG, logging::function_name_t(__FUNCTION__), "wrote", fn);
      }
   };

   auto output_waters = [] (std::vector<std::vector<water_info_t> > water_clusters) {

      std::cout << "DEBUG:: output_waters() " << water_clusters.size() << std::endl;
      nlohmann::json j = nlohmann::json::array();
      for (unsigned int iclust=0; iclust<water_clusters.size(); iclust++) {
         std::vector<water_info_t> &v = water_clusters[iclust];
         nlohmann::json j_cluster;
         j_cluster["cluster_index"] = iclust;
         nlohmann::json j_waters = nlohmann::json::array();
         for (unsigned int i=0; i<v.size(); i++) {
            const water_info_t &wi = v[i];
            nlohmann::json j_water;
            j_water["imol"] = wi.imol;
            j_water["chain_id"] = wi.residue_spec.chain_id;
            j_water["res_no"] = wi.residue_spec.res_no;
            j_water["pos"] = {wi.pos.x, wi.pos.y, wi.pos.z};
            j_waters.push_back(j_water);
         }
         j_cluster["waters"] = j_waters;
         j.push_back(j_cluster);
      }
      std::string sjson = j.dump(4);
      std::string fn = "waters.json";
      std::ofstream f(fn.c_str());
      if (!f) {
         logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                    "output_clusters() failed to open", fn);
      } else {
         f << sjson;
         f.close();
         logger.log(log_t::DEBUG, logging::function_name_t(__FUNCTION__), "wrote", fn);
      }
   };

   auto output_waters_plain = [] (const std::vector<std::vector<water_info_t> > &water_clusters) {

      if (false)
         std::cout << "DEBUG:: output_waters_plain() " << water_clusters.size() << " water cluster"
                   << std::endl;
      for (unsigned int iclust=0; iclust<water_clusters.size(); iclust++) {
         const std::vector<water_info_t> &v = water_clusters.at(iclust);
         if (false)
            std::cout << "debug:: water cluster iclust " << iclust << std::endl;
         for (unsigned int i=0; i<v.size(); i++) {
            const water_info_t &wi = v[i];
            if (false)
               std::cout << "  debug:: output_waters_plain(): water iclust: " << iclust << " at "
                         << wi.pos.x << " " << wi.pos.y << " " << wi.pos.z << std::endl;
         }
      }
   };

   auto feature_is_close_to_a_protein_atom = [] (const feature_info_t &fi,
                                                 const std::vector<mmdb::Residue *> &residues_near_ligand,
                                                 float dist_crit) {

      for (unsigned int i=0; i<residues_near_ligand.size(); i++) {
         mmdb::Residue *residue_p = residues_near_ligand[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               float delta_x = at->x - fi.pos.x;
               float delta_y = at->y - fi.pos.y;
               float delta_z = at->z - fi.pos.z;
               float dd = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
               if (dd < dist_crit * dist_crit)
                  return true;
            }
         }
      }
      return false;
   };

   // --------------------- main line -------------------

   std::vector<feature_info_t> feature_infos;
   std::vector<water_info_t> water_infos;

   // the feature factory is made once, rather than for every ligand
   std::unique_ptr<RDKit::MolChemicalFeatureFactory> factory;
   std::string feature_data_file_name = find_feature_definitions_file();
   if (feature_data_file_name.empty()) {
      logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                 "BaseFeatures.fdef not found (looked in the current directory, RDKIT_BASEFEATURES and",
                 coot::rdkit_package_data_dir() + "/Data)");
   } else {
      std::string feature_data_string = file_to_string(feature_data_file_name);
      try {
         factory.reset(RDKit::buildFeatureFactory(feature_data_string));
      }
      catch (const std::exception &e) {
         logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                    {"failed to read", feature_data_file_name, e.what()});
      }
   }

   if (! mol_infos.empty() && factory) {

      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      int imol_ref = mol_infos[0].imol;
      mmdb::Manager *mol_ref = mol_infos[0].mol;

      // if the molecules are not superposed, they will be superposed onto the reference
      // (the first molecule) using the residues near the centre - or, if we don't have
      // a centre, near the ligand in the reference molecule.
      std::vector<coot::lsq_range_match_info_t> lsq_matchers;
      if (! options.molecules_are_superposed) {
         if (mol_ref) {
            bool have_pt = false;
            clipper::Coord_orth pt = options.centre;
            if (options.use_centre) {
               have_pt = true;
            } else {
               mmdb::Residue *ligand_ref = find_ligand_residue(mol_ref, mol_infos[0].res_name, false, pt);
               if (ligand_ref) {
                  auto rc = coot::util::get_residue_centre(ligand_ref);
                  if (rc.first) {
                     pt = rc.second;
                     have_pt = true;
                  }
               }
            }
            if (have_pt)
               lsq_matchers = make_lsq_matchers(mol_ref, pt, options.superposition_radius);
         }
         if (lsq_matchers.empty())
            logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                       "no residues for superposition onto reference", imol_ref);
      }

      for (unsigned int i=0; i<mol_infos.size(); i++) {

         mmdb::Manager *mol = mol_infos[i].mol;
         int imol = mol_infos[i].imol;
         if (! mol) {
            logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                       "no mol found for imol", imol);
            continue;
         }

         // rtop takes this molecule onto the reference molecule
         clipper::RTop_orth rtop = clipper::RTop_orth::identity();
         if (i > 0 && ! options.molecules_are_superposed) {
            bool superposed = false;
            if (! lsq_matchers.empty()) {
               std::cout << "info:: superposing ref: " << imol_ref << " mov: " << imol << std::endl;
               bool summary_to_screen = false;
               std::pair<short int, clipper::RTop_orth> rtop_info =
                  coot::util::get_lsq_matrix(mol_ref, mol, lsq_matchers, 1, summary_to_screen);
               if (rtop_info.first) {
                  rtop = rtop_info.second;
                  superposed = true;
               }
            }
            if (! superposed) {
               logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
                          "failed to superpose - ignoring imol", imol);
               continue;
            }
         }

         // residues

         clipper::Coord_orth centre_in_mol = rtop.inverse() * options.centre;
         const std::string &res_name = mol_infos[i].res_name;
         mmdb::Residue *residue_p = find_ligand_residue(mol, res_name, options.use_centre, centre_in_mol);
         if (residue_p) {
            std::vector<mmdb::Residue *> residues_near_residue = coot::residues_near_residue(residue_p, mol, 4.3);
            auto new_chemical_features = chemical_features(residue_p, imol, imol_enc, geom, factory.get());
            for (auto &fi : new_chemical_features) {
               float dist_crit = 4.3;
               // if (feature_is_close_to_a_protein_atom(fi, residues_near_residue, dist_crit)) {
               if (true) {
                  fi.pos = transform_point(rtop, fi.pos);
                  feature_infos.push_back(fi);
               }
            }

            // waters

            auto rc = coot::util::get_residue_centre(residue_p);
            if (rc.first) {
               std::vector<water_info_t> waters = waters_near_ligand_site(imol, mol, rc.second);
               for (auto &wi : waters) {
                  // a molecule may be listed more than once (with different ligands)
                  // so don't add the same water twice
                  auto same_water = [&wi] (const water_info_t &w) {
                     return w.imol == wi.imol && w.residue_spec == wi.residue_spec; };
                  if (std::find_if(water_infos.begin(), water_infos.end(), same_water) != water_infos.end())
                     continue;
                  wi.pos = transform_point(rtop, wi.pos);
                  water_infos.push_back(wi);
               }
            }
         } else {
            logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__),
                       {"no ligand residue", res_name, "found in imol", imol});
         }
      }
   }

   unsigned int n_features = feature_infos.size();
   if (! feature_infos.empty()) {
      std::cout << "DEBUG:: :::::::::::::::::::::::::::::::::: found a total of "
                << n_features << " features" << std::endl;
      logger.log(log_t::DEBUG, logging::function_name_t(__FUNCTION__),
                 "Found a total of", n_features, "features");

      std::vector<typed_cluster_t> typed_clusters = cluster_features(feature_infos);
      std::vector<std::vector<water_info_t> > water_clusters_dir; // = cluster_waters_dirichlet(water_infos);

      std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
      std::vector<std::vector<water_info_t> > water_clusters_gmm = cluster_waters_gmm(water_infos);
      std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;

      std::vector<std::vector<water_info_t> > water_clusters = water_clusters_gmm;

      auto water_sorter = +[] (const std::vector<water_info_t> &wi1,
                               const std::vector<water_info_t> &wi2) {
         return wi2.size() < wi1.size();
      };
      std::sort(water_clusters.begin(), water_clusters.end(), water_sorter);

      // std::cout << "Found " << water_clusters.size() << " water clusters from "
      //           << water_infos.size() << " water" << std::endl;
      logger.log(log_t::DEBUG, logging::function_name_t(__FUNCTION__),
                 {"Found", water_clusters.size(), "water clusters from",
                  water_infos.size(), "waters"});

      output_clusters(typed_clusters, mol_infos);

      output_waters(water_clusters);
      output_waters_plain(water_clusters); // to terminal

      if (false) {
         std::cout << "--------- at end of chemical_feature_clustering() --- "
                   << std::endl;
         for (unsigned int i=0; i<water_clusters.size(); i++) {
            std::set<int> imols_in_cluster;
            const auto &wi = water_clusters[i];
            for (unsigned int jj=0; jj<wi.size(); jj++) {
               int imol = wi[jj].imol;
               imols_in_cluster.insert(imol);
            }
            std::cout << "debug:: water_clusters[" << i << "] has contributions from "
                      << imols_in_cluster.size() << " molecules"
                      << std::endl;
         }
         std::cout << "------------" << std::endl;
      }

      return std::make_pair(typed_clusters, water_clusters);

   } else {
      logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__), "empty feature_infos");
   }

   // return dummy
   std::vector<typed_cluster_t> d1;
   std::vector<std::vector<water_info_t> > d2;
   return std::make_pair(d1, d2);
}

#endif // MAKE_ENHANCED_LIGAND_FEATURES

