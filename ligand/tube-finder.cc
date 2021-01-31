
#include <fstream>
#include <chrono>

#include "analysis/stats.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/merge-molecules.hh"
#include "coot-utils/merge-atom-selections.hh"
#include "coot-utils/fib-sphere.hh"
#include "helix-placement.hh"
#include "tube-finder.hh"


coot::tube_finder_t::tube_finder_t(const clipper::Xmap<float> &xmap) {

   float b_factor = 200;
   unsigned n_test_points = 1500;
   clipper::Xmap<float> blurred_xmap = util::sharpen_blur_map(xmap, b_factor);

   std::pair<float, float> mv = util::mean_and_variance(blurred_xmap);
   float sd = sqrt(mv.second);
   float density_level_crit = 2.0 * sd;

   util::map_molecule_centre_info_t mmci = util::map_molecule_centre(xmap);
   if (mmci.success) {
      clipper::Cell_descr cell = xmap.cell();
      float quarter_x = 0.25 * cell.a();
      float quarter_y = 0.25 * cell.b();
      float quarter_z = 0.25 * cell.c();
      float half_x = 0.5 * cell.a();
      float half_y = 0.5 * cell.b();
      float half_z = 0.5 * cell.c();

      clipper::Coord_orth mmm = mmci.updated_centre; // middle of the molecule blob
      std::cout << "mmm: " << mmm.format() << std::endl;
      float fudge = 0.5;
      std::pair<float, float> x_limits(mmm.x() - fudge * quarter_x, mmm.x() + fudge * quarter_x);
      std::pair<float, float> y_limits(mmm.y() - fudge * quarter_y, mmm.y() + fudge * quarter_y);
      std::pair<float, float> z_limits(mmm.z() - fudge * quarter_z, mmm.z() + fudge * quarter_z);

      std::vector<clipper::Coord_orth> test_points;

      float rmi = 1.0/float(RAND_MAX);
      for (unsigned int i=0; i<n_test_points; i++) {
         float rx = x_limits.first + random() * rmi * half_x * fudge;
         float ry = y_limits.first + random() * rmi * half_y * fudge;
         float rz = z_limits.first + random() * rmi * half_z * fudge;
         // rx = x_limits.first + quarter_x;
         // ry = y_limits.first + quarter_y;
         // rz = z_limits.first + quarter_z;
         clipper::Coord_orth pt(rx, ry, rz);
         test_points.push_back(pt);
      }

      // points are along the X axis -8, -4, 0, 4, 8
      std::vector<clipper::Mat33<double> > orientations;
      clipper::Mat33<double> r0(1,0,0,  0,1,0,  0,0,1);
      clipper::Mat33<double> r1(0,1,0, -1,0,0,  0,0,1);
      clipper::Mat33<double> r2(0,0,1,  0,1,0, -1,0,0);
      orientations.push_back(r0);
      // orientations.push_back(r1);
      // orientations.push_back(r2);
      // more here

      std::vector<clipper::Coord_orth> potential_centres;
      auto tp_0 = std::chrono::high_resolution_clock::now();
      for (unsigned int i=0; i<n_test_points; i++) {
         for (unsigned int i_ori=0; i_ori<orientations.size(); i_ori++) {
            //std::cout << "constructor test_point i " << i << " " << test_points[i].format() << std::endl;
            std::pair<bool, clipper::Coord_orth> fitted_point =
               fit_stick_in_map(test_points[i], orientations[i_ori], blurred_xmap, density_level_crit);
            if (fitted_point.first) {
               std::cout << i << " " << i_ori << std::endl;
               std::cout << fitted_point.second.format();
               clipper::Coord_orth pt = fitted_point.second;
               std::ofstream f("potential-tube-centres.table", std::ofstream::app);
               f << pt.x() << " ";
               f << pt.y() << " ";
               f << pt.z() << "\n";
               potential_centres.push_back(pt);
            }
         }
      }

      // debugging:
      // potential_centres.clear();
      // potential_centres.push_back(clipper::Coord_orth(148.8, 152.6, 133));
      
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << " timings: " << d10 << " ms, which is " << d10/n_test_points << " per point"
                << std::endl;
      // variance scored
      unsigned  int n_c = potential_centres.size();
      std::vector<clipper::Coord_orth> sphere_points = fibonacci_sphere(50);
      float radius = 4.6;
      std::vector<std::pair<float, clipper::Coord_orth> > scored_potential_centres(n_c);
      for (unsigned int i=0; i<potential_centres.size(); i++) {
         const clipper::Coord_orth &pt = potential_centres[i];
         float var = sphere_variance(pt, sphere_points, radius, blurred_xmap);
         scored_potential_centres[i] = std::pair<float, clipper::Coord_orth>(var, pt);
      }

      auto score_func = [] (const std::pair<float, clipper::Coord_orth> &pt1,
                            const std::pair<float, clipper::Coord_orth> &pt2) {
                           return pt2.first < pt1.first;
                        };
      std::sort(scored_potential_centres.begin(),
                scored_potential_centres.end(), score_func);

      // debug
      for (unsigned int i=0; i<scored_potential_centres.size(); i++)
         std::cout << i << " "
                   << scored_potential_centres[i].first << " "
                   << scored_potential_centres[i].second.format() << std::endl;

      // maybe not so great density fit though!
      std::vector<helix_placement_info_t> successfully_placed_helices;

      coot::helix_placement p(xmap);
      float min_density_limit = 8.5 * sd;
      float high_density_turning_point = 20 * sd;
      float bf = 40;
      float map_rmsd = sd;
      std::cout << "INFO:: density_level for trim " << min_density_limit << std::endl;
      std::cout << "INFO:: density_level for score turning point " << high_density_turning_point << std::endl;
      for (unsigned int i=0; i<scored_potential_centres.size(); i++) {
         const clipper::Coord_orth &pt = scored_potential_centres[i].second;
         int n_helix_residues = 20;
         // map_rmsd is used for step size in rigid-body fitting.
         helix_placement_info_t n =
            p.place_alpha_helix_near_kc_version(pt, n_helix_residues,
                                                min_density_limit, high_density_turning_point, bf, map_rmsd);
         if (n.success) {
            successfully_placed_helices.push_back(n);
         }
      }

      auto score_helix_position = [] (const coot::minimol::molecule &m, const clipper::Xmap<float> &xmap) {
                                     float score = 0.0;
                                     for (unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
                                        for (int ires=m[ifrag].min_res_no(); ires<m[ifrag].max_residue_number(); ires++) {
                                           for (unsigned int iat=0; iat<m[ifrag][ires].atoms.size(); iat++) {
                                              score += util::density_at_point(xmap, m[ifrag][ires][iat].pos) *
                                                 m[ifrag][ires][iat].occupancy;
                                           }
                                        }
                                     }
                                     return score; 
                                  };
      
      std::vector<scored_helix_info_t> scored_helices(successfully_placed_helices.size());
      for (unsigned int i=0; i<successfully_placed_helices.size(); i++) {
         const minimol::molecule &mol = successfully_placed_helices[i].mol[0];
         float score = score_helix_position(mol, xmap);
         scored_helices[i].score = score;
         scored_helices[i].mol = mol;
      }

      auto helix_sorter_func = [] (const scored_helix_info_t &h1,
                                   const scored_helix_info_t &h2) {
                                  return h2.score < h1.score;
                               };

      std::sort(scored_helices.begin(), scored_helices.end(), helix_sorter_func);

      if (false)  {
         for (unsigned int i=0; i<scored_helices.size(); i++) {
            std::cout << "helix_index " << i << " score " << scored_helices[i].score << std::endl;
            std::string file_name = "sorted-placed-helix-" + std::to_string(i) + ".pdb";
            float bf = 40;
            scored_helices[i].mol.write_file(file_name, bf);
         }
      }

      unsigned int n_top = 20;
      if (scored_helices.size() > n_top)
         scored_helices.resize(n_top);
      else
         n_top = scored_helices.size();

      for (unsigned int i=0; i<n_top; i++) {
         std::cout << "helix_index " << i << " score " << scored_helices[i].score << std::endl;
         std::string file_name = "sorted-placed-helix-" + std::to_string(i) + ".pdb";
         float bf = 40;
         scored_helices[i].mol.write_file(file_name, bf);
      }
      if (n_top > 0) {
         mmdb::Manager *mol_first = scored_helices[0].mol.pcmmdbmanager();
         std::vector<mmdb::Manager *> mol_others;
         for (unsigned int i=1; i<scored_helices.size(); i++) {
            const minimol::molecule &m = scored_helices[i].mol;
            mmdb::Manager *mol =m.pcmmdbmanager();
            mol_others.push_back(mol);
         }

         merge_molecules(mol_first, mol_others);

         // check that there are the right number of chains in mol_first
         mmdb::Model *model_p = mol_first->GetModel(1);
         int n_chains = model_p->GetNumberOfChains();
         std::cout << "################################ " << n_chains << " chains ###########################"
                   << std::endl;
         mol_first->WritePDBASCII("pre-merge-fragments.pdb");
         merge_atom_selections(mol_first);

         mol_first->WritePDBASCII("merged.pdb");
         mol_first->WriteCIFASCII("merged.cif");

      }
   }
}


// If you edit this function again, get rid of it.
std::pair<bool, clipper::Coord_orth>
coot::tube_finder_t::fit_stick_in_map(const clipper::Coord_orth &test_point,
                                      const clipper::Mat33<double> &orientation,
                                      const clipper::Xmap<float> &xmap,
                                      float density_level_crit) const {

   std::pair<bool, clipper::Coord_orth> fitted_point =
      fit_to_map_by_simplex_rigid(test_point, orientation, xmap, density_level_crit);
   return fitted_point;
}

// static
clipper::RTop_orth
coot::tube_finder_t::construct_matrix(const gsl_vector *v) {

   clipper::Coord_orth trans(gsl_vector_get(v, 0),
			     gsl_vector_get(v, 1),
			     gsl_vector_get(v, 2));
   
   double sin_t;
   double cos_t;
   double f = 0.2;

   sin_t = sin(-clipper::Util::d2rad(f * gsl_vector_get(v, 3)));
   cos_t = cos(-clipper::Util::d2rad(f * gsl_vector_get(v, 3)));
   clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

   sin_t = sin(-clipper::Util::d2rad(f * gsl_vector_get(v, 4)));
   cos_t = cos(-clipper::Util::d2rad(f * gsl_vector_get(v, 4)));
   clipper::Mat33<double> y_mat(cos_t,0,sin_t, 0,1,0, -sin_t,0,cos_t);

   sin_t = sin(-clipper::Util::d2rad(f * gsl_vector_get(v, 5)));
   cos_t = cos(-clipper::Util::d2rad(f * gsl_vector_get(v, 5)));
   clipper::Mat33<double> z_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);

   clipper::Mat33<double> angle_mat = x_mat * y_mat * z_mat;

   clipper::RTop_orth rtop(angle_mat, trans);

   return rtop;
}


// static
double
coot::tube_finder_t::my_f_simplex_rigid_internal(const gsl_vector *v, void *params) {

   simplex_param_t *p = static_cast<simplex_param_t *>(params);

   double score = 0;

   if (false)
      std::cout << "in my_f_simplex_rigid_internal: v: "
                << gsl_vector_get(v, 0) << " "
                << gsl_vector_get(v, 1) << " "
                << gsl_vector_get(v, 2) << " "
                << gsl_vector_get(v, 3) << " "
                << gsl_vector_get(v, 4) << " "
                << gsl_vector_get(v, 5) << " "
                << std::endl;

   clipper::Coord_orth point;
   unsigned int n_points_in_stick = 5;
   clipper::RTop_orth rtop = construct_matrix(v);

   for (unsigned int i=0; i<n_points_in_stick; i++) {
      const clipper::Coord_orth &orig_pt = p->original_positions[i];
      // ::cout << "orig_pt: " << orig_pt.format() << std::endl;
      point = p->test_point_centre_orig + (orig_pt - p->test_point_centre_orig).transform(rtop);
      // std::cout << "point for density test:   " << point.format() << std::endl;
      // we are trying to minimize, don't forget:
      score -= util::density_at_point(*(p->xmap), point);
   }

   // std::cout << "my_f_simplex_rigid_internal() score: " << score << std::endl;
   return score;
}

std::pair<bool, clipper::Coord_orth>
coot::tube_finder_t::fit_to_map_by_simplex_rigid(const clipper::Coord_orth &test_point,
                                                 const clipper::Mat33<double> &orientation,
                                                 const clipper::Xmap<float> &xmap,
                                                 float density_level_crit) const {

   bool fit_status = false;
   clipper::Coord_orth fitted_point(0,0,0);
   // std::cout << "fit_to_map_by_simplex_rigid:: test_point " << test_point.format() << std::endl;

   const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
   gsl_multimin_fminimizer *s = NULL;
   gsl_multimin_function minex_func;

   int rval = GSL_CONTINUE;
   int status = GSL_SUCCESS;
   simplex_param_t par;
   par.test_point_centre_orig = test_point;
   int num_params  = 6; // 3 translation, 3 rotation

   gsl_vector *ss = gsl_vector_alloc(num_params); // step size

   if (ss == NULL) {
      std::cout << "ERROR:: fit_to_map_by_simplex_rigid() no memory" << std::endl;
      return std::pair<bool, clipper::Coord_orth> (false, clipper::Coord_orth(0,0,0));
   }

   gsl_vector *x = gsl_vector_alloc(num_params);  // for the "atoms" on the stick

   par.xmap = &xmap;
   // now fill x with positions
   clipper::Coord_orth zero(0,0,0);
   for (int i= -2; i<3; i++) {
      clipper::Coord_orth pt_delta(static_cast<float>(i)*3.8, 0, 0);
      clipper::Coord_orth rot_pt = pt_delta.transform(clipper::RTop_orth(orientation, zero));
      clipper::Coord_orth position = test_point + rot_pt;
      par.original_positions.push_back(position);
   }

   // debuggging
   par.original_positions.clear();
   par.original_positions.push_back(test_point);

   gsl_vector_set_all(ss, 0.2); // step size
   gsl_vector_set_all(x,  0.0); // no rotations or translations initially.

   // setup_simplex_x_internal(x, atom_selection, n_selected_atoms);

   minex_func.f = my_f_simplex_rigid_internal;
   minex_func.n = num_params;
   minex_func.params = static_cast<void *>(&par); // par stays in scope for this function

   s = gsl_multimin_fminimizer_alloc(T, num_params);
   gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

   size_t iter = 0;
   while (rval == GSL_CONTINUE && iter < 2000) {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
	 break;
      rval = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), 1e-2);
      double ssval = gsl_multimin_fminimizer_size(s);

      // std::cout << "ssval " << ssval << " rval " << rval << std::endl;

      if (rval == GSL_SUCCESS) {
	 std::cout << "converged at minimum\n";
	 status = 1;
         // now apply the parameters to the positions
	 fitted_point = apply_shifts_rigid_internal(s->x, par);
         if (util::density_at_point(xmap, fitted_point) > density_level_crit)
            fit_status = true;
      }

      if (false) {
         std::cout << "iterating " << iter << ": "
                   << gsl_vector_get(s->x, 0) << " "
                   << gsl_vector_get(s->x, 1) << " "
                   << gsl_vector_get(s->x, 2) << " "
                   << gsl_vector_get(s->x, 3) << " "
                   << gsl_vector_get(s->x, 4) << " "
                   << gsl_vector_get(s->x, 5) << " ";
         std::cout << "f() " << s->fval << " step-size " << ssval << "\n";
      }
   }


   gsl_vector_free(x);
   gsl_vector_free(ss);
   gsl_multimin_fminimizer_free(s);

   // maybe return success here, in a pair.
   if (! status) std::cout << "fail" << std::endl;
   return std::pair<bool,clipper::Coord_orth> (fit_status, fitted_point);
}

std::vector<clipper::Coord_orth>
coot::tube_finder_t::get_positions() const {
   return positions;  
}


clipper::Coord_orth
coot::tube_finder_t::apply_shifts_rigid_internal(gsl_vector *x, const simplex_param_t &par) const {

   clipper::Coord_orth position_delta(gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2));
   const clipper::Coord_orth &orig_pt = par.test_point_centre_orig;
   clipper::Coord_orth pt = orig_pt + position_delta;

   // std::cout << "apply_shifts_rigid_internal() orig_pt: " << orig_pt.format() << std::endl;
   //  std::cout << "apply_shifts_rigid_internal() point:   " << pt.format() << std::endl;

   return pt;
}


// static
float
coot::tube_finder_t::sphere_variance(const clipper::Coord_orth &centre_point,
                                     const std::vector<clipper::Coord_orth> &points,
                                     float radius,
                                     const clipper::Xmap<float> &xmap) {

   float var = 0;

   stats::single data;
   for (unsigned int i=0; i<points.size(); i++) {
      const clipper::Coord_orth &pt = points[i];
      float d = util::density_at_point(xmap, radius * pt + centre_point);
      data.add(d);
   }

   var = data.variance();

   return var;
}
