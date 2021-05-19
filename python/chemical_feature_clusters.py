#  Copyright 2016 by Medical Research Council
#  Author: Paul Emsley
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or (at
#  your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#  02110-1301, USA

from __future__ import print_function
import math
import itertools
import coot

import numpy as np
from scipy import linalg
from sklearn import mixture

import res_spec_utils as rsu


def cluster_star_obj(obj, pos, thick, v_0):
    
    delta_0 = v_0
    delta_r2 = delta_0/math.sqrt(2)
    delta_r3 = delta_0/math.sqrt(3)
    delta = delta_r3

    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]+delta, pos[1]+delta, pos[2]+delta)
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]-delta, pos[1]-delta, pos[2]-delta)
    delta = delta_0
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]-delta, pos[1],       pos[2])
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0],       pos[1]-delta, pos[2])
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0],       pos[1],       pos[2]-delta)
    delta = delta_r2
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0],       pos[1]-delta, pos[2]-delta)
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]-delta, pos[1]-delta, pos[2])
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]-delta, pos[1],       pos[2]-delta)

    delta = delta_0
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]+delta, pos[1],       pos[2])
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0],       pos[1]+delta, pos[2])
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0],       pos[1],       pos[2]+delta)
    delta = delta_r2
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0],       pos[1]+delta, pos[2]+delta)
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]+delta, pos[1]+delta, pos[2])
    coot.to_generic_object_add_line(obj, "yellowtint", thick,
        pos[0],       pos[1],       pos[2],
        pos[0]+delta, pos[1],       pos[2]+delta)


def get_residue_name(imol, res_spec):
    return coot.residue_name(imol,
                             rsu.residue_spec_to_chain_id(res_spec),
                             rsu.residue_spec_to_res_no(res_spec),
                             "")

# return the dpgmm, so that we can get weights, means, covar and predictions
#
def cluster_and_display_waters(site_number, w_positions_np):

   def optimize_n(positions_np, n_data):

       """return the optimal value for n (that maximizes bic)"""

       bic = {}
       for n in [x+1 for x in range(50)]:
          if n < len(positions_np):
              # gmm = mixture.GMM(n_components=n, covariance_type='spherical', n_iter=20) # old version of sklearn
             gmm = mixture.GaussianMixture(n_components=n, covariance_type='spherical', max_iter=200)
             gmm.fit(positions_np)
             # score = sum(gmm.score(positions_np))
             score = gmm.score(positions_np)
             lambda_c = 15 # 3 too few
             lambda_c = 0.06
             lambda_c = 0.02
             bic_l = score - lambda_c * 0.5 * math.log(n_data) * n
             print("debug:: optimize_n: ", n, score, bic_l)
             bic[n] = bic_l

       for key in bic:
           print("   debug:: water bic", key, bic[key])

       key,value = max(iter(bic.items()), key=lambda x:x[1])
       return key


   print("cluster_and_display_waters(): ########### size of w_positions_np: ", len(w_positions_np))

   n_components = optimize_n(w_positions_np, len(w_positions_np))
   print("optimize_n for water::::::::::::: n_components: ", n_components)
   # dpgmm = mixture.GMM(n_components, covariance_type='spherical', n_iter=40)  # old
   dpgmm = mixture.GaussianMixture(n_components, covariance_type='spherical', max_iter=40)
   dpgmm.fit(w_positions_np)

   cluster_assignments = dpgmm.predict(w_positions_np)
   
   color_list=['green', 'greentint', "sea", 'yellow', "yellowtint", "aquamarine", "forestgreen",
	       "goldenrod", "orangered", "orange", "cyan", 'red', "blue"]
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)
   color_list.extend(color_list)

   means   = dpgmm.means_
   # cvs     = dpgmm._get_covars()
   cvs     = dpgmm.covariances_
   weights = dpgmm.weights_

   print("debug:: cluster_and_display_waters(): ################ means:", means)
   print("debug:: cluster_and_display_waters(): ################ weights:", weights)
   print("debug:: cluster_and_display_waters(): ################ cvs:", cvs)

   obj = coot.new_generic_object_number("CFC Site " + str(site_number) + " selected waters")
   for i,pos in enumerate(w_positions_np):
       mean = means[cluster_assignments[i]]
       # reject spheres at the origin - (from DPGMM strangeness)
       d = mean[0]*mean[0] + mean[1]*mean[1] + mean[2]*mean[2]
       if d > 1.0:
	   col = color_list[cluster_assignments[i]]
	   coot.to_generic_object_add_point(obj, col, 10, pos[0], pos[1], pos[2])
       else:
	   print("reject prediction", i, "for cluster", cluster_assignments[i])

   # set_display_generic_object(obj, 1)

   obj = coot.new_generic_object_number("CFC Site " + str(site_number) + " water cluster means")

   for i,cv in enumerate(cvs):
       print("debug:: cluster_and_display_waters(): ################ cvs: i", i, "cv", cv)

   for i,cv in enumerate(cvs):

      mean = means[i]
      d = mean[0]*mean[0] + mean[1]*mean[1] + mean[2]*mean[2]

      # v,w = linalg.eigh(cv)
      # print "mean  ", mean
      # print "weight", weights[i], "prec", precs[i]
      # print "weight", weights[i]
      # print "v", v
      
      # if d > 1.0:

      #     pos = mean
      #     thick = 2
      #     cluster_star_obj(obj, pos, thick, v[0])

      # else:
      #     print("reject", mean, v)

      if d > 1.0:
          pos = mean
          cluster_star_obj(obj, pos, 2, 1.0)

   coot.set_display_generic_object(obj, 1)

   cluster_assignments_as_list = [int(x) for x in cluster_assignments]

   return (dpgmm, cluster_assignments_as_list)


def cluster_and_display_chemical_features(site_number, type, chemical_features_list):


   def optimize_n(type, positions_np, n_data):

      print("cluster_and_display_chemical_features.optimize_n called " \
             "with n_data = ", n_data)

      bic = {}
      for n in [x+1 for x in range(10)]:
         if n < n_data:
            # gmm = mixture.GMM(n_components=n, covariance_type='spherical', n_iter=20)
            gmm = mixture.GaussianMixture(n_components=n, covariance_type='spherical', max_iter=20)
            gmm.fit(positions_np)
            # score = sum(gmm.score(positions_np))
            score = gmm.score(positions_np)
            lambda_c = 15
            if type == 'Aromatic':
               lambda_c = 20
            bic_l = score - lambda_c * 0.5 * math.log(n_data) * n
            bic[n] = bic_l

      if len(bic) > 1:
          key, value = max(iter(bic.items()), key=lambda x:x[1])
          return key
      else:
          return 1

   def analyse_bic(type, positions_np, n_data):

      for n in [x+1 for x in range(14)]:
         # gmm = mixture.GMM(n_components=n, covariance_type='spherical', n_iter=20)
         gmm = mixture.GaussianMixture(n_components=n, covariance_type='spherical', max_iter=20)
         gmm.fit(positions_np)
         # score = sum(gmm.score(positions_np))
         score = gmm.score(positions_np)
         lambda_c = 3
         if type == 'Aromatic':
             lambda_c = 3000
         bic = score - lambda_c * 0.5 * n_data * n
         print(type, len(positions_np), n, "converged?", gmm.converged_, "score:", score, "bic", bic)

   def get_cfc_col(type):
       if type == "Donor":
           return "blue"
       if type == "Acceptor":
           return "red"
       if type == "Hydrophobe":
           return "yellow"
       if type == "Aromatic":
           return "orange"
       return "grey"

   # --- main line ----

   # no fake points
   # positions_np = np.array([item[0] for item in chemical_features_list])

   ext_chemical_features_list = [item[0] for item in chemical_features_list]
   
   for item_b in chemical_features_list:
      delta = 0.25
      item = item_b[0]
      p1  = [item[0],       item[1],       item[2]+delta]
      p2  = [item[0],       item[1],       item[2]-delta]
      p3  = [item[0],       item[1]+delta, item[2]      ]
      p4  = [item[0],       item[1]-delta, item[2]      ]
      p5  = [item[0]+delta, item[1],       item[2]      ]
      p6  = [item[0]-delta, item[1],       item[2]      ]
      ext_chemical_features_list.append(p1)
      ext_chemical_features_list.append(p2)
      ext_chemical_features_list.append(p3)
      ext_chemical_features_list.append(p4)
      ext_chemical_features_list.append(p5)
      ext_chemical_features_list.append(p6)

   positions_np = np.array(ext_chemical_features_list)

   # analyse_bic(type, positions_np, len(chemical_features_list))

   n_data = len(chemical_features_list)
   n = 1
   if n_data > 1:
       n = optimize_n(type, positions_np, n_data)
       
   if n <= len(chemical_features_list):
       # gmm = mixture.GMM(n_components=n, covariance_type='spherical', n_iter=20)
       gmm = mixture.GaussianMixture(n_components=n, covariance_type='spherical', max_iter=20)
       gmm.fit(positions_np)
       # score = sum(gmm.score(positions_np))
       score = gmm.score(positions_np)
       print(type, len(positions_np), n, "converged? ", gmm.converged_, "score:", score)

       cluster_assignments = gmm.predict(positions_np)

       features = []
       for i,cf in enumerate(chemical_features_list):
           # print "     ", cf, cluster_assignments[i]
           features.append([cf, int(cluster_assignments[i])])

       means = gmm.means_
       means_as_list = [[x[0], x[1], x[2]] for x in means]

       obj_name = "CFC Site " + str(site_number) + " " + type + " pharmacophore-clusters"
       cfc_obj = coot.new_generic_object_number(obj_name)
       cfc_col = get_cfc_col(type)
       for mean in means_as_list:
           # coot.to_generic_object_add_dodecahedron(cfc_obj, cfc_col, 0.2, mean[0], mean[1], mean[2])
           coot.to_generic_object_add_pentakis_dodecahedron(cfc_obj, cfc_col, 2.3, 0.1, mean[0], mean[1], mean[2])
       coot.set_display_generic_object(cfc_obj, 1)
           
       return [type, features, means_as_list]

   # oops too many parameters for the model
   return False


def make_ball_and_stick_by_spec(imol, ligand_spec):

    # this function no longer worked with a prefixed ligand spec. So we assume 3 elements

    print("debug:: in make_ball_and_stick_by_spec: imol", imol, "spec:", ligand_spec)

    if len(ligand_spec) == 3:
        s = "//" + ligand_spec[0] + "/" + str(ligand_spec[1])
        coot.make_ball_and_stick(imol, s, 0.125, 0.125, 1)



# find the chemical feature and water clusters of the give site.
# call this in turn - with a ranking based on the number of structures per site.
# (The first time this is called should be the site that has the most structures)
#
def cfc_process_site(site_number, wrapped_imol_ligand_specs, imol_first, first_ligand_spec):

   imol_ligand_specs = wrapped_imol_ligand_specs[0] # unwrap

   print("debug:: in cfc_process_site() with wrapped_imol_ligand_specs", wrapped_imol_ligand_specs)
   print("debug:: in cfc_process_site() with imol_ligand_specs", imol_ligand_specs)
   print("debug:: in cfc_process_site() with non-first imol_ligand_specs", imol_ligand_specs[1:])
   print("debug:: in cfc_process_site() with imol_first", imol_first)

   for ls in imol_ligand_specs:
       print("debug:: in cfc_process_site() with ligand", ls)

   print("debug:: in cfc_process_site() calling residues_near_residue_py() with args", imol_first, first_ligand_spec)
   env_residue_specs = coot.residues_near_residue_py(imol_first, first_ligand_spec, 6)
   # print("env_residue_specs", env_residue_specs)
   protein_res_specs = [r for r in env_residue_specs if get_residue_name(imol_first, r) != "HOH"]

   # only lsq the first (0th) one - that one has the most ligands in the site
   #
   if site_number == 0:
       # print("protein_res_specs (for lsqing):")
       # for spec in protein_res_specs:
       #     print("   ", spec, get_residue_name(imol_first, spec))

       print("debug:: ================= in cfc_process_site() count the protein_res_specs: ", len(protein_res_specs))

       for res_spec in protein_res_specs:
	   chain_id = rsu.residue_spec_to_chain_id(res_spec)
	   res_no   = rsu.residue_spec_to_res_no(res_spec)
	   coot.add_lsq_match(res_no, res_no, chain_id, res_no, res_no, chain_id, 1)
           print("debug:: ================ in cfc_process_site() adding protein residue", res_spec)

       for imol_and_spec in imol_ligand_specs[1:]:  # lsq fit others to the first in the list
	   print('============================ lsq-match ', imol_first, imol_and_spec, imol_and_spec[0])
           imol,spec = imol_and_spec
	   # coot.apply_lsq_matches_py(imol_first, imol_and_spec[0])
	   coot.apply_lsq_matches_py(imol_first, imol)
           # 4-element (i.e. prefixed) residue specs are need to be phased out, but for now, test and fix
           if len(spec) == 4:
               spec = spec[1:]
           make_ball_and_stick_by_spec(imol, spec)
	   # pass

       print("Here with first_ligand_spec:", first_ligand_spec)
       ligand_centre = coot.residue_centre_py(imol_first,
                                              rsu.residue_spec_to_chain_id(first_ligand_spec),
                                              rsu.residue_spec_to_res_no(first_ligand_spec), '')
       coot.set_go_to_atom_molecule(imol_first)
       coot.set_rotation_centre(*ligand_centre)
                                    

   combo_list = []
   try: 

       # we have a large radius for the water selection          
       radius = 10  # water must be within radius of it's own ligand
       radius_2 = 5 # water must be with radius_2 of any ligand atom (not just its own)
       
       combo_list = coot.chemical_feature_clusters_py(env_residue_specs,
                                                      imol_ligand_specs,
                                                      radius, radius_2)
       print("debug:: ########### in cfc_process_site() env_residue_specs is ", env_residue_specs)
       print("debug:: ########### in cfc_process_site() imol_ligand_specs is ", imol_ligand_specs)
       print("debug:: ########### in cfc_process_site() radius is ", radius)
       print("debug:: ########### in cfc_process_site() radius_2 is ", radius_2)

   except TypeError as e:
       print(e)

       # the rest is unlikely to work if we get here

   if True:

       water_position_list      = combo_list[0]
       chemical_feature_list    = combo_list[1]
       # residues_sidechains_list = combo_list[1]

       # ----------- handle waters -----------

       w_positions_list = []

       for item in [wat[2] for wat in water_position_list]:
           w_positions_list.append(item)

       for w in w_positions_list:
           print("debug:: ########### in cfc_process_site() water start position list", w)

       for item in [wat[2] for wat in water_position_list]:
           delta = 0.1
           p1  = [item[0],       item[1],       item[2]+delta]
           p2  = [item[0],       item[1],       item[2]-delta]
           p3  = [item[0],       item[1]+delta, item[2]      ]
           p4  = [item[0],       item[1]-delta, item[2]      ]
           p5  = [item[0]+delta, item[1],       item[2]      ]
           p6  = [item[0]-delta, item[1],       item[2]      ]

           w_positions_list.append(p1)
           w_positions_list.append(p2)
           w_positions_list.append(p3)
           w_positions_list.append(p4)
           w_positions_list.append(p5)
           w_positions_list.append(p6)

       w_positions_np = np.array(w_positions_list)

       if True: # debug water positions
           f = open("debug-waters-positions.table", "w")
           for pos in water_position_list:
               f.write("water_position: ")
               f.write(str(pos))
               f.write("\n")
           f.close()

       print("debug:: ########### in cfc_process_site() water_position_list size is ", len(water_position_list))
       print("debug:: ########### in cfc_process_site() w_positions_np size is ", len(w_positions_np))

       # move these to the origin
       # w_positions_np = w_positions_np_at_ligand
       # for pos in w_positions_np:
       #     pos -= np.array(ligand_centre)

       # dpgmm = mixture.DPGMM(n_components=25, covariance_type='spherical', alpha=1.101,
       #                       n_iter=40000, params='wmc', init_params='wmc', tol=1e-4,
       #                       verbose=0)
       # 
       # the number of clusters is highly related to the dist_cutoff (the
       # distance of an accepted water atom to any any atom in any of the
       # ligands = currently 4.2)
       #
       gmm, cluster_assignments = cluster_and_display_waters(site_number, w_positions_np)

       means   = gmm.means_
       # cvs     = gmm._get_covars()
       cvs     = gmm.covariances_ # not a square matrix now
       weights = gmm.weights_

       print("water means:")
       for mean in means:
           print("   ", mean)

       # each water has been assigned a cluster, that is the cluster_assignments
       #   
       # need to convert the array cluster_assignments to a list of items:
       #   [imol water_residue_spec cluster_number]
       #
       water_cluster_info_for_input = []
       for i, water_pos in enumerate(water_position_list):
           # print water_pos, cluster_assignments[i]
           item = [water_pos[0], water_pos[1], cluster_assignments[i]]
           water_cluster_info_for_input.append(item)

       # cluster_info is a list of
       #  list of water cluster info
       #      list of [position, weight, length]  where length is the eigenvalue v[0],
       #              (same as v[1], v[2] - all the same for spherical model)
       #      list of cluster predictions for then input positions
       #
       fake_radii = [0.03 for i in range(len(weights))]
       print("debug:: weights", weights)
       # ci = list(zip([[l[0],l[1],l[2]] for l in means], weights, [cv[0][0] for cv in cvs]))
       ci = list(zip([[l[0],l[1],l[2]] for l in means], weights, fake_radii))
       water_cluster_info = [ci, water_cluster_info_for_input]
       # give those results back to c++ so that we can use them for display
       #

       coot.set_display_generic_objects_as_solid(1)

       # ----------- handle chemical features -----------

       # make a dictionary from the list of chemical features
       chemical_features_dict = {}
       for item in chemical_feature_list:
           for type in ['Donor', 'Acceptor', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe']:
               if item[0] == type:
                   try:
                       chemical_features_dict[type].append(item[1:])
                   except KeyError:
                       chemical_features_dict[type] = [item[1:]]

       chemical_feature_clusters_info = []
       for key in chemical_features_dict:
           # list of [type, features-annotated-by-cluster-number, cluster_means]
           clusters = cluster_and_display_chemical_features(site_number, key, chemical_features_dict[key])
           chemical_feature_clusters_info.append(clusters)

       # print 'water_cluster_info'
       # for wc in water_cluster_info:
       #    print wc

       cluster_info = [water_cluster_info, chemical_feature_clusters_info]

       coot.chemical_feature_clusters_accept_info_py(site_number, protein_res_specs,
                                                     imol_ligand_specs, cluster_info)

def cfc_process_sites(sites):

    # what is a site?  Good question!

    for site in sites:
       print("debug:: #### in cfc_process_sites() Here's a site:", site)

    # this adds a static empty non-null widget
    #
    coot.chemical_feature_clusters_setup_dialog()

    # does this (always) do what I want? (I want the site (list of specs) with the most members on top)
    #
    sorted_sites = sorted(sites)
    # this is debugging/testing/checking the above sort:
    for i,site in enumerate(sorted_sites):
       print("sorted site #{} site: {}".format(i, site))

    # imol_first = sites[0][0]
    imol_first = sites[0][0][0]

    print("debug:: in cfc_process_sites sites is ", sites)
    print("debug:: in cfc_process_sites sites[0] is ", sites[0])
    print("debug:: in cfc_process_sites sites[0][0] is ", sites[0][0])
    print("debug:: in cfc_process_sites imol_first is ", imol_first)

    # first_spec = sites[0][1]
    first_spec = sites[0][0][1]
    print("debug:: ########################## in cfc_proces_sites: first_spec:", first_spec)
    for i,site in enumerate(sorted_sites):
       print("   debug site-idx", i, "site:", site)

    site_number = 0
    cfc_process_site(site_number, sites, imol_first, first_spec)


class cfc_ligand_sites:

   def optimize_n(self, positions_np, n_data):

      bic = {}
      for n in [x+1 for x in range(5)]:
         if n < len(positions_np):
            # gmm = mixture.GMM(n_components=n, covariance_type='full', n_iter=20)
             gmm = mixture.GaussianMixture(n_components=n, covariance_type='full', max_iter=20)
             gmm.fit(positions_np)
             score = gmm.score(positions_np)
             # score = sum(score_list)
             lambda_c = 15 # 20 # 17 is enough (to result in 1 cluster) for 5en*
             bic_l = score - lambda_c * 0.5 * math.log(n_data) * n
             bic[n] = bic_l

      for key in bic:
          print("bic", key, bic[key])
      
      key, value = max(bic.iteritems(), key=lambda x:x[1])
      return key


   def merge_clusters(self, cluster_assignments, merge_map):

       continue_merging = True
       while continue_merging:
           key_longest_list = 'Unset'
           longest_list_len = 0;
           for key in merge_map:
               l = len(merge_map[key])
               if l > longest_list_len:
                   key_longest_list = key
                   longest_list_len = l
           if key_longest_list == 'Unset':
               continue_merging = False
           else:
               merged_something = False
               mergeable_list = merge_map[key_longest_list]
               for i in range(len(cluster_assignments)):
                   for m_c in mergeable_list:
                       cluster_idx = cluster_assignments[i]
                       if cluster_idx == m_c:
                           cluster_assignments[i] = key_longest_list
                           merged_something = True

               # remove key_longest_list and its data from map
               # print "del merge_map[", key_longest_list, "]"
               del merge_map[key_longest_list]
               if not(merged_something):
                   continue_merging = False
                             
       return cluster_assignments
  
    
   # cluster the ligands - i.e. give us the sites
   #
   def find_the_sites(self, file_name_comp_id_list):

      # main line
      #
      coords_with_spec = []

      for fn_comp_id in file_name_comp_id_list:
          fn      = fn_comp_id[0]
          comp_id = fn_comp_id[1]
          imol = coot.handle_read_draw_molecule_with_recentre(fn_comp_id[0], 0)
          # what are the residue specs for the given comp_ids?
          residue_specs = coot.get_residue_specs_in_mol_py(imol, comp_id)
          print(fn, residue_specs)

          for spec in residue_specs:
              # centre = residue_centre_from_spec_py(imol, spec)
              chain_id = rsu.residue_spec_to_chain_id(spec)
              res_no   = rsu.residue_spec_to_res_no(spec)
              ins_code = ''

              res_info = coot.residue_info_py(imol, chain_id, res_no, ins_code)

              for atom in res_info:
                  coords_with_spec.append([rsu.residue_atom_to_position(atom), imol, spec])

      # print coords_with_spec

      # now cluster coords. There will be 1 (usually), maybe 2 possibly 3 sites

      if len(coords_with_spec) < 3:

          return False

      else:

          coords = [x[0] for x in coords_with_spec]
          positions_np = np.array(coords)
          n_components = self.optimize_n(positions_np, len(positions_np))
          print( "optimize_n for sites::::::::::::", n_components)
          # dpgmm = mixture.GMM(n_components, covariance_type='full', n_iter=40)
          dpgmm = mixture.GaussianMixture(n_components, covariance_type='full', max_iter=40)
          dpgmm.fit(positions_np)

          cluster_assignments = dpgmm.predict(positions_np)
          means   = dpgmm.means_
          weights = dpgmm.weights_

          print(cluster_assignments)
          print(means)
          print(weights)

          print("cluster_assignments", cluster_assignments)

          merge_map = self.find_mergeable_clusters(means, weights)
          # which key (i.e. cluster index) has the most number of other clusters
          # that can be merged in?
          #
	  # convert to a list of ints (not <type 'numpy.int64'>) (because, on decoding Python->C++ object
	  # we do a PyInt_Check for the site_idx (and a <type 'numpy.int64'> fails that test)
	  #
          new_cluster_assignments = [int(x) for x in self.merge_clusters(cluster_assignments, merge_map)]
          print("new cluster_assignments", new_cluster_assignments)
 
          specs = [x[1:] for x in coords_with_spec]
          cluster_assignments_with_specs = zip(new_cluster_assignments, specs)

          sites = coot.chemical_feature_clusters_accept_site_clusters_info_py(cluster_assignments_with_specs)

          print("debug:: ########## in cfc_ligand_sites find_the_sites() chemical_feature_clusters_accept_site_clusters_info_py() returned", sites)
          for idx,s in enumerate(sites):
              print("site debug", idx, s)

          # show me them
          if True:  # debug
              o = coot.new_generic_object_number("site clusters")
              for mean in means:
                  cluster_star_obj(o, mean, 2, 2)
              # coot.set_display_generic_object(o, 1) this is for debugging

          self.sites = sites


   # if any site cluster is within 5A of any other site cluster, then
   # those site clusters are the same
   #
   def find_mergeable_clusters(self, means, weights):
       cluster_dist_crit = 5 # Bradley et al.
       merge_map = {}
       for i,mean_i in enumerate(means):
           for j, mean_j in enumerate(means):
               if j > i:
                   d = np.linalg.norm(mean_i-mean_j)
                   # print(mean_i - mean_j, d)
                   if (d < cluster_dist_crit):
                      try:
                         merge_map[i].append(j)
                      except KeyError as e:
                         merge_map[i] = [j]
                      try:
                         merge_map[j].append(i)
                      except KeyError as e:
                         merge_map[j] = [i]
       return merge_map
   

   def __init__(self, file_name_comp_id_list):
      self.sites = []
      self.find_the_sites(file_name_comp_id_list)


   def get_sites(self):
       return self.sites

