
#include <iostream>
#include <chrono>
#include <thread>

#include "utils/coot-utils.hh"
#include "utils/split-indices.hh"
#include "contacts-by-bricks.hh"

coot::contacts_by_bricks::contacts_by_bricks(mmdb::PAtom *atoms_in, int n_atoms_in, const std::set<unsigned int> &fixed_atom_indices) {

   atoms = atoms_in;
   n_atoms = n_atoms_in;
   brick_size = 20.0;
   dist_nbc_max = 8.0;
   only_between_different_residues_flag = false;

   range[0] = 0; range[1] = 0; range[2] = 0;
   set_lower_left_and_range(atoms_in, n_atoms_in);

   int n_bricks = range[0] * range[1] * range[2];
   atoms_in_bricks.resize(n_bricks);

   fill_the_bricks();

   fixed_flags.resize(n_atoms, false);
   std::set<unsigned int>::const_iterator it;
   for(it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
      fixed_flags[*it] = true;

   unsigned int n_threads = get_max_number_of_threads();
   unsigned int n_thread_sets = n_threads -1;
   if (n_thread_sets < 1)
      n_thread_sets = 1;
   split_indices(&thread_index_sets, n_bricks, n_thread_sets);

}

void
coot::contacts_by_bricks::set_dist_max(float f) {

   dist_nbc_max = f;
}

void
coot::contacts_by_bricks::fill_the_bricks() {

   float inv_brick_size = 1.0/brick_size;
   for(int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atoms[i];
      int idx_3d[3];
      // beware when copying this later - when atoms move?
      idx_3d[0] = static_cast<int> ((at->x - lower_left[0]) * inv_brick_size);
      idx_3d[1] = static_cast<int> ((at->y - lower_left[1]) * inv_brick_size);
      idx_3d[2] = static_cast<int> ((at->z - lower_left[2]) * inv_brick_size);
      unsigned int idx_1d = idx_3d_to_idx_1d(idx_3d);
      // atoms that fly over the edge don't have NBCs :-)
      unsigned int n_bricks = atoms_in_bricks.size();
      if (idx_1d < n_bricks) {
	 std::set<unsigned int> &ss = atoms_in_bricks.at(idx_1d);
	 // ss.contains()?
	 if (ss.find(i) == ss.end()) {
	    ss.insert(i);
	    // now delete i from the neighbour brick sets
	    // ..
	    if (true) { // this block ~5ms for 400 residues
	       for (int ix=idx_3d[0]-1; ix<=idx_3d[0]+1; ix++) {
		  if (ix >= 0) {
		     if (ix < range[0]) {
			for (int iy=idx_3d[1]-1; iy<=idx_3d[1]+1; iy++) {
			   if (iy >= 0) {
			      if (iy < range[1]) {
				 for (int iz=idx_3d[2]-1; iz<=idx_3d[2]+1; iz++) {
				    if (iz >= 0) {
				       if (iz < range[2]) {
					  int idx_3d_neighb[3];
					  if (! ((ix==0) && (iy==0) && (iz==0))) {
					     idx_3d_neighb[0] = ix;
					     idx_3d_neighb[1] = iy;
					     idx_3d_neighb[2] = iz;
					     unsigned int idx_neighb(idx_1d + idx_3d_to_idx_1d(idx_3d_neighb));
					     if (idx_neighb < n_bricks)
						atoms_in_bricks.at(idx_neighb).erase(i);
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::contacts_by_bricks::set_lower_left_and_range(mmdb::PAtom *atoms_in, int n_atoms_in) {

   for (int i=0; i<3; i++)
      lower_left[i] = 9999.9;
   //
   unsigned int n_atoms = static_cast<unsigned int>(n_atoms_in);
   for(unsigned int i=0; i<n_atoms; i++) {
      mmdb::Atom *atom = atoms_in[i];
      float pos[3];
      pos[0] = atom->x; pos[1] = atom->y; pos[2] = atom->z;
      for (int j=0; j<3; j++)
	 if (pos[j] < lower_left[j])
	    lower_left[j] = pos[j];
   }

   for (int i=0; i<3; i++)
      lower_left[i] -= 6.0; // say, maybe more needed

   float inv_brick_size = 1.0/brick_size;
   for(unsigned int i=0; i<n_atoms; i++) {
      mmdb::Atom *atom = atoms_in[i];
      float pos[3];
      pos[0] = atom->x; pos[1] = atom->y; pos[2] = atom->z;
      for (int j=0; j<3; j++) {
	 if (false)
	    std::cout
	       << pos[0] << " " << lower_left[0] << " "
	       << pos[1] << " " << lower_left[1] << " "
	       << pos[2] << " " << lower_left[2] << " "
	       << std::endl;
	 float delta = (pos[j] - lower_left[j]);
	 float f = (pos[j] - lower_left[j]) * inv_brick_size;
	 int brick_idx = static_cast<int>(f);
	 if (brick_idx > range[j]) {
	    range[j] = brick_idx;
	 }
      }
   }

   for (int i=0; i<3; i++)
      range[i] += 1;  // maybe more needed

   // std::cout << "ranges: " << range[0] << " " << range[1] << " " << range[2] << std::endl;
}

unsigned int
coot::contacts_by_bricks::idx_3d_to_idx_1d(int idx_3d[3]) const {

   unsigned int idx = range[0] * range[1] * idx_3d[2] + range[1] * idx_3d[1] + idx_3d[0];
   return idx;
}


void
coot::contacts_by_bricks::find_the_contacts(std::vector<std::set<unsigned int> > *vec_p,
					    bool only_between_different_residues_flag) {

   vec_p->resize(n_atoms); // can't reserve std::sets
   fill_the_bricks();
   find_the_contacts_in_bricks(vec_p, only_between_different_residues_flag);
   find_the_contacts_between_bricks(vec_p, only_between_different_residues_flag);
}

void
coot::contacts_by_bricks::find_the_contacts_in_bricks(std::vector<std::set<unsigned int> > *vec,
						      bool only_between_different_residues_flag) const {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   unsigned int n_in_brick = 0;
   int n_bricks = atoms_in_bricks.size();
   float dist_max_sqrd = dist_nbc_max * dist_nbc_max;
   for (int ib=0; ib<n_bricks; ib++) {
      const std::set<unsigned int> &brick_base = atoms_in_bricks[ib];
      std::set<unsigned int>::const_iterator it_base;
      std::set<unsigned int>::const_iterator it_neighb;
      for (it_base=brick_base.begin(); it_base!=brick_base.end(); it_base++) {
	 if (!fixed_flags[*it_base]) {
	    mmdb::Atom *at_1 = atoms[*it_base];
	    for (it_neighb=brick_base.begin(); it_neighb!=brick_base.end(); it_neighb++) {
	       if (it_neighb != it_base) {
		  mmdb::Atom *at_2 = atoms[*it_neighb];
		  if (only_between_different_residues_flag)
		     if (at_2->residue == at_1->residue)
			continue;
		  float d_x(at_1->x - at_2->x);
		  float d_y(at_1->y - at_2->y);
		  float d_z(at_1->z - at_2->z);
		  float dd = d_x * d_x + d_y * d_y + d_z * d_z;
		  if (dd < dist_max_sqrd) {
		     vec->at(*it_base).insert(*it_neighb);
		     n_in_brick++;
		  }
	       }
	    }
	 }
      }
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   // std::cout << "------- contacts_by_bricks(): in bricks: " << d10 << " milliseconds " << std::endl;

   // std::cout << "Found n_in_brick " << n_in_brick << std::endl;
}

void
coot::contacts_by_bricks::find_the_contacts_between_bricks(std::vector<std::set<unsigned int> > *vec_p,
							   bool only_between_different_residues_flag) const {

   // find_the_contacts_between_bricks_simple(vec_p);

   find_the_contacts_between_bricks_multi_thread(vec_p, only_between_different_residues_flag);

   if (false) {
      for (std::size_t ii=0; ii<vec_p->size(); ii++) {
	 const std::set<unsigned int> &ss = vec_p->at(ii);
	 std::cout << "Atom " << ii << " : ";
	 std::set<unsigned int>::const_iterator it;
	 for (it=ss.begin(); it!=ss.end(); it++) {
	    std::cout << *it << " ";
	 }
	 std::cout << "\n";
      }
   }

}

void
coot::contacts_by_bricks::find_the_contacts_between_bricks_multi_thread(std::vector<std::set<unsigned int> > *vec_p,
									bool only_between_different_residues_flag) const {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   unsigned int n_btwn_bricks = 0;

   // ~20 microseconds for 400 residues
   // auto tp_s1 = std::chrono::high_resolution_clock::now();
   // auto ds10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_s1 - tp_0).count();
   // std::cout << "------- brick set splitting takes : " << ds10 << " microseconds " << std::endl;

   int brick_index_max = range[0] * range[1] * range[2];
   float dist_max_sqrd = dist_nbc_max * dist_nbc_max;

   std::vector<std::thread> threads;

   for (std::size_t ii=0; ii<thread_index_sets.size(); ii++) {
      const std::vector<unsigned int> &index_set = thread_index_sets[ii];
      threads.push_back(std::thread(find_the_contacts_between_bricks_multi_thread_workpackage,
      				    vec_p, std::cref(index_set), std::cref(atoms_in_bricks),
				    std::cref(fixed_flags), range, atoms, brick_index_max, dist_nbc_max,
				    only_between_different_residues_flag));
   }
   for (std::size_t ii=0; ii<thread_index_sets.size(); ii++)
      threads[ii].join();

   if (false) {
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "------- contacts_by_bricks(): between_brick_multi: " << d10 << " milliseconds " << std::endl;
   }

}

// the function for the thread:
//
// We need to be sure that there is only one index of each atom across all of the bricks, otherwise race
void
coot::contacts_by_bricks::find_the_contacts_between_bricks_multi_thread_workpackage(std::vector<std::set<unsigned int> > *vec,
										    const std::vector<unsigned int> &index_set,
										    const std::vector<std::set<unsigned int> > &atoms_in_bricks,
										    const std::vector<bool> &fixed_flags,
										    const int brick_range[3],
										    mmdb::PAtom *atoms,
										    int brick_index_max,
										    float dist_max,
										    bool only_between_different_residues_flag) {

   float dist_max_sqrd = dist_max * dist_max;
   for (std::size_t ii=0; ii<index_set.size(); ii++) {

      int ib = index_set[ii];
      const std::set<unsigned int> &brick_base = atoms_in_bricks[ib];

      if (brick_base.size() > 0) {
	 for (int iz=-1; iz<2; iz++) { // or do I mean ix?
	    for (int iy= -1; iy<2; iy++) {
	       for (int ix= -1; ix<2; ix++) {
		  int ib_neighb = ib + ix + iy * brick_range[0] + iz * brick_range[0] * brick_range[1];
		  if ((ib_neighb >= 0) && (ib_neighb != ib)) {
		     if (ib_neighb < brick_index_max) {
			const std::set<unsigned int> &brick_neighb = atoms_in_bricks[ib_neighb];
			std::set<unsigned int>::const_iterator it_base;
			std::set<unsigned int>::const_iterator it_neighb;
			for (it_base=brick_base.begin(); it_base!=brick_base.end(); it_base++) {
			   if (!fixed_flags[*it_base]) {
			      mmdb::Atom *at_1 = atoms[*it_base];
			      for (it_neighb=brick_neighb.begin(); it_neighb!=brick_neighb.end(); it_neighb++) {
				 mmdb::Atom *at_2 = atoms[*it_neighb];
				 if (only_between_different_residues_flag)
				    if (at_2->residue == at_1->residue)
				       continue;
				 float d_x(at_1->x - at_2->x);
				 float d_y(at_1->y - at_2->y);
				 float d_z(at_1->z - at_2->z);
				 float dd(d_x * d_x + d_y * d_y + d_z * d_z);
				 // std::cout << "MP " << *it_base << " " << *it_neighb << " sqrt(dd) " << sqrt(dd) << std::endl;
				 if (dd < dist_max_sqrd) {
				    // If this is not the first time around, it's probably already there.
				    // This is not a "go-faster" test, it's here to cut down the output.
				    if (vec->at(*it_base).find(*it_neighb) == vec->at(*it_base).end()) {
				       vec->at(*it_base).insert(*it_neighb);
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}


void
coot::contacts_by_bricks::find_the_contacts_between_bricks_simple(std::vector<std::set<unsigned int> > *vec,
								  bool only_between_different_residues_flag) const {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   unsigned int n_btwn_bricks = 0;

   int n_bricks = atoms_in_bricks.size();

   // std::vector<std::vector<unsigned int> > index_sets;
   // split_indices(&index_sets, n_bricks, n_sets);

   // ~20 microseconds for 400 residues
   // auto tp_s1 = std::chrono::high_resolution_clock::now();
   // auto ds10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_s1 - tp_0).count();
   // std::cout << "------- brick set splitting takes : " << ds10 << " microseconds " << std::endl;

   int brick_index_max = range[0] * range[1] * range[2];
   float dist_max_sqrd = dist_nbc_max * dist_nbc_max;

   for (int ib=0; ib<n_bricks; ib++) {
      const std::set<unsigned int> &brick_base = atoms_in_bricks[ib];
      if (brick_base.size() > 0) {
	 for (int iz=-1; iz<2; iz++) { // or do I mean ix?
	    for (int iy= -1; iy<2; iy++) {
	       for (int ix= -1; ix<2; ix++) {
		  int ib_neighb = ib + ix + iy * range[0] + iz * range[0] * range[1];
		  if ((ib_neighb >= 0) && (ib_neighb != ib)) {
		     if (ib_neighb < brick_index_max) {
			const std::set<unsigned int> &brick_neighb = atoms_in_bricks[ib_neighb];
			std::set<unsigned int>::const_iterator it_base;
			std::set<unsigned int>::const_iterator it_neighb;
			for (it_base=brick_base.begin(); it_base!=brick_base.end(); it_base++) {
			   if (!fixed_flags[*it_base]) {
			      mmdb::Atom *at_1 = atoms[*it_base];
			      for (it_neighb=brick_neighb.begin(); it_neighb!=brick_neighb.end(); it_neighb++) {
				 mmdb::Atom *at_2 = atoms[*it_neighb];
				 if (only_between_different_residues_flag)
				    if (at_2->residue == at_1->residue)
				       continue;
				 float d_x(at_1->x - at_2->x);
				 float d_y(at_1->y - at_2->y);
				 float d_z(at_1->z - at_2->z);
				 float dd(d_x * d_x + d_y * d_y + d_z * d_z);
				 if (dd < dist_max_sqrd) {
				    vec->at(*it_base).insert(*it_neighb);
				    n_btwn_bricks++;
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   // std::cout << "------- between bricks: " << d10 << " milliseconds " << std::endl;

   // std::cout << "Found n_btwn_bricks " << n_btwn_bricks << std::endl;
}
