
#ifndef COOT_SRC_INTERMEDIATE_ATOM_HH
#define COOT_SRC_INTERMEDIATE_ATOM_HH

#include <mmdb2/mmdb_manager.h>
#include "coords/Cartesian.hh"

class intermediate_atom_distance_t {
   coot::Cartesian static_position;
   mmdb::Atom *dynamic_atom;
   bool static_pos_filled_flag;

public:
   intermediate_atom_distance_t() {
      dynamic_atom = 0;
      static_pos_filled_flag = 0;
   }
   explicit intermediate_atom_distance_t(const coot::Cartesian &pt) : static_position(pt) {
      dynamic_atom = 0;
      static_pos_filled_flag = 1;
   }
   explicit intermediate_atom_distance_t(mmdb::Atom *at) : dynamic_atom(at) {
      static_pos_filled_flag = 0;
   }
   void draw_dynamic_distance() const;
   bool static_position_is_filled() const { return static_pos_filled_flag; }
   bool atom_is_filled() const {
      return (dynamic_atom != 0);
   }
   void add_atom(mmdb::Atom *at) {
      dynamic_atom = at;
   }
   void add_static_point(coot::Cartesian &pt) {
      static_position = pt;
      static_pos_filled_flag = 1;
   }
   bool filled() const {
      return (static_pos_filled_flag && dynamic_atom);
   }
};


#endif // COOT_SRC_INTERMEDIATE_ATOM_HH
