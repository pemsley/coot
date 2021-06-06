
#ifndef COOT_MAP_INDEX_HH
#define COOT_MAP_INDEX_HH

namespace coot {
   class map_index_t {
      int index_;
   public:
      enum index_type { UNASSIGNED = -1 };
      map_index_t() { index_ = UNASSIGNED; }
      explicit map_index_t(int i) { index_ = i; }
      int index() const { return index_; }
      bool is_assigned() const { return (index_ != UNASSIGNED); }
      bool operator==(const map_index_t &ti) const {
	 return (ti.index() == index_);
      }
   };

}


#endif // MAP_INDEX_HH
