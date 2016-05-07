
#ifndef BOND_TABLE_RECORD_T_HH
#define BOND_TABLE_RECORD_T_HH

#include <string>

#include "cod-atom-type-t.hh"

namespace cod {
   
   class bond_table_record_t {
   public:
      // std::string cod_type_1;
      // std::string cod_type_2;
      atom_type_t cod_type_1;
      atom_type_t cod_type_2;
      double mean;
      double std_dev;
      unsigned int count;
      unsigned int approx_level;

      bond_table_record_t() {}
      bond_table_record_t(const atom_type_t &cod_type_1_in,
			  const atom_type_t &cod_type_2_in,
			  const double &mean_in,
			  const double &std_dev_in,
			  unsigned int count_in,
			  unsigned int al=0) {
	 cod_type_1 = cod_type_1_in;
	 cod_type_2 = cod_type_2_in;
	 mean = mean_in;
	 std_dev = std_dev_in;
	 count = count_in;
	 approx_level = al;
      }
      bool operator<(const bond_table_record_t &btr) const {

	 if (cod_type_1 == btr.cod_type_1) {
	    return (cod_type_2 < btr.cod_type_2);
	 } else {
	    return (cod_type_1 < btr.cod_type_1);
	 }
      }
      void write(std::ostream &s) const;
      void write(std::ostream &s,
		 unsigned int type_index_1,
		 unsigned int type_index_2) const;
      friend std::ostream &operator<<(std::ostream &s, const bond_table_record_t &btr);
   };
   std::ostream &operator<<(std::ostream &s, const bond_table_record_t &brc);

}

#endif // BOND_TABLE-RECORD_T_HH
