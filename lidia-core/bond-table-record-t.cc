
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

// some-header

#include <iostream>
#include <iomanip>

#include "bond-table-record-t.hh"

void
cod::bond_table_record_t::write(std::ostream &s) const { 

   s << std::setw(10) << mean;
   s << std::setw(10) << std_dev;
   s << std::setw(6) << count;

   std::string::size_type s1 = cod_type_1.level_4.length();
   std::string::size_type s2 = cod_type_2.level_4.length();

   s << std::setw(4) << s1;
   s << std::setw(4) << s2;

   s << " ";
   s << cod_type_1.level_4;
   s << " ";
   s << cod_type_2.level_4;
   s << "\n";
}

void
cod::bond_table_record_t::write(std::ostream &s,
				unsigned int type_index_1,
				unsigned int type_index_2) const {

      s << std::setw(10) << mean;
      s << std::setw(10) << std_dev;
      s << std::setw(6) << count;

      s << " ";
      
      s << std::setw(7) << type_index_1 << " ";
      s << std::setw(7) << type_index_2 <<"\n";
}


std::ostream &
cod::operator<<(std::ostream &s, const cod::bond_table_record_t &btr) {

   s << btr.cod_type_1.level_4 << " " << btr.cod_type_2.level_4 << " " << btr.mean
     << " " << btr.std_dev << " " << btr.count;
   return s;
}


#endif // MAKE_ENHANCED_LIGAND_TOOLS
