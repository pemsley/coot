
#include "zo-rama.hh"

int main(int argc, char **argv) {

   zo::rama_table_set rts;
   std::map<std::string, zo::rama_table>::const_iterator it;
   for (it=rts.table_map.begin(); it!=rts.table_map.end(); it++) {
      const zo::rama_table &rt = it->second;
      std::string zo_residue_type = it->first;
      std::cout << zo_residue_type << " " << rt.rama_vec.size() << std::endl;
      std::string file_name = zo_residue_type + ".png";
      rt.make_a_png(600, file_name);
   }
   return 0;

}
