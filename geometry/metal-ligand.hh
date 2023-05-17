
class coordinated_atom_mode_t {
public:
   float mode;
   float std;
};

class coordinated_atom_t {
public:
   int coordination_number;
   float median;
   float mad;
   float mean;
   float std;
   unsigned int count;
   std::vector<coordinated_atom_mode_t> modes;
   coordinated_atom_t(int cn, float m1, float m2, float m3, float s, int c) :
      coordination_number(cn), median(m1), mad(m2), mean(m3), std(s), count(c) {}
};

class metal_ligand_t {
public:
   std::string element;
   std::vector<coordinated_atom_t> coordinated_atoms;
   explicit metal_ligand_t(const std::string &e) : element(e) {}
   void add(const coordinated_atom_t &ca) { coordinated_atoms.push_back(ca); }
};

