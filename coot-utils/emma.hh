
namespace coot { 
   namespace util {
      
      class emma {
	 void sfs_from_boxed_molecule(CMMDBManager *mol, float border);
	 double f(double r) const;
      public:
	 emma(CMMDBManager *mol, float border) {
	    sfs_from_boxed_molecule(mol, border);
	 } 
	 clipper::Spacegroup spacegroup;
	 clipper::Cell cell;
	 clipper::Resolution reso;
	 clipper::HKL_info hkl_info;
	 clipper::HKL_data<clipper::data32::F_phi> fc;
	 void integrate(const clipper::Xmap<float> &xmap) const;
	 void test() const;
      };
   }
}
  
