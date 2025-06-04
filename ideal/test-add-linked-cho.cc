
#include "clipper/core/xmap.h"
#include "clipper/core/map_utils.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats

#include "coot-utils/atom-selection-container.hh"
#include "add-linked-cho.hh"


// copied from mini-rsr. Consider making this a util function
std::pair<bool, clipper::Xmap<float> >
map_from_mtz(std::string mtz_file_name,
	     std::string f_col,
	     std::string phi_col,
	     std::string weight_col,
	     int use_weights,
	     bool is_debug_mode) {


   bool status = 0; // not filled

   clipper::HKL_info myhkl;
   clipper::MTZdataset myset;
   clipper::MTZcrystal myxtl;
   clipper::Xmap<float> xmap;

   try {
      std::cout << "reading mtz file..." << std::endl;
      clipper::CCP4MTZfile mtzin;
      mtzin.open_read( mtz_file_name );       // open new file
      mtzin.import_hkl_info( myhkl );         // read sg, cell, reso, hkls
      clipper::HKL_data< clipper::datatypes::F_sigF<float> >   f_sigf_data(myhkl, myxtl);
      clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data(myhkl, myxtl);
      clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphidata(myhkl, myxtl);

      std::string mol_name = mtz_file_name + " ";
      mol_name += f_col;
      mol_name += " ";
      mol_name += phi_col;

      if (use_weights) {
	 mol_name += " ";
	 mol_name += weight_col;
      }

      if (use_weights) {
	 clipper::String dataname = "/*/*/[" + f_col + " " + f_col + "]";
	 std::cout << dataname << "\n";
	 mtzin.import_hkl_data(  f_sigf_data, myset, myxtl, dataname );
	 dataname = "/*/*/[" + phi_col + " " + weight_col + "]";
	 std::cout << dataname << "\n";
	 mtzin.import_hkl_data( phi_fom_data, myset, myxtl, dataname );
	 mtzin.close_read();
	 std::cout << "We should use the weights: " << weight_col << std::endl;
	 // it seems to me that we should make 2 data types, an F_sigF and a phi fom
	 // and then combine them using a Convert_fsigf_phifom_to_fphi();

	 fphidata.compute(f_sigf_data, phi_fom_data,
			  clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

      } else {
	 clipper::String dataname = "/*/*/[" + f_col + " " + phi_col + "]";
	 mtzin.import_hkl_data(     fphidata, myset, myxtl, dataname );
	 mtzin.close_read();
      }
      std::cout << "Number of reflections: " << myhkl.num_reflections() << "\n";
      xmap.init( myhkl.spacegroup(), myhkl.cell(),
		 clipper::Grid_sampling( myhkl.spacegroup(),
					 myhkl.cell(),
					 myhkl.resolution()) );
      std::cout << "Grid..." << xmap.grid_sampling().format() << "\n";
      std::cout << "doing fft..." << std::endl;

      if (is_debug_mode) {
	 int count = 0;
	 clipper::HKL_info::HKL_reference_index hri;
	 for (hri=fphidata.first(); !hri.last(); hri.next()) {
	    if (count == 10)
	       break;
	    std::cout << "sample data " << " "
		      << hri.hkl().h() << " "
		      << hri.hkl().k() << " "
		      << hri.hkl().l() << " : "
		      << fphidata[hri].f() << " " << fphidata[hri].phi()*180/M_PI << std::endl;
	    count++;
	 }
      }
      xmap.fft_from(fphidata);                  // generate map
      std::cout << "done fft..." << std::endl;
      status = 1;
   }
   catch (const clipper::Message_base &exc) {  // "exception" is a protected word, it seems.
      std::cout << "Failed to read mtz file " << mtz_file_name << std::endl;
   }

   return std::pair<bool, clipper::Xmap<float> > (status, xmap);
}




int main(int argc, char **argv) {

   int  status = 0;

   std::string pdb_file_name = "2qc1-sans-cho.pdb";
   std::string mtz_file_name = "2qc1_map.mtz";
   std::string asn_chain_id = "B";
   int asn_res_no = 141;

   pdb_file_name = "pdb8zwp-sans-cho.pdb";
   mtz_file_name = "8zwp_map.mtz";
   asn_res_no = 174;

   int imol = 0;
   bool use_gemmi = false;
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, use_gemmi);
   if (asc.read_success) {

      std::pair<bool, clipper::Xmap<float> > xmap_pair =
         map_from_mtz(mtz_file_name, "FWT", "PHWT", "W", false, false);
      if (xmap_pair.first) {
         const clipper::Xmap<float> &xmap = xmap_pair.second;
         coot::protein_geometry geom;
         geom.init_standard();
         coot::cho::add_named_glyco_tree("NAG-NAG-BMA", &asc, imol, xmap, &geom, asn_chain_id, asn_res_no);
         asc.mol->WritePDBASCII("done.pdb");
      }
   }
   return status;
}
