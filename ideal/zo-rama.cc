
#include "zo-rama.hh"


void
zo::rama_table::test_analytical_derivs() const {

   int width  = 10;
   int height = 10;
   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 zo::realtype phi = M_PI *  2.0 * zo::realtype(i-0.5*height)/zo::realtype(height);  // -pi to +pi
	 zo::realtype psi = M_PI * -2.0 * zo::realtype(j-0.5*width )/zo::realtype(width);
	 zo::realtype sum = 0;
	 std::pair<zo::realtype, zo::realtype> df_n = df_numerical(phi, psi);
	 std::pair<zo::realtype, zo::realtype> df_a = df          (phi, psi);
      }
   }
}

void
zo::rama_table::make_a_png(int width, const std::string &file_name) {

   
   int height = width; // this is so for a ramachandran plot

   png_bytep *row_pointers = static_cast<png_bytep *>(malloc(sizeof(png_bytep) * height));
   for (int i=0; i<height; i++)
      row_pointers[i] = static_cast<png_byte*>(malloc(width));
   std::vector<std::vector<zo::realtype> > v(height);
   for (int i=0; i<height; i++)
      v[i] = std::vector<zo::realtype>(width, 0);

   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 zo::realtype phi = M_PI *  2.0 * double(i-0.5*height)/double(height);
	 zo::realtype psi = M_PI * -2.0 * double(j-0.5*width )/double(width); // -pi to +pi
	 zo::realtype sum = value(phi,psi);
	 zo::realtype d = expf(sum);
	 v[j][i] = d;
      }
   }

   // replace by gradient for testing
   if (false) {
      for (int j=0; j<width; j++) {
	 for (int i=0; i<height; i++) {
	    zo::realtype phi = M_PI *  2.0 * double(i-0.5*height)/double(height);
	    zo::realtype psi = M_PI * -2.0 * double(j-0.5*width )/double(width); // -pi to +pi
	    std::pair<zo::realtype, zo::realtype> grads = df_numerical(phi,psi);
	    // zo::realtype d = expf(sqrt(grads.first*grads.first + grads.second*grads.second));
	    zo::realtype d = sqrt(grads.first*grads.first + grads.second*grads.second);
	    v[j][i] = d;
	 }
      }
   }

   zo::realtype min_v = 9.99e12;
   zo::realtype max_v = 0.0;
   zo::realtype sum = 0;
   int n = 0;
   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 sum += v[j][i];
	 n++;
	 if (v[j][i] < min_v) min_v = v[j][i];
	 if (v[j][i] > max_v) max_v = v[j][i];
      }
   }
   zo::realtype mean = sum/zo::realtype(n);
   zo::realtype sf = 0.1 / mean;

   for (int j=0; j<width; j++) {
      for (int i=0; i<height; i++) {
	 zo::realtype vv = v[j][i] * sf;
	 vv = 255 * (1 - vv);
	 if (vv < 0) vv = 0;
	 // int pixel_value = std::lround(vv);
	 int pixel_value = int(vv+0.5);
	 row_pointers[j][i] = pixel_value;
      }
   }

   write_png_file(width, height, row_pointers, file_name.c_str());

   // now clean up
   for (int y=0; y<height; y++)
      free(row_pointers[y]);
   free(row_pointers);

}

