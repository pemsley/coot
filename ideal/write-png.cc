
#include <iostream>
#include "zo-rama.hh"

void
zo::write_png_file(int width, int height, png_bytep *row_pointers,
		   const std::string &file_name) {

   /* create file */

   png_byte bit_depth = 8;
   png_byte color_type = PNG_COLOR_TYPE_GRAY; // or PNG_COLOR_TYPE_GRAY_ALPHA

   FILE *fp = fopen(file_name.c_str(), "wb");
   if (!fp) {
      std::cout << "[write_png_file] File " << file_name
		<< " could not be opened for writing." << std::endl;
   } else {

      /* initialize stuff */
      png_structp png_ptr;
      png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

      if (!png_ptr) {
	 std::cout << "[write_png_file] png_create_write_struct failed\n";
	 return;
      }

      png_infop info_ptr;
      info_ptr = png_create_info_struct(png_ptr);
      if (!info_ptr) {
	 std::cout << "[write_png_file] png_create_info_struct failed\n";
	 return;
      }

      if (setjmp(png_jmpbuf(png_ptr))) {
	 std::cout << "[write_png_file] Error during init_io\n";
	 return;
      }

      png_init_io(png_ptr, fp);


      /* write header */
      if (setjmp(png_jmpbuf(png_ptr))) {
	 std::cout << "[write_png_file] Error during writing header" << std::endl;
	 return;
      }

      png_set_IHDR(png_ptr, info_ptr, width, height,
		   bit_depth, color_type, PNG_INTERLACE_NONE,
		   PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

      png_write_info(png_ptr, info_ptr);


      /* write bytes */
      if (setjmp(png_jmpbuf(png_ptr))) {
	 std::cout << "[write_png_file] Error during writing bytes\n";
	 return;
      }

      png_write_image(png_ptr, row_pointers);


      /* end write */
      if (setjmp(png_jmpbuf(png_ptr))) {
	 std::cout << "[write_png_file] Error during end of write\n";
	 return;
      }

      png_write_end(png_ptr, NULL);

      fclose(fp);
   }
}
