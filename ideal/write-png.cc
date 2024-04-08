/*
 * idea/write-png.cc
 *
 * Copyright 2017 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


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
