/*
 * src/molecule-class-info-string.cc
 *
 * Copyright 2020 by Medical Research Council
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


// Eugene Krissinel wrote this.
//
// I just put a coot molecule wrapper around it

#include "molecule-class-info.h"


std::string
molecule_class_info_t::pdb_string() const {

   std::string s;

   mmdb::io::PFile f = new mmdb::io::File();   // file class instance
   mmdb::pstr buf = new char[100001];  // initial buffer for writing

   f->assign ( 100000,  // initial buffer size, not the difference with allocation_size
               10000,   // if buffer is not sufficient, this is the increment_size
               buf      // initial buffer
               );

   /*  I think that this would also work:
       f->assign ( 0,  // initial buffer size, not the difference with allocation size
       10000,          //  buffer increment size
       null            // no initial buffer
       );
   */

   f->rewrite();   // opens class for writing

   // then use normal MMDB write function, e.g.

   // atom_sel.mol->WritePDBASCII( &(*f) ); // the argument should cast to &File

   atom_sel.mol->WritePDBASCII(*f); // the argument should cast to &File

   // now extract the string with content: note that this may differ from
   // 'buf' due to memory reallocations:

   mmdb::pstr stringFile;
   mmdb::word stringFileLength;
   f->takeFilePool ( stringFile,        // should come 0-terminated if text I/O
                     stringFileLength   // useful in case of binary I/O
                     );

   s = std::string(stringFile, stringFileLength);
   f->shut();  // 'close' file

   delete [] buf;
   delete f;

   return s;

}
