
#ifndef COOT_UTILS_G_TRIANGLE_HH
#define COOT_UTILS_G_TRIANGLE_HH

#include <ostream>

//! g_triangle is the container for the vertices of a triangles (indexing into the vertices vector).
class g_triangle {
public:
   //! constructor
   g_triangle(const unsigned int &a0,
              const unsigned int &a1,
              const unsigned int &a2) {
      point_id[0] = a0;
      point_id[1] = a1;
      point_id[2] = a2;
      // colour_index = -1;
   }
   g_triangle() {} // for resize
   unsigned int point_id[3];
   //! int colour_index;
   unsigned int &operator[] (const unsigned int &i) { return point_id[i]; }
   const unsigned int &operator[] (const unsigned int &i) const { return point_id[i]; }
   //! use ``rebase()`` when adding more vertices and triangles into a mesh.
   void rebase(const unsigned int &idx_base) {
      for (unsigned int i=0; i<3; i++) {
         point_id[i] += idx_base;
      }
   }
   //! in-place reverse the winding of this triangle
   void reverse_winding() {
      std::swap(point_id[0], point_id[1]);
   }
   friend std::ostream& operator <<(std::ostream &s, const g_triangle &t);
};

std::ostream& operator <<(std::ostream &s, const g_triangle &t);


// This is the g_triangle class that should be used for blender - and
// recently used in webassembly.
//! a triangle with a colour index - helpful, I think when transfering
//! colour/material information to Blender.
class g_triangle_with_colour_index : public g_triangle {
   public:
   int colour_index;
   g_triangle_with_colour_index(const unsigned int &a0,
                                const unsigned int &a1,
                                const unsigned int &a2) : g_triangle(a0, a1, a2) {
      colour_index = -1;
   }
};

#endif // G_TRIANGLE_HH

