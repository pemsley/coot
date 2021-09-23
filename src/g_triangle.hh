
#ifndef G_TRIANGLE_HH
#define G_TRIANGLE_HH

#include <ostream>

class g_triangle {
public:
   g_triangle(const unsigned int &a0,
              const unsigned int &a1,
              const unsigned int &a2) {
      point_id[0] = a0;
      point_id[1] = a1;
      point_id[2] = a2;
   }
   g_triangle() {} // for resize
   unsigned int point_id[3];
   unsigned int &operator[] (const unsigned int &i) { return point_id[i]; }
   const unsigned int &operator[] (const unsigned int &i) const { return point_id[i]; }
   void rebase(const unsigned int &idx_base) {
      for (unsigned int i=0; i<3; i++) {
         point_id[i] += idx_base;
      }
   }
   friend std::ostream& operator <<(std::ostream &s, const g_triangle &t);
};
std::ostream& operator <<(std::ostream &s, const g_triangle &t);


#endif // G_TRIANGLE_HH

