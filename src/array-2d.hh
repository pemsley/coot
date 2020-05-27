
#ifndef ARRAY_2D_HH
#define ARRAY_2D_HH


class array_2d {
public:
   std::vector<float> v;
   unsigned int nx;
   array_2d(unsigned int nx, unsigned int ny) : nx(nx) {
      v.resize(nx * ny);
   }
   void set(const unsigned int &i, const unsigned int &j, const float &f) {
      v[j * nx + i] = f;
   }
   float get(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i];
   }

};

class coord_array_2d {
public:
   std::vector<std::pair<clipper::Coord_orth, float> > v;
   unsigned int nx;
   coord_array_2d(unsigned int nx, unsigned int ny) : nx(nx) {
      v.resize(nx * ny);
   }
   void set(const unsigned int &i, const unsigned int &j, const clipper::Coord_orth &co, const float &f) {
      v[j * nx + i] = std::pair<clipper::Coord_orth, float> (co, f);
   }
   std::pair<clipper::Coord_orth, float>  get(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i];
   }
   clipper::Coord_orth get_co(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i].first;
   }
   float get_f(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i].second;
   }

};

#endif // ARRAY_2D

