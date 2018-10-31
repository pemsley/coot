

#ifndef DENSITY_RESULTS_CONTAINER_T
#define DENSITY_RESULTS_CONTAINER_T

class density_results_t {
   public:

   // angle in degrees
   density_results_t(const clipper::Coord_orth &p, double angle_in, float density_value_in) :
         position(p), angle(angle_in), density_value(density_value_in)  {}
   clipper::Coord_orth position;
   float angle; // degrees
   float density_value;
};

class density_results_container_t {

   public:
     std::vector<density_results_t> scored_points;

};
#endif // DENSITY_RESULTS_CONTAINER_T
