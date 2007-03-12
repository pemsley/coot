
#ifndef XMAP_UTILS_H
#define XMAP_UTILS_H

template<class T>
class mean_and_variance { 

 public:
  T mean; 
  T variance; 
  T range; 
  T max_density;
  T min_density;

  // Fix this in a spare moment:  To complex (or I am too stupid
  // to do it while rushing now). 
  // 
  // friend ostream& operator<<(ostream&, mean_and_variance<T>);


};

  
template<class T>
mean_and_variance<T>
map_density_distribution(const clipper::Xmap<T> &map, short int write_output_flag); 

// I suppose that this could be templated to (i.e. an int is also useful)
// 
float nearest_step(float val, float step); 

int 
nint(float val);


class Xmap_plus_bits {
 public:
   clipper::Xmap<float> xmap;
   std::string name;
   float min, max, mean, std_dev; 
}; 


#endif // XMAP_UTILS_H
