

class data_pair_remover {
   public:
   double r_lim;
   data_pair_remover(const double &r_in) { r_lim = r_in; }
   bool operator()(const std::pair<double,double> &rp) const {
      double r1 = coot::util::random()/float(RAND_MAX);
      // std::cout << "   " << r1 << " " << r_lim << "\n";
      return (r1 > r_lim);
   }
};
