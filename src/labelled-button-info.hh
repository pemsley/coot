#ifndef LABELLED_BUTTON_INFO_HH
#define LABELLED_BUTTON_INFO_HH

#include <clipper/core/coords.h>

class labelled_button_info_t {
public:
   std::string label;
   clipper::Coord_orth position;
   labelled_button_info_t(const std::string &l, const clipper::Coord_orth &co) : label(l), position(co) {}
};

#endif // LABELLED_BUTTON_INFO_HH
