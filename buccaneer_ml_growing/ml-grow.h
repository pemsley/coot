#pragma once

#include <vector>
#include "buccaneer-prot.h"
#include "k2c_tensor_include.h"

namespace ml_grow
{
    Ca_group next_ca_group(const Ca_group &ca_group, const clipper::Xmap<float> &xmap);
    Ca_group prev_ca_group(const Ca_group &ca_group, const clipper::Xmap<float> &xmap);
    void _zscore(k2c_tensor &tensor);
}
