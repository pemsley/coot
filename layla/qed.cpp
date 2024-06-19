#include "qed.hpp"

namespace coot::layla::RDKit {

const QEDproperties QED::WEIGHT_MAX = QEDproperties({0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00});
const QEDproperties QED::WEIGHT_MEAN = QEDproperties({0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95});
const QEDproperties QED::WEIGHT_NONE = QEDproperties({1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00});
// RDKit::MolFromSmarts()

}