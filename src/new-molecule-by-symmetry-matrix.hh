/*  ----------------------------------------------------------------------- */
/*                  Abstraction of New molecule by symmetry functions       */
/*  ----------------------------------------------------------------------- */

#include <mmdb2/mmdb_manager.h>

mmdb::Manager *new_molecule_by_symmetry_matrix_from_molecule(mmdb::Manager *mol,
                                                            double m11, double m12, double m13,
                                                            double m21, double m22, double m23,
                                                            double m31, double m32, double m33,
                                                            double tx, double ty, double tz,
                                                            int pre_shift_to_origin_na,
                                                            int pre_shift_to_origin_nb,
                                                            int pre_shift_to_origin_nc);
