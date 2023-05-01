

void gensurf_and_add_vecs_threaded_workpackage(const clipper::Xmap<float> *xmap_p,
					       float contour_level, float dy_radius,
					       coot::Cartesian centre,
					       int iream_start, int n_reams,
					       int isample_step,
					       bool is_em_map,
					       std::vector<coot::density_contour_triangles_container_t> *draw_vector_sets_p);

// std::vector<std::pair<const coot::CartesianPair *,int> > draw_vector_sets;
