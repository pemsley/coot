// Additional C++ Doxygen Documentation for cc-interface.hh
// These are enhanced documentation blocks for commonly-used functions
// that currently have minimal or no documentation.
//
// To use: Copy these documentation blocks into the appropriate locations
// in cc-interface.hh, replacing or enhancing existing minimal comments.

// ============================================================================
// Navigation and "Go To Atom" Functions
// ============================================================================

//! \brief Set the rotation centre to a specific 3D coordinate
//!
//! This function sets the view center to the specified Cartesian coordinates.
//! The camera will rotate around this point.
//!
//! \param pos Cartesian coordinates (Coord_orth object) for the new rotation center
//!
//! \note The view will not automatically move to this position - use go_to_atom
//!       functions to also move the camera
//!
//! Example:
//! \code{.cpp}
//! clipper::Coord_orth new_center(10.0, 20.0, 30.0);
//! set_rotation_centre(new_center);
//! \endcode
void set_rotation_centre(const clipper::Coord_orth &pos);

#ifdef USE_GUILE
//! \brief Navigate to the next atom in the sequence (Guile interface)
//!
//! Given the current atom position, return the specification for the next atom
//! in the sequence. This traverses atoms in order: within a residue, then to
//! the next residue, then to the next chain.
//!
//! \param chain_id Current chain identifier
//! \param resno Current residue number
//! \param ins_code Current insertion code
//! \param atom_name Current atom name
//!
//! \return SCM - List containing the next atom specification [chain_id, resno, ins_code, atom_name]
//!
//! \note Returns the current position if already at the last atom
//!
//! Example usage:
//! \code{.scm}
//! (define next-atom (goto-next-atom-maybe "A" 42 "" "CA"))
//! ;; next-atom will be something like ("A" 42 "" "C") or ("A" 43 "" "N")
//! \endcode
SCM goto_next_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);

//! \brief Navigate to the previous atom in the sequence (Guile interface)
//!
//! Given the current atom position, return the specification for the previous atom
//! in the sequence. This traverses atoms in reverse order: within a residue, then to
//! the previous residue, then to the previous chain.
//!
//! \param chain_id Current chain identifier
//! \param resno Current residue number
//! \param ins_code Current insertion code
//! \param atom_name Current atom name
//!
//! \return SCM - List containing the previous atom specification [chain_id, resno, ins_code, atom_name]
//!
//! \note Returns the current position if already at the first atom
SCM goto_prev_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

#ifdef USE_PYTHON
//! \brief Navigate to the next atom in the sequence (Python interface)
//!
//! Given the current atom position, return the specification for the next atom
//! in the sequence. Useful for iterating through atoms programmatically.
//!
//! \param chain_id Current chain identifier
//! \param resno Current residue number
//! \param ins_code Current insertion code (use "" if none)
//! \param atom_name Current atom name
//!
//! \return PyObject* - List [chain_id, resno, ins_code, atom_name] for next atom
//!
//! Example usage:
//! \code{.py}
//! # Navigate forward through atoms
//! current = ["A", 42, "", "CA"]
//! next_atom = goto_next_atom_maybe_py(*current)
//! print(f"Next atom: {next_atom}")
//! \endcode
PyObject *goto_next_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);

//! \brief Navigate to the previous atom in the sequence (Python interface)
//!
//! Given the current atom position, return the specification for the previous atom.
//! Useful for iterating backwards through atoms programmatically.
//!
//! \param chain_id Current chain identifier
//! \param resno Current residue number
//! \param ins_code Current insertion code (use "" if none)
//! \param atom_name Current atom name
//!
//! \return PyObject* - List [chain_id, resno, ins_code, atom_name] for previous atom
PyObject *goto_prev_atom_maybe_py(const char *chain_id, int resno, const char *ins_code, const char *atom_name);
#endif

//! \brief Center the view on an atom specified by atom_spec
//!
//! Sets the rotation center and moves the view to the specified atom.
//!
//! \param atom_spec Complete atom specification including molecule, chain, residue, and atom
//!
//! \return 1 on success, 0 on failure (atom not found)
//!
//! Example:
//! \code{.cpp}
//! coot::atom_spec_t spec("A", 42, "", "CA", "");
//! int success = set_go_to_atom_from_spec(spec);
//! \endcode
int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec);

//! \brief Center the view on a residue specified by residue_spec
//!
//! Sets the rotation center and moves the view to the CA atom (or first atom)
//! of the specified residue.
//!
//! \param spec Residue specification (chain, residue number, insertion code)
//!
//! \return 1 on success, 0 on failure (residue not found)
//!
//! \note This is useful when you want to center on a residue without
//!       specifying which atom within that residue
//!
//! Example:
//! \code{.cpp}
//! coot::residue_spec_t spec("A", 42, "");
//! set_go_to_atom_from_res_spec(spec);
//! \endcode
int set_go_to_atom_from_res_spec(const coot::residue_spec_t &spec);

#ifdef USE_PYTHON
//! \brief Center the view on a residue (Python interface)
//!
//! Moves the view to center on the specified residue. The view will center
//! on the CA atom if present, otherwise the first atom of the residue.
//!
//! \param residue_spec Python list [chain_id, resno, ins_code]
//!
//! \return 1 on success, 0 if residue not found
//!
//! Example usage:
//! \code{.py}
//! # Center on residue A 42
//! set_go_to_atom_from_res_spec_py(["A", 42, ""])
//!
//! # Center on residue with insertion code
//! set_go_to_atom_from_res_spec_py(["B", 100, "A"])
//! \endcode
int set_go_to_atom_from_res_spec_py(PyObject *residue_spec);

//! \brief Center the view on an atom (Python interface)
//!
//! Moves the view to center on the specified atom.
//!
//! \param atom_spec Python list [chain_id, resno, ins_code, atom_name, alt_conf]
//!
//! \return 1 on success, 0 if atom not found
//!
//! Example usage:
//! \code{.py}
//! # Center on CA atom of residue A 42
//! set_go_to_atom_from_atom_spec_py(["A", 42, "", "CA", ""])
//!
//! # Center on alternate conformation
//! set_go_to_atom_from_atom_spec_py(["A", 50, "", "CA", "A"])
//! \endcode
int set_go_to_atom_from_atom_spec_py(PyObject *atom_spec);
#endif

//! \brief Get the currently active (centered) atom
//!
//! Returns information about the atom currently at the rotation center.
//! Useful for determining what the user is currently looking at.
//!
//! \return Pair where:
//!         - first.first: bool indicating if an atom was found
//!         - first.second: pair of (molecule_number, atom_spec)
//!
//! Example:
//! \code{.cpp}
//! auto result = active_atom_spec();
//! if (result.first) {
//!     int imol = result.second.first;
//!     coot::atom_spec_t spec = result.second.second;
//!     std::cout << "Active atom: " << spec << " in molecule " << imol << std::endl;
//! } else {
//!     std::cout << "No active atom" << std::endl;
//! }
//! \endcode
std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec();

#ifdef USE_PYTHON
//! \brief Get the currently active atom (Python interface)
//!
//! Returns the atom specification for the atom at the current rotation center.
//!
//! \return PyObject* - Tuple (found, (molecule_number, atom_spec)) where:
//!         - found: Boolean indicating if an atom exists at center
//!         - molecule_number: Integer molecule ID
//!         - atom_spec: List [chain_id, resno, ins_code, atom_name, alt_conf]
//!
//! Example usage:
//! \code{.py}
//! found, (imol, atom_spec) = active_atom_spec_py()
//! if found:
//!     chain, resno, ins, atom, alt = atom_spec
//!     print(f"Active: {chain} {resno} {atom} in molecule {imol}")
//! else:
//!     print("No active atom")
//! \endcode
PyObject *active_atom_spec_py();
#endif

// ============================================================================
// Density Scoring Functions
// ============================================================================

#ifdef USE_GUILE
//! \brief Calculate density score for a residue (Guile interface)
//!
//! Computes a numerical score indicating how well a residue fits into
//! the electron density map. Higher scores indicate better fit.
//!
//! \param imol Model molecule number
//! \param residue_spec Scheme list [chain_id, resno, ins_code]
//! \param imol_map Map molecule number
//!
//! \return Float score (typically 0.0 to 1.0+, higher is better)
//!
//! \note This uses a simplified scoring compared to map_to_model_correlation
//! \note Useful for quick assessment but map_to_model_correlation_* functions
//!       provide more detailed analysis
//!
//! Example usage:
//! \code{.scm}
//! (define score (density-score-residue-scm 1 '("A" 42 "") 2))
//! (format #t "Density score: ~a~%" score)
//! \endcode
float density_score_residue_scm(int imol, SCM residue_spec, int imol_map);
#endif

#ifdef USE_PYTHON
//! \brief Calculate density score for a residue (Python interface)
//!
//! Computes how well a residue fits into the electron density map.
//! This is a simpler alternative to the more comprehensive
//! map_to_model_correlation functions.
//!
//! \param imol Model molecule number
//! \param residue_spec Python list [chain_id, resno, ins_code]
//! \param imol_map Map molecule number
//!
//! \return Float score - higher values indicate better fit to density
//!
//! \note For comprehensive density validation, use
//!       map_to_model_correlation_stats_per_residue_range_py() instead
//!
//! Example usage:
//! \code{.py}
//! # Score a single residue
//! score = density_score_residue_py(1, ["A", 42, ""], 2)
//! print(f"Density fit score: {score:.3f}")
//!
//! # Find residues with poor density fit
//! for resno in range(1, 100):
//!     score = density_score_residue_py(1, ["A", resno, ""], 2)
//!     if score < 0.5:
//!         print(f"Poor fit: A {resno}, score = {score:.3f}")
//! \endcode
float density_score_residue_py(int imol, PyObject *residue_spec, int imol_map);
#endif

//! \brief Simple density score for given residue (C++ interface)
//!
//! Calculates a basic density fit score for the specified residue.
//! This function provides a quick assessment of how well atoms fit density.
//!
//! \param imol Model molecule number
//! \param chain_id Chain identifier
//! \param res_no Residue number
//! \param ins_code Insertion code (use "" if none)
//! \param imol_map Map molecule number to score against
//!
//! \return Float score indicating density fit quality (higher is better)
//!
//! \note This is a simplified scoring function. For detailed validation
//!       including correlation statistics, use map_to_model_correlation
//!       functions instead.
//!
//! Example:
//! \code{.cpp}
//! float score = density_score_residue(1, "A", 42, "", 2);
//! if (score < 0.5) {
//!     std::cout << "Residue A 42 has poor density fit: " << score << std::endl;
//! }
//! \endcode
float density_score_residue(int imol, const char *chain_id, int res_no, const char *ins_code, int imol_map);

// ============================================================================
// Map Statistics Functions
// ============================================================================

#ifdef USE_GUILE
//! \brief Get the mean value of a map (Guile interface)
//!
//! Returns the mean (average) density value across all grid points in the map.
//!
//! \param imol Map molecule number
//!
//! \return SCM - Mean value as a number, or #f if imol is not a valid map
//!
//! \note Useful for understanding map scale and detecting data problems
//!
//! Example usage:
//! \code{.scm}
//! (define mean (map-mean-scm 2))
//! (format #t "Map mean: ~a~%" mean)
//! \endcode
SCM map_mean_scm(int imol);

//! \brief Get the standard deviation (sigma) of a map (Guile interface)
//!
//! Returns the standard deviation of density values in the map.
//! This is the "sigma" used for contouring at "N sigma" levels.
//!
//! \param imol Map molecule number
//!
//! \return SCM - Standard deviation as a number, or #f if invalid map
//!
//! \note The contouring level "1.5 sigma" means 1.5 times this value above the mean
//!
//! Example usage:
//! \code{.scm}
//! (define sigma (map-sigma-scm 2))
//! (format #t "Contour at 1.5 sigma = ~a~%" (* 1.5 sigma))
//! \endcode
SCM map_sigma_scm(int imol);

//! \brief Get comprehensive map statistics (Guile interface)
//!
//! Returns detailed statistical measures for the map including mean,
//! standard deviation, skew, and kurtosis.
//!
//! \param imol Map molecule number
//!
//! \return SCM - List (mean std-dev skew kurtosis) or #f if invalid map
//!
//! \note Skew and kurtosis help identify if the map has unusual distributions
//!       that might indicate problems with the data
//!
//! Example usage:
//! \code{.scm}
//! (define stats (map-statistics-scm 2))
//! (if stats
//!     (let ((mean (list-ref stats 0))
//!           (sigma (list-ref stats 1))
//!           (skew (list-ref stats 2))
//!           (kurtosis (list-ref stats 3)))
//!       (format #t "Mean: ~a, Sigma: ~a, Skew: ~a, Kurtosis: ~a~%"
//!               mean sigma skew kurtosis)))
//! \endcode
SCM map_statistics_scm(int imol);
#endif

#ifdef USE_PYTHON
//! \brief Get the mean value of a map (Python interface)
//!
//! Returns the mean (average) density value for all grid points in the map.
//!
//! \param imol Map molecule number
//!
//! \return PyObject* - Float mean value, or False if imol is not a valid map
//!
//! Example usage:
//! \code{.py}
//! mean = map_mean_py(2)
//! if mean is not False:
//!     print(f"Map mean: {mean}")
//! \endcode
PyObject *map_mean_py(int imol);

//! \brief Get the standard deviation (sigma) of a map (Python interface)
//!
//! Returns the standard deviation of density values. This is the "sigma"
//! value used when you set contouring to "1.5 sigma".
//!
//! \param imol Map molecule number
//!
//! \return PyObject* - Float sigma value, or False if invalid map
//!
//! Example usage:
//! \code{.py}
//! sigma = map_sigma_py(2)
//! if sigma is not False:
//!     print(f"1.5 sigma contour level: {1.5 * sigma}")
//! \endcode
PyObject *map_sigma_py(int imol);

//! \brief Get comprehensive map statistics (Python interface)
//!
//! Returns detailed statistical information about the map distribution.
//!
//! \param imol Map molecule number
//!
//! \return PyObject* - List [mean, std_dev, skew, kurtosis] or False if invalid
//!         - mean: Average density value
//!         - std_dev: Standard deviation (sigma)
//!         - skew: Asymmetry of the distribution
//!         - kurtosis: "Tailedness" of the distribution
//!
//! \note Normal distributions have skew≈0 and kurtosis≈3
//! \note Large deviations may indicate data problems
//!
//! Example usage:
//! \code{.py}
//! stats = map_statistics_py(2)
//! if stats is not False:
//!     mean, sigma, skew, kurtosis = stats
//!     print(f"Map statistics:")
//!     print(f"  Mean: {mean:.3f}")
//!     print(f"  Sigma: {sigma:.3f}")
//!     print(f"  Skew: {skew:.3f}")
//!     print(f"  Kurtosis: {kurtosis:.3f}")
//!     
//!     if abs(skew) > 1.0:
//!         print("  Warning: Unusual skew detected")
//! \endcode
PyObject *map_statistics_py(int imol);
#endif

// ============================================================================
// Rotamer Scoring Functions
// ============================================================================

//! \brief Score all possible rotamers for a residue
//!
//! Evaluates all rotameric conformations for the specified residue and
//! returns them ranked by quality. Scoring considers both rotamer probability
//! and fit to density (if a map is provided).
//!
//! \param imol Model molecule number
//! \param chain_id Chain identifier
//! \param res_no Residue number
//! \param ins_code Insertion code (use "" if none)
//! \param alt_conf Alternate conformation identifier (use "" for default)
//! \param imol_map Map molecule number for density scoring (use -1 to skip)
//! \param clash_flag Check for clashes: 1=check, 0=don't check
//! \param lowest_probability Minimum rotamer probability to consider (e.g., 0.01 for 1%)
//!
//! \return Vector of named_rotamer_score objects, each containing:
//!         - Rotamer name (e.g., "mt-85", "tp175")
//!         - Probability score
//!         - Density fit score (if map provided)
//!         - Richardson rotamer name
//!
//! \note Results are typically sorted by overall quality (probability × density_fit)
//! \note Use lowest_probability to filter out very unlikely rotamers
//!
//! Example:
//! \code{.cpp}
//! // Score rotamers for residue A 42, considering density and clashes
//! auto rotamers = score_rotamers(1, "A", 42, "", "", 2, 1, 0.01);
//! 
//! std::cout << "Found " << rotamers.size() << " rotamers above 1% probability" << std::endl;
//! for (const auto &rot : rotamers) {
//!     std::cout << rot.name << ": prob=" << rot.probability 
//!               << ", fit=" << rot.density_score << std::endl;
//! }
//! \endcode
std::vector<coot::named_rotamer_score> score_rotamers(int imol,
                                                      const char *chain_id,
                                                      int res_no,
                                                      const char *ins_code,
                                                      const char *alt_conf,
                                                      int imol_map,
                                                      int clash_flag,
                                                      float lowest_probability);

#ifdef USE_GUILE
//! \brief Score rotamers for a residue (Guile interface)
//!
//! Returns a list of possible rotamer conformations with their scores.
//! Each rotamer is scored based on rotamer library probability and
//! (optionally) density fit.
//!
//! \param imol Model molecule number
//! \param chain_id Chain identifier
//! \param res_no Residue number
//! \param ins_code Insertion code
//! \param alt_conf Alternate conformation
//! \param imol_map Map for density scoring (-1 to skip)
//! \param clash_flag 1 to check clashes, 0 to skip
//! \param lowest_probability Minimum probability threshold (0.0-1.0)
//!
//! \return SCM - List of rotamer descriptions, each containing rotamer name,
//!         probability, and density score. Empty list if residue not found
//!         or no rotamers above threshold.
//!
//! \note The density score is only meaningful if imol_map is a valid map
//!
//! Example usage:
//! \code{.scm}
//! ;; Score rotamers for LEU 42 in chain A
//! (define rotamers (score-rotamers-scm 1 "A" 42 "" "" 2 1 0.01))
//! (for-each
//!   (lambda (rot)
//!     (format #t "Rotamer: ~a, Probability: ~a, Fit: ~a~%"
//!             (list-ref rot 0)  ; name
//!             (list-ref rot 1)  ; probability
//!             (list-ref rot 2))) ; density fit
//!   rotamers)
//! \endcode
SCM score_rotamers_scm(int imol,
                       const char *chain_id,
                       int res_no,
                       const char *ins_code,
                       const char *alt_conf,
                       int imol_map,
                       int clash_flag,
                       float lowest_probability);
#endif

#ifdef USE_PYTHON
//! \brief Score all rotamers for a residue (Python interface)
//!
//! **USEFUL FOR FIXING BAD ROTAMERS**
//!
//! Evaluates all possible rotamer conformations and returns them with scores.
//! This is the function to call before using auto_fit_best_rotamer.
//!
//! \param imol Model molecule number
//! \param chain_id Chain identifier
//! \param res_no Residue number
//! \param ins_code Insertion code (use "" if none)
//! \param alt_conf Alternate conformation (use "" for default)
//! \param imol_map Map molecule for density scoring (use -1 to ignore density)
//! \param clash_flag 1 to check for clashes with other atoms, 0 to skip
//! \param lowest_probability Filter: only return rotamers above this probability
//!
//! \return PyObject* - List of rotamer dictionaries, each containing:
//!         - 'name': Rotamer name (e.g., "mt-85")
//!         - 'probability': Rotamer library probability (0.0-1.0)
//!         - 'density_score': Fit to density (if map provided)
//!         - 'richardson_name': Rotamer name in Richardson notation
//!         Empty list if no suitable rotamers found.
//!
//! \note Rotamers are ranked by combined probability and density fit
//! \note Use clash_flag=1 to avoid rotamers that clash with nearby atoms
//!
//! Example usage:
//! \code{.py}
//! # Score rotamers for LEU 42, considering density and clashes
//! rotamers = score_rotamers_py(
//!     imol=1,
//!     chain_id="A",
//!     res_no=42,
//!     ins_code="",
//!     alt_conf="",
//!     imol_map=2,           # Use map 2 for density scoring
//!     clash_flag=1,         # Check for clashes
//!     lowest_probability=0.01  # Only show rotamers >1% probability
//! )
//!
//! print(f"Found {len(rotamers)} possible rotamers")
//! for rot in rotamers:
//!     print(f"{rot['name']}: "
//!           f"prob={rot['probability']:.1%}, "
//!           f"fit={rot['density_score']:.3f}")
//!
//! # The best rotamer is typically first in the list
//! if rotamers:
//!     best = rotamers[0]
//!     print(f"Best rotamer: {best['name']}")
//! \endcode
PyObject *score_rotamers_py(int imol,
                            const char *chain_id,
                            int res_no,
                            const char *ins_code,
                            const char *alt_conf,
                            int imol_map,
                            int clash_flag,
                            float lowest_probability);
#endif

// ============================================================================
// Utility Functions
// ============================================================================

//! \brief Get the Git commit hash of the Coot build
//!
//! Returns the Git commit identifier for the version of Coot currently running.
//! Useful for bug reports and version tracking.
//!
//! \return String containing the Git commit hash (e.g., "a3f2c1d")
//!
//! Example:
//! \code{.cpp}
//! std::string commit = git_commit();
//! std::cout << "Coot version: commit " << commit << std::endl;
//! \endcode
std::string git_commit();

//! \brief Filter files by glob pattern
//!
//! Returns a list of files in the specified directory that match the given
//! glob pattern, useful for file selection dialogs.
//!
//! \param pre_directory Directory path to search
//! \param data_type Type of data to filter (specific values TBD)
//!
//! \return Vector of filename strings matching the criteria
std::vector<std::string> filtered_by_glob(const std::string &pre_directory, int data_type);

//! \brief Check if a string exists in a vector
//!
//! Searches for an exact match of the search string in the provided list.
//!
//! \param search String to search for
//! \param list Vector of strings to search within
//!
//! \return 1 if found, 0 if not found
//!
//! Example:
//! \code{.cpp}
//! std::vector<std::string> chains = {"A", "B", "C"};
//! if (string_member("B", chains)) {
//!     std::cout << "Chain B exists" << std::endl;
//! }
//! \endcode
short int string_member(const std::string &search, const std::vector<std::string> &list);

//! \brief Compare two strings
//!
//! Performs string comparison for sorting purposes.
//!
//! \param a First string
//! \param b Second string
//!
//! \return true if a < b in lexicographic order
//!
//! \note Useful as a comparator function for std::sort
bool compare_strings(const std::string &a, const std::string &b);

// ============================================================================
// CaBLAM Validation
// ============================================================================

//! \brief Add CaBLAM validation markup to a model
//!
//! CaBLAM (CA B-factor, Length, Angle, Measurement) is a protein backbone
//! validation tool. This function reads CaBLAM output and adds visual markup
//! to highlight problematic regions.
//!
//! \param imol Model molecule number
//! \param cablam_file_name Path to CaBLAM output file
//!
//! \return Vector of pairs (residue_spec, score) for all validated residues
//!
//! \note Lower CaBLAM scores indicate poorer backbone geometry
//! \note Run CaBLAM externally and provide the output file to this function
//!
//! Example:
//! \code{.cpp}
//! auto results = add_cablam_markup(1, "cablam_output.txt");
//! for (const auto &[spec, score] : results) {
//!     if (score < 0.5) {
//!         std::cout << "Poor CaBLAM score: " << spec << " = " << score << std::endl;
//!     }
//! }
//! \endcode
std::vector<std::pair<coot::residue_spec_t, double> >
add_cablam_markup(int imol, const std::string &cablam_file_name);

#ifdef USE_PYTHON
//! \brief Add CaBLAM validation markup (Python interface)
//!
//! Reads CaBLAM output and adds colored markup to the model showing
//! backbone validation results.
//!
//! \param imol Model molecule number  
//! \param cablam_log_file_name Path to CaBLAM output file
//!
//! \return PyObject* - List of tuples [(residue_spec, score), ...]
//!         where residue_spec is [chain_id, resno, ins_code]
//!         and score is the CaBLAM validation score
//!
//! Example usage:
//! \code{.py}
//! # Read CaBLAM results and add markup
//! results = add_cablam_markup_py(1, "cablam_results.txt")
//!
//! # Find worst CaBLAM outliers
//! outliers = [r for r in results if r[1] < 0.5]
//! print(f"Found {len(outliers)} CaBLAM outliers")
//!
//! for (chain, resno, ins), score in outliers:
//!     print(f"  {chain} {resno}: score = {score:.3f}")
//! \endcode
PyObject *add_cablam_markup_py(int imol, const std::string &cablam_log_file_name);
#endif

