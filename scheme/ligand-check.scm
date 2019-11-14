;;;; Copyright 2015 by Medical Research Council

;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


(use-modules (json reader)
	     (ice-9 format)
	     (ice-9 q))


;; we need a function that takes 2 molecules, 
;; a protein molecule and a ligand molecule,
;;    merges them to create a new molecule 
;;    then run probe on the complex - but only looking at
;;    protein->ligand score (actually that might mean that we give 2
;;    files to probe - look at isolate-ligand.scm) but we need the
;;    hydrogenated molecules of course - and the hydrogens should be
;;    added in the complex - so we *will* have to do the merge (and a
;;    possible split afterwards) need to generate a score from the
;;    probe result (how does clashscore do it? - should be easy)

;; (load "cache-or-net-get.scm")


(load "contact-score-isolated-ligand.scm")
(load "run-mogul.scm")


;; These are fei's column labels
;;
;; (define *refmac-data-f-col*     "FP")
;; (define *refmac-data-sigf-col*  "SIGFP")
;; (define *refmac-data-rfree-col* "FREE")

(define mogul-results mogul-results-scm)

;; test:
;; (using-active-atom (get-metrics-for-ligand aa-imol aa-chain-id aa-res-no aa-ins-code "rnasa-1.8-all_refmac1.mtz" "FGMP18" "SIGFGMP18" "FreeR_flag" "coot-refmac"))
;; 
;; refmac-dir: the dir where input/output refmac files should be written.
;; 
;; return (list correlation mogul bumps difference-map-stats b-factor-info)
;; where:
;;     mogul is an improper pair (mogul-z-worst . mogul-out-file-name) or an error status
;;     difference-map-stats:
;;       the output of map-to-model-correlation-stats-scm
;;       (list correlation var_x var_y n sum_x sum_y D D2 (based on mean/sd of the map at the ligand) 
;;             map_mean map_mean_at_ligand map_sd map_sd_at_ligand SCM_UNDEFINED)
;; and b-factor-info
;;       list of (median-ratio median-ligand median-env ks-test-result)
;; 
(define (get-metrics-for-ligand imol chain-id res-no ins-code refmac-input-mtz-file-name fobs-col sig-fobs-col rfree-col refmac-dir)

  ;; get-correlation and get-ligand-difference-map-stats are very similar
  ;; Let's factor out the similarities.
  ;; 
  ;; success-function takes args: refmac-out-sfs-file-name
  ;;                              f-col phi-col
  ;;                              is-difference-map-flag
  ;;                              ligand-spec
  ;;                              neighbs
  ;; is that true?
  ;; 
  (define (local-refmac stub-name)

    (let* ((ligand-spec (list chain-id res-no ins-code))
	   (neighbs (residues-near-residue imol ligand-spec 4)))

      (let* ((rn (residue-name imol chain-id res-no ins-code)))

	(if (not (string? rn))
	    (begin
	      (format #t "No residue name imol: ~s chain-id ~s res-no ~s~%" imol chain-id res-no)
	      'fail-no-residue-name)

	    (let ((n-ligand-atoms (het-group-n-atoms rn)))

	      (if (not (number? n-ligand-atoms))
		  'fail-n-ligand-atoms-not-a-number

		  ;; 
		  (begin
		    ;; (delete-residue imol chain-id res-no ins-code)
		    (let ((refmac-out-sfs-file-name (append-dir-file
						     refmac-dir
						     (string-append
						      stub-name
						      "-with-ligand-refmac.mtz")))
			  (with-ligand-pdb-file-name (append-dir-file
						      refmac-dir
						      (string-append
						       stub-name "-with-ligand.pdb"))))

		      (make-directory-maybe refmac-dir)
		      (make-directory-maybe "coot-refmac") ;; XYZOUT goes here
		      (write-pdb-file imol with-ligand-pdb-file-name)
		      (let ((r (refmac-calc-sfs-make-mtz-with-columns with-ligand-pdb-file-name
								      refmac-input-mtz-file-name
								      refmac-out-sfs-file-name
								      fobs-col sig-fobs-col rfree-col)))
			(if (eq? r #f)

			    (begin 
			      ;; test if refmac-input-mtz-file-name has intensites
			      'fail-problem-calculating-sfs-using-refmac)

			    ;; happy path
			    refmac-out-sfs-file-name))))))))))

  ;; get-correlation at the ligand for the direct (FWT) map 
  ;;
  (define (get-correlation stub-name)

    (set! refmac-extra-params (list (string-append
				     "refine exclude all from " 
				     (number->string res-no)
				     " " 
				     chain-id
				     " to " 
				     (number->string res-no)
				     " " 
				     chain-id
				     )))

    (format #t "in get-correlation(): refmac-extra-params: ~s~%" refmac-extra-params)

    (let ((ligand-spec (list chain-id res-no ins-code))
          (refmac-out-sfs-file-name (local-refmac stub-name)))

      (if (not (string? refmac-out-sfs-file-name))

	  refmac-out-sfs-file-name

	  ;; happy path
	  (let ((imol-map (make-and-draw-map refmac-out-sfs-file-name
					     "FWT" "PHWT" "" 0 0))
	        (neighbs (residues-near-residue imol ligand-spec 4)))

	    (let ((c (map-to-model-correlation imol (list ligand-spec)
					       neighbs 0 imol-map)))
	      (close-molecule imol-map)
	      c)))))

  ;; return a failure symbol or a list of stats.
  ;; 
  (define (get-ligand-difference-map-stats stub-name)
    (set! refmac-extra-params #f) ;; reset - no ligand exclusion
    
    (let* ((refmac-out-sfs-file-name 
            (local-refmac (string-append stub-name "-for-ligand-diff-map")))
           (ligand-spec (list chain-id res-no ins-code))
           (neighbs (residues-near-residue imol ligand-spec 4)))

      (if (not (string? refmac-out-sfs-file-name))

	  refmac-out-sfs-file-name

	  ;; happy path
	  (let ((imol-map (make-and-draw-map refmac-out-sfs-file-name
					     "DELFWT" "PHDELWT" "" 0 1)))
	    ;; now do some stats on the map at the ligand site

            (let ((c (map-to-model-correlation-stats-scm imol (list ligand-spec) neighbs 10 imol-map)))
                (format #t "### residue ~s density statistics ~s~%~!" ligand-spec c)
                c)
	    ))))

  ;; Return an error status (not a list) or a improper pair
  ;; (mogul-z-worst . mogul-out-file-name) on success.
  ;; 
  ;; (we want the filename so that we can show mogul markup on the ligand when we activate this 
  ;;  function from the gui)
  (define (get-mogul-score use-cache?)

    ;; if run-result is a string, then it is the mogul-output file name
    ;; if it is a symbol something went wrong.
    ;; 
    (let ((run-result (run-mogul 'bonds-and-angles imol chain-id res-no ins-code "ligand-check" use-cache?)))
      (format #t "run-result (mogul): ~s ~s ~s~%" chain-id res-no run-result)

      (if (not (string? run-result))

	  run-result

	  ;; happy path
	  (let* ((mogul-results-list (mogul-results run-result))
		 (mogul-score (if (not (list? mogul-results-list))
				  mogul-results-list
				  (if (= (length mogul-results-list) 0)
				      'mogul-no-stats
				      (apply max mogul-results-list)))))
	    (cons mogul-score run-result)))))

  ;; return a list: n_bad_overlaps n_hydrogen_bonds n_small_overlaps n_close_contacts n_wide_contacts
  ;; 
  (define (get-bump-score)
    (let* ((ligand-spec (list chain-id res-no ins-code))
	   (cs (contact-score-ligand imol ligand-spec)))
      ;; (format #t "debug:: contact score cs: ~s~%" cs)
      (graphics-draw)
      cs))

  ;; return a list of (median-ratio median-ligand median-env ks-test-result)
  ;; stub-name is only passed so that we can write diagnostics to standard-out.
  ;; 
  (define (get-b-factor-distribution-metrics stub-name)

    (define (median-ratio ligand-b-factors env-b-factors)
      (/ (median ligand-b-factors) (median (apply append env-b-factors))))

    (define (filter-out-waters imol env-residues)
      (filter (lambda (residue-item)
		(let ((rn (residue-name imol
					(residue-spec->chain-id residue-item)
					(residue-spec->res-no   residue-item)
					(residue-spec->ins-code residue-item))))
		  (not (or (string-match "HOH" rn)
			   (string-match "WAT" rn)))))
	      env-residues))

    (define (median number-list)
      (print-var number-list)
      (if (not (list? number-list))
         1.0
         (let ((numbers
	        (let loop ((number-list number-list)
			(nums '()))
	          (cond
	             ((null? number-list) nums)
                     ((not (number? (car number-list)))
                        (loop (cdr number-list) nums))
                  (else
                     (loop (cdr number-list) (cons (car number-list) nums)))))))
           (print-var numbers)
	   (let* ((sorted-nums (sort-list numbers >))
	          (n (length sorted-nums))
	          (mid (/ n 2))
                  (nov (print-var sorted-nums))
	          (median (if (= (remainder n 2) 0)
			   (begin
			     ;; (format #t "averaging mid: ~s  n: ~s ~%" mid n)
			     ;; (format #t "idx: ~s  idx: ~s ~%" (/ n 2) (- (/ n 2) 1))
			     (/ (+ (list-ref sorted-nums (/ n 2))
				   (list-ref sorted-nums (- (/ n 2) 1))) 2))
			   (begin
			     ;; (format #t "simple take: ~s~%" mid)
			     (list-ref sorted-nums (/ (- n 1) 2))))))
	     median))))


    ;; Return a list of length 2: 
    ;; The car is a list of atoms for the residue specified by ligand-spec
    ;; The car may be #f, in which case this function failed to return a result
    ;; If the car is not #f, the cdr is a list list of atom b-factors
    ;; 
    (define (ligand-environment-temperature-factors imol ligand-spec radius)
      (let ((atoms (residue-info imol
				 (residue-spec->chain-id ligand-spec)
				 (residue-spec->res-no   ligand-spec)
				 (residue-spec->ins-code ligand-spec))))
	(let* ((env-residues (residues-near-residue imol ligand-spec radius))
	       (non-water-env-residues (filter-out-waters imol env-residues))
	       (env-atoms (map (lambda (res-spec)
				 (residue-info imol 
					       (residue-spec->chain-id res-spec)
					       (residue-spec->res-no   res-spec)
					       (residue-spec->ins-code res-spec))) non-water-env-residues)))
	  (list (if (list? atoms)
		    (map (lambda (atom)
			   (let* ((occ (cadr (cadr atom))))
			     (if (list? occ) (car occ) occ)))
			 atoms)
		    #f)
		(map (lambda (atom-list)
		       (map (lambda (atom)
			      (let* ((occ (cadr (cadr atom))))
				(if (list? occ) (car occ) occ)))
			    atom-list))
		     env-atoms)))))

    ;; main line of get-b-factor-distribution-metrics
    ;;
    (let ((ligand-spec (list chain-id res-no ins-code)))
      (print-var ligand-spec)

      (let* ((lig-env-temp-factors (ligand-environment-temperature-factors imol ligand-spec 5)))
         (format #t "$$$$$$$$$$$$\n")
         (print-var lig-env-temp-factors)
         (if (not (list? lig-env-temp-factors))
             (begin
                (format #t  "Ligand env temp factors not a list\n")
                #f)

             (begin
                (if (= (length (list-ref lig-env-temp-factors 1)) 0)
                   (begin
                      (format #t  "No values in Ligand env temp factors\n")
                      #f)

                   (let ((mr (apply median-ratio lig-env-temp-factors)))
                      (let ((temp-factor-median-ratio mr))

                         (let ((v1 (car lig-env-temp-factors))
	                       (v2 (cadr lig-env-temp-factors)))
	                    (format #t "b-factor kolmogorov-smirnov lig: ~s ~s ~s~%~!" stub-name ligand-spec (car lig-env-temp-factors))
	                    (format #t "b-factor kolmogorov-smirnov env: ~s ~s ~s~%~!" stub-name ligand-spec (cadr lig-env-temp-factors))
	                    (let ((kolmogorov-smirnov-result (kolmogorov-smirnov v1 (apply append v2))))
	                      (list temp-factor-median-ratio 
		                    (median (car lig-env-temp-factors))
		                    (median (apply append (cadr lig-env-temp-factors)))
		                    kolmogorov-smirnov-result)))))))))))


;   pre-20150803-PE  
;   ;; main line of get-metrics-for-ligand
;   ;; 
;   (let* ((stub-name (molecule-name-stub imol 0)))
;     (let ((cor (get-correlation stub-name)))
;       (if (number? cor)
; 	  (let ((mog (get-mogul-score #f))) ;; use the cache for the ligand? - testing only!
; 	    (if (number? mog)
; 		(let ((bmp (get-bump-score)))
; 		  (if (number? bmp)
; 		      (list cor mog bmp)
; 		      bmp)) ;; error symbol/string
; 		mog)) ;; error symbol/string
; 	  cor)))) ;; error symbol/string


  ;; main line of get-metrics-for-ligand
  ;; 
  (let* ((stub-name (molecule-name-stub imol 0)))

    (let ((b-factor-info (get-b-factor-distribution-metrics stub-name)))

      ;; add error checking to this 
      ;; 
      (let ((cor (get-correlation stub-name)))
	(if (number? cor)
	    (let ((dms (get-ligand-difference-map-stats stub-name)))
	      (if (not (list? dms))
		  dms ;; error symbol
		  (let ((mog (get-mogul-score #f))) ;; use the cache for the ligand? - testing only!

		    (if (pair? mog)
			(let ((bmp (get-bump-score)))

			  ;; (format #t "------------- bmp: ~s~%" bmp)
			  (if (list? bmp)
			      (list cor mog bmp dms b-factor-info)
			      bmp)) ;; error symbol/string
			mog)))) ;; error symbol/string
	    cor))))) ;; error symbol/string


;; only look at ligands in maps with resolution worse than this:
;; (to look at everything, set it to 0.1)
;; 
(define *ligand-check-resolution-limit* 1.2)


;; remove residues that are waters from env-residues
;; 
(define (filter-out-waters imol env-residues)
  (filter (lambda (residue-item)
	    (let ((rn (residue-name imol
				    (residue-spec->chain-id residue-item)
				    (residue-spec->res-no   residue-item)
				    (residue-spec->ins-code residue-item))))
	      (not (or (string-match "HOH" rn)
		       (string-match "WAT" rn)))))
	  env-residues))
