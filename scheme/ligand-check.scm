
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
(define (get-metrics-for-ligand imol chain-id res-no ins-code 
				refmac-input-mtz-file-name fobs-col sig-fobs-col rfree-col
				refmac-dir)

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
  (define (local-refmac stub-name refmac-out-sfs-file-name success-function)
    #f)


  ;; get-correlation at the ligand for the direct (FWT) map 
  ;;
  (define (get-correlation stub-name)
    (let* ((ligand-spec (list chain-id res-no ins-code))
	   (neighbs (residues-near-residue imol ligand-spec 4)))

      (let* ((rn (residue-name imol chain-id res-no ins-code))
	     (n-ligand-atoms (het-group-n-atoms rn)))

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

		      'fail-problem-calculating-sfs-using-refmac

		      ;; happy path
		      (let ((imol-map (make-and-draw-map refmac-out-sfs-file-name
							 "FWT" "PHWT" "" 0 0)))
			(let ((c (map-to-model-correlation imol (list ligand-spec)
							   neighbs 0 imol-map)))
			  (close-molecule imol-map)
			  c))))))))))

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
	    mogul-score))))

  (define (get-bump-score)
    (let* ((ligand-spec (list chain-id res-no ins-code))
	   (cs (contact-score-ligand imol ligand-spec)))
      (format #t "debug:: contact score cs: ~s~%" cs)
      (graphics-draw)
      (if (list? cs)
	  (car cs)
	  cs)))
  
  ;; main line of get-metrics-for-ligand
  ;; 
  (let* ((stub-name (molecule-name-stub imol 0)))
    (let ((cor (get-correlation stub-name)))
      (if (number? cor)
	  (let ((mog (get-mogul-score #f))) ;; use the cache for the ligand? - testing only!
	    (if (number? mog)
		(let ((bmp (get-bump-score)))
		  (if (number? bmp)
		      (list cor mog bmp)
		      bmp)) ;; error symbol/string
		mog)) ;; error symbol/string
	  cor)))) ;; error symbol/string


;; only look at ligands in maps with resolution worse than this:
;; (to look at everything, set it to 0.1)
;; 
(define *ligand-check-resolution-limit* 1.2)


;; remove residues that are waters from env-residues
(define (filter-out-waters imol env-residues)
  (filter (lambda (residue-item)
	    (let ((rn (residue-name imol
				    (residue-spec->chain-id residue-item)
				    (residue-spec->res-no   residue-item)
				    (residue-spec->ins-code residue-item))))
	      (not (or (string-match "HOH" rn)
		       (string-match "WAT" rn)))))
	  env-residues))


