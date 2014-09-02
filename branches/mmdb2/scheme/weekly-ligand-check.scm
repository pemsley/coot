
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

(load "cache-or-net-get.scm")
(load "contact-score-isolated-ligand.scm")
(load "run-mogul.scm")

(define *min-n-atoms-for-ligand* 7) ;; 7 excludes GOL


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

  (define (get-correlation stub-name)
    (let* ((ligand-spec (list chain-id res-no ins-code))
	   (neighbs (residues-near-residue imol ligand-spec 4)))

      (let* ((rn (residue-name imol chain-id res-no ins-code))
	     (n-ligand-atoms (het-group-n-atoms rn)))

	(if (not (number? n-ligand-atoms))
	    'n-ligand-atoms-not-a-number

	    ;; 
	    (let ((imol-new (copy-molecule imol)))
	      (delete-residue imol-new chain-id res-no ins-code)
	      (let ((refmac-out-sfs-file-name (append-dir-file
					       refmac-dir
					       (string-append
						stub-name
						"-sans-ligand-refmac.mtz")))
		    (sans-ligand-pdb-file-name (append-dir-file
						refmac-dir
						(string-append
						 stub-name "-sans-ligand.pdb"))))

		(make-directory-maybe refmac-dir)
		(write-pdb-file imol-new sans-ligand-pdb-file-name)
		(set! *refmac-data-f-col* fobs-col)
		(set! *refmac-data-sigf-col* sig-fobs-col)
		(set! *refmac-data-rfree-col* rfree-col)
		(let ((r (refmac-calc-sfs-make-mtz sans-ligand-pdb-file-name
						   refmac-input-mtz-file-name
						   refmac-out-sfs-file-name)))
		  (if (eq? r #f)

		      'problem-calculating-sfs-using-refmac

		      ;; happy path
		      (let ((imol-map (make-and-draw-map refmac-out-sfs-file-name
							 "FWT" "PHWT" "" 0 0)))
			(map-to-model-correlation
			 imol 
			 (list ligand-spec)
			 neighbs
			 0 imol-map))))))))))
  

  (define (get-mogul-score use-cache?)
    ;; if run-result is a string, then it is the mogul-output file name
    ;; if it is a symbol something went wrong.
    (let ((run-result (run-mogul 'bonds-and-angles imol chain-id res-no ins-code "vvv" use-cache?)))
      (if (not (string? run-result))
	  run-result
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
	  (let ((mog (get-mogul-score #t))) ;; use the cache for the ligand - testing only!
	    (if (number? mog)
		(let ((bmp (get-bump-score)))
		  (if (number? bmp)
		      (list cor mog bmp)
		      #f))
		#f))
	  #f))))


(define (score-ligand accession-code imol imol-map ligand-spec res-name)

  ;; get the correlation of the ligand
  (let* ((neighbs (residues-near-residue imol ligand-spec 4))
	 (c (map-to-model-correlation
	     imol 
	     (list ligand-spec)
	     neighbs
	     0 imol-map)))

    (format #t "---------------- correlation c : ~s for imol ~s specs-list ~s with imol-map ~s ~%" 
	    c imol (list ligand-spec) imol-map)
    
    (let* ((mogul-run (run-mogul-show-results imol accession-code 
					      'bonds-and-angles
					      (residue-spec->chain-id ligand-spec)
					      (residue-spec->res-no   ligand-spec)
					      (residue-spec->ins-code ligand-spec)))
	   (nov (format #t "----- mogul-run: ~s~%" mogul-run))
	   (mogul-results-list (if (not (string? mogul-run))
				   mogul-run
				   (mogul-results mogul-run)))
	   ;; average z score
	   ;; 
	   (mogul-score-1 (if (not (list? mogul-results-list))
			      mogul-results-list
			      (let ((n (length mogul-results-list)))
				(if (= n 0)
				    'mogul-no-stats
				    (/ (apply + mogul-results-list) n)))))
	   ;; worse z score
	   ;; 
	   (mogul-score (if (not (list? mogul-results-list))
			    mogul-results-list
			    (if (= (length mogul-results-list) 0)
				'mogul-no-stats
				(apply max mogul-results-list))))

	   (lig-env-temp-factors (ligand-environment-temperature-factors imol ligand-spec 5))

	   (temp-factor-median-ratio (apply median-ratio lig-env-temp-factors))
					
	   (kolmogorov-smirnov-result (let ((v1 (car lig-env-temp-factors))
					    (v2 (cadr lig-env-temp-factors)))
					(kolmogorov-smirnov v1 (apply append v2)))))
      
      ;; get contact score (probe and reduce)
      (let ((cs (contact-score-ligand imol ligand-spec)))
	(format #t "::::::::::: ---- code ~s ---~%" accession-code)
	(format #t "::::::::::: density correl: ~5f for ~s ~s ~%" c ligand-spec res-name)
	(format #t "::::::::::: clash score: ~s~%" cs)
	(format #t "::::::::::: mogul score : ~s~%" mogul-score)
	(format #t "::::::::::: mogul results : ~s~%" mogul-results-list)
	
	(list accession-code c cs mogul-results-list temp-factor-median-ratio ligand-spec res-name 
        ;;         0         1  2      3                         4                5          6

	      kolmogorov-smirnov-result 
        ;;              7
         )))))

;; Is there a LINK record for the given ligand spec?  (any model).
;; 
(define (linked-ligand? ligand-spec imol)

  (define (atom-spec->residue-spec atom-spec)
    (list-head atom-spec 3))

  (define (specs-match? ligand-spec link-res-spec-1)
    (if (not (string=? (residue-spec->chain-id ligand-spec)
		       (residue-spec->chain-id link-res-spec-1)))
	#f
	(if (not (= (residue-spec->res-no ligand-spec)
		    (residue-spec->res-no link-res-spec-1)))
	    #f
	    (if (not (string=? (residue-spec->ins-code ligand-spec)
			       (residue-spec->ins-code link-res-spec-1)))
		#f
		#t))))

  ;; main line
  (let ((links (link-info imol)))
    (let loop ((link-infos links))
      (cond
       ((null? link-infos) #f)
       (else 
	(let* ((link-info (car link-infos))
	       (link-atom-spec-1 (cdr (list-ref link-info 1)))
	       (link-atom-spec-2 (cdr (list-ref link-info 2)))
	       (link-res-spec-1 (atom-spec->residue-spec link-atom-spec-1))
	       (link-res-spec-2 (atom-spec->residue-spec link-atom-spec-2)))

	  (if (specs-match? ligand-spec link-res-spec-1)
	      #t
	      (if (specs-match? ligand-spec link-res-spec-2)
	      #t
	      (loop (cdr link-infos))))))))))



;; we already have the files on disk
(define (get-ligand-metrics-for-code accession-code)

  (format #t "--------------------------- accession-code: ~s --------------------------~%" 
	  accession-code)

  ;; we should test if these files are cached first
  ;;
  (let* ((pdb-file-name (append-dir-file "coot-download" (string-append accession-code ".pdb")))
	 (cif-file-name (append-dir-file "coot-download" (string-append "r" accession-code 
								      "sf.cif"))))

    (if (not (file-exists? pdb-file-name))
	(begin
	  (format #t "pdb file does not exist ~s~%" pdb-file-name))

	(begin
	  ;; happy path
	  (if (not (file-exists? cif-file-name))
	      (begin
		(format #t "cif file does not exist ~s~%" cif-file-name))
	  
	      
	      ;; OK so we have the files 
	      (let* ((sfs-mtz-file-name (append-dir-file "coot-download" 
							(string-append "r" 
								       accession-code
								       "sfs.mtz")))
		     (convert-status (if (file-exists? sfs-mtz-file-name) 
					 1
					 (mmcif-sfs-to-mtz cif-file-name sfs-mtz-file-name))))

		(if (not (= convert-status 1))
		    ;; 
		    (begin
		      (format #t "problem in converting~s~%" cif-file-name))

		    ;; OK, happy path
		    (begin
		      (let* ((imol (read-pdb pdb-file-name))
			     ;; returns a list of the new molecule number and the spec of the 
			     ;; removed residue (or scheme false).
			     (imol-new-info (new-molecule-sans-biggest-ligand-scm imol)))
			(if (not (list? imol-new-info))
			    (begin
			      (format #t "No (interesting) ligand ~s ~s~%" 
				      pdb-file-name imol-new-info)
			      'no-interesting-ligand-here)

			    ;; happy path
			    (begin

			      (let ((refmac-out-sfs-file-name (append-dir-file
							       "coot-download"
							       (string-append  
								"r" 
								accession-code
								"-sans-ligand-refmac.mtz"))))
				(let* ((ligand-spec (cdr (cadr imol-new-info)))
				       ;; (nov (format #t ":::::::::::::::: ligand-spec: ~s~%" 
				       ;; ligand-spec))
				       (imol-new (car imol-new-info))
				       (rn (residue-name imol 
							 (residue-spec->chain-id ligand-spec)
							 (residue-spec->res-no   ligand-spec)
							 (residue-spec->ins-code ligand-spec)))
				       (n-ligand-atoms (het-group-n-atoms rn)))

				  ;; hack in env check.  should be somewhere else
				  ;; (let ((b-factors-set 
				  ;;        (ligand-environment-temperature-factors imol ligand-spec 4)))
				  ;;    (format #t ":::::::::::::: b-factors-set: ~s~%" b-factors-set))

				  (if (not (number? n-ligand-atoms))
				      'n-ligand-atoms-not-a-number

				      (if (not (>= n-ligand-atoms *min-n-atoms-for-ligand*))
					  (begin 
					    (format #t "INFO:: insufficient atoms in dict for type ~s, found ~s, need ~s~%"
						    rn n-ligand-atoms *min-n-atoms-for-ligand*)
					    'n-ligand-atoms-insufficient)

					  (if (linked-ligand? ligand-spec imol)
					      
					      'linked-ligand
					      
					      (if (string-member? rn (list "MSE" "OCS" "SEP" "LE1" "HEM" "HEC"))

						  'non-interesting-het-group

						  ;; now we have a modified molecule, let's write it out so
						  ;; that we can use it with refmac.
						  ;; 
						  (let ((sans-ligand-pdb-file-name
							 (append-dir-file
							  "coot-download"
							  (string-append
							   accession-code
							   "-sans-ligand.pdb"))))
						    (write-pdb-file imol-new sans-ligand-pdb-file-name)

						    ;; let's calc sfs (or check the cache)
						    ;; 
						    (let ((r (if (file-exists? refmac-out-sfs-file-name) 
								 #t
								 (refmac-calc-sfs-make-mtz sans-ligand-pdb-file-name
											   sfs-mtz-file-name
											   refmac-out-sfs-file-name))))
						      (if (eq? #f r)
							  (begin
							    (format #t "problem calculating sfs with refmac ~s ~s~%"
								    sans-ligand-pdb-file-name
								    sfs-mtz-file-name)
							    'no-refmac-output-mtz)
							  
							  ;; happy path
							  (begin
							    (let ((imol-map (make-and-draw-map
									     refmac-out-sfs-file-name
									     "FWT" "PHWT" "" 0 0)))
							      (if (not (valid-map-molecule? imol-map))
								  (begin
								    (format #t "problem making map from ~s~%" 
									    refmac-out-sfs-file-name)
								    'problem-making-map-from-refmac-mtz)

								  (score-ligand accession-code imol imol-map ligand-spec rn))))))))))))))))))))))))

;; 
(define (get-pdbe-files-for-code code)
  (let* ((pdb-file-name (append-dir-file "coot-download" (string-append code ".pdb")))
	 (cif-file-name (append-dir-file "coot-download" (string-append "r" code "sf.cif"))))
    (let* ((coords-type ".pdb") ;; can/will be ".cif"
	   (pdb-url (string-append 
		     "http://www.ebi.ac.uk/pdbe-srv/view/files/"
		     code coords-type))
	   (cif-url (string-append
		     "http://www.ebi.ac.uk/pdbe-srv/view/files/r"
		     code "sf.ent")))
      (let ((status-1 #t)
	    (status-2 #t))
	(if (not (file-exists? pdb-file-name))
	    (set! status-1 (coot-get-url pdb-url pdb-file-name)))
	(if (not (file-exists? cif-file-name))
	    (set! status-2 (coot-get-url cif-url cif-file-name)))
	(values status-1 status-2)))))

;; only look at ligands in maps with resolution worse than this:
;; (to look at everything, set it to 0.1)
;; 
(define *ligand-check-resolution-limit* 2.5)

(define (do-everything)

  (define (handle-pdb-entry-entity e)
    ;; (format #t "::: e: ~s~%" e)
    (if (hash-table? e)
	(let* ((entry-id   (hash-ref e "EntryID"))
	       (title      (hash-ref e "Title"))
	       (contains-ligands (hash-ref e "ContainsLigands"))
	       (resolution (hash-ref e "Resolution"))
	       (method     (hash-ref e "Method")))
	  (if (and (string? entry-id)
		   (> (length contains-ligands) 0)
		   (and (string? method)
			(string=? method "x-ray diffraction")
			(if (not (string? resolution))
			    #f
			    (let ((n (string->number resolution)))
			      (if (not (number? n))
				  #f
				  (>= n *ligand-check-resolution-limit*))))))
	      entry-id
	      #f))))


  (let ((url "http://www.ebi.ac.uk/pdbe-apps/jsonizer/latest/released/")
	(json-file-name "latest-releases.json"))
    (if (not (file-exists? json-file-name))
	(coot-get-url url json-file-name))

    (let ((t (json-object-from-file json-file-name)))
      ;; main line of recent-structure-browser
      (if (not (hash-table? t))
	  (begin
	    #f ;; in silence please
	    )

	  (begin
	    (let ((a (hash-ref t "ResultSet")))
	      (if (not (hash-table? a))
	      #f
	      (let ((b (hash-ref a "Result")))
		(if (not (list? b))
		    #f
		    (let ((codes-from-xrays-with-ligands 
			   (filter string? (map handle-pdb-entry-entity b))))
		      (format #t "codes-from-xrays-with-ligands: ~s~%" 
			      codes-from-xrays-with-ligands)
		      
		      (for-each (lambda (code)
				  (receive (status-1 status-2)
					   (get-pdbe-files-for-code code)
					   (if (and status-1 status-2)
					       (get-correlation-of-ligand code))))
				codes-from-xrays-with-ligands)))))))))))



;; get files (mtz and pdb) for each of the accession codes
;; 
;; put them in coot-download
;;
;; copy mtz files from fei, and get pdb files from pdbe
;; 
(define (get-download-files-for-code code fei-mtz-dir)

  (let* ((pdb-file-name (append-dir-file "coot-download" (string-append code ".pdb")))
	 (coords-type ".pdb") ;; can/will be ".cif"
	 (coords-ext  ".ent")
	 (pdb-url (string-append 
		   "http://www.ebi.ac.uk/pdbe/entry-files/download/r/"
		   code coords-ext))
	 (sfs-mtz-file-name (append-dir-file "coot-download" 
					     (string-append "r" code "sfs.mtz")))
	 (cif-file-name (append-dir-file "coot-download" (string-append "r" code 
									"sf.cif")))
	 (cif-url (string-append "http://www.ebi.ac.uk/pdbe-srv/view/files/r"
				 code "sf.ent")))

  (if (not (file-exists? cif-file-name))
      (begin
	(format #t "getting cif-file-name: ~s~%" cif-file-name)
	(coot-get-url cif-url cif-file-name)))


  (if (not (file-exists? pdb-file-name))
      (let ((status-1 (coot-get-url pdb-url pdb-file-name)))
	(if (not status-1)
	    #f
	    (let ((fei-mtz-file-name 
		   (append-dir-file fei-mtz-dir
				    (string-append code ".mtz"))))
	      (if (not (file-exists? fei-mtz-file-name))
		  #f
		  (begin
                    (format #t "copying file ~s ~s ~%" fei-mtz-file-name sfs-mtz-file-name)
		    (copy-file fei-mtz-file-name sfs-mtz-file-name)))))))))




(define (process-accession-codes code-list results-out-file-name)

  (define fei-mtz-dir "/teraraid/garibgrp/database/mtz")

  (define (print-results results port)
    (for-each (lambda (result)
		(if (list? result)
		    (let* ((mogul-z-list (list-ref result 3))
			   (mogul-score (if (not (list? mogul-z-list))
					    mogul-z-list
					    (apply max mogul-z-list)))
			  (cs (list-ref result 2)))
		      (format port 
			      "               code: ~s ~s n-atoms: ~s spec: ~s contact-score: ~s ~s ~s ~s cor: ~6f mogul: ~s med-ratio ~s ks-test: ~s ~%" 
			      (list-ref result 0) ;; accession code
			      (list-ref result 6) ;; ligand comp-id
			      (het-group-n-atoms (list-ref result 6))
			      (list-ref result 5)  ;; ligand spec
			      (list-ref cs 0)  ;; clash score short contacts
			      (list-ref cs 1)  ;; clash score h-bonds 
			      (list-ref cs 2)  ;; clash score ... 
			      (list-ref cs 3)  ;; clash score ... 
			      (list-ref result 1)  ;; correlation
			      mogul-score ;; mogul-score (worst z)
			      (list-ref result 4) ;; b-factor median ratio
			      (list-ref result 7) ;; ks distance
			      ))
		    (format port "          bad result ~s~%" result)))
	      results))

  ;; main line
  ;;
  ;; get files for each of the accession codes
  (map (lambda (code)
	 (get-download-files-for-code code fei-mtz-dir))
       code-list)

  (let ((results (map get-correlation-of-ligand code-list)))
    (print-results results #t)
    (if (string? results-out-file-name)
	(call-with-output-file results-out-file-name
	  (lambda (port)
	    (print-results results port))))))

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


;; return a list:
;; (list list-of-bs-in-ligand (list list-of-bs-in-env-residue))
;; So to get the mean of the environment, you will need to do an append to the cadr.
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

(define (median number-list)

  (let ((numbers 
	 (let loop ((number-list number-list)
		    (nums '()))
	   (cond
	    ((null? number-list) nums)
	    ((not (number? (car number-list)))
	     (loop (cdr number-list) nums))
	    (else 
	     (loop (cdr number-list) (cons (car number-list) nums)))))))
    (let* ((sorted-nums (sort-list numbers >))
	   (n (length sorted-nums))
	   (mid (/ n 2))
	   (median (if (= (remainder n 2) 0)
		       (begin
			 ;; (format #t "averaging mid: ~s  n: ~s ~%" mid n)
			 ;; (format #t "idx: ~s  idx: ~s ~%" (/ n 2) (- (/ n 2) 1))
			 (/ (+ (list-ref sorted-nums (/ n 2))
			       (list-ref sorted-nums (- (/ n 2) 1))) 2))
		       (begin
			 ;; (format #t "simple take: ~s~%" mid)
			 (list-ref sorted-nums (/ (- n 1) 2))))))
      median)))
	   

(define (median-ratio ligand-b-factors env-b-factors)

  (/ (median ligand-b-factors) (median (apply append env-b-factors))))

  


;; test
(define (test-thing) 
  ;; not  "4feg", because, as yet, we can't run refmac with mtz file only intensities
  (let ((code-list
	 ;; (list "3vff" "3tuc" "3t64" "3tud" "3vk2" "4eg6")
	 ;; (list "3vk2")
	 ;; (list "3vk2")
	 ;; (list "3t64")
	 ;; "4ekg" no dictionary
	 ;; (list "3tyg"  "4b5m" "4gdc" "4gns" "4gdd" "4b5g" "4hfr")
	 ;;  (list "4gdc")
	 ;; (list "3s05") ;; caused a crash.
	 ;; (list "3u08")
	 ;; (list "3rfe") ;; not a great build
	 ;; (list "3rkr")
	 ;; (list "3tfj")
	 (list "3rmz")
	 ))
    (process-accession-codes code-list #f)))

;; (define this-weeks-file  "/lmb/home/pemsley/fei/All_pdb/pdb_Dec_1_2010.list")
;; (define this-weeks-file  "/lmb/home/pemsley/fei/tiny-sample-bit")
(define this-weeks-file  "/lmb/home/pemsley/fei/sample-bit")

(define (fei-test)
  (if (not (file-exists? this-weeks-file))
      (begin
	(format #t "file does not exist: ~s ~%" this-weeks-file)
	#f)
      (let ((code-list 
	     (call-with-input-file this-weeks-file
	       (lambda (port)
		 (let loop ((code (read-line port))
			    (codes '()))
		   (if (not (eof-object? code))
		       (loop (read-line port) (cons code codes))
		       codes)))))
	    (results-file-name 
	     (string-append "ligand-of-the-week-"
			    (strip-path this-weeks-file)
			    ".tab")))

	(format #t "code list: ~s ~%" code-list)
	(process-accession-codes code-list results-file-name))))
	
