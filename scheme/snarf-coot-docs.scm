
(use-modules (ice-9 popen)
             (ice-9 string-fun)
	     (ice-9 regex))

(use-modules (os process))

(define is-valid-model-molecule 1) ;; hack to get past problem loading
				   ;; coot-utils on bubbles

(define generic-object-name-scm #f)
(define additional-representation-info-scm #f)
(define missing-atom-info-scm #f)
(define drag-intermediate-atom-scm #t)
(define mark-atom-as-fixed-scm #t)
(define ncs-chain-ids-scm #t)
(define ncs-chain-differences-scm #t)
(define refmac-parameters-scm #t)
(define water-chain-scm #f)
(define water-chain-from-shelx-ins-scm #f)
(define key-sym-code-scm #f)
(define ncs-ghosts-scm #f)
(define inverse-rtop-scm #f)
(define coot-has-python-p #f)
(define map-sigma-scm #f)
(define pucker-info-scm #f)
(define map-parameters-scm #f)
(define map-cell-scm #f)
(define cell-scm #f)
(define ccp4i-projects-scm #f)
(define get-rotamer-name-scm #f)
(define test-internal-scm #f)
(define atom-info-string-scm #f)
(define get-refmac-sad-atom-info-scm #f)
(define set-find-hydrogen-torsions #f)
(define find-terminal-residue-type #f)
(define residues-near-position-scm #f)
(define non-standard-residue-names-scm #f)
(define refine-residues-scm #f)
(define refine-residues-with-alt-conf-scm #f)
(define regularize-residues-scm #f)
(define map-peaks-scm #f)
(define map-peaks-near-point-scm #f)
(define add-dipole-scm #f)
(define add-dipole-for-residues-scm #f)
(define get-torsion-scm #f)
(define set-torsion-scm #f)
(define test-internal-single-scm #f)
(define user-defined-click-scm #f)
(define add-lsq-atom-pair-scm #f)
(define coot-sys-build-type-scm #f)
(define add-alt-conf-scm #f)
(define origin-pre-shift-scm #f)
(define alignment-mismatches-scm #f)
(define rigid-body-refine-by-residue-ranges-scm #f)
(define average-map-scm #f)
(define symmetry-operators-scm #f)
(define symmetry-operators-to-xHM-scm #f)
(define user-mods-scm #f)
(define refine-zone-with-full-residue-spec-scm #f)
(define get-pkgdatadir-scm #f)
(define matching-compound-names-from-sbase-scm #f)
(define chain-id-scm #f)
(define highly-coordinated-waters-scm #f)
(define handle-pisa-interfaces-scm #f)
(define space-group-scm #f)
(define nearest-residue-by-sequence-scm #f)
(define list-extra-restraints-scm #f)
(define delete-extra-restraint-scm #f)
(define do-clipped-surface-sc #f)

(load "filter.scm")
(load "coot-utils.scm")

(define is-empty-file?
  (lambda (file-name)

    (let ((strs (run-command/strings "wc" (list "-l" file-name) '())))
      (if (list? strs)
	  (let ((ls (string->list-of-strings (car strs))))
	    (if (not (null? ls))
		(let ((n (string->number (car ls))))
		  (if (number? n)
		      (= 0 n)
		      #f))
		#f))
	  #f))))



;; texi2html on user-manual.texi with coot-scheme-functions.texi included gives:
;; 
;; *** Duplicate node found: mutate (in ./../../coot/scheme/coot-scheme-functions.texi l. 42)
;; *** Duplicate node found: povray (in ./../../coot/scheme/coot-scheme-functions.texi l. 341)
;;*** Duplicate node found: raster3d (in ./../../coot/scheme/coot-scheme-functions.texi l. 361)
;; 
(let* ((source-files (glob "*.scm" ".")))

  (map (lambda (source-file)
     
     (let* ((doc-file
	     (string-append 
	      (strip-extension source-file) ".texi"))
	    (tmp-file-1 (string-append doc-file ".tmp")))
       
       ;; problems snarfing coot docs?  check that the function-name
       ;; and args are define all on one line.
       (let ((args (list "doc-snarf"
			    "--texinfo"
			    "-o"
			    tmp-file-1
			    source-file)))
	 (goosh-command "guile-tools" args '() "snarf-log" #t)
	 (goosh-command "grep" (list "-v" "^" tmp-file-1) '() doc-file #f))

       (if (file-exists? tmp-file-1)
	   (begin
	     (delete-file tmp-file-1))
	   (begin
	     (format #t "OOPS: Problem snarfing ~s~%" source-file)
	     (format #t "OOPS: ~s does not exist - not deleting~%" tmp-file-1)))))

       source-files))
  
(define add-to-list-section-texis
  (lambda (texi-list)

    (let f ((ls texi-list)
	    (new-list '()))
      
      (cond 
       ((null? ls) new-list)
       (else 
	(let* ((stub (strip-extension (car ls)))
	       (section-texi-name (string-append stub "-section.texi")))
	  (f (cdr ls)
	     (append (list section-texi-name (car ls)) new-list))))))))

(define delete-section-texi-files
  (lambda ()

    (let ((section-texi-files (glob "*-section.texi" ".")))
      (map delete-file section-texi-files))))
    

(delete-section-texi-files)
(let ((all-raw-texi-files (glob "*.texi" ".")))
  
  (let ((non-empty-texis
	 (let f ((ls all-raw-texi-files)
		 (non-null-list '()))
	   
	   (cond 
	    ((null? ls) non-null-list)
	    ;; we also want to reject coot-scheme.texi, the master
	    ;; texi file (we don't want to cat that) and the catted file
	    ((string=? (car ls) "coot-scheme.texi")
	     (f (cdr ls) non-null-list))
	    ((string=? (car ls) "coot-scheme-functions.texi")
	     (f (cdr ls) non-null-list))
	    ((is-empty-file? (car ls)) 
	     (delete-file (car ls))
	     (f (cdr ls) non-null-list))
	    (else 
	     (f (cdr ls) (cons (car ls) non-null-list)))))))

    ;; we want to add section/node info before each of the texi files
    ;; so let's create a section-texi file that has that info in.
    ;;
    ;; make the section texis here
    (call-with-output-file "menu-section.texi"
      (lambda (menu-port)
	(format menu-port "@menu~%")
	(for-each (lambda (file)
		    (format #t "processing ~s~%" file)
		    (let* ((stub (strip-extension file))
			   (node-name (cond
				  ((string=? stub "mutate") "mutate-in-scheme")
				  ((string=? stub "povray") "mutate-from-scheme")
				  ((string=? stub "raster3d") "raster3d-from-scheme")
				  (else
				   stub)))
			   (section-name (string-append stub "-section.texi")))

		      (format menu-port "* ~a::~%" node-name)
		      
		      (call-with-output-file section-name
			(lambda (port)
			  
			  (format port "@node    ~a ~%" node-name)
			  (format port "@section ~a ~%" node-name)
			  ;; (format port "@cindex  ~a ~%" stub) not cindex
			  ))))
		  non-empty-texis)
	(format menu-port "@end menu~%")))
	
    
    (let ((all-useful-texis (cons "menu-section.texi"
				  (add-to-list-section-texis non-empty-texis))))

      (format #f "all-useful-texis: ~s~%" all-useful-texis)
      (goosh-command "cat" all-useful-texis '() "coot-scheme-functions.texi" #f))))


