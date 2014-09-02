
;;;; Copyright 2006 Joel Bard
;;;; Copyright 2006 Paul Emsley
    
  
;; Read in cns coeff-data (given filenames) and a pdb molecule filename to make
;; maps.
;; 
(define (cns->coot 2fofc-coeffs fofc-coeffs model-pdb)


  ;; remove any trailing spaces of string
  (define (trim-trailing-spaces s)
    
    (let f ((ls (reverse (string->list s))))
      
      (cond 
       ((null? ls) "")
       ((eq? #\space (car ls)) (f (cdr ls)))
       (else
	(list->string (reverse ls))))))


  ;; Return '() on bad, else return a 6 membered list.
  ;; 
  (define (get-cell-from-pdb pdb-file)
    
    (let ((s (run-command/strings "awk" (list "/CRYST1/ {printf \"%.3f %.3f %.3f %.3f %.3f %.3f\\n\", $2, $3, $4, $5, $6, $7}" pdb-file) '())))
      
      (if (not (= (length s) 1))
	  '()
	  (car s))))

  ;; return "" on bad, else return a string that is the space group.
  ;; 
  ;; #have to go through all sorts of hoops because sometimes we might
  ;; #get P 21 and sometimes we might get P 1 21 1 for example
  ;; #sftools needs P21, can't handle P1211
  ;; #that being said - there must be a better way
  ;; set SYMM1=`awk 'BEGIN { FIELDWIDTHS="55 11 4" } /CRYST1/ {gsub(" ","",$2); printf "%s\n", $2}' $3`
  ;; 
  (define (get-symm-from-pdb-file pdb-file)
    
;    (let ((s (run-command/strings "awk" (list "BEGIN { FIELDWIDTHS=\"55 11 4\" } /CRYST1/ {gsub(\" \",\"\",$2); printf \"%s\\n\", $2}" pdb-file) '())))

    (let* ((a-str "/CRYST1/ {s=substr($0,56,11);gsub(\" \",\"\",s); print s}")
	   (s (run-command/strings "awk" (list a-str pdb-file) '())))
      (if (not (= (length s) 1))
	  ""
	  (car s))))

  (define (get-symm-from-pdb-file-2 pdb-file)
    
    (let* ((a-str "/CRYST1/ {s=substr($0,56,11);gsub(\" \",\"\",s); print s}")
	  (s (run-command/strings "awk" (list a-str pdb-file) '())))
      (if (not (= (length s) 1))
	  ""
	  ;; because this is not shell sript, we have to be careful
	  ;; to remove trailing spaces in the variable.  Csh does
	  ;; this implicitly.
	  (trim-trailing-spaces (car s)))))

  ;; Return a symmetry string.
  ;; on a problem, return a #f.  Use that to exit (elsewhere).
  ;; 
  (define (get-symm-from-symop-lib symm2)

    (let ((e (getenv "CCP4")))
      (if (not (string? e))
	  (begin 
	    (format #t "CCP4 not set up.~%")
	    #f)
	  (let ((s (run-command/strings "awk" (list (string-append "/" symm2 "/ { print $4}")
						    (append-dir-file 
						     (append-dir-dir (append-dir-dir e "lib") "data")
						     "symop.lib")) 
					'())))
	    (if (= (length s) 0)
		#f ;; no hits, oops
		(begin 
		  (car s))))))) ;; pick the first one ! (yikes?)

    ;; main body
    (let* ((map1-prefix (strip-extension 2fofc-coeffs))
	   (map2-prefix (strip-extension  fofc-coeffs))
	   (map1-tmp (string-append map1-prefix "_tmp.pdb"))
	   (map2-tmp (string-append map2-prefix "_tmp.pdb"))
	   (cell (get-cell-from-pdb model-pdb)))
      
      (format #t "map1-prefix: ~s~%" map1-prefix)
      (format #t "map2-prefix: ~s~%" map2-prefix)
      (format #t "cell: ~s~%" cell)
      
      (let ((symm1 (get-symm-from-pdb-file model-pdb))
	    (pdbset-log "cns2coot-pdbset-tmp.log"))
	
	(format #t "symm1: ~s~%" symm1)
	(goosh-command "pdbset" (list "XYZIN" model-pdb "XYZOUT" map1-tmp) 
		       (list (string-append "CELL " cell)
			     (string-append "SPACEGROUP " symm1))
		       pdbset-log #t)
	(let* ((symm2 (get-symm-from-pdb-file-2 map1-tmp))
	       (nov (format #t "symm2: ~s~%" symm2))
	       (symm (get-symm-from-symop-lib symm2)))
	  
	  (if (eq? #f symm)
	      (format #t "Failed to find symm in symop.lib~%")
	      (begin 
		(format #t "INFO:: SYMM is ~s~%" symm)
		(let ((map1-mtz (string-append map1-prefix ".mtz"))
		      (map2-mtz (string-append map2-prefix ".mtz"))
		      (map-coot-mtz (string-append map1-prefix "-coot.mtz"))
		      (cad-log (string-append "cns2coot-cad-tmp.log"))
		      (sftools-1-log  (string-append "cns2coot-sftools-1-tmp.log"))
		      (sftools-2-log  (string-append "cns2coot-sftools-2-tmp.log")))
		  
		  (if (file-exists? map1-mtz)
		      (delete-file map1-mtz))
		  (if (file-exists? map2-mtz)
		      (delete-file map2-mtz))
		  
		  (goosh-command "sftools" '() (list (string-append "read " 2fofc-coeffs)
						     "cns"
						     cell
						     symm
						     "END"
						     "W"
						     "P"
						     "R"
						     "SET LABELS"
						     "FOM"
						     "PHIC"
						     "SCALE"
						     "FWT"
						     "PHWT"
						     (string-append "WRITE " map1-mtz))
				 sftools-1-log #t)
		  
		  (goosh-command "sftools" '() (list (string-append "read " fofc-coeffs)
						     "cns"
						     cell
						     symm
						     "END"
						     "W"
						     "P"
						     "R"
						     "SET LABELS"
						     "FOM"
						     "PHIC"
						     "SCALE"
						     "DELFWT"
						     "PHDELWT"
						     (string-append "WRITE " map2-mtz))
				 sftools-2-log #t)
		  
		  (goosh-command "cad" (list "HKLIN1" map1-mtz 
					     "HKLIN2" map2-mtz 
					     "HKLOUT" map-coot-mtz)
				 (list "LABIN FILE_NUMBER 1 E1=FOM E2=PHIC E3=FWT E4=PHWT"
				       "LABIN FILE_NUMBER 2 E1=DELFWT E2=PHDELWT"
				       "END") cad-log #t)
		  
		  (if (file-exists? map1-mtz)
		      (delete-file map1-mtz))
		  (if (file-exists? map2-mtz)
		      (delete-file map2-mtz))

		  ;; now load them into Coot
		  ;; 
		  (read-pdb model-pdb)
		  (make-and-draw-map-with-reso-with-refmac-params map-coot-mtz "FWT" "PHWT" "" 0 0 0 
								  "Fobs:None-specified" "SigF:None-specified" 
								  "RFree:None-specified" 0 0 0 -1.00 -1.00)

		  (make-and-draw-map-with-reso-with-refmac-params map-coot-mtz "DELFWT" "PHDELWT" "" 0 1 0
								  "Fobs:None-specified" "SigF:None-specified" 
								  "RFree:None-specified" 0 0 0 -1.00 -1.00)
		  )))))))
    

