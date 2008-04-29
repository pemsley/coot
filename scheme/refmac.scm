;;;; refmac.scm
;;;;
;;;; Copyright 2004, 2005, 2006, 2007 by The University of York
;;;; Author: Paul Emsley
;;;; 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 2, or (at your option)
;;;; any later version.
;;;; 
;;;; This program is distributed in the hope that it will be useful,
;;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;;; GNU General Public License for more details.
;;;; 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this software; see the file COPYING.  If not, write to
;;;; the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
;;;; Boston, MA 02111-1307 USA
;;;;

;;; Extra parameters can be passed to refmac using either a file
;;; "refmac-extra-params" or by setting the variable
;;; refmac-extra-params.  The LABIN line should not be part those
;;; extra parameters, of course - Coot takes care of that.

;; This is the default refmac version, it is presumed to be in the
;; path.  It can be overridden using a re-definition either at the
;; scripting interface or in one's ~/.coot file. E.g.:
;; (define refmac-exec "/y/programs/xtal/refmac-latest/bin/refmac5-3-dec-2004")
(define refmac-exe "refmac5") 

;; Set this to a list of parameter strings:
;; 
;; If refmac-extra-params is a list of strings, it is used in
;; preference to the "refmac-extra-params" file (should it exist).
;; e.g. (set! refmac-extra-params (list "WEIGHT 0.2" 
;;                                      "NCYC 10" 
;;                                      "REFI BREF ISO"
;;                                      "REFI METH CGMAT"
;;                                      "REFI TYPE REST RESO 20 1.64"))
;; 
(define refmac-extra-params #f)

(define refmac-count 0)

;; If phase-combine-flag is 1 then we should do a phase combination
;; with the coefficients that were used to make the map, specifically
;; PHIB and FOM, that are the result of (say) a DM run.  It is best
;; of course to use HL coefficients, if DM (say) wrote them out.
;; 
;; If phase-combine-flag is 0, then phib-fom-pair is ignored.  If
;; phase-combine-flag is 1 phib-fom-pair is presumed to be an
;; improper list of the PHIB and FOM column labels (as strings).
;; 
;; Is there a case where you want phase recombination, but do not
;; have FOMs?  Don't know, I presume not.  So we required that if
;; phase-combine-flag is 1, then PHIB and FOM are present.  Else we
;; don't run refmac.
;; 
;; If ccp4i project dir is "" then we ignore it and write the refmac
;; log file here.  If it is set to something, it is the prefix of the
;; refmac log file.  In that case it should end in a "/".  This is
;; not tested for here!  The pdb-in-filename etc are dealt with in
;; the c++ execute_refmac function.
;; 
;; if swap_map_colours_post_refmac? is not 1 then imol-mtz-molecule
;; is ignored.
;;
;; It is possible that we get here with long refmac column labels,
;; e.g.  "/RNASE3GMP/COMPLEX/FGMP18" "/RNASE3GMP/COMPLEX/SIGFGMP18"
;; "/RNASE/NATIVE/FreeR_flag".
;; 
;; So we only use text after the last /.
;;
;; 
(define (run-refmac-by-filename pdb-in-filename pdb-out-filename mtz-in-filename mtz-out-filename extra-cif-lib-filename imol-refmac-count swap-map-colours-post-refmac? imol-mtz-molecule show-diff-map-flag phase-combine-flag phib-fom-pair force-n-cycles ccp4i-project-dir f-col sig-f-col . r-free-col) 

    (format #t "got args: pdb-in-filename: ~s, pdb-out-filename: ~s, mtz-in-filename: ~s, mtz-out-filename: ~s, imol-refmac-count: ~s, show-diff-map-flag: ~s, phase-combine-flag: ~s, phib-fom-pair: ~s, force-n-cycles: ~s, f-col: ~s, sig-f-col: ~s, r-free-col: ~s~%"
	    pdb-in-filename pdb-out-filename
	    mtz-in-filename mtz-out-filename
	    imol-refmac-count show-diff-map-flag
	    phase-combine-flag
	    phib-fom-pair
	    force-n-cycles
	    f-col sig-f-col r-free-col)
	
    (let* ((local-r-free-col (if (null? r-free-col) '() (car r-free-col)))
;	   (labin-string
	   (labin-string (if (string=? f-col "") ""
			     (apply string-append (append (list "LABIN" " "
								"FP=" (strip-path f-col) " "
								"SIGFP=" (strip-path sig-f-col))
							  (if (null? local-r-free-col) 
							      '()
;					     (list " FREE=" (strip-path local-r-free-col))))))
							      (list " FREE=" (strip-path local-r-free-col)))))))
	  (command-line-args
	   (append
	    (list 
	     "XYZIN"  pdb-in-filename
	     "XYZOUT" pdb-out-filename
	     "HKLIN"  mtz-in-filename
	     "HKLOUT" mtz-out-filename)
	    (if (string=? extra-cif-lib-filename "")
		(begin
		  (format #t "Not Passing LIBIN to refmac LIBIN~%")
		  (list)) ; nothing
		(begin
		  (format #t "Passing to refmac LIBIN ~s~%" extra-cif-lib-filename)
		  (list "LIBIN" extra-cif-lib-filename)))))
      
	  (data-lines (let* ((std-lines
			      (list 
			       "MAKE HYDROGENS NO" ; Garib's suggestion 8 Sept 2003
			       (if (number? force-n-cycles)
				   (if (>= force-n-cycles 0)
				       (string-append
					"NCYCLES " 
					(number->string force-n-cycles))
				       "")
				   "")
			       ))
			     (extra-params (get-refmac-extra-params)))

			(if (extra-params-include-weight? extra-params)
			    (append std-lines
				    extra-params
				    (list labin-string))
			    (append std-lines
				    (list "WEIGHT AUTO")
				    extra-params
				    (list labin-string)))))

	  (nov (format #t "DEBUG:: refmac-extra-params returns ~s~%"
		       (get-refmac-extra-params)))
	  ;; this should be a database filename:
	  ;; 
	  (refmac-log-file-name (string-append 
				 (if (> (string-length ccp4i-project-dir) 0)
				     ccp4i-project-dir
				     "")
				 "refmac-from-coot-" 
				 (number->string refmac-count) ".log")))

      (set! refmac-count (+ imol-refmac-count 1))
      (format #t "INFO:: Running refmac with these command line args: ~s~%"
	      command-line-args)
      (format #t "INFO:: Running refmac with these data lines: ~s~%"
	      data-lines)
      (format #t "environment variable:  SYMOP: ~s~%" (getenv "SYMOP"))
      (format #t "environment variable: ATOMSF: ~s~%" (getenv "ATOMSF"))
      (format #t "environment variable:  CLIBD: ~s~%" (getenv "CLIBD"))
      (format #t "environment variable:   CLIB: ~s~%" (getenv "CLIB"))

      ;; first check if refmac exists?
      (format #t "INFO:: now checking for refmac exe - [it should give status 1...]~%" )
      (let ((refmac-status  (goosh-command refmac-exe '() (list "END")
					   refmac-log-file-name #f)))

	;; potentially a no-clobber problem here, I think.

	(format #t "INFO:: refmac-status: ~s~%" refmac-status)

	(if (not (number? refmac-status))
	    -3
	    (if (not (= refmac-status 1))
		
		;; problem finding refmac executable
		(begin 
		  (format #t "refmac failed (no executable)")
		  (format #t " - no new map and molecule available~%"))

		;; OK, we found the executable, this should be OK then...
		(let ((status (goosh-command refmac-exe 
					     command-line-args 
					     data-lines 
					     refmac-log-file-name #t))) ; to screen too
		  
		  (if (file-exists? refmac-log-file-name)
		      (run-concurrently "loggraph" refmac-log-file-name))
		  
		  (if (and (number? status) (= status 0)) ; refmac ran OK
		      (begin
					; now let's read in those
					; newly-created coordinates
					; and phases for a map:
			
			(let* ((r-free-bit (if (null? r-free-col)
					       (list  "" 0)
					       (list local-r-free-col 1)))
			       (args 
				(append
					; numbers: use-weights? is-diff-map? have-refmac-params?
				 (list mtz-out-filename "FWT" "PHWT" "" 0 0 1 f-col sig-f-col)
				 r-free-bit))
			       (recentre-status (recentre-on-read-pdb))
			       (novalue (set-recentre-on-read-pdb 0))
			       (novalue2 (format #t "DEBUG:: recentre status: ~s~%" recentre-status))
			       (imol (handle-read-draw-molecule pdb-out-filename)))
			  
			  (if recentre-status
			      (set-recentre-on-read-pdb 1))
			  (set-refmac-counter imol (+ imol-refmac-count 1))
			  
			  (let ((new-map-id (apply make-and-draw-map-with-refmac-params args)))
			    
			    (if (= swap-map-colours-post-refmac? 1)
				(swap-map-colours imol-mtz-molecule new-map-id)))
			  
			  (if (= 1 show-diff-map-flag) ; flag was set
			      (let ((args (append
					   (list mtz-out-filename "DELFWT" "PHDELWT" "" 0 1 1
						 f-col sig-f-col)
					   r-free-bit)))
				
				(apply make-and-draw-map-with-refmac-params args)))))

		      (begin 
			(format #t "refmac failed - no new map and molecule available~%")))))))))

;; Return #t if the list of strings @var{params-list} contains a
;; string beginning with "WEIGHT".  If not return #f
;; 
(define (extra-params-include-weight? params-list)

  ;; This tests all words in the string.  We should more rigourously
  ;; test only the string either at the beginning or after "".
  ;; 
  (define has-weight-word?
    (lambda (string-list)
      
      (cond
       ((null? string-list) #f)
       ((string-ci=? (car string-list) "WEIGHT") #t)
       (else 
	(has-weight-word? (cdr string-list))))))
  
  ;;main body
  (if (not (list? params-list))
      #f
      (let f ((params-list params-list))

	(cond
	 ((null? params-list) #f)
	 ((has-weight-word? (string->list-of-strings (car params-list))) #t)
	 (else 
	  (f (cdr params-list)))))))


;; If refmac-extra-params is defined (as a list of strings), then
;; return that, else read the file "refmac-extra-params".
;; 
;; Return a list a list of strings.
;; 
(define get-refmac-extra-params
  (lambda ()

    (if (list-of-strings? refmac-extra-params)

	refmac-extra-params

	(let ((extras-file-name
	       "refmac-extra-params")) 
	  
	  (if (not (file-exists? extras-file-name))
	      '()
	      
	      (call-with-input-file extras-file-name
		(lambda (port)
		  
		  (let ((r-list 
			 (let f ((line-list '())
				 (line (read-line port)))
			   
			   ; (format #t "line-list: ~s~%" line-list)
			   
			   (cond
			    ((null? line) line-list)
			    ((eof-object? line) line-list)
			    (else (f (cons line line-list)
				     (read-line port)))))))
		    
		    (reverse r-list)))))))))

;; 
(define run-refmac-for-phases
  (lambda (imol mtz-file-name f-col sig-f-col)
    
    (if (file-exists? mtz-file-name)
	(if (valid-model-molecule? imol)
	      (let ((dir-state (make-directory-maybe "coot-refmac")))
		(if (not (= 0))
		    (format #t "Failed to make coot-refmac directory\n")
		    (let* ((stub (string-append "coot-refmac/refmac-for-phases"))
			   (pdb-in  (string-append stub ".pdb"))
			   (pdb-out (string-append stub "-tmp.pdb"))
			   (mtz-out (string-append stub ".mtz"))
			   (cif-lib-filename ""))
		      (write-pdb-file imol pdb-in)
		      (run-refmac-by-filename pdb-in pdb-out
					      mtz-file-name mtz-out
					      cif-lib-filename 0 0 -1
					      1 0 '() 0 "" 
					      f-col sig-f-col))))))))
	

(define refmac-for-phases-and-make-map
  (lambda (mtz-file-name f-col sig-f-col)

    (if (file-exists? mtz-file-name)
	(molecule-chooser-gui 
	 "  Choose a molecule from which to calculate Structure factors:  "  
	 ; a lambda function that accepts the choose imol as its arg:
	 (lambda (imol)
	   (run-refmac-for-phases imol mtz-file-name f-col sig-f-col))))))
