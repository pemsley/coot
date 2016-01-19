
;; This is a wrapper for the python function ligand_validation_metrics_gui() 
;; in python/ligand_validation_sliders.py.

(define (refmac-columns f-list sigf-list rfree-list)

  (if (not (> (vector-length f-list) 0))
      #f
      (if (not (> (vector-length sigf-list) 0))
	  #f
	  (if (not (> (vector-length rfree-list) 0))
	      #f
	      
	      ;; happy path
	      
	      ;; Using the first sigf (there should only be one typically)
	      ;; get the F label (Fx) from x/y/SIGFPx
	      ;; 
	      (let ((field (vector-ref sigf-list 0)))
		(format #t "debug:: f-list: ~s~%" f-list)
		(format #t "debug:: sigf-list: ~s~%" sigf-list)
		(format #t "debug:: field: ~s~%" field)
		(let ((ls (split-after-char-last #\/ field list)))
		  (format #t "debug:: ls: ~s~%" ls)
		  (let ((sigf-col (cadr ls)))
		    (format #t "debug:: sigf-col: ~s~%" sigf-col)
		    (let ((sm (string-match "[Ss][Ii][Gg]F" sigf-col)))
		      (format #t "debug:: sm: ~s~%" sm)
		      ;; (let ((f-col (substring sigf-col (cdr (vector-ref sm 1))))) ;; what was I thinking?

		      (let ((f-col (vector-ref f-list 0))) ;; hacketty hack!
			   
			(format #t "debug:: f-col: ~s~%" f-col)
			(list 
			 (string-append (car ls) "F" f-col)
			 field
			 (vector-ref rfree-list 0)))))))))))


(define (get-refmac-columns-for-map imol-map)
  (let ((m (mtz-file-name imol-map)))
    (format #t "refmac-mtz-filename: ~s~%" m)
    (if (not (> (string-length m) 0))
	#f
	(let ((f-cols (get-f-cols m))
	      (sigf-cols (get-sigf-cols m))
	      (free-cols (get-r-free-cols m)))
	  (format #t "f      columns: ~s~%" f-cols)
	  (format #t "sigf   columns: ~s~%" sigf-cols)
	  (format #t "r-free columns: ~s~%" free-cols)
	  (let ((rc (refmac-columns f-cols sigf-cols free-cols)))
	    (format #t "---------- here in get-refmac-columns-for-map #1 ----------~%~!")
	    (print-var rc)
	    rc)))))



(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Ligand")))

      (add-simple-coot-menu-menuitem
       menu "Ligand Validation Tests"
       (lambda ()
	 (using-active-atom

	  (format #t "---------- here #1 ----------~%~!")

	  (let ((imol-map (imol-refinement-map)))

	    (let ((refmac-input-mtz-file-name (mtz-file-name imol-map))
		  (refmac-cols (get-refmac-columns-for-map imol-map)))

	      (format #t "---------- here #2 ----------~%~!")
	      (print-var refmac-cols)

	      (if (not (list? refmac-cols))
		  
		  (begin
		    (format #t "Couldn't get columns from mtz file: ~s~%~!" refmac-input-mtz-file-name))

		  ;; happy path
		  (let ((fp-col    (list-ref refmac-cols 0))
			(sigfp-col (list-ref refmac-cols 1))
			(rfree-col (list-ref refmac-cols 2)))
		    
		    (let ((refmac-dir (get-directory "coot-refmac")))

		      (format #t ":::::::::: (get-metrics-for-ligand ~s ~s ~s ~s ~s ~s ~s ~s ~s)~%~!"
			      aa-imol aa-chain-id aa-res-no aa-ins-code 
			      refmac-input-mtz-file-name 
			      fp-col sigfp-col rfree-col refmac-dir)

		      (let ((m (get-metrics-for-ligand aa-imol aa-chain-id aa-res-no aa-ins-code 
						       refmac-input-mtz-file-name 
						       fp-col sigfp-col rfree-col refmac-dir)))
			(format #t "get-metrics-for-ligand returns ~s~%~!" m)
			(if (list? m)
			    't
			    ))))))))))))




      
		    


      

