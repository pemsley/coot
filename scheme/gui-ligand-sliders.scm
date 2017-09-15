;;;; Copyright 2016 by Medical Research Council

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


(define *cache-ligand-metrics* #f)


(define (to-python-arg-string ls)

  (define (item-to-string item)
    (cond
     ((string? item) (if (string=? item ",")
			 item
			 (string-append "'" item "'")))
     ((number? item) (number->string item))
     ((equal? item #t) "True")
     ((equal? item #f) "False")
     (else
      "False")))

  (let ((string-parts (map item-to-string ls)))
    (let ((r (string-concatenate string-parts)))
      (string-append "[" r "]"))))



;; return 3 values 
;; 
(define (mtz-file-name->refinement-column-labels file-name)

  (let ((columns (get-mtz-columns file-name)))
    (let ((read-success (mtz-column-types-info-t-read-success-get columns)))
      
      (if (not (= read-success 1))
	  
	  ;; unhappy path
	  (begin 
	    (format #t "Failed to read columms from file ~s for map molecule ~s~%" file-name imol-map)
	    (values #f #f #f))
	  
	  ;; happy path
	  (let ((f-cols (mtz-column-types-info-t-f-cols-get columns))
		(sigf-cols (mtz-column-types-info-t-sigf-cols-get columns))
		(rfree-cols (mtz-column-types-info-t-r-free-cols-get columns)))

	    (let ((l1 (vector-mtz-type-label-length     f-cols))
		  (l2 (vector-mtz-type-label-length  sigf-cols))
		  (l3 (vector-mtz-type-label-length rfree-cols)))

	      (if (not (all-true? (list
				   (> l1 0)
				   (> l2 0)
				   (> l3 0))))
		  
		  ;; unhappy path
		  (begin
		    (format #t "Failed to find columns of the necessary types from ~s : ~s ~s ~s~%~!"
			    file-name l1 l2 l3)
		    (values #f #f #f))

		  ;; happy path
		  (let ((f-col      (vector-mtz-type-label-ref     f-cols 0))
			(sigf-col   (vector-mtz-type-label-ref  sigf-cols 0))
			(r-free-col (vector-mtz-type-label-ref rfree-cols 0)))
		    
		    (let ((f-col-label      (mtz-type-label-column-label-get      f-col))
			  (sigf-col-label   (mtz-type-label-column-label-get   sigf-col))
			  (r-free-col-label (mtz-type-label-column-label-get r-free-col)))

		      (values f-col-label sigf-col-label  r-free-col-label)
		      
		      )))))))))



(define (ligand-validation-metrics-gui-list-wrapper imol chain-id res-no ins-code refmac-input-mtz-file-name
						    fp-col sigfp-col rfree-col refmac-dir)


  (let ((m (get-metrics-for-ligand imol chain-id res-no ins-code refmac-input-mtz-file-name 
				   fp-col sigfp-col rfree-col refmac-dir)))

    (if *cache-ligand-metrics*
	(set! m (list 0.935339450836182 10.6290054321289 (list 0 10 18 70 143) 
		      (list -0.0306930494488489 0.386785203839611 0.0146206432429434 10865.0 
			    5373.01106146351 -27.3189968762454 0.0173681300165985 0.0280039555722363 
			    -8.86567375069092e-10 -0.0025144037621947 0.117837198078632 0.120915851909265) 
		      (list 0.952432048573266 13.3150000572205 13.9800000190735 0.176323987538941))))

    (if (list? m)
	
	(let* ((diff-d 0.05)
	       (d            (list-ref m 0))
	       (mwz          (list-ref m 1))  ;; mogul Z
	       (contact-info (list-ref m 2))  ;; number of bad contacts
	       (bc (list-ref contact-info 0))
	       (low-is-good 1)
	       (high-is-good 0))

	  ;; if mogul ran OK, then we can display the mogul markup
	  ;; 
	  (if (pair? mwz) ;; mwz an improper pair (mogul-z-worst . mogul-out-file-name) or an error status
	      (let ((mogul-out-file-name (cdr mwz)))
		(mogul-markup imol chain-id res-no ins-code mogul-out-file-name)))
	      

	  (let ((percentile-d      (get-ligand-percentile "density_correlation" d high-is-good))
		(percentile-diff-d (get-ligand-percentile "coot_diff_map_correlation" diff-d low-is-good))
		(percentile-mwz    (get-ligand-percentile "mogul_z_worst" (car mwz) low-is-good))
		(percentile-bc     (get-ligand-percentile "bumps_1" bc low-is-good))
		(c ",")) ;; comma
	    
	    (let ((input-to-sliders 
		   (string-append
		    "["
		    (to-python-arg-string (list "Direct map density correl." c percentile-d c d))
		    ","
		    (to-python-arg-string (list " Diff map density correl." c percentile-diff-d c diff-d))
		    ","
		    (to-python-arg-string (list "            Mogul Z-worst" c percentile-mwz c mwz))
		    ","
		    (to-python-arg-string (list "             Bad contacts" c percentile-bc  c bc))
		    "]")))
	      
	      (run-python-command (string-append "ligand_validation_metrics_gui_list_wrapper(" 
						 input-to-sliders
						 ")"))))))))

  

(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Ligand")))

      (add-simple-coot-menu-menuitem
       menu "Ligand Metric Sliders"
       (lambda ()

	 (let ((imol-map (imol-refinement-map)))

	   (if (not (valid-map-molecule? imol-map))
	       (add-status-bar-text "No valid refinement map molecule")

	       (let ((refmac-input-mtz-file-name 
		      (let ((l (refmac-parameters imol-map)))
			(if (null? l)
			    (mtz-file-name imol-map)
			    (car l))))
		     (refmac-dir (get-directory "coot-refmac")))

		 (receive (f-col-label sigf-col-label r-free-col-label)
			  (mtz-file-name->refinement-column-labels refmac-input-mtz-file-name)
			  (format #t "    f-col-label: ~s~%"      f-col-label)
			  (format #t " sigf-col-label: ~s~%"   sigf-col-label)
			  (format #t "rfree-col-label: ~s~%" r-free-col-label)

			  (using-active-atom
			   (ligand-validation-metrics-gui-list-wrapper 
			    aa-imol aa-chain-id aa-res-no aa-ins-code 
			    refmac-input-mtz-file-name 
			    f-col-label sigf-col-label r-free-col-label refmac-dir))))))))))




