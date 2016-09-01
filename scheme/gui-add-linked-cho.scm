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


(define (add-pyranose-pseudo-ring-plane-restraints comp-id)

  (define (filter-out plane-name-sub-string plane-restraints)
    plane-restraints)

  (let ((restraints (monomer-restraints comp-id)))

    (if (not (list? restraints))
	(begin
	  (format #t "failed to get ~s restraints~%" comp-id))
	(begin
	  (let ((new-restraints 
		 (let f ((restraints restraints))
		   (cond
		    ((null? restraints) '())
		    ((string=? (car (car restraints)) "_chem_comp_plane_atom")
		     (let* ((plane-restraints (cdr (car restraints)))
			    (new-plane-restraints
			     (cons "_chem_comp_plane_atom"
				   (append 
				    (list
				     (list "pseudo-ring-1" (list " C1 " " C2 " " C4 " " C5 ") 0.01)
				     (list "pseudo-ring-2" (list " C2 " " C3 " " C5 " " O5 ") 0.01)
				     (list "pseudo-ring-3" (list " C3 " " C4 " " O5 " " C1 ") 0.01)
				     )
				    (filter-out "pseudo-ring-3" plane-restraints)))))
		       (cons new-plane-restraints (f (cdr restraints)))))
		    (else 
		     (cons (car restraints)
			   (f (cdr restraints))))))))
	    ;; (format #t "new-restraints: ~s ~%" new-restraints)
	    (set-monomer-restraints comp-id new-restraints))))))


(define (add-synthetic-carbohydrate-planes)
  (add-pyranose-pseudo-ring-plane-restraints "NAG")
  (add-pyranose-pseudo-ring-plane-restraints "BMA")
  (add-pyranose-pseudo-ring-plane-restraints "MAN")
  (add-pyranose-pseudo-ring-plane-restraints "GAL")
  (add-pyranose-pseudo-ring-plane-restraints "GLC"))



(define (multi-add-linked-residue imol res-spec residues-to-add)

  (set-go-to-atom-molecule imol)
  (set-matrix 8)

  (let ((current-refinement-rate (dragged-refinement-steps-per-frame)))
  
    (let loop ((residues-to-add residues-to-add)
	       (current-res-spec res-spec))
	       

      (set-go-to-atom-from-res-spec current-res-spec)
      (set-dragged-refinement-steps-per-frame 300)

      (if (not (list? current-res-spec))
	  (begin
	    (format #t "OOps not a proper res-spec ~s with residues-to-add: ~s~%"
		    current-res-spec residues-to-add)
	    #f)
	  (cond
	   ((null? residues-to-add) 'done)
	   (else 
	    (let* ((new-res-pair (car residues-to-add)))
	      (if (not (pair? new-res-pair))
		  (begin
		    (format #t "Oops - not a residue-link string pair when adding new-res-pair~%"))
		  
		  (let ((new-res  (car new-res-pair))
			(new-link (cdr new-res-pair)))

		    (let ((new-res-spec (add-linked-residue imol 
							    (res-spec->chain-id current-res-spec)
							    (res-spec->res-no   current-res-spec)
							    (res-spec->ins-code current-res-spec)
							    new-res
							    new-link)))

		      ;; residues-near-residue takes a 3-part spec and makes 3-part specs.
		      (let* ((ls (residues-near-residue imol current-res-spec 1.9)))
			(with-auto-accept
			 (refine-residues imol (cons current-res-spec ls))))
		      
		      (loop (cdr residues-to-add) (cdr new-res-spec)))))))))

    ;; restore refinement mode
    (set-dragged-refinement-steps-per-frame current-refinement-rate))))


(define oligomannose-tree '(("NAG-ASN" . "NAG")
			    (("BETA1-4" . "NAG")
			     (("BETA1-4" . "BMA")
			      (("ALPHA1-6" . "MAN")
			       (("ALPHA1-6" . "MAN")
				(("ALPHA1-2" . "MAN")))
			       (("ALPHA1-3" . "MAN")
				(("ALPHA1-2" . "MAN"))))
			      (("ALPHA1-3" . "MAN")
			       (("ALPHA1-2" . "MAN")
				(("ALPHA1-2" . "MAN"))))))))

;(define oligomannose-tree '(("NAG-ASN" . "NAG")
;			    (("BETA1-4" . "NAG")
;			     (("BETA1-4" . "BMA")
;			      (("ALPHA1-6" . "MAN"))
;			      (("ALPHA1-3" . "MAN")
;			       (("ALPHA1-2" . "MAN")
;				(("ALPHA1-2" . "MAN"))))))))

;(define oligomannose-tree '(("NAG-ASN" . "NAG")
;			    (("BETA1-4" . "NAG")
;			     (("BETA1-4" . "BMA")
;			      (("ALPHA1-6" . "MAN"))
;			      (("ALPHA1-3" . "MAN")
;			       (("ALPHA1-2" . "MAN")))))))

			       

;;; testing tree
;(define oligomannose-tree '(("NAG-ASN" . "NAG")
;			    (("BETA1-4" . "NAG")
;			     (("BETA1-4" . "BMA")
;			      (("ALPHA1-6" . "MAN"))))))

;; (set-add-linked-residue-do-fit-and-refine 0)

;; now users can set this
(define *add-linked-residue-tree-correlation-cut-off* 0.6)

(define (add-linked-residue-tree imol parent tree)

  (define (centre-view-on-residue-centre res-spec)
    (let ((res-centre (residue-centre imol 
				      (residue-spec->chain-id res-spec)
				      (residue-spec->res-no res-spec)
				      (residue-spec->ins-code res-spec))))
      (if (list? res-centre)
	  (apply set-rotation-centre res-centre))))
	  
  
  
  (define func-test
    (let ((count 10))
      (lambda (parent res-pair)
	(let ((rv (list "B" count "")))
	  (set! count (+ count 1))
	  (format #t "process res-pair ~s with parent ~s producing ~s ~%" res-pair parent rv)
	  rv))))

  (define (well-fitting? res-spec)
    (using-active-atom
     (let ((neighbs (residues-near-residue aa-imol res-spec 4)))
       (let ((c (map-to-model-correlation imol (list res-spec) neighbs 0 (imol-refinement-map))))
	 (format #t "############# new residue ~s density correlation: ~s~%" res-spec c)
	 (if (not (number? c))
	     #f
	     (> c *add-linked-residue-tree-correlation-cut-off*))))))

  (define (delete-residue-by-spec spec)
    (delete-residue imol
		    (residue-spec->chain-id spec)
		    (residue-spec->res-no   spec)
		    (residue-spec->ins-code spec)))
		    

  (define (func parent res-pair)
    (if (not (list? parent))
        (begin
          (format #t "WARNING:: Oops not a proper res-spec ~s with residues-to-add: ~s~%"
                  parent res-pair)
	  #f)
	(if (not (pair? res-pair))
	    (begin
	      (format #t "Oops - not a residue-link string pair when adding res-pair~%" res-pair)
	      #f)
	    ;; OK! go!
	    (let ((new-link     (car res-pair))
		  (new-res-type (cdr res-pair)))

	      ;; (set-go-to-atom-from-res-spec parent)
	      (centre-view-on-residue-centre parent)

;	      (format #t "================= calling add-linked-residue with args ~s ~s ~s ~s ~s ~s~%"
;		      imol 
;		      (res-spec->chain-id parent)
;		      (res-spec->res-no   parent)
;		      (res-spec->ins-code parent)
;		      new-res-type
;		      new-link)

	      (let ((new-res-spec (add-linked-residue imol 
						      (res-spec->chain-id parent)
						      (res-spec->res-no   parent)
						      (res-spec->ins-code parent)
						      new-res-type
						      new-link)))
		(let* ((ls (residues-near-residue imol parent 1.9))
		       (local-ls (cons parent ls)))
		  (with-auto-accept (refine-residues imol local-ls))
		  (if (list? new-res-spec)
		      (begin
			(let ((preped-new-res-spec (cdr new-res-spec))) ;; strip off leading result
			  (if (well-fitting? preped-new-res-spec)
			      (begin
				preped-new-res-spec)
			      (begin
				(format #t "------------ That was not well-fitting. Deleting ~s: ~%"
					preped-new-res-spec)
				(delete-residue-by-spec preped-new-res-spec)
				(with-auto-accept (refine-residues imol local-ls))
				#f))))
		      #f))))))) ;; oops, something bad...

  (define (process-tree parent tree proc-func)
    (cond 
     ((null? tree) '())
     ((list? (car tree))
      (let ((part-1 (process-tree parent (car tree) proc-func))
	    (part-2 (process-tree parent (cdr tree) proc-func)))
	(cons part-1 part-2)))
     (else 
      (let ((new-res (proc-func parent (car tree))))
	(cons new-res
	      (process-tree new-res (cdr tree) proc-func))))))

  ;; main line of add-linked-residue
  ;; 
  (add-synthetic-carbohydrate-planes)
  (set-residue-selection-flash-frames-number 1)
  (set-go-to-atom-molecule imol)
  (let* ((previous-m (default-new-atoms-b-factor))
	(m (median-temperature-factor imol))
	(new-m (* m 1.55)))

    (if (number? m)
	(set-default-temperature-factor-for-new-atoms new-m))

    (process-tree parent tree func)
    (set-default-temperature-factor-for-new-atoms previous-m)))

  

(define (add-module-carbohydrate) 

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "Glyco")))

	(add-simple-coot-menu-menuitem
	 menu "Add a ASN-NAG NAG"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "NAG" "NAG-ASN"))))

	(add-simple-coot-menu-menuitem
	 menu "Add a BETA1-4 NAG"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "NAG" "BETA1-4"))))

	(add-simple-coot-menu-menuitem
	 menu "Add a BETA1-4 BMA"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "BMA" "BETA1-4"))))

	(add-simple-coot-menu-menuitem
	 menu "Add an ALPHA1-2 MAN"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "MAN" "ALPHA1-2"))))

	(add-simple-coot-menu-menuitem
	 menu "Add an ALPHA1-3 MAN"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "MAN" "ALPHA1-3"))))

	(add-simple-coot-menu-menuitem
	 menu "Add an ALPHA2-3 MAN"
	 (lambda () 
	   (set-matrix 8)
	   ;; we should do this only if we are sitting on an SIA.
	   ;; Attaching a SIA to a MAN (i.e. reverse order) would be a
	   ;; good test too...
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "MAN" "ALPHA2-3"))))

	(add-simple-coot-menu-menuitem
	 menu "Add an ALPHA2-3 GAL"
	 (lambda () 
	   (set-matrix 8)
	   ;; we should do this only if we are sitting on an SIA.
	   ;; Attaching a SIA to a MAN (i.e. reverse order) would be a
	   ;; good test too...
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "GAL" "ALPHA2-3"))))

	(add-simple-coot-menu-menuitem
	 menu "Add an ALPHA1-6 MAN"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "MAN" "ALPHA1-6"))))

	(add-simple-coot-menu-menuitem
	 menu "Add an ALPHA1-6 FUC"
	 (lambda () 
	   (set-matrix 8)
	   (using-active-atom
	    (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code "FUC" "ALPHA1-6"))))

	(add-simple-coot-menu-menuitem
	 menu "N-link add NAG, NAG, BMA"
	 (lambda ()
	   (using-active-atom
	    (multi-add-linked-residue 
	     aa-imol
	     (list aa-chain-id aa-res-no aa-ins-code)
	     (list 
	      (cons "NAG" "NAG-ASN")
	      (cons "NAG" "BETA1-4")
	      (cons "BMA" "BETA1-4")
	      )))))

	(add-simple-coot-menu-menuitem
	 menu "Auto Fit & Refine on Link Addition"
	 (lambda ()
	   (set-add-linked-residue-do-fit-and-refine 1)))

	(add-simple-coot-menu-menuitem
	 menu "Auto Fit & Refinement Off for Link Addition"
	 (lambda ()
	   (set-add-linked-residue-do-fit-and-refine 0)))

	(add-simple-coot-menu-menuitem
	 menu "Add Oligomannose"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    ;; (with-no-backups aa-imol
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     oligomannose-tree))))

	(add-simple-coot-menu-menuitem
	 menu "Torsion Fit this residue"
	 (lambda ()
	   (using-active-atom
	    (multi-residue-torsion-fit aa-imol (list (list aa-chain-id aa-res-no aa-ins-code)) 30000))))

	(add-simple-coot-menu-menuitem 
	 menu "Torsion Fit & Refine this residue"
	 (lambda ()
	   (using-active-atom
	    (let ((centre-residue (list aa-chain-id aa-res-no aa-ins-code)))
	      (multi-residue-torsion-fit aa-imol (list centre-residue ) 30000)
	      (with-auto-accept
	       (refine-residues aa-imol (list centre-residue)))))))

	(add-simple-coot-menu-menuitem
	 menu "Add synthetic pyranose plane restraints"
	 (lambda () 
	   (add-synthetic-carbohydrate-planes)))

	)))

