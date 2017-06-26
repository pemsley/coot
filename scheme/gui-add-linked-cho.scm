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
    (format #t "DEBUG:: in filter-out: remove ~s from ~s~%"
	    plane-name-sub-string plane-restraints)
    ;; (filter (lambda (s) (string-match? plane-name-sub-string ...?)
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
				    (filter-out "pseudo-ring-" plane-restraints)))))
		       (cons new-plane-restraints (f (cdr restraints)))))
		    (else 
		     (cons (car restraints)
			   (f (cdr restraints))))))))
	    ;; (format #t "new-restraints: ~s ~%" new-restraints)
	    (set-monomer-restraints comp-id new-restraints))))))


(define (add-synthetic-pyranose-planes)
  (for-each (lambda(comp-id)
	      (add-pyranose-pseudo-ring-plane-restraints comp-id))
	    (list "NAG" "BMA" "MAN" "GAL" "GLC" "FUC" "XYP")))

(define (use-unimodal-pyranose-ring-torsions)
  (for-each (lambda(tlc)
	      (use-unimodal-ring-torsion-restraints tlc))
	    (list "NAG" "BMA" "MAN" "GAL" "GLC" "FUC" "XYP")))

;; fill this later
(define (add-cho-restraints-for-residue imol new-res-spec)
  #f)

(define (multi-add-linked-residue imol res-spec residues-to-add)

  (format #t "---------------- multi-add-linked-residue ~s ~s ~%~!" imol res-spec)
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
							    new-res new-link 2))) ;; add and fit mode

		      ;; residues-near-residue takes a 3-part spec and makes 3-part specs.
		      (let* ((ls (residues-near-residue imol current-res-spec 1.9)))

			(add-cho-restraints-for-residue imol new-res-spec)

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

(define paucimannose-tree '(("NAG-ASN" . "NAG")
			    (("ALPHA1-3" . "FUC"))
			    (("BETA1-4" . "NAG")
			     (("BETA1-4" . "BMA")
			      (("ALPHA1-6" . "MAN"))
			      (("ALPHA1-3" . "MAN"))
			      (("XYP-BMA"  . "XYP"))))))

(define complex-tree '(("NAG-ASN" . "NAG")
		       (("BETA1-6" . "FUL"))
		       (("BETA1-4" . "NAG")
			(("BETA1-4" . "BMA")
			 (("ALPHA1-6" . "MAN")
			  (("BETA1-2" . "NAG")))))))
;			   (("BETA1-4" . "GAL"))))
;			 (("ALPHA1-3" . "MAN")
;			  (("BETA1-2"  . "NAG")
;			   (("BETA1-4" . "GAL"))))))))


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

	      (let* ((tree-residues (glyco-tree-residues imol parent))
		     (imol-save (new-molecule-by-residue-specs imol tree-residues))
		     (new-res-spec (add-linked-residue imol
						       (res-spec->chain-id parent)
						       (res-spec->res-no   parent)
						       (res-spec->ins-code parent)
						       new-res-type
						       new-link 2))) ;; add and link mode

		(set-mol-displayed imol-save 0)
		(set-mol-active    imol-save 0)
		(let* ((ls (residues-near-residue imol parent 1.9))
		       (local-ls (cons parent ls)))
		  (add-cho-restraints-for-residue imol new-res-spec)
		  (rotate-y-scene 100 0.5)
		  (with-auto-accept (refine-residues imol local-ls))
		  (if (list? new-res-spec)
		      (begin
			(let ((preped-new-res-spec (cdr new-res-spec))) ;; strip off leading result
			  (if (well-fitting? preped-new-res-spec)
			      (begin
				preped-new-res-spec)
			      (begin
				;; ------------ bad fit -----------------
				;; delete residue and restore others
				(format #t "------------ That was not well-fitting. Deleting ~s: ~%"
					preped-new-res-spec)
				(delete-extra-restraints-for-residue-spec imol preped-new-res-spec)
				(delete-residue-by-spec preped-new-res-spec)
				;; restore glyco-tree residues from imol-save
				(replace-fragment imol imol-save "//")
				;; (with-auto-accept (refine-residues imol local-ls))
				#f))))
		      #f))))))
    ) ;; oops, something bad...

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

  ;; main line of add-linked-residue-tree
  ;;
  (add-synthetic-pyranose-planes)
  (use-unimodal-pyranose-ring-torsions)
  (set-refine-with-torsion-restraints 1)
  (set-matrix 8)
  (set-residue-selection-flash-frames-number 1)
  (set-go-to-atom-molecule imol)
  (let* ((previous-m (default-new-atoms-b-factor))
	 (m (median-temperature-factor imol))
	 (new-m (* m 1.55)))

    (set-default-temperature-factor-for-new-atoms previous-m)

    (if (number? m)
	(set-default-temperature-factor-for-new-atoms new-m))

    (let ((start-pos-view (add-view-here "Glyo Tree Start Pos")))
      (process-tree parent tree func)
      (go-to-view-number start-pos-view 0)
      ;; add a test here that the tree here (centre of screen) matches a known tree.
      ;; 
      ;; and that each is 4C1 (or 1C4 for FUC?) (XYP?)
      )))


(define (add-linked-residue-with-extra-restraints-to-active-residue new-res-type link-type)
  (set-matrix 8)
  (set-refine-with-torsion-restraints 1)
  (set-add-linked-residue-do-fit-and-refine 0)
  (using-active-atom
   (let ((new-res-spec (add-linked-residue aa-imol aa-chain-id aa-res-no aa-ins-code new-res-type link-type 2)))
     (if (list? new-res-spec)
	 (begin
	   (add-cho-restraints-for-residue aa-imol new-res-spec)
	   ;; refine that
	   (with-auto-accept
	    (let ((residues (cons aa-res-spec (residues-near-residue aa-imol aa-res-spec 1.9))))
	      (refine-residues aa-imol residues))))))))

(define (delete-all-cho)
  (let ((delete-cho-list '()))
    (using-active-atom
     (if (valid-model-molecule? aa-imol)
	 (begin
	   (for-each (lambda (chain-id)
		       (for-each (lambda (res-serial)
				   (let ((res-no (seqnum-from-serial-number aa-imol chain-id res-serial)))
				     (let ((rn (residue-name aa-imol chain-id res-no "")))
				       (if (string? rn)
					   (if (or (string=? "NAG" rn) (string=? "MAN" rn) (string=? "BMA" rn) (string=? "FUL" rn)
						   (string=? "FUC" rn) (string=? "XYP" rn) (string=? "SIA" rn) (string=? "GAL" rn))
					       (let* ((residue-spec (list chain-id res-no "")))
						 (set! delete-cho-list (cons (list chain-id res-no "") delete-cho-list))))))))
				 (range (chain-n-residues chain-id aa-imol))))
		     (chain-ids aa-imol))

	   ;; now we have delete-residues, we don't need to delete them one by one
	   ;;(for-each (lambda(cho-res-spec)
	   ;; (delete-residue aa-imol (residue-spec->chain-id cho-res-spec) (residue-spec->res-no cho-res-spec) ""))
	   ;;   delete-cho-list)))))))
	   ;;
	   (delete-residues aa-imol delete-cho-list))))))

(define (interactive-add-cho-dialog)

  ;; (define (update-for-current-residue-inner vbox))
  ;; (define (update-for-current-residue) (update-for-current-residue-inner vbox))
  (add-synthetic-pyranose-planes)
  (use-unimodal-pyranose-ring-torsions)
  (let ((buttons (list ;; (list label func)
		  (list "Update for Current Residue" (lambda () (format #t "dummy\n")))

		  (list "Add a NAG-ASN NAG"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "NAG" "NAG-ASN")))
		  (list "Add a BETA1-4 NAG"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "NAG" "BETA1-4")))
		  (list "Add a BETA1-4 BMA"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "BMA" "BETA1-4")))
		  (list "Add an ALPHA1-2 MAN"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "MAN" "ALPHA1-2")))
		  (list "Add an ALPHA1-3 MAN"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "MAN" "ALPHA1-3")))
		  (list "Add an ALPHA2-3 MAN"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "MAN" "ALPHA2-3")))
		  (list "Add an ALPHA2-3 GAL"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "GAL" "ALPHA2-3")))
		  (list "Add an ALPHA1-6 MAN"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "MAN" "ALPHA1-6")))
		  (list "Add a BETA1-2 NAG"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "NAG" "BETA1-2")))
		  (list "Add a BETA1-4 GAL"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "GAL" "BETA1-4")))
		  (list "Add an ALPHA1-3 FUC"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "FUC" "ALPHA1-3" )))
		  (list "Add an ALPHA1-6 FUC"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "FUC" "ALPHA1-6")))
		  (list "Add an BETA1-6 FUL"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "FUL" "BETA1-6")))
		  (list "Add an XYP-BMA XYP"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "XYP" "XYP-BMA"))))))
    (let ((vbox (dialog-box-of-buttons "Add N-linked Glycan" (cons 360 520) buttons "Close")))
      (gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox)
      ;; set the callback on the first button
      (let ((children  (gtk-container-children vbox)))
	(if (list? children)
	    (if (> (length children) 0)
		(let ((first-button (car children)))
		  (gtk-signal-connect first-button
				      "clicked"
				      (lambda ()
					(gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox)))))))
      ;; add a widget to allow the user to choose the tree type
      (let* ((hbox-1 (gtk-hbox-new #f 2))
	     (hbox-2 (gtk-hbox-new #f 2))
	     (butt-1 (gtk-radio-button-new-with-label #f "Oligomannose"))
	     (butt-2 (gtk-radio-button-new-with-label butt-1 "Hybrid"))
	     (butt-3 (gtk-radio-button-new-with-label butt-1 "Complex"))
	     (butt-4 (gtk-radio-button-new-with-label butt-1 "Expert Mode")))

	(gtk-box-pack-start hbox-1 butt-1 #f #f 2)
	(gtk-box-pack-start hbox-1 butt-2 #f #f 2)
	(gtk-box-pack-start hbox-2 butt-3 #f #f 2)
	(gtk-box-pack-start hbox-2 butt-4 #f #f 2)

	(gtk-widget-show butt-1)
	(gtk-widget-show butt-2)
	(gtk-widget-show butt-3)
	(gtk-widget-show butt-4)
	(gtk-widget-show hbox-1)
	(gtk-widget-show hbox-2)
	(gtk-box-pack-start vbox hbox-1 #f #f 2)
	(gtk-box-pack-start vbox hbox-2 #f #f 2)
	(gtk-box-set-homogeneous hbox-1 #t)
	(gtk-box-set-homogeneous hbox-2 #t)
	(gtk-box-reorder-child vbox hbox-1 0)
	(gtk-box-reorder-child vbox hbox-2 1))

      ;; global var post-set-rotation-centre-hook
      (set! post-set-rotation-centre-hook
	    (lambda ()
	      (gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox))))))

;;
(define (glyco-tree-dialog-set-button-active-state button glyco-id tree-type)

  (define (get-sensitive-button-list glyco-id tree-type)

    (if (not (list? glyco-id))
	'()
	(let ((level-number        (list-ref glyco-id 0))
	      (residue-type        (list-ref glyco-id 1))
	      (link-type           (list-ref glyco-id 2))
	      (parent-residue-type (list-ref glyco-id 3))
	      (residue-spec        (list-ref glyco-id 4)))
	    (let ((active-button-label-list '()))

	      (if (eq? tree-type 'expert-mode)
		  (set! active-button-label-list 'expert-mode))

	      (if (eq? tree-type 'oligomannose)

		  (begin
		    (if (= level-number 0)
			(if (string=? residue-type "ASN")
			    (set! active-button-label-list (list "Add a NAG-ASN NAG"))))

		    (if (= level-number 1)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 NAG"))))
		    (if (= level-number 2)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 BMA"))))
		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"))))
		    (if (= level-number 4)
			(if (string=? residue-type "MAN")
			    (begin
			      (if (string=? link-type "ALPHA1-3")
				  (set! active-button-label-list (list "Add an ALPHA1-2 MAN")))
			      (if (string=? link-type "ALPHA1-6")
				  (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								       "Add an ALPHA1-6 MAN"))))))
		    (if (= level-number 5)
			(if (string=? residue-type "MAN")
			    (begin
			      (if (string=? link-type "ALPHA1-2")
				  (set! active-button-label-list (list "Add an ALPHA1-2 MAN")))
			      (if (string=? link-type "ALPHA1-6")
				  (set! active-button-label-list (list "Add an ALPHA1-2 MAN")))
			      (if (string=? link-type "ALPHA1-3")
				  (set! active-button-label-list (list "Add an ALPHA1-2 MAN"))))))

		    (if (= level-number 6)
			(if (string=? residue-type "MAN")
			    (begin
			      (if (string=? link-type "ALPHA1-2")
				  (set! active-button-label-list (list "Add an ALPHA1-3 GLC"))))))

		    (if (= level-number 7)
			(if (string=? residue-type "GLC")
			    (begin
			      (if (string=? link-type "ALPHA1-2")
				  (set! active-button-label-list (list "Add an ALPHA1-3 GLC"))))))

		    (if (= level-number 8)
			(if (string=? residue-type "GLC")
			    (begin
			      (if (string=? link-type "ALPHA1-2")
				  (set! active-button-label-list (list "Add an ALPHA1-2 GLC"))))))))

	      ;; hybrid
	      (if (eq? tree-type 'paucimannose)

		  (begin

		    (if (= level-number 0)
			(if (string=? residue-type "ASN")
			    (set! active-button-label-list (list "Add a NAG-ASN NAG"))))

		    (if (= level-number 1)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 NAG"
								 "Add an ALPHA1-3 FUC"))))

		    (if (= level-number 2)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 BMA"))))

		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"))))
		    (if (= level-number 4)
			(if (string=? residue-type "MAN")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add an XYP-BMA XYP"))))
		    ))

	      ;; complex
	      ;;
	      (if (eq? tree-type 'complex)

		  (begin
		    (if (= level-number 0)
			(if (string=? residue-type "ASN")
			    (set! active-button-label-list (list "Add a NAG-ASN NAG"))))

		    (if (= level-number 1)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 NAG"))))

		    (if (= level-number 2)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 BMA"))))

		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add an XYP-BMA XYP"))))))


	      active-button-label-list))))

  ;; main line
  ;;
  ;; (format #t "here in glyco-tree-dialog-set-button-active-state ~s ~s ~s ~%" button glyco-id tree-type)
  (let ((l (gtk-button-get-label button)))
    (let ((active-button-label-list (get-sensitive-button-list glyco-id tree-type)))
      ;; (format #t "active-button-label-list: ~s~%" active-button-label-list)
      (if (eq? active-button-label-list 'expert-mode)
	  (gtk-widget-set-sensitive button #t)
	  (if (not (string=? l "Update for Current Residue"))
	      (gtk-widget-set-sensitive button (string-member? l active-button-label-list)))))))

;; return a value
(define (gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox)

  (define (get-tree-type)
    (let ((tree-type 'oligomannose))
      (let ((children (gtk-container-children vbox)))
	(for-each (lambda (child)
		    ;; (format #t "child: ~s~%" child)
		    (if (gtk-box? child)
			(begin
			  (for-each (lambda (box-child)
				      (if (gtk-radio-button? box-child)
					  (if (gtk-toggle-button-get-active box-child)
					      (let ((l (gtk-button-get-label box-child)))
						(if (string=? l "Oligomannose") (set! tree-type 'oligomannose))
						(if (string=? l "Hybrid")       (set! tree-type 'hybrid))
						(if (string=? l "Expert Mode")  (set! tree-type 'expert-mode))
						(if (string=? l "Complex")      (set! tree-type 'complex))))))
				    (gtk-container-children child)))))
		  children)
	tree-type)))

  (using-active-atom
   (let ((glyco-id (glyco-tree-residue-id aa-imol aa-res-spec)))
     ;; (format #t "glyco-id (first): ~s~%" glyco-id)
     ;; if it was an ASP create a level-0 glyco-id for that (glyco-tree-residue-id doesn't
     ;; do that (not sure why)).
     (if (not glyco-id)
	 (let ((rn (residue-name aa-imol aa-chain-id aa-res-no aa-ins-code)))
	   (if (string? rn)
	       (if (string=? rn "ASN")
		   (set! glyco-id (list 0 "ASN" "" "" aa-res-spec))))))
     ;; (format #t "glyco-id (second): ~s~%" glyco-id)
     (if (list? glyco-id)
	 (let ((tree-type (get-tree-type)))
	   (let ((children (gtk-container-children vbox)))
	     (for-each (lambda (child)
			 (if (gtk-button? child)
			     (glyco-tree-dialog-set-button-active-state child glyco-id tree-type)))
		       children)
	     #t)
	   #f)))))

(define (add-module-carbohydrate) 

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "Glyco")))

	(add-simple-coot-menu-menuitem
	 menu "N-linked Glycan Addition Dialog..."
	 (lambda ()
	   (interactive-add-cho-dialog)))

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

;  the mode in the function call now takes take of this
;
;	(add-simple-coot-menu-menuitem
;	 menu "Auto Fit & Refine on Link Addition"
;	 (lambda ()
;	   (set-add-linked-residue-do-fit-and-refine 1)))

;	(add-simple-coot-menu-menuitem
;	 menu "Auto Fit & Refinement Off for Link Addition"
;	 (lambda ()
;	   (set-add-linked-residue-do-fit-and-refine 0)))

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
	 menu "Add Paucimannose"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    ;; (with-no-backups aa-imol
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     paucimannose-tree))))

	(add-simple-coot-menu-menuitem
	 menu "Add Complex Tree"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     complex-tree))))

	(add-simple-coot-menu-menuitem
	 menu "Delete All Carbohydrate"
	 (lambda()
	   (delete-all-cho)))

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
	   (add-synthetic-pyranose-planes)))

	(add-simple-coot-menu-menuitem
	 menu "Use Unimodal ring torsion restraints"
	 (lambda ()
	   (use-unimodal-pyranose-ring-torsions)))

	)))
