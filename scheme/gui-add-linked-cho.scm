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



(define (interactive-add-cho-dialog)

  ;; (define (update-for-current-residue-inner vbox))
  ;; (define (update-for-current-residue) (update-for-current-residue-inner vbox))
  (add-synthetic-pyranose-planes)
  (use-unimodal-pyranose-ring-torsions)
  (let ((buttons (list ;; (list label func)
		  (list "Update for Current Residue" (lambda () (format #t "dummy\n")))

		  (list "Refine Tree" (lambda ()
					(with-auto-accept
					 (using-active-atom
					  (refine-residues aa-imol
							   (glyco-tree-residues aa-imol aa-res-spec))))))

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
		  (list "Add an ALPHA1-2 FUC"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "FUC" "ALPHA1-2" )))
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
			  (add-linked-residue-with-extra-restraints-to-active-residue "XYP" "XYP-BMA")))
		  (list "Add an ALPHA2-3 SIA"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "SIA" "ALPHA2-3"))) ;; do these work?
		  (list "Add an ALPHA2-6 SIA"
			(lambda ()
			  (add-linked-residue-with-extra-restraints-to-active-residue "SIA" "ALPHA2-6")))
		  )))

    (let ((vbox (car (dialog-box-of-buttons "Add N-linked Glycan" (cons 420 600) buttons "Close"))))
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
      (let* ((hbox-1 (gtk-hbox-new #t 2))
	     (hbox-2 (gtk-hbox-new #t 2))
	     (table (gtk-table-new 3 2 #f))
	     (butt-1 (gtk-radio-button-new-with-label #f     "High Mannose"))
	     (butt-2 (gtk-radio-button-new-with-label butt-1 "Hybrid (Mammal)"))
	     (butt-4 (gtk-radio-button-new-with-label butt-1 "Complex (Mammal)"))
	     (butt-5 (gtk-radio-button-new-with-label butt-1 "Complex (Plant) "))
	     (butt-6 (gtk-radio-button-new-with-label butt-1 "Expert User Mode")))

	(gtk-widget-show butt-1)
	(gtk-widget-show butt-2)
	(gtk-widget-show butt-4)
	(gtk-widget-show butt-5)
	(gtk-widget-show butt-6)

; table

	;; add buttons for nice(?) layout/order
	(gtk-table-attach table butt-1 0 1 0 1 '(expand fill) '(expand fill) 0 0) ;; high mannose
	(gtk-table-attach table butt-4 1 2 0 1 '(expand fill) '(expand fill) 0 0) ;; complex mammal
	(gtk-table-attach table butt-6 0 1 1 2 '(expand fill) '(expand fill) 0 0) ;; Expert
	(gtk-table-attach table butt-5 1 2 1 2 '(expand fill) '(expand fill) 0 0) ;; complex plant
	(gtk-table-attach table butt-2 2 3 0 1 '(expand fill) '(expand fill) 0 0) ;; hybrid mammal

	(gtk-box-pack-start vbox table #t #t 2)
	(gtk-widget-show table)
	(gtk-box-reorder-child vbox table 0)

	(for-each (lambda (butt)
		    (gtk-signal-connect butt "toggled"
					(lambda ()
					  (gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox))))
		  (list butt-1 butt-2 butt-4 butt-5 butt-6)))

      ;; global var post-set-rotation-centre-hook
      (set! post-set-rotation-centre-hook
	    (lambda ()
	      (gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox))))))

;;
(define (glyco-tree-dialog-set-button-active-state button glyco-id tree-type)

  (define (glyco-id->level-number glyco-id)
    (list-ref glyco-id 0))

  (define (glyco-id->prime-arm-sym glyco-id)
    (list-ref glyco-id 1))

  (define (glyco-id->residue-type glyco-id)
    (list-ref glyco-id 2))

  (define (glyco-id->link-type glyco-id)
    (list-ref glyco-id 3))

  (define (glyco-id->parent-residue-type glyco-id)
    (list-ref glyco-id 4))

  (define (glyco-id->residue-spec glyco-id)
    (list-ref glyco-id 5))

  ;; can return a list (of label strings) or 'expert-user-mode (because there is no test)
  ;; in that case
  (define (get-sensitive-button-list glyco-id tree-type)

    (if (not (list? glyco-id))
	'()
	(let ((level-number        (glyco-id->level-number        glyco-id))
	      (prime-arm-sym       (glyco-id->prime-arm-sym       glyco-id))
	      (residue-type        (glyco-id->residue-type        glyco-id))
	      (link-type           (glyco-id->link-type           glyco-id))
	      (parent-residue-type (glyco-id->parent-residue-type glyco-id))
	      (residue-spec        (glyco-id->residue-spec        glyco-id)))
	    (let ((active-button-label-list '()))

	      (if (eq? tree-type 'expert-user-mode)
		  (set! active-button-label-list 'expert-user-mode))

	      ;; ------------------------------------------------------------------------------------
	      ;; Note the trees tested here match those from (get-tree-type) which examines the button
	      ;; label of the active radio button in the dialog (these are not the auto-build trees)
	      ;; ------------------------------------------------------------------------------------
	      ;;
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

	      ;; hybrid mammal
	      (if (eq? tree-type 'hybrid-mammal)

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

		    (if (= level-number 2)
			(if (string=? residue-type "FUC")
			    (set! active-button-label-list (list "Add a BETA1-4 GAL"))))

		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"))))

		    (if (= level-number 3)
			(if (string=? residue-type "GAL")
			    (set! active-button-label-list (list "Add an ALPHA1-2 FUC"))))

		    (if (= level-number 4)
			(if (string=? residue-type "MAN")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"))))

		    (if (= level-number 5)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 GAL"))))

		    (if (= level-number 5)
			(if (string=? residue-type "GAL")
			    (set! active-button-label-list (list "Add an ALPHA2-3 SIA"
								 "Add an ALPHA2-6 SIA"
								 ))))
		    ))

	      ;; hybrid plant
	      (if (eq? tree-type 'hybrid-plant)

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

		    (if (= level-number 2)
			(if (string=? residue-type "FUC")
			    (set! active-button-label-list (list "Add a BETA1-4 GAL"))))

		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add an XYP-BMA XYP"))))

		    (if (= level-number 3)
			(if (string=? residue-type "GAL")
			    (set! active-button-label-list (list "Add an ALPHA1-2 FUC"))))

		    ;; note that level-number is not enough for complete disambiguation,
		    ;; we need to know we are 4 or 4' (Vliegenthart et al 1983 nomenclature).
		    ;;
		    (if (= level-number 4)
			(if (string=? residue-type "MAN")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add a BETA1-2 NAG"))))

		    (if (= level-number 5)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 GAL"))))

		    (if (= level-number 5)
			(if (string=? residue-type "GAL")
			    (set! active-button-label-list (list "Add an ALPHA2-3 SIA"
								 "Add an ALPHA2-6 SIA"
								 ))))
		    ))

	      ;; complex-mammal
	      ;;
	      (if (eq? tree-type 'complex-mammal)

		  (begin
		    (if (= level-number 0)
			(if (string=? residue-type "ASN")
			    (set! active-button-label-list (list "Add a NAG-ASN NAG"))))

		    (if (= level-number 1)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 NAG"
								 "Add an ALPHA1-6 FUC"))))
		    (if (= level-number 2)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 BMA"))))

		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add a BETA1-4 NAG"))))

		    ;; note that level-number is not enough for complete disambiguation, 
		    ;; we need to know we are 4 or 4' (Vl...? nomenclature).
		    ;;
		    (if (= level-number 4)
			(if (string=? residue-type "MAN")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add a BETA1-2 NAG"))))

		    (if (= level-number 5)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 GAL"))))

		    (if (= level-number 6)
			(if (string=? residue-type "GAL")
			    (set! active-button-label-list (list "Add an ALPHA2-3 SIA"
								 "Add an ALPHA2-6 SIA"
								 ))))
		    ))


	      ;; complex-plant
	      ;;
	      (if (eq? tree-type 'complex-plant)

		  (begin
		    (if (= level-number 0)
			(if (string=? residue-type "ASN")
			    (set! active-button-label-list (list "Add a NAG-ASN NAG"))))

		    (if (= level-number 1)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 NAG"
								 "Add an ALPHA1-3 FUC"
								 "Add an ALPHA1-6 FUC"))))
		    (if (= level-number 2)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 BMA"))))

		    (if (= level-number 3)
			(if (string=? residue-type "BMA")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add an XYP-BMA XYP"
								 "Add an BETA1-4 NAG"))))

		    (if (= level-number 4)
			(if (string=? residue-type "MAN")
			    (set! active-button-label-list (list "Add an ALPHA1-3 MAN"
								 "Add an ALPHA1-6 MAN"
								 "Add a BETA1-2 NAG"))))

		    (if (= level-number 5)
			(if (string=? residue-type "NAG")
			    (set! active-button-label-list (list "Add a BETA1-4 GAL"))))

		    (if (= level-number 6)
			(if (string=? residue-type "GAL")
			    (set! active-button-label-list (list "Add an ALPHA2-3 SIA"
								 "Add an ALPHA2-6 SIA"))))
		    ))
	      active-button-label-list))))

  ;; main line
  ;;
  ;; (format #t "here in glyco-tree-dialog-set-button-active-state ~s ~s ~s ~%"
  ;;        button glyco-id tree-type)
  ;;
  (let ((l (gtk-button-get-label button)))
    (let ((active-button-label-list (get-sensitive-button-list glyco-id tree-type)))
      ;; (format #t "active-button-label-list: ~s~%" active-button-label-list)
      (if (eq? active-button-label-list 'expert-user-mode)
	  (gtk-widget-set-sensitive button #t)
	  (if (and (not (string=? l "Update for Current Residue"))
                   (not (string=? l "Refine Tree")))
	      ;; don't do this for active-button-label-list is 'expert-user-mode
	      (gtk-widget-set-sensitive button (string-member? l active-button-label-list)))))))


;; vbox is the vbox of the dialog box of buttons. One of the children of the vbox
;; is the table that contains the buttons
;;
(define (gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook vbox)

  (define (get-tree-type)
    (let ((tree-type 'oligomannose))
      (let ((children (gtk-container-children vbox)))
	(for-each (lambda (child)
		    (if (gtk-table? child)
			(begin
			  (for-each (lambda (table-child)
				      (if (gtk-radio-button? table-child)
					  (if (gtk-toggle-button-get-active table-child)
					      (let ((l (gtk-button-get-label table-child)))
						;; this is a bit ugly because we are testing that these strings
						;; match button labels (set in interactive-add-cho-dialog)
						(if (string=? l "High Mannose") (set! tree-type 'oligomannose))
						(if (string=? l "Hybrid (Mammal)")  (set! tree-type 'hybrid-mammal))
						(if (string=? l "Hybrid (Plant)")   (set! tree-type 'hybrid-plant))
						(if (string=? l "Expert User Mode") (set! tree-type 'expert-user-mode))
						(if (string=? l "Complex (Mammal)") (set! tree-type 'complex-mammal))
						(if (string=? l "Complex (Plant)")  (set! tree-type 'complex-plant))
						))))
				    (gtk-container-children child)))))
		  children)
	tree-type)))

  (using-active-atom
   (let ((glyco-id (glyco-tree-residue-id aa-imol aa-res-spec)))
     
     ;; debug
     ;; (format #t "gui-add-linked-cho-dialog-vbox-set-rotation-centre-hook: glyco-id: ~s~%" glyco-id)
     ;; if it was an ASP create a level-0 glyco-id for that (glyco-tree-residue-id doesn't
     ;; do that (not sure why)).
     (if (not glyco-id)
	 (let ((rn (residue-name aa-imol aa-chain-id aa-res-no aa-ins-code)))
	   (if (string? rn)
	       (if (string=? rn "ASN")
		   (set! glyco-id (list 0 'unset "ASN" "" "" aa-res-spec))))))
     (if (list? glyco-id)
	 (let ((tree-type (get-tree-type)))
	   (let ((children (gtk-container-children vbox)))
	     (for-each (lambda (child)
			 (if (gtk-button? child)
			     (glyco-tree-dialog-set-button-active-state child glyco-id tree-type)))
		       children)
	     #t)
	   #f)))))

(define (add-module-carbohydrate-gui)

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "Glyco")))

	(add-simple-coot-menu-menuitem
	 menu "N-linked Glycan Addition..."
	 (lambda ()
	   (interactive-add-cho-dialog)))

	(add-simple-coot-menu-menuitem
	 menu "Set Default N-linked CHO Atoms B-factor"
	 (lambda ()
	   (using-active-atom
	    (let ((residues (residues-near-residue aa-imol aa-res-spec 10)))
	      (let ((imol-region (new-molecule-by-residue-specs aa-imol residues)))
		(let ((m (median-temperature-factor imol-region)))
		  (close-molecule imol-region)
		  (if (number? m)
		      (let ((new-m  (* m 1.55)))
			(set-default-temperature-factor-for-new-atoms new-m)
			(let ((s (format #f "New Temperature Factor set to ~6,2f" new-m)))
			  (info-dialog s))))))))))

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
	 menu "Add High Mannose"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    ;; (with-no-backups aa-imol
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     oligomannose-tree))))

	(add-simple-coot-menu-menuitem
	 menu "Add Hybrid (Mammal)"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    ;; (with-no-backups aa-imol
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     hybrid-mammal-tree))))

; in practice, no one will be doing this.
; 	(add-simple-coot-menu-menuitem
; 	 menu "Add Hybrid (Plant)"
; 	 (lambda ()
; 	   (using-active-atom
; 	    (make-backup aa-imol)
; 	    ;; (with-no-backups aa-imol
; 	    (add-linked-residue-tree aa-imol
; 				     (list aa-chain-id aa-res-no aa-ins-code)
; 				     hybrid-plant-derived-tree))))

	(add-simple-coot-menu-menuitem
	 menu "Add Complex (Mammal)"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     complex-mammal-tree))))

	(add-simple-coot-menu-menuitem
	 menu "Add Complex (Plant)"
	 (lambda ()
	   (using-active-atom
	    (make-backup aa-imol)
	    (add-linked-residue-tree aa-imol
				     (list aa-chain-id aa-res-no aa-ins-code)
				     complex-plant-tree))))

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

	(add-simple-coot-menu-menuitem
	 menu "Extract this Tree"
	 (lambda()
	   (new-molecule-from-this-glyco-tree)))

	)))
