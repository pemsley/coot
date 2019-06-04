
(define (add-pyranose-pseudo-ring-plane-restraints comp-id)

  (define (filter-out plane-name-sub-string plane-restraints)
    (filter (lambda (s)
	      (let ((v (regexp-match? (string-match plane-name-sub-string (car s)))))
		(not v))) ;; filter function keeps items that for func(s) is true
	    plane-restraints))

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
;; (define (add-cho-restraints-for-residue imol new-res-spec)
;; #f)

(define (multi-add-linked-residue imol res-spec residues-to-add)

  (format #t "---------------- multi-add-linked-residue ~s ~s ~%~!" imol res-spec)
  (set-go-to-atom-molecule imol)
  (let ((wm (matrix-state)))
    (set-matrix (/ wm 4))

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
	(set-dragged-refinement-steps-per-frame current-refinement-rate)
	(set-matrix wm)))))

;; return the new molecule number
(define (new-molecule-from-this-glyco-tree)
  (using-active-atom
   (let* ((tree-residues (glyco-tree-residues aa-imol aa-res-spec)))
     (new-molecule-by-residue-specs aa-imol tree-residues))))


;; also "precursor"
;; high mannose was used for human
;; high mannose is now used for human too
;; call this "High Mannose" in the GUI
;;
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
				(("ALPHA1-2" . "MAN")
				 (("ALPHA1-3" . "GLC")
				  (("ALPHA1-3" . "GLC")
				   (("ALPHA1-2" . "GLC")))))))))))

;; Hybrid is for any system
;;
;; Plant Hybrid also allows an alpha1-3 FUC
;;

;; hybrid mammal
;;
(define hybrid-mammal-tree  '(("NAG-ASN" . "NAG")
			     (("BETA1-4" . "NAG")
			      (("BETA1-4" . "BMA")
			       (("ALPHA1-6" . "MAN")
				(("ALPHA1-6" . "MAN"))
				(("ALPHA1-3" . "MAN")))
			       (("ALPHA1-3" . "MAN")
				(("BETA1-2" . "NAG")
				 (("BETA1-4" . "GAL")
				  (("ALPHA2-3" . "SIA"))
				  (("ALPHA2-6" . "SIA")))))
			       (("BETA1-4" "NAG"))))
			     ("ALPHA1-6" . "FUC")))


;; hybrid plant
;;
(define hybrid-plant-derived-tree  '(("NAG-ASN" . "NAG")
				     (("BETA1-4" . "NAG")
				      (("BETA1-4" . "BMA")
				       (("ALPHA1-6" . "MAN")
					(("ALPHA1-6" . "MAN"))
					(("ALPHA1-3" . "MAN")))
				       (("ALPHA1-3" . "MAN")
					(("BETA1-2" . "NAG")
					 (("BETA1-4" . "GAL")
					  (("ALPHA2-3" . "SIA"))
					  (("ALPHA2-6" . "SIA")))))
				       (("XYP-BMA"  . "XYP"))
				       (("BETA1-4" "NAG"))))
				     ("ALPHA1-6" . "FUC")
				     ("ALPHA1-3" . "FUC")))


;; complex mammal
;; bianntennary mammal
(define complex-mammal-tree  '(("NAG-ASN" . "NAG")
			       (("BETA1-4" . "NAG")
				(("BETA1-4" . "BMA")
				 (("ALPHA1-6" . "MAN")
				  (("BETA1-2" . "NAG")
				   (("BETA1-4" . "GAL")
				    (("ALPHA2-3" . "SIA"))
				    (("ALPHA2-6" . "SIA")))))
				 (("ALPHA1-3" . "MAN")
				  (("BETA1-2" . "NAG")
				   (("BETA1-4" . "GAL")
				    (("ALPHA2-3" . "SIA"))
				    (("ALPHA2-6" . "SIA")))))
				 (("BETA1-4" "NAG"))))
			       ("ALPHA1-6" . "FUC")))


;; complex plant
;; plant bianntennary
(define complex-plant-tree  '(("NAG-ASN" . "NAG")
			      (("BETA1-4" . "NAG")
			       (("BETA1-4" . "BMA")
				(("ALPHA1-6" . "MAN")
				 (("BETA1-2" . "NAG")
				  (("BETA1-4" . "GAL")
				   (("ALPHA2-3" . "SIA"))
				   (("ALPHA2-6" . "SIA")))))
				(("ALPHA1-3" . "MAN")
				 (("BETA1-2" . "NAG")
				  (("BETA1-4" . "GAL")
				   (("ALPHA2-3" . "SIA"))
				   (("ALPHA2-6" . "SIA")))))
				(("BETA1-4" "NAG"))
				(("XYP-BMA" "XYP"))))
			      (("ALPHA1-6" . "FUC"))
			      (("ALPHA1-3" . "FUC"))))


;; (set-add-linked-residue-do-fit-and-refine 0)

;; now users can set this
(define *add-linked-residue-tree-correlation-cut-off* 0.50)

(define (add-linked-residue-add-cho-function imol parent res-pair)

  (define (well-fitting? res-spec)
    (using-active-atom
     (let ((neighbs (residues-near-residue aa-imol res-spec 4)))
       (let ((c (map-to-model-correlation imol (list res-spec) neighbs 0 (imol-refinement-map))))
	 (format #t "############# new residue ~s density correlation: ~s~%" res-spec c)
	 (if (not (number? c))
	     #f
	     (if (> c *add-linked-residue-tree-correlation-cut-off*)
		 (let ((symm-clash (clashes-with-symmetry imol
							  (residue-spec->chain-id res-spec)
							  (residue-spec->res-no   res-spec)
							  (residue-spec->ins-code res-spec) 2.0)))
		   (if (= symm-clash 1) #f #t))
		 #f))))))

  (define (centre-view-on-residue-centre res-spec)
    (let ((res-centre (residue-centre imol
				      (residue-spec->chain-id res-spec)
				      (residue-spec->res-no res-spec)
				      (residue-spec->ins-code res-spec))))
      (if (list? res-centre)
	  (apply set-rotation-centre res-centre))))

  (define (delete-residue-by-spec spec)
    (delete-residue imol
		    (residue-spec->chain-id spec)
		    (residue-spec->res-no   spec)
		    (residue-spec->ins-code spec)))

  ;; main line
  ;;
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


(define (add-linked-residue-tree imol parent tree)

  (define func-test
    (let ((count 10))
      (lambda (parent res-pair)
	(let ((rv (list "B" count "")))
	  (set! count (+ count 1))
	  (format #t "process res-pair ~s with parent ~s producing ~s ~%" res-pair parent rv)
	  rv))))

  (define (process-tree parent tree proc-func)
    (cond
     ((null? tree) '())
     ((list? (car tree))
      (let ((part-1 (process-tree parent (car tree) proc-func))
	    (part-2 (process-tree parent (cdr tree) proc-func)))
	(cons part-1 part-2)))
     (else
      (let ((new-res (proc-func imol parent (car tree))))
	(cons new-res
	      (process-tree new-res (cdr tree) proc-func))))))

  (define (is-just-an-ASN? imol glyco-tree)

    (format #t "glyco-tree: ~s~%" glyco-tree)
    (if (not (list? glyco-tree))
	#f
	(let ((l (length glyco-tree)))
	  (format #t "l: ~s~%" l)
	  (if (not (= l 1))
	      #f
	      (using-active-atom
	       (let ((res-spec aa-res-spec))
		 (format #t "res-spec: ~s~%" res-spec)
		 (let ((rn (residue-name imol
					 (residue-spec->chain-id res-spec)
					 (residue-spec->res-no   res-spec)
					 (residue-spec->ins-code res-spec))))
		   (format #t "--------- Here with rn: ~s~%" rn)
		   (if (not (string? rn))
		       #f
		       (string=? rn "ASN")))))))))


  ;; main line of add-linked-residue-tree
  ;;
  (add-synthetic-pyranose-planes)
  (use-unimodal-pyranose-ring-torsions)
  (set-refine-with-torsion-restraints 1)
  (let ((wm (matrix-state)))
    (set-matrix (/ wm 4))
    (set-residue-selection-flash-frames-number 1)
    (set-go-to-atom-molecule imol)
    (set-go-to-atom-from-res-spec parent)
    (let* ((previous-m (default-new-atoms-b-factor))
	   (m (median-temperature-factor imol))
	   (new-m (* m 1.55)))

      (set-default-temperature-factor-for-new-atoms previous-m)

      (if (number? m)
	  (set-default-temperature-factor-for-new-atoms new-m))

      ;; prevent tree building if there is already a partial tree here
      ;; (only proceed with the one ASN)
      ;;
      (using-active-atom
       (let ((start-tree (glyco-tree-residues aa-imol aa-res-spec)))
	 
	 ;; why do I need both test here? 
	 ;; 5n11 needs the is-just-an-ASN? test
	 ;; 5n09/5wzy needs the null test. 
	 ;; Hmmm.
	 (if (not (or (null? start-tree)
		      (is-just-an-ASN? aa-imol start-tree)))

	     (begin
	       (info-dialog "Must start on Single ASN")
	       (format #t "start-tree: ~s~%" start-tree))

	     ;; OK, continue
	     (let ((start-pos-view (add-view-here "Glyo Tree Start Pos")))
	       (process-tree parent tree add-linked-residue-add-cho-function)
	       (go-to-view-number start-pos-view 0)
	       (with-auto-accept (refine-residues aa-imol (glyco-tree-residues aa-imol aa-res-spec)))
	       ;; add a test here that the tree here (centre of screen) matches a known tree.
	       ;;
	       ;; and that each is 4C1 (or 1C4 for FUC?) (XYP?)
	       )))))
    (set-matrix wm)))


(define (add-linked-residue-with-extra-restraints-to-active-residue new-res-type link-type)
  (let ((wm (matrix-state)))
    (set-matrix (/ wm 8))
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
		(refine-residues aa-imol residues)))))))
    (set-matrix wm)))

(define (delete-all-cho)
  (let ((delete-cho-list '()))
    (using-active-atom
     (if (valid-model-molecule? aa-imol)
	 (begin
	   (for-each (lambda (chain-id)
		       (for-each (lambda (res-serial)
				   (let ((res-no (seqnum-from-serial-number aa-imol chain-id res-serial))
					 (ins-code (insertion-code-from-serial-number aa-imol chain-id res-serial)))
				     (let ((rn (residue-name aa-imol chain-id res-no ins-code)))
				       (if (string? rn)
					   ;; a better test is to find all the hetgroups and look at the _chem_comp group or type
					   (if (or (string=? "NAG" rn) (string=? "MAN" rn) (string=? "BMA" rn) (string=? "FUL" rn)
						   (string=? "FUC" rn) (string=? "XYP" rn) (string=? "SIA" rn) (string=? "GAL" rn)
						   (string=? "NDG" rn) (string=? "BGC" rn) (string=? "A2G" rn))
					       (let* ((residue-spec (list chain-id res-no ins-code)))
						 (set! delete-cho-list (cons residue-spec delete-cho-list))))))))
				 (range (chain-n-residues chain-id aa-imol))))
		     (chain-ids aa-imol))

	   ;; now we have delete-residues, we don't need to delete them one by one
	   ;;(for-each (lambda(cho-res-spec)
	   ;; (delete-residue aa-imol (residue-spec->chain-id cho-res-spec) (residue-spec->res-no cho-res-spec) ""))
	   ;;   delete-cho-list)))))))
	   ;;
	   (delete-residues aa-imol delete-cho-list))))))
