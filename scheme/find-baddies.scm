
;; find bad things in the structure - rama, C-beta, rotamer, atom clashes baddies
;;
(define (validation-outliers-dialog imol imol-map)

  (let ((dialog-vbox #f)
	(window #f)
	(missing-sidechains-checkbutton #f)
	(cg-torsion-diff-checkbutton #f)
	(poor-density-checkbutton #f))

    (define (find-rama-baddies)
      (let ((rs (all-molecule-ramachandran-score imol)))
	;; (format #t "rs: ~s~%" rs)
	(let ((scored-residues (list-ref rs 5)))
	  (let ((interesting (filter (lambda(item)
				       (let ((phi-psi      (list-ref item 0))
					     (residue-spec (list-ref item 1))
					     (pr           (list-ref item 2))
					     (res-names    (list-ref item 3)))
					 (< pr 0.02)))
				     scored-residues)))

	    ;; remove phi-psi and res-names from the return value
	    (let ((munged (map (lambda(item)
				 (list (list-ref item 1)
				       (list-ref item 2)))
			       interesting)))
	      (list "Ramachandran Improbables" munged))))))

    ;; a list of atom specs
    (define (find-chiral-volume-baddies)
      (let ((r (chiral-volume-errors imol)))
	(if (not (list? r))
	    '()
	    r)))

    (define (make-window-title n)
      (string-append "Coot Interesting/Outliers/Problems: "
		     (number->string (- n 1))))

    (define (find-c-beta-baddies)
      (if (defined? 'c-beta-deviations)
	  (c-beta-deviations imol)
	  '()))

    (define (find-em-ringer-baddies)
      (if (not (defined? 'CG-spin-search))
	  '()
	  (let ((scored-residues (CG-spin-search imol imol-map)))
	    (if (not (list? scored-residues))
		'()
		(let ((interesting (filter (lambda(item)
					     (let ((delta (cadr item)))
					       (or (< delta -30)
						   (> delta 30))))
					   scored-residues)))
		  ;; return a list of residue specs - not what is expected?
		  interesting)))))

    (define (rotamer-score-residues imol)

      (let ((residues (all-residues-sans-water imol)))
	(map (lambda(residue-spec)
	       (let ((alt-conf ""))
		 (let ((score (rotamer-score imol
					     (residue-spec->chain-id residue-spec)
					     (residue-spec->res-no   residue-spec)
					     (residue-spec->ins-code residue-spec)
					     alt-conf)))
		   (list residue-spec score))))
	     residues)))

    (define (filter-rotamer-baddies baddies)

      (let ((het-groups-in-mol (het-group-residues imol)))

	(filter (lambda(baddie)
		  (let ((spec (car baddie))
			(score (cadr baddie)))
		    (let ((res-name (residue-name imol
						  (residue-spec->chain-id spec)
						  (residue-spec->res-no   spec)
						  (residue-spec->ins-code spec))))

		      ;; (format #t "filter-rotamers testing baddie: ~s~%" baddie)

		      (cond
		       ((string=? res-name "ALA") #f)
		       ((string=? res-name "GLY") #f)
		       ((string=? res-name "UNK") #f)
		       ((string=? res-name "HOH") #f)
		       ((and (= score 0.0) (not (ok-to-have-missing-sidechain-buttons?))) #f)
		       (else
			;; if spec is a het-group then no rotamers for that (return #f)
			(let ((is-het (find (lambda(item)
					      (residue-specs-match? item spec)) het-groups-in-mol)))
			  ;; (format #t "spec: ~s res-name ~s het-groups-in-mol ~s is-het ~s ~%"
			  ;;            spec res-name het-groups-in-mol is-het)
			  (if is-het
			      #f
			      (if (string=? res-name "LEU")
				  (begin
				    ;; there is something strange with LEU rotamers, D 427?
				    ;; (format #t "------------- here with spec: ~s res-name ~s score ~s ~%"
				    ;;         res-name spec score)
				    (< score 1.1))
				  (< score 0.021)))))))))
		baddies)))

    (define (molecule-atom-overlap-baddies)
      (molecule-atom-overlaps imol))

    (define (filter-molecule-atom-overlap-baddies mao-baddies)
      (let ((baddie-limit 2.2)) ;; more than this is marked as a baddie
	(let ((fn (lambda(mao-item)
		    (let ((atom-spec-1 (list-ref mao-item 0))
			  (atom-spec-2 (list-ref mao-item 1))
			  (overlap     (list-ref mao-item 4)))
		      (> overlap baddie-limit)))))
	  (filter fn mao-baddies))))

    (define (non-pro-cis-peptide-baddies) ;; do the filter here - just for consistency
      (let ((cis-peps (cis-peptides imol)))
	(filter (lambda(peptide)
		    (format #t "cis peptide: ~s~%" peptide)
		    (let ((spec-1 (list-ref peptide 0))
			  (spec-2 (list-ref peptide 1))
			  (omega  (list-ref peptide 2)))
		      (let ((rn (residue-spec->residue-name imol spec-2)))
			(not (string=? rn "PRO")))))
		  cis-peps)))

    (define (destroy-buttons-with-label label-fragment-txt dialog-vbox)
      (let ((current-buttons (gtk-container-children dialog-vbox)))
	(for-each (lambda (button)
		    (let ((label (gtk-button-get-label button)))
		      (if (string-match label-fragment-txt label)
			  (gtk-widget-destroy button))))
		  current-buttons)))

    (define (regenerate-button-fn)

      (if dialog-vbox
	  (begin
	    (let ((buttons (make-buttons)))
	      (let ((old-buttons (gtk-container-children dialog-vbox)))
		(for-each (lambda(button-spec)
			    (let ((button (gtk-button-new-with-label
					   (car button-spec))))
			      (gtk-signal-connect
			       button "clicked"
			       (cadr button-spec))
			      (gtk-box-pack-start
			       dialog-vbox button #f #f 2)
			      (gtk-widget-show button)))
			  buttons)
		(if window
		    (gtk-window-set-title
		     window
		     (make-window-title (length buttons))))
		(for-each gtk-widget-destroy old-buttons))))
	  (validation-outliers-dialog imol imol-map)))

    (define (ok-to-do-density-correlations?)
      (if poor-density-checkbutton
	  (gtk-toggle-button-get-active poor-density-checkbutton)
	  #t))

    (define (ok-to-have-missing-sidechain-buttons?)
      (if missing-sidechains-checkbutton
	  (gtk-toggle-button-get-active missing-sidechains-checkbutton)
	  #t))

    (define (ok-to-do-CG-torsion-diffs?)
      (if (not (defined? 'CG-spin-search))
	  #f
	  (gtk-toggle-button-get-active cg-torsion-diff-checkbutton)))

    (define (make-buttons)

      (let ((frb (find-rama-baddies))
	    (fcbb (find-c-beta-baddies))
	    (filtered-mao-baddies (filter-molecule-atom-overlap-baddies (molecule-atom-overlap-baddies)))
	    (residue-correlations
	     (if (not (ok-to-do-density-correlations?))
		 '()
		 (map-to-model-correlation-per-residue imol (all-residues-sans-water imol) 0 imol-map))))

	(let* (

	       ;; rama
	       (rama-filter-fn (lambda(baddie)
			    (let ((spec (car baddie))
				  (rama-prob (cadr baddie)))
			      (< rama-prob 0.002)))) ;; tested
	       (baddies (filter rama-filter-fn (cadr frb)))
	       (sorted-filtered-rama-baddies (sort baddies (lambda(ele-1 ele-2)
							     (let ((rama-prob-1 (list-ref ele-1 1))
								   (rama-prob-2 (list-ref ele-2 1)))
							       (< rama-prob-1 rama-prob-2)))))

	       ;; c-beta
	       (c-beta-filter-fn (lambda(baddie)
				   (let ((spec (car baddie))
					 (alt-conf-map (cadr baddie))
					 (score (cadr (car (cadr baddie))))) ;; only the first score (lazy?)
				     (> score 0.25)))) ;; pretty close
	       (c-beta-baddies (filter c-beta-filter-fn fcbb))
	       (sorted-filtered-c-beta-baddies (reverse (sort c-beta-baddies
							      (lambda(ele-1 ele-2)
								(let ((spec-1 (car ele-1))
								      (spec-2 (car ele-2))
								      (score-1 (cadr (car (cadr ele-1))))
								      (score-2 (cadr (car (cadr ele-2)))))
								  (< score-1 score-2))))))

	       ;; Density correlations
	       ;;
	       (density-baddies-filter-fn (lambda(baddie)
					    (let ((correlation (cadr baddie)))
					      (< correlation 0.8))))
	       (density-baddies (filter density-baddies-filter-fn residue-correlations))

	       ;; CG Torsion
	       ;;
	       (cg-torsion-baddies (find-em-ringer-baddies))

	       ;; Rotamers
	       ;;
	       (rotamer-baddies (rotamer-score-residues imol))
	       (filtered-rotamer-baddies (filter-rotamer-baddies rotamer-baddies))
	       (sorted-filtered-rotamer-baddies (sort filtered-rotamer-baddies
						      (lambda (ele-1 ele-2)
							(let ((score-1 (cadr ele-1))
							      (score-2 (cadr ele-2)))
							  ;; score things with score 0.0 (meaning missing sidechain)
							  ;; as if they are better than low probability outliers
							  ;;
							  (cond
							   ((and (= score-1 0.0)
								 (> score-2 0.0)) #f)
							   ((and (= score-2 0.0)
								 (> score-1 0.0)) #t)
							   (else
							    (< score-1 score-2))))))))

	  (let ((rama-buttons (map (lambda (baddie)
				     (let ((spec (car baddie))
					   (rama-prob (cadr baddie)))
				       (let* ((score-string
					       (format #f "~5f %" (* 100 rama-prob)))
					      (button-label (string-append
							     "Ramachandran Outlier "
							     (residue-spec->chain-id spec)
							     " "
							     (number->string (residue-spec->res-no spec))
							     (residue-spec->ins-code spec)
							     " "
							     (residue-spec->residue-name imol spec)
							     " "
							     score-string))
					      (fn (lambda()
						    (set-go-to-atom-molecule imol)
						    (set-go-to-atom-from-res-spec spec))))
					 (list button-label fn))))
				   sorted-filtered-rama-baddies))

		(c-beta-buttons (map (lambda(baddie)
				       (let ((spec (car baddie))
					     (score (cadr (car (cadr baddie))))) ;; only the first score
					 (let ((score-string (format #f "~5f" score)))
					   (let ((button-label
						  (string-append
						   "C-beta deviant "
						   (residue-spec->string spec)
						   "  "
						   (residue-spec->residue-name imol spec)
						   " "
						   score-string "Å"))
						 (fn (lambda()
						       (set-go-to-atom-molecule imol)
						       (set-go-to-atom-from-res-spec spec))))
					     (list button-label fn)))))
				     sorted-filtered-c-beta-baddies))

		(non-pro-cis-peptide-buttons (map (lambda (baddie)
						    (let ((spec-1 (list-ref baddie 0))
							  (spec-2 (list-ref baddie 1))
							  (omega  (list-ref baddie 2)))
						      (let ((button-label (string-append
									   "Non-PRO cis-peptide "
									   (residue-spec->string spec-1)
									   " - "
									   (residue-spec->string spec-2)))
							    (fn (lambda()
								  (set-go-to-atom-molecule imol)
								  (set-go-to-atom-from-res-spec spec-1))))
							(list button-label fn))))
						  (non-pro-cis-peptide-baddies)))

		(rota-buttons (map (lambda(baddie)
				     (let ((spec (car baddie))
					   (score (cadr baddie)))

				       ;; I am not sure that I like a score of 0.0 meaning "Missing sidechain"
				       ;; we have lost some information on the way
				       ;;
				       (let ((score-string (format #f "~5f %" (* score 100)))
					     (ms-string (if (= score 0.0) "Missing Sidechain" "Rotamer Outlier"))
					     (rot-name (get-rotamer-name imol
									 (residue-spec->chain-id spec)
									 (residue-spec->res-no   spec)
									 (residue-spec->ins-code spec))))
					 (let ((button-label
						(string-append
						 ms-string " "
						 (residue-spec->string spec)
						 " "
						 (residue-spec->residue-name imol spec)
						 " "
						 (if (string? rot-name) rot-name)
						 "  "
						 (if (= score 0.0) ""
						     score-string)))
					       (fn (lambda()
						     (set-go-to-atom-molecule imol)
						     (set-go-to-atom-from-res-spec spec))))
					   (list button-label fn)))))
				   sorted-filtered-rotamer-baddies))

		(density-baddies-buttons (map (lambda(baddie)
						(let ((spec (car baddie))
						      (score (cadr baddie)))
						  (let ((button-label (string-append
								       "Poor Density Fit "
								       (residue-spec->string spec)
								       " "
								       (format #f "~5f" score))))
						    (let ((fn (lambda()
								(set-go-to-atom-molecule imol)
								(set-go-to-atom-from-res-spec spec))))
						      (list button-label fn)))))
					      density-baddies))

		(cg-torsion-buttons (map (lambda(baddie)
					   (let ((spec (car baddie))
						 (score (cadr baddie)))
					     ;; (format #t "cg score: ~s~%" score)
					     (let ((button-label (string-append
								  "CG Torsion Diff. "
								  (residue-spec->string spec)
								  " "
								  (format #f "~5f" score))))
					       (let ((fn (lambda()
							   (set-go-to-atom-molecule imol)
							   (set-go-to-atom-from-res-spec spec))))
						 (list button-label fn)))))
					 cg-torsion-baddies))

		(chiral-volume-buttons (map (lambda (baddie-atom-spec)
					      (let ((button-label
						     (string-append "Chiral Volume Error "
								    (atom-spec->string baddie-atom-spec)))
						    (fn (lambda ()
							   (set-go-to-atom-molecule imol)
							   (set-go-to-atom-from-atom-spec baddie-atom-spec))))
						(list button-label fn)))
					    (find-chiral-volume-baddies)))

		(atom-overlap-buttons (map (lambda(baddie)
					     (let ((atom-spec-1 (cdr (list-ref baddie 0))) ;; unprefix
						   (atom-spec-2 (cdr (list-ref baddie 1))) ;; ditto
						   (overlap     (list-ref baddie 4)))
					       (let ((buton-label
						      (string-append
						       "Atom Overlap "
						       (atom-spec->string atom-spec-1)
						       " on "
						       (atom-spec->string atom-spec-2)
						       " OV: "
						       (format #f "~5f" overlap)))
						     (fn (lambda()
							   (set-go-to-atom-molecule imol)
							   (set-go-to-atom-from-atom-spec atom-spec-1))))
						 (list buton-label fn))))
					   filtered-mao-baddies)))

	    (let ((buttons (append chiral-volume-buttons
				   rama-buttons
				   rota-buttons
				   non-pro-cis-peptide-buttons
				   density-baddies-buttons
				   c-beta-buttons
				   cg-torsion-buttons
				   atom-overlap-buttons)))
	      buttons)))))

  ;; --- main line ---

    (let* ((buttons (make-buttons)))

      (let ((p (dialog-box-of-buttons (make-window-title (length buttons)) (cons 350 400) buttons " Close ")))
      (set! dialog-vbox (car p))
      (set! window (cadr p))

      (let ((window-bits (gtk-container-children window)))
	(let ((vbox-outer (car window-bits)))
	  (let ((control-button-vbox-1 (gtk-hbox-new #f 2)))
	    (let ((missing-sidechains-checkbutton-local (gtk-check-button-new-with-label "Missing Sidechains"))
		  (poor-density-checkbutton-local (gtk-check-button-new-with-label "Poor Density Fit"))
		  (cg-torsion-diff-checkbutton-local (gtk-check-button-new-with-label "CG Torsion Diff."))
		  (regenerate-button-local (gtk-button-new-with-label "Update")))

	      (set! missing-sidechains-checkbutton missing-sidechains-checkbutton-local)
	      (set! poor-density-checkbutton poor-density-checkbutton-local)
	      (set! cg-torsion-diff-checkbutton cg-torsion-diff-checkbutton-local)

	      (gtk-toggle-button-set-active missing-sidechains-checkbutton #t)
	      (gtk-toggle-button-set-active poor-density-checkbutton       #t)
	      (gtk-toggle-button-set-active cg-torsion-diff-checkbutton    #t)
	      (gtk-box-pack-start control-button-vbox-1 missing-sidechains-checkbutton #f #f 2)
	      (gtk-box-pack-start control-button-vbox-1 poor-density-checkbutton #f #f 2)
	      (gtk-box-pack-start vbox-outer regenerate-button-local #f #f 6)

	      (gtk-box-pack-start vbox-outer control-button-vbox-1 #f #f 2)

	      (gtk-signal-connect regenerate-button-local "clicked" regenerate-button-fn)

	      (gtk-signal-connect missing-sidechains-checkbutton "toggled"
				  (lambda()
				    (let ((state (gtk-toggle-button-get-active missing-sidechains-checkbutton)))
				      (if (not state) ;; i.e. no buttons with "Missing Sidechain"
					  (destroy-buttons-with-label "Missing Sidechain" dialog-vbox)))))

	      (gtk-signal-connect poor-density-checkbutton "toggled"
				  (lambda()
				    (let ((state (gtk-toggle-button-get-active poor-density-checkbutton)))
				      (if (not state)
					  (destroy-buttons-with-label "Poor Density" dialog-vbox)))))

	      ;; 20190102-PE depends on the version of coot that we are using
	      ;;
	      (if (defined? 'CG-spin-search)
		  (begin
		    (gtk-box-pack-start control-button-vbox-1 cg-torsion-diff-checkbutton #f #f 2)
		    (gtk-widget-show cg-torsion-diff-checkbutton)
		    (gtk-signal-connect cg-torsion-diff-checkbutton "toggled"
				  (lambda()
				    (let ((state (gtk-toggle-button-get-active cg-torsion-diff-checkbutton)))
				      (if (not state)
					  (destroy-buttons-with-label "CG Torsion" dialog-vbox)))))))

	      (gtk-widget-show control-button-vbox-1)
	      (gtk-widget-show missing-sidechains-checkbutton)
	      (gtk-widget-show poor-density-checkbutton)
	      (gtk-widget-show regenerate-button-local)
	      ))))

      ))))

