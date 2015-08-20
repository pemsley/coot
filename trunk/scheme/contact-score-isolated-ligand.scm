

(define (deactivate-molecules-except imol)
  (for-each (lambda (i)
	      (if (not (= i imol))
		  (begin
		    (set-mol-active i 0))))
	    (model-molecule-list)))


;; This will hydrogenate the active residue, not imol
;; 
(define (contact-score-ligand imol res-spec)

  (deactivate-molecules-except imol)

  (let* ((chain-id (res-spec->chain-id res-spec))
	 (res-no   (res-spec->res-no   res-spec))
	 (ins-code (res-spec->ins-code res-spec))
	 (ss (string-append "//" chain-id "/" (number->string res-no)))
	 (imol-selection (new-molecule-by-atom-selection imol ss))
	 (coot-molprobity-dir (get-directory "coot-molprobity"))
	 (ligand-selection-pdb  (append-dir-file coot-molprobity-dir (string-append "tmp-selected-ligand-for-probe-" 
										  (number->string imol)
										  ".pdb")))
	 (protein-selection-pdb (append-dir-file coot-molprobity-dir (string-append "tmp-protein-for-probe-" 
										  (number->string imol)
										  ".pdb")))
	 (dots-file-name (append-dir-file coot-molprobity-dir
					  (string-append "probe"
							 "-"
							 chain-id
							 "-"
							 (number->string res-no)
							 ".dots"))))
    
    (set-mol-active    imol-selection 0)
    (set-mol-displayed imol-selection 0)
    (set-go-to-atom-molecule imol)
    (let ((rc (residue-centre imol chain-id res-no ins-code)))
      (apply set-rotation-centre rc)
      (hydrogenate-region 6)

      (write-pdb-file imol-selection ligand-selection-pdb)
      (write-pdb-file imol protein-selection-pdb)
      (goosh-command
       *probe-command*
       (list "-q" "-u" "-once" ;; -once or -both
             ;; first pattern
             (string-append "CHAIN" chain-id " "(number->string res-no))
	     ;; second pattern
             ;; (string-append "not " (number->string res-no))
             (string-append "not " (number->string res-no))
	     ;; consider og33 (occ > 0.33)
	     "-density30"
	     ligand-selection-pdb protein-selection-pdb)
       '()
       dots-file-name
       #f)


      ;; debugging
      (handle-read-draw-probe-dots-unformatted dots-file-name imol 0)
      
      (probe-clash-score-scm dots-file-name))))

(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Ligand")))
      (add-simple-coot-menu-menuitem 
       menu "Isolated dots for this ligand"
       (lambda ()
        (using-active-atom
	   (contact-score-ligand aa-imol (list aa-chain-id aa-res-no aa-ins-code)))))))




; ;; test it
; (let* ((imol (read-pdb "1x8b-bumpy.pdb"))
;        (result (contact-score-ligand imol (list "A" 901 ""))))

;   (format #t "result: ~s~%" result))

