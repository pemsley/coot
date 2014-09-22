
(define *ligand-check-refmac-mtz*         #f)
(define *ligand-check-refmac-fobs-col*    #f)
(define *ligand-check-refmac-sigfobs-col* #f)
(define *ligand-check-refmac-r-free-col*  #f)

(define (ligand-check-refmac-columns f-list sigf-list rfree-list)

  (cond

   ((all-true?
     (string? *ligand-check-refmac-fobs-col*)
     (string? *ligand-check-refmac-sigfobs-col*)
     (string? *ligand-check-refmac-r-free-col*))
    
    (list *ligand-check-refmac-fobs-col*
	  *ligand-check-refmac-sigfobs-col*
	  *ligand-check-refmac-r-free-col*))

   ((not (> (vector-length f-list) 0))     #f)
   ((not (> (vector-length sigf-list) 0))  #f)
   ((not (> (vector-length rfree-list) 0)) #f)

   (else 
    ;; happy path
    
    ;; Using the first sigf (there should only be one typically)
    ;; get the F label (Fx) from x/y/SIGFPx
    ;; 
    (let ((field (vector-ref sigf-list 0)))
      (let ((ls (split-after-char-last #\/ field list)))
	(let ((sigf-col (cadr ls)))
	  (let ((sm (string-match "SIGF" sigf-col)))
	    (let ((f-col (substring sigf-col (cdr (vector-ref sm 1)))))
	      (list 
	       (string-append (car ls) "F" f-col)
	       field
	       (vector-ref rfree-list 0))))))))))


(if (enhanced-ligand-coot?)

    (begin

      (set! import-from-3d-generator-from-mdl 
	    (lambda (mdl-file-name comp-id)

	      (let ((status (goosh-command
			     "pyrogen"
			     (list "-m" mdl-file-name "-t" comp-id)
			     '()
			     (string-append "pyrogen-" comp-id ".log")
			     #f)))
		(format #t "goosh status for pyrogen: ~s~%" status)
		(if (ok-goosh-status? status)
		    (let ((pdb-name (string-append comp-id "-pyrogen.pdb"))
			  (cif-name (string-append comp-id "-pyrogen.cif")))
		      
		      (read-pdb pdb-name)
		      (read-cif-dictionary cif-name))))))

      (if (defined? 'coot-main-menubar)
	  (let ((menu (coot-menubar-menu "Ligand")))
	    (add-simple-coot-menu-menuitem 
	     menu "SMILES -> 2D"
	     (lambda ()
	       (generic-single-entry "SMILES string" "" " Send to 2D Viewer " 
				     (lambda (text)
				       (smiles-to-ligand-builder text)))))

	    (add-simple-coot-menu-menuitem 
	     menu "Residue -> 2D"
	     (lambda ()
	       (using-active-atom 
		(residue-to-ligand-builder aa-imol aa-chain-id aa-res-no aa-ins-code 0.015))))

	    (add-simple-coot-menu-menuitem
	     menu "FLEV this residue"
	     (lambda () 
	       (using-active-atom 
		(fle-view-with-rdkit aa-imol aa-chain-id aa-res-no aa-ins-code 4.2)
		(set-flev-idle-ligand-interactions 1))))
	    
	    (add-simple-coot-menu-menuitem
	     menu "Toggle FLEV Ligand Interactions" 
	     (lambda ()
	       (toggle-flev-idle-ligand-interactions)))

	    (add-simple-coot-menu-menuitem
	     menu "Show Chemical Features"
	     (lambda () 
	       (set-display-generic-objects-as-solid 1) ;; there may be consequences...
	       (using-active-atom
		(show-feats aa-imol aa-chain-id aa-res-no aa-ins-code))))))



))


