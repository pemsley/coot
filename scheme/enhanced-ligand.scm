;;;; Copyright 2013 by Medical Research Council

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


;;(if (enhanced-ligand-coot?)
;;    (format #t "----------------------------- with enhanced-ligand! -------------------~%~!")
;;    (format #t "----------------------------- without enhanced-ligand! -------- -------~%~!"))


(if (enhanced-ligand-coot?)

    (begin

      ;; Use pyrogen if we have mogul
      ;; 
      (if *use-mogul*
	  (set! import-from-3d-generator-from-mdl 
		(lambda (mdl-file-name comp-id)

		  (let ((status (goosh-command
				 "pyrogen"
				 (list "--mol" mdl-file-name "--residue-type" comp-id)
				 '()
				 (string-append "pyrogen-" comp-id ".log")
				 #f)))
		    (format #t "goosh status for pyrogen: ~s~%" status)
		    (if (ok-goosh-status? status)
			(let ((pdb-name (string-append comp-id "-pyrogen.pdb"))
			      (cif-name (string-append comp-id "-pyrogen.cif")))

			  (read-pdb pdb-name)
			  (read-cif-dictionary cif-name)))))))

      (if (defined? 'coot-main-menubar)
	  (let ((menu (coot-menubar-menu "Ligand")))

	    (add-simple-coot-menu-menuitem
	     menu "Find Ligands..."
	     do-find-ligands-dialog)

	    (add-simple-coot-menu-menuitem 
	     menu "Jiggle-Fit Ligand"
	     (lambda ()
	       (using-active-atom
		(fit-to-map-by-random-jiggle
		 aa-imol aa-chain-id aa-res-no aa-ins-code 200 1.5))))

	    (add-simple-coot-menu-menuitem
	     menu "Hydrogenate region"
	     (lambda ()
	       (hydrogenate-region 6)))

            (add-simple-coot-menu-menuitem
             menu "Contact Dots for Ligand"
             (lambda ()
               (using-active-atom
                (coot-contact-dots-for-ligand-scm aa-imol (list aa-chain-id aa-res-no aa-ins-code)))))

	    (add-simple-coot-menu-menuitem
	     menu "SMILES → 2D"
	     (lambda ()
	       (generic-single-entry "SMILES string" "" " Send to 2D Viewer "
				     (lambda (text)
				       (smiles-to-ligand-builder text)))))

	    (add-simple-coot-menu-menuitem
	     menu "SMILES → Simple 3D"
	     (lambda ()
	       (generic-double-entry "Residue Name" "SMILES string  " "LIG" "" #f #f "Import Molecule"
				     (lambda (text-1 text-2 dum)
				       (import-rdkit-mol-from-smiles text-2 text-1)))))

	    (add-simple-coot-menu-menuitem 
	     menu "Residue → 2D"
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
	     menu "Solid Generic Objects"
	     (lambda ()
	       (set-display-generic-objects-as-solid 1)
	       (graphics-draw)))

	    (add-simple-coot-menu-menuitem 
	     menu "Unsolid Generic Objects"
	     (lambda ()
	       (set-display-generic-objects-as-solid 0)
	       (graphics-draw)))

	    (add-simple-coot-menu-menuitem
	     menu "Show Chemical Features"
	     (lambda () 
	       (set-display-generic-objects-as-solid 1) ;; there may be consequences...
	       (using-active-atom
		(show-feats aa-imol aa-chain-id aa-res-no aa-ins-code))))

	    (add-simple-coot-menu-menuitem 
	     menu "Tabulate (on terminal) Ligand Distortions"
	     (lambda ()
	       (using-active-atom
		(print-residue-distortions aa-imol aa-chain-id aa-res-no aa-ins-code))))

	    (add-simple-coot-menu-menuitem 
	     menu "Display Ligand Distortions"
	     (lambda ()
	       (set-display-generic-objects-as-solid 1) ;; there may be consequences...
	       (using-active-atom
		(display-residue-distortions aa-imol aa-chain-id aa-res-no aa-ins-code))))

; 	    (add-simple-coot-menu-menuitem
; 	     menu "write sdf file" 
; 	     (lambda ()
; 	       (using-active-atom
; 		(let* ((rn (residue-name aa-imol aa-chain-id aa-res-no aa-ins-code))
; 		       (file-name (string-append rn ".sdf")))
; 		  (residue-to-sdf-file aa-imol aa-chain-id aa-res-no aa-ins-code file-name)))))

;; this is not interesting for the normal user
; 	    (add-simple-coot-menu-menuitem 
; 	     menu "Density Score Ligand"
; 	     (lambda ()
; 	       (using-active-atom
; 		(let ((spec (list aa-chain-id aa-res-no aa-ins-code)))
; 		  (let ((r (density-score-residue aa-imol spec (imol-refinement-map))))
; 		    (format #t "density at ligand atoms: ~s~%" r))))))


; 	    (add-simple-coot-menu-menuitem
; 	     menu "### [Fetch ligand description & generate restraints]"
; 	     (lambda ()
; 	       (using-active-atom
; 		(let ((comp-id (residue-name aa-imol aa-chain-id aa-res-no aa-ins-code)))
; 		  (format #t "here with residue name ~s~%" comp-id)
; 		  (let ((s (get-SMILES-for-comp-id-from-pdbe comp-id)))
; 		    (if (string? s)
; 			(let* ((pdbe-cif-file-name (append-dir-file "coot-download"
; 								    (string-append "PDBe-" comp-id ".cif"))))
; 			  (import-from-3d-generator-from-mdl pdbe-cif-file-name comp-id))))))))


;; this function is in contact-score-isolated-ligand.scm now.
;;
; 	    (add-simple-coot-menu-menuitem 
; 	     menu "Isolate and probe ligand"
; 	     (lambda ()
	       
; 	       (using-active-atom
; 		(let* ((ss (string-append "//" aa-chain-id "/" (number->string aa-res-no)))
; 		       (imol-selection (new-molecule-by-atom-selection aa-imol ss))
; 		       (work-dir (get-directory "coot-molprobity"))
; 		       (tmp-selected-ligand-for-probe-pdb 
; 			(append-dir-file work-dir "tmp-selected-ligand-for-probe.pdb"))
; 		       (tmp-protein-for-probe-pdb
; 			(append-dir-file work-dir "tmp-protein-for-probe.pdb"))
; 		       (probe-dots-file-name
; 			(append-dir-file work-dir "probe.dots")))

; 		  (set-mol-displayed imol-selection 0)
; 		  (set-mol-active    imol-selection 0)
; 		  (write-pdb-file imol-selection tmp-selected-ligand-for-probe-pdb)
; 		  (write-pdb-file aa-imol tmp-protein-for-probe-pdb)
; 		  (goosh-command 
; 		   *probe-command* 
; 		   (list "-u" "-once" (number->string aa-res-no)  ;; -once or -both
; 			 (string-append "not " (number->string aa-res-no))
; 			 "-density60"
; 			 tmp-selected-ligand-for-probe-pdb
; 			 tmp-protein-for-probe-pdb)
; 		   '()
; 		   probe-dots-file-name
; 		   #f)
; 		  (handle-read-draw-probe-dots-unformatted probe-dots-file-name aa-imol 0)
; 		  (graphics-draw)))))


	    ))))


