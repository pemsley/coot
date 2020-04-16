
;;;; Copyright 2007, 2008 by The University of Oxford
;;;; Copyright 2014, 2015 by Medical Research Council
;;;; Author: Paul Emsley

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
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;;; 02110-1301, USA

(greg-testcase "Get monomer test" #t 
   (lambda ()

     (let ((imol (get-monomer "3GP")))
       (if (valid-model-molecule? imol)
	   (begin
	     (set! imol-ligand imol) ; for use in next test
	     (delete-residue-hydrogens imol "A" 1 "" "")
	     #t)
	   (begin
	     (format #t "   No ligand molecule - monomer test~%")
	     (throw 'untested)
	     #f)))))


(greg-testcase "Delete all-molecule Hydrogens" #t
   (lambda ()

     (let ((imol (greg-pdb "monomer-3GP.pdb")))
       (if (not (valid-model-molecule? imol))
	   (throw 'fail))
       (format #t "here 1 ~%")
       (let ((n (delete-hydrogens imol)))
	 (> n 0)))))



(greg-testcase "Non-Autoloads" #t
   (lambda ()

     (define (get-ccp4-version)
       (let ((s (string-split (shell-command-to-string "cad -i") #\newline)))
	 (let loop ((s s))
	   (cond 
	    ((null? s) #f)
	    ((string-match "CCP4 software suite: patch level" (car s))
	     (let ((sl (string->list-of-strings (car s))))
	       (car (reverse sl))))
	    (else (loop (cdr s)))))))

     (define (old-ccp4-restraints?)
       (if (not have-ccp4?)
           #f
	   (let ((v (get-ccp4-version)))
	     (if (not (string? v))
		 #t ;; bizarre
		 (string<? v "6.2")))))

     (let ((r-1 (monomer-restraints "LIG"))
	   (o (old-ccp4-restraints?)))
       (greg-pdb "test-LIG.pdb")
       (let ((r-2 (monomer-restraints "LIG")))
	 (remove-non-auto-load-residue-name "LIG")
	 (greg-pdb "test-LIG.pdb")
	 (let ((r-3 (monomer-restraints "LIG")))
	   (delete-restraints "LIG")
	   (add-non-auto-load-residue-name "LIG")
	   (let ((r-4 (monomer-restraints "LIG")))

	     ;; r-1, r-2, r-4 should be #f, r-3 should be filled.
	     (if (all-true? 
		  (list 
		   (eq? #f r-1)
		   (eq? #f r-2)
		   (eq? #f r-4)
		   ;; just to be clear here: either this is modern
		   ;; restraints and we should get a list in r-3... or
		   ;; this is old restraints. Either of those should
		   ;; result in a PASS.
                   (or 
                    (and (not (old-ccp4-restraints?))
                         (list? r-3))
                    (old-ccp4-restraints?))))
		 #t ;; good return value
		 (begin
		   (format #t "restraints: r-1 ~s~%" r-1)
		   (format #t "restraints: r-2 ~s~%" r-2)
		   (format #t "restraints: r-3 ~s~%" r-3)
		   (format #t "restraints: r-4 ~s~%" r-4)
		   #f))))))))




(greg-testcase "Merge molecules of a ligand with a spec" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb"))
	   (imol-lig (get-monomer "3GP")))

       (let ((spec (list "L" 1 "")))

	 (set-merge-molecules-ligand-spec spec)
	 (merge-molecules (list imol-lig) imol)

	 ;; now check that L1 exists in imol

	 ;; residue-spec->residue-name expects a 4-ele spec. Hmm.	 
	 (let ((rn (residue-spec->residue-name imol (cons #t spec))))
	   (if (not (string? rn))
	       #f
	       (string=? rn "3GP")))))))



(greg-testcase "Move and Refine Ligand test" #t 
   (lambda ()

     (let ((imol (greg-pdb "monomer-3GP.pdb")))
       (if (not (valid-model-molecule? imol))
	   #f
	   (let ((new-rc (list 55.3 9.1 20.6)))
	     ;; set the view
	     (let ((view-number (add-view (list    54.5698 8.7148 20.5308)
					  (list 0.046229 -0.157139 -0.805581 0.569395)
					  19.8858
					  "ligand-view")))
	       (go-to-view-number view-number 1))
				
	     ;; updates the map:
	     (apply set-rotation-centre new-rc)
	     (move-molecule-here imol)
	     (let ((backup-mode (backup-state imol))
		   (alt-conf "")
		   (replacement-state (refinement-immediate-replacement-state)))
	       
	       (turn-off-backup imol)
	       (set-refinement-immediate-replacement 1)
	       (refine-zone imol "A" 1 1 alt-conf)
	       (accept-regularizement)
	       (rotate-y-scene (rotate-n-frames 600) 0.1)
	       (if (= replacement-state 0)
		   (set-refinement-immediate-replacement 0))
	       (if (= backup-mode 1)
		   (turn-on-backup imol))
	       #t ; success
	       ))))))

(greg-testcase "Many Molecules - Ligand Fitting" #t 
   (lambda ()

     (let* ((npo-pdb  (append-dir-file greg-data-dir "monomer-NPO.pdb"))
	    (43ca-pdb (append-dir-file greg-data-dir "pdb43ca-sans-NPO-refmaced.pdb"))
	    (43ca-mtz (append-dir-file greg-data-dir "pdb43ca-sans-NPO-refmaced.mtz"))
	    (imol-npo (handle-read-draw-molecule-with-recentre npo-pdb 0)))
       
       (if (not (valid-model-molecule? imol-npo))
	   #f
	   (begin 
	     (let loop ((count 0))
	       (cond
		((= count 5) 'done)
		(else 
		 (let ((imol-copy (copy-molecule imol-npo)))
		   (set-mol-displayed imol-copy 0))
		 (loop (+ count 1)))))

	     (let* ((imol-protein (read-pdb 43ca-pdb))
		    (imol-maps (auto-read-make-and-draw-maps 43ca-mtz))
		    (imol-map-1 (vector-ref imol-maps 0))
		    (imol-map-2 (vector-ref imol-maps 1)))
	     
	       (add-ligand-clear-ligands)
	       (set-ligand-search-protein-molecule imol-protein)
	       (set-ligand-search-map-molecule imol-map-1)
	       (add-ligand-search-ligand-molecule imol-npo)

	       (let ((solutions (execute-ligand-search))) ; crash
		 (format #t "   Fitting NPO gave these results: ~s~%" solutions)
		 (set-map-displayed imol-map-1 0)
		 (set-map-displayed imol-map-2 0)
		 (set-mol-displayed imol-protein 0)
		 (set-mol-displayed imol-npo 0)
		 #t)))))))




(greg-testcase "flip residue (around eigen vectors)" #t 
   (lambda ()

     (let* ((imol-orig (greg-pdb "monomer-3GP.pdb"))
	    (imol-copy (copy-molecule imol-orig)))
       (if (not (valid-model-molecule? imol-copy))
	   #f

	   (begin

	     ;; we need this, otherwise active-atom is (accidentally) the wrong molecule
	     (set-go-to-atom-molecule imol-copy)
	     (set-go-to-atom-chain-residue-atom-name "A" 1 " C8 ")

	     (let ((active-atom (active-residue)))
	       (if (not active-atom)
		   (begin
		     (format #t "No active atom~%")
		     #f)
		   (let ((imol      (list-ref active-atom 0))
			 (chain-id  (list-ref active-atom 1))
			 (res-no    (list-ref active-atom 2))
			 (ins-code  (list-ref active-atom 3))
			 (atom-name (list-ref active-atom 4))
			 (alt-conf  (list-ref active-atom 5)))
		     (if (= imol imol-orig)
			 (begin
			   (format #t "oops - didn't pick the copy for active res~%")
			   #f)
			 (begin
			   (flip-ligand imol chain-id res-no)
			   (let ((atom-orig-1 (get-atom imol-orig "A" 1 "" " C8 "))
				 (atom-move-1 (get-atom imol      "A" 1 "" " C8 ")))

			     (if (not (list? atom-orig-1))
				 (begin
				   (format #t "atom-orig-1 not found~%")
				   (throw 'fail)))
			     
			     (if (not (list? atom-move-1))
				 (begin
				   (format #t "atom-move-1 not found~%")
				   (throw 'fail)))

			     (let ((d (bond-length (list-ref atom-orig-1 2)
						   (list-ref atom-move-1 2))))
			       (format #t "distance: ~s~%" d)
			       (if (not (> d 2.1))
				   (begin
				     (format #t "fail to move test atom d1~%"))
				   (begin
				     (flip-ligand imol chain-id res-no)
				     (flip-ligand imol chain-id res-no)
				     (flip-ligand imol chain-id res-no)
				     ;; having flipped it round the axes 4
				     ;; times, we should be back where we
				     ;; started.
				     (let ((atom-orig-1 (get-atom imol-orig "A" 1 "" " C8 "))
					   (atom-move-1 (get-atom imol      "A" 1 "" " C8 ")))
				       (let ((d2 (bond-length (list-ref atom-orig-1 2)
							      (list-ref atom-move-1 2))))
					 (format #t "distance d2: ~s~%" d2)
					 (if (not (< d2 0.001))
					     (begin
					       (format #t "fail to move atom back to start d2~%"))
					     #t)))))))))))))))))

;; no partial charges in the dictionary now.
;;
;(greg-testcase "Test dipole" #t
;   (lambda ()

;     (let ((imol (greg-pdb "dipole-residues.pdb")))

;       (if (not (valid-model-molecule? imol))
;           (begin
;             (format #t "dipole-residues.pdb not found~%")
;             #f)

;           (let* ((residue-specs
;                   (list
;                    (list "A" 1 "")
;                    (list "A" 2 "")
;                    (list "A" 3 "")))
;                  (dipole (add-dipole-for-residues imol residue-specs)))

;             (if (not dipole)
;                 (begin
;                   (format #t "bad dipole ~s~%" dipole)
;                   #f)
;                 (let ((d (car dipole))
;                       (dip (cadr dipole)))

;                   (let ((dip-x (list-ref dip 0))
;                         (dip-y (list-ref dip 1))
;                         (dip-z (list-ref dip 2)))

;                     (format #t "info:: dipole components ~s ~%" dip)

;                     (if (not (and (close-float? dip-y 0)
;                                   (close-float? dip-z 0)))
;                         (begin
;                           (format #t "bad dipole y z components ~s ~s~%"
;                                   dip-y dip-z)
;                           #f)
			 
;			 ;; dipole points in the negative x direction
;                         (if (and (< dip-x 0) 
;                                  (> dip-x -20))

;                             #t
;                             #f))))))))))


(greg-testcase "Reading new dictionary restraints replaces" ;; not adds
   #t
   (lambda ()

     (define (get-torsions r)
       ; (format #t "r: ~s~%" r)
       (cdr (assoc "_chem_comp_tor" r)))

     (read-cif-dictionary (append-dir-file greg-data-dir "libcheck_3GP.cif"))
     (read-cif-dictionary (append-dir-file greg-data-dir "libcheck_3GP.cif"))
     (read-cif-dictionary (append-dir-file greg-data-dir "libcheck_3GP.cif"))

     (let* ((r (monomer-restraints "3GP"))
	    (t (get-torsions r)))

       ;; 
       (if (not (< (length t) 26)) ;; 22 in new dictionary, it seems.
	   (begin
	     (format #t "torsions: ~s ~s~%" (length t) t)
	     (throw 'fail))
	   #t))))


(if #f
(greg-testcase "Pyrogen Runs OK?"  #t
   (lambda ()

     (if (not (enhanced-ligand-coot?))
	 #t ;; fake a pass

	 ;; bad things may well happen if we run the wrong version of pyrogen.
	 ;; so force pyrogen to be the one that is installed alongside this version of coot 
	 ;; that we are running. We do that by manipulation of pkglibdir (which points to xxx/lib/guile)
	 ;;
	 (let* ((p (cdr (assoc 'pkglibdir %guile-build-info)))
		(lib-dir (file-name-directory p))
		(prefix-dir (append-dir-dir lib-dir ".."))
		(bin-dir (append-dir-dir prefix-dir "bin"))
		(pyrogen-bin (append-dir-file bin-dir "pyrogen")))

	   (let ((smiles "C1CNC1")
		 (tlc-text "XXX")
		 (log-file-name "pyrogen.log"))

             ;; no need to do this test twice.
	     (if (not (enhanced-ligand-coot?))

		 #t ;; don't test pyrogen

		 (let ((goosh-status
			(goosh-command
			 pyrogen-bin
			 (if *use-mogul*
			     (list "--residue-type" tlc-text smiles)
			     (if (command-in-path? "mogul")
				 (list                   "--residue-type" tlc-text smiles)
				 (list "--no-mogul" "-M" "--residue-type" tlc-text smiles)))
			 '() log-file-name #t)))
		   (if (not (ok-goosh-status? goosh-status))
		       (begin
			 (format #t "WARNING:: OOps - pyrogen exited with status: ~s~%~!" goosh-status)
			 #f)
		       (begin
			 (let* ((pdb-file-name (string-append tlc-text "-pyrogen.pdb"))
				(cif-file-name (string-append tlc-text "-pyrogen.cif"))
				(imol (handle-read-draw-molecule-with-recentre pdb-file-name 0)))

			   (format #t "INFO:: pyrogen test will try to read ~s~%~!" pdb-file-name)
			   ;; add test for chirality in the dictionary here 
			   (valid-model-molecule? imol)))))))))))
)

(if #f
(greg-testcase "pyrogen dictionary does not make double-quoted atom names" #t
   (lambda ()

     ;; make sure that you are running the correct pyrogen

     (if (not (enhanced-ligand-coot?))
	 #t ;; fake a pass
	 (begin
	   (if (file-exists? "UVP-pyrogen.cif")
	       (delete-file  "UVP-pyrogen.cif"))

	   (let ((goosh-status (goosh-command
				"pyrogen"
				(list "-nM" "-r" "UVP"
				      "CO[C@@H]1[C@H](O)[C@H](O[C@H]1[n+]1ccc(O)nc1O)\\C=C\\P(O)(O)=O")
				'()
				"pyrogen.log" #f)))
	     (if (not (ok-goosh-status? goosh-status))
		 (begin
		   (format #t "Fail to correctly run pyrogen\n")
		   #f)
		 (begin
		   (read-cif-dictionary "UVP-pyrogen.cif")
		   (let ((imol (get-monomer "UVP")))
		     (if (not (valid-model-molecule? imol))
			 (begin
			   (format #t "Fail to load molecule from pyrogen dictionary\n"))
			 (let ((atoms-info (residue-info imol "A" 1 ""))
			       (passes #t))
			   (for-each (lambda(atom)
				       (let ((atom-name (residue-atom->atom-name atom)))
					 ;; (format #t "Atom name ~s~%" atom-name)
					 (if (string-match "\"" atom-name)
					     (begin
					       (format #t "Atom name quote fail ~s~%" atom-name)
					       (set! passes #f)))))
				     atoms-info)
			   passes))))))))))
)

;; FLEV will not make a PNG if it is not compiled with 
;; C++-11 - and that is OK for 0.8.9.x.
;;
; (greg-testcase "FLEV makes a PNG" #t
;    (lambda ()

;      (let ((fn "test-flev-greg-testcase.png"))

;        (if (file-exists? fn)
; 	   (delete-file fn))

;        (let ((imol (greg-pdb "tutorial-modern.pdb"))
; 	     (imol-ligand (get-monomer "3GP")))
	 
; 	 (set-rotation-centre 54 10 20)
; 	 (move-molecule-to-screen-centre imol-ligand)
; 	 (set-merge-molecules-ligand-spec (list "L" 1 ""))
; 	 (merge-molecules (list imol-ligand) imol)
; 	 (fle-view-with-rdkit-to-png imol "L" 1 "" 4.8 fn)
; 	 (file-exists? fn)))))
