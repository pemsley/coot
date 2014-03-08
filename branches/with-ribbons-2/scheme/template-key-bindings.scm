
(add-key-binding "Refine Active Residue" "r" (lambda () (manual-refine-residues 0)))
(add-key-binding "Refine Active Residue AA" "x" (lambda () (refine-active-residue)))
(add-key-binding "Triple Refine" "t" (lambda () (manual-refine-residues 1)))
(add-key-binding "Autofit Rotamer" "j" (lambda () (auto-fit-rotamer-active-residue)))
(add-key-binding "Pepflip" "q" (lambda () (pepflip-active-residue)))
(add-key-binding "Go To Blob" "g" (lambda () (blob-under-pointer-to-screen-centre)))
(add-key-binding "Eigen-flip Ligand" "e" (lambda() (flip-active-ligand)))
(add-key-binding "Add Water" "w" (lambda () (place-typed-atom-at-pointer "Water")))
(add-key-binding "Add Water +"  "W"
   (lambda ()
     (blob-under-pointer-to-screen-centre)
     (place-typed-atom-at-pointer "Water")
     (refine-active-residue)))

(add-key-binding "Add terminal residue" "y" 
   (lambda ()
     (using-active-atom
      (add-terminal-residue aa-imol aa-chain-id aa-res-no "auto" 1))))

(add-key-binding "Fill Partial" "k" 
   (lambda ()
     (using-active-atom
      (fill-partial-residue aa-imol aa-chain-id aa-res-no aa-ins-code))))

(add-key-binding "Rotamers dialog for Active Residue" "Q"
   (lambda () 
     (using-active-atom
      (show-rotamers-dialog aa-imol aa-chain-id aa-res-no aa-ins-code aa-alt-conf))))

(add-key-binding "Rotamer name in Status Bar" "~" 
   (lambda ()
     (using-active-atom
      (let ((name (get-rotamer-name aa-imol aa-chain-id aa-res-no aa-ins-code)))
	(if (not name)
	    (add-status-bar-text "No name found")
	    (if (string=? "" name)
		(add-status-bar-text "No name for this")
		(add-status-bar-text (string-append "Rotamer name: " name))))))))

(add-key-binding "Sphere Refinement" "R"
   (lambda ()

     (if (not (valid-map-molecule? (imol-refinement-map)))
	 (info-dialog "Must set the refinement map"))

     (using-active-atom
      
      (let* ((rc-spec (list aa-chain-id aa-res-no aa-ins-code))
	     (ls (residues-near-residue aa-imol rc-spec 4.5)))
	(refine-residues aa-imol (cons rc-spec ls))))))

(add-key-binding "Neighbours Refine" "h"
   (lambda ()

     (if (not (valid-map-molecule? (imol-refinement-map)))
         (info-dialog "Must set the refinement map"))

     (using-active-atom
      
      (let* ((rc-spec (list aa-chain-id aa-res-no aa-ins-code))
             (ls (residues-near-residue aa-imol rc-spec 1.9)))
        (with-auto-accept
         (refine-residues aa-imol (cons rc-spec ls)))))))

(add-key-binding "Sphere Regularization" "B"
   (lambda ()

     (using-active-atom
      
      (let* ((rc-spec (list aa-chain-id aa-res-no aa-ins-code))
	     (ls (residues-near-residue aa-imol rc-spec 3.0)))
	(regularize-residues aa-imol (cons rc-spec ls))))))
	
(add-key-binding "Edit Chi Angles" "X"
   (lambda ()
     (using-active-atom
      (edit-chi-angles aa-imol aa-chain-id aa-res-no aa-ins-code aa-alt-conf))))

(add-key-binding "Jiggle Fit" "J" (lambda () 
				    (using-active-atom 
				     (fit-to-map-by-random-jiggle 
				      aa-imol aa-chain-id aa-res-no aa-ins-code 100 1.0))))

(add-key-binding "Step scrollable map number" "M" 
   (lambda()
     (let ((maps    (map-molecule-list))
	   (current (scroll-wheel-map)))
       (if (not (null? maps))
	   (let ((l (memq current maps)))
	     (set-scroll-wheel-map
	     (if (> (length l) 1) 
		 (cadr l)
		 (car maps))))))))

(add-key-binding "Delete Residue Hydrogens" "P"
   (lambda ()
     (using-active-atom (delete-residue-hydrogens aa-imol aa-chain-id aa-res-no aa-ins-code aa-alt-conf))))


(add-key-binding "ball-and-stickify residue" "$"
   (lambda ()
     (using-active-atom
      (additional-representation-by-attributes aa-imol aa-chain-id aa-res-no aa-res-no aa-ins-code
					       2 2 0.15 1))))

(add-key-binding "Undo Symmetry View" "V" undo-symmetry-view)

(add-key-binding "accept baton position" "A" accept-baton-position)


;; not sure about this one :-)
;; (add-key-binding "Delete this water" "D" (lambda () (apply delete-atom (active-residue))))

