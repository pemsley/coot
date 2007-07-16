

(greg-testcase "Get monomer test" #t 
   (lambda ()
     (if have-ccp4?
	 (let ((imol (monomer-molecule-from-3-let-code "3GP" "")))
	   (if (valid-model-molecule? imol)
	       (begin
		 (set! imol-ligand imol) ; for use in next test
		 (delete-residue-hydrogens imol "A" 1 "" "")
		 #t)))
	 (begin
	   (format #t "CCP4 not set up - skipping 3GP test~%")
	   (throw 'untested)))))

(greg-testcase "Set Bond thickness" #t 
   (lambda ()	       
     (if (valid-model-molecule? imol-ligand)
	 (begin
	   (set-bond-thickness imol-ligand 5)
	   #t)
	 (begin
	   (format #t "No ligand molecule - Skipping bond thickness test~%")
	   (throw 'untested)))))
	 

(greg-testcase "Move and Refine Ligand test" #t 
   (lambda ()
     (let ((new-rc (list 55.3 9.1 20.6)))
       (if (not (valid-model-molecule? imol-ligand))
	   (throw 'untested)
	   (begin
	     (apply spin-zoom-trans (append (list 2 200 0.2 0.2)
					    (map (lambda (rc nc) (- nc rc))
						 (rotation-centre) new-rc)))
	     ;; updates the map:
	     (apply set-rotation-centre new-rc)
	     (move-molecule-here imol-ligand)
	     (let ((backup-mode (backup-state imol-ligand))
		   (alt-conf "")
		   (replacement-state (refinement-immediate-replacement-state)))
	       
	       (turn-off-backup imol-ligand)
	       (set-refinement-immediate-replacement 1)
	       (refine-zone imol-ligand "A" 1 1 alt-conf)
	       (accept-regularizement)
	       (rotate-y-scene 1000 0.1)
	       (if (= replacement-state 0)
		   (set-refinement-immediate-replacement 0))
	       (if (= backup-mode 1)
		   (turn-on-backup imol-ligand))
	       #t ; pah
	       ))))))
