
(define rnase-pdb "data/rnase/tutorial-modern.pdb")
(define rnase-mtz "data/rnase/rnasa-1.8-all_refmac1.mtz")
(define terminal-residue-test-pdb "data/tutorial-add-terminal-0-test.pdb")
(define base-imol (graphics-n-molecules))

(define have-ccp4? #f)
(define imol-rnase -1) 
(define imol-ligand -1) 
(define imol-terminal-residue-test -1)

;; CCP4 is set up? If so, set have-ccp4? #t

(let ((ccp4-master (getenv "CCP4_MASTER")))
  (if (string? ccp4-master)
      (set! have-ccp4? #t)))

(greg-testcase "Close bad molecule" #t
   (lambda ()
     (close-molecule -2)
     #t))

(greg-testcase "Read coordinates test" #t 
   (lambda ()
     (let ((imol (read-pdb rnase-pdb)))
       (set! imol-rnase imol)
       (valid-model-molecule? imol))))

(greg-testcase "Read MTZ test" #t 
   (lambda ()
     (let ((imol-map (make-and-draw-map 
		      rnase-mtz "FWT" "PHWT" "" 0 0)))
       (change-contour-level 0)
       (change-contour-level 0)
       (change-contour-level 0)
       (set-imol-refinement-map imol-map)
       (valid-map-molecule? imol-map))))

(greg-testcase "Another Level test" #t 
   (lambda ()
     (let ((imol-map-2 (another-level)))
       (valid-map-molecule? imol-map-2))))

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
	   #t))))

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

(greg-testcase "Add Terminal Residue Test" #t 
   (lambda ()

     (if (not (string? terminal-residue-test-pdb))
	 (begin 
	   (format #t "~s does not exist - skipping test~%" 
		   terminal-residue-test-pdb)
	   (throw 'untested))
	 (begin
	   (if (file-exists? terminal-residue-test-pdb)
	       (let ((imol (read-pdb terminal-residue-test-pdb)))
		 (if (not (valid-model-molecule? imol))
		     (begin
		       (format #t "~s bad pdb read - skipping test~%" 
			       terminal-residue-test-pdb)
		       (throw 'untested))
		     (begin
		       (add-terminal-residue imol "A" 1 "ALA" 1)
		       (write-pdb-file imol "regression-test-terminal-residue.pdb")
		       ;; where did that extra residue go?
		       ;; 
		       ;; Let's copy a fragment from imol and go to
		       ;; the centre of that fragment, and check where
		       ;; we are and compare it to where we expected
		       ;; to be.
		       ;; 
		       (let ((new-mol (new-molecule-by-atom-selection
				       imol "//A/0")))
			 (if (not (valid-model-molecule? new-mol))
			     (throw 'fail)
			     (move-molecule-here new-mol)
			     (let ((rc (rotation-centre)))
			       (let ((r (apply + (map (lambda (rc1 x1)
							(- rc1 x1))
						      rc (list 45.6 15.8 11.8)))))
			       (if (> r 0.66)
				   (begin 
				     (format #t "Bad placement of terminal residue~%")
				     #f)
				   #t)))))))))))))

				      
				      

	       
