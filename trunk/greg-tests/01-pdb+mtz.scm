
(define rnase-pdb "greg-data/tutorial-modern.pdb")
(define rnase-mtz "greg-data/rnasa-1.8-all_refmac1.mtz")
(define terminal-residue-test-pdb "greg-data/tutorial-add-terminal-0-test.pdb")
(define base-imol (graphics-n-molecules))

(define have-ccp4? #f)
(define imol-rnase -1) 
(define imol-ligand -1) 
(define imol-terminal-residue-test -1)

;; CCP4 is set up? If so, set have-ccp4? #t


;(greg-gui 'register "Close bad molecule")
;(greg-gui 'register "Read coordinates test")
;(greg-gui 'register "Read MTZ test")
;(greg-gui 'register "Another Level test")
;(greg-gui 'register "Add Terminal Residue Test")
;(greg-gui 'register "Select by Sphere")

;(greg-gui 'make-gui)


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

(greg-testcase "Another Level Test" #t 
   (lambda ()
     (let ((imol-map-2 (another-level)))
       (valid-map-molecule? imol-map-2))))

(greg-testcase "Set Atom Attribute Test" #t
	       (lambda ()
		 (set-atom-attribute imol-rnase "A" 11 "" " CA " "" "x" 64.5) ; an Angstrom or so
		 (let ((atom-ls (residue-info imol-rnase "A" 11 "")))
		   (let f ((atom-ls atom-ls))
		     (cond
		      ((null? atom-ls) (throw 'fail))
		      (else 
		       (let* ((atom (car atom-ls))
			      (compound-name (car atom))
			      (atom-name (car compound-name)))
			 (if (string=? atom-name " CA ")
			     (let* ((xyz (car (cdr (cdr atom))))
				    (x (car xyz)))
			       (close-float? x 64.5))
			     (f (cdr atom-ls))))))))))


(greg-testcase "Add Terminal Residue Test" #t 
   (lambda ()

     (if (not (string? terminal-residue-test-pdb))
	 (begin 
	   (format #t "~s does not exist - skipping test~%" 
		   terminal-residue-test-pdb)
	   (throw 'untested))
	 (begin
	   (if (not (file-exists? terminal-residue-test-pdb))
	       (begin 
		 (format #t "~s does not exist - skipping test~%" 
			 terminal-residue-test-pdb)
		 (throw 'untested))
	       
	       ;; OK, file exists
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
			     (begin
			       (move-molecule-here new-mol)
			       (let ((rc (rotation-centre)))
				 (let ((r (apply + (map (lambda (rc1 x1)
							  (- rc1 x1))
							rc (list 45.6 15.8 11.8)))))
				   (if (> r 0.66)
				       (begin 
					 (format #t "Bad placement of terminal residue~%")
					 #f)
				       #t))))))))))))))



(greg-testcase "Select by Sphere" #t
   (lambda ()

     (let ((imol-sphere (new-molecule-by-sphere-selection 
			 imol-rnase 
			 24.6114959716797 24.8355808258057 7.43978214263916
			 3.6)))

       (if (not (valid-model-molecule? imol-sphere))
	   (begin
	     (format #t "Bad sphere molecule~%")
	     #f)
	   
	   (let ((n-atoms 
		  (apply + (map (lambda (chain-id)
				  ;; return the number of atoms in this chain
				  (let ((n-residues (chain-n-residues chain-id imol-sphere)))
				    (format #t "Sphere mol: there are ~s residues in chain ~s~%" 
					    n-residues chain-id)

				    (let loop ((residue-list (number-list 0 (- n-residues 1)))
					       (chain-n-atoms 0))

				      (cond 
				       ((null? residue-list) 
					(format #t "chain-n-atoms for chain ~s is ~s~%" chain-id chain-n-atoms)
					chain-n-atoms)
				       (else 
					(let ((serial-number (car residue-list)))
					  (let ((res-name (resname-from-serial-number imol-sphere chain-id serial-number))
						(res-no   (seqnum-from-serial-number  imol-sphere chain-id serial-number))
						(ins-code (insertion-code-from-serial-number imol-sphere chain-id serial-number)))
					    (let ((residue-atoms-info (residue-info imol-sphere chain-id res-no ins-code)))
					      (loop (cdr residue-list) (+ (length residue-atoms-info) chain-n-atoms))))))))))
				(chain-ids imol-sphere)))))
	     
	     (format #t "Found ~s sphere atoms ~%" n-atoms)

	     (= n-atoms 20))))))

(greg-testcase "Test Views" #t 
	       (lambda ()
		 (let ((view-number 
			(add-view (list    32.0488 21.6531 13.7343)
				  (list -0.12784 -0.491866 -0.702983 -0.497535)
				  20.3661
				  "B11 View")))
		   (go-to-view-number view-number 1)
		   #t)))


(greg-testcase "Label Atoms and Delete" #t 
   (lambda ()

     (let ((imol-frag (new-molecule-by-atom-selection imol-rnase "//B/10-12")))
       (set-rotation-centre 31.464  21.413  14.824)
       (map (lambda (n)
	      (label-all-atoms-in-residue imol-frag "B" n ""))
	    '(10 11 12))
       (rotate-y-scene 100 0.1)
       (delete-residue imol-frag "B" 10 "")
       (delete-residue imol-frag "B" 11 "")
       (delete-residue imol-frag "B" 12 "")
       (rotate-y-scene 100 0.1)
       #t))) ; it didn't crash - good :)


(greg-testcase "Rotamer outliers" #t

   (lambda ()

     ;; pre-setup so that residues 1 and 2 are not rotamers "t" but 3
     ;; to 8 are (1 and 2 are more than 40 degrees away). 9-16 are
     ;; other residues and rotamers.
     (let ((imol-rotamers (read-pdb (append-dir-file greg-data-dir 
						     "rotamer-test-fragment.pdb"))))

       (let ((rotamer-anal (rotamer-graphs imol-rotamers)))
	 (if (not (= (length rotamer-anal) 14))
	     (throw 'fail)
	     (let ((a-1 (list-ref rotamer-anal 0))
		   (a-2 (list-ref rotamer-anal 1))
		   (a-last (car (reverse rotamer-anal))))

	       (let ((anal-str-a1 (car (reverse a-1)))
		     (anal-str-a2 (car (reverse a-2)))
		     (anal-str-a3 (car (reverse a-last))))

		 (if (and (string=? anal-str-a1 "Rotamer not recognised")
			  (string=? anal-str-a2 "Rotamer not recognised")
			  (string=? anal-str-a3 "Missing Atoms"))
		     #t 
		     (begin 
		       (format #t "  failure rotamer test: ~s ~s ~s~%"
			       a-1 a-2 a-last)
		       (throw 'fail))))))))))


;; Don't reset the occupancies of the other parts of the residue
;; on regularization using alt confs
;; 
(greg-testcase "Alt Conf Occ Sum Reset" #t

  (lambda () 

     (define (get-occ-sum imol-frag)

       (define (occ att-atom-name att-alt-conf atom-ls)

	 (let ((atom-ls atom-ls)
	       (occ-sum 0))

	   (let f ((atom-ls atom-ls))
	     (cond 
	      ((null? atom-ls) (throw 'fail))
	      (else 
	       (let* ((atom (car atom-ls))
		      (compound-name (car atom))
		      (atom-name (car compound-name))
		      (alt-conf (car (cdr compound-name)))
		      (occ (car (car (cdr atom)))))
		 (if (string=? atom-name att-atom-name)
		     (if (string=? alt-conf att-alt-conf)
			 (begin 
			   (format #t "for atom ~s ~s returning occ ~s~%" 
				   att-atom-name att-alt-conf occ)
			   occ)
			 (f (cdr atom-ls)))
		     (f (cdr atom-ls)))))))))

       ;; get-occ-sum body
       (let ((atom-ls (residue-info imol-frag "X" 15 "")))
	 (if (not (list? atom-ls))
	     (throw 'fail)
	     (+ (occ " CE " "A" atom-ls) (occ " CE " "B" atom-ls)
		(occ " NZ " "A" atom-ls) (occ " NZ " "B" atom-ls)))))

	 

     ;; main body 
     ;; 
     (let ((imol-frag (read-pdb (append-dir-file greg-data-dir "res098.pdb"))))

       (let ((occ-sum-pre (get-occ-sum imol-frag)))

	 (let ((replace-state (refinement-immediate-replacement-state)))
	   (set-refinement-immediate-replacement 1)
	   (regularize-zone imol-frag "X" 15 15 "A")
	   (accept-regularizement)
	   (if (= replace-state 0)
	       (set-refinement-immediate-replacement 0))
	  
	   (let ((occ-sum-post (get-occ-sum imol-frag)))
	     
	     (format #t "test for closeness: ~s ~s~%" occ-sum-pre occ-sum-post)	    
	     (close-float? occ-sum-pre occ-sum-post)))))))


;; This test we expect to fail until the CISPEP correction code is in
;; place (using mmdb-1.10+).
;; 
(greg-testcase "Correction of CISPEP test" #f
   (lambda ()	      

;; In this test the cis-pep-12A has indeed a CIS pep and has been
;; refined with refmac and is correctly annotated.  There was 4 on
;; read.  There should be 3 after we fix it and write it out.
;; 
;; Here we sythesise user actions:
;; 1) CIS->TRANS on the residue 
;; 2) convert leading residue to a ALA
;; 3) mutate the leading residue and auto-fit it to correct the chiral
;;    volume error :)

      (let ((chain-id "A")
	    (resno 11)
	    (ins-code ""))

	(let ((cis-pep-mol (read-pdb 
			    (append-dir-file greg-data-dir 
					     "tutorial-modern-cis-pep-12A_refmac0.pdb"))))

	  (let ((view-number (add-view 
			      (list    63.455 11.764 1.268)
			      (list -0.760536 -0.0910907 0.118259 0.631906)
			      15.7374
			      "CIS-TRANS cispep View")))
	    (go-to-view-number view-number 1))

	  (cis-trans-convert cis-pep-mol chain-id resno ins-code)
	  (pepflip resno chain-id cis-pep-mol)
	  (let ((res-type (residue-name cis-pep-mol chain-id resno ins-code)))
	    (mutate resno chain-id cis-pep-mol "GLY")
	    (with-auto-accept
	     (refine-zone cis-pep-mol chain-id resno (+ resno 1) "")
	     (accept-regularizement)
	     (mutate resno chain-id cis-pep-mol res-type)
	     (auto-fit-best-rotamer resno "" ins-code chain-id cis-pep-mol 
				    (imol-refinement-map) 1 1)
	     (refine-zone cis-pep-mol chain-id resno (+ resno 1) "")
	     (accept-regularizement))
	     ;; should be repared now.  Write it out

	    (let ((tmp-file "tmp-fixed-cis.pdb"))
	      (write-pdb-file cis-pep-mol tmp-file)
	      (let ((o (run-command/strings "grep" (list"-c" "CISPEP" tmp-file) '())))
		(if (not (list? o))
		    (throw 'fail)
		    (if (not (= (length o) 1))
			(throw 'fail)
			(let ((parts (split-after-char-last #\: (car o) list)))
			  (format #t "CISPEPs: ~s~%" (car (cdr parts)))
			  (if (not (string=? "3" (car (cdr parts))))
			      (throw 'fail)
			      #t)))))))))))


(greg-testcase "Rigid Body Refine Alt Conf Waters" #t
   (lambda ()

     (let ((imol-alt-conf-waters (read-pdb (append-dir-file greg-data-dir
							    "alt-conf-waters.pdb"))))
       
       (let ((rep-state (refinement-immediate-replacement-state)))
	 (set-refinement-immediate-replacement 1)
	 (refine-zone imol-alt-conf-waters "D" 71 71 "A")
	 (accept-regularizement)
	 (refine-zone imol-alt-conf-waters "D" 2 2 "")
	 (accept-regularizement)
	 (set-refinement-immediate-replacement rep-state)
	 #t))))  ;; good it didn't crash.




