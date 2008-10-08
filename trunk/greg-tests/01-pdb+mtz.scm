;;;; Copyright 2007, 2008 by The University of Oxford
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

;; provides read-line:
(use-modules (ice-9 rdelim))

(define rnase-pdb (append-dir-file greg-data-dir "tutorial-modern.pdb"))
(define rnase-mtz (append-dir-file greg-data-dir "rnasa-1.8-all_refmac1.mtz"))
(define terminal-residue-test-pdb (append-dir-file greg-data-dir "tutorial-add-terminal-1-test.pdb"))
(define base-imol (graphics-n-molecules))
(define rnase-seq     (append-dir-file greg-data-dir "rnase.seq"))


(define have-ccp4? #f)
(define imol-rnase -1) 
(define imol-rnase-map -1) 
(define imol-ligand -1) 
(define imol-terminal-residue-test -1)

(define horne-works-cif (append-dir-file greg-data-dir "lib-B3A.cif"))
(define horne-cif       (append-dir-file greg-data-dir "lib-both.cif"))
(define horne-pdb       (append-dir-file greg-data-dir "coords-B3A.pdb"))
(define ins-code-frag-pdb (append-dir-file greg-data-dir "ins-code-fragment-pre.pdb"))


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
      (begin
	(format "==== It seems CCP4 is setup ===~%")
	(format "==== CCP4_MASTER: ~s~%" ccp4-master)
	(set! have-ccp4? #t))))


(greg-testcase "Post Go To Atom no molecules"  #t
   (lambda () 
     (post-go-to-atom-window) ; crashes pre r820
     #t)) 
     

(greg-testcase "Close bad molecule" #t
   (lambda ()
     (close-molecule -2)
     #t))

(greg-testcase "Read coordinates test" #t 
   (lambda ()
     (let ((imol (read-pdb rnase-pdb)))
       (set! imol-rnase imol)
       (valid-model-molecule? imol))))

(greg-testcase "New molecule from bogus molecule" #t
    (lambda ()
      (let* ((pre-n-molecules (graphics-n-molecules))
	     (new-mol (new-molecule-by-atom-selection -5 "//A/0"))
	     (post-n-molecules (graphics-n-molecules)))
	(if (not (= new-mol -1))
	    (begin
	      (format #t "  fail on non-matching n molecules.")
	      (throw 'fail))
	    (if (not (= pre-n-molecules post-n-molecules))
		(begin 
		  (format #t "  fail on non-matching n molecules.")
		  (throw 'fail))
		#t)))))

(greg-testcase "New molecule from bogus atom selection" #t
    (lambda ()
      (let* ((pre-n-molecules (graphics-n-molecules))
	     (new-mol (new-molecule-by-atom-selection  imol-rnase "//A/100"))
	     (post-n-molecules (graphics-n-molecules)))
	;; what should happen if there are no atoms in the new-mol?
	(format #t "   INFO:: pre-n-molecules ~s   post-n-molecules ~s~%"
		pre-n-molecules post-n-molecules)
	#t)))



(greg-testcase "ins code change and Goto atom over an ins code break" #t
   (lambda () 	       

     (define (matches-attributes atts-obs atts-ref)
       (equal? atts-obs atts-ref))
       

     ;; main line
     (if (not (file-exists? ins-code-frag-pdb))
	 (begin
	   (format #t "  WARNING:: file not found: ~s~%" ins-code-frag-pdb)
	   (throw 'fail)))
     (let ((frag-pdb (handle-read-draw-molecule-with-recentre 
		      ins-code-frag-pdb 0)))
       (set-go-to-atom-molecule frag-pdb)
       (set-go-to-atom-chain-residue-atom-name "A" 68 " CA ")
       (let* ((ar-1 (active-residue))
	      (ins-1 (list-ref ar-1 3)))
	 (change-residue-number frag-pdb "A" 68 "" 68 "A")
	 (change-residue-number frag-pdb "A" 69 "" 68 "B")
	 (change-residue-number frag-pdb "A" 67 "" 68 "")
	 (let* ((ar-2 (active-residue))
		(ins-2 (list-ref ar-2 3)))
	   (format #t "pre and post ins codes: ~s ~s~%" ins-1 ins-2)
	   (if (not (string=? ins-2 "A"))
	       (begin
		 (format #t " Fail ins code set: ~s is not \"A\"~%" ins-2)
		 (throw 'fail)))

	   (let ((test-expected-results
		  (list 
		   (cons (goto-next-atom-maybe "A" 67 ""  " CA ") (list "A" 68 ""  " CA "))
		   (cons (goto-next-atom-maybe "A" 68 "A" " CA ") (list "A" 68 "B" " CA "))
		   (cons (goto-next-atom-maybe "A" 68 "B" " CA ") (list "A" 70 ""  " CA "))
		   (cons (goto-next-atom-maybe "D" 10 ""  " O  ") (list "A" 62  "" " CA "))
		   (cons (goto-prev-atom-maybe "A" 70 ""  " CA ") (list "A" 68 "B" " CA "))
		   (cons (goto-prev-atom-maybe "A" 68 "B" " CA ") (list "A" 68 "A" " CA "))
		   (cons (goto-prev-atom-maybe "A" 68 "A" " CA ") (list "A" 68  "" " CA "))
		   (cons (goto-prev-atom-maybe "A" 68 ""  " CA ") (list "A" 66  "" " CA ")))))

	     (let loop ((real-results (map car test-expected-results))
			(expected-results (map cdr test-expected-results)))

	       (cond 
		((null? real-results) #t) ; pass
		((equal? (car real-results) (car expected-results))
		 (format #t "   pass: ~s ~%" (car expected-results))
		 (loop (cdr real-results) (cdr expected-results)))
		(else 
		 (format #t "   fail: real: ~s expected: ~s ~%" 
			 (car real-results) (car expected-results))
		 #f)))))))))


(greg-testcase "Read a bogus map" #t
    (lambda ()
      (let ((pre-n-molecules (graphics-n-molecules))
	    (imol (handle-read-ccp4-map "bogus.map" 0)))
	(if (not (= -1 imol))
	    (begin
	      (format #t "bogus ccp4 map returns wrong molecule number~%")
	      (throw 'fail))
	    (let ((now-n-molecules (graphics-n-molecules)))
	      (if (not (= now-n-molecules pre-n-molecules))
		  (begin
		    (format #t "bogus ccp4 map creates extra map ~s ~s~%"
			    pre-n-molecules now-n-molecules)
		    (throw 'fail))
		  #t))))))


(greg-testcase "Read MTZ test" #t 
   (lambda ()

     ;; bogus map test
     (let ((pre-n-molecules (graphics-n-molecules))
	   (imol-map (make-and-draw-map "bogus.mtz" "FWT" "PHWT" "" 5 6)))
       (if (not (= -1 imol-map))
	   (begin
	     (format #t "   bogus MTZ returns wrong molecule number~%")
	     (throw 'fail))
	   (let ((now-n-molecules (graphics-n-molecules)))
	     (if (not (= now-n-molecules pre-n-molecules))
		 (begin
		   (format #t "   bogus MTZ creates extra map ~s ~s~%"
			   pre-n-molecules now-n-molecules)
		   (throw 'fail))))))

     ;; correct mtz test
     (let ((imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))
       (change-contour-level 0)
       (change-contour-level 0)
       (change-contour-level 0)
       (set-imol-refinement-map imol-map)
       (set! imol-rnase-map imol-map)
       (valid-map-molecule? imol-map))))


(greg-testcase "Map Sigma " #t 
   (lambda ()

     (if (not (valid-map-molecule? imol-rnase-map))
	 #f
	 (let ((v (map-sigma imol-rnase-map)))
	   (if (not (and (> v 0.2)
			 (< v 1.0)))
	       #f
	       (let ((v2 (map-sigma imol-rnase)))
		 (format #t "INFO:: map sigmas ~s ~s~%" v v2)
		 (eq? v2 #f)))))))
		     
(greg-testcase "Another Level Test" #t 
   (lambda ()
     (let ((imol-map-2 (another-level)))
       (valid-map-molecule? imol-map-2))))

(greg-testcase "NCS maps test" #t 
   (lambda ()

     (if (not (valid-model-molecule? imol-rnase))
	 (begin
	   (format #t "imol-rnase not valid~%")
	   (throw 'fail)))

     (if (not (valid-map-molecule? imol-rnase-map))
	 (begin
	   (format #t "imol-rnase-map not valid~%")
	   (throw 'fail)))

     (let ((n-mols (graphics-n-molecules)))
       ;; try to make it trip up by doing it twice:
       (let ((imol-map-2 (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))
	 (make-dynamically-transformed-ncs-maps imol-rnase imol-rnase-map)
	 (make-dynamically-transformed-ncs-maps imol-rnase imol-map-2)
	 ;; 2*2 + 1 new maps should have been made
	 (let ((n-new (graphics-n-molecules)))
	   (if (not (=  n-new (+ n-mols 5)))
	       (begin
		 (print-molecule-names)
		 (format #t "no match in number of molecules ~s ~s~%"
			 n-mols n-new)
		 (throw 'fail))
	       #t))))))



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
     
     (if (= (recentre-on-read-pdb) 0)
	 (set-recentre-on-read-pdb 1))

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
				    (format #t "   Sphere mol: there are ~s residues in chain ~s~%" 
					    n-residues chain-id)

				    (let loop ((residue-list (number-list 0 (- n-residues 1)))
					       (chain-n-atoms 0))

				      (cond 
				       ((null? residue-list) 
					chain-n-atoms)
				       (else 
					(let ((serial-number (car residue-list)))
					  (let ((res-name (resname-from-serial-number imol-sphere chain-id serial-number))
						(res-no   (seqnum-from-serial-number  imol-sphere chain-id serial-number))
						(ins-code (insertion-code-from-serial-number imol-sphere chain-id serial-number)))
					    (let ((residue-atoms-info (residue-info imol-sphere chain-id res-no ins-code)))
					      (loop (cdr residue-list) (+ (length residue-atoms-info) chain-n-atoms))))))))))
				(chain-ids imol-sphere)))))
	     
	     (format #t "   Found ~s sphere atoms ~%" n-atoms)

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
       (rotate-y-scene (rotate-n-frames 200) 0.1)
       (delete-residue imol-frag "B" 10 "")
       (delete-residue imol-frag "B" 11 "")
       (delete-residue imol-frag "B" 12 "")
       (rotate-y-scene (rotate-n-frames 200) 0.1)
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

		 ; we can only do this if we have graphics, must test
		 ; first.
		 ; (run-gtk-pending-events)
		 
		 ;; with new Molprobity rotamer probabilities,
		 ;; residues 1 and 2 are no longer "not recognised"
		 ;; they are in fact, just low probabilites.
		 ;; 
		 (if (and (string=? anal-str-a1 "VAL")
			  (string=? anal-str-a2 "VAL")
			  (string=? anal-str-a3 "Missing Atoms"))

		     ;; Now we test that the probabilites of the
		     ;; rotamer is correct:
		     (let ((pr-1 (list-ref (list-ref rotamer-anal 0) 3))
			   (pr-2 (list-ref (list-ref rotamer-anal 1) 3)))
		       (and (< pr-1 0.3)
			    (< pr-2 0.3)
			    (> pr-1 0.0)
			    (> pr-2 0.0)))
		     
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
			   (format #t "   For atom ~s ~s returning occ ~s~%" 
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

       (if (not (valid-model-molecule? imol-frag))
	   (begin
	     (format #t "bad molecule for reading coords in Alt Conf Occ test~%")
	     (throw 'fail))
	   (let ((occ-sum-pre (get-occ-sum imol-frag)))
	     
	     (let ((replace-state (refinement-immediate-replacement-state)))
	       (set-refinement-immediate-replacement 1)
	       (regularize-zone imol-frag "X" 15 15 "A")
	       (accept-regularizement)
	       (if (= replace-state 0)
		   (set-refinement-immediate-replacement 0))
	       
	       (let ((occ-sum-post (get-occ-sum imol-frag)))
		 
		 (format #t "   test for closeness: ~s ~s~%" occ-sum-pre occ-sum-post)	    
		 (close-float? occ-sum-pre occ-sum-post))))))))



;; This test we expect to fail until the CISPEP correction code is in
;; place (using mmdb-1.10+).
;; 
;; 20080812: OK, so we have mmdb-1.11, we expect this to pass (note
;; that an unexpected pass causes an overall fail - and the tar file
;; is not built).
;; 
(greg-testcase "Correction of CISPEP test" #t
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
	  (pepflip cis-pep-mol chain-id resno ins-code) 
	  (let ((res-type (residue-name cis-pep-mol chain-id resno ins-code)))
	    (if (not (string? res-type))
		(throw 'fail)
		(begin
		  (mutate cis-pep-mol chain-id resno "" "GLY")
		  (with-auto-accept
		   (refine-zone cis-pep-mol chain-id resno (+ resno 1) "")
		   (accept-regularizement)
		   (mutate cis-pep-mol chain-id resno "" res-type)
		   (auto-fit-best-rotamer resno "" ins-code chain-id cis-pep-mol 
					  (imol-refinement-map) 1 1)
		   (refine-zone cis-pep-mol chain-id resno (+ resno 1) "")
		   (accept-regularizement))
		  ;; should be repaired now.  Write it out.
		  
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
				    #t)))))))))))))



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




(greg-testcase "Setting multiple atom attributes" #t
   (lambda ()

     (if (not (valid-model-molecule? imol-rnase))
	 (begin 
	   (format #t "   Error invalid imol-rnase~%")
	   (throw 'fail))
	      
         (let* ((x-test-val 2.1)
		(y-test-val 2.2)
		(z-test-val 2.3)
		(o-test-val 0.5)
		(b-test-val 44.4)
		(ls (list (list imol-rnase "A" 2 "" " CA " "" "x" x-test-val)
			  (list imol-rnase "A" 2 "" " CA " "" "y" y-test-val)
			  (list imol-rnase "A" 2 "" " CA " "" "z" z-test-val)
			  (list imol-rnase "A" 2 "" " CA " "" "occ" o-test-val)
			  (list imol-rnase "A" 2 "" " CA " "" "b" b-test-val))))

           (set-atom-attributes ls)
	   (let ((atom-ls (residue-info imol-rnase "A" 2 "")))
	     (let f ((atom-ls atom-ls))
	       (cond 
		((null? atom-ls) (throw 'fail))
		(else 
		 (let* ((atom (car atom-ls))
			(compound-name (car atom))
			(atom-name (car compound-name)))
		   (if (string=? atom-name " CA ")
		       (let* ((xyz (car (cdr (cdr atom))))
			      (x (car xyz))
			      (y (car (cdr xyz)))
			      (z (car (cdr (cdr xyz))))
			      (occ-b-ele (car (cdr atom)))
			      (occ (car occ-b-ele))
			      (b (car (cdr occ-b-ele))))
			      
			 (if 
			  (and 
			   (close-float? x   x-test-val)
			   (close-float? y   y-test-val)
			   (close-float? z   z-test-val)
			   (close-float? occ o-test-val)
			   (close-float? b   b-test-val))
			  #t ; success
			  (begin
			    (format #t "   Error in setting multiple atom attributes~%")
			    (format #t "   These are not close: ~%")
			    (format #t "      ~s ~s~%"   x x-test-val)
			    (format #t "      ~s ~s~%"   y y-test-val)
			    (format #t "      ~s ~s~%"   z z-test-val)
			    (format #t "      ~s ~s~%" occ o-test-val)
			    (format #t "      ~s ~s~%"   b b-test-val)
			    #f)))
		       (f (cdr atom-ls))))))))))))



;; 
(greg-testcase "Tweak Alt Confs on Active Residue" #t 
    (lambda ()

      ;; did it reset it?
      (define (matches-alt-conf? imol chain-id resno inscode atom-name-ref alt-conf-ref)
	(let ((atom-ls (residue-info imol chain-id resno inscode)))
	  (if (not (list? atom-ls))
	      (begin
		(format #t "No atom list found - failing.") 
		(throw 'fail))
	      (let loop ((atom-ls atom-ls))
		(cond
		 ((null? atom-ls) #f)
		 (else 
		  (let* ((atom (car atom-ls))
			 (compound-name (car atom))
			 (atom-name (car compound-name))
			 (alt-conf (car (cdr compound-name))))
		    ; (format #t "DEBUG:: atom: ~s~%" atom)
		    (if (not (string=? atom-name atom-name-ref))
			(loop (cdr atom-ls))
			(string=? alt-conf alt-conf-ref)))))))))

      ;; main line
      (let ((chain-id "B")
	    (resno 58)
	    (inscode "")
	    (atom-name " CB ")
	    (new-alt-conf-id "B"))
	(set-go-to-atom-molecule imol-rnase)
	(set-go-to-atom-chain-residue-atom-name chain-id resno " CA ")
	(set-atom-string-attribute imol-rnase chain-id resno inscode
				   atom-name "" "alt-conf" new-alt-conf-id)
	(if (not (matches-alt-conf? imol-rnase chain-id resno inscode 
				    atom-name new-alt-conf-id))
	    (begin
	      (format #t "   No matching pre CB altconfed - failing.~%")
	      (throw 'fail)))
	(sanitise-alt-confs-in-residue imol-rnase chain-id resno inscode)
	(if (not (matches-alt-conf? imol-rnase chain-id resno inscode atom-name ""))
	    (begin
	      (format #t "   No matching post CB (unaltconfed) - failing.~%")
	      (throw 'fail)))
	#t)))

	   

(greg-testcase "Libcif horne" #t
   (lambda ()

     (if (not (file-exists? horne-pdb))
	 (begin
	   (format #t "file ~s not found - skipping test~%" horne-pdb)
	   (throw 'untested))

	 (if (not (file-exists? horne-cif))
	     (begin
	       (format #t "file ~s not found - skipping test~%" horne-cif)
	       (throw 'untested))

	     (let ((imol (read-pdb horne-pdb)))
	       (if (not (valid-model-molecule? imol))
		   (begin
		     (format #t "bad molecule number~%" imol)
		     (throw 'fail))
		   
		   (begin
		     (with-auto-accept
		      (regularize-zone imol "A" 1 1 "")) ; refine, no dictionary.  Used 
                                                         ; to hang here with no graphics.
		     (read-cif-dictionary horne-works-cif)
		     (with-auto-accept
		      (regularize-zone imol "A" 1 1 ""))
		     (newline) (newline)
		     (newline) (newline)
		     (read-cif-dictionary horne-cif)
		     (with-auto-accept
		      (regularize-zone imol "A" 1 1 ""))
		     (let ((cent (molecule-centre imol)))
		       (if (not (< (car cent) 2000))
			   (begin 
			     (format #t "Position fail 2000 test: ~s in ~s~%"
				     (car cent) cent)
			     (throw 'fail))
			   #t ; pass
			   )))))))))


(greg-testcase "Refmac Parameters Storage" #t
   (lambda ()

     (let* ((arg-list (list rnase-mtz "/RNASE3GMP/COMPLEX/FWT" "/RNASE3GMP/COMPLEX/PHWT" "" 0 0 1 
			    "/RNASE3GMP/COMPLEX/FGMP18" "/RNASE3GMP/COMPLEX/SIGFGMP18" 
			    "/RNASE/NATIVE/FreeR_flag" 1))
	    (imol (apply make-and-draw-map-with-refmac-params arg-list)))

       (if (not (valid-map-molecule? imol))
	   (begin
	     (format #t "  Can't get valid refmac map molecule from ~s~%" rnase-mtz)
	     (throw 'fail)))
       
       (let ((refmac-params (refmac-parameters imol)))
	 
;; 	 (format #t "          comparing: ~s~% with refmac params: ~s~%" arg-list refmac-params)
	 (if (not (equal? refmac-params arg-list))
	     (begin
	       (format #t "        non matching refmac params~%")
	       (throw 'fail))
	     #t)))))



(greg-testcase "The position of the oxygen after a mutation" #t
   (lambda ()

     ;; Return the o-pos (can be #f) and the number of atoms read.
     ;; 
     (define (get-o-pos hist-pdb-name)
       (call-with-input-file hist-pdb-name
	 (lambda (port)
	   
	   (let loop ((line (read-line port))
		      (o-pos #f)
		      (n-atoms 0))
	     (cond
	      ((eof-object? line) (values o-pos n-atoms))
	      (else 
	       (if (< (string-length line) 20)
		   (loop (read-line port) o-pos n-atoms)
		   (let ((atom-name (substring line 12 16))
			 (new-n-atoms (if (string=? (substring line 0 4) "ATOM")
					  (+ n-atoms 1)
					  n-atoms)))
		     (if (string=? atom-name " O  ")
			 (loop (read-line port) new-n-atoms new-n-atoms)
			 (loop (read-line port) o-pos new-n-atoms))))))))))
	    
     ;; main body
     ;; 
     (let ((hist-pdb-name "his.pdb"))
       
       (let ((imol (read-pdb (append-dir-file greg-data-dir "val.pdb"))))
	 (if (not (valid-model-molecule? imol))
	     (begin
	       (format #t "   failed to read file val.pdb~%")
	       (throw 'fail))
	     (let ((mutate-state (mutate imol "C" 3 "" "HIS")))
	       (if (not (= 1 mutate-state))
		   (begin
		     (format #t "   failure in mutate function~%")
		     (throw 'fail))
		   (begin
		     (write-pdb-file imol hist-pdb-name)
		     (if (not (file-exists? hist-pdb-name))
			 (begin
			   (format #t "   file not found: ~s~%" hist-pdb-name)
			   (throw 'fail))
			 (call-with-values
			     (lambda () (get-o-pos hist-pdb-name))
			   (lambda (o-pos n-atoms)
			     
			     (if (not (= n-atoms 10))
				 (begin
				   (format #t "   found ~s atoms (not 10)~%" n-atoms)
				   (throw 'fail))
				 (if (not (number? o-pos))
				     (begin
				       (format #t "   Oxygen O position not found~%")
				       (throw 'fail))
				     (if (not (= o-pos 4))
					 (begin
					   (format #t "   found O atom at ~s (not 4)~%" o-pos)
					   (throw 'fail))
					 #t))))))))))))))


(greg-testcase "Deleting (non-existing) Alt conf and Go To Atom [JED]" #t
   (lambda ()

     ;; alt conf "A" does not exist in this residue:
     (delete-residue-with-altconf imol-rnase "A" 88 "" "A")
     ;; to activate the bug, we need to search over all atoms
     (active-residue) ; crash
     #t))


(greg-testcase "Mask and difference map" #t 
  (lambda ()

    (define (get-ca-coords imol-model resno-1 resno-2)
      (let loop ((resno resno-1)
	    (coords '()))
	(cond 
	 ((> resno resno-2) coords)
	 (else 
	  (let ((atom-info (atom-specs imol-model "A" resno "" " CA " "")))
	    (loop (+ resno 1)
		  (cons 
		   (cdr (cdr (cdr atom-info)))
		   coords)))))))

    (let ((d-1 (difference-map -1 2 -9))
	  (d-2 (difference-map 2 -1 -9)))

      (if (not (= d-1 -1))
	  (begin
	    (format #t "failed on bad d-1~%")
	    (throw 'fail)))

      (if (not (= d-2 -1))
	  (begin
	    (format #t "failed on bad d-2~%")
	    (throw 'fail))))

    (let* ((imol-map-nrml-res (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0))
	   (prev-sampling-rate (get-map-sampling-rate))
	   (nov-1 (set-map-sampling-rate 2.2))
	   (imol-map-high-res (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0))
	   (nov-2 (set-map-sampling-rate prev-sampling-rate))
	   (imol-model (handle-read-draw-molecule-with-recentre rnase-pdb 0)))
      
      (let ((imol-masked (mask-map-by-atom-selection imol-map-nrml-res imol-model
						     "//A/1-10" 0)))
	(if (not (valid-map-molecule? imol-map-nrml-res))
	    (throw 'fail))

	;; check that imol-map-nrml-res is good here by selecting
	;; regions where the density should be positive. The masked
	;; region should be 0.0
	;; 
	(let ((high-pts (get-ca-coords imol-model 20 30))
	      (low-pts  (get-ca-coords imol-model  1 10)))

	(let ((high-values (map (lambda(pt) 
				  (apply density-at-point imol-masked pt)) high-pts))
	      (low-values (map (lambda(pt) 
				 (apply density-at-point imol-masked pt)) low-pts)))

	  (format #t "high-values: ~s  low values: ~s~%" 
		  high-values low-values)


	  (if (not (< (apply + low-values) 0.000001))
	      (begin
		(format #t "Bad low values~%")
		(throw 'fail)))

	  (if (not (> (apply + high-values) 5))
	      (begin
		(format #t "Bad high values~%")
		(throw 'fail)))

	  (let ((diff-map (difference-map imol-masked imol-map-high-res 1.0)))
	    (if (not (valid-map-molecule? diff-map))
		(begin
		  (format #t "failure to make good difference map")
		  (throw 'fail)))

	    ;; now in diff-map low pt should have high values and high
	    ;; pts should have low values
	    
	    (let ((diff-high-values (map (lambda(pt) 
					   (apply density-at-point diff-map pt)) high-pts))
		  (diff-low-values (map (lambda(pt) 
					  (apply density-at-point diff-map pt)) low-pts)))
	   
	      (format #t "diff-high-values: ~s  diff-low-values: ~s~%" 
		      diff-high-values diff-low-values)

	      (if (not (< (apply + diff-high-values) 0.03))
		  (begin
		    (format #t "Bad diff high values: value: ~s target: ~s~%"
			    (apply + diff-high-values)  0.03)
		    (throw 'fail)))
	      
	      (if (not (< (apply + diff-low-values) -5))
		  (begin
		    (format #t "Bad diff low values ~s ~%" (apply + diff-low-values))
		    (throw 'fail)))
	      
	      #t))))))))


(greg-testcase "Make a glycosidic linkage" #t 
   (lambda ()

     (let* ((carbo "multi-carbo-coot-2.pdb")
	    (imol (greg-pdb carbo)))

       (if (not (valid-model-molecule? imol))
	   (begin 
	     (format #t "file not found: ~s~%" carbo)
	     #f)
	   
	   (let ((atom-1 (get-atom imol "A" 1 " O4 "))
		 (atom-2 (get-atom imol "A" 2 " C1 ")))
	     
	     (format #t "bond-length: ~s: ~%"
		     (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))

	     (let ((s (dragged-refinement-steps-per-frame)))
	       (set-dragged-refinement-steps-per-frame 300)
	       (with-auto-accept
		(regularize-zone imol "A" 1 2 ""))
	       (set-dragged-refinement-steps-per-frame s))

	     (let ((atom-1 (get-atom imol "A" 1 " O4 "))
		   (atom-2 (get-atom imol "A" 2 " C1 ")))

	       (format #t "bond-length: ~s: ~%"
		       (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))
	       
	       (bond-length-within-tolerance? atom-1 atom-2 1.439 0.04)))))))




(greg-testcase "Test for flying hydrogens on undo" #t 
   (lambda ()

     (let ((imol (greg-pdb "monomer-VAL.pdb")))
       
       (if (not (valid-model-molecule? imol))
	   (begin 
	     (format "   Failure to read monomer-VAL.pdb~%")
	     (throw 'fail)))

       (with-auto-accept (regularize-zone imol "A" 1 1 ""))
       (set-undo-molecule imol)
       (apply-undo)
       (with-auto-accept (regularize-zone imol "A" 1 1 ""))
       
       (let ((atom-1 (get-atom imol "A" 1 "HG11"))
	     (atom-2 (get-atom imol "A" 1 " CG1")))

	 (if (bond-length-within-tolerance? atom-1 atom-2 0.96 0.02)
	     #t
	     (begin
	       (if (and (list? atom-1)
			(list? atom-2))
		   (format "   flying hydrogen failure, bond length ~s, should be 0.96~%"
			   (bond-length-from-atoms atom-1 atom-2))
		   (format "   flying hydrogen failure, atoms: ~s ~s~%" atom-1 atom-2))
	       #f))))))



(greg-testcase "Test for mangling of hydrogen names from a PDB v 3.0" #t
   (lambda ()

     (let ((imol (greg-pdb "3ins-6B-3.0.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "Bad read of greg test pdb: 3ins-6B-3.0.pdb~%")
	     (throw 'fail)))

       (with-auto-accept (regularize-zone imol "B" 6 6 ""))
       (let ((atom-pairs 
	      (list 
	       (cons "HD11" " CD1")
	       (cons "HD12" " CD1")
	       (cons "HD13" " CD1")
	       (cons "HD21" " CD2")
	       (cons "HD22" " CD2")
	       (cons "HD23" " CD2"))))

	 (all-true?
	  (map
	   (lambda (pair)
	     (let ((atom-1 (get-atom imol "B" 6 (car pair)))
		   (atom-2 (get-atom imol "B" 6 (cdr pair))))
	       (if (bond-length-within-tolerance? atom-1 atom-2 0.96 0.02)
		   #t 
		   (begin
		     (format #t "Hydrogen names mangled from PDB ~s ~s~%"
			     atom-1 atom-2)
		     #f))))
	   atom-pairs))))))
	       

(greg-testcase "update monomer restraints" #t 
   (lambda () 

     (let ((atom-pair (list " CB " " CG "))
	   (m (monomer-restraints "TYR")))

       (if (not m)
	   (begin 
	     (format #t "update bond restraints - no momomer restraints~%")
	     (throw 'fail)))

       (let ((n (strip-bond-from-restraints atom-pair m)))
	 (set-monomer-restraints "TYR" n)
	 
	 (let ((imol (new-molecule-by-atom-selection imol-rnase "//A/30")))
	   
	   (with-auto-accept
	    (refine-zone imol "A" 30 30 ""))
	   
	   (let ((atom-1 (get-atom imol "A" 30 " CB "))
		 (atom-2 (get-atom imol "A" 30 " CG ")))
	     
	     (format #t "   Bond-length: ~s: ~%"
		     (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))

	     (if (not (bond-length-within-tolerance? atom-1 atom-2 2.8 0.6))
		 (begin
		   (format #t "   Fail 2.8 tolerance test~%")
		   #f)
		 (begin 
		   (format #t "pass intermediate 2.8 tolerance test~%")
		   (set-monomer-restraints "TYR" m)

		   (with-auto-accept
		    (refine-zone imol "A" 30 30 ""))
		 
		   (let ((atom-1 (get-atom imol "A" 30 " CB "))
			 (atom-2 (get-atom imol "A" 30 " CG ")))

		     (let ((post-set-plane-restraints (assoc-ref (monomer-restraints "TYR")
								 "_chem_comp_plane_atom")))
		       
		       ;; (format #t "plane-restraint-pre:  ~s~%" (assoc-ref m "_chem_comp_plane_atom"))
		       ;; (format #t "plane-restraint-post: ~s~%" post-set-plane-restraints)
		       
		       (let ((atom (car (car (cdr (car post-set-plane-restraints))))))
			 (if (string=? atom " CB ")
			     (format #t "  OK plane atom ~s~%" atom)
			     (begin
			       (format #t "FAIL plane atom ~s~%" atom)
			       (throw 'fail)))))
		     
		     (format #t "   Bond-length: ~s: ~%"
			     (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))
		     (bond-length-within-tolerance? atom-1 atom-2 1.512 0.04))))))))))


(greg-testcase "Change Chain IDs and Chain Sorting" #t 
   (lambda () 
     (let ((imol (copy-molecule imol-rnase)))
       (change-chain-id imol "A" "D" 0  0  0)
       (change-chain-id imol "B" "E" 1 80 90)  
       (change-chain-id imol "B" "F" 1 70 80)
       (change-chain-id imol "B" "G" 1 60 70)
       (change-chain-id imol "B" "J" 1 50 59)
       (change-chain-id imol "B" "L" 1 40 49)
       (change-chain-id imol "B" "N" 1 30 38)
       (change-chain-id imol "B" "Z" 1 20 28)

;        (sort-chains imol)

;       (let ((c (chain-ids imol)))
;	 (chains-in-order? c)))))

       #t)))


(greg-testcase "Replace Fragment" #t
   (lambda ()

     (let* ((atom-sel-str "//A70-80")
	    (imol-rnase-copy (copy-molecule imol-rnase))
	    (imol (new-molecule-by-atom-selection imol-rnase atom-sel-str)))
       (if (not (valid-model-molecule? imol))
	   (throw 'fail))
		
       (translate-molecule-by imol 11 12 13)
       (let ((reference-res (residue-info imol-rnase "A" 75 ""))
	     (moved-res     (residue-info imol       "A" 75 "")))

	 (replace-fragment imol-rnase-copy imol atom-sel-str)
	 (let ((replaced-res (residue-info imol-rnase-copy "A" 75 "")))
	   
	   ;; now replace-res should match reference-res.
	   ;; the atoms of moved-res should be 20+ A away from both.
	   
;	   (format #t "reference-res: ~s~%" reference-res) 
;	   (format #t "    moved-res: ~s~%"     moved-res) 
;	   (format #t " replaced-res: ~s~%"  replaced-res) 

	   (if (not (all-true? (map atoms-match? moved-res replaced-res)))
	       (begin
		 (format #t "   moved-res and replaced-res do not match~%")
		 (throw 'fail)))

	   (if (all-true? (map atoms-match? moved-res reference-res))
	       (begin
		 (format #t "   fail - moved-res and replaced-res Match!~%")
		 (throw 'fail)))

	   (format #t "distances: ~s~%" (map atom-distance reference-res replaced-res))
	   (all-true? (map (lambda (d) (> d 20)) (map atom-distance reference-res replaced-res))))))))



(greg-testcase "Residues in Region of Residue" #t 
   (lambda ()

     (let ((rs (residues-near-residue imol-rnase (list "A" 40 "") 4)))
       (if (not (= (length rs) 6))
           (begin
             (format #t "wrong number of neighbours ~s ~s~%" (length rs) rs)
             #f)
           (begin
             (format #t "found ~s neighbours ~s~%" (length rs) rs)
             #t)))))


(greg-testcase "Empty molecule on type selection" #t 
   (lambda () 
     (let ((imol1 (new-molecule-by-residue-type-selection imol-rnase "TRP"))
	   (imol2 (new-molecule-by-residue-type-selection imol-rnase "TRP")))
       (if (not (= imol1 -1))
	  (begin
	    (format #t "failed on empty selection 1 gives not imol -1 ~%")
	    (throw 'fail)))
       (if (not (= imol2 -1))
	  (begin
	    (format #t "failed on empty selection 2 gives not imol -1 ~%")
	    (throw 'fail)))
       #t)))


(greg-testcase "Set Rotamer" #t 
   (lambda ()

     (let ((chain-id "A")
	   (resno 51))

       (let ((n-rot (n-rotamers -1 "ZZ" 45 "")))
	 (if (not (= n-rot -1))
	     (throw 'fail)))
       
       (let ((n-rot (n-rotamers imol-rnase "Z" 45 "")))
	 (if (not (= n-rot -1))
	     (throw 'fail)))

       (let ((residue-pre (residue-info imol-rnase chain-id resno "")))

	 ;; note that the rotamer number is 0-indexed (unlike the rotamer
	 ;; number dialog)
	 (set-residue-to-rotamer-number imol-rnase chain-id resno "" 1)
	 
	 (let ((residue-post (residue-info imol-rnase chain-id resno "")))
	   
	   (if (not (= (length residue-pre) (length residue-post)))
	       (throw 'fail))
	   
	   ;; average dist should be > 0.1 and < 0.3.
	   ;; 
	   (let ((dists 
		  (map (lambda (atom-number)
			 (let ((atom-pre (list-ref residue-pre atom-number))
			       (atom-post (list-ref residue-post atom-number)))
			   (let ((d (atom-distance atom-pre atom-post)))
			     d)))
		       (range (length residue-pre)))))
	     (map (lambda (d)
		    (if (> d 0.6)
			(throw 'fail))
		    (if (< d 0.0)
			(throw 'fail)))
		  dists)
	     #t)))))) ;; everything OK


		   

(greg-testcase "Align and mutate a model with deletions" #t 
   (lambda ()

     (define (residue-in-molecule? imol chain-id resno ins-code)
       (let ((r (residue-info imol chain-id resno ins-code)))
;;	 (format #t "::: res-info: ~s ~s ~s ~s -> ~s~%" imol chain-id resno ins-code r)
	 (if r
	     #t
	     #f)))

     ;; in this PDB file 60 and 61 have been deleted. Relative to the
     ;; where we want to be (tutorial-modern.pdb, say) 62 to 93 have
     ;; been moved to 60 to 91
     ;; 
     (let ((imol (greg-pdb "rnase-A-needs-an-insertion.pdb"))
	   (rnase-seq-string (file->string rnase-seq)))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "   Missing file rnase-A-needs-an-insertion.pdb~%")
	     (throw 'fail))
	   (begin

	     (align-and-mutate imol "A" rnase-seq-string)
	     (write-pdb-file imol "mutated.pdb")

	     (all-true? 
	      (map 
	       (lambda (residue-info)
		 (let ((residue-spec (cdr residue-info))
		       (expected-status (car residue-info)))
		   (format #t "::::: ~s ~s ~s~%" residue-spec 
			   (apply residue-in-molecule? residue-spec)
			   expected-status)
		   (eq? 
		    (apply residue-in-molecule? residue-spec)
		    expected-status)))
	       (list 
		(list #f imol "A"  1 "")
		(list #t imol "A"  4 "")
		(list #t imol "A" 59 "")
		(list #f imol "A" 60 "")
		(list #f imol "A" 61 "")
		(list #t imol "A" 92 "")
		(list #f imol "A" 94 "")
		))))))))
