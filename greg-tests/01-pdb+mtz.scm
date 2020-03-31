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


(set-map-radius 4.5) ;; faster


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


(greg-testcase "Don't crash on empty NCS from mmCIF file" #t
   (lambda ()

     (let ((imol (greg-pdb "2WF6.cif")))
       (format #t "   closing molecule number ~s~%" imol)
       (close-molecule imol)
       #t)))


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
	   (format #t "   pre and post ins codes: ~s ~s~%" ins-1 ins-2)
	   (if (not (string=? ins-2 "A"))
	       (begin
		 (format #t "Fail ins code set: ~s is not \"A\"~%" ins-2)
		 (throw 'fail)))

	   (write-pdb-file frag-pdb "post-ins-change.pdb")

	   (let ((test-expected-results
		  (list 

		   ;; note: 68B -> 70 doesn't happen unless we are on 68B (rc distance check)
		   ;; 
		   (cons (goto-next-atom-maybe "A" 67 ""  " CA ") (list "A" 68 ""  " CA "))
		   (cons (goto-next-atom-maybe "A" 68 "A" " CA ") (list "A" 68 "B" " CA "))
		   (cons (goto-next-atom-maybe "A" 68 "B" " CA ") (list "A" 68 "B" " CA ")) 
		   (cons (goto-prev-atom-maybe "A" 70 ""  " CA ") (list "A" 68 "B" " CA "))
		   (cons (goto-prev-atom-maybe "A" 68 "B" " CA ") (list "A" 68 "A" " CA "))
		   (cons (goto-prev-atom-maybe "A" 68 "A" " CA ") (list "A" 68  "" " CA "))
		   (cons (goto-prev-atom-maybe "A" 68 ""  " CA ") (list "A" 66  "" " CA ")))))

;; Note:: these are taken out beause goto-next-atom-maybe now checks
;; where the screen centre is before the move, so the previously
;; expected results no longer apply.  We could recentre and in these
;; tests - (todo...).
;; 		   (cons (goto-next-atom-maybe "A" 68 "B" " CA ") (list "A" 70 ""  " CA "))
;; 		   (cons (goto-next-atom-maybe "D" 10 ""  " O  ") (list "A" 62  "" " CA "))

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

;; remove this test because the PTR in the dictionary has the wrong chirality. And
;; now is not the time to fix it.
;;
(if #f
(greg-testcase "Replace Residue gets correct residue number" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (mutate-by-overlap imol "A" 86 "PTR")
       (let ((rn (residue-name imol "A" 86 "")))
	 (if (not (string? rn))
	     #f
	     (if (not (string=? rn "PTR"))
		 #f
		 ;; OK, did the the refinement run OK? Check the C-N distance
		 (let ((N-atom (get-atom imol "A" 86 "" " N  " ""))
		       (C-atom (get-atom imol "A" 85 "" " C  " "")))

		   (let ((dd (bond-length-from-atoms N-atom C-atom)))
		     (if (> dd 1.4)
			 #f
			 (if (< dd 1.25)
			     #f
			     (begin
			       (format #t "C-N dist good enough: ~s~%" dd)
			       #t)))))))))))
)


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


(greg-testcase "Auto-read bad MTZ test" #t 
   (lambda ()

     (let ((mtz-list (list "xx-missing.mtz" 
			   (append-dir-file greg-data-dir "broken.mtz"))))

       (map (lambda (file-name)
	      (let ((r (auto-read-make-and-draw-maps file-name)))
		
		(format #t "   got status: ~s~%" r)))
	    mtz-list))
     #t ;; no crash
))



(greg-testcase "Map Sigma " #t 
   (lambda ()

     (if (not (valid-map-molecule? imol-rnase-map))
	 #f
	 (let ((v (map-sigma imol-rnase-map)))
	   (if (not (and (> v 0.2)
			 (< v 1.0)))
	       #f
	       (let ((v2 (map-sigma imol-rnase)))
		 (format #t "   INFO:: map sigmas ~s ~s~%" v v2)
		 (eq? v2 #f)))))))
		     
(greg-testcase "Another Level Test" #t 
   (lambda ()
     (let ((imol-map-2 (another-level)))
       (valid-map-molecule? imol-map-2))))


(greg-testcase "Sharpen map from map" #t
   (lambda ()
     ;; don't crash
     (let* ((mtz-file-name (append-dir-file greg-data-dir "3hfl_sigmaa.mtz"))
	    (imol-map (make-and-draw-map mtz-file-name
					 "2FOFCWT" "PH2FOFCWT" "" 0 0)))
       (if (not (valid-map-molecule? imol-map))
	   (begin
	     (format #t "fail to get map from 3hfl_sigmaa.mtz~%")
	     (throw 'fail)))

       (export-map imol-map "test-3hfl.map")
       (let ((new-mol (handle-read-ccp4-map "test-3hfl.map" 0)))
	 (if (not (valid-map-molecule? new-mol))
	     (begin
	       (format #t "fail to get map from 3hfl map~%")
	       (throw 'fail)))
	 (sharpen new-mol 5.0)
	 #t ; didn't crash
	 ))))



(greg-testcase "db-main makes mainchain" #t
   (lambda ()	      

     (read-pdb ".")

     (let ((imol (read-pdb rnase-pdb)))
       (db-mainchain imol"A" 10 20 "forward")
       #t ;; didn't hang
       )))


(greg-testcase "Negative Residues in db-mainchain don't cause a crash" #t 
   (lambda ()

     ;; Oliver Clarke spotted this bug

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (if (not (valid-model-molecule? imol))
	   #f
	   (begin
	     (renumber-residue-range imol "A" 1 50 -20)
	     (let ((imol-mc-1 (db-mainchain imol "A" -19  6 "forwards"))
		   (imol-mc-2 (db-mainchain imol "B"  10 30 "backwards")))
	       ;; (write-pdb-file imol-mc-1 "dbmc-1.pdb")
	       ;; (write-pdb-file imol-mc-2 "dbmc-2.pdb")
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
		       (set-default-temperature-factor-for-new-atoms 45)
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

				       ;; now test that the new atoms have the correct
				       ;; B factor.
				       (let ((new-atoms (residue-info imol "A" 0 "")))
					 (if (not (> (length new-atoms) 4))
					     (begin
					       (format #t "Not enough new atoms ~s~%" new-atoms)
					       (throw 'fail))
					     (if (not 
						  (all-true? 
						   (map (lambda (atom)
							  (close-float? 45
									(list-ref (list-ref atom 1) 1)))
							new-atoms)))
						 (begin 
						   (format #t "Fail b-factor test ~s~%" new-atoms)
						   #f)
						 
						 #t)))))))))))))))))


(greg-testcase "Adding residue by phi psi, no crash" #t
   (lambda ()

     (let ((imol (greg-pdb "frag-2wot.pdb")))
       (if (not (valid-model-molecule? imol))
	   (throw 'fail))
       (let ((v1 (add-terminal-residue-using-phi-psi imol "A" 275 "ALA" -60 -60)))
	 (if (not (= v1 1))
	     (throw 'fail)))
       (let ((v2 (add-terminal-residue-using-phi-psi imol "A" 276 "ALA" -60 -60)))
	 (if (not (= v2 1))
	     (throw 'fail)))
       (let ((v3 (add-terminal-residue-using-phi-psi imol "XX" 276 "ALA" -60 -60)))
	 (if (not (= v3 0))
	     (throw 'fail)))
       #t)))

(greg-testcase "Add Terminal Residue O Position" #t
   (lambda ()

     (let* ((imol (greg-pdb "tutorial-modern.pdb"))
	    (mtz-file-name (append-dir-file
			    greg-data-dir
			    "rnasa-1.8-all_refmac1.mtz"))
	    (imol-map (make-and-draw-map mtz-file-name "FWT" "PHWT" "" 0 0)))

       ;; move the O close to where the N will end up - a bad place.
       ;; (we check that it moves there)
       ;;
       (let ((attribs (list (list imol "A" 93 "" " O  " "" "x" 58.5)
			    (list imol "A" 93 "" " O  " "" "y"  2.9)
			    (list imol "A" 93 "" " O  " "" "z" -1.9))))
	 (set-atom-attributes attribs)
	 (let ((O-atom-o (get-atom imol "A" 93 "" " O  " "")))
	   (with-no-backups imol
			    (let ((dummy 'making-a-thunk))
			      (add-terminal-residue imol "A" 93 "ALA" 1)
			      (let ((N-atom-n (get-atom imol "A" 94 "" " N  " "")))
			      (let ((O-atom-n (get-atom imol "A" 93 "" " O  " "")))
				(let ((dd-1 (bond-length-from-atoms O-atom-o N-atom-n))
				      (dd-2 (bond-length-from-atoms O-atom-n N-atom-n)))
				  (format #t "Add terminal residue bond check dd-1: ~s~%" dd-1)
				  (format #t "Add terminal residue bond check dd-2: ~s~%" dd-2)
				  ;; the new N will not always go into the same place
				  ;; allow a bit more generosity in position difference
				  (if (not (< dd-1 0.4)) ;; 0.4 should be generous enough
				      #f
				      (if (not (> dd-2 1.9))
					  #f
					  #t ; hooray
					  ))))))))))))

				       
(greg-testcase "Select by Sphere" #t
   (lambda ()

     (let ((imol-sphere (new-molecule-by-sphere-selection 
			 imol-rnase 
			 24.6114959716797 24.8355808258057 7.43978214263916
			 3.6 1))) ;; allow partial residues in output

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



(greg-testcase "HIS with unusual atom order rotates correct fragment for 180 sidechain flip" #t
   (lambda ()

     (let ((imol (greg-pdb "eleanor-HIS.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "BAD imol for 180 sidechain flip test imol: ~s~%" imol)
	     (throw 'fail))

	   (let ((N-atom-o   (get-atom imol "A" 111 "" " N  " ""))
		 (ND1-atom-o (get-atom imol "A" 111 "" " ND1" "")))
	     (with-no-backups imol
	      (do-180-degree-side-chain-flip imol "A" 111 "" ""))
	     (let ((N-atom-n   (get-atom imol "A" 111 "" " N  " ""))
		   (ND1-atom-n (get-atom imol "A" 111 "" " ND1" "")))

	       ;; the N-atom stays still
	       ;; the ND1 atom moves by > 1A.

	       (let ((dd-1 (bond-length-from-atoms N-atom-o N-atom-n))
		     (dd-2 (bond-length-from-atoms ND1-atom-o ND1-atom-n)))

		 (format #t "dd-1: ~s dd-2: ~s~%" dd-1 dd-2)

		 (if (> dd-1 0.01)
		     (begin
		       (format #t "N atom moved - fail\n")
		       (throw 'fail))

		     (if (< dd-2 1.01)
			 (begin
			   (format #t "ND1 atom did not move enough - fail\n")
			   (throw 'fail))

			 #t)))))))))

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


(greg-testcase "Correct occupancies after auto-fit rotamer on alt-confed residue" #t
   (lambda ()

     (let* ((imol (greg-pdb "tutorial-modern.pdb"))
	    (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0))
	    (new-alt-conf (add-alt-conf imol "A" 93 "" "" 0)))

       (accept-regularizement) ;; Presses the OK button for the alt conf

       (auto-fit-best-rotamer 93 "A" "" "A" imol imol-map 1 0.01) 
       
       ;; Now test the ocupancies.  The problem was that we had been
       ;; ending up with atoms in the residue with occupancies of 0.0.
       ;; (They should be 0.8 and 0.2 - for me at least).  So test
       ;; that the atoms have occupancies of greater than 0.1 (I
       ;; suppose also test for occs > 0.85 would be good too).
       ;; 
       (let ((atoms (residue-info imol "A" 93 "")))
	 (map (lambda (atom)
		(let ((occupancy (car (list-ref atom 1)))
		      (alt-conf  (cadr (car atom))))
		  ;;(format #t "alt-conf: ~s occupancy: ~s   in atom ~s~%" 
		  ;; alt-conf occupancy atom)
		  (if (< occupancy 0.1)
		      (begin
			(format #t "bad occupancy in atom: ~s~%" atom)
			(throw 'fail)))
		  (if (> occupancy 0.85)
		      (begin
			(format #t "bad occupancy in atom: ~s~%" atom)
			(throw 'fail)))))
	      atoms))
       #t ; no fail -> success
       )))
	 

(greg-testcase "Rotamers work on MSE" #t
   (lambda ()

     (let ((imol (greg-pdb "pdb3knw.ent")))

       (let ((se-1 (get-atom imol "A" 89 "" "SE  ")))

	 (set-residue-to-rotamer-number imol "A" 89 "" "" 3)

	 (let ((se-2 (get-atom imol "A" 89 "" "SE  ")))

	   (format #t "    se-1: ~s~%" se-1)
	   (format #t "    se-2: ~s~%" se-2)

 	   (if (atoms-match? se-1 se-2) 
 	       #f
 	       #t ;; the SE moved, test passes
 	       ))))))




(greg-testcase "Hs are correctly swapped on a TYR" #t
   (lambda ()

     (let ((imol (greg-pdb "pdb1py3.ent")))
       (if (not (valid-model-molecule? imol))
	   (throw 'missing-or-bad-pdb1py3)
	   
	   (begin

	     ;; the pdb file contains hydrogens and nomenclature
	     ;; errors, lets see if we can fix them
	     ;; 
	     (fix-nomenclature-errors imol)
	     (let* ((atoms (residue-info imol "C" 54 "")))
	       (if (not (list? atoms))
		   (throw 'atoms-not-a-list)
		   (let ((cd1 (get-atom-from-residue " CD1" atoms ""))
			 (cd2 (get-atom-from-residue " CD2" atoms ""))
			 (hd1 (get-atom-from-residue " HD1" atoms ""))
			 (hd2 (get-atom-from-residue " HD2" atoms ""))
			 (ce1 (get-atom-from-residue " CE1" atoms ""))
			 (ce2 (get-atom-from-residue " CE2" atoms ""))
			 (he1 (get-atom-from-residue " HE1" atoms ""))
			 (he2 (get-atom-from-residue " HE2" atoms "")))
		     (let ((bonded-atoms (list (list cd1 hd1)
					       (list cd2 hd2)
					       (list ce1 he1)
					       (list ce2 he2))))
		       (let ((results (map (lambda (atom-1 atom-2) 
					     (bond-length-within-tolerance? atom-1 atom-2 0.93 0.02))
					   (map car  bonded-atoms)
					   (map cadr bonded-atoms))))
			 (format #t "results: ~s~%" results)
			 (all-true? results)))))))))))


(greg-testcase "Splitting residue leaves no atoms with negative occupancy" #t 
   (lambda ()

     ;; return #f if there are negative occupancies
     ;; 
     (define (check-for-negative-occs occs)
       (if (not (list? occs))
	   #f
	   (let loop ((occs occs))
	     (cond
	      ((null? occs) #t)
	      ((< (car occs) 0.0) #f)
	      (else
	       (loop (cdr occs)))))))
	       
     (let* ((imol (greg-pdb "tutorial-modern.pdb"))
	    (mtz-file-name (append-dir-file
			    greg-data-dir
			    "rnasa-1.8-all_refmac1.mtz"))
	   (imol-map (make-and-draw-map mtz-file-name "FWT" "PHWT" "" 0 0)))

       (zero-occupancy-residue-range imol "A" 37 37)
       (let* ((new-alt-conf (add-alt-conf imol "A" 37 "" "" 0))
	      (atoms (residue-info imol "A" 37 ""))
	      (occs (map (lambda (atom) (car (list-ref atom 1))) atoms)))
	 
	 (if (< (length occs) 5)
	     #f ;; too few atoms

	     (let ((occs-ok-status (check-for-negative-occs occs)))
	       (if (not occs-ok-status)
		   (format #t "Ooops: bad occupancies: ~s~%" occs))
	       (close-molecule imol)
	       (close-molecule imol-map)
	       occs-ok-status))))))




(greg-testcase "Pepflip flips the correct alt confed atoms" #t
   (lambda () 

     (let ((imol (greg-pdb "alt-conf-pepflip-test.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "   not found greg test pdb alt-conf-pepflip-test.pdb~%")
	     (throw 'fail)))
       
       ;; get the original coords 
       (let ((c-atom-A-o (get-atom imol "A" 65 "" " C  " "A"))
	     (o-atom-A-o (get-atom imol "A" 65 "" " O  " "A"))
	     (n-atom-A-o (get-atom imol "A" 66 "" " N  " "A"))
	     (c-atom-B-o (get-atom imol "A" 65 "" " C  " "B"))
	     (o-atom-B-o (get-atom imol "A" 65 "" " O  " "B"))
	     (n-atom-B-o (get-atom imol "A" 66 "" " N  " "B")))
	   
	 (pepflip imol "A" 65 "" "B")

	 ;; get the new coords 
	 (let ((c-atom-A-n (get-atom imol "A" 65 "" " C  " "A"))
	       (o-atom-A-n (get-atom imol "A" 65 "" " O  " "A"))
	       (n-atom-A-n (get-atom imol "A" 66 "" " N  " "A"))
	       (c-atom-B-n (get-atom imol "A" 65 "" " C  " "B"))
	       (o-atom-B-n (get-atom imol "A" 65 "" " O  " "B"))
	       (n-atom-B-n (get-atom imol "A" 66 "" " N  " "B")))

	   ;; now, the *A-n atoms should match the position of the
	   ;; *A-o atoms: 

	   (let ((b1 (bond-length-from-atoms c-atom-A-o c-atom-A-n))
		 (b2 (bond-length-from-atoms o-atom-A-o o-atom-A-n))
		 (b3 (bond-length-from-atoms n-atom-A-o n-atom-A-n))
		 (b4 (bond-length-from-atoms c-atom-B-o c-atom-B-n))
		 (b5 (bond-length-from-atoms o-atom-B-o o-atom-B-n))
		 (b6 (bond-length-from-atoms n-atom-B-o n-atom-B-n)))

	     (if (not (and (< b1 0.001)
			   (< b2 0.001)
			   (< b3 0.001)))
		 (begin
		   (format #t "   bad! A conf moved ~s ~s ~s~%" b1 b2 b3)
		   (throw 'fail)))
	     
	     (if (not (and (> b4 0.8)
			   (> b5 2.0)
			   (> b6 0.5)))
		 (begin
		   (format #t "   bad! B conf moved too little ~s ~s ~s~%" b3 b4 b5)
		   (throw 'fail)))

	     #t ;; success then

	     ))))))



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
	  (pepflip cis-pep-mol chain-id resno ins-code "")
	  (let ((res-type (residue-name cis-pep-mol chain-id resno ins-code)))
	    (if (not (string? res-type))
		(throw 'fail)
		(begin
		  (mutate cis-pep-mol chain-id resno "" "GLY")
		  (with-auto-accept
		   (refine-zone cis-pep-mol chain-id resno (+ resno 1) ""))
		  (mutate cis-pep-mol chain-id resno "" res-type)
		  (auto-fit-best-rotamer resno "" ins-code chain-id cis-pep-mol 
					 (imol-refinement-map) 1 1)
		  (with-auto-accept
		   (refine-zone cis-pep-mol chain-id resno (+ resno 1) ""))
		  ;; should be repaired now.  Write it out.

		  (let ((tmp-file "tmp-fixed-cis.pdb"))
		    (write-pdb-file cis-pep-mol tmp-file)

		    (call-with-input-file tmp-file
		      (lambda (port)

			(let ((n-cispeps
			       (let loop ((line (read-line port))
					  (n-cispeps 0))
				 (cond
				  ((eof-object? line) n-cispeps)
				  ((string=? (substring line 0 4) "CISP") (loop (read-line port) (+ n-cispeps 1)))
				  ((string=? (substring line 0 4) "ATOM") n-cispeps)
				  (else
				   (loop (read-line port) n-cispeps))))))

			  (= n-cispeps 3))))))))))))


(greg-testcase "H on a N moves on cis-trans convert" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (let ((imol-2 (new-molecule-by-atom-selection imol "//A/1-10")))
	 (coot-reduce imol-2)
	 (let ((H-atom-o (get-atom imol-2 "A" 6 "" " H  " "")))
	     (with-no-backups imol-2 (cis-trans-convert imol-2 "A" 5 "")) ;; 5-6 peptide
	     (let ((H-atom-n   (get-atom imol-2 "A" 6 "" " H  " "")))
	       (let ((dd (bond-length-from-atoms H-atom-o H-atom-n)))
		 (close-molecule imol)
		 (close-molecule imol-2)
		 (format #t "dd: ~s~%" dd)
		 (> dd 1.4))))))))


(greg-testcase "HA on a ALA exists after mutation to GLY" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (let ((imol-2 (new-molecule-by-atom-selection imol "//A/5-11")))
	 (coot-reduce imol-2)
	 (let ((H-atom-o (get-atom imol-2 "A" 10 "" " HA " "")))
	   (if (not (list? H-atom-o))
	       (throw 'fail)

	       (begin
		 (mutate imol-2 "A" 10 "" "GLY")
		 (let ((H-atom-n (get-atom imol-2 "A" 10 "" " HA " "")))

		   (if (list? H-atom-n)
		       (begin
			 (format #t "atom still exists ~s~%" H-atom-n)
			 #f)

		       #t)))))))))



(greg-testcase "Refine Zone with Alt conf" #t 
   (lambda ()

     (let* ((imol (greg-pdb "tutorial-modern.pdb"))
	    (mtz-file-name (append-dir-file
			    greg-data-dir
			    "rnasa-1.8-all_refmac1.mtz"))
	   (imol-map (make-and-draw-map mtz-file-name "FWT" "PHWT" "" 0 0)))

       (set-imol-refinement-map imol-map)
       (let ((at-1 (get-atom imol "B" 72 "" " SG " "B")))
	 (with-auto-accept
	  (refine-zone imol "B" 72 72 "B"))
	 (let ((at-2 (get-atom imol "B" 72 "" " SG " "B")))
	   (let ((d (bond-length-from-atoms at-1 at-2)))
	     ;; the atom should move in refinement
	     (format #t "   refined moved: d=~s~%" d)

	     (if (< d 0.09) ;; 20180606-PE New style refinement (wrapping residue 
                           ;; vec means that refinement has changed) atoms are still
		           ;; still moving though. Hmm. Perhaps the A atoms should not
		           ;; move and we should test for that.
                           ;;
                           ;; 20120110 new-style NBCs means that the
                           ;; atoms move less here
		           ;; 20160608 - they move still less  (not sure
                           ;; why this time)
		 (begin
		   (format #t "   refined atom failed to move: d=~s~%" d)
		   (throw 'fail))
		 #t)))))))


(greg-testcase "Sphere Refine" #t
   (lambda () 

     (define (sphere-refine-here)
       (let ((active-atom (active-residue)))
	 (if (not (list? active-atom))
	     (begin 
	       (format #t "No active atom~%")
	       (throw 'fail)))
	 (let* ((centred-residue (list-head (cdr active-atom) 3))
		(imol (car active-atom))
		(other-residues (residues-near-residue imol 
						       centred-residue 3.2))
		(all-residues (if (list? other-residues)
				  (cons centred-residue other-residues)
				  (list centred-residue))))
	   
	   ;; (format #t " ===== in sphere-refine-here  imol: ~s residues: ~s~%" imol all-residues)
	   (with-auto-accept
	    (refine-residues imol all-residues)))))


     ;; main line

     ;; refine against something, temporary for single test
     ;;
     ;; (let ((imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))
     ;; (valid-map-molecule? imol-map))


     (let ((imol (greg-pdb "tutorial-add-terminal-1-test.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "   molecule pdb not found~%")
	     (throw 'fail)))

       (let ((new-alt-conf (add-alt-conf imol "A" 93 "" "" 0)))
	 
	 (set-go-to-atom-molecule imol)
	 (let ((success (set-go-to-atom-chain-residue-atom-name  "A" 93 " CA ")))
	   (if (= success 0)
	       (begin
		 (format #t "   failed to go to A93 CA atom~%")
		 (throw 'fail)))

	   (sphere-refine-here)

	   (set-go-to-atom-chain-residue-atom-name "A" 41 " CA ")
	   (sphere-refine-here)

	   ;; these should be a standard peptide bond length - however, if we don't find
	   ;; the peptide link, they get pushed apart.  Test for that. 
	   (let* ((atom-1 (get-atom imol "A" 40 "" " C  " ""))
		  (atom-2 (get-atom imol "A" 41 "" " N  " ""))
		  (bl (bond-length-from-atoms atom-1 atom-2)))
	     (format #t "======= got bond length ~s~%" bl)
	     (< bl 1.4)))))))




(greg-testcase "Refinement gives useful results" #t 
   (lambda () 

     ;; 
     (define (no-non-bonded ls)
       (cond 
	((null? ls) '())
	((string=? (car (car ls)) "Non-bonded")
	 (no-non-bonded (cdr ls)))
	(else 
	 (cons (car ls) (no-non-bonded (cdr ls))))))
	

     ;; return #f or a number, which is the current overweighting of
     ;; the density terms.  When not overweighted, the geometric chi
     ;; squareds will be 1.0.
     ;; 
     (define (weight-scale-from-refinement-results rr)
       (if (not (list? rr))
	   #f
	   (let* ((nnb-list (no-non-bonded (list-ref rr 2)))
		  (chi-squares (map (lambda (x) (list-ref x 2)) nnb-list))
		  (n (length chi-squares))
		  (sum (apply + chi-squares)))
	     (/ sum n))))
	     
		
     (let ((imol (greg-pdb "tutorial-modern.pdb"))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))
       (set-imol-refinement-map imol-map)

       (let loop ((results (refine-zone-with-full-residue-spec imol "A" 40 "" 43 "" "")))
	 (format #t "   refinement results: ~s~%" results)
	 (let ((ow (weight-scale-from-refinement-results results)))
	   (format #t "   ow factor: ~s~%" ow)
	   (if (and (< ow 1.1)
		    (> ow 0.9))
	       #t
	       (let ((new-weight (/ (matrix-state) (* ow ow))))
		 (format #t "   INFO:: setting refinement weight to ~s~%" new-weight)
		 (set-matrix new-weight)
		 (loop (refine-zone-with-full-residue-spec imol "A" 40 "" 43 "" "")))))))))



(greg-testcase "Neighbour-Refine doesn't destroy disulfide bonds" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))

       (set-imol-refinement-map imol-rnase-map)
       (set-go-to-atom-molecule imol)
       (set-go-to-atom-chain-residue-atom-name "B" 7 " CA ")

	(let* ((rc-spec (list "B" 7 ""))
	       (ls (residues-near-residue imol rc-spec 2.2))
	       (residues-list (cons rc-spec ls)))

	  ;; (format #t "debug::::: neighbour-refine test: residues-list: ~s~%" residues-list)

	  (with-auto-accept (refine-residues imol residues-list))
	  (let ((at-spec-1 (list "B"  7 "" " SG " ""))
		(at-spec-2 (list "B" 96 "" " SG " "")))
	    (let ((at-1 (get-atom-from-spec imol at-spec-1))
		  (at-2 (get-atom-from-spec imol at-spec-2)))
	      (let ((bl (bond-length-from-atoms at-1 at-2)))

		;; (format #t "bl-1: ~s~%" bl)

		(let ((state (bond-length-within-tolerance? at-1 at-2 2.0 0.05)))

		  ;; (format #t "bond-length within tolerance?: ~s~%" state)

		  (if (not state)
		      #f

		      ;; do it again
		      (begin
			(with-auto-accept (refine-residues imol residues-list))
			(let ((at-spec-1 (list "B"  7 "" " SG " ""))
			      (at-spec-2 (list "B" 96 "" " SG " "")))
			  (let ((at-1 (get-atom-from-spec imol at-spec-1))
				(at-2 (get-atom-from-spec imol at-spec-2)))
			    (let ((bl (bond-length-from-atoms at-1 at-2)))

			      ;; (format #t "bl-2: ~s~%" bl)

			      (let ((state (bond-length-within-tolerance? at-1 at-2 2.0 0.05)))

				state))))))))))))))





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


       
(greg-testcase "Backrub rotamer" #t
   (lambda ()

     ;; 
     (define (get-mover-list imol-1 imol-2 res-no)
       (let ((atoms-1 (residue-info imol-1 "A" res-no ""))
	     (atoms-2 (residue-info imol-2 "A" res-no "")))
	 (let loop ((atoms-1 atoms-1)
		    (atoms-2 atoms-2)
		    (mover-list '()))
	   (cond 
	    ((null? atoms-1) (sort mover-list string<?))
	    (else 
	     (let ((pos-1 (list-ref (car atoms-1) 2))
		   (pos-2 (list-ref (car atoms-2) 2)))
	       (if (all-true? (map close-float? pos-1 pos-2))
		   ;; the atom did not move
		   (loop (cdr atoms-1) (cdr atoms-2) mover-list)
		   
		   ;; did move
		   (loop (cdr atoms-1) (cdr atoms-2) (cons (car (car (car atoms-1)))
							   mover-list)))))))))

     ;; main line
     (set-imol-refinement-map imol-rnase-map)

     (let* ((imol (handle-read-draw-molecule-with-recentre 
		   (append-dir-file greg-data-dir "backrub-fragment.pdb") 0))
	    (imol-copy (copy-molecule imol)))

       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "backrub-rotamer backrub-fragment.pdb not found it seems~%")
	     (throw 'fail)))
		  
       (let ((status (backrub-rotamer imol-copy "A" 37 "" "")))

	 (if (not (= status 1))
	     (begin
	       (format #t "backrub-rotamer status failure~%")
	       (throw 'fail)))

	 ;; Now, N should have moved in 38 (nothing else)
	 ;; 
	 ;; C and O should have moved in 36 (nothing else)
	 ;; 
	 ;; all atoms in 37.
	 ;; 

	 (let* ((atoms-37 (residue-info imol "A" 37 ""))
		(calc-37-movers-pre (map car (map car atoms-37)))
		(calc-37-movers (sort calc-37-movers-pre string<?))
		(calc-38-movers (list " N  "))
		(calc-36-movers (list " C  " " O  ")))
	   
	   (let ((36-movers (get-mover-list imol imol-copy 36))
		 (37-movers (get-mover-list imol imol-copy 37))
		 (38-movers (get-mover-list imol imol-copy 38)))

;	     (format #t "calc-37-movers: ~s~%" calc-36-movers)
;	     (format #t " obs-37-movers: ~s~%"      36-movers)
;	     (format #t "calc-38-movers: ~s~%" calc-37-movers)
;	     (format #t " obs-38-movers: ~s~%"      37-movers)
;	     (format #t "calc-39-movers: ~s~%" calc-38-movers)
;	     (format #t " obs-39-movers: ~s~%"      38-movers)

	     (if (not (equal? 36-movers calc-36-movers))
		 (begin
		   (format "  fail on 36 movers similar~%")
		   (throw 'fail)))
	     
	     (if (not (equal? 37-movers calc-37-movers))
		 (begin
		   (format "  fail on 37 movers similar~%")
		   (throw 'fail)))
	     
	     (if (not (equal? 38-movers calc-38-movers))
		 (begin
		   (format "  fail on 38 movers similar~%")
		   (throw 'fail)))
	     
	     #t))))))

	   

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


(greg-testcase "OXT is added before TER record - add only one"  #t
   (lambda () 

     (let ((imol (greg-pdb "test-TER-OXT.pdb"))
	   (opdb "test-TER-OXT-added.pdb"))

       (if (not (valid-model-molecule? imol))
	   (begin 
	     (format #t "Failed to read test-TER-OXT.pdb")
	     (throw 'fail)))
	   
       (let ((add-status-1 (add-OXT-to-residue imol 14 "" "A"))
	     (add-status-2 (add-OXT-to-residue imol 14 "" "A")))

	 ;; The second add should be ignored, only 1 OXT allowed
	 ;; per residue.

	 (if (not (= add-status-1 1))
	     (begin 
	       (format #t "   add-status-1 not success - fail~%")
	       (throw 'fail)))

	 (if (not (= add-status-2 0))
	     (begin 
	       (format #t "   add-status-2 was success - fail~%")
	       (throw 'fail)))


	 (write-pdb-file imol opdb)
	 (call-with-input-file opdb
	   (lambda (port)
	     (let loop ((line (read-line port))
			(count 0)
			(ter-line #f)
			(oxt-line #f))

	       (cond
		((eof-object? line) #f) ; should not have got here!
		((string=? (substring line 0 3) "END")
		 #f)  ; should not have got here!
		((string=? (substring line 0 3) "TER")
		 (format #t "   found TER ~s~%" line)
		 (if (number? oxt-line)
		     #t; TER happens after OXT line
		     (loop (read-line port) (+ count 1) count oxt-line)))
		((string=? (substring line 13 16) "OXT")
		 (if (number? ter-line)
		     #f ; fail because TER has already happened!
		     (if (number? oxt-line)
			 (begin
			   (format #t "   Encountered another OXT! - fail~%")
			   #f) ; Fail because there an OXT was encountered before this one!
			 (loop (read-line port) (+ count 1) ter-line count))))
		(else 
		 (loop (read-line port) (+ count 1) ter-line oxt-line))))))))))




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


(greg-testcase "TER is at the end of a nucleotide after mutation" #t 
   (lambda ()

     ;; Before this test, if we mutated a nucleotide at the end of a
     ;; chain, then the TER record appeared in the PDB file before the
     ;; additional new atoms.  Wrongness.  James Parker bug.

     (let ((imol (greg-pdb "2yie-frag.pdb")))
       
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "failed to read 2yie-frag.pdb~%")
	     (throw 'fail)))

       (let ((status (mutate-base imol "X" 54 "" "C")))

	 (if (not (= status 1))
	   (begin
	     (format #t "failed to mutate 2yie-frag.pdb~%")
	     (throw 'fail)))

	 (write-pdb-file imol "2yie-mutated.pdb")
	 
	 (call-with-input-file "2yie-mutated.pdb"
	   (lambda (port)

	     (let loop ((line (read-line port))
			(ter-status #f))
	       (cond
		((eof-object? line) 
		 (throw 'fail))
		((string-contains line "END") 
		 #t) ;; good, we got to the end without seeing ATOM records 
		     ;; after a TER record.
		((string-contains line "TER") (loop (read-line port) #t))
		((string-contains line "ATOM")
		 (if ter-status 
		     #f ;; fail if we see an ATOM record after a TER record.
		     (loop (read-line port) #f)))
		(else 
		 (loop (read-line port) ter-status))))))))))


(greg-testcase "C7 is removed on mutation from a DC" #t
   (lambda ()

     (let ((imol (greg-pdb "4f8g.pdb")))

       (if (not (valid-model-molecule? imol))
	   (throw 'missing-4f8g-pdb))

       (let ((status (mutate-base imol "A" 19 "" "DC")))
	 (if (not (= status 1))
	   (begin
	     (format #t "failed to mutate 4f8g.pdb~%")
	     (throw 'fail-mutate-1)))

	 (let ((atom-list (residue-info imol "A" 19 "")))
	   (if (not (list? atom-list))
	       (throw 'bad-atom-list-1))
	   
	   (if (not (> (length atom-list) 10))
	       (throw 'bad-atom-list-need-more-atoms-1))

	   (let ((atoms (map residue-atom->atom-name atom-list)))
	     (if (member? " C7 " atoms)
		 (throw 'C7-still-present))

	     (write-pdb-file imol "test-1.pdb")

	     #t))))))


(greg-testcase "C7 is added on mutation to a DC" #t
   (lambda ()

     (let ((imol (greg-pdb "4f8g.pdb")))
       ;; now try mutate 10 to a T
		 
       (let ((status (mutate-base imol "A" 10 "" "DT")))
	 (if (not (= status 1))
	   (begin
	     (format #t "failed to mutate back 4f8g.pdb~%")
	     (throw 'fail-mutate-2)))
	 (let ((atom-list (residue-info imol "A" 10 "")))
	   (if (not (list? atom-list))
	       (throw 'bad-atom-list-2))

	   (let ((atoms (map residue-atom->atom-name atom-list)))
	     (print-var atoms)
	     (if (not (member? " C7 " atoms))
		 (throw 'C7-not-present-in-T)))))

       (write-pdb-file imol "test-2.pdb")

       #t)))




;; Restore this when the delete-residue-with-full-spec makes it to trunk.
;; 

; (greg-testcase "Deleting (non-existing) Alt conf and Go To Atom [JED]" #t
;    (lambda ()

;      ;; alt conf "A" does not exist in this residue:
;      (delete-residue-with-altconf imol-rnase "A" 88 "" "A")
;      ;; to activate the bug, we need to search over all atoms
;      (active-residue) ; crash
;      #t))


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

	  (format #t "   high-values: ~s  low values: ~s~%~!" 
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
	   
	      (format #t "   diff-high-values: ~s  diff-low-values: ~s~%~!" 
		      diff-high-values diff-low-values)

	      (if (not (< (apply + diff-high-values) 0.06)) 
		  (begin
		    (format #t "Bad diff high values: value: ~s target: ~s~%"
			    (apply + diff-high-values)  0.06)
		    (throw 'fail)))
	      
	      (if (not (< (apply + diff-low-values) -5))
		  (begin
		    (format #t "Bad diff low values ~s ~%" (apply + diff-low-values))
		    (throw 'fail)))
	      
	      #t))))))))

(greg-testcase "Skeletonize a map" #t
   (lambda ()

     (let ((imol (read-pdb rnase-pdb))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))

       (skeletonize-map  1 imol-map)
       (skeletonize-map  0 imol-map)
       (skeletonize-map -1 -1)
       (skeletonize-map 0 0)
       (close-molecule imol)
       (close-molecule imol-map)
       #t ;; don't crash
       )))


(greg-testcase "Simple Averaged maps" #t 
   (lambda ()

     (let ((imol-map-1 (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0))
	   (novalue-1 (set-map-sampling-rate 2.5))
	   (imol-map-2 (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0))
	   (novalue-2 (set-map-sampling-rate 1.5))) ;; reset it

       (let ((new-map (average-map (list (list imol-map-1 1)
					 (list imol-map-2 1)))))
	 (if (not (valid-map-molecule? new-map))
	     (begin 
	       (format #t " average map fail~%")
	       (throw 'fail)))

	 ;; the difference map should be nearly flat (0.0)
	 ;; 
	 (let ((diff-map (difference-map imol-map-1 new-map 1.0)))
	   
	   (let ((rms-1 (map-sigma  new-map))
		 (rms-2 (map-sigma diff-map)))

	     (format #t "  INFO:: map sigmas: normal ~s and diff-map: ~s~%" 
		     rms-1 rms-2)

	     (if (not (> rms-1 0.3))
		 (begin 
		   (format #t " map sigma 1 fail~%")
		   (throw 'fail)))

	     (if (not (< rms-2 0.003))
		 (begin 
		   (format #t " map sigma for diff average map fail~%")
		   (throw 'fail)))

	     #t))))))




(greg-testcase "Make a glycosidic linkage" #t 
   (lambda ()

     (let* ((carbo "multi-carbo-coot-3.pdb")
	    (imol (greg-pdb carbo)))

       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "file not found: ~s~%" carbo)
	     #f)
	   
	   (let ((atom-1 (get-atom imol "A" 1 "" " O4 "))
		 (atom-2 (get-atom imol "A" 2 "" " C1 ")))
	     
	     (format #t "   bond-length: ~s: ~%"
		     (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))

	     (let ((s (dragged-refinement-steps-per-frame)))
	       (set-dragged-refinement-steps-per-frame 300)
	       (with-auto-accept
		(regularize-zone imol "A" 1 2 ""))
	       (set-dragged-refinement-steps-per-frame s))

	     (let ((atom-1 (get-atom imol "A" 1 "" " O4 "))
		   (atom-2 (get-atom imol "A" 2 "" " C1 ")))

	       (format #t "   bond-length: ~s: ~%"
		       (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))
	       
	       (bond-length-within-tolerance? atom-1 atom-2 1.439 0.04)))))))


(greg-testcase "Refine an NAG-ASN Link" #t
   (lambda ()
     
     (let ((imol (greg-pdb "pdb2qc1-sans-cho.pdb"))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))

       (if (not (valid-model-molecule? imol))
	   (throw 'failed-to-find-pdb2qc1-sans-cho))
       
       (let ((status (add-linked-residue imol "B" 141 "" "NAG" "NAG-ASN" 3000)))

	 (with-auto-accept (refine-residues imol (list (list "B" 141 "") (list "B" 464 ""))))

	 (let ((atom-1 (get-atom imol "B" 141 "" " ND2"))
	       (atom-2 (get-atom imol "B" 464 "" " C1 ")))

	   (bond-length-within-tolerance? atom-1 atom-2 1.43 0.2))))))



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
       
       (let ((atom-1 (get-atom imol "A" 1 "" "HG11"))
	     (atom-2 (get-atom imol "A" 1 "" " CG1")))

	 (if (bond-length-within-tolerance? atom-1 atom-2 0.96 0.02)
	     #t
	     (begin
	       (if (and (list? atom-1)
			(list? atom-2))
		   (format "   flying hydrogen failure, bond length ~s, should be 0.96~%"
			   (bond-length-from-atoms atom-1 atom-2))
		   (format "   flying hydrogen failure, atoms: ~s ~s~%" atom-1 atom-2))
	       #f))))))




(greg-testcase "Test for regularization and mangling of hydrogen names from a PDB v 3.0" #t
   (lambda ()

     ;; Note that it seems to me that the bonds are not within
     ;; tolerance before the regularization.
     ;; 
     ;; Also note that the D atom has been removed because now we have
     ;; a test for the atom names matching (hydrogens are not checked
     ;; currently, but this atom does not appear to be a hydrogen
     ;; (element is " D").

     (let ((imol (greg-pdb "3ins-6B-3.0-no-peptide-D.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "Bad read of greg test pdb: 3ins-6B-3.0-no-peptide-D.pdb~%")
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
	     (let ((atom-1 (get-atom imol "B" 6 "" (car pair)))
		   (atom-2 (get-atom imol "B" 6 "" (cdr pair))))
	       (if (bond-length-within-tolerance? atom-1 atom-2 0.96 0.02)
		   #t
		   (begin
		     (let ((bl (bond-length-from-atoms atom-1 atom-2)))
		       (format #t "  Oops! bond length not within tolerance: ~s~%" bl)
		       (format #t "  Atom names:~%   ~s~%   ~s~%" atom-1 atom-2)
		       #f)))))
	   atom-pairs))))))


(greg-testcase "correct matching dictionary names from test name" #t
   (lambda ()

     ;; new dictionary
     (let ((ls (matching-compound-names-from-dictionary "gua" 0)))
       (if 

	;; This is too fragile
	;;
	;; (or (= (length ls) 153)    ;; new dictionary
	;;     (= (length ls)  63)    ;; 6.1.3 dictionary
	;;     (= (length ls)  61))   ;; 6.0.2 dictionary

	(> (length ls) 60)

	   #t ;; good

	   (begin
	     (format #t "   found ~s matching names: ~s~%" (length ls) ls)
	     #f)))))



(greg-testcase "update monomer restraints" #t 
   (lambda () 

     ;; this test needs a refinment map

     (let ((atom-pair ; (list " CB " " CG ")) ;; round the wrong way in the new dictionary
	                (list " CG " " CB "))
	   (m (monomer-restraints "TYR")))

       (if (not m)
	   (begin 
	     (format #t "   update bond restraints - no momomer restraints~%")
	     (throw 'fail)))

       ;; let's delete both ways
       ;;
       (let* ((n-t (strip-bond-from-restraints atom-pair m))
	      (n   (strip-bond-from-restraints (reverse atom-pair) n-t)))
	 (set-monomer-restraints "TYR" n)
	 
	 (let ((imol (new-molecule-by-atom-selection imol-rnase "//A/30")))
	   
	   (with-auto-accept
	    (refine-zone imol "A" 30 30 ""))
	   
	   (let ((atom-1 (get-atom imol "A" 30 "" " CB "))
		 (atom-2 (get-atom imol "A" 30 "" " CG ")))

	     (format #t "   Bond-length: ~s: ~%"
		     (bond-length (list-ref atom-1 2) (list-ref atom-2 2)))

             ;; with the new refinement system, I am not sure that I expect the
             ;; atoms to interact only with NBC term.
             ;; So remove the test.
	     (if #t
		 (begin 
		   (set-monomer-restraints "TYR" m)

		   (with-auto-accept
		    (refine-zone imol "A" 30 30 ""))
		 
		   (let ((atom-1 (get-atom imol "A" 30 "" " CB "))
			 (atom-2 (get-atom imol "A" 30 "" " CG ")))

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


(greg-testcase "Write mmCIF restraints correctly" #t
   (lambda ()

     ;; not a very good test, but better than nothing and it would
     ;; fail on cif writing prior to today 20100223.

     (if (file-exists? "coot-test-ala.cif")
	 (delete-file "coot-test-ala.cif"))

     (write-restraints-cif-dictionary "ALA" "coot-test-ala.cif")
     (if (not (file-exists? "coot-test-ala.cif"))
	 #f
	 (call-with-input-file "coot-test-ala.cif"
	   (lambda (port)
	     (let loop ((line (read-line port))
			(n-found 0))
	       (cond
		((= n-found 2) #t)
		((eof-object? line) #f) ; failed if we get here
		((string=? line "data_comp_list")
		 (loop (read-line port) 1))
		((string=? line "data_comp_ALA")
		 (loop (read-line port) (+ n-found 1)))
		(else
		 (loop (read-line port) n-found)))))))))



(greg-testcase "Refinement OK with zero bond esd" #t
   (lambda ()
     
     ;; return the restraints as passed (in r-list) except for a bond
     ;; restraint atom-1 to atom-2 set the esd to 0.0 (new-dist is
     ;; passed but not useful).
     (define (zero-bond-restraint r-list atom-1 atom-2 new-dist)
     
       (let loop ((r-list r-list))
	 
	 (cond
	  ((null? r-list) '())
	  ((string=? "_chem_comp_bond" (car (car r-list)))

	   (let ((v 
		  (let ((bond-list (cdr (car r-list))))

		    (let loop-2 ((bond-list bond-list))
		      
		      (cond 
		       ((null? bond-list) '())
		       ((and (string=? (car (car bond-list)) atom-1)
			     (string=? (car (cdr (car bond-list))) atom-2))
			(cons 
			 (list atom-1 atom-2 new-dist 0.0)
			 (loop-2 (cdr bond-list))))
		       (else
			(cons (car bond-list) (loop-2 (cdr bond-list)))))))))

	     (cons (cons "_chem_comp_bond" v) (loop (cdr r-list)))))
	   
	  (else
	   (cons (car r-list) (loop (cdr r-list)))))))
     
     ;; main line
     ;; 
     (let ((imol (greg-pdb "monomer-ACT.pdb")))
       (read-cif-dictionary (append-dir-file greg-data-dir "libcheck_ACT.cif"))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "   bad molecule from ACT from greg data dir~%")
	     (throw 'fail)))
       
       (let* ((r (monomer-restraints "ACT")))
	 
	 (case (list? r)
	   ((#f)
	    (format #t "   ACT restraints are #f~%")
	    (throw 'fail)))

	 (let ((r-2 (zero-bond-restraint r " O  " " C  " 1.25)))
	   
	   ;; (format #t "r: ~s~%=======~%" r)
	   ;; (format #t "r-2: ~s~%" r-2)
	   
	   (case r-2
	     (('())
	      ((format #t "   null modified restraints~%")
	      (throw 'fail))))
	   
	   (let ((mr-set (set-monomer-restraints "ACT" r-2)))
	     
	     (case mr-set
	       ((#f)
		(format #t "   set restraints fail~%")
		(throw 'fail)))
	     
	     (let ((r-3 (monomer-restraints "ACT")))
	       ;; (format #t "r-3: ~s~%" r-3)
	       
	       (case r-3
		 ((#f)
		  (format #t "   get modified restraints fail~%")
		  (throw 'fail)))
	       
	       (with-auto-accept
		(refine-zone imol "A" 1 1 ""))
	       
	       #t))))))) ;; it didn't crash

;; enable when we are in 07
;; 
;;; simple-minded Hydrogen-based filtered (not topological).
;;; 
;(greg-testcase "filter out multi-hydrogen chiral centres" #t
;   (lambda ()

;     (let ((imol (greg-pdb "chiral-test-1.pdb")))
;       (if (not (valid-model-molecule? imol))
;	   (throw "chiral-test-1.pdb not found"))
;       (read-cif-dictionary (append-dir-file greg-data-dir "chiral-test-1-dict.cif"))
;       (let ((r (monomer-restraints "DRG")))
;	 (let ((chirals (assoc-ref r "_chem_comp_chir")))
;	   (= (length chirals) 2))))))



(greg-testcase "Change Chain IDs and Chain Sorting" #t 
   (lambda () 

     (define (chains-in-order? chain-list)
       
       (let loop ((ls chain-list)
		  (ref-chain-id ""))
	 (cond
	  ((null? ls) #t)
	  ((string<? (car ls) ref-chain-id)
	   (begin
	     (format #t "ERROR:: ~s was less than ~s in ~s ~%" 
		     (car ls) ref-chain-id chain-list)
	     #f))
	  (else 
	   (loop (cdr ls) (car ls))))))
	    
     ;; 
     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (change-chain-id imol "A" "D" 0  0  0)
       (change-chain-id imol "B" "N" 1 30 38)
       (change-chain-id imol "B" "L" 1 40 49)
       (change-chain-id imol "B" "J" 1 50 59)
       (change-chain-id imol "B" "Z" 1 20 28)
       (change-chain-id imol "B" "G" 1 60 70)
       (change-chain-id imol "B" "F" 1 70 80)
       (change-chain-id imol "B" "E" 1 80 90)  

       (sort-chains imol)

       (let ((c (chain-ids imol)))
	 (chains-in-order? c)))))


(greg-testcase "Chain-ids in links change also on change chain id" #t 
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))

       (let ((spec-1 (list "B" 42 "" " N  " ""))
	     (spec-2 (list "B" 62 "" " O  " "")))

	 (make-link imol spec-1 spec-2 "test-link" 2.2)

	 (let ((li-1 (link-info imol)))

	   (change-chain-id imol "B" "C" 0 0 0)

	   (let ((li-2 (link-info imol)))

	     ;; li-1 should not contain C and should contain B
	     ;; li-2 should not contain B and should contain C

	     ;; (format #t "li-1: ~s~%" li-1)
	     ;; (format #t "li-2: ~s~%" li-2)

	     (let ((ch-B-1 (list-ref (list-ref (list-ref li-1 0) 1) 1)) ;; before
		   (ch-B-2 (list-ref (list-ref (list-ref li-1 0) 2) 1))
		   (ch-A-1 (list-ref (list-ref (list-ref li-2 0) 1) 1)) ;; after
		   (ch-A-2 (list-ref (list-ref (list-ref li-2 0) 2) 1)))

	       ;(format #t "ch-B-1: ~s~%" ch-B-1)
	       ;(format #t "ch-B-2: ~s~%" ch-B-2)
	       ;(format #t "ch-A-1: ~s~%" ch-A-1)
	       ;(format #t "ch-A-2: ~s~%" ch-A-2)

	       (all-true? (list
			   (string=? ch-B-1 "B")
			   (string=? ch-B-2 "B")
			   (string=? ch-A-1 "C")
			   (string=? ch-A-2 "C"))))))))))



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

	   (format #t "   distances: ~s~%" (map atom-distance reference-res replaced-res))
	   (all-true? (map (lambda (d) (> d 20)) (map atom-distance reference-res replaced-res))))))))



(greg-testcase "Residues in Region of Residue" #t 
   (lambda ()

     (all-true?
      (map (lambda (dist n-neighbours)
	     
	     (let ((rs (residues-near-residue imol-rnase (list "A" 40 "") dist)))
	       (if (not (= (length rs) n-neighbours))
		   (begin
		     (format #t "wrong number of neighbours ~s ~s~%" (length rs) rs)
		     #f)
		   (begin
		     (format #t "   found ~s neighbours ~s~%" (length rs) rs)
		     #t))))
	   (list 4 0) (list 6 0)))))


(greg-testcase "Residues in region of a point" #t 
   (lambda () 

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (let* ((pt (list 43.838   0.734  13.811)) ; CA 47 A
	      (residues (residues-near-position imol pt 2)))
	 (if (not (= (length residues) 1))
	     (begin
	       (format #t "  Fail, got residues: ~s~%" residues)
	       (throw 'fail))
	     (if (not (equal? (car residues)
			      (list #t "A" 47 "")))
		 (begin
		   (format #t "  Fail 2, got residues: ~s~%" residues)
		   (throw 'fail))
		 
		 (let ((residues-2 (residues-near-position imol pt 4)))
		   (if (not (= (length residues-2) 3)) ;; its neighbours too.
		       (begin
			 (format #t "  Fail 3, got residues-2: ~s~%" residues-2)
			 (throw 'fail))
		       #t))))))))



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
	 (set-residue-to-rotamer-number imol-rnase chain-id resno "" "" 1)
	 
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


(greg-testcase "Rotamer names and scores are correct" #t 
   (lambda () 

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (with-no-backups 
	imol 
	
	(let* ((residue-attributes (list "A" 28 ""))
	       (alt-conf "")
	       (residue-attributes-with-alt-conf (append residue-attributes (list alt-conf)))
	       (results 
	  
		(map
		 (lambda (rotamer-number correct-name correct-prob)
		   (apply set-residue-to-rotamer-number imol 
			  (append residue-attributes (list "" rotamer-number)))
		   (let ((rotamer-name (apply get-rotamer-name imol residue-attributes))
			 (rotamer-prob (apply rotamer-score imol residue-attributes-with-alt-conf)))
		     (format #t "   Rotamer ~s : ~s ~s ~%" 
			     rotamer-number rotamer-name rotamer-prob)
		     (if (not (close-float? correct-prob rotamer-prob))
			 (begin
			   (format #t "fail on rotamer probability: result: ~s should be: ~s ~%"
				   rotamer-prob correct-prob)
			   (throw 'fail)))
		     (string=? rotamer-name correct-name)))
		 (range 0 5)
		 ;; the spaces at the end of the name are in the Complete_rotamer_lib.csv
		 ;; and hence in richardson-rotamers.cc.  C'est la vie.
		 (list "m-85" "t80" "p90" "m -30 " "m -30 ")
		 (list 100 90.16684 50.707787 21.423154 21.423154))))
	  (all-true? results))))))


(greg-testcase "Align and mutate a model with deletions" #t 
   (lambda ()

     (define (residue-in-molecule? imol chain-id resno ins-code)
       (let ((r (residue-info imol chain-id resno ins-code)))
	 (if r #t #f)))

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
	   (let ((renumber? 1))

	     (set-alignment-gap-and-space-penalty -3.0 -0.5)
	     (align-and-mutate imol "A" rnase-seq-string renumber?)
	     (write-pdb-file imol "mutated.pdb")

	     (let ((results
		    (map
		     (lambda (residue-info)
		       (let ((residue-spec (cdr residue-info))
			     (expected-status (car residue-info)))
			 (format #t "    ::::: ~s ~s ~s~%" residue-spec
				 (apply residue-in-molecule? residue-spec)
				 expected-status)
			 (eq?
			  (apply residue-in-molecule? residue-spec)
			  expected-status)))
		     (list
		      (list #f imol "A"  1 "")
		      (list #t imol "A"  4 "")
		      (list #t imol "A" 57 "")
		      (list #f imol "A" 60 "")
		      (list #f imol "A" 61 "")
		      (list #t imol "A" 92 "")
		      (list #f imol "A" 94 "")
		      ))))

	       (format #t "results: ~s~%" results)

	       (all-true? results)))))))


(greg-testcase "renumbered residues should be in seqnum order" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (renumber-residue-range imol "A" 10 20 100)
       (let ((iser (seqnum-from-serial-number imol "A" 90)))
	 (if (not (= iser 118)) 
	     (throw 'fail))))

     ;; and again this time adding residue to the front of the molecule
     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (renumber-residue-range imol "A" 10 20 -100)
       (let ((iser (seqnum-from-serial-number imol "A" 0)))
	 (if (not (= iser -90))
	     (throw 'fail))))

     #t))



;; 
(greg-testcase "Autofit Rotamer on Residues with Insertion codes" #t 
   (lambda ()

     ;; what's the plan?  
     
     ;; we need to check that H52 LEU and H53 GLY do not move and H52A does move

   (define (centre-atoms mat)
     (map (lambda (ls)
	    (/ (apply + ls) (length ls)))
	  (transpose-mat (map (lambda (x) (list-ref x 2)) mat))))
     
   ;; 
   (let* ((imol (greg-pdb "pdb3hfl.ent"))
	  (mtz-file-name (append-dir-file greg-data-dir "3hfl_sigmaa.mtz"))
	  (imol-map (make-and-draw-map mtz-file-name
				       "2FOFCWT" "PH2FOFCWT" "" 0 0)))

     (if (not (valid-model-molecule? imol))
	 (begin
	   (format #t " ERROR:: pdb3hfl.ent nto found~%")
	   (throw 'fail)))

     (if (not (valid-map-molecule? imol-map))
	 (begin
	   (format #t "   no map from 3hfl_sigmaa.mtz~%")
	   (throw 'fail))

	 (let ((leu-atoms-1 (residue-info imol "H" 52 ""))
	       (leu-resname (residue-name imol "H" 52 ""))
	       (gly-atoms-1 (residue-info imol "H" 53 ""))
	       (gly-resname (residue-name imol "H" 53 ""))
	       (pro-atoms-1 (residue-info imol "H" 52 "A"))
	       (pro-resname (residue-name imol "H" 52 "A")))

	   ;; First check that the residue names are correct
	   (if (not
		(and 
		 (string=? leu-resname "LEU")
		 (string=? gly-resname "GLY")
		 (string=? pro-resname "PRO")))
	       (begin
		 (format #t "  failure of residues names: ~s ~s ~s~%"
			 leu-resname gly-resname pro-resname)
		 (throw 'fail)))
	   

	   ;; OK, what are the centre points of these residues?
	   (let ((leu-centre-1 (centre-atoms leu-atoms-1))
		 (gly-centre-1 (centre-atoms gly-atoms-1))
		 (pro-centre-1 (centre-atoms pro-atoms-1)))

	     (auto-fit-best-rotamer 52 "" "A" "H" imol imol-map 0 0)

	     ;; OK, what are the centre points of these residues?
	     (let ((leu-atoms-2 (residue-info imol "H" 52 ""))
		   (gly-atoms-2 (residue-info imol "H" 53 ""))
		   (pro-atoms-2 (residue-info imol "H" 52 "A")))
	       
	       (let ((leu-centre-2 (centre-atoms leu-atoms-2))
		     (gly-centre-2 (centre-atoms gly-atoms-2))
		     (pro-centre-2 (centre-atoms pro-atoms-2)))
		 
		 (let ((d-leu (pos-difference leu-centre-1 leu-centre-2))
		       (d-gly (pos-difference gly-centre-1 gly-centre-2))
		       (d-pro (pos-difference pro-centre-1 pro-centre-2)))
		   
		   (if (not (close-float? d-leu 0.0))
		       (begin 
			 (format #t "   Failure: LEU 52 moved~%")
			 (throw 'fail)))

		   (if (not (close-float? d-gly 0.0))
		       (begin 
			 (format #t "   Failure: GLY 53 moved~%")
			 (throw 'fail)))

		   (if (< d-pro 0.05)
		       (begin 
			 (format #t "   Failure: PRO 52A not moved enough: ~s~%" d-pro)
			 (map (lambda (a) (format #t "   PRO-atoms 1: ~s~%" a)) pro-atoms-1)
			 (map (lambda (a) (format #t "   PRO-atoms 2: ~s~%" a)) pro-atoms-2)
			 (throw 'fail)))
		   
		   ;; rotamer 4 is out of range.
		   (set-residue-to-rotamer-number imol "H" 52 "A" "" 4) ;; crash
		   #t
		   )))))))))


;; new tests

(greg-testcase "RNA base has correct residue type after mutation" #t 
   (lambda ()

     (define (test-vs rna-mol base-name)

       (let ((previous-name (residue-name rna-mol "A" 2 ""))
	     (success (mutate-base rna-mol "A" 2 "" base-name)))
	 (if (not (= success 1))
	     (begin
	       (format #t "  mutation fail!~%")
	       (throw 'fail)))
	 (let ((rn (residue-name rna-mol "A" 2 "")))
	   (format #t "  mutated base to type ~s - was ~s ~%" rn previous-name)
	   (string=? rn base-name))))


     ;; main line
     ;; 
     (let ((rna-mol (ideal-nucleic-acid "RNA" "A" 0 "GACUCUAG")))
       (let ((res-1 (test-vs rna-mol "C")))
	 (if (not res-1)
	     (begin
	       (format #t "  incorrect base! (default names) ~%")
	       (throw 'fail))))
	   
       (set-convert-to-v2-atom-names 1)

       (let ((rna-mol-old-names (ideal-nucleic-acid "RNA" "A" 0 "GACUCUAG")))
	 (let ((res-2 (test-vs rna-mol-old-names "Cr")))

	   ;; back to normal
	   (set-convert-to-v2-atom-names 0)
	   
	   (if (not res-2)
	       (begin
		 (format #t "  incorrect base! (old names) ~%")
		 (throw 'fail))))
       
	 #t))))


;; Wolfram Tempel bug
;; 
(greg-testcase "resname from serial number doesnt crash on silly input" #t 
   (lambda ()
     (resname-from-serial-number 0 "DNA" 1)
     (resname-from-serial-number 0 "A"  -1)
     (resname-from-serial-number imol-rnase "A"  -1)
     (resname-from-serial-number imol-rnase "A"  65)
     #t))



(greg-testcase "DNA bases are the correct residue type after mutation" #t  
   (lambda ()


     (define (correct-base-type? rna-mol target-base-type)
       (let ((success (mutate-base rna-mol "A" 2 "" target-base-type)))
	 (if (not (= success 1))
	     (begin 
	       (format #t "  DNA base mutation fail!~%")
	       (throw 'fail)))
	 (let ((rn (residue-name rna-mol "A" 2 "")))
	   (format #t "  mutated base to type ~s~%" rn)
	   (string=? rn target-base-type))))

     ;; main line
     (let ((rna-mol (ideal-nucleic-acid "DNA" "A" 0 "GACTCTAG")))

       (if (not (all-true?
		 (map (lambda (base)
			(correct-base-type? rna-mol base))
		      (list "DC" "DG" "DA" "DT"))))

	   (begin
	     (format #t "Fail in DNA~%")
	     (throw 'fail))))


     (set-convert-to-v2-atom-names 1)
     (let ((rna-mol (ideal-nucleic-acid "DNA" "A" 0 "GACTCTAG")))
       
       (if (not (all-true?
		 (map (lambda (base)
			(correct-base-type? rna-mol base))
		      (list "Cd" "Gd" "Ad" "Td"))))
	   
	   (begin
	     (set-convert-to-v2-atom-names 0)
	     (format #t "Fail in DNA~%")
	     (throw 'fail))))

     (set-convert-to-v2-atom-names 0)
     #t))



(greg-testcase "SegIDs are correct after mutate" #t
   (lambda () 

     ;; main line
     (let ((imol (copy-molecule imol-rnase)))
       
       (if (not (valid-model-molecule? imol))
	   (throw 'fail))

       (let ((atoms (residue-info imol "A" 32 "")))
	 
	 (if (not (atoms-have-correct-seg-id? atoms ""))
	     (begin
	       (format #t "wrong seg-id ~s should be ~s~%" atoms "")
	       (throw 'fail))))

       ;; now convert that residue to segid "A"
       ;; 
       (let ((attribs (map (lambda (atom)
			     (list imol "A" 32 ""
				   (car (car atom))
				   (car (cdr (car atom)))
				   "segid" "A"))
			   (residue-info imol "A" 32 ""))))
	 
	 (set-atom-attributes attribs))

       (let ((atoms (residue-info imol "A" 32 "")))
	 (if (not (atoms-have-correct-seg-id? atoms "A"))
	     (begin
	       (format #t "wrong seg-id ~s should be ~s~%" atoms "A")
	       (throw 'fail))))

       ;; now let's do the mutation
       (mutate imol "A" 32 "" "LYS")

       (let ((rn (residue-name imol "A" 32 "")))
	 (if (not (string=? rn "LYS"))
	     (begin 
	       (format #t "  Wrong residue name after mutate ~s~%" rn)
	       (throw 'fail))))

       (let ((atoms (residue-info imol "A" 32 "")))
	 (if (not (atoms-have-correct-seg-id? atoms "A"))
	     (begin
	       (format #t "wrong seg-id ~s should be ~s~%" atoms "A")
	       (throw 'fail))))

       #t)))


(greg-testcase "TER on water chain is removed on adding a water by hand" #t 
   (lambda ()

     (let ((imol (greg-pdb "some-waters-with-ter.pdb")))
       
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "bad read of some-waters-with-ter.pdb~%")
	     (throw 'fail)))
       
       (set-rotation-centre 3 4 5)
       (set-pointer-atom-molecule imol)
       (place-typed-atom-at-pointer "Water")

       ;; OK let's write out that molecule now as a PDB file and look
       ;; to see if the PDB file contains a TER record - it should
       ;; not.
       ;; 
       (let ((opdb  "tmp-with-new-water.pdb"))
	 (write-pdb-file imol opdb)
	 (call-with-input-file opdb
	   (lambda (port)

	     (let loop ((obj (read-line port)))
	       (cond
		((eof-object? obj) #t) ; test did not fail before now,
				       ; so it passes.
		((string-match "TER" obj)
		 (format #t "   TER card found: ~%    ~s~%" obj)
		 #f ;; fail!
		 )
		(else (loop (read-line port)))))))))))


(greg-testcase "TER on water chain is removed on adding waters automatically" #t 
   (lambda ()
     
     (let ((imol-model (greg-pdb "tm+some-waters.pdb")))
       
       (if (not (valid-model-molecule? imol-model))
	   (begin
	     (format #t "tm+some-waters.pdb not found~%")
	     (throw 'fail)))

       (find-waters imol-rnase-map imol-model 0 2.0 0)
       (write-pdb-file imol-model "auto-waters.pdb")

       ;; OK let's write out that molecule now as a PDB file and look
       ;; to see if the PDB file contains a TER record - it should
       ;; not.
       ;; 
       (let ((opdb  "tmp-with-new-water.pdb"))
	 (write-pdb-file imol-model opdb)
	 (call-with-input-file opdb
	   (lambda (port)

	     (let loop ((obj (read-line port)))
	       (cond
		((eof-object? obj) #t) ; test did not fail before now,
				       ; so it passes.
		((string-match "TER" obj)
		 (format #t "   TER card found: ~%    ~s~%" obj)
		 #f ;; fail!
		 )
		(else (loop (read-line port)))))))))))


(greg-testcase "Adding atoms to Many-Chained Molecule" #t 
   (lambda () 

     (let ((imol (read-pdb rnase-pdb)))
       (set-pointer-atom-molecule imol)
       (for-each 
	(lambda (i)
	  (place-typed-atom-at-pointer "Mg"))
	(range 100))
       #t))) ;; doesn't crash :)


(greg-testcase "Arrange waters round protein" #t
   (lambda ()

     (let ((imol (greg-pdb "water-test-no-cell.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "ERROR:: water-test-no-cell not found~%")
	     (throw 'fail)))

       (let ((status (move-waters-to-around-protein imol)))
	 (if (not (= status 0))
	   (begin
	     (format #t "ERROR:: failure with water-test-no-cell~%")
	     (throw 'fail)))))

     (let ((imol (greg-pdb "pathological-water-test.pdb")))
       (if (not (valid-model-molecule? imol))
	   (begin
	     (format #t "ERROR:: pathological-waters pdb not found~%")
	     (throw 'fail)))

       (let ((status (move-waters-to-around-protein imol)))
	 (write-pdb-file imol "waters-moved-failure.pdb")
	 (if (not (= status 181))
	   (begin
	     (format #t "ERROR:: failure with pathological-waters moved ~s~%" status)
	     (throw 'fail)))

	 (let ((v (max-water-distance imol)))

	   (if (not (< v 5.0))
	       (begin
		 (format #t "ERROR:: failure to move waters close ~s~%" v)
		 (throw 'fail))

	       #t))))))

	     

(greg-testcase "Correct Segid After Add Terminal Residue" #t 
   (lambda () 

     (let ((imol (copy-molecule imol-rnase))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))

       ;; now convert that residue to segid "A"
       ;; 
       (let ((attribs (map (lambda (atom)
			     (list imol "A" 93 ""
				   (car (car atom))
				   (car (cdr (car atom)))
				   "segid" "A"))
			   (residue-info imol "A" 93 ""))))
	 (set-atom-attributes attribs))

       (add-terminal-residue imol "A" 93 "ALA" 1)

       (let ((atoms (residue-info imol "A" 94 "")))
	 (if (not (atoms-have-correct-seg-id? atoms "A"))
	     (begin
	       (format #t "wrong seg-id ~s should be ~s~%" atoms "A")
	       (throw 'fail))))

       #t)))


(greg-testcase "Correct Segid after NCS residue range copy" #t 
   (lambda () 

     (define (convert-residue-seg-ids imol chain-id resno seg-id)
       (let ((attribs (map (lambda (atom)
			     (list imol chain-id resno ""
				   (car (car atom))
				   (car (cdr (car atom)))
				   "segid" seg-id))
			   (residue-info imol chain-id resno ""))))
	 (set-atom-attributes attribs)))
       

     ;; main line
     (let ((imol (copy-molecule imol-rnase)))
       
       ;; convert a residue range
       (for-each (lambda (resno)
		  (convert-residue-seg-ids imol "A" resno "X"))
		(range 20 30))

       ;; NCS copy
       (make-ncs-ghosts-maybe imol)
       (copy-residue-range-from-ncs-master-to-others imol "A" 20 30)

       ;; does B have segid X for residues 20-30?
       (for-each (lambda (resno)
		  (let ((atoms (residue-info imol "A" resno "")))
		    (if (not (atoms-have-correct-seg-id? atoms "X"))
			(begin
			  (format #t "wrong seg-id ~s should be ~s~%" atoms "X")
			  (throw 'fail)))))
		(range 20 30))

       #t)))


(greg-testcase "Merge Water Chains" #t
   (lambda ()

     (define (create-water-chain imol from-chain-id chain-id n-waters offset prev-offset)
       (for-each (lambda (n)
		   (place-typed-atom-at-pointer "Water")
		   ;; move the centre of the screen
		   (let ((rc (rotation-centre)))
		     (apply set-rotation-centre 
			    (+ (car rc) 2)
			    (cdr rc))))
		 (range n-waters))
       (change-chain-id imol from-chain-id chain-id 1 
			(+ 1 prev-offset) 
			(+ n-waters prev-offset))
       (renumber-residue-range imol chain-id (+ 1 prev-offset) (+ 5 prev-offset) offset))
     


     ;; main line
     (let ((imol (greg-pdb "tutorial-modern.pdb")))

       (with-no-backups imol
			(let ((imol imol))
			  (set-pointer-atom-molecule imol)
			  
			  ;; first create a few chains of waters
			  (create-water-chain imol "C" "D" 5 10 0)
			  (create-water-chain imol "D" "E" 5 5 10)
			  (create-water-chain imol "D" "F" 5 2 15)
			  
			  ;; OK, we have set up a molecule, now let's test it:
			  ;; 
			  (merge-solvent-chains imol)
			  ))
       
       ;; Test the result:
       ;; 
       (let ((nc (n-chains imol)))
	 (if (not (= nc 3))
	     (begin
	       (format #t "  wrong number of chains ~s~%" nc)
	       (throw 'fail))
	     ;; There should be 15 waters in the last chain
	     (let ((solvent-chain-id (chain-id imol 2)))
	       (let ((n-res (chain-n-residues solvent-chain-id imol)))
		 (if (not (= n-res 15))
		     (begin
		       (format #t "  wrong number of residues ~s~%" n-res)
		       (throw 'fail))
		     (let ((r1  (seqnum-from-serial-number imol "D" 0))
			   (r15 (seqnum-from-serial-number imol "D" 14)))
		       (if (not (= r1 1))
			   (begin
			     (format #t "  wrong residue number r1 ~s~%" r1)
			     (throw 'fail))
			   (if (not (= r15 15))
			       (begin
				 (format #t "  wrong residue number r15 ~s~%" r15)
				 (throw 'fail))
			       #t))))))))))) ;; passes test

;; tested merge by consolidation
;; 
(greg-testcase "Consolidated merge" #t
   (lambda () 

     (let ((imol (greg-pdb "pdb1hvv.ent"))
	   (imol-lig-1 (greg-pdb "monomer-ACT.pdb"))	   
	   (imol-lig-2 (greg-pdb "monomer-NPO.pdb")))


       (format #t "-------- starting chain list ----------- ~%")
       (print-var (chain-ids imol))


       (merge-molecules (list imol-lig-1 imol-lig-2) imol)

       (let ((imol-symm-copy (new-molecule-by-symop imol "-X,-X+Y,-Z+1/3" 0 0 0)))
	 
	 (if (not (valid-model-molecule? imol-symm-copy))
	     (throw 'symm-molecule-problem))
	 
	 (merge-molecules (list imol-symm-copy) imol)

         (let ((chain-list (chain-ids imol)))

           (print-var chain-list)

           (write-pdb-file imol "sym-merged.pdb")
   
           (equal? chain-list (list "A" "B" "C" "D" ""  ;; original
				    "E" "F"  ;; merged ligs
				    "G" "H" "I" "J" "K" ;; protein chain copies
				    )))))))


(greg-testcase "Test for good chain ids after a merge" #t
   (lambda ()

     (let ((imol (greg-pdb "tutorial-modern.pdb")))

       (change-chain-id imol "A" "AAA" 0 0 0)

       (let ((imol-new (new-molecule-by-atom-selection imol "//B/1-90")))

	 (change-chain-id imol-new "B" "B-chain" 0 0 0)

	 (merge-molecules (list imol-new) imol)

	 (let ((chids (chain-ids imol)))

	   (format #t "chain-ids: ~s~%" chids)

	   ;; should be '("AAA" "B" "B-chain")

	   (if (not (= (length chids) 3))
	       (throw 'fail))

	   (if (not (string=? (list-ref chids 2) "B-chain"))
	       (throw 'fail))

	   ;; now Wolfram Tempel test: multi-char chain matcher needs
	   ;; prefix, not new single letter

	   (change-chain-id imol-new "B-chain" "AAA" 0 0 0)

	   (merge-molecules (list imol-new) imol)

	   (let ((chids-2 (chain-ids imol)))

	     (format #t "--- chain-ids: ~s~%" chids-2)

	     (if (not (= (length chids-2) 4))
		 (throw 'fail))

	     (if (not (string=? (list-ref chids-2 3) "AAA2"))
		 (throw 'fail))

	     #t))))))


	
(greg-testcase "LSQ by atom" #t
   (lambda ()

     (define (make-spec-ref atom-name)
       (list "A" 35 "" atom-name ""))

     (define (make-spec-mov atom-name)
       (list "B" 35 "" atom-name ""))

     (clear-lsq-matches)
     (let ((imol-1 (copy-molecule imol-rnase))
	   (imol-2 (copy-molecule imol-rnase)))

       (let ((spec-refs (map make-spec-ref (list " CG2" " CG1" " CB " " CA ")))
	     (spec-movs (map make-spec-mov (list " CG2" " CG1" " CB " " CA "))))

	 (map add-lsq-atom-pair spec-refs spec-movs)

	 (let ((result (apply-lsq-matches imol-1 imol-2)))

	   (if (not result) 
	       (begin
		 (format #t "Bad match~%")
		 (throw 'fail)))
	   
	   (let* ((c-1 (get-atom imol-1 "A" 35 "" " C  "))
		  (c-2 (get-atom imol-2 "B" 35 "" " C  "))
		  (b (bond-length-from-atoms c-1 c-2)))
	     
	     (bond-length-within-tolerance? c-1 c-2 0.0 0.2)))))))


(greg-testcase "LSQing changes the space-group and cell to that of the reference molecule" #t
   (lambda ()

     (let ((imol-mov (greg-pdb "tutorial-modern.pdb"))
	   (imol-ref (greg-pdb "pdb1py3.ent")))

       (let ((sg-mov-orig (show-spacegroup imol-mov))
	     (sg-ref-orig (show-spacegroup imol-ref))
	     (cell-mov-orig (cell imol-mov))
	     (cell-ref-orig (cell imol-ref)))

       (clear-lsq-matches)
       (add-lsq-match 10 50 "A" 8 48 "B" 1)
       (let ((rtop (apply-lsq-matches imol-ref imol-mov)))
	 (let ((sg-mov-curr (show-spacegroup imol-mov))
	       (cell-mov-curr (cell imol-mov)))

	   (if (not (string=? sg-mov-curr sg-ref-orig))
	       (begin
		 (format #t "   fail on matching spacegroups: ~s and ~s~%" 
			 sg-mov-curr sg-ref-orig)
		 (throw 'fail)))

	   (let ((r (all-true? (map close-float? cell-ref-orig cell-mov-curr))))
	     (if (not r)
	       (begin
		 (format #t "   fail on matching cells: ~s and ~s~%" 
			 cell-ref-orig cell-mov-curr)
		 (throw 'fail)))

	     #t)))))))
		 
	   


(greg-testcase "set-residue-name sets the correct residue" #t
   (lambda () 
     (let ((imol (greg-pdb "tutorial-modern.pdb")))
       (set-residue-name imol "A" 37 "" "FRE")

       ;; There should be only one residue with that residue type and it
       ;; should be A37.
       ;; 
       (let* ((specs (fit-protein-make-specs imol 'all-chains))
	      (residues-with-name (get-residues-in-molecule-of-type imol "FRE")))
	 (if (not (= (length residues-with-name) 1))
	     #f 
	     (specs-match? (car residues-with-name) (list 0 "A" 37 "")))))))



(greg-testcase "fit-protein-make-specs makes all specs" #t 
   (lambda () 
     (let* ((imol (greg-pdb "tutorial-modern.pdb"))
	    (specs (fit-protein-make-specs imol 'all-chains)))
       (format #t "   specs: ~s ~s~%" (length specs) specs)
       (= (length specs) 189))))





(greg-testcase "Phosphate distance in pucker analysis is sane" #t 
   (lambda () 

     (let ((imol (greg-pdb "2goz-manip.pdb")))
       (let ((pi (pucker-info imol (list "B" 14 "") 0)))
	 (let ((phosphate-distance (car pi)))
	   (if (> phosphate-distance 0.2)
	       (begin
		 (format #t "  Bad phosphate distance on 14 ~s~%" phosphate-distance)
		 (throw 'fail)))))

       (let ((pi (pucker-info imol (list "B" 15 "") 0)))
	 (let ((phosphate-distance (car pi)))
	   (if (< phosphate-distance 2.0)
	       (begin
		 (format #t "  Bad phosphate distance on 15 ~s~%" phosphate-distance)
		 (throw 'fail)))))

       #t)))



(greg-testcase "Fix for Oliver Clarke fit by atom selection bug"  #t
   (lambda ()

     ;; Atoms of the A chain (up to residue 71) moves when the B chain is refined
     ;;

     (let ((imol-rnase (greg-pdb "tutorial-modern.pdb"))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))

       (set-imol-refinement-map imol-map)

       (with-auto-accept
	(rigid-body-refine-by-atom-selection imol-rnase "//B"))

       ;; did we get silly superposition?
       ;; test that the A-chain atom is close to where it should be

       (let ((atom (get-atom imol-rnase "A" 43 "" " CA " "")))
	 (if (not (list? atom))
	     (begin
	       (format #t "Failure to extract atom~%~!")
	       #f)
	     (let ((atom-pos (list-ref atom 2)))
	       (let ((bl (bond-length atom-pos (list 46.4 11.6 12.1))))
		 (format #t "bl: ~s~%" bl)
		 (if (> bl 1.0)
		     (begin
		       (format #t  "Fail: moved atom ~s~%" bl)
		       #f)
		     #t ; OK, it didn't move much
		     ))))))))

