
;; NOTE: hof.fcf has  _refln_A_calc and _refln_B_calc, not fcalc and phase.
;; This will not work.  Why is it here?
;; 
(define hof-fcf (append-dir-file greg-data-dir "hof.fcf"))
(define hof-res (append-dir-file greg-data-dir "HOF.RES"))

(define imol-hof-res #f)

(define insulin-fcf (append-dir-file greg-data-dir "insulin.fcf"))
(define insulin-res (append-dir-file greg-data-dir "insulin.res"))
(define hollander-ins (append-dir-file greg-data-dir "hollander.ins"))
(define imol-insulin-res -1) ; set later


(greg-testcase "Read small molecule .res file" #t
   (lambda ()
     (if (string? hof-res)
	 (if (file-exists? hof-res)
	     (let ((imol (read-pdb hof-res)))
	       (let ((success (valid-model-molecule? imol)))
		 (if success
		     (set! imol-hof-res imol))
		 success))
	     (begin
	       (format #t "~s does not exist - skipping test~%" hof-res)
	       (throw 'untested)))
	 (begin
	   (format #t "hof-res not defined - skipping test~%")
	   (throw 'untested)))))


(greg-testcase "Read hollander small molecule .res file" #t
   (lambda ()

     (if (not (file-exists? hollander-ins))
	 (begin
	   (format #t "~s does not exist - skipping test~%" hollander-ins)
	   (throw 'untested))
	 (let ((imol (read-pdb hollander-ins)))
	   (if (not (valid-model-molecule? imol))
	       (begin
		 (format #t "   fail: bad molecule for ~%" hollander-ins)
		 (throw 'fail))
	       (let ((spg (show-spacegroup imol)))
		 (if (not (string=? spg "I 41 2 2"))
		     (begin
		       (format #t "   fail: wrong spacegroup for ~s ~s ~%" 
			       hollander-ins spg)
		       (throw 'fail))
		     #t)))))))


(greg-testcase "read shelx insulin with fcf" #t 
   (lambda ()

     (let ((imol-insulin-res-local
	    (handle-read-draw-molecule-with-recentre insulin-res 1)))
       (if (not (valid-model-molecule? imol-insulin-res-local))
	   (begin
	     (format #t "   Bad insulin.res: ~s for ~s~%" 
		     insulin-res imol-insulin-res-local)
	     (throw 'fail)))

       (set! imol-insulin-res imol-insulin-res-local) ; used in water addition test
       (let ((imol (handle-shelx-fcf-file insulin-fcf)))
	 (if (not (valid-map-molecule? imol))
	     (begin
	       (format #t "   Bad read of ~s ~s~%" insulin-fcf imol)
	       (throw 'fail))
	     (begin
	       (let ((name (molecule-name imol))
		     (cif-name (string-append insulin-fcf ".cif SigmaA")))
		 (if (not (string=? name cif-name))
		     (begin
		       (format #t "   Bad name match ~s != ~s~%" name cif-name)
		       (throw 'fail))))

	       (if (not (string=? (show-spacegroup imol)
				  (show-spacegroup imol-insulin-res-local)))
		   (begin
		     (format #t "   Mismatch spacegroups ~s ~s~%" 
			     (show-spacegroup imol)
			     (show-spacegroup imol-insulin-res-local))
		     (throw 'fail))
		   
		   (if (not (string=? (show-spacegroup imol)
				      "I 21 3"))
		       (begin
			 (format #t "   Bad spacegroups ~s~%" 
				 (show-spacegroup imol))
			 (throw 'fail))

		       ;; good then:
		       (begin
			 (rotate-y-scene (rotate-n-frames 200) 0.1)
			 (close-molecule imol)
			 ;; (close-molecule imol-insulin-res-local) ; needed later
			 #t)))))))))


;; The "SHELX Woe" problem
;; 
(greg-testcase "Write a INS from PDB test" #t 
   (lambda ()	      

     (if (not (valid-model-molecule? imol-rnase))
         (begin
            (format #t "   imol-rnase not valid.~%")
	    (throw 'fail))
	 (let* ((rnase-ins "rnase.ins")
		(status (write-shelx-ins-file imol-rnase rnase-ins)))
	   (if (not (= status 1)) 
	       (begin
		 (format #t "   failure to write INS file ~s from PDB~%" rnase-ins)
		 (throw 'fail))
	       #t)))))



(greg-testcase "Add water to SHELX molecule" #t
   (lambda ()

     (set-pointer-atom-molecule imol-insulin-res)
     (set-rotation-centre 3 -1 60)
     (place-typed-atom-at-pointer "Water")
     ;; test is to have to occupancy of the new HOH to be 11.0
     (let ((chain-id (water-chain imol-insulin-res)))
       (let ((n-residues (chain-n-residues chain-id imol-insulin-res)))
	 (let ((serial-number (- n-residues 1)))
	   (let ((res-name (resname-from-serial-number imol-insulin-res chain-id serial-number))
		 (res-no   (seqnum-from-serial-number  imol-insulin-res chain-id serial-number))
		 (ins-code (insertion-code-from-serial-number imol-insulin-res chain-id serial-number)))
		  
	     (if (string=? res-name "HOH")
		 (let ((atom-list (residue-info imol-insulin-res chain-id res-no ins-code)))
		   (for-each 
		    (lambda (atom)
		      (let ((occ (car (car (cdr atom)))))
			(if (not (close-float? occ 11.0))
			    (begin
			      (format #t "  bad occupancy in SHELXL molecule ~s~%" atom)
			      (throw 'fail)))))
		    atom-list)))))
	 #t))))


;(greg-testcase "Aniso Shelx Atoms - Mitch Miller" #t
;   (lambda ()

;     (let* ((bad-anis-res (append-dir-file greg-data-dir "mmiller.res"))
;	    (aniso-bad-mol (read-pdb bad-anis-res)))

;       (if (not (valid-model-molecule? aniso-bad-mol))
;	   (begin
;	     (format #t "   No Mitch Miller molecule yet~%")
;	     (throw 'untested))
;	   (begin
;	     (set-show-aniso 1) ; should (used to) crash
;	     #t)))))

