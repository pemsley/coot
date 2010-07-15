
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
(define imol-insulin-map -1) ; set later
(define m-miller-res (append-dir-file greg-data-dir "miller/shelx-test4-NPD-mini.res"))


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
			 (set! imol-insulin-map imol)
			 (rotate-y-scene (rotate-n-frames 200) 0.1)
			 ;; (close-molecule imol) ; needed later
			 ;; (close-molecule imol-insulin-res-local) ; needed later
			 #t)))))))))


;; The "SHELX Woe" problem
;; 
(greg-testcase "Write an INS from PDB test" #t 
   (lambda ()	      

     ;; First, check return status on a bogus molecule
     ;; 
     (let ((status (write-shelx-ins-file 204050 "no-molecule.ins")))
       (if (not (= status 0))
	   (begin
	     (format #t "bad exit status from write-shelx-ins-file on bogus molecule ~s~%"
		     status)
	     (throw 'fail))))

     ;; Now check return status on a bogus molecule.  The Happy Path.
     ;; 
     (if (not (valid-model-molecule? imol-rnase))
         (begin
            (format #t "   imol-rnase not valid.~%")
	    (throw 'fail))
	 (let* ((rnase-ins "rnase.ins")
		(status (write-shelx-ins-file imol-rnase rnase-ins)))
	   (if (not (= status 1)) 
	       (begin
		 (format #t "   failure to write INS file ~s from PDB: status ~s~%"
			 rnase-ins status)
		 (throw 'fail))
	       #t)))))



(greg-testcase "new molecule by atom selection inherits shelx molecule flag"  #t 
   (lambda ()

     (let ((insulin-frag (new-molecule-by-atom-selection imol-insulin-res "//B/2010-2020")))
       (if (not (valid-model-molecule? insulin-frag))
	   (begin
	     (format #t " bad fragment of insulin res molecule~%")
	     (throw 'fail)))

       (if (not (shelx-molecule? insulin-frag))
	   (begin 
	     (format #t " bad shelx flag from insulin-frag~%")
	     (throw 'fail)))
       
       #t)))
     


(greg-testcase "Addition of Terminal Residue on SHELX molecule has correct occupancy" #t
   (lambda ()

     (let ((insulin-frag (new-molecule-by-atom-selection imol-insulin-res "//B/2010-2020")))
       (if (not (valid-model-molecule? insulin-frag))
	   (begin
	     (format #t " bad fragment of insulin res molecule~%")
	     (throw 'fail)))

       (set-imol-refinement-map imol-insulin-map)
       (add-terminal-residue insulin-frag "B" 2020 "ALA" 1)

       (let ((res-atoms (residue-info insulin-frag "B" 2021 "")))

	 (if (null? res-atoms)
	   (begin
	     (format #t " bad residue info after add terminal residue~%")
	     (throw 'fail)))

	 (if (not 
	      (all-true? 
	       (map (lambda (atom)
		      (let ((b-factor (car (list-ref atom 1))))
			(close-float? b-factor 11.0)))
		    res-atoms)))
	     (begin
	       (format #t " bad occupancides in new residue ~%")
	       (for-each (lambda (at) (format #t "~s~%" at)) res-atoms)
	       (throw 'fail)))

	 #t ;; correct then
	 
	 ))))


(greg-testcase "Add water to SHELX molecule" #t
   (lambda ()

     (set-pointer-atom-molecule imol-insulin-res)
     (set-rotation-centre 3 -1 60)
     (place-typed-atom-at-pointer "Water")
     ;; test is to have to occupancy of the new HOH to be 11.0
     (shelx-waters-all-good-occ? imol-insulin-res)))


;; Tobias Beck test.
(greg-testcase "Find Waters for a SHELXL molecule" #t
   (lambda () 
     
     (let ((n-chains-pre (n-chains imol-insulin-res)))
       (find-waters imol-insulin-map imol-insulin-res 0 0.6 1)
       (let ((n-chains-post (n-chains imol-insulin-res)))
	 (if (not (= n-chains-pre n-chains-post))
	     (begin
	       (format #t "Find waters on a shelx molecule created a new chain ~s ~s~%" 
		       n-chains-pre n-chains-post)
	       #f)
	     (shelx-waters-all-good-occ? imol-insulin-res))))))



;; non positive definite anistropic atom (reported by Mitch Miller)
;; crash test
(greg-testcase "NPD Anisotropic Atom [Mitch Miller]" #t 
   (lambda ()

     (let ((imol-miller (handle-read-draw-molecule-with-recentre 
                         m-miller-res 1)))
       (if (not (valid-model-molecule? imol-miller))
	   (begin
	     (format #t "Bad read of miller test molecule~%")
	     (throw 'fail)))
       (set-show-aniso 1) ; crash
       (rotate-y-scene (rotate-n-frames 20) 0.1)
       (close-molecule imol-miller)
       (set-show-aniso 0)
       #t)))


;; cheesy test, close the shelx molecules
;; 
(greg-testcase "close shelx molecules" #t
   (lambda () 
	       
     (close-molecule imol-insulin-map)
     (close-molecule imol-insulin-res)

     (if (valid-model-molecule? imol-insulin-res)
	 #f
	 (if (valid-map-molecule? imol-insulin-map)
	     #f
	     #t))))


(greg-testcase "Aniso Bs in P21" #t
   (lambda () 

     (define (aniso-b-from-atom atom)
       (let ((occ-etc (list-ref atom 1)))
	 (list-ref occ-etc 1)))

     (let ((imol (greg-pdb "horma-p21.res")))
       (write-shelx-ins-file imol "new-horma.ins")
       (let ((imol-2 (read-pdb "new-horma.ins")))
	 (let ((at-1 (get-atom imol   "A" 4 "" " N1 "))
	       (at-2 (get-atom imol-2 "A" 4 "" " N1 ")))
	   (let ((b-1 (aniso-b-from-atom at-1))
		 (b-2 (aniso-b-from-atom at-2)))
	     
	     (format #t "b-1: ~s ~%" b-1)
	     (format #t "b-2: ~s ~%" b-2)

	     (if (not (list? b-1))
		 (throw 'fail))

	     (if (not (list? b-2))
		 (throw 'fail))

	     (all-true? (map close-float? b-1 b-2))))))))




(greg-testcase "Don't crash on reading a strange HAT file" #t 
   (lambda ()

     (handle-read-draw-molecule-with-recentre (append-dir-file greg-data-dir "crash.hat") 0)
     #t))

