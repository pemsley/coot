
(greg-testcase "RNA NCS Ghosts" #t 
    (lambda ()

     (define (jiggle-random)
        (- (/ (random 10000) 10000) 0.25)) ; between -0.5 to 0.5
     
     ;; a bit of ugliness here, we set! additions to
     ;; atom-attribute-settings, rather than consing them with
     ;; recursion over chain-id, serial-number and atom-list.  
     ;; Hmm...  Come back to improve one day?
     ;; 
     (define (jiggle-atoms-of-mol imol)
       (let ((atom-attribute-settings '()))
	 (map (lambda (chain-id)
		(if (not (is-solvent-chain? imol chain-id))
		    (let ((n-residues (chain-n-residues chain-id imol)))
		      (format #t "   There are ~s residues in chain ~s~%" n-residues chain-id)
		      
		      (for-each 
		       (lambda (serial-number)
			 
			 (let ((res-name (resname-from-serial-number imol chain-id serial-number))
			       (res-no   (seqnum-from-serial-number  imol chain-id serial-number))
			       (ins-code (insertion-code-from-serial-number imol chain-id serial-number)))
			   (let ((atom-ls (residue-info imol chain-id res-no ins-code)))
			     (let f ((atom-ls atom-ls))
			       (cond 
				((null? atom-ls) 'done)
				(else
				 (let* ((atom (car atom-ls))
					(compound-name (car atom))
					(atom-name (car compound-name))
					(alt-conf (car (cdr compound-name)))
					(xyz (car (cdr (cdr atom))))
					(x (car xyz))
					(y (car (cdr xyz)))
					(z (car (cdr (cdr xyz)))))
				   (set! atom-attribute-settings
					 (cons (list imol chain-id res-no ins-code atom-name alt-conf
						     "x" (+ 24 x (* 0.3 (jiggle-random))))
					       atom-attribute-settings))
				   (set! atom-attribute-settings
					 (cons (list imol chain-id res-no ins-code atom-name alt-conf
						     "y" (+ 24 y (* 0.3 (jiggle-random))))
					       atom-attribute-settings))
				   (set! atom-attribute-settings
					 (cons (list imol chain-id res-no ins-code atom-name alt-conf
						     "z" (+ 24 z (* 0.3 (jiggle-random))))
					       atom-attribute-settings))
				   (f (cdr atom-ls)))))))))
		       (number-list 0 (- n-residues 1))))))
	      (chain-ids imol))
	 (set-atom-attributes atom-attribute-settings)))

     ;; main body
     (let* ((rna-mol (ideal-nucleic-acid "RNA" "A" 0 "GACUCUAG"))
	    (copy-rna-mol (copy-molecule rna-mol)))
       
       ;; move the view over a bit so we can see the atoms being jiggled
       (let ((rc (rotation-centre)))
	 (set-rotation-centre (+ 12 (car rc))
			      (+ 3 (car (cdr rc)))
			      (car (cdr (cdr rc)))))
       (let ((view-number (add-view (list    74.7079 10.6267 24.3308)
				    (list -0.713385 -0.0433099 -0.105865 -0.691373)
				    70.3919
				    "RNA-builder-view")))
	 (go-to-view-number view-number 1))

       ;; now jiggle the atoms of copy-rna-mol
       (jiggle-atoms-of-mol copy-rna-mol)
       (merge-molecules (list copy-rna-mol) rna-mol)

       (let ((imol-copy (copy-molecule rna-mol)))
	 (clear-lsq-matches)
	 (add-lsq-match 1 6 "A" 1 6 "C" 0) ; ref mov - all atoms
	 (let ((rtop (apply-lsq-matches imol-copy imol-copy)))
	   (close-molecule imol-copy)
	   (close-molecule copy-rna-mol)
; 	   (format #t "Here is my NCS matrix: ~s~%" rtop) 
	   (if (not rtop)
	       (begin 
		 (format #t "Failed to get matching matrix~%")
		 #f)
	       (begin
		 (set-draw-ncs-ghosts rna-mol 1)
		 (apply add-ncs-matrix rna-mol "C" "A" (apply append rtop))
		 (let ((view-number (add-view (list    72.3306 10.6899 24.073)
					      (list -0.240736 -0.674651 -0.690658 -0.0994136)
					      14.9021
					      "RNA-ghots-view")))
		   (go-to-view-number view-number 1)
		   (rotate-y-scene (rotate-n-frames 200) 1)
		   #t))))))))

