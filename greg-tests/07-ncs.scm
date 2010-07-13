
(define insulin-res (append-dir-file greg-data-dir "insulin.res"))
(define imol-insulin (read-pdb insulin-res))

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
	 (make-dynamically-transformed-ncs-maps imol-rnase imol-rnase-map 0)
	 (make-dynamically-transformed-ncs-maps imol-rnase imol-map-2 0)
	 ;; 2*2 + 1 new maps should have been made
	 (let ((n-new (graphics-n-molecules)))
	   (if (not (=  n-new (+ n-mols 5)))
	       (begin
		 (print-molecule-names)
		 (format #t "no match in number of molecules ~s ~s~%"
			 n-mols n-new)
		 (throw 'fail))
	       #t))))))



(greg-testcase "NCS chains info" #t
   (lambda ()

     ;; should return #f
     (let ((ncs-chain-info (ncs-chain-ids -1)))
       (if ncs-chain-info
	   (begin
	     (format #t "   Fail: ncs-chains returns ~s, should be #f~%" ncs-chain-info)
	     (throw 'fail))))

     ;; a normal case
     (make-ncs-ghosts-maybe imol-rnase)
     (let ((ncs-chain-info (ncs-chain-ids imol-rnase)))
       (if (not ncs-chain-info)
	   (begin
	     (format #t "   Fail: ncs-chain-ids returns #f~%")
	     (throw 'fail))
	   (if (not (> (length ncs-chain-info) 0))
	       (begin
		 (format #t "   Fail: ncs-chains returns ~%" ncs-chain-info)
		 (throw 'fail))
	       (let ((first-ghost (car ncs-chain-info)))
		 (if (not (> (length first-ghost) 1))
		     (begin
		       (format #t "   Fail: first-ghost ~%" first-ghost)
		       (throw 'fail))

		     (if (not (and (string? (list-ref first-ghost 0))
				   (string? (list-ref first-ghost 1))))
			 (begin
			   (format #t "   Fail: not strings first-ghost ~%" first-ghost)
			   (throw 'fail))
			 (begin
			   (format #t "   NCS info: ~s~%" ncs-chain-info)
			   #t)))))))))

	       
     
(greg-testcase "NCS deviation info" #t
   (lambda ()

     ;; should return #f
     (let ((ncs-chain-info (ncs-chain-differences -1 "XX")))
       (if ncs-chain-info
	   (begin
	     (format #t "   Fail: ncs-chains returns ~s, should be #f~%" ncs-chain-info)
	     (throw 'fail))))

     ;; should return #f for insulin
     (let ((ncs-chain-info (ncs-chain-differences imol-insulin "A")))
       (if ncs-chain-info
	   (begin
	     (format #t "   Fail: ncs-chains for insulin returns ~s, should be #f~%" 
		     ncs-chain-info)
	     (throw 'fail))))

     ;; a normal case
     (make-ncs-ghosts-maybe imol-rnase)
     (let ((ncs-chain-info (ncs-chain-differences imol-rnase "A")))
       (if (not ncs-chain-info)
	   (begin
	     (format #t "   Fail: ncs-chain-differences returns #f~%")
	     (throw 'fail))
	   (if (not (= (length ncs-chain-info) 3))
	       (begin
		 (format #t "   Fail on length: length ncs-chain-differences should be 3 is ~s~%" 
			 (length ncs-chain-info))
		 (format #t "        ncs-chain-differences returns ~s~%" ncs-chain-info)
		 (throw 'fail))
	       (begin
		 ; (format #t "   ncs-chain-differences returns ~s~%" ncs-chain-info)
		 #t))))))



;; Phil spotted this bug (it wouldn't refine across the residue 3->4
;; peptide bond - because out out of order residues
;;
(greg-testcase "NCS Residue Range copy" #t
    (lambda ()

     ;; ls is a list of numbers.  Is it in ascending order?  return #t or #f.
     (define (ascending-order? ls)
       
       (cond
	((null? ls) #t)
	((= (length ls) 1) #t)
	((> (car ls) (car (cdr ls))) #f)
	(else
	 (ascending-order? (cdr ls)))))

     ;; prepare the input
     (let ((imol (read-pdb rnase-pdb)))
       (map (lambda (r)
	      (delete-residue imol "B" r ""))
	    (range 1 4))
       ;; make the ghosts
       (make-ncs-ghosts-maybe imol)
       ;; evaluate the function to be tested
       (copy-residue-range-from-ncs-master-to-others imol "A" 1 6)
       ;; check the result.
       (let* ((chain-id "B")
	      (n-residues  (chain-n-residues chain-id imol))
	      (seqnum-order (map (lambda (serial-number)
				   (seqnum-from-serial-number imol chain-id serial-number))
				 (range n-residues))))

	 ;; return a boolean value 
	 (if (not (ascending-order? seqnum-order))
	     (begin
	       (display seqnum-order)
	       (newline)
	       (display "fail")
	       (newline)
	       #f)

	     ;; Now we can test a mutation.  
	     ;;
	     ;; check that the mutate function returns success.
	     ;; 
	     ;; check that the new residue type is indeed a TRP
	     ;; 
	     ;; make the NCS copy
	     ;;
	     ;; check that the NCS target residue has had it's type
	     ;; updated to TRP.
	     ;; 
	     (let ((mutate-success (mutate imol "A" 2 "" "TRP")))
	       (if (not (= mutate-success 1))
		   (begin
		     (format #t "Mutate fails.~%")
		     (throw 'fail))
		   (let ((rname (residue-name imol "A" 2 "")))
		     (if (not (string=? rname "TRP"))
			 (begin
			   (format #t "Mutate fails - master not a TRP. ~s~%" rname)
			   (throw 'fail))
			 (begin
			   (copy-residue-range-from-ncs-master-to-others imol "A" 1 3)
			   (let ((rname (residue-name imol "B" 2 "")))
			     (if (not (string=? rname "TRP"))
				 (begin
				   (format #t "Mutate fails - peer not a TRP (~s)~%" rname)
				   (throw 'fail))
				 #t))))))))))))


;; Perhaps combine this with previous test, it tests the same function.
;; 
(greg-testcase "NCS Residue Range edit to all chains" #t 
   (lambda ()


     (let ((imol (greg-pdb "pdb1t6q.ent")))

       (mutate imol "A" 50 "" "ASP")
       (skip-to-next-ncs-chain 'forward) ;; generate the ghosts
       
       (copy-residue-range-from-ncs-master-to-others imol "A" 50 50)
       
       ;; did it apply?
       (let ((result 
	      (map (lambda (chain-id-peer)
		     (let ((resname (residue-name imol chain-id-peer 50 "")))
		       (string=? resname "ASP")))
		   (list "B" "C"))))

	 (format #t "result: ~s~%" result)
	 (all-true? result)))))


;; This excercises a failure reported by Engin Ozkan 20081209.  Oh
;; D'oh, I'd hard-coded "A" as the master chain id into
;; manual-ncs-ghosts, I should have used the beginning of the list of
;; chain-ids instead.
;; 
(greg-testcase "Manual NCS ghosts generates correct NCS chain ids" #t 
   (lambda ()

     (let ((imol (greg-pdb "pdb1hvv.ent")))

       (set-draw-ncs-ghosts imol 1)
       (ncs-control-change-ncs-master-to-chain-id imol "B")
       (make-ncs-ghosts-maybe imol)
       (let ((ncs-ghost-chains-1 (ncs-chain-ids imol)))
	 (manual-ncs-ghosts imol 220 230 (list "B" "A" "C" "D"))
	 (let ((ncs-ghost-chains-2 (ncs-chain-ids imol)))

	   (format #t "   NCS ghost chain IDs pre:  ~S~%" ncs-ghost-chains-1)
	   (format #t "   NCS ghost chain IDs post: ~S~%" ncs-ghost-chains-2)
	   
	   (if (not (equal? ncs-ghost-chains-1 (list (list "B" "A" "C" "D"))))
	       (throw 'fail))

	   (if (not (equal? ncs-ghost-chains-2 (list (list "B" "A" "C" "D"))))
	       (throw 'fail))

	   #t)))))



(greg-testcase "NCS maps overwrite existing maps" #t 
   (lambda ()

     ;; first close all maps that have "NCS found" in the name:
     (for-each 
      (lambda (imol)
	(if (string-match "NCS found" (molecule-name imol))
	    (close-molecule imol)))
      (molecule-number-list))

     (let ((imol (greg-pdb "pdb1hvv.ent"))
	   (imol-map (make-and-draw-map 
		      (append-dir-file greg-data-dir "1hvv_sigmaa.mtz")
		      "2FOFCWT" "PH2FOFCWT" "" 
		      0 0)))
       
       (make-dynamically-transformed-ncs-maps imol imol-map 0)
       (make-dynamically-transformed-ncs-maps imol imol-map 0)
       (make-dynamically-transformed-ncs-maps imol imol-map 1)

       (all-true?
	(map (lambda (chain-id)
	       (let ((test-name (string-append 
				 "Map "
				 (number->string imol-map) " "
				 "NCS found from matching Chain "
				 chain-id
				 " onto Chain A")))
		 (let loop ((molecule-names (map molecule-name (molecule-number-list)))
			    (n-matchers 0))
		   (cond
		    ((null? molecule-names)
		     (format #t "==== test-name: ~s   n-matchers: ~s~%" test-name n-matchers)
		     (if (= n-matchers 2)
			 #t ;; all good
			 (begin ;; bad
			   (format #t "   failed to match 2 molecule to ~s given names: ~%"
				   test-name)
			   (print-molecule-names)
			   #f)))
		    ((string=? test-name (car molecule-names))
		     (loop (cdr molecule-names) (+ n-matchers 1)))
		    (else
		     (loop (cdr molecule-names) n-matchers))))))
	     (list "B" "C" "D"))))))

