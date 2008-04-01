
(define insulin-res (append-dir-file greg-data-dir "insulin.res"))
(define imol-insulin (read-pdb insulin-res))

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
	 (if (ascending-order? seqnum-order)
	     #t
	     (begin
	       (display seqnum-order)
	       (newline)
	       (display "fail")
	       (newline)
	       #f))))))

       
       
