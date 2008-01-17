
(define insulin-res (append-dir-file greg-data-dir "insulin.res"))
(define imol-insulin (read-pdb insulin-res))

(greg-testcase "NCS chains info" #t
   (lambda ()

     ;; should return #f
     (let ((ncs-chain-info (ncs-chains -1)))
       (if ncs-chain-info
	   (begin
	     (format #t "   Fail: ncs-chains returns ~s, should be #f~%" ncs-chain-info)
	     (throw 'fail))))

     ;; a normal case
     (make-ncs-ghosts-maybe imol-rnase)
     (let ((ncs-chain-info (ncs-chains imol-rnase)))
       (if (not ncs-chain-info)
	   (begin
	     (format #t "   Fail: ncs-chains returns #f~%")
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

