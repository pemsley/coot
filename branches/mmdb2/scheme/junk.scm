	    (lambda (port)
	      (call-with-input-string s
	       (lambda (str-port)
		(read-line str-port)
		(read-line str-port)
		(let f ((str (read-line str-port)))
		  (if (not (eof-object? str))
		      (begin
			(display str port)
			(newline port)
			(f (read-line str-port)))))))))


	  (data-lines (append (list "WEIGHT AUTO" 
				    "MAKE HYDROGENS NO" ; Garib's suggestion 8 Sept 2003
				    )
			      (get-refmac-extra-params)
			      (list labin-string)))


(define (atom-specs imol chain-id resno ins-code atom-name alt-conf)

  (eval 
   (call-with-input-string (atom-info-string imol chain-id resno ins-code atom-name alt-conf) 
			   (lambda (port) (read port)))
   (interaction-environment)))

  (lambda (m00 m01 m02
	   m10 m11 m12
	   m20 m21 m22)


  (if (not (char-ready? soc))
      (begin
;   	(format #t "nothing on the line...~%")
	(usleep 10000)
	#t)
      (let f ((c (read-char soc))
	      (read-bits '()))
	(cond
	 ((eof-object? c) (format #t "server gone\n"))
	 (else 
	  (let ((ready-char-flag (char-ready? soc)))
	   ; (format #t "DEBUG: ready-char-flag: ~s~%" ready-char-flag)
	    (if ready-char-flag
		(begin
		  (f (read-char soc) (cons c read-bits)))
		(begin
		  ;; was the end? Let's wait a bit to see if the
		  ;; buffer gets refilled - and if it still is not
		  ;; ready, then let's evaluate what we have...
		  (usleep 500000)
		  (format #t "waiting for new char...\n")
		  (let ((ready-char-2-flag (char-ready? soc)))
		    (if ready-char-2-flag 
			(begin
			  (f (read-char soc) (cons c read-bits)))
			(begin (evaluate-char-list read-bits))))))))))))
			  


	 (generic-chooser-and-entry "Molecule for refinement:"
				    "HKL data filename (leave blank for default)"
				    "" 
				    (lambda (imol text)
				      (if (= (string-length text) 0)
					  (shelxl-refine imol)
					  (shelxl-refine imol text))))))
