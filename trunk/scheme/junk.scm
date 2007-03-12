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
