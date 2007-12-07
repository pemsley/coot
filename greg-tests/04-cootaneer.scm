
(define rnase-pir     (append-dir-file greg-data-dir "rnase.pir"))
(define poly-ala-frag (append-dir-file greg-data-dir "crashes_on_cootaneering.pdb"))

(greg-testcase "Cootaneer Beta Strand" #t 
   (lambda ()
     (let ((imol-model (read-pdb poly-ala-frag))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))

       (if (not (file-exists? rnase-pir))
	   (begin
	     (format #t "missing rnase pir file~%")
	     #f)
	   (let ((seq-text 
		  (call-with-input-file rnase-pir
		    (lambda (port)
		      (let loop  ((lines '())
				  (line (read-line port)))
			(cond
			 ((eof-object? line) 
			  (string-append-with-string (reverse lines) "\n"))
			 (else
			  (loop (cons line lines) (read-line port)))))))))
	     
	     (assign-pir-sequence imol-model "A" seq-text)
	     
	     (set-rotation-centre 64.271 7.036 14.42)
	     
	     (let ((n-atom (closest-atom imol-model)))
	       (if (not n-atom)
		   (begin 
		     (format #t "missing closest atom~%")
		     #f)
		   (let ((imol     (list-ref n-atom 0))
			 (chain-id (list-ref n-atom 1))
			 (resno    (list-ref n-atom 2))
			 (inscode  (list-ref n-atom 3))
			 (at-name  (list-ref n-atom 4))
			 (alt-conf (list-ref n-atom 5)))
		     (format #t "   Cootaneering: imol ~s chain-id ~s resno ~s inscode ~s at-name ~s alt-conf ~s~%"
			     imol chain-id resno inscode at-name alt-conf)
		     (cootaneer imol-map imol (list chain-id resno inscode 
						    at-name alt-conf))
		     #t))))))))


