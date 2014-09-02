
(define rnase-pir     (append-dir-file greg-data-dir "rnase.pir"))
(define poly-ala-frag (append-dir-file greg-data-dir "crashes_on_cootaneering-v2.pdb"))


(greg-testcase "Assignment of new PIR sequence overwrites old assignment" #t
   (lambda ()

     (let* ((imol (greg-pdb "tutorial-modern.pdb"))
	    (seq-1 "ACDEFGHIKLMNPQ")
	    (seq-2 "ACDEFGHIKLMNPQRST")
	    (pir-seq-1 (string-append ">test\n\n" seq-1 "*"))
	    (pir-seq-2 (string-append ">test\n\n" seq-2 "*")))
	   
       (assign-pir-sequence imol "A" pir-seq-1)
       (let ((si (sequence-info imol)))
	 (let ((seq (assoc "A" si)))

	   (format #t "debug seq: ~s~%" seq)

	   (if (not seq)
	       (begin
		 (format #t "Bad sequence assignment~%")
		 (throw 'fail)))

	   (if (not (string=? (cdr seq) seq-1))
	       (begin
		 (format #t "bad sequence - not matched ~s vs ~s with assoc:~s~%" 
			 (cdr seq) seq-1 seq)
		 (throw 'fail)))

	   (assign-pir-sequence imol "A" pir-seq-2)
	   (let ((si (sequence-info imol)))
	     (let ((seq (assoc "A" si)))
	       (if (not seq)
		   (begin
		     (format #t "Bad sequence assignment - 2 ")
		     (throw 'fail)))

	       (if (not (string=? (cdr seq) seq-2))
		   (begin
		     (format #t "bad sequence - not matched ~s vs ~s with assoc:~s~%" 
			     (cdr seq) seq-2 seq)
		     (throw 'fail)))

	       #t)))))))


(greg-testcase "Cootaneer Beta Strand" #t 
   (lambda ()
     (let ((imol-model (read-pdb poly-ala-frag))
	   (imol-map (make-and-draw-map rnase-mtz "FWT" "PHWT" "" 0 0)))

       (if (not (valid-model-molecule? imol-model))
	   (begin
	     (format #t "bad imol-model: ~s from file ~s ~%" imol-model poly-ala-frag)
	     (throw 'fail)))

       (if (not (file-exists? rnase-pir))
	   (begin
	     (format #t "missing rnase pir file~%")
	     #f)
	   (let ((seq-text (file->string rnase-pir)))
	     
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


