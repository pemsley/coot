
;; NOTE: hof.fcf has  _refln_A_calc and _refln_B_calc, not fcalc and phase.
;; This will not work.  Why is it here?
;; 
(define hof-fcf (append-dir-file greg-data-dir "hof.fcf"))
(define hof-res (append-dir-file greg-data-dir "HOF.RES"))

(define imol-hof-res #f)

(define insulin-fcf (append-dir-file greg-data-dir "insulin.fcf"))
(define insulin-res (append-dir-file greg-data-dir "insulin.res"))

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


(greg-testcase "read shelx insulin with fcf" #t 
   (lambda ()

     (let ((imol-insulin-res
	    (handle-read-draw-molecule-with-recentre insulin-res 1)))
       (if (not (valid-model-molecule? imol-insulin-res))
	   (begin
	     (format #t "Bad insulin.res~%")
	     (throw 'fail)))

       (let ((imol (handle-shelx-fcf-file insulin-fcf)))
	 (if (not (valid-map-molecule? imol))
	     (begin
	       (format #t "Bad read of insulin.fcf~%")
	       (throw 'fail))
	     (begin
	       (rotate-y-scene (rotate-n-frames 200) 0.1)
	       (close-molecule imol)
	       (close-molecule imol-insulin-res)
	       #t))))))



