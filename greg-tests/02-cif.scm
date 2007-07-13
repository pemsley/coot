
;; NOTE: hof.fcf has  _refln_A_calc and _refln_B_calc, not fcalc and phase.
;; This will not work.  Why is it here?
;; 
(define hof-fcf "data/shelx/egs/hof.fcf")
(define hof-res "data/shelx/egs/HOF.RES")

(define imol-hof-res #f)

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


;(greg-testcase "read shelx fcf" #t 
;   (lambda ()

;     ;; close imol-hof-res too if it was set properly.
;     (if (not (string? hof-fcf))
;	 (begin
;	   (format #t "hof-fcf not defined - skipping test~%")
;	   (throw 'untested))
;	 (begin
;	   (if (not (file-exists? hof-fcf))
;	       (begin
;		 (format #t "~s does not exist - skipping test~%" hof-fcf)
;		 (throw 'untested))
	       
;	       (let ((imol (handle-shelx-fcf-file hof-fcf)))
;		 (let ((success (valid-map-molecule? imol)))
;		   (if success
;		       (begin
;			 (rotate-y-scene 1000 0.1)
;			 (if success
;			     (close-molecule imol))))
;		   (if (valid-model-molecule? imol-hof-res)
;		       (close-molecule imol-hof-res))
;		   success)))))))

