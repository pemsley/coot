
;; convert a view matrix to a view quaternion to set Coot view internals.
;; 
(define matrix->quaternion (lambda (m00 m10 m20 m01 m11 m21 m02 m12 m22)

    ;; From an idea by "Christian" at euclidianspace.com.  The
    ;; rotation matrix is special orthogonal, so (1 + trace) > 0. So
    ;; we can simply do a sqrt on the 4 trace variants.  Most of the
    ;; code here is simply recovering the sign.
    
    ;; return x with the sign of y
    (define convert-sign 
      (lambda (x y)
	
	(cond
	 ((and (> x 0)
	       (> y 0)) x)
	 ((and (< x 0)
	       (> y 0)) (- x))
	 ((and (> x 0)
	       (< y 0)) (- x))
	 (else x))))

    (let ((pw (+ 1 m00 m11 m22))
	  (px (+ 1 m00 (- m11) (- m22)))
	  (py (+ 1 (- m00) m11 (- m22)))
	  (pz (+ 1 (- m00) (- m11) m22)))

      (let ((pr (map (lambda (v)
		       (let ((v1 (if (< v 0) 0 v)))
			 (/ (sqrt v1) 2)))
		     (list pw px py pz))))

	(format #t "prs: ~s~%" pr)

	(append (map (lambda (v signed)
		       (convert-sign v signed))
		     (cdr pr) (list (- m21 m12)
				    (- m02 m20)
				    (- m10 m01)))
	      (list (car pr)))))))
  

; e.g (quat-convert 0.0347695872187614 0.773433089256287 0.632923781871796
;                   0.774806916713715 0.379149734973907 -0.505885183811188
;                  -0.631241261959076 0.507983148097992 -0.586078405380249)
; -> 
; (-0.55715757608 -0.694704711 -7.549694273e-4 0.45492890477)

;(matrix->quaternion 0.0347695 0.7734330  0.632923
;		    0.7748069 0.3791497 -0.50588
;		    -0.631241 0.5079831 -0.58607)

; (matrix->quaternion 0.034 0.77 0.63 0.774  0.379 -0.505 -0.631 0.507 -0.586)
; (matrix->quaternion 1 0 0 0 1 0 0 0 1)

;; Set the view matrix
;; 
(define (set-view-matrix m00 m10 m20 m01 m11 m21 m02 m12 m22)

  (apply set-view-quaternion (matrix->quaternion m00 m10 m20
						 m01 m11 m21
						 m02 m12 m22)))



