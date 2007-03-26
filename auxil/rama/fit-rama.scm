
(use-modules (ice-9 format))

(define dimensions 72)
(define base-array (make-array 0.0 dimensions dimensions))

(define read-rama-data 
  (lambda (filename)

    (let ((arr (make-array 0.0 dimensions dimensions))
	  (fac (/ 360.0 dimensions))) ; degrees/step

      (call-with-input-file filename
	(lambda (port)
	  
	  (let f ((s (read-line port)))
	    (if (not (eof-object? s))
		
		; so a line has 3 numbers
		(let ((parts (string-split s #\space)))
		  
		  (if (= 3 (length parts))
		      
		      (let ((phi (string->number (list-ref parts 0)))
			    (psi (string->number (list-ref parts 1)))
			    (val (string->number (list-ref parts 2))))

			(let ((phi-index (inexact->exact (/ phi fac)))
			      (psi-index (inexact->exact (/ psi fac))))
			  
			  (array-set! arr val phi-index psi-index))))

		  (f (read-line port)))))))

      arr)))

;; return the indices of the peak in arr, look for negative peaks too.
;
(define indices-for-peak
  (lambda (arr)

    (let f ((peak-pair (cons -1 -1))
	    (max-value    -1)
	    (phi-index 0)
	    (psi-index 0))
      (cond 
       ((= phi-index dimensions) peak-pair) ; return value
       ((= psi-index dimensions) (f peak-pair max-value (+ phi-index 1) 0))
       ((> (abs (array-ref arr phi-index psi-index)) max-value)
	(f (cons phi-index psi-index)
	   (abs (array-ref arr phi-index psi-index))
	   phi-index (+ psi-index 1)))
       (else (f peak-pair max-value phi-index (+ psi-index 1)))))))


;; What is the value of the "background" of the gaussian? Estimate it 
;; by averaging points about the 3 sigma mark.
;; 
;; return a number
(define gaussian-base-value 
  (lambda (arr peak-pair variance)

    (let ((three-sigma (* (sqrt variance) 3.0)))

    4.0))
    
    
;; 
(define gaussian-value
  (lambda (centre height variance point-x point-y) 

    ; exp[-2(x-xo)^2/s^2]
    (let* ((degrees/step (/ 360 dimensions))
	   (diff-y (* (- (cdr centre) point-y) degrees/step))
	   (diff-x (* (- (car centre) point-x) degrees/step))
	   (diff-y2 (* diff-y diff-y))
	   (diff-x2 (* diff-x diff-x))
	   (z2 (/ (+ diff-y2 diff-x2) variance)))

      (* height (exp (* -0.5 z2))))))


(define remove-gaussian
  (lambda (arr peak-pair height variance)

    ; given a sigma of 10, typically, we need 6 grid points (for
    ; 360/dimension (=5 degree spacing)) eitheir side, so a
    ; 12x12 grid around peak.
    ; 3 is for 3 sigma 
    (let* ((sigma-guess 10)
	   (grid-size (inexact->exact (* 3 (/ sigma-guess (/ 360 dimensions))))))
      
      (format #t "choosing grid size ~s from sigma ~s ~%" grid-size (sqrt variance))
      
      (let loop ((phi-index (- (car peak-pair) grid-size))
		 (psi-index (- (cdr peak-pair) grid-size)))
	
	(cond
	 ((> phi-index (+ (car peak-pair) grid-size)) 'done)
	 ((> psi-index (+ (cdr peak-pair) grid-size)) (loop (+ phi-index 1) 
							    (- (cdr peak-pair) grid-size)))
	 (else 
	  (let ((g (gaussian-value peak-pair height variance phi-index psi-index)))

; 	    (let ((a (sqrt (+ (* (- phi-index (car peak-pair))
; 				 (- phi-index (car peak-pair)))
; 			      (* (- psi-index (cdr peak-pair))
; 				 (- psi-index (cdr peak-pair)))))))
	      
; 	      (format #t "debug: gaussian-value ~s from height ~s at dist: ~s nsigma: ~s~%"
; 		      g height a (* (/ 360 dimensions) (/ a (sqrt variance)))))
	    
	    (let* ((current-pair (cons phi-index psi-index))
		   (current-value (arr-value arr current-pair)))
	      
; 	      (format #t "removing ~s from ~s at ~s~%" 
; 		      g current-value (cons phi-index psi-index))
	      
	      (arr-remove base-array current-pair (- g)) ; so add it
	      (arr-remove arr current-pair g)
	      (loop phi-index (+ psi-index 1))))))))))


		    
;; the engine of all this:
(define gaussian-fit 
  (lambda (arr peak-pair gaussian-height-scale gaussian-variance-scale)

    (define variance-est
      (lambda (peak-value neighbours x-dist)
	
	(let ((y (average (map (lambda (x)
				 (arr-value arr x))
			       neighbours))))

	  (if (< (* y peak-value) 0.0)
	      (begin
		(format #t "bad neighb avearge: peak: ~s average ~s neigb values: ~s~%"
			peak-value y (map (lambda (x)
					    (arr-value arr x))
					  neighbours))

		;; return a variance equal to the grid spacing
		;; (squared because it is variance)
		(let ((degrees/step (/ 360 dimensions)))
		  (* 0.8 degrees/step degrees/step)))
	  
	      ; normal case, we have neighbours average of the same
	      ; sign as the peak:
	      (- (/ (* x-dist x-dist ) (* 2 (log (/ y peak-value)))))))))

    (define neighbours-avearge
      (lambda (neighbours)

	(average (map (lambda (x)
			(arr-value arr x))
		      neighbours))))

    ;; so we have to estimate the height and the width, 

    ;; Let's start off with 10 degrees width and height is peak height
    ;; take that away from the points...

    ;; aim: minimize Sigma d_i_2, where
    ;; d_i = y_i - g(x_i)
    ;; and g(x) = A\exp((x-x_o)^2)
    ;; 
    ;; we have a set of y_i

    ;; what is the sum of the nearest neighbours?
    ;; what is the sum of the next neest neighbours?
    ;; 
    ;; up to 2 cells away, we have 5 different types of neighbours.
    
    ;; So some strategy then...
    ;; the centre of the gaussian goes at peak-pair.
    ;; we then get the neighbours, the 5 different kinds, and average
    ;; each kind of neighbour.
    ;; 
    ;; each of these can give us an estimate of the standard
    ;; deviation.  Let's average these estimates to assign sigma.
    ;; 
    
    (let* ((phi-index (car peak-pair))
	   (psi-index (cdr peak-pair))
	   (neighbours-1 (list (cons (+ phi-index 1) psi-index)
			       (cons (- phi-index 1) psi-index)
			       (cons phi-index (+ psi-index 1))
			       (cons phi-index (- psi-index 1))))
	   (neighbours-2 (list (cons (+ phi-index 1) (+ psi-index 1))
			       (cons (+ phi-index 1) (- psi-index 1))
			       (cons (- phi-index 1) (+ psi-index 1))
			       (cons (- phi-index 1) (- psi-index 1))))
	   (neighbours-3 (list (cons (+ phi-index 2) psi-index)
			       (cons (- phi-index 2) psi-index)
			       (cons phi-index (+ psi-index 2))
			       (cons phi-index (- psi-index 2))))
	   (neighbours-4 (list (cons (+ phi-index 1) (+ psi-index 2))
			       (cons (+ phi-index 1) (- psi-index 2))
			       (cons (- phi-index 1) (+ psi-index 2))
			       (cons (- phi-index 1) (- psi-index 2))
			       (cons (+ phi-index 2) (+ psi-index 1))
			       (cons (+ phi-index 2) (- psi-index 1))
			       (cons (- phi-index 2) (+ psi-index 1))
			       (cons (- phi-index 2) (- psi-index 1))))
	   (neighbours-5 (list (cons (+ phi-index 2) (+ psi-index 2))
			       (cons (+ phi-index 2) (- psi-index 2))
			       (cons (- phi-index 2) (+ psi-index 2))
			       (cons (- phi-index 2) (- psi-index 2))))
	   (n-1-dist (* 1.0 (/ 360 dimensions)))
	   (n-2-dist (* (sqrt 2.0) (/ 360 dimensions)))
	   (n-3-dist (* 2.0 (/ 360 dimensions)))
	   (n-4-dist (* (sqrt 5.0) (/ 360 dimensions)))
	   (n-5-dist (* (sqrt 8.0) (/ 360 dimensions))))

;       (format #t "peat at ~s has neighbours: ~s~%"
; 	      peak-pair neighbours-1)

       (format #t "peak: ~s at ~s with neibs average: ~s ~%"
 	      (arr-value arr peak-pair)
	      peak-pair
	      (neighbours-avearge neighbours-1))

      (let* ((peak-value (arr-value arr peak-pair))
	     (variance-estimate-1 (variance-est peak-value neighbours-1 n-1-dist))
	     (variance-estimate-2 (variance-est peak-value neighbours-2 n-2-dist))
	     (variance-estimate-3 (variance-est peak-value neighbours-3 n-3-dist))
	     (variance-estimate-4 (variance-est peak-value neighbours-4 n-4-dist))
	     (variance-estimate-5 (variance-est peak-value neighbours-5 n-5-dist)))

;  	(format #t "variance estimates: ~s ~s ~s ~s ~s~%" 
;  		variance-estimate-1 variance-estimate-2 variance-estimate-3
;  		variance-estimate-4 variance-estimate-5)
	
	(let ((average-var (/ (+ (* variance-estimate-1 (length neighbours-1))
				 (* variance-estimate-2 (length neighbours-2))
				 (* variance-estimate-3 (length neighbours-3))
				 (* variance-estimate-4 (length neighbours-4))
				 (* variance-estimate-5 (length neighbours-5)))
			      (+ (length neighbours-1)
				 (length neighbours-2)
				 (length neighbours-3)
				 (length neighbours-4)
				 (length neighbours-5)))))

;  	  (format #t "average variance: ~s from variance estimates: ~%   ~s~%    ~s~%   ~s~%   ~s~%   ~s~%"
; 		  average-var variance-estimate-1 variance-estimate-2 
; 		  variance-estimate-3 variance-estimate-4 variance-estimate-5)

	  ; Now we have a position, a height and a variance, 
	  ; let's remove that gaussian from the plot 
	  (remove-gaussian arr peak-pair 
			   (* gaussian-height-scale (arr-value arr peak-pair))
			   (* gaussian-variance-scale average-var)))))))
	    

;; 
(define wrapped-phi-psi-pair
  (lambda (phi-psi-pair)

    (let ((phi (car phi-psi-pair))
	  (psi (cdr phi-psi-pair)))
      (cond
       ((and (>= phi 0) (< phi dimensions)
	     (>= psi 0) (< psi dimensions))
	phi-psi-pair) ; no change
       ((and (>= phi 0) (< phi dimensions))
					; psi needs changing
	(if (< psi 0)
	    (cons phi (+ psi dimensions))
	    (if (>= psi dimensions)
		(cons phi (- psi dimensions))
		'this-cannot-happen-psi)))
       ((and (>= psi 0) (< psi dimensions))
					; phi needs changing
	(if (< phi 0) 
	    (cons (+ phi dimensions) psi)
	    (if (>= phi dimensions)
		(cons (- phi dimensions) psi)
		'this-cannot-happen-phi)))
       
					; both need changing then:
       ((and (< phi 0) (< psi 0))
	(cons (+ phi dimensions) (+ psi dimensions)))
       ((and (< phi 0) (>= dimensions psi))
	(cons (+ phi dimensions) (- psi dimensions)))
       ((and (>= dimensions phi) (< psi 0))
	(cons (- phi dimensions) (+ psi dimensions)))
       ((and (>= dimensions phi) (>= dimensions psi))
	(cons (- phi dimensions) (- psi dimensions)))
       (else 
	'this-cannont-happen-fall-through)))))
      
;; 
(define average 
  (lambda (lsi) 

    (if (null? lsi)
	0.0
	(let loop ((sum 0.0)
		   (ls lsi))
	  (cond
	   ((null? ls) (/ sum (length lsi)))
	   (else 
	    (loop (+ sum (car ls)) (cdr ls))))))))
      

;; does wrapping, phi and psi indices can be negative or = to
;; dimension (or more).
;; 
(define arr-value
  (lambda (arr phi-psi-pair)

    (let ((w (wrapped-phi-psi-pair phi-psi-pair)))
      (array-ref arr (car w) (cdr w)))))

;; 
(define arr-remove 
  (lambda (arr phi-psi-pair val)

    (let* ((w (wrapped-phi-psi-pair phi-psi-pair))
	   (current-value (arr-value arr w)))

      (array-set! arr (- current-value val) (car w) (cdr w)))))
    
;; 
(define unfitted-sum 
  (lambda (arr)

    (let f ((sum 0.0)
	    (phi-index 0)
	    (psi-index 0))
      (cond 
       ((= phi-index dimensions) (sqrt sum)) ; return value
       ((= psi-index dimensions) (f sum (+ phi-index 1) 0))
       (else
	(let ((v (array-ref arr phi-index psi-index)))
	  (f (+ sum (* v v)) phi-index (+ psi-index 1))))))))

(define print-arr 
  (lambda (arr)

    (let f ((phi-index 0)
	    (psi-index 0))
      (cond 
       ((= phi-index dimensions) 'done) ; return value
       ((= psi-index dimensions) 
	(format #t "~%")
	(f (+ phi-index 1) 0))
       (else 
	(format #t "~5f " (array-ref arr phi-index psi-index))
	(f phi-index (+ psi-index 1)))))))

(define print-arr-for-gnuplot
  (lambda (arr filename) 

    (call-with-output-file filename
      (lambda (port)
	(let f ((phi-index 0)
		(psi-index 0)
		(degrees/step (/ 360 dimensions)))
	  (cond 
	   ((= phi-index dimensions) 'done) ; return value
	   ((= psi-index dimensions) 
	    (format port "~%")
	    (f (+ phi-index 1) 0 degrees/step))
	   (else 
	    (let ((v (array-ref arr phi-index psi-index)))
	      (format port "~s ~s ~s~%" 
		      (* phi-index degrees/step)
		      (* psi-index degrees/step) 
		      (if (> v 0)
			  v v))
	      (f phi-index (+ psi-index 1) degrees/step)))))))))


;; first read in the data into arr

(define fit-rama
  (lambda (rounds height-scale var-scale) 
    (let* ((arr (read-rama-data "rama-data.tab"))
	   (initial-sum (unfitted-sum arr)))
      
      (let loop ((n 0))
	
	(cond
	 ((= n rounds ) (* 100.0 (/ (unfitted-sum arr) initial-sum)))
	 (else 
	  (let ((peak-pair (indices-for-peak arr)))
    	    ; let us know if that chose a negative peak:
	    (if (< (arr-value arr peak-pair) 0.0)
		(format #t "Fitted a negative peak (~s) at: ~s~%"
			(arr-value arr peak-pair) peak-pair))
	    
	    (let ((us (unfitted-sum arr)))
	      (format #t "----- unfitted sum: ~s (~s%)~%" us (* 100.0 (/ us initial-sum)))
	      (gaussian-fit arr peak-pair height-scale var-scale) ; and remove it
	      (loop (+ n 1)))))))
;       arr)))
      (* 100.0 (/ (unfitted-sum arr) initial-sum)))))


; (let* ((rounds 500)
;        (arr (fit-rama rounds 0.7 0.55)))
;   (print-arr-for-gnuplot arr (string-append "post-fit-" (number->string rounds)
; 					    ".tab"))
;   (print-arr-for-gnuplot base-array (string-append "base-fit-" (number->string rounds)
;  						   ".tab")))

			 

;; now that I have corrected the estimate of variance (with a factor
;; of 2), both the optimal height-scale and the optimal variance scale
;; are 1.0

(let loop ((height-scale 0.4))
  
   (cond 
    ((> height-scale 2.0) 'done)
    (else 
     (let ((r (fit-rama 30 height-scale 1.0)))
     
       (format #t "height: ~s final remainder: ~s~%" height-scale r)
       (loop (+ height-scale 0.1))))))
    

  

  