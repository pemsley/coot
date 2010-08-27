
;;; A silly (but somewhat amusing) demo that changes the background colour
;;; over a range of colours.
;;;
;;; The scene is spun too, for extra eye-candy (so have a molecule
;;; loaded before you run this).
;;;
;;; What's the point?  It's to show that we can display the result of
;;; an arbitarily complex computation, i.e, we have a real extension
;;; language, not just a list of commands (like some other molecular
;;; graphics programs).

; n-steps is divided into 3 sets:
(define n-steps 300)

;; flash the background different colours in some uninteresting way.
(define background-demo
  (lambda ()

    ; debugging
    (define b
      (lambda (arg)
	(format #t "~s~%" arg)))

    ; return a list of 3 colour components
    ;
    (define rgb-colours
      (lambda (n-local)

	(let ((r-step (exact->inexact (/ 1 n-steps))))

	  (let ((red (cond 
		      ((< n-local (/ n-steps 3))        (* 3 r-step n-local))
		      ((< n-local (/ (* 2 n-steps) 3))  (- 1 (* 3 r-step (- n-local (/ n-steps 3)))))
		      (else 0.0)))
		(green (cond
		      ((< n-local (/ n-steps 3))        0.0)
		      ((< n-local (/ (* 2 n-steps) 3))  (* 3 r-step (- n-local (/ n-steps 3))))
		      (else (- 1 (* 3 r-step (- n-local (* 2 (/ n-steps 3))))))))
		(blue 0.0))

	    (list red green blue)))))


    (let f ((n-local 0)) ; get incremented
      (if (not (= n-local n-steps))
	    (let ((r-col (rgb-colours n-local)))
; 	      (format #t "colours ~s~%" r-col)
	      (apply set-background-colour r-col)
	      (rotate-y-scene 10 0.1) ; n-frames frame-interval(degrees)
	      (f (+ n-local 1))
	      )))))

