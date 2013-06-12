;;;; Copyright 1999 by Paul Emsley

;;; These 3 functions not uninfluenced by SICP, it must be said
;;;
(define accumulate 
  (lambda (op initial sequence)

    (if (null? sequence)
	initial
	(op (car sequence)
	    (accumulate op initial (cdr sequence))))))

;;;
(define accumulate-n
  (lambda (op init seqs)
    
    (if (null? (car seqs))
	'()
	(cons (accumulate   op init (map car seqs))
	      (accumulate-n op init (map cdr seqs))))))

;;; transform
;;; Just a (eg 3x3) (no vector) matrix
;;; The matrix does not have to be symmetric.
;;;
(define transpose-simple-matrix
  (lambda (mat)

    (accumulate-n cons '() mat)))

