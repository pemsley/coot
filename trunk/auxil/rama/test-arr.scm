

(define fiddle-arr 
  (lambda (arr)

    (array-set! arr 1 1 1)))

(let ((arr (make-array 0 5 5)))

  (format #t "arr: ~s~%" arr)
  (fiddle-arr arr)
  (format #t "arr: ~s~%" arr))



	