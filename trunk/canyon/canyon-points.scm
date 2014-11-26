

(define (get-canyon-points file-name)
  
  (let ((n (new-generic-object-number "canyon-grid"))
	(g (new-generic-object-number "canyon-grid"))
	(m (new-generic-object-number "path-model"))
	(u (new-generic-object-number "unit-vectors")))

    (call-with-input-file grid-data-file-name
      (lambda (port)

	(let loop ((line (read-line port)))
	  (if (not (eof-object? line))
	      (let ((bits (string->list-of-strings line)))

		(cond 
		 ;; bump point test
		 ;; 
		 ((= (length bits) 9)
		  (let ((x (list-ref bits 5))
			(y (list-ref bits 6))
			(z (list-ref bits 7))
			(col (list-ref bits 8)))
		    (to-generic-object-add-point n col 3
						 (string->number x)
						 (string->number y)
						 (string->number z))
		    (loop (read-line port))))
		 (else 
		  (loop (read-line port)))))))))
    (set-display-generic-object n 1)

    (call-with-input-file file-name
      (lambda (port)

		 ;; line test
		 ;; 
		 ((and (= (length bits) 8)
		       (string=? (car bits) "line"))

		  (let ((x1  (string->number (list-ref bits 1)))
			(y1  (string->number (list-ref bits 2)))
			(z1  (string->number (list-ref bits 3)))
			(x2  (string->number (list-ref bits 4)))
			(y2  (string->number (list-ref bits 5)))
			(z2  (string->number (list-ref bits 6)))
			(col (list-ref bits 7)))
		    (to-generic-object-add-line g col 1 x1 y1 z1 x2 y2 z2))
		  (loop (read-line port)))

		 ;; unit vector
		 ((string=? (car bits) "line-uv")
		  (let ((x1  (string->number (list-ref bits 1)))
			(y1  (string->number (list-ref bits 2)))
			(z1  (string->number (list-ref bits 3)))
			(x2  (string->number (list-ref bits 4)))
			(y2  (string->number (list-ref bits 5)))
			(z2  (string->number (list-ref bits 6)))
			(col (list-ref bits 7)))
		    (to-generic-object-add-line u col 2 x1 y1 z1 x2 y2 z2))
		  (loop (read-line port)))

		 (else 
		  (loop (read-line port)))))

    (set-display-generic-object g 1)))

    ))


(define (canyon-path-points file-name col-1 col-2)
  (call-with-input-file file-name
    (lambda (port)

      (let ((pp (new-generic-object-number "actual-path-points"))
	    (mp (new-generic-object-number "model-path-points")))

      (let loop ((line (read-line port)))
	
	(cond

	 ((eof-object? line) 'done)
	 (else 
	  (let ((bits (string->list-of-strings line)))
	    (if (not (= (length bits) 6))
		
		(begin
		  (loop (read-line port)))
		
		(begin
		  (to-generic-object-add-point pp col-1 4 
					       (string->number (list-ref bits 0))
					       (string->number (list-ref bits 1))
					       (string->number (list-ref bits 2)))
		  (to-generic-object-add-point mp col-2 6
					       (string->number (list-ref bits 3))
					       (string->number (list-ref bits 4))
					       (string->number (list-ref bits 5)))
		  (loop (read-line port))))))))))))



(define (ugo)
  (for-each (lambda (i) (set-display-generic-object i 0)) (range (number-of-generic-objects))))


(get-canyon-points "../../canyon/points")
(canyon-path-points  "../../canyon/canyon-path-fit-initial.tab" "#607080" "#4070aa")
(canyon-path-points  "../../canyon/canyon-path-fit-refined.tab" "#CC7080" "#4070cc")


