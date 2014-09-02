
(define a-molecule
  (list ; a list of models
   (list ; a list of chains
    (list "A"
	  (list 
	   (list 1081 "" "ASN"
		 (list 
		  (list '(" O  " "")
			'(1.0 22.87 " O")
			'(0.287 20.822 40.742))
		  (list '(" C  " "")
			'(1.0 22.17 " C")
			'(0.303 20.8649 39.5379))
		  (list '(" ND2" "")
			'(1.0 29.17 " N")
			'(-1.363 16.5869 37.908))
		  (list '(" OD1" "")
			'(1.0 32.810001373291 " O")
			'(-2.877 18.049 38.597))
		  (list '(" CG " "")
			'(1.0 27.4 " C")
			'(-1.711 17.676 38.529))
		  (list '(" CB " "")
			'(1.0 24.58 " C")
			'(-0.547 18.528 39.102))
		  (list '(" CA " "")
			'(1.0 23.18 " C")
			'(-0.657 19.99 38.7))
		  (list '(" N  " "")
			'(1.0 22.87 " N")
			'(-0.531 20.0909 37.2929))))

	   (list 1082 "" "GLN"
		 (list
		  (list '(" O  " "")
			'(1.0 19.223 " O")
			'(0.419 24.068 39.868))
		  (list '(" C  " "")
			'(1.0 20.54 " C")
			'(1.431 23.549 40.3401))
		  (list '(" NE2" "")
			'(1.0 14.77 " N")
			'(5.637 22.1359 39.458))
		  (list '(" OE1" "")
			'(1.0 25.2299 " O")
			'(5.2859 23.3729 41.34))
		  (list '(" CD " "")
			'(1.0 25.7 " C")
			'(5.08 23.136 40.01))
		  (list '(" CG " "")
			'(1.0 19.98 " C")
			'(4.176 24.007 39.253))
		  (list '(" CB " "")
			'(1.0 21.63 " C")
			'(3.038 23.178 38.604))
		  (list '(" CA " "")
			'(1.0 20.56 " C")
			'(2.118 22.47 39.593))
		  (list '(" N  " "")
			'(1.0 21.99 " N")
			'(1.156 21.65 38.894)))))))))

;; Return a scm-mol
;;
;; 
;; 
(define (jiggled-mol reference-mol current-mol traj-frac)

  (define (jiggle-random)
    (- (/ (random 10000) 10000) 0.5))

  (define (jiggled-pos ref-pos current-pos)
    (if (< traj-frac 0) ;; magic value
	;; make a starting set of coords
	(map (lambda (x1)
	       (+ x1 (* 1.0 (jiggle-random))))
	     ref-pos)
	(map (lambda (x x-ref)
	       (let ((q (- 1 traj-frac)))
		 (+ (* traj-frac x-ref) (* q x)
		    (* 0.4 q (jiggle-random)))))
	     current-pos ref-pos)))

  (define (jiggled-atom ref-atom current-atom)
    (let ((ref-pos (car (cdr (cdr ref-atom))))
	  (cur-pos (car (cdr (cdr current-atom)))))
    (list (car ref-atom)
	  (car (cdr ref-atom))
	  (jiggled-pos ref-pos cur-pos))))

  (define (jiggled-residue ref-res cur-res)
    (list (car ref-res)
	  (car (cdr ref-res))
	  (car (cdr (cdr ref-res)))
	  (map jiggled-atom
	       (car (cdr (cdr (cdr ref-res))))
	       (car (cdr (cdr (cdr cur-res)))))))

  (define (jiggled-chain ref-chain cur-chain)
    (list (car ref-chain)
	  (map jiggled-residue 
	       (car (cdr ref-chain))
	       (car (cdr cur-chain)))))

  (define (jiggled-model ref-model cur-model)
    (map jiggled-chain ref-model cur-model))

  (map jiggled-model reference-mol current-mol))


;;   
(define (disrupt reference-mol biggness)

  (jiggled-mol reference-mol reference-mol -1))

;;
(define max-count 5000)

(let ((mol-no (add-molecule a-molecule "test molecule")))
  (if (not (= mol-no -1))
      (begin
	(apply set-rotation-centre (centre-of-mass mol-no))
	(let loop ((count 0)
		   (current-mol (disrupt a-molecule 0.8)))
	  (cond 
	   ((= count max-count) 'done)
	   (else
	    (let ((new-mol (jiggled-mol a-molecule current-mol
					(/ count max-count))))
	      ;; (format #t "cycle ~s ~s~%" count (/ count max-count))
	      (clear-and-update-molecule mol-no new-mol)
	      (if (gtk-events-pending)
		  (gtk-main-iteration))
	      (loop (+ count 1) new-mol))))))))

