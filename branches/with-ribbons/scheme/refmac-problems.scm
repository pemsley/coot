

;The following would be easy from the coot side:

;(list 
; (problem-description place-to-go action-description action-set)
;)

;the place to go can either be an atom, between an atom pair or simply 3D
;coordinates

;place description-type is 'atom 'atom-pair or '3d-coords 
;(place-description-type args ...)

;e.g.
;'('atom molecule-number chain-id residue-number ins-code atom-name
;alt-conf)
;or 
;'('atom-pair molecule-number chain-id-1 residue-number-1 ins-code-1
;atom-name-1 alt-conf-1 chain-id-2 residue-number-2 ins-code-2
;atom-name-2 alt-conf-2)
;or 
;'(45.3 12.1 -9.4)

;action-set can be #f (no Coot commands to fix this problem automatically) or something like

;(lambda ()
;   (flip-peptide molecule-number chain-id residue-number ins-code)
;   (refine-residue-range molecule-number chain-id residue-number ins-code))

(define problem-list
  (list 
   (list "Clash: A 45 CG2 -> B 54 CA" 
	 (list 'atom-pair 
	       (list "A" 45 "" " CG " "")	       
	       (list "A" 45 "" " CG " ""))
	 "Re-fit Rotamer A 45"
	 (lambda (imol)
	   (autofit-best-rotamer 45 "" "" "A" imol 1 0.1)))
   (list "Extra Density near A 21 OG"
	 (list 'atom
	       (list "C" 21 "" " OG " ""))
	 "Add Alt Conf for A 21"
	 (lambda () 
	   (add-alt-conf imol "A" 21 "")))
   (list "High Anisotropy D 43"
	 (list 'atom
	       (list "D" 43 "" " CA " ""))
	 "Add Alt Conf for D 43" 
	 (lambda () 
	   (add-alt-conf imol "D" 43 "")))))

;; The gui for the refmac problems
;; 
(define (refmac-problems-gui imol problem-list)

  (define go-to-position-function
    (lambda (position-description)

      (if (list? position-description)
	  (let ((description-type (car position-description)))
	    
	    (cond
	     ((eq? description-type 'atom-pair) 
	      (lambda (imol)
		(set-go-to-atom-molecule imol)
		(format #t "  chain: ~s~%" (list-ref (car (cdr position-description)) 0))
		(format #t "  resno: ~s~%" (list-ref (car (cdr position-description)) 1))
		(format #t "atomname: ~s~%" (list-ref (car (cdr position-description)) 3))
		(set-go-to-atom-chain-residue-atom-name
		 (list-ref (car (cdr position-description)) 0)
		 (list-ref (car (cdr position-description)) 1)
		 (list-ref (car (cdr position-description)) 3))))
	     ((eq? description-type 'atom) 
	      (lambda (imol)
		(set-go-to-atom-molecule imol)
		(set-go-to-atom-chain-residue-atom-name 
		 (list-ref (car (cdr position-description)) 0)
		 (list-ref (car (cdr position-description)) 1)
		 (list-ref (car (cdr position-description)) 3))))
	     ((eq? description-type '3d-coords)
	      (lambda (imol)
		(apply set-rotation-centre (car (cdr problem-description)))))
	     (else 
	      #f))))))
		       


  (let ((button-ls 
	 (map (lambda (problem)
		(if (not (list? problem))
		    #f
		    (if (not (= (length problem) 4))
			#f
			(let ((button-label1 (list-ref problem 0))
			      (button-label2 (list-ref problem 2))
			      (action-button-2 (list-ref problem 3))
			      (action-button-1 (go-to-position-function
						(list-ref problem 1))))

			  (list button-label1
				action-button-1
				button-label2
				action-button-2)))))
	      problem-list)))

    ;; check that the refinement map has been set, if not, give us a
    ;; message and a chooser (if that's possible)
    (let ((imol-map (imol-refinement-map)))
      (if (= imol-map -1)
	  ;; need to set it then
	  (let ((s "It's a good idea to set the map for refinement now,\n"
		   "since the \"Fix\" actions will use it"))
	    (info-dialog s))))


    (format #t "button list: ~s~%" button-ls)
    (dialog-box-of-pairs-of-buttons imol 
				    "Refmac Noted Problems"
				    (cons 300 200)
				    button-ls "  Close  ")))

(refmac-problems-gui 0 problem-list)


