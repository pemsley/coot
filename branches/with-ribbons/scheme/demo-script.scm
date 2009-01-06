
; For Glasgow Demo:
;(define pdb-file "tetc/11nov.final.sialic.pdb")
;(define mtz-file "tetc/11novsia.5.mtz")
;(define f-col    "2FOFCWT")

; For York Demo:
;
(define pdb-file "29nov/29nov01.cit.ref.pdb")
(define mtz-file "29nov/m72.compl.freer.scaled_refmac1.mtz")
(define f-col "FWT")

(let ((mol (handle-read-draw-molecule pdb-file)))

  (make-and-draw-map mtz-file f-col "PHWT" "no_weights" 0 0 1)
  (set-go-to-atom-chain-residue-atom-name "A" 1130 "CA")
  (scale-zoom 0.3) ; small numbers zoom in
  (set-smooth-scroll-steps 30) ; depends on the computer/graphics
  (set-go-to-atom-molecule mol); 
  
  (let f ((res-no 1130))
    (if (< res-no 1181)
	(begin

	  ; interestingly different CA specs, see
	  ; graphics_info_t::find_atom_index

	  (add-atom-label    mol "A" res-no       " CA ")
	  (remove-atom-label mol "A" (- res-no 3) " CA ")
	  (set-go-to-atom-chain-residue-atom-name "A" res-no "CA")
	  (if (= (modulo res-no 5) 0)
	      (let ((adj (if (even? res-no) 1 0)))
		(for-each (lambda (x)
			    (change-contour-level adj))
			  (list 1 2 3 4 5))))
	  (rotate-y-scene 200 0.3) ; n-frames frame-interval(degrees)
	  (f (+ res-no 1))))))

(toggle-idle-function)
