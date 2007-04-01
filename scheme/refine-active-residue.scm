
;; This totally ignores insertion codes.  A clever algorithm would
;; need a re-write, I think.  Well, we'd have at this end a function
;; that took a chain-id res-no-1 ins-code-1 res-no-2 ins-code-2 
;; 
;; And refine-zone would need to be re-written too, of course.  So
;; let's save that for a rainy day (days... (weeks)).
;; 
(define (refine-active-residue-generic side-residue-offset)

  (let ((active-atom (active-residue)))
    
    (if (not active-atom)
	(format #t "No active atom~%")
	(let ((imol      (list-ref active-atom 0))
	      (chain-id  (list-ref active-atom 1))
	      (res-no    (list-ref active-atom 2))
	      (ins-code  (list-ref active-atom 3))
	      (atom-name (list-ref active-atom 4))
	      (alt-conf  (list-ref active-atom 5)))

	  (format #t "active-atom: ~s~%" active-atom)
	  (let ((backup-mode (backup-state imol))
		(imol-map (imol-refinement-map))
		(replacement-state (refinement-immediate-replacement-state)))
	    
	    (if (= imol-map -1)
		(info-dialog "Oops.  Must Select Map to fit to!")
		
		(begin
		  (turn-off-backup imol)
		  (set-refinement-immediate-replacement 1)
		  (refine-zone imol chain-id 
			       (- res-no side-residue-offset)
			       (+ res-no side-residue-offset)
			       alt-conf)
		  (accept-regularizement)))
		  
	    (if (= replacement-state 0)
		(set-refinement-immediate-replacement 0))
	    (if (= backup-mode 1)
		(turn-on-backup imol)))))))
		  

(define (refine-active-residue)
  (refine-active-residue-generic 0))

(define (refine-active-residue-triple)
  (refine-active-residue-generic 1))


(define (auto-fit-rotamer-active-residue)

  (let ((active-atom (active-residue)))
    (if (not active-atom)
	(format #t "No active atom~%")
	(let ((imol      (list-ref active-atom 0))
	      (chain-id  (list-ref active-atom 1))
	      (res-no    (list-ref active-atom 2))
	      (ins-code  (list-ref active-atom 3))
	      (atom-name (list-ref active-atom 4))
	      (alt-conf  (list-ref active-atom 5)))

	  (format #t "active-atom: ~s~%" active-atom)
	  (let ((backup-mode (backup-state imol))
		(imol-map (imol-refinement-map)))
	    
	    (if (= imol-map -1)
		(info-dialog "Oops.  Must Select Map to fit to!")
		
		(begin
		  (turn-off-backup imol)
		  (auto-fit-best-rotamer res-no alt-conf ins-code chain-id imol imol-map 1 0.1)))
		  
	    (if (= backup-mode 1)
		(turn-on-backup imol)))))))

