
(define (refine-active-residue)

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
		(add-status-bar-text "Oops.  Must set a map to fit")
		
		(begin
		  (turn-off-backup imol)
		  (set-refinement-immediate-replacement 1)
		  (refine-zone imol chain-id res-no res-no alt-conf)
		  (accept-regularizement)))
		  
	    (if (= replacement-state 0)
		(set-refinement-immediate-replacement 0))
	    (if (= backup-mode 1)
		(turn-on-backup imol)))))))


		  

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
		(add-status-bar-text "Oops.  Must set a map to fit")
		
		(begin
		  (turn-off-backup imol)
		  (auto-fit-best-rotamer res-no alt-conf ins-code chain-id imol imol-map 1 0.1)))
		  
	    (if (= backup-mode 1)
		(turn-on-backup imol)))))))

