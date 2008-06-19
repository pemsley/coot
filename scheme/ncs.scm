
(define find-first-model-molecule
  (lambda ()
    
    (for-each 
     (lambda (molecule)
       (if (valid-model-molecule? molecules)
	   (break molecule)))
     (molecule-number-list))))

;; Skip the residue in the next chain (typically of a molecule with
;; NCS) with the same residue number.  If on the last chain, then wrap
;; to beginning.  If it can't find anything then don't move (and put a
;; message in the status bar)
;; 
(define (skip-to-next-ncs-chain)

  ;; Given a chain-id and a list of chain-ids, return the chain-id of
  ;; the next chain to be jumped to (use wrapping).  If the list of
  ;; chain-ids is less then length 2, return #f.
  ;; 
  (define (skip-to-chain-internal this-chain-id chain-id-list)
;    (format #t "DEBUG this-chain-id: ~s~%" this-chain-id)
    (if (< (length chain-id-list) 2)
	#f 
	(let loop ((local-chain-id-list chain-id-list))
;	  (format #t "considering local-chain-id-list: ~s~%"
;		  local-chain-id-list)
	  (cond
	   ((null? local-chain-id-list) (car chain-id-list))
	   ((string=? this-chain-id (car local-chain-id-list))
	    (if (null? (cdr local-chain-id-list))
		(car chain-id-list)
		(car (cdr local-chain-id-list))))
	   (else 
	    (loop (cdr local-chain-id-list)))))))


  (define (skip-to-chain imol this-chain-id chain-id-list)
    
    ;; Given a chain-id and a list of chain-ids, return the chain-id of
    ;; the next chain to be jumped to (use wrapping).  If the list of
    ;; chain-ids is less then length 2, return #f.
    ;; 
    (let ((chain-guess (skip-to-chain-internal this-chain-id chain-id-list)))
      
;      (format #t "DEBUG chain-guess: ~s~%" chain-guess)
      (if (not (string? chain-guess))
	  chain-guess
	  (if (not (is-solvent-chain? imol chain-guess))	  
	      chain-guess
	      (skip-to-chain imol chain-guess chain-id-list)))))

  (define (get-chain-id-list imol this-chain-id)
    (let ((att (ncs-chain-ids imol)))
      (if (not att)
	  (chain-ids imol)
	  (let loop ((attempts att))
	    (cond
	     ((null? attempts) (chain-ids imol))
	     ((string-member? this-chain-id (car attempts))
	      (car attempts))
	     (else loop (cdr attempts)))))))
	       

  ;; First, what is imol?  imol is the go to atom molecule
  (let* ((imol (go-to-atom-molecule-number))
	 (this-chain-id (go-to-atom-chain-id))
	 (chains (get-chain-id-list imol this-chain-id))
	 ; (nov (format #t "Got chains: ~s~%" chains))
	 (next-chain (skip-to-chain imol this-chain-id chains)))

;     (format #t "next-chain: ~s~%" next-chain)
     (make-ncs-ghosts-maybe imol)
    
    (let loop ((try-next-chain next-chain)) 
      
      ;; OK, stop trying for next chain if we have looped round
      ;; (as it were) so that we are back at the starting chain:
      ;; e.g. consider the case: ["A" is protein, "B" is water,
      ;; "C" is ligand]
      ;; 
      (if (not next-chain)
	  (add-status-bar-text "No \"NCS Next Chain\" found")
	  
	  (if (not (string=? try-next-chain this-chain-id))
	      
	      (let ((found-atom-state (set-go-to-atom-chain-residue-atom-name-no-redraw
				       try-next-chain 
				       (go-to-atom-residue-number)
				       (go-to-atom-atom-name)
				       0)))

		;; now, did that set-go-to-atom function work (was there a
		;; real atom)?  If not, then that could have been the ligand
		;; chain or the water chain that we tried to go to.  We want
		;; to try again, and we shbould keep trying again until we get
		;; back to this-chain-id - in which case we have a "No NCS
		;; Next Chain atom" status-bar message.
		
; 		(format #t "DEBUG:: gone to : ~s ~s s~%" 
;			try-next-chain
;			(go-to-atom-residue-number)
;			(go-to-atom-atom-name))

; 		(format #t "DEBUG:: found atom state: ~s~%" found-atom-state)

		(if (= found-atom-state 0)
		    ;; then we did *not* find the atom, e.g. next-chain was
		    ;; the water chain
		    (loop (skip-to-chain imol try-next-chain chains))
		    
		    ;; otherwise all was hunky-dorey
		    (begin
		      ; set the orientation:
;		      (format #t "DEBUG:: apply-ncs-to-view-orientation-and-screen-centre ~s ~s ~s~%"
;			      imol this-chain-id try-next-chain)
		      (apply-ncs-to-view-orientation-and-screen-centre
		       imol this-chain-id try-next-chain)
		      #t))))))))

	 
   
;; A function inspired by a question from Bill Scott.  He wanted to
;; RNA ghosts.  Because RNA does not work with SSM, we need to define
;; the matrix manually.  Let's make a copy of given rna-mol and get
;; the rtop from that.  Typical usage (manual-ncs-ghosts 0 1 10 "A" "C")
;; 
(define (manual-ncs-ghosts rna-mol resno-start resno-end ref-chain peer-chain)

  (let ((imol-copy (copy-molecule rna-mol)))
    (clear-lsq-matches)
    (add-lsq-match resno-start resno-end ref-chain resno-start resno-end peer-chain 0) ; ref mov - all atoms
    (let ((rtop (apply-lsq-matches imol-copy imol-copy)))
      (close-molecule imol-copy)
      (if (not rtop)
	  (format #t "Failed to get matching matrix~%")
	  (begin
	    (clear-ncs-ghost-matrices rna-mol)
	    (set-draw-ncs-ghosts rna-mol 1)
	    (apply add-ncs-matrix rna-mol peer-chain ref-chain (apply append rtop)))))))


(define (ncs-ligand imol-protein ncs-master-chain-id imol-ligand chain-id-ligand resno-ligand-start resno-ligand-stop)

  ;; find ghost in ghosts that has a chain-id matching
  ;; chain-id-protein and get its rtop.  Return #f on not finding the
  ;; ghost
  (define (rtop-from-ghost-with-chain-id ghosts chain-id)
    (let loop ((ghosts ghosts))
      (cond 
       ((null? ghosts) #f)
       ((string=? (list-ref (car ghosts) 1) chain-id)
	(list-ref (car ghosts) 3))
       (else 
	(loop (cdr ghosts))))))

  (let ((chain-ids-from-ncs (ncs-chain-ids imol-protein)))
    (if (list? chain-ids-from-ncs)
	(if (> (length chain-ids-from-ncs) 0)
	    (let* ((ligand-selection-string 
		    (string-append "//" chain-id-ligand "/" (number->string resno-ligand-start)
				   "-" (number->string resno-ligand-stop)))
		   (imol-ligand-fragment (new-molecule-by-atom-selection imol-ligand ligand-selection-string))
		   (ghosts (ncs-ghosts imol-protein)))
	      (let loop ((chain-ids-from-ncs chain-ids-from-ncs))
		(cond
		 ((null? chain-ids-from-ncs) 'failure-to-find-protein-ncs-master-chain)
		 ((string=? (car (car chain-ids-from-ncs)) ncs-master-chain-id)
		  (let ((peer-chains (cdr (car chain-ids-from-ncs)))
			(candidate-name "Candidate NCS-related ligand"))
		    (map (lambda (chain-id-protein count)
			   (let ((rtop (rtop-from-ghost-with-chain-id ghosts chain-id-protein)))
			     (if (not rtop)
				 (begin
				   (format #t "Opps - ncs-ligand: Missing ghost rt-op!~%")
				   (info-dialog "Opps - ncs-ligand: Missing ghost rt-op!~%"))
				 (let ((new-lig-mol (copy-molecule imol-ligand-fragment)))
				   (transform-coords-molecule new-lig-mol (inverse-rtop rtop))
				   (set-molecule-name new-lig-mol (string-append 
								   (number->string new-lig-mol)
								   ": "
								   candidate-name
								   " to protein chain "
								   chain-id-protein))))))
			 peer-chains (number-list 1 (length peer-chains)))
		    (molecules-matching-criteria 
		     (lambda (imol)
		       (if (not (valid-model-molecule? imol))
			   #f
			   (let ((name (molecule-name imol)))
			     (if (string-match candidate-name name)
				 (cons name (molecule-centre imol))
				 #f)))))))
		 (else
		  (loop (cdr chain-ids-from-ncs))))))))))
