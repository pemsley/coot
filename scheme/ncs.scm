;;;; Copyright 2012 by the University of Oxford
;;;; Copyright 2013 by Medical Research Council

;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


(define find-first-model-molecule
  (lambda ()
    
    (for-each 
     (lambda (molecule)
       (if (valid-model-molecule? molecule)
	   (break molecule)))
     (molecule-number-list))))

;; Skip the residue in the next chain (typically of a molecule with
;; NCS) with the same residue number.  If on the last chain, then wrap
;; to beginning.  If it can't find anything then don't move (and put a
;; message in the status bar)
;; 
(define (skip-to-next-ncs-chain direction)

  ;; Given a chain-id and a list of chain-ids, return the chain-id of
  ;; the next chain to be jumped to (use wrapping).  If the list of
  ;; chain-ids is less then length 2, return #f.
  ;; 
  (define (skip-to-next-chain-id-internal this-chain-id chain-id-list)
    ;; (format #t "DEBUG next this-chain-id: ~s chain-id-list: ~s~%" this-chain-id chain-id-list)
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

  ;; and the previous chain id:
  (define (skip-to-previous-chain-id-internal this-chain-id chain-id-list)
    ;; (format #t "DEBUG prev this-chain-id: ~s chain-id-list: ~s~%" this-chain-id chain-id-list)
    (if (< (length chain-id-list) 2)
	#f 
	(let ((l (length chain-id-list)))
	  (call-with-current-continuation 
	   (lambda (break)
	     (for-each (lambda (ch-id-index)
			 (let ((ch-id (list-ref chain-id-list ch-id-index)))
			   (if (string=? ch-id this-chain-id)
			       (break 
				(if (> ch-id-index 0)
				    (list-ref chain-id-list (- ch-id-index 1))
				    (list-ref chain-id-list (- l 1)))))))
		       (range l)))))))
			    

  ;; 
  (define (skip-to-chain imol this-chain-id chain-id-list)
    
    ;; Given a chain-id and a list of chain-ids, return the chain-id of
    ;; the next chain to be jumped to (use wrapping).  If the list of
    ;; chain-ids is less then length 2, return #f.
    ;; 
    (let ((chain-guess 
	   (if (eq? direction 'forward)
	       (skip-to-next-chain-id-internal     this-chain-id chain-id-list)
	       (skip-to-previous-chain-id-internal this-chain-id chain-id-list))))
      
      ;; (format #t "DEBUG chain-guess: ~s~%" chain-guess)

      (if (not (string? chain-guess))
	  chain-guess
	  (if (not (is-solvent-chain? imol chain-guess))	  
	      chain-guess
	      (skip-to-chain imol chain-guess chain-id-list)))))


  (define (get-chain-id-list imol this-chain-id)
    ;; (format #t "get-chain-id-list given ~s ~s~%" imol this-chain-id)
    (let ((att (ncs-chain-ids imol)))
      (if (not att)
	  (chain-ids imol)
	  (let loop ((attempts att))
	    ;; (format #t "attempts: ~s testing for ~s~%" attempts this-chain-id)
	    (cond
	     ((null? attempts) (chain-ids imol))
	     ((string-member? this-chain-id (car attempts))
	      (car attempts))
	     (else (loop (cdr attempts))))))))
	       

  ;; main line 
  ;; 
  ;; First, what is imol?  imol is the go to atom molecule
  ;; 
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
		      (let ((forward-flag (if (eq? direction 'forward) 1 0)))
			(apply-ncs-to-view-orientation-and-screen-centre
			 imol this-chain-id try-next-chain forward-flag)
			#t)))))))))

   
;; A function inspired by a question from Bill Scott.  He wanted to
;; RNA ghosts.  Because RNA does not work with SSM, we need to define
;; the matrix manually.  Let's make a copy of given imol and get
;; the rtop from that.  Typical usage (manual-ncs-ghosts 0 1 10 "A" "C")
;;
;; We can only see one peer at a time with this (each time we do a
;; clear-ncs-ghost-matrices).
;; 
(define (single-manual-ncs-ghost imol resno-start resno-end ref-chain peer-chain)

  (let ((imol-copy (copy-molecule imol)))
    (clear-lsq-matches)
    (add-lsq-match resno-start resno-end ref-chain resno-start resno-end peer-chain 0) ; ref mov - all atoms
    (let ((rtop (apply-lsq-matches imol-copy imol-copy)))
      (close-molecule imol-copy)
      (if (not rtop)
	  (format #t "Failed to get matching matrix~%")
	  (begin
	    (clear-ncs-ghost-matrices imol)
	    (set-draw-ncs-ghosts imol 1)
	    (apply add-ncs-matrix imol peer-chain ref-chain (apply append rtop)))))))


;; chain-id-list is (list "A" "B" "C" "D")), i.e. the
;; reference/target/master chain-id first and then the peers.  This
;; allows us to add many peers at the same time (unlike above
;; function).
;; 
(define (manual-ncs-ghosts imol resno-start resno-end chain-id-list)

  (if (list? chain-id-list)
      (if (> (length chain-id-list) 1)
	  (begin
	    ;; OK, to the standard SSM-based NCS matrices are bad for this
	    ;; molecule, lets use LSQ.
	    (clear-ncs-ghost-matrices imol)
	    (for-each 
	     (lambda (chain-id)
	       (let ((imol-copy (copy-molecule imol)))
		 (clear-lsq-matches)
		 (add-lsq-match resno-start resno-end (car chain-id-list)
				resno-start resno-end chain-id 1)
		 (let ((rtop (apply-lsq-matches imol imol-copy)))
		   (close-molecule imol-copy)
		   (if (not rtop)
		       (format #t "Failed to get LSQ matching matrix ~s~%" chain-id)
		       (begin
			 (let ((master (car chain-id-list)))
			   (if (string? master)
			       (apply add-ncs-matrix imol chain-id master (apply append rtop)))))))))
	     (cdr chain-id-list))))))




;; Update NCS ghosts based on local environment (residues within 6A of
;; (and including) the active residue).
;; 
;; Typically one would bind this function to a key.
;; 
(define (update-ncs-ghosts-by-local-sphere)
  (using-active-atom 
   ;; (clear-ncs-ghost-matrices aa-imol)
   (let* ((ghost-ncs-chain-ids (ncs-chain-ids aa-imol)))
     ;; (print-var ghost-ncs-chain-ids)
     (if (list? ghost-ncs-chain-ids)
	 (let ((ghost-chain-id-list (car ghost-ncs-chain-ids)))
	   ;; (print-var ghost-chain-id-list)
	   (if (list? ghost-chain-id-list)
	       (begin
		 (clear-ncs-ghost-matrices aa-imol)
		 (for-each 
		  (lambda (chain-id)
		    (let ((imol-copy (copy-molecule aa-imol))) ;; temporary
		      (set-mol-displayed imol-copy 0)
		      (set-mol-active    imol-copy 0)
		      (clear-lsq-matches)
		      (let* ((active-res-spec (list aa-chain-id aa-res-no aa-ins-code))
			     (near-residues (residues-near-residue aa-imol active-res-spec 6))
			     (sphere-residues (cons active-res-spec near-residues)))
			(for-each (lambda(residue-spec)
				    ;; (print-var residue-spec)
				    (let ((res-no  (residue-spec->res-no residue-spec)))
				      (add-lsq-match res-no res-no (car ghost-chain-id-list)
						     res-no res-no chain-id 1)))
				  sphere-residues)

			(let ((rtop (apply-lsq-matches aa-imol imol-copy)))
			  (close-molecule imol-copy)
			  (if (not rtop)
			      (begin
				(format #t "Failed to get LSQ matching matrix ~s~%" chain-id))
			      (begin
				(let ((master (car ghost-chain-id-list)))
				  (format #t "chain-id ~s master ~s  rtop: ~s~%" chain-id master rtop)
				  (if (string? master)
				      (apply add-ncs-matrix aa-imol chain-id master (apply append rtop))))))))))
		  (cdr ghost-chain-id-list)))))))))




;; Return the first master chain id (usually there is only one of course) or #f.
;; 
(define (ncs-master-chain-id imol)

  (let ((cids-ls (ncs-chain-ids imol)))
    (if (not (list? cids-ls))
	#f
	(if (= (length cids-ls) 0)
	    #f
	    (let ((cids (car cids-ls)))
	      (if (= (length cids) 0)
		  #f 
		  (car cids)))))))


;; This was designed to create an NCS copy of a ligand (or range of
;; residues) in the active site of one chain to the as yet unoccupied
;; active site of another, i.e. it makes a NCS ligand "D"1 that is a NCS
;; copy of ligand "C"1 using an NCS operator that maps protein chain "A"
;; onto chain "B".
;; 
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
