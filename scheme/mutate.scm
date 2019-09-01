;;;; Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
 
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
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, 
;;;; Boston, MA 02110-1301, USA

;; Mutate chain-id of molecule number imol given sequence.  This
;; presumes a protein sequence.
;; 
;; The number of residues in chain-id must match the length of sequence.
;; 
(define (mutate-chain imol chain-id sequence)

  (let ((sequence-list (string->list sequence)))

    (if (not (= (length sequence-list) (chain-n-residues chain-id imol)))

	(format #t "sequence mismatch: molecule: ~s new seq: ~s~%"
		(chain-n-residues chain-id imol) 
		(length sequence-list))

	;; do it then
	(begin 
	  
	  (make-backup imol) ; do a backup first
	  (let ((backup-mode (backup-state imol)))
					; turn off backups for imol
	    (turn-off-backup imol) 

	    (let f ((ires 0) (seq sequence-list) (baddies 0))

	      (if (null? seq)
		  
		  (format #t "multi-mutate of ~s residues had ~s ~a.\n"
			  ires baddies (if (= 1 baddies) "error" "errors"))
		  
		  (begin
					; tedious in the real world (mutate in
					; molecule-class-info says more or less the same
					; thing 
					; (format #t "mutating ~s to ~s~%" ires (car seq))
		    
		    (let ((res (mutate-single-residue-by-serial-number
				ires chain-id imol (car seq))))
		      
					;add a baddy if result was 0 (fail) or -1 (fail)
		      (let ((bad-count-addition? (if (not (= res 1)) 1 0)))
			
			(f (+ 1 ires) (cdr seq) (+ baddies bad-count-addition?)))))))

	    (set-have-unsaved-changes imol)
	    (if (= backup-mode 1)
		(turn-on-backup imol))
	    (update-go-to-atom-window-on-changed-mol imol)
	    (graphics-draw))))))

;; An internal function of mutate.  This presumes a protein sequence.
(define multi-mutate
  (lambda (mutate-function imol start-res-no chain-id sequence-list)

    (let f ((ires start-res-no) (seq sequence-list) (baddies 0))
	    
      (if (null? seq)
	  
	  (format #t "multi-mutate of ~s residues had ~s ~a.\n"
		  ires baddies (if (= 1 baddies) "error" "errors"))
	  
	  (begin
	    ; (format #t "mutating ~s to ~s~%" ires (car seq)) ; done elsewhere
	    
	    (let ((result (mutate-function ires "" chain-id imol (car seq))))
	      
              ;add a baddy if result was 0 (fail)
	      (let ((bad-count-addition? (if (= result 0) 1 0)))
		
		(f (+ 1 ires) (cdr seq) (+ baddies bad-count-addition?)))))))))



;; The stop-res-no is inclusive, so usage e.g. (mutate-residue-range 0 "A" 1 2 "AC")
;; 
;; This presumes a protein sequence (not nucleic acid).
;; 
(define (mutate-residue-range imol chain-id start-res-no stop-res-no sequence)

  (if (is-nucleotide-chain? imol chain-id)
      (mutate-nucleotide-range imol chain-id start-res-no stop-res-no sequence)

      ;; protein
      (begin
	(make-backup imol)

	(let ((sequence-list (string->list sequence))
	      (n-residues (+ (- stop-res-no start-res-no) 1)))

	  (if (not (= (length sequence-list) n-residues))
	      (format #t "sequence length mismatch: ~s ~s~%"
		      (length sequence-list) n-residues)

	      (begin
		(let ((backup-mode (backup-state imol)))
		  (turn-off-backup imol)
		  (multi-mutate mutate-single-residue-by-seqno
				imol
				start-res-no
				chain-id
				sequence-list)

		  (set-have-unsaved-changes imol)
		  (if (= backup-mode 1)
		      (turn-on-backup imol))
		  (update-go-to-atom-window-on-changed-mol imol)
		  (graphics-draw))))))))


;; mutate and auto fit a residue range.
;; 
;; This presumes a protein sequence (not nucleic acid).
;; 
;; The sequence is a string of one letter codes
;; 
(define (mutate-and-autofit-residue-range imol chain-id start-res-no stop-res-no sequence)

  (mutate-residue-range imol chain-id start-res-no stop-res-no sequence) ;; does nucleic acids now too
  (let ((mol-for-map (imol-refinement-map))) 
    
    (if (number? mol-for-map)
	(if (>= mol-for-map 0)
	    
	    (let ((backup-mode (backup-state imol)))
	      (turn-off-backup imol)
	      (for-each 
	       (lambda (resno)
		 (let ((clash 1)
		       (altloc "")
		       (inscode ""))
		   (format #t "auto-fit-best-rotamer ~s ~s ~s ~s ~s ~s ~s 0.5~%"
			   resno altloc inscode chain-id imol mol-for-map clash)
		   (let ((score (auto-fit-best-rotamer resno altloc inscode 
						       chain-id imol mol-for-map
						       clash 0.5)))
		     (format #t "   Best score: ~s~%" score))))
	       (number-list start-res-no stop-res-no))
	      (if (= backup-mode 1)
		  (turn-on-backup imol)))
	    
	    (format #t "WARNING:: no map set for refinement.  Can't fit~%")))))

	  
;; mutate and autofit whole chain
;; 
;; This presumes a protein sequence (not nucleic acid).
;; 
(define (mutate-and-auto-fit residue-number chain-id mol mol-for-map residue-type)

  (mutate mol chain-id residue-number "" residue-type)
  (auto-fit-best-rotamer residue-number "" "" chain-id mol mol-for-map 0 0.5))

;; a short-hand for mutate-and-auto-fit
;; 
(define maf mutate-and-auto-fit)

;; return a char, return #\A for unknown residue-type
;; 
(define (3-letter-code->single-letter residue-type)
  (cond 
   ((string=? residue-type "ALA" ) #\A)
   ((string=? residue-type "ARG" ) #\R)
   ((string=? residue-type "ASN" ) #\N)
   ((string=? residue-type "ASP" ) #\D)
   ((string=? residue-type "CYS" ) #\C)
   ((string=? residue-type "GLN" ) #\Q)
   ((string=? residue-type "GLU" ) #\E)
   ((string=? residue-type "GLY" ) #\G)
   ((string=? residue-type "HIS" ) #\H)
   ((string=? residue-type "ILE" ) #\I)
   ((string=? residue-type "LEU" ) #\L)
   ((string=? residue-type "LYS" ) #\K)
   ((string=? residue-type "MET" ) #\M)
   ((string=? residue-type "MSE" ) #\M)
   ((string=? residue-type "PHE" ) #\F)
   ((string=? residue-type "PRO" ) #\P)
   ((string=? residue-type "SER" ) #\S)
   ((string=? residue-type "THR" ) #\T)
   ((string=? residue-type "TRP" ) #\W)
   ((string=? residue-type "TYR" ) #\Y)
   ((string=? residue-type "VAL" ) #\V)
   (else 
    #f)))

;; Mutate a residue range of nucleotides
;; where sequence is a string (for example: "atgccgta")
;;
(define (mutate-nucleotide-range imol chain-id resno-start resno-end sequence)

  ;; Return #f on unknown letter.  We don't need to be specific about
  ;; RNA/DNA here because mutate-base uses canonical_base_name() to
  ;; set the base name according to the residue type before it was
  ;; mutated.
  ;; 
  (define (nucleotide-letter->3-letter-code letter)
    (cond 
     ((not (char? letter)) #f)
     ((eq? letter #\a) "A")
     ((eq? letter #\c) "C")
     ((eq? letter #\g) "G")
     ((eq? letter #\t) "T")
     ((eq? letter #\u) "U")
     ((eq? letter #\A) "A")
     ((eq? letter #\C) "C")
     ((eq? letter #\G) "G")
     ((eq? letter #\T) "T")
     ((eq? letter #\U) "U")
     (else 
      #f)))

  ;; main line
  (if (not (string? sequence))
      (let ((s "sequence must be a string"))
	(info-dialog s)
	(format #t "~s~%" s))
      (let ((residue-range-length (+ (- resno-end resno-start) 1))
	    (seq-length (string-length sequence)))

	(if (not (= residue-range-length seq-length))
	    (let ((s "residue range must equal sequence length"))
	      (info-dialog s)
	      (format #t "~s~%" s))
	    (let ((residue-types 
		   (map nucleotide-letter->3-letter-code 
			(string->list sequence))))

	      (map (lambda (res-no res-type)
		     (if (string? res-type) 
			 (begin
			   (mutate-base imol chain-id res-no "" res-type))))
		     (range resno-start (+ resno-end 1))
		     residue-types))))))



;; a wrapper for mutate-single-residue-by-seqno (which uses slightly
;; inconvenient single letter code)
;; 
;; Here residue-type is the 3-letter code
;; 
(define mutate-residue-redundant
  (lambda (residue-number chain-id imol residue-type . insertion-code)

    (let ((ins-code (if (null? insertion-code) "" (car insertion-code)))
	  (res-type-1lc
	   (cond 
	    ((string=? residue-type "ALA" ) #\A)
	    ((string=? residue-type "ARG" ) #\R)
	    ((string=? residue-type "ASN" ) #\N)
	    ((string=? residue-type "ASP" ) #\D)
	    ((string=? residue-type "CYS" ) #\C)
	    ((string=? residue-type "GLN" ) #\Q)
	    ((string=? residue-type "GLU" ) #\E)
	    ((string=? residue-type "GLY" ) #\G)
	    ((string=? residue-type "HIS" ) #\H)
	    ((string=? residue-type "ILE" ) #\I)
	    ((string=? residue-type "LEU" ) #\L)
	    ((string=? residue-type "LYS" ) #\K)
	    ((string=? residue-type "MET" ) #\M)
	    ((string=? residue-type "PHE" ) #\F)
	    ((string=? residue-type "PRO" ) #\P)
	    ((string=? residue-type "SER" ) #\S)
	    ((string=? residue-type "THR" ) #\T)
	    ((string=? residue-type "TRP" ) #\W)
	    ((string=? residue-type "TYR" ) #\Y)
	    ((string=? residue-type "VAL" ) #\V)
	    (else #\A))))
      (mutate imol chain-id residue-number "" res-type-1lc))))

;;; Prompted by Tim Gruene's email to CCP4bb 20060201.  
;;; Turn all residues (including GLY) of imol to ALA.
;;; 
;;; type is an optional argument.  if type is 'SER then build polySer,
;;; if type is 'GLY, build polyGly.
;;; 
(define (poly-ala imol . type) 

  (if (valid-model-molecule? imol)
      (begin
	(make-backup imol)
	(let ((backup-mode (backup-state imol))
	      (imol-map (imol-refinement-map)))

	  (turn-off-backup imol)

	  (let ((single-letter-code 
		 (cond
		  ((null? type) #\A)
		  ((eq? (car type) 'SER) #\S)
		  ((eq? (car type) 'Ser) #\S)
		  ((eq? (car type) 'GLY) #\G)
		  ((eq? (car type) 'Gly) #\G)
		  (else 
		   #\A))))

	    (map (lambda (chain-id)
		   (let ((n-residues (chain-n-residues chain-id imol)))
		     
		     (for-each 
		      (lambda (serial-number)
			
			(mutate-single-residue-by-serial-number 
			 serial-number chain-id imol single-letter-code))
		      
		      (number-list 0 (- n-residues 1)))))
		 (chain-ids imol))
	    
	    (if (= backup-mode 1)
		(turn-on-backup imol)))))))




;; Delete (back to the CB stub) the side change in the range
;; resno-start to resno-end
;; 
(define delete-sidechain-range
   (lambda (imol chain-id resno-start resno-end)
     
  (if (valid-model-molecule? imol)
      (begin
	(make-backup imol)
	(let ((backup-mode (backup-state imol))
	      (imol-map (imol-refinement-map)))

	  (turn-off-backup imol)

	  (map (lambda (resno)
		 (delete-residue-sidechain imol chain-id resno "" 0))
	       (number-list resno-start resno-end))

	  (if (= backup-mode 1)
	      (turn-on-backup imol)))))))

