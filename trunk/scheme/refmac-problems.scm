

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
	       (list "A" 45 "" " CG2" "")	       
	       (list "A" 45 "" " CG2" ""))
	 "Correct Rotamer A 45?"
	 (lambda ()
	   (auto-fit-rotamer-by-atom-spec 45 "" "" "A" imol 1 0.1)))
   (list "Extra Density near C 21 OG"
	 (list 'atom
	       (list "C" 21 "" " OG " ""))
	 "Add Alt Conf"
	 (lambda () 
	   (add-alt-conf imol "C" 21 "")))
   (list "High Anisotropy D 43"
	 (list 'atom
	       (list "D" 43 "" " CA " ""))
	 "Add Alt Conf" 
	 (lambda () 
	   (add-alt-conf imol "D" 43 "")))))

;; The gui for the refmac problems
;; 
(define (refmac-problems-gui imol problem-list)

  (let ((buton-ls 
	 (map (lambda (problem)
		
		problem))))
		
    (dialog-box-of-pairs-of-buttons "Refmac Noted Problems"
				    (cons 300 200)
				    button-ls "  Close  ")))



