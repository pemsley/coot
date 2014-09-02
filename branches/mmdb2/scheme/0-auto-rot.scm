
(define mol 0) ; molecule number of the coordinates
(define fragment-chain-id "B")
(define mol-for-map 1) ; molecule number of the map


;;;
;;; Make list of integers, a to b: eg (2 3 4 5)
(define number-list
  (lambda (a b)

    (cond
     ((= a b) (list a))
     ((> a b) '())
     (else
      (cons a (number-list (+ a 1) b))))))

(set-go-to-atom-molecule mol)

(for-each 
 (lambda (residue-number)

   ; args: resno inscode chain-id imol-coords imol-map clash-flag lowest-prob
   (add-atom-label    mol fragment-chain-id residue-number       " CA ")
   (remove-atom-label mol fragment-chain-id (- residue-number 3) " CA ")
   (set-go-to-atom-chain-residue-atom-name fragment-chain-id residue-number "CA")
   (rotate-y-scene 40 1.0) ; n-frames frame-interval(degrees)
   (auto-fit-best-rotamer residue-number "" fragment-chain-id mol mol-for-map 0 1.0)
   )

 (number-list 2 25))
