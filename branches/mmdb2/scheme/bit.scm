
(handle-read-draw-molecule "29nov/29nov01.cit.ref.pdb")
(make-and-draw-map "29nov/m72.compl.freer.scaled_refmac1.mtz"
		   "FWT" "PHWT" "no_weights" 0 0 1)
(set-go-to-atom-chain-residue-atom-name "A" 1130 "CA")
(scale-zoom 0.3) ; small numbers zoom in
(toggle-idle-function)

(let f ((res-no 1130))
  (if (< res-no 1150)
      (begin
	(set-go-to-atom-chain-residue-atom-name "A" res-no "CA")
	(f (+ res-no 1)))))
