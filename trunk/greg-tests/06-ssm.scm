
(greg-testcase "SSM - Frank von Delft's Example" #t 
    (lambda ()

      (let ((imol-a (handle-read-draw-molecule-with-recentre (append-dir-file greg-data-dir "1wly.pdb") 0))
	    (imol-b (handle-read-draw-molecule-with-recentre (append-dir-file greg-data-dir "1yb5.pdb") 1)))
	
	(if (not (and (valid-model-molecule? imol-a)
		      (valid-model-molecule? imol-b)))
	    #f
	    (begin 
	      (graphics-to-ca-plus-ligands-representation imol-a)
	      (graphics-to-ca-plus-ligands-representation imol-b)
	      
	      (superpose-with-atom-selection imol-a imol-b  "A/2-111" "A/6-115" 0)
	      (set-rotation-centre 65.65 -3 -4)
	      (let ((view-number (add-view (list    49.7269 7.69693 3.93221)
					   (list -0.772277 0.277494 0.292497 0.490948)
					   98.9608 "SSM View")))
		(go-to-view-number view-number 1)
		(rotate-y-scene (rotate-n-frames 100) 0.1)
		(set-mol-displayed imol-a 0)
		(set-mol-displayed imol-b 0)
		
		#t)))))) ; didn't crash.  That's success!   Thanks, Frank.

(greg-testcase "SSM - Alice Dawson's Example" #t 
    (lambda ()
      (let ((imol-s (handle-read-draw-molecule-with-recentre (append-dir-file greg-data-dir "1pyd.pdb") 0)))

	(graphics-to-ca-plus-ligands-representation imol-s)
	(set-graphics-window-size 687 452)
	
	(superpose-with-atom-selection imol-s imol-s "A/100-400" "B/50-450" 1)
	(let ((imol-copy (- (graphics-n-molecules) 1)))
	  (graphics-to-ca-plus-ligands-representation imol-copy)
	  (rotate-y-scene (rotate-n-frames 100) 0.1)
	  #t)))) ; didn't crash!  Thanks Alice Dawson

