
(define hole-ify 
  (let ((start-pos #f)
	(end-pos #f))
    (lambda ()

      (define (status-bar-pos position pos-type)
	(let ((s (apply format #f "Hole ~a point set: (~6,2f ~6,2f ~6,2f)" pos-type position)))
	  (add-status-bar-text s)))


      (let* ((window (gtk-window-new 'toplevel))
	     (vbox (gtk-vbox-new #f 0))
	     (hbox (gtk-hbox-new #f 0))
	     (hbox-pos-buttons (gtk-hbox-new #f 0))
	     (hbox-calc-cancel-buttons (gtk-hbox-new #f 0))
	     (start-button (gtk-button-new-with-label "Set Start Point"))
	     (end-button   (gtk-button-new-with-label "Set End Point"))
	     (calculate-button (gtk-button-new-with-label "Calculate"))
	     (cancel-button    (gtk-button-new-with-label "Cancel")))
	(let ((option-menu-and-model-mol-list (generic-molecule-chooser hbox "HOLE-ify molecule: ")))
	  
	  (gtk-container-add window vbox)
	  (gtk-box-pack-start hbox-pos-buttons start-button #t #t 6)
	  (gtk-box-pack-start hbox-pos-buttons   end-button #t #t 6)
	  (gtk-box-pack-start hbox-calc-cancel-buttons calculate-button #t #t 6)
	  (gtk-box-pack-start hbox-calc-cancel-buttons    cancel-button #t #t 6)
	  (gtk-box-pack-start vbox hbox #t #t 6)
	  (gtk-box-pack-start vbox hbox-pos-buttons #t #t 6)
	  (gtk-box-pack-start vbox hbox-calc-cancel-buttons #t #t 6)
	  
	  (gtk-signal-connect start-button "clicked"
			      (lambda ()
				(set! start-pos (rotation-centre))
				(format #t "Start pos set to: ~s~%" start-pos)
				(status-bar-pos start-pos "start")))


	  (gtk-signal-connect end-button "clicked"
			      (lambda ()
				(set! end-pos (rotation-centre))
				(format #t "End pos set to: ~s~%" end-pos)
				(status-bar-pos end-pos "end")))

	  (gtk-signal-connect calculate-button "clicked"
			      (lambda ()
				(let ((imol (apply get-option-menu-active-molecule 
						   option-menu-and-model-mol-list)))
				  (if (number? imol)
				      (begin
					(format #t "~s ~s~%" start-pos end-pos)
					(if (not (list? start-pos))
					    (begin
					      (add-status-bar-text "Start position not set")
					      (format #t "start pos not set"))
					    (if (not (list? end-pos))
						(begin
						  (add-status-bar-text "End position not set")
						  (format #t "end pos not set"))
						(begin
						  (format #t "hole ~s ~s ~s~%" imol start-pos end-pos)
						  (let ((colour-map-multiplier 1)
							(colour-map-offset 0))
						    (apply hole (append (list imol) start-pos end-pos
									(list colour-map-multiplier colour-map-offset 1 1)))
						    (gtk-widget-destroy window))))))))))
							 

	  (gtk-signal-connect cancel-button "clicked"
			      (lambda ()
				(gtk-widget-destroy window)))
	  
	  (gtk-widget-show-all window))))))




(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Test Hole")))
      (add-simple-coot-menu-menuitem 
       menu "Test" (lambda ()
		     (hole-ify)))))

