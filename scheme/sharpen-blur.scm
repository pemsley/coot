
(define (sharpen-blur-map-gui)

  (let* ((window (gtk-window-new 'toplevel))
	 (chooser-label "Map for Sharpen/Blur")
	 (entry-hint-text-1 "Sharpen/Blur:")
	 (entry-hint-text-2 "Factor:")
	 (check-button-label "Resample")
	 (default-entry-1-text "20")
	 (default-entry-2-text "1.3")
	 (label (gtk-label-new chooser-label))
	 (vbox (gtk-vbox-new #f 2))
	 (hbox-for-sharpen  (gtk-hbox-new #f 0))
	 (hbox-for-resample (gtk-hbox-new #f 0))
	 (entry-1 (gtk-entry-new))
	 (entry-2 (gtk-entry-new))
	 (entry-label-1 (gtk-label-new entry-hint-text-1))
	 (entry-label-2 (gtk-label-new entry-hint-text-2))
	 (hbox-buttons (gtk-hbox-new #t 2))
	 (menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new))
	 (ok-button (gtk-button-new-with-label "  Make Map  "))
	 (cancel-button (gtk-button-new-with-label " Close "))
	 (h-sep (gtk-hseparator-new))
	 (map-mol-list (fill-option-menu-with-map-mol-options menu))
	 (check-button (gtk-check-button-new-with-label check-button-label))
	 (check-button-for-refinement-map (gtk-check-button-new-with-label "Make the new map the Refinement Map")))

    (gtk-window-set-title window "Coot: Sharpen/Blur Map")
    (gtk-window-set-default-size window 400 100)
    (gtk-container-add window vbox)
    (gtk-widget-set-usize entry-1 50 -1)
    (gtk-widget-set-usize entry-2 50 -1)

    (gtk-box-pack-start vbox label #f #f 5)
    (gtk-box-pack-start vbox option-menu #t #t 0)
    (gtk-box-pack-start vbox hbox-for-sharpen  #f #f 5)
    (gtk-box-pack-start vbox hbox-for-resample #f #f 5)
    (gtk-box-pack-start vbox check-button-for-refinement-map #f #f 5)
    (gtk-box-pack-start vbox h-sep #t #f 2)
    (gtk-box-pack-start vbox hbox-buttons #f #f 5)

    (gtk-box-pack-start hbox-for-resample check-button #f #f 2)
    (gtk-box-pack-start hbox-buttons ok-button #t #f 5)
    (gtk-box-pack-start hbox-buttons cancel-button #f #f 5)
    (gtk-box-pack-start hbox-for-sharpen entry-label-1 #f #f 4)
    (gtk-box-pack-start hbox-for-sharpen entry-1 #f #f 4)
    (gtk-box-pack-start hbox-for-resample entry-label-2 #f #f 4)
    (gtk-box-pack-start hbox-for-resample entry-2 #f #f 4)
    (gtk-entry-set-text entry-1 default-entry-1-text)
    (gtk-entry-set-text entry-2 default-entry-2-text)

    (gtk-toggle-button-set-active check-button-for-refinement-map #t)
    (gtk-widget-set-sensitive entry-2 #f)
    (gtk-widget-set-sensitive entry-label-2 #f)

    (gtk-option-menu-set-menu option-menu menu)

    ;; check button callback
    (gtk-signal-connect check-button "toggled"
			(lambda ()
			  (format #t "toggled\n")
			  (if (gtk-toggle-button-get-active check-button)
			      (begin
				(gtk-widget-set-sensitive entry-2 #t)
				(gtk-widget-set-sensitive entry-label-2 #t))
			      (begin
				(gtk-widget-set-sensitive entry-2 #f)
				(gtk-widget-set-sensitive entry-label-2 #f)))))

    ;; button callbacks:
    (gtk-signal-connect ok-button "clicked"
			(lambda args
			  ;; what is the molecule number of the option menu?
			  (let ((active-mol-no (get-option-menu-active-molecule 
						option-menu map-mol-list)))

			    (if (number? active-mol-no)
				(begin
				  (let* ((text-1 (gtk-entry-get-text entry-1))
					 (text-2 (gtk-entry-get-text entry-2))
					 (state (gtk-toggle-button-get-active check-button)))

				    (if state
					(let ((imol (sharpen-blur-map-with-resampling active-mol-no
										      (string->number text-1)
										      (string->number text-2))))
					  (if (gtk-toggle-button-get-active check-button-for-refinement-map)
					      (set-imol-refinement-map imol)))
					  
					(let ((imol (sharpen-blur-map active-mol-no (string->number text-1))))
					  (if (gtk-toggle-button-get-active check-button-for-refinement-map)
					      (set-imol-refinement-map imol))))))))))


    (gtk-signal-connect cancel-button "clicked"
			(lambda args
			  (gtk-widget-destroy window)))

    (gtk-widget-show-all window)))


(if #f
    (if (defined? 'coot-main-menubar)
	(let ((menu (coot-menubar-menu "Cryo-EM")))

	  (add-simple-coot-menu-menuitem
	   menu "Sharpen/Blur Map"
	   (lambda ()
	     (sharpen-blur-map-gui))))))

