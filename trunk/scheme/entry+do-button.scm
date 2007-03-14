
(use-modules (gtk gdk)
             (gtk gtk)
             (gui event-loop)
             (gui entry-port)
             (gui text-output-port))

;; vbox is the vbox to which this compound widget should be added.
;; button-press-func is the lambda function called on pressing return
;; or the button, which takes one argument, the entry.
;;
(define entry+do-button
  (lambda (vbox button-label button-press-func)

    (let ((hbox (gtk-hbox-new #f 0))
	  (entry (gtk-entry-new))
	  (button (gtk-button-new-with-label button-label)))
      
      (gtk-box-pack-start hbox entry #t #f 2)
      (gtk-box-pack-start hbox button #f #f 2)
      (gtk-box-pack-start vbox hbox)
      (gtk-signal-connect button "clicked" 
			  (lambda ()
			    (button-press-func entry))))))


;; test that function:
(define (teb)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (go-button (gtk-button-new-with-label "Done")))

    (gtk-container-add window vbox)
    (entry+do-button vbox "Apply" 
		     (lambda (entry)
		       (let ((entry-text (gtk-entry-get-text entry)))
			 (format #t "got text: ~s~%" entry-text))))
    (gtk-box-pack-start vbox go-button)
    (gtk-widget-show-all window)))


;; test that function:
(define (teb2)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (go-button (gtk-button-new-with-label "Done")))

    (gtk-container-add window vbox)
    (entry+do-button vbox "  File...  "
		     (lambda (entry)
		       (let* ((window (gtk-file-selection-new "file selection")))
			 (gtk-signal-connect (gtk-file-selection-ok-button window)
					     "clicked" 
					     (lambda ()
					       (let ((t (gtk-file-selection-get-filename window)))
						 (format #t "~s~%" t)
						 (gtk-entry-set-text entry t))))
			 (gtk-signal-connect (gtk-file-selection-cancel-button window)
					     "clicked" 
					     (lambda ()
					       (gtk-widget-destroy window)))
			 (gtk-widget-show window))))
			  
    (gtk-box-pack-start vbox go-button)
    (gtk-widget-show-all window)))


			
; (teb2)
; (gtk-main)

