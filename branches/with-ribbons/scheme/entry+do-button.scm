
(use-modules (gtk gdk)
             (gtk gtk)
             (gui event-loop)
             (gui entry-port)
             (gui text-output-port))


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
	 (go-button (gtk-button-new-with-label "  Done   ")))

    (gtk-container-add window vbox)
    (file-selector-entry vbox "Hint Text: ")

    (gtk-box-pack-start vbox go-button #f #f)
    (gtk-widget-show-all window)))


			
(teb2)
(gtk-main)

