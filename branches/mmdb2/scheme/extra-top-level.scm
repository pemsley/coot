
(use-modules (gtk gdk)
             (gtk gtk))

;; This is what to do when the button is pressed.
;; It can be any guile or coot function.
;; 
(define my-button-callback
  (lambda ()

    (format #t "button pressed!~%")))

;; define a simple window and put a button in it
;; 
(define (my-top-level)

  (let* ((window (gtk-window-new 'toplevel))
	 (button (gtk-button-new-with-label "my button")))

    (gtk-signal-connect button "clicked" my-button-callback)
    (gtk-container-add window button)
    (gtk-container-border-width window 10)
    (gtk-widget-show button)
    (gtk-widget-show window)))

;; make it happen.
(my-top-level)
