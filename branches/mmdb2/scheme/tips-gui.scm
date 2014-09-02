;;;; Copyright 2006 by The University of York
;;;; Author: Paul Emsley
 
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
;;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

(use-modules (gtk gdk)
             (gtk gtk)
             (gui event-loop)
             (gui text-output-port))

; (load "tips.scm")
; make the random seed depend on the time.
(set! *random-state* (seed->random-state (current-time)))
(define *coot-tip-number* 0)

;; given a number and a gtk text widget @var{text}, put tip number
;; @var{n} into the widget.
(define (show-coot-tip-from-list n text)
  (let ((tip (string-append "\n  " (list-ref tip-list n))))
    (gtk-text-insert text #f "black" "#bfe6bf" tip -1)))

;; increment the tip number when the user sees a tip
(define increment-coot-tip-number
  (lambda ()
    (set! *coot-tip-number* (+ *coot-tip-number* 1))
    (if (= *coot-tip-number* (length tip-list))
	(set! *coot-tip-number* 0))))

;; decrement the tip number when the user sees a tip
(define decrease-coot-tip-number
  (lambda ()
    (set! *coot-tip-number* (- *coot-tip-number* 1))
    (if (= *coot-tip-number* -1)
	(set! *coot-tip-number* (- (length tip-list) 1)))))

;; run the tips gui.
(define (tips-gui)
  (if *do-coot-tips-flag*
      (let* ((window (gtk-window-new 'toplevel))
	     (text (gtk-text-new #f #f))
	     (scrolled-win (gtk-scrolled-window-new))
	     (cancel-button (gtk-button-new-with-label "  Close  "))
	     (next-tip-button (gtk-button-new-with-label "  Next Tip  "))
	     (prev-tip-button (gtk-button-new-with-label " Previous Tip "))
	     (hbox (gtk-hbox-new #f 0))
	     (vbox (gtk-vbox-new #f 0)))

	(gtk-window-set-default-size window 470 120)
	(gtk-container-add window vbox)
	(gtk-window-set-title window "Tips")
	(gtk-container-border-width vbox 10)

	(gtk-container-add scrolled-win text)
	(gtk-container-add vbox scrolled-win)
	(gtk-container-add vbox hbox)
	(gtk-container-add hbox prev-tip-button)
	(gtk-container-add hbox next-tip-button)
	(gtk-container-add hbox cancel-button)
	(gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

	(gtk-signal-connect cancel-button "clicked"
			    (lambda args
			      (gtk-widget-destroy window)))

	(gtk-signal-connect next-tip-button "clicked"
			    (lambda args
			      (increment-coot-tip-number)
			      (let ((nchars (gtk-text-get-length text)))
				(gtk-text-backward-delete text nchars)
				(show-coot-tip-from-list *coot-tip-number* text))))

	(gtk-signal-connect prev-tip-button "clicked"
			    (lambda args
			      (decrease-coot-tip-number)
			      (let ((nchars (gtk-text-get-length text)))
				(gtk-text-backward-delete text nchars)
				(show-coot-tip-from-list *coot-tip-number* text))))

	(set! *coot-tip-number* (random (length tip-list)))
	(show-coot-tip-from-list *coot-tip-number* text)

	(gtk-widget-show-all window))))

;(tips-gui)
;(gtk-main)



