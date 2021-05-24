;;;; Copyright 2007 by Paul Emsley
;;;; 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3 of the License, or (at
;;;; your option) any later version.
;;;; 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;;; 02110-1301, USA

;;;; coot-gui is based on Neil Jerram's guile-gui, I modified the
;;;; widget somewhat to make it more pretty.
;;;; 
;;;; Wouldn't it be nice to have command completion?

(use-modules (gtk gdk)
             (gtk gtk)
             (gui event-loop)
             (gui entry-port)
             (gui text-output-port))
(use-modules (ice-9 threads))

;; 
(define (run-gtk-pending-events)
  (if (gtk-events-pending)
      (gtk-main-iteration)))


;; Fire up the coot scripting gui.  This function is called from the
;; main C++ code of coot.  Not much use if you don't have a gui to
;; type functions in to start with.
;; 
(define (coot-gui)
  (let* ((window (gtk-window-new 'toplevel))
         (text (gtk-text-new #f #f))
         (scrolled-win (gtk-scrolled-window-new))
         (entry (gtk-entry-new))
	 (close-button (gtk-button-new-with-label "  Close  "))
         (vbox (gtk-vbox-new #f 0))
	 (hbox (gtk-hbox-new #f 0))
	 (label (gtk-label-new "Command: "))
         (port (make-entry-port entry #:paren-matching-style '(highlight-region)))
         (oport (make-text-output-port text))
         (ifont (gdk-font-load "fixed"))
         (icolour (gdk-color-parse "darkblue")))
    (gtk-window-set-policy window #t #t #f)
    (gtk-window-set-default-size window 550 250)
    ; (gtk-box-pack-start vbox text #t #t 5)
    (gtk-container-add window vbox)
    (gtk-window-set-title window "Coot Scheme Scripting")
    (gtk-container-border-width vbox 5)

    (gtk-box-pack-start hbox label #f #f 0) ; expand fill padding
    (gtk-box-pack-start hbox entry #t #t 0)
    (gtk-box-pack-start vbox hbox #f #f 5)
    (gtk-container-add scrolled-win text)
    (gtk-container-add vbox scrolled-win)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    (gtk-box-pack-end vbox close-button #f #f 5)

    (gtk-entry-set-text entry "(inactive)")

    (add-hook! (entry-read-complete-hook entry)
               (lambda (str)
                 (with-fluids (; (text-output-font ifont)
                               (text-output-colour icolour))
                              (display str))
                 (newline)))

    (gtk-signal-connect text "key-press-event"
			(lambda args
			  ;; (format #t "text key press ~s~%" args)))
			  (format #t "Don't type here.  Type in the white entry bar.~%")
			  #f))

    (gtk-signal-connect close-button "clicked" 
			(lambda args 
			  (gtk-widget-destroy window)))

    (gtk-widget-show-all window)
    (set-current-input-port port)
    (set-current-output-port oport)
    (set-current-error-port oport)
    (top-repl)))


(set-repl-prompt! "coot> ")

; (guile-gui)

;; The callback from pressing the Go button in the smiles widget, an
;; interface to run libcheck.
;; 
(define (handle-smiles-go tlc-entry smiles-entry)
  
  (let ((tlc-text (gtk-entry-get-text tlc-entry))
	(smiles-text (gtk-entry-get-text smiles-entry)))
    (new-molecule-by-smiles-string tlc-text smiles-text)))


;; smiles GUI
(define (smiles-gui)

  (define (smiles-gui-internal)
    (let* ((window (gtk-window-new 'toplevel))
	   (vbox (gtk-vbox-new #f 0))
	   (hbox1 (gtk-hbox-new #f 0))
	   (hbox2 (gtk-hbox-new #f 0))
	   (tlc-label (gtk-label-new "  3-letter code "))
	   (tlc-entry (gtk-entry-new))
	   (smiles-label (gtk-label-new "SMILES string "))
	   (smiles-entry (gtk-entry-new))
	   (text (gtk-label-new
		  (if (enhanced-ligand-coot?)
		      "  [SMILES interface works by using Pyrogen]  "
		      "  [SMILES interface works by using CCP4's LIBCHECK]  ")))
	   (go-button (gtk-button-new-with-label "  Go  ")))

      (gtk-box-pack-start vbox hbox1 #f #f 0)
      (gtk-box-pack-start vbox hbox2 #f #f 4)
      (gtk-box-pack-start vbox text #f #f 2)
      (gtk-box-pack-start vbox go-button #f #f 6)
      (gtk-box-pack-start hbox1 tlc-label #f 0)
      (gtk-box-pack-start hbox1 tlc-entry #f 0)
      (gtk-box-pack-start hbox2 smiles-label #f 0)
      (gtk-box-pack-start hbox2 smiles-entry #t 4)
      (gtk-container-add window vbox)
      (gtk-container-border-width vbox 6)

      (gtk-signal-connect go-button "clicked"
			  (lambda args
			    (handle-smiles-go tlc-entry smiles-entry)
			    (gtk-widget-destroy window)))
      
      (gtk-signal-connect smiles-entry "key-press-event"
			  (lambda (event)
			    (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
				(begin
				  (handle-smiles-go tlc-entry smiles-entry)
				  (gtk-widget-destroy window)))
			    #f))
      
      (gtk-widget-show-all window)))

  ;; first check that libcheck is available... if not put up and info
  ;; dialog.
  (if (command-in-path? libcheck-exe)
      (smiles-gui-internal)
      (info-dialog "You need to setup CCP4 (specifically LIBCHECK) first.")))



;; Generic single entry widget
;; 
;; Pass the hint labels of the entries and a function that gets called
;; when user hits "Go".  The handle-go-function accepts one argument
;; that is the entry text when the go button is pressed.
;; 
(define (generic-single-entry function-label entry-1-default-text go-button-label handle-go-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox1 (gtk-hbox-new #f 0))
	 (hbox2 (gtk-hbox-new #f 0))
	 (hbox3 (gtk-hbox-new #f 0))
	 (function-label (gtk-label-new function-label))
	 (smiles-entry (gtk-entry-new))
	 (cancel-button (gtk-button-new-with-label "  Cancel  "))
	 (go-button (gtk-button-new-with-label go-button-label)))

    (gtk-box-pack-start vbox hbox1 #f #f 0)
    (gtk-box-pack-start vbox hbox2 #f #f 4)
    (gtk-box-pack-start vbox hbox3 #f #f 4)
    (gtk-box-pack-start hbox3 go-button #t #t 6)
    (gtk-box-pack-start hbox3 cancel-button #t #t 6)
    (gtk-box-pack-start hbox1 function-label #f #f 0)
    (gtk-box-pack-start hbox2 smiles-entry #t #t 0)
    (gtk-container-add window vbox)
    (gtk-container-border-width vbox 6)

    (if (string? entry-1-default-text)
	(gtk-entry-set-text smiles-entry entry-1-default-text))

    (gtk-signal-connect cancel-button "clicked"
			(lambda args
			  (gtk-widget-destroy window)))
    
    (gtk-signal-connect go-button "clicked"
			(lambda args
			  (handle-go-function (gtk-entry-get-text smiles-entry))
			  (gtk-widget-destroy window)))
    
    (gtk-signal-connect smiles-entry "key-press-event"
			(lambda (event)
			  (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
			      (begin
				(handle-go-function (gtk-entry-get-text smiles-entry))
				(gtk-widget-destroy window)))
			  #f))

    (gtk-widget-show-all window)))

;; generic double entry widget, now with a check button
;; 
;; pass a the hint labels of the entries and a function
;; (handle-go-function) that gets called when user hits "Go" (which
;; takes two string aguments and the active-state of the check button
;; (either #t of #f).
;;
;; handle-go-function takes 3 arguments, the third of which is the
;; state of the check button.
;; 
;; if check-button-label not a string, then we don't display (or
;; create, even) the check-button.  If it *is* a string, create a
;; check button and add the callback handle-check-button-function
;; which takes as an argument the active-state of the the checkbutton.
;; 
(define (generic-double-entry label-1 label-2 entry-1-default-text entry-2-default-text check-button-label handle-check-button-function go-button-label handle-go-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox1 (gtk-hbox-new #f 0))
	 (hbox2 (gtk-hbox-new #f 0))
	 (hbox3 (gtk-hbox-new #f 0))
	 (tlc-label (gtk-label-new label-1))
	 (tlc-entry (gtk-entry-new))
	 (smiles-label (gtk-label-new label-2))
	 (smiles-entry (gtk-entry-new))
	 (h-sep (gtk-hseparator-new))
	 (cancel-button (gtk-button-new-with-label "  Cancel  "))
	 (go-button (gtk-button-new-with-label go-button-label)))
    
    (gtk-box-pack-start vbox hbox1 #f #f 0)
    (gtk-box-pack-start vbox hbox2 #f #f 0)
    (gtk-box-pack-start hbox3 go-button #f #f 6)
    (gtk-box-pack-start hbox3 cancel-button #t #f 6)
    (gtk-box-pack-start hbox1 tlc-label #f #f 0)
    (gtk-box-pack-start hbox1 tlc-entry #f #f 0)
    (gtk-box-pack-start hbox2 smiles-label #f #f 0)
    (gtk-box-pack-start hbox2 smiles-entry #t #t 0)

    (let ((check-button 
	   (if (string? check-button-label)
	       (let ((c-button (gtk-check-button-new-with-label check-button-label)))
		 (gtk-box-pack-start vbox c-button #f #f 2)
		 (gtk-signal-connect c-button "toggled"
				     (lambda ()
				       (let ((active-state (gtk-toggle-button-get-active 
							    c-button)))
					 (handle-check-button-function active-state))))
		 c-button)
	       #f))) ; the check-button when we don't want to see it
      
      (gtk-box-pack-start vbox h-sep #t #f 3)
      (gtk-box-pack-start vbox hbox3 #f #f 0)
      (gtk-container-add window vbox)
      (gtk-container-border-width vbox 6)
      
      (if (string? entry-1-default-text) 
	  (gtk-entry-set-text tlc-entry entry-1-default-text))
      
      (if (string? entry-2-default-text )
	  (gtk-entry-set-text smiles-entry entry-2-default-text))

      (gtk-signal-connect cancel-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window)))
      
      (gtk-signal-connect go-button "clicked"
			  (lambda args
			    (handle-go-function (gtk-entry-get-text tlc-entry)
						(gtk-entry-get-text smiles-entry)
						(if check-button
						    (gtk-toggle-button-get-active check-button)
						    'no-check-button))
			    (gtk-widget-destroy window)))
      
      (gtk-signal-connect smiles-entry "key-press-event"
			  (lambda (event)
			    (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
				(begin
				  (handle-go-function (gtk-entry-get-text tlc-entry)
						      (gtk-entry-get-text smiles-entry)
						      (if check-button
							  (gtk-toggle-button-get-active check-button)
							  'no-check-button))
				  (gtk-widget-destroy window)))
			    #f))

      (gtk-window-set-default-size window 400 100)
      (gtk-widget-show-all window))))

;; generic double entry widget, now with a check button
;; 
;; 
;; 
;; OLD:
;; 
;; pass a the hint labels of the entries and a function
;; (handle-go-function) that gets called when user hits "Go" (which
;; takes two string aguments and the active-state of the check button
;; (either #t of #f).
;; 
;; if check-button-label not a string, then we don't display (or
;; create, even) the check-button.  If it *is* a string, create a
;; check button and add the callback handle-check-button-function
;; which takes as an argument the active-state of the the checkbutton.
;; 
(define (generic-multiple-entries-with-check-button entry-info-list check-button-info go-button-label handle-go-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox3 (gtk-hbox-new #f 0))
	 (h-sep (gtk-hseparator-new))
	 (cancel-button (gtk-button-new-with-label "  Cancel  "))
	 (go-button (gtk-button-new-with-label go-button-label)))

      
    ;; all the labelled entries
    ;; 
    (let ((entries (map (lambda (entry-info)
			  
			  (let ((entry-1-hint-text (car entry-info))
				(entry-1-default-text (car (cdr entry-info)))
				(hbox1 (gtk-hbox-new #f 0)))
			    
			    (let ((label (gtk-label-new entry-1-hint-text))
				  (entry (gtk-entry-new)))
			      
			      (if (string? entry-1-default-text) 
				  (gtk-entry-set-text entry entry-1-default-text))
			      
			      (gtk-box-pack-start hbox1 label #f 0)
			      (gtk-box-pack-start hbox1 entry #f 0)
			      (gtk-box-pack-start vbox hbox1 #f #f 0)
			      entry)))
			entry-info-list)))
      
      (let ((check-button 
	     (if (not (and (list? entry-info-list)
			   (= (length check-button-info) 2)))
		 (begin
		   (format #t "check-button-info failed list and length test~%")
		   #f)
		 (if (string? (car check-button-info))
		     (let ((c-button (gtk-check-button-new-with-label (car check-button-info))))
		       (gtk-box-pack-start vbox c-button #f #f 2)
		       (gtk-signal-connect c-button "toggled"
					   (lambda ()
					     (let ((active-state (gtk-toggle-button-get-active 
								  c-button)))
					       ((car (cdr check-button-info)) active-state))))
		       c-button)
		     #f)))) ; the check-button when we don't want to see it
	
	(gtk-box-pack-start vbox h-sep #t #f 3)
	(gtk-box-pack-start vbox hbox3 #f #f 0)
	(gtk-container-add window vbox)
	(gtk-container-border-width vbox 6)
	
	(gtk-box-pack-start hbox3 go-button #t #f 6)
	(gtk-box-pack-start hbox3 cancel-button #t #f 6)
	
	(gtk-signal-connect cancel-button "clicked"
			    (lambda args
			      (gtk-widget-destroy window)))
	
	(gtk-signal-connect go-button "clicked"
			    (lambda ()
			      (handle-go-function (map gtk-entry-get-text entries)
						  (if check-button
						      (gtk-toggle-button-get-active check-button)
						      #f))
			      (gtk-widget-destroy window)))
	
	(gtk-widget-show-all window)))))



;; A demo gui to move about to molecules.
;; 
(define (molecule-centres-gui)

  ;; first, we create a window and a frame to be put into it.
  ;; 
  ;; we create a vbox (a vertical box container) that will contain the
  ;; buttons for each of the coordinate molecules
  ;; 
  (let* ((window (gtk-window-new 'toplevel))
	 (frame (gtk-frame-new "Molecule Centres"))
	 (vbox (gtk-vbox-new #f 3)))

    ;; add the frame to the window and the vbox to the frame
    ;; 
    (gtk-container-add window frame)
    (gtk-container-add frame vbox)
    (gtk-container-border-width vbox 6)

    ;; for each molecule, test if this is a molecule that has a
    ;; molecule (it might have been closed or contain a map).  The we
    ;; construct a button label that is the molecule number and the
    ;; molecule name and create a button with that label.
    ;; 
    ;; then we attach a signal to the "clicked" event on the newly
    ;; created button.  In the function that subsequently happen (on a
    ;; click) we add a text message to the status bar and set the
    ;; graphics rotation centre to the centre of the molecule.  Each
    ;; of these buttons is packed into the vbox (the #f #f means no
    ;; expansion and no filling).  The padding round each button is 3
    ;; pixels.
    ;; 
    (for-each 
     (lambda (molecule-number)

       (if (valid-model-molecule? molecule-number)
	   (let* ((name (molecule-name molecule-number))
		  (label (string-append (number->string molecule-number) " " name))
		  (button (gtk-button-new-with-label label)))
	     (gtk-signal-connect button "clicked"
				 (lambda args
				   (let ((s (format #f "Centred on ~a" label)))
				     (add-status-bar-text s)
				     (apply set-rotation-centre (molecule-centre molecule-number)))))

	     (gtk-box-pack-start vbox button #f #f 3)
	     (gtk-widget-show button))))
       
       (molecule-number-list))

    (gtk-widget-show-all window)))



;; old coot test
(define (old-coot?)

  (if (= (random 10) 0)
      ;; when making a new date, recall that localtime has months that are
      ;; 0-indexed
      ;; 
      ;; use (strftime "%d %m %Y" (localtime xxx))
      (let* (;; (new-release-time 1200000000) ; 10 Jan 2008
	     ;; (new-release-time 1205678900) ; 16 Mar 2008 : 0.3.3
	     ;; (new-release-time 1222222222) ; 24 Jul 2008 : 0.4
	     ;; (new-release-time 1237270000) ; 17 Mar 2009   
	     ;; (new-release-time 1250000000) ; 11 Aug 2009 : 0.5
	     ;; (new-release-time 1280000000) ; 24 Jul 2010 : --
	     ;; (new-release-time 1310000000) ;  7 Jul 2011 : 0.6
	     ;; (new-release-time 1330000000) ; 23 Jan 2012
	     (new-release-time 1330000000) ; 12 Nov 2012
	     (time-diff (- (current-time) new-release-time)))
	(if (> time-diff 0)
	    (let ((s (if (> time-diff 8600000) ;; 100 days
			 "You've got an Old Coot!\n\nIt's time to upgrade."
			 (if (= (random 10) 0)
			     ;; Jorge Garcia:
			     "(Nothing says \"patriotism\" like an Ireland shirt...)\n"
			     "You've got an Old Coot!\n\nIt's time to upgrade."))))
	      (info-dialog s))))))

;; (old-coot?)

;; We can either go to a place (in which case the element is a list of
;; button label (string) and 3 numbers that represent the x y z
;; coordinates) or an atom (in which case the element is a list of a
;; button label (string) followed by the molecule-number chain-id
;; residue-number ins-code atom-name altconf)
;; 
;; e.g. (interesting-things-gui 
;;       "Bad things by Analysis X" 
;;       (list (list "Bad Chiral" 0 "A" 23 "" "CA" "A") 
;;             (list "Bad Density Fit" 0 "B" 65 "" "CA" "") 
;;             (list "Interesting blob" 45.6 46.7 87.5)))
;; 
(define (interesting-things-gui title baddie-list)

  (interesting-things-with-fix-maybe title baddie-list))

;; In this case, each baddie can have a function at the end which is
;; called when the fix button is clicked.
;; 
(define (interesting-things-with-fix-maybe title baddie-list)

  ;; does this baddie have a fix at the end?.  If yes, return the
  ;; func, if not, return #f.
  (define baddie-had-fix?
    (lambda (baddie)

      (let ((l (length baddie)))
	(if (< l 5)
	    #f
	    (let ((func-maybe (list-ref baddie (- l 1))))
	      (if (procedure? func-maybe)
		  func-maybe
		  #f))))))

  ;; main body	
  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-vbox (gtk-vbox-new #f 0)))
    
    (gtk-window-set-default-size window 250 250)
    (gtk-window-set-title window title)
    (gtk-container-border-width inside-vbox 4)
    
    (gtk-container-add window outside-vbox)
    (gtk-container-add outside-vbox scrolled-win)
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    
    (let loop ((baddie-items baddie-list))

      (cond 
       ((null? baddie-items) 'done)
       ((not (list? (car baddie-items))) ; skip this badly formatted baddie
	(loop (cdr baddie-items)))
       (else
	(let* ((hbox (gtk-hbox-new #f 0)) ; hbox to contain Baddie button and 
					; Fix it button
	       (label (car (car baddie-items)))
	       (button (gtk-button-new-with-label label)))

	  ;; (gtk-box-pack-start inside-vbox button #f #f 2) no buttons hbox
	  (gtk-box-pack-start inside-vbox hbox #f #f 2)
	  (gtk-box-pack-start hbox button #f #f 2)
	  
	  ;; add the a button for the fix func if it exists.  Add
	  ;; the callback.
	  (let ((fix-func (baddie-had-fix? (car baddie-items))))
	    (if (procedure? fix-func)
		(let ((fix-button (gtk-button-new-with-label "  Fix  ")))
		  (gtk-box-pack-start hbox fix-button #f #f 2)
		  (gtk-widget-show fix-button)
		  (gtk-signal-connect fix-button "clicked" (lambda () 
							     (fix-func))))))
	  
	  (if (= (length (car baddie-items)) 4) ; "e.g ("blob" 1 2 3)

					; we go to a place
	      (gtk-signal-connect button "clicked" 
				  (lambda args 
				    (apply set-rotation-centre
					   (cdr (car baddie-items)))))

					; else we go to an atom 
	      (let ((mol-no (car (cdr (car baddie-items))))
					; current practice, 
					; (atom-info (cdr (cdr (car baddie-items)))))
					; but we need to kludge out just the chain-id resno 
					; and atom name, because we use 
					; set-go-to-atom-chain-residue-atom-name and it doesn't 
					; take args for altconf and ins-code
		    (atom-info (list (list-ref (car baddie-items) 2)
				     (list-ref (car baddie-items) 3)
				     (list-ref (car baddie-items) 5))))

		(gtk-signal-connect button "clicked"
				    (lambda args 
				      (format #t "Attempt to go to chain: ~s resno: ~s atom-name: ~s~%" 
					      (car atom-info) (cadr atom-info) (caddr atom-info))
				      (set-go-to-atom-molecule mol-no)
				      (let ((success (apply set-go-to-atom-chain-residue-atom-name
							    atom-info)))
					(if (= 0 success) ; failed
					    (let* ((new-name (unmangle-hydrogen-name (caddr atom-info)))
						   (success2 (set-go-to-atom-chain-residue-atom-name
							      (car atom-info) (cadr atom-info) new-name)))
					      (if (= 0 success2) 
						  (begin 
						    (format #t "Failed to centre on demangled name: ~s~%"
							    new-name)
						    (set-go-to-atom-chain-residue-atom-name
						     (car atom-info) (cadr atom-info) " CA "))))))))))

	  (loop (cdr baddie-items))))))

    (gtk-container-border-width outside-vbox 4)
    (let ((ok-button (gtk-button-new-with-label "  OK  ")))
      (gtk-box-pack-start outside-vbox ok-button #f #f 6)
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window))))

    (gtk-widget-show-all window)))

;; Fill an option menu with the "right type" of molecules.  If
;; filter-function returns #t then add it.  Typical value of
;; filter-function is valid-model-molecule?
;; 
(define (fill-option-menu-with-mol-options menu filter-function)

      (let ((n-molecules (graphics-n-molecules)))

	(let loop ((mol-no-ls (molecule-number-list))
		   (rlist '()))
	  
	  (cond 
	   ((null? mol-no-ls) (reverse rlist))
	   ((filter-function (car mol-no-ls))
	    (let ((label-str (molecule-name (car mol-no-ls))))
	      (if (string? label-str)
		  (begin
		    (let* ((mlabel-str (string-append 
					(number->string (car mol-no-ls)) " " label-str))
			   (menuitem (gtk-menu-item-new-with-label mlabel-str)))
		      (gtk-menu-append menu menuitem)
		      (gtk-widget-show menuitem))
		    (loop (cdr mol-no-ls) (cons (car mol-no-ls) rlist)))
		  
		  (begin
		    (format #t "OOps molecule name for molecule ~s is ~s~%" 
			    (car mol-no-ls) label-str)
		    (loop (cdr mol-no-ls) rlist)))))
	   (else 
	    (loop (cdr mol-no-ls) rlist))))))
    
(define (fill-option-menu-with-map-mol-options menu)

  (fill-option-menu-with-mol-options menu valid-map-molecule?))

;; Helper function for molecule chooser.  Not really for users.
;; 
;; Return a list of models, corresponding to the menu items of the
;; option menu.
;; 
;; The returned list will not contain references to map or closed
;; molecules.
;; 
(define (fill-option-menu-with-coordinates-mol-options menu)

  (fill-option-menu-with-mol-options menu valid-model-molecule?))

;; 
(define (fill-option-menu-with-number-options menu number-list default-option-value)

  (let loop ((number-list number-list)
	     (count 0))
    (cond 
     ((null? number-list) '())
     (else 
      (let* ((mlabel-str (number->string (car number-list)))
	     (menuitem (gtk-menu-item-new-with-label mlabel-str)))
	(gtk-menu-append menu menuitem)
	(if (= (car number-list) default-option-value)
	    (begin
	      (gtk-menu-set-active menu count)
	      (format #t "setting menu active ~s ~s~%" 
		      default-option-value count)))
	(loop (cdr number-list)
	      (+ count 1)))))))

    
;; Helper function for molecule chooser.  Not really for users.
;; 
;; return the molecule number of the active item in the option menu,
;; or return #f if there was a problem (e.g. closed molecule)
;; 
(define (get-option-menu-active-molecule option-menu model-mol-list)

      (let* ((menu (gtk-option-menu-get-menu option-menu))
	     (active-item (gtk-menu-get-active menu))
	     (children (gtk-container-children menu)))
	
	(if (= (length children) (length model-mol-list))
	    (begin 
	      (let loop ((children children)
			 (model-mol-list model-mol-list))
		
		(cond
		 ((null? children) #f)
		 ((eqv? active-item (car children))
		  (if (or (valid-model-molecule? (car model-mol-list))
			  (valid-map-molecule?   (car model-mol-list)))
		      (car model-mol-list)
		      #f))
		 (else
		  (loop (cdr children) (cdr model-mol-list))))))
	    
	    (begin
	      (format #t "Failed children length test : ~s ~s~%" children model-mol-list)
	      #f))))

;; Here we return the active item in an option menu of generic items
;; 
(define (get-option-menu-active-item option-menu item-list)

      (let* ((menu (gtk-option-menu-get-menu option-menu))
	     (active-item (gtk-menu-get-active menu))
	     (children (gtk-container-children menu)))

	(if (= (length children) (length item-list))
	    (begin 
	      (let loop ((children children)
			 (item-list item-list))
		
		(cond
		 ((null? children) #f)
		 ((eqv? active-item (car children))
		  (car item-list))
		 (else
		  (loop (cdr children) (cdr item-list))))))
	    
	    (begin
	      (format #t "Failed children length test : ~s ~s~%" 
		      children item-list)
	      #f))))



;; Typically option-menu-fill-function is fill-option-menu-with-coordinates-mol-options
;; 
(define (molecule-chooser-gui-generic chooser-label callback-function option-menu-fill-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (label (gtk-label-new chooser-label))
	 (vbox (gtk-vbox-new #f 6))
	 (hbox-buttons (gtk-hbox-new #f 5))
	 (menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new))
	 (ok-button (gtk-button-new-with-label "  OK  "))
	 (cancel-button (gtk-button-new-with-label " Cancel "))
	 (h-sep (gtk-hseparator-new))
	 (model-mol-list (option-menu-fill-function menu)))
    (gtk-window-set-default-size window 370 100)
    (gtk-container-add window vbox)
    (gtk-box-pack-start vbox label #f #f 5)
    (gtk-box-pack-start vbox option-menu #t #t 0)
    (gtk-box-pack-start vbox h-sep #t #f 2)
    (gtk-box-pack-start vbox hbox-buttons #f #f 5)
    (gtk-box-pack-start hbox-buttons ok-button #t #f 5)
    (gtk-box-pack-start hbox-buttons cancel-button #t #f 5)
    
    (gtk-option-menu-set-menu option-menu menu)
    
    ;; button callbacks:
    ;;
    (gtk-signal-connect ok-button "clicked"
			(lambda args
			  ;; what is the molecule number of the option menu?
			  (let ((active-mol-no (get-option-menu-active-molecule 
						option-menu
						model-mol-list)))

			    (if (not (number? active-mol-no))
				(begin ;; something bad happend
				  (format #t "Failed to get a (molecule) number~%"))

				(begin ;; good...
				  (callback-function active-mol-no))))
			  
			  (gtk-widget-destroy window)))

    (gtk-signal-connect cancel-button "clicked"
			(lambda args
			  (gtk-widget-destroy window)))
    
    (gtk-widget-show-all window)))

;; Fire up a coordinates/model molecule chooser dialog, with a given
;; label and on OK we call the call-back-fuction with an argument of
;; the chosen molecule number.
;; 
;; chooser-label is a directive to the user such as "Choose a Molecule"
;; 
;; callback-function is a function that takes a molecule number as an
;; argument.
;; 
(define molecule-chooser-gui
  (lambda (chooser-label callback-function)
    (molecule-chooser-gui-generic chooser-label callback-function 
				  fill-option-menu-with-coordinates-mol-options)))
;; 
;; Fire up a map molecule chooser dialog, with a given label and on OK we
;; call the call-back-fuction with an argument of the chosen molecule
;; number. 
;; 
;; chooser-label is a directive to the user such as "Choose a Molecule"
;; 
;; callback-function is a function that takes a molecule number as an
;; argument.
;; 
(define map-molecule-chooser-gui
  (lambda (chooser-label callback-function)
    (molecule-chooser-gui-generic chooser-label callback-function 
				  fill-option-menu-with-map-mol-options)))

;; A pair of widgets, a molecule chooser and an entry.  The
;; callback-function is a function that takes a molecule number and a
;; text string.
;;
;; if always-dismiss-on-ok-clicked? is false then the dialog is not dismissed if 
;; callback-function returns #f.
;; 
(define (generic-chooser-and-entry chooser-label entry-hint-text defaut-entry-text callback-function . always-dismiss-on-ok-clicked?)

  (format #t "---- deal with this: always-dismiss-on-ok-clicked?: ~s~%" always-dismiss-on-ok-clicked?)
  (generic-chooser-and-entry-and-checkbutton chooser-label entry-hint-text defaut-entry-text #f callback-function always-dismiss-on-ok-clicked?))

;; as above, plus we also have a check-button and an additional argument in the callback
;; If check-button-label is false, then don't create a check-button.
(define (generic-chooser-and-entry-and-checkbutton chooser-label entry-hint-text defaut-entry-text check-button-label callback-function . always-dismiss-on-ok-clicked?)

  ;; (format #t "always-dismiss-on-ok-clicked? is ~s~% " always-dismiss-on-ok-clicked?)
  ;; if this function is called withouth the always-dismiss-on-ok-clicked? then
  ;; always-dismiss-on-ok-clicked? here is ().
  ;; otherwise it is a list containing the passed always-dismiss-on-ok-clicked? value.

  (let* ((window (gtk-window-new 'toplevel))
	 (label (gtk-label-new chooser-label))
	 (vbox (gtk-vbox-new #f 2))
	 (hbox-for-entry (gtk-hbox-new #f 0))
	 (entry (gtk-entry-new))
	 (entry-label (gtk-label-new entry-hint-text))
	 (hbox-buttons (gtk-hbox-new #t 2))
	 (menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new))
	 (ok-button (gtk-button-new-with-label "  OK  "))
	 (cancel-button (gtk-button-new-with-label " Cancel "))
	 (h-sep (gtk-hseparator-new))
	 (model-mol-list (fill-option-menu-with-coordinates-mol-options menu)))
    
    (gtk-window-set-default-size window 400 100)
    (gtk-container-add window vbox)
    (gtk-box-pack-start vbox label #f #f 5)
    (gtk-box-pack-start vbox option-menu #t #t 0)
    (gtk-box-pack-start vbox hbox-for-entry #f #f 5)
    (let ((check-button 
	   (if (string? check-button-label)
	       (let ((c-button (gtk-check-button-new-with-label check-button-label)))
		 (gtk-box-pack-start vbox c-button #f #f 2)
		 c-button)
	       #f))) ; the check-button when we don't want to see it

      (gtk-box-pack-start vbox h-sep #t #f 2)
      (gtk-box-pack-start vbox hbox-buttons #f #f 5)
      (gtk-box-pack-start hbox-buttons ok-button #t #f 5)
      (gtk-box-pack-start hbox-buttons cancel-button #f #f 5)
      (gtk-box-pack-start hbox-for-entry entry-label #f #f 4)
      (gtk-box-pack-start hbox-for-entry entry #t #t 4)
      (gtk-entry-set-text entry defaut-entry-text)
      
      (gtk-option-menu-set-menu option-menu menu)
      
      ;; button callbacks:
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    ;; what is the molecule number of the option menu?
			    (let ((active-mol-no (get-option-menu-active-molecule 
						  option-menu
						  model-mol-list)))

			      (if (number? active-mol-no)
				  (begin
				    (let* ((text (gtk-entry-get-text entry))
					   (button-state (if check-button
							     (gtk-toggle-button-get-active check-button)
							     'no-check-button)))

				      (let ((func-return-value
					     (if (eq? button-state 'no-check-button)
						 (callback-function active-mol-no text)
						 (callback-function active-mol-no text button-state))))

					(format #t "func-return-value: ~s~%" func-return-value)

					(if (null? always-dismiss-on-ok-clicked?) ;; default
					    (gtk-widget-destroy window)
					    (if func-return-value
						(gtk-widget-destroy window)
						(begin
						  ;; (format #t "not going anywhere\n")
						  #t))))))))))

      (gtk-signal-connect cancel-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window)))

      (gtk-widget-show-all window))))

;; works with below function generic-chooser-entry-and-file-selector
;;
(define *generic-chooser-entry-and-file-selector-file-entry-default-text* "")

;; Create a window
;; 
;; Return a pair of widgets, a molecule chooser and an entry.  The
;; callback-function is a function that takes a molecule number and 2
;; text strings (e.g chain-id and file-name)
;; 
;; chooser-filter is typically valid-map-molecule? or valid-model-molecule?
;; 
;; If file-entry-default-text is passed, then *generic-chooser-entry-and-file-selector-file-entry-default-text* is set
;; to be the contents of the file selection entry on "OK" button click.
;;
(define (generic-chooser-entry-and-file-selector chooser-label chooser-filter entry-hint-text default-entry-text file-selector-hint callback-function . file-entry-default-text)

  (let* ((window (gtk-window-new 'toplevel))
	 (label (gtk-label-new chooser-label))
	 (vbox (gtk-vbox-new #f 2))
	 (hbox-for-entry (gtk-hbox-new #f 0))
	 (entry (gtk-entry-new))
	 (entry-label (gtk-label-new entry-hint-text))
	 (hbox-buttons (gtk-hbox-new #t 2))
	 (menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new))
	 (ok-button (gtk-button-new-with-label "  OK  "))
	 (cancel-button (gtk-button-new-with-label " Cancel "))
	 (h-sep (gtk-hseparator-new))
	 (model-mol-list (fill-option-menu-with-mol-options menu chooser-filter)))

    ;; (fill-option-menu-with-coordinates-mol-options menu)))
    
    (gtk-window-set-default-size window 400 100)
    (gtk-container-add window vbox)
    (gtk-box-pack-start vbox label #f #f 5)
    (gtk-box-pack-start vbox option-menu #t #t 2)
    (gtk-box-pack-start vbox hbox-for-entry #f #f 5)
    (gtk-box-pack-start hbox-buttons ok-button #t #f 5)
    (gtk-box-pack-start hbox-buttons cancel-button #t #f 5)
    (gtk-box-pack-start hbox-for-entry entry-label #t #f 5)
    (gtk-box-pack-start hbox-for-entry entry #t #t 4)
    (gtk-entry-set-text entry default-entry-text)

    (let ((file-sel-entry (file-selector-entry vbox file-selector-hint)))
      (gtk-box-pack-start vbox h-sep #t #f 2)
      (gtk-box-pack-start vbox hbox-buttons #f #f 5)
    
      (gtk-option-menu-set-menu option-menu menu)

      (if (not (null? file-entry-default-text))
	 (let ((file-name (car file-entry-default-text)))
            (gtk-entry-set-text file-sel-entry file-name)))
      
      ;; button callbacks:
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    ;; what is the molecule number of the option menu?
			    (let ((active-mol-no (get-option-menu-active-molecule 
						  option-menu
						  model-mol-list)))
			      
			      (if (number? active-mol-no)
				  (begin
				    (let ((text (gtk-entry-get-text entry))
					  (file-sel-text (gtk-entry-get-text file-sel-entry)))
                                      (if (not (null? file-entry-default-text))
                                          (set! *generic-chooser-entry-and-file-selector-file-entry-default-text* file-sel-text))
				      (callback-function active-mol-no text file-sel-text)))))

			    (gtk-widget-destroy window)))
      
      (gtk-signal-connect cancel-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window)))
      
      (gtk-widget-show-all window))))

;; Create a window.
;; 
;; Return a pair of widgets, a molecule chooser and an entry.
;; callback-function is a function that takes a molecule number and a
;; file-name
;; 
;; chooser-filter is typically valid-map-molecule? or valid-model-molecule?
;; 
(define (generic-chooser-and-file-selector chooser-label chooser-filter file-selector-hint default-file-name callback-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (label (gtk-label-new chooser-label))
	 (vbox (gtk-vbox-new #f 2))
	 (hbox-for-entry (gtk-hbox-new #f 0))
	 (hbox-buttons (gtk-hbox-new #t 2))
	 (menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new))
	 (ok-button (gtk-button-new-with-label "  OK  "))
	 (cancel-button (gtk-button-new-with-label " Cancel "))
	 (h-sep (gtk-hseparator-new))
	 (model-mol-list (fill-option-menu-with-mol-options menu chooser-filter)))
    
    (gtk-window-set-default-size window 400 100)
    (gtk-container-add window vbox)
    (gtk-box-pack-start vbox label #f #f 5)
    (gtk-box-pack-start vbox option-menu #t #t 2)
    (gtk-box-pack-start hbox-buttons ok-button #t #f 5)
    (gtk-box-pack-start hbox-buttons cancel-button #t #f 5)

    (let ((file-sel-entry (file-selector-entry vbox file-selector-hint)))
      (gtk-box-pack-start vbox h-sep #t #f 2)
      (gtk-box-pack-start vbox hbox-buttons #f #f 5)
    
      (gtk-option-menu-set-menu option-menu menu)
      (gtk-entry-set-text file-sel-entry default-file-name)

      (let ((go-func
	     (lambda ()
	       ;; what is the molecule number of the option menu?
	       (let ((active-mol-no (get-option-menu-active-molecule 
				     option-menu
				     model-mol-list)))
		 
		 (if (not (number? active-mol-no))
		     (begin 
		       (format #t "Ooops ~s is not a valid molecule~%" active-mol-no))
		     (begin
		       (let ((file-sel-text (gtk-entry-get-text file-sel-entry)))
			 (callback-function active-mol-no file-sel-text))))))))

	;; button callbacks:
	(gtk-signal-connect ok-button "clicked"
			    (lambda args
			      (go-func)
			      (gtk-widget-destroy window)))

	(gtk-signal-connect file-sel-entry "key-press-event"
			    (lambda (event)
			      (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
				  (begin
				    (go-func)
				    (gtk-widget-destroy window)))
			      #f))
				    
	
	(gtk-signal-connect cancel-button "clicked"
			    (lambda args
			      (gtk-widget-destroy window)))
	
	(gtk-widget-show-all window)))))


;; If a menu with label menu-label is not found in the coot main
;; menubar, then create it and return it. 
;; If it does exist, simply return it.
;;  
(define (coot-menubar-menu menu-label)
  
  (define menu-bar-label-list
    (lambda () 
      (map 
       (lambda (menu-child)
	 (let* ((ac-lab-ls (gtk-container-children menu-child))
		;; ac-lab-ls is a GtkAccelLabel in a list
		(ac-lab (car ac-lab-ls))
		;; ac-lab is a simple GtkAccelLabel
		(lab (gtk-label-get ac-lab)))
	   (cons lab (gtk-menu-item-submenu menu-child)))) ; improper list
       (gtk-container-children (coot-main-menubar)))))

  ;; main body
  ;; 
  (let ((found-menu
	 (let f ((ls (menu-bar-label-list)))
	   (cond 
	    ((null? ls) #f)
	    ((string=? menu-label (car (car ls)))
	     (cdr (car ls)))
	    (else 
	     (f (cdr ls)))))))
    (if found-menu
	found-menu
	(let ((menu (gtk-menu-new))
	      (menuitem (gtk-menu-item-new-with-label menu-label)))
	  (gtk-menu-item-set-submenu menuitem menu)
	  (gtk-menu-bar-append (coot-main-menubar) menuitem)
	  (gtk-widget-show menuitem)
	  menu))))

;; test script
;(let ((menu (coot-menubar-menu "File")))
;  (add-simple-coot-menu-menuitem
;   menu "XXX addition..."
;   (lambda ()
;     (molecule-chooser-gui "Find and Fill residues with missing atoms"
;			   (lambda (imol)
;			     (fill-partial-residues imol))))))


;; Given that we have a menu (e.g. one called "Extensions") provide a
;; cleaner interface to adding something to it:
;; 
;; activate-function is a thunk.
;; 
(define (add-simple-coot-menu-menuitem menu menu-item-label activate-function)

  (let ((submenu (gtk-menu-new))
	(sub-menuitem (gtk-menu-item-new-with-label menu-item-label)))
    
    (gtk-menu-append menu sub-menuitem)
    (gtk-widget-show sub-menuitem)
    
    (gtk-signal-connect sub-menuitem "activate"
			activate-function)))



;;; Make an interesting things GUI for residues of molecule number
;;; imol that have alternate conformations.
;;;
(define (alt-confs-gui imol)

  (interesting-residues-gui imol "Residues with Alt Confs"
			    (residues-with-alt-confs imol)))

;; Make an interesting things GUI for residues with missing atoms
;; 
(define (missing-atoms-gui imol)

  (interesting-residues-gui imol "Residues with missing atoms"
			    (map (lambda (v) (cons #t v))
				 (missing-atom-info imol))))


;;; Make an interesting things GUI for residues of molecule number
;;; imol for the given imol.   A generalization of alt-confs gui
;;;
(define (interesting-residues-gui imol title interesting-residues)

  (if (valid-model-molecule? imol)

      (let* ((residues interesting-residues)
	     (centre-atoms (map (lambda (spec)
				  (if (list? spec)
				      (apply residue-spec->atom-for-centre
					     (cons imol (cdr spec)))
				      #f))
				residues)))
	(interesting-things-gui
	 title
	 (map (lambda (residue-cpmd centre-atom)
		;; see (from above) that centre-atom can be false
		(if centre-atom 
		    (let* ((residue (cdr residue-cpmd))
			   (label (string-append 
				   (car residue) " "
				   (number->string (car (cdr residue))) " "
				   (car (cdr (cdr residue))) " "
				   (car centre-atom) " "
				   (car (cdr centre-atom)))))
		      (append (list label imol) residue centre-atom))

		    ;; protection from centre-atom being false.
		    ;; 
		    ;; make a synthetic baddie
		    (list "[oops - why did this happen?]" 0 0 0)))
	      residues centre-atoms)))))


;; button-list is a list of pairs (improper list) the first item of
;; which is the button label text the second item is a lambda
;; function, what to do when the button is pressed.
;; 
(define (generic-buttons-dialog dialog-name button-list)

  ;; main body	
  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-vbox (gtk-vbox-new #f 0)))
  
    (gtk-window-set-default-size window 250 250)
    (gtk-window-set-title window dialog-name)
    (gtk-container-border-width inside-vbox 4)
    
    (gtk-container-add window outside-vbox)
    (gtk-container-add outside-vbox scrolled-win)
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    
    (let loop ((button-items button-list))
      (cond 
       ((null? button-items) 'done)
       (else 
	(if (pair? (car button-items))
	    (let ((button-label (car (car button-items)))
		  (action (cdr (car button-items))))
	      
	      (let ((button (gtk-button-new-with-label button-label)))
		(gtk-box-pack-start inside-vbox button #f #f 2)
		(gtk-signal-connect button "clicked" action)
		(gtk-widget-show button))))
	(loop (cdr button-items)))))

    (gtk-container-border-width outside-vbox 4)
    (let ((ok-button (gtk-button-new-with-label "  OK  ")))
      (gtk-box-pack-start outside-vbox ok-button #f #f 6)
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window))))

    (gtk-widget-show-all window)))
 

;; Generic interesting things gui: user passes a function that takes 4
;; args: the chain-id, resno, inscode and residue-serial-number
;; (should it be needed) and returns either #f or something
;; interesting (e.g. a label/value).  It is the residue-test-func of
;; the residue-matching-criteria function.
;; 
(define (generic-interesting-things imol gui-title-string residue-test-func)
  
  (if (valid-model-molecule? imol)

      (let* ((interesting-residues (residues-matching-criteria imol residue-test-func))
	     (centre-atoms (map (lambda (spec)
				  (if (list? spec)
				      (apply residue-spec->atom-for-centre 
					     (cons imol (cdr spec)))
				      #f))
				interesting-residues)))
	(interesting-things-gui
	 gui-title-string
	 (map (lambda (interesting-residue centre-atom)
		(if (not (list? centre-atom))
		    "Atom in residue name failure"
		    (let ((label (string-append (list-ref interesting-residue 0) " "
						(list-ref interesting-residue 1) " "
						(number->string 
						 (list-ref interesting-residue 2)) " "
						(list-ref interesting-residue 3)
						(car centre-atom) " "
						(car (cdr centre-atom)))))
		      (append (list label imol) (cdr interesting-residue) centre-atom))))
	      interesting-residues centre-atoms)))))

;; A gui that makes a generic number chooser the go function is a
;; lambda function that takes the value of the active menu item - as a
;; number.
;; 
(define (generic-number-chooser number-list default-option-value hint-text go-button-label go-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox1 (gtk-hbox-new #f 0))
	 (hbox2 (gtk-hbox-new #t 0)) ;; for Go and Cancel
	 (function-label (gtk-label-new hint-text))
	 (h-sep (gtk-hseparator-new))
	 (go-button (gtk-button-new-with-label go-button-label))
	 (cancel-button (gtk-button-new-with-label "  Cancel  "))
	 (menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new)))

    (fill-option-menu-with-number-options menu number-list default-option-value)

    (gtk-box-pack-start vbox hbox1 #t #f 0)
    (gtk-box-pack-start vbox function-label #f 0)
    (gtk-box-pack-start vbox option-menu #t 0)
    (gtk-box-pack-start vbox h-sep)
    (gtk-box-pack-start vbox hbox2 #f #f 0)
    (gtk-box-pack-start hbox2 go-button #f #f 6)
    (gtk-box-pack-start hbox2 cancel-button #f #f 6)
    (gtk-container-add window vbox)
    (gtk-container-border-width vbox 6)
    (gtk-container-border-width hbox1 6)
    (gtk-container-border-width hbox2 6)
    (gtk-signal-connect go-button "clicked" 
			(lambda () 
			  (let ((active-number
				 (get-option-menu-active-item
				  option-menu number-list)))
			    (go-function active-number)
			    (gtk-widget-destroy window))))
    (gtk-signal-connect cancel-button "clicked" 
			(lambda()
			  (gtk-widget-destroy window)))

    (gtk-option-menu-set-menu option-menu menu)

    (gtk-widget-show-all window)))

;; vbox is the vbox to which this compound widget should be added.
;; button-press-func is the lambda function called on pressing return
;; or the button, which takes one argument, the entry.
;;
;; Add this widget to the vbox here.
;;
(define (entry+do-button vbox hint-text button-label button-press-func)

  (let ((hbox (gtk-hbox-new #f 0))
	(entry (gtk-entry-new))
	(button (gtk-button-new-with-label button-label))
	(label (gtk-label-new hint-text)))
    
    (gtk-box-pack-start hbox label #f #f 2)
    (gtk-box-pack-start hbox entry #t #t 2)
    (gtk-box-pack-start hbox button #f #f 2)
    (gtk-signal-connect button "clicked" 
			(lambda ()
			  (button-press-func entry)))

    (gtk-widget-show label)
    (gtk-widget-show button)
    (gtk-widget-show entry)
    (gtk-widget-show hbox)
    (gtk-box-pack-start vbox hbox #t #f)
    entry)) ; return the entry so that we can ask what it says elsewhere

;; pack a hint text and a molecule chooser option menu into the given vbox.
;; 
;; return the option-menu and model molecule list:
(define (generic-molecule-chooser hbox hint-text)

  (let* ((menu (gtk-menu-new))
	 (option-menu (gtk-option-menu-new))
	 (label (gtk-label-new hint-text))
	 (model-mol-list (fill-option-menu-with-coordinates-mol-options menu)))
    
    (gtk-box-pack-start hbox label #f #f 2)
    (gtk-box-pack-start hbox option-menu #t #t 2)
    
    (gtk-option-menu-set-menu option-menu menu)
    (list option-menu model-mol-list)))
    

;; Return an entry, the widget is inserted into the hbox passed to
;; this function.
;; 
(define (file-selector-entry hbox hint-text)

  (let* ((vbox (gtk-vbox-new #f 0))
	 (entry 
	  (entry+do-button vbox hint-text "  File...  "
			   (lambda (entry)
			     (let* ((fs-window (gtk-file-selection-new "file selection")))
			       (gtk-signal-connect 
				(gtk-file-selection-ok-button fs-window)
				"clicked" 
				(lambda ()
				  (let ((t (gtk-file-selection-get-filename fs-window)))
				    (format #t "~s~%" t)
				    (gtk-entry-set-text entry t)
				    (gtk-widget-destroy fs-window))))
			       (gtk-signal-connect 
				(gtk-file-selection-cancel-button fs-window)
				"clicked" 
				(lambda ()
				  (gtk-widget-destroy fs-window)))
			       (gtk-widget-show fs-window))))))

    (gtk-box-pack-start hbox vbox #f #f 2)
    (gtk-widget-show vbox)
    entry))


;; The gui for the strand placement function
;;
(define place-strand-here-gui
  (lambda ()

    (generic-number-chooser (number-list 4 12) 7
			    " Estimated number of residues in strand "
			    "  Go  "
			    (lambda (n)
			      (place-strand-here n 15)))))

;; Cootaneer gui
(define (cootaneer-gui imol)

  (define (add-text-to-text-box text-box description)
    (gtk-text-insert text-box #f "black" "#c0e6c0" description -1))

  ;; return the (entry . text-box)
  ;; 
  (define (entry-text-pair-frame seq-info)

    (let ((frame (gtk-frame-new ""))
	  (vbox (gtk-vbox-new #f 3))
	  (entry (gtk-entry-new))
	  (text-box (gtk-text-new #f #f))
	  (chain-id-label (gtk-label-new "Chain ID"))
	  (sequence-label (gtk-label-new "Sequence")))

      (gtk-container-add frame vbox)
      (gtk-box-pack-start vbox chain-id-label #f #f 2)
      (gtk-box-pack-start vbox entry #f #f 2)
      (gtk-box-pack-start vbox sequence-label #f #f 2)
      (gtk-box-pack-start vbox text-box #f #f 2)
      (add-text-to-text-box text-box (cdr seq-info))
      (gtk-entry-set-text entry (car seq-info))
      (list frame entry text-box)))
      
		 
  ;; main body
  (let ((imol-map (imol-refinement-map)))
    (if (= imol-map -1)
	(show-select-map-dialog)
    
	(let* ((window (gtk-window-new 'toplevel))
	       (outside-vbox (gtk-vbox-new #f 2))
	       (inside-vbox (gtk-vbox-new #f 2))
	       (h-sep (gtk-hseparator-new))
	       (buttons-hbox (gtk-hbox-new #t 2))
	       (go-button (gtk-button-new-with-label "  Dock Sequence  "))
	       (cancel-button (gtk-button-new-with-label "  Cancel  ")))

	  (let ((seq-info-ls (sequence-info imol)))
	    
	    (if (not seq-info-ls)
		(begin
		  (let ((s (format #f "No sequence assigned for molecule number ~s" imol)))
		    (format #t "~s~%"s)
		    (info-dialog s)))
		(begin
		  (map (lambda (seq-info)
			 (let ((seq-widgets (entry-text-pair-frame seq-info)))
			   (gtk-box-pack-start inside-vbox (car seq-widgets) #f #f 2)))
		       seq-info-ls)

		  (gtk-box-pack-start outside-vbox inside-vbox #f #f 2)
		  (gtk-box-pack-start outside-vbox h-sep #f #f 2)
		  (gtk-box-pack-start outside-vbox buttons-hbox #t #f 2)
		  (gtk-box-pack-start buttons-hbox go-button #t #f 6)
		  (gtk-box-pack-start buttons-hbox cancel-button #t #f 6)

		  (gtk-signal-connect cancel-button "clicked"
				      (lambda ()
					(gtk-widget-destroy window)))
		  
		  (gtk-signal-connect 
		   go-button "clicked"
		   (lambda ()

		     ;; no active atom won't do.  We need
		     ;; to find the nearest atom in imol to (rotation-centre).
		     ;;
		     ;; if it is too far away, give a
		     ;; warning and do't do anything.

		     (let ((n-atom (closest-atom imol)))
		       (if n-atom
			   (let ((imol     (list-ref n-atom 0))
				 (chain-id (list-ref n-atom 1))
				 (resno    (list-ref n-atom 2))
				 (inscode  (list-ref n-atom 3))
				 (at-name  (list-ref n-atom 4))
				 (alt-conf (list-ref n-atom 5)))
			     (let ((cootaneer-results
				    (cootaneer imol-map imol (list chain-id resno inscode 
								   at-name alt-conf))))
			       (format #t "Cootaneering status: ~s~%" cootaneer-results)
			       (if (= cootaneer-results 0)
				   (let ((s (string-append
					     "Insufficiently confident in alignment to make a fit."
					     "\n"
					     "Perhaps you could improve or extend this fragment.")))
				     (gtk-widget-destroy window)
				     (info-dialog s)))))))))
		     
		  (gtk-container-add window outside-vbox)
		  (gtk-widget-show-all window))))))))
    

;; The gui for saving views
(define view-saver-gui
  (lambda ()

    (define (local-view-name)
      (let loop ((view-count 0))
        (let ((str (string-append "View"
                                  (if (> view-count 0)
                                      (string-append "-" (number->string view-count))
                                      ""))))
          ;; now is a view already called str?
          (let iloop ((jview 0))
            (let ((jview-name (view-name jview)))

              (cond
               ((>= jview (n-views)) str)
               ((eq? #f jview-name) str)
               ((string=? str jview-name) (loop (+ view-count 1)))
               (else
                (iloop (+ jview 1)))))))))

    ;; main line
    ;; 
    (let ((view-name (local-view-name)))
      (generic-single-entry "View Name: " view-name " Add View " 
                            (lambda (text)
                              (let ((new-view-number (add-view-here text)))
                                (add-view-to-views-panel text new-view-number)))))))
				

;; 
(define (add-view-to-views-panel view-name view-number)

  (if *views-dialog-vbox*
      (let ((button (gtk-button-new-with-label view-name)))
        (gtk-signal-connect button "clicked" (lambda () (go-to-view-number view-number 0)))
        (gtk-box-pack-start *views-dialog-vbox* button #f #f 2)
        (gtk-widget-show button))))


;; geometry is an improper list of ints.
;; 
;; return a list of (h-box-buttons window)
;;
;; a button is a list of (label callback-thunk text-description)
;; where text-description is an option arg (can be omitted)
;;
(define (dialog-box-of-buttons window-name geometry buttons close-button-label)
  (dialog-box-of-buttons-with-check-button window-name geometry 
					   buttons close-button-label
					   #f #f #f))

;; geometry is an improper list of ints.
;; 
;; return a list of (h-box-buttons window)
;;
;; a button is a list of (label callback-thunk text-description)
;;
;; If check-button-label is #f, don't make one, otherwise create with with
;; the given label and "on" state.
;; 
(define (dialog-box-of-buttons-with-check-button window-name geometry buttons close-button-label check-button-label check-button-func check-button-is-initially-on-flag)

  (define (add-text-to-text-widget text-box description)
    (gtk-text-insert text-box #f "black" "#c0e6c0" description -1))

  ;; main line
  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (h-sep (gtk-hseparator-new))
	 (inside-vbox (gtk-vbox-new #f 0)))
    
    (gtk-window-set-default-size window (car geometry) (cdr geometry))
    (gtk-window-set-title window window-name)
    (gtk-container-border-width inside-vbox 2)
    (gtk-container-add window outside-vbox)

    (if (string? check-button-label)
	(begin
	  (let ((check-button (gtk-check-button-new-with-label 
			       check-button-label)))
	    (gtk-signal-connect check-button "toggled" 
				(lambda ()
				  (check-button-func check-button inside-vbox)))
	    (if check-button-is-initially-on-flag
		(gtk-toggle-button-set-active check-button #t))
	    (gtk-box-pack-start outside-vbox check-button #f #f 2))))

    (gtk-box-pack-start outside-vbox scrolled-win #t #t 0) ; expand fill padding
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

    (map (lambda (button-info)
	   (add-button-info-to-box-of-buttons-vbox button-info inside-vbox))
	 buttons)

    (gtk-container-border-width outside-vbox 2)
    (gtk-box-pack-start outside-vbox h-sep #f #f 2)
    (let ((ok-button (gtk-button-new-with-label close-button-label)))
      (gtk-box-pack-end outside-vbox ok-button #f #f 0)

      ;; Note to self: the setting of inside-vbox should not be done
      ;; in this generic dialog-box-of-buttons, but passed as an
      ;; argment to the function (e.g. destroy-window-extra-func)
      ;; 
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window)
			    (set! inside-vbox #f)))) ; redundant?
      (gtk-signal-connect window "destroy"
			  (lambda args
			    (set! inside-vbox #f)))
      (gtk-widget-show-all window)
      (list inside-vbox window)))


(define (dialog-box-of-buttons-from-specs window-name geometry imol specs)

  (let ((buttons (map (lambda(spec)
                        (let ((label (residue-spec->string spec))
                              (cbf (lambda ()
                                     (set-go-to-atom-molecule imol)
                                     (set-go-to-atom-chain-residue-atom-name
                                      (residue-spec->chain-id spec)
                                      (residue-spec->res-no spec) " C  "))))
                          (list label cbf)))
                      specs)))
    (dialog-box-of-buttons window-name geometry buttons " Close ")))


;; This is exported outside of the box-of-buttons gui because the
;; clear-and-add-back function (e.g. from using the check button)
;; also needs to add buttons - let's not duplicate that code.
;;
(define (add-button-info-to-box-of-buttons-vbox button-info vbox)

  (define (add-text-to-text-widget text-box description)
    (gtk-text-insert text-box #f "black" "#c0e6c0" description -1))

  (let* ((buton-label (car button-info))
	 (callback (car (cdr button-info)))
	 (description (if (= (length button-info) 2)
			  #f ; it doesn't have one
			  (list-ref button-info 2)))
	 (button (gtk-button-new-with-label buton-label)))
    (gtk-signal-connect button "clicked" callback)

    (if (string? description)
	(let ((text-box (gtk-text-new #f #f)))
	  (add-text-to-text-widget text-box description)
	  (gtk-box-pack-start vbox text-box #f #f 2)

	  (gtk-widget-realize text-box)
	  (gtk-text-thaw text-box)))

    ;; "description" is a bad name!  We run a function and pass the button as the argument
    (if (procedure? description)
	(description button))

    (gtk-box-pack-start vbox button #f #f 2)
    (gtk-widget-show button)))


;; geometry is an improper list of ints
;; buttons is a list of: (list (list button-1-label button-1-action
;;                                   button-2-label button-2-action))
;; The button-1-action function takes as an argument the imol
;; The button-2-action function takes as an argument the imol
;; 
(define (dialog-box-of-pairs-of-buttons imol window-name geometry buttons close-button-label)

  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-vbox (gtk-vbox-new #f 0)))
    
    (gtk-window-set-default-size window (car geometry) (cdr geometry))
    (gtk-window-set-title window window-name)
    (gtk-container-border-width inside-vbox 2)
    (gtk-container-add window outside-vbox)
    (gtk-box-pack-start outside-vbox scrolled-win #t #t 0) ; expand fill padding
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

    (map (lambda (buttons-info)
;	   (format #t "buttons-info ~s~%" buttons-info) 
	   (if (list? buttons-info)
	       (let* ((button-label-1 (car buttons-info))
		      (callback-1  (car (cdr buttons-info)))

		      (button-label-2 (car (cdr (cdr buttons-info))))
		      (callback-2  (car (cdr (cdr (cdr buttons-info)))))

		      (button-1 (gtk-button-new-with-label button-label-1))
		      (h-box (gtk-hbox-new #f 2)))

;		 (format #t "button-label-1 ~s~%" button-label-1) 
;		 (format #t "callback-1 ~s~%" callback-1) 
;		 (format #t "buton-label-2 ~s~%" button-label-2) 
;		 (format #t "callback-2 ~s~%" callback-2) 

		 (gtk-signal-connect button-1 "clicked" 
				     (lambda ()
				       (callback-1 imol)))
		 (gtk-box-pack-start h-box button-1 #f #f 2)

		 (if callback-2 
		     (let ((button-2 (gtk-button-new-with-label button-label-2)))
		       (gtk-signal-connect button-2 "clicked" 
					   (lambda ()
					     (callback-2 imol)))
		       (gtk-box-pack-start h-box button-2 #f #f 2)))

		 (gtk-box-pack-start inside-vbox h-box #f #f 2))))
	 buttons)

    (gtk-container-border-width outside-vbox 2)
    (let ((ok-button (gtk-button-new-with-label close-button-label)))
      (gtk-box-pack-end outside-vbox ok-button #f #f 0)
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window))))
    (gtk-widget-show-all window)))

(define *views-dialog-vbox* #f)

;; A gui showing views:
(define views-panel-gui
  (lambda ()

    (let* ((number-of-views (n-views))
           (buttons

            (let loop ((button-number 0)
                       (ls '()))

              (cond
               ((= number-of-views button-number) (reverse ls))
               (else
                (let ((button-label (view-name button-number)))
                  (if button-label
                      (let ((desciption (view-description button-number)))
                        (if (eq? #f desciption)
                            (loop (+ button-number 1)
                                  (cons
                                   (list button-label
                                         (lambda ()
                                           (go-to-view-number button-number 0)))
                                   ls))
                            (loop (+ button-number 1)
                                  (cons
                                   (list button-label
                                         (lambda ()
                                           (go-to-view-number button-number 0))
                                         desciption)
                                   ls))))

                      (loop (+ button-number 1) ls))))))))
      (let ((all-buttons
             (if (> (length buttons) 1)
                 (let ((view-button (list "  Play Views "
                                          (lambda ()
                                            (go-to-first-view 1)
                                            (usleep 100000)
                                            (play-views)))))
                   (reverse (cons view-button (reverse buttons))))
                 buttons)))

        (let ((views-vbox
               (dialog-box-of-buttons " Views " (cons 200 140) all-buttons "  Close  ")))

          (set! *views-dialog-vbox* views-vbox))))))

   
;   (* nudge (apply + (map * (list-ref tvm 0) ori)))
;   (* nudge (apply + (map * (list-ref tvm 1) ori)))
;   (* nudge (apply + (map * (list-ref tvm 2) ori)))))
  

; (display  
;(screen-coords-nudge (transpose-simple-matrix (view-matrix))
;		     1.0
;		     (list 1 0 0)))
;(newline)
     
;; nudge screen centre box.  Useful when Ctrl left-mouse has been
;; taken over by another function.
;; 
(define (nudge-screen-centre-gui)

  (define (screen-coords-nudge tvm nudge ori)
    (map (lambda (e) (* nudge (apply + (map * e ori)))) tvm))

  (let* ((zsc 0.02)
	 (f (lambda (axes)
	      (let ((nudge (* (zoom-factor) zsc))
		    (rc (rotation-centre))
		    (tvm (transpose-simple-matrix (view-matrix))))
		;; new-x new-y new-z:
		;; rc + tvm*nudge
		(apply set-rotation-centre 
		       (map + rc (screen-coords-nudge tvm nudge axes))))))
	 (buttons
	  (list 
	   (list "Nudge +X" (lambda () (f (list  1 0 0))))
	   (list "Nudge -X" (lambda () (f (list -1 0 0))))
	   (list "Nudge +Y" (lambda () (f (list 0  1 0))))
	   (list "Nudge -Y" (lambda () (f (list 0 -1 0))))
	   (list "Nudge +Z" (lambda () (f (list 0 0  1))))
	   (list "Nudge -Z" (lambda () (f (list 0 0 -1)))))))

    (dialog-box-of-buttons "Nudge Screen Centre" 
			   (cons 200 190) buttons "  Close ")))


;; A gui to make a difference map (from arbitrarily gridded maps
;; (that's it's advantage))
;; 
(define (make-difference-map-gui)

  (let ((window (gtk-window-new 'toplevel))
	(diff-map-vbox (gtk-vbox-new #f 2))
	(h-sep (gtk-hseparator-new))
	(title (gtk-label-new "Make a Difference Map:"))
	(ref-label (gtk-label-new "Reference Map:"))
	(sec-label (gtk-label-new "Subtract this map:"))
	(second-map-hbox (gtk-hbox-new #f 2))
	(buttons-hbox (gtk-hbox-new #t 6))
	(option-menu-ref-mol (gtk-option-menu-new))
	(option-menu-sec-mol (gtk-option-menu-new))
	(menu-ref (gtk-menu-new))
	(menu-sec (gtk-menu-new))
	(scale-label (gtk-label-new "Scale"))
	(scale-entry (gtk-entry-new))
	(ok-button (gtk-button-new-with-label "   OK   "))
	(cancel-button (gtk-button-new-with-label " Cancel ")))

    (let ((map-molecule-list-ref (fill-option-menu-with-map-mol-options
				  menu-ref))
	  (map-molecule-list-sec (fill-option-menu-with-map-mol-options
				  menu-sec)))

      (gtk-option-menu-set-menu option-menu-ref-mol menu-ref)
      (gtk-option-menu-set-menu option-menu-sec-mol menu-sec)

      (gtk-container-add window diff-map-vbox)
      (gtk-box-pack-start diff-map-vbox title #f #f 2)
      (gtk-box-pack-start diff-map-vbox ref-label #f #f 2)
      (gtk-box-pack-start diff-map-vbox option-menu-ref-mol #t #t 2)

      (gtk-box-pack-start diff-map-vbox sec-label #f #f 2)
      (gtk-box-pack-start diff-map-vbox second-map-hbox #f #f 2)

      (gtk-box-pack-start second-map-hbox option-menu-sec-mol #t #t 2)
      (gtk-box-pack-start second-map-hbox scale-label #f #f 2)
      (gtk-box-pack-start second-map-hbox scale-entry #f #f 2)

      (gtk-box-pack-start diff-map-vbox h-sep #t #f 2)
      (gtk-box-pack-start diff-map-vbox buttons-hbox #t #f 2)
      (gtk-box-pack-start buttons-hbox ok-button #t #f 2)
      (gtk-box-pack-start buttons-hbox cancel-button #t #f 2)
      (gtk-entry-set-text scale-entry "1.0")

      (gtk-signal-connect ok-button "clicked"
			  (lambda()
			    (format #t "make diff map here~%")
			    (let* ((active-mol-no-ref
				   (get-option-menu-active-molecule 
				    option-menu-ref-mol
				    map-molecule-list-ref))
				   (active-mol-no-sec
				    (get-option-menu-active-molecule 
				     option-menu-sec-mol
				     map-molecule-list-sec))
				   (scale-text (gtk-entry-get-text scale-entry))
				   (scale (string->number scale-text)))

			      (if (not (number? scale))
				  (format #t "can't decode scale ~s~%" scale-text)
				  (difference-map active-mol-no-ref
						  active-mol-no-sec
						  scale))
			      (gtk-widget-destroy window))))
      
      (gtk-signal-connect cancel-button "clicked"
			  (lambda()
			    (gtk-widget-destroy window)))

      (gtk-widget-show-all window))))

;; A GUI to display all the CIS peptides and navigate to them.
;; 
(define (cis-peptides-gui imol)

  (define (get-ca atom-list)
    
    (let loop ((atom-list atom-list))
      (cond
       ((null? atom-list) #f)
       (else 
	(let ((atom-name (car (car (car atom-list)))))
	  (if (string=? atom-name " CA ")
	      (car atom-list)
	      (loop (cdr atom-list))))))))

  
  (let ((cis-peps (cis-peptides imol)))
    
    (if (null? cis-peps)
	(info-dialog "No Cis Peptides found")
	(interesting-things-gui
	 "Cis Peptides:"
	 (map (lambda (cis-pep-spec)

		(let ((r1 (list-ref cis-pep-spec 0))
		      (r2 (list-ref cis-pep-spec 1))
		      (omega (list-ref cis-pep-spec 2)))
		  (let ((atom-list-r1 (apply residue-info 
					     imol (cdr r1)))
			(atom-list-r2 (apply residue-info 
					     imol (cdr r2))))
		    (let ((ca-1 (get-ca atom-list-r1))
			  (ca-2 (get-ca atom-list-r2)))

		      (if (not (and ca-1 ca-2))
			  (list (string-append "Cis Pep (atom failure) " 
					       (car (cdr r1))
					       " "
					       (number->string (car (cdr (cdr r1)))))
				imol chain-id 
				(car (cdr r1))
				(car (cdr (cdr r1)))
				"CA" "")
			  (let ((p-1 (car (cdr (cdr ca-1))))
				(p-2 (car (cdr (cdr ca-2)))))
			    (let* ((pos (map (lambda (x) (/ x 2)) 
					    (map + p-1 p-2)))
				   (tors-s1 (number->string (list-ref cis-pep-spec 2)))
				   (tors-string (if (< (string-length tors-s1) 6)
						    tors-s1
						    (substring tors-s1 0 6)))
				   (mess (string-append 
					  "Cis Pep: "
					  (car (cdr r1)) ; "A"
					  " "
					  (number->string (car (cdr (cdr r1))))
					  " " 
					  (apply residue-name imol (cdr r1))
					  " - "
					  (number->string (car (cdr (cdr r2))))
					  " " 
					  (apply residue-name imol (cdr r2))
					  "   "
					  tors-string)))
							      
			      (cons mess pos))))))))
	      cis-peps)))))


(define (transform-map-using-lsq-matrix-gui)

  ;; atom-sel-type is either 'Reference or 'Moving
  ;; 
  ;; return the (list frame option-menu model-mol-list)
  (define (atom-sel-frame atom-sel-type)
    (let* ((frame (gtk-frame-new (symbol->string atom-sel-type)))
	   (option-menu (gtk-option-menu-new))
	   (menu (gtk-menu-new))
	   (model-mol-list (fill-option-menu-with-coordinates-mol-options menu))
	   (atom-sel-vbox (gtk-vbox-new #f 2))
	   (atom-sel-hbox (gtk-hbox-new #f 2))
	   (chain-id-label (gtk-label-new " Chain ID "))
	   (resno-1-label (gtk-label-new " Resno Start "))
	   (resno-2-label (gtk-label-new " Resno End "))
	   (chain-id-entry (gtk-entry-new))
	   (resno-1-entry (gtk-entry-new))
	   (resno-2-entry (gtk-entry-new)))
      (gtk-container-add frame atom-sel-vbox)
      (gtk-box-pack-start atom-sel-vbox   option-menu  #f #f 2)
      (gtk-box-pack-start atom-sel-vbox atom-sel-hbox  #f #f 2)
      (gtk-box-pack-start atom-sel-hbox chain-id-label #f #f 2)
      (gtk-box-pack-start atom-sel-hbox chain-id-entry #f #f 2)
      (gtk-box-pack-start atom-sel-hbox resno-1-label  #f #f 2)
      (gtk-box-pack-start atom-sel-hbox resno-1-entry  #f #f 2)
      (gtk-box-pack-start atom-sel-hbox resno-2-label  #f #f 2)
      (gtk-box-pack-start atom-sel-hbox resno-2-entry  #f #f 2)
      (gtk-option-menu-set-menu option-menu menu)
      (list frame option-menu model-mol-list chain-id-entry resno-1-entry resno-2-entry)))
	  
  (let* ((window (gtk-window-new 'toplevel))
	 (dialog-name "Map Transformation")
	 (main-vbox (gtk-vbox-new #f 2))
	 (buttons-hbox (gtk-hbox-new #f 2))
	 (cancel-button (gtk-button-new-with-label "  Cancel  "))
	 (ok-button (gtk-button-new-with-label "  Transform  "))
	 (usage (string-append "Note that this will transform the current refinement map "
			       "to around the screen centre"))
	 (usage-label (gtk-label-new usage))
	 (h-sep (gtk-hseparator-new))
	 (frame-info-ref (atom-sel-frame 'Reference))
	 (frame-info-mov (atom-sel-frame 'Moving))
	 (radius-hbox (gtk-hbox-new #f 2))
	 (radius-label (gtk-label-new "  Radius "))
	 (radius-entry (gtk-entry-new)))

    (gtk-box-pack-start radius-hbox radius-label #f #f 2)
    (gtk-box-pack-start radius-hbox radius-entry #f #f 2)

    (gtk-box-pack-start buttons-hbox     ok-button #f #f 4)
    (gtk-box-pack-start buttons-hbox cancel-button #f #f 4)

    (gtk-container-add window main-vbox)
    (gtk-box-pack-start main-vbox (car frame-info-ref) #f #f 2)
    (gtk-box-pack-start main-vbox (car frame-info-mov) #f #f 2)
    (gtk-box-pack-start main-vbox radius-hbox #f #f 2)
    (gtk-box-pack-start main-vbox usage-label #f #f 4)
    (gtk-box-pack-start main-vbox h-sep #f #f 2)
    (gtk-box-pack-start main-vbox buttons-hbox #f #f 6)
    
    (gtk-entry-set-text (list-ref frame-info-ref 3) "A")
    (gtk-entry-set-text (list-ref frame-info-ref 4) "1")
    (gtk-entry-set-text (list-ref frame-info-ref 5) "10")
    (gtk-entry-set-text (list-ref frame-info-mov 3) "B")
    (gtk-entry-set-text (list-ref frame-info-mov 4) "1")
    (gtk-entry-set-text (list-ref frame-info-mov 5) "10")

    (gtk-entry-set-text radius-entry "8.0")

    (gtk-signal-connect cancel-button "clicked"
			(lambda ()
			  (gtk-widget-destroy window)))

    (gtk-signal-connect ok-button "clicked"
			(lambda ()

			  (let ((active-mol-ref (apply get-option-menu-active-molecule 
						       (cdr (first-n 3 frame-info-ref))))
				(active-mol-mov (apply get-option-menu-active-molecule 
						       (cdr (first-n 3 frame-info-mov))))
				(chain-id-ref     (gtk-entry-get-text (list-ref frame-info-ref 3)))
				(resno-1-ref-text (gtk-entry-get-text (list-ref frame-info-ref 4)))
				(resno-2-ref-text (gtk-entry-get-text (list-ref frame-info-ref 5)))

				(chain-id-mov     (gtk-entry-get-text (list-ref frame-info-mov 3)))
				(resno-1-mov-text (gtk-entry-get-text (list-ref frame-info-mov 4)))
				(resno-2-mov-text (gtk-entry-get-text (list-ref frame-info-mov 5)))
				
				(radius-text (gtk-entry-get-text radius-entry)))

			    (let ((imol-map (imol-refinement-map))
				  (resno-1-ref (string->number resno-1-ref-text))
				  (resno-2-ref (string->number resno-2-ref-text))
				  (resno-1-mov (string->number resno-1-mov-text))
				  (resno-2-mov (string->number resno-2-mov-text))
				  (radius (string->number radius-text)))

			      (if (and (number? resno-1-ref)
				       (number? resno-2-ref)
				       (number? resno-1-mov)
				       (number? resno-2-mov)
				       (number? radius))

				  (if (not (valid-map-molecule? imol-map))
				      (format #t "Must set the refinement map~%")

				      (let ((imol-copy (copy-molecule active-mol-mov)))
					(let ((new-map-number 
					       (transform-map-using-lsq-matrix 
						active-mol-ref chain-id-ref resno-1-ref resno-2-ref
						imol-copy chain-id-mov resno-1-mov resno-2-mov
						imol-map (rotation-centre) radius)))
					  
					  (set-molecule-name imol-copy
							     (string-append "Transformed copy of " 
									    (strip-path
									     (molecule-name active-mol-mov))))
					  (set-molecule-name new-map-number
							     (string-append
							      "Transformed map: from map "
							      (number->string imol-map)
							      " by matrix that created coords " 
							      (number->string imol-copy)))
					  (set-mol-active    imol-copy 0)
					  (set-mol-displayed imol-copy 0)))))
			  
			      (gtk-widget-destroy window)))))
    
    (gtk-widget-show-all window)
    (if (not (valid-map-molecule? (imol-refinement-map)))
	(show-select-map-dialog))))



		      

(define (ncs-ligand-gui) 
  
  (let ((window (gtk-window-new 'toplevel))
	(ncs-ligands-vbox (gtk-vbox-new #f 2))
	(title (gtk-label-new "Find NCS-Related Ligands"))
	(ref-label (gtk-label-new "Protein with NCS"))
	(ref-chain-hbox (gtk-hbox-new #f 2))
	(chain-id-ref-label (gtk-label-new "NCS Master Chain"))
	(chain-id-ref-entry (gtk-entry-new))
	(lig-label (gtk-label-new "Molecule containing ligand"))
	(specs-hbox (gtk-hbox-new #f 2))
	(h-sep (gtk-hseparator-new))
	(buttons-hbox (gtk-hbox-new #f 6))
	(chain-id-lig-label (gtk-label-new "Chain ID: "))
	(resno-start-label (gtk-label-new " Residue Number "))
	(to-label (gtk-label-new "  to  "))
	(chain-id-lig-entry (gtk-entry-new))
	(resno-start-entry (gtk-entry-new))
	(resno-end-entry (gtk-entry-new))
	(ok-button (gtk-button-new-with-label "   Find Candidate Positions  "))
	(cancel-button (gtk-button-new-with-label "    Cancel    "))
	(menu-ref (gtk-menu-new))
	(menu-lig (gtk-menu-new))
	(option-menu-ref-mol (gtk-option-menu-new))
	(option-menu-lig-mol (gtk-option-menu-new)))
	
    (let ((molecule-list-ref (fill-option-menu-with-coordinates-mol-options
				  menu-ref))
	  (molecule-list-lig (fill-option-menu-with-coordinates-mol-options
				  menu-lig)))

      (gtk-option-menu-set-menu option-menu-ref-mol menu-ref)
      (gtk-option-menu-set-menu option-menu-lig-mol menu-lig)

      (gtk-container-add window ncs-ligands-vbox)
      (gtk-box-pack-start ncs-ligands-vbox title #f #f 6)
      (gtk-box-pack-start ncs-ligands-vbox ref-label #f #f 2)
      (gtk-box-pack-start ncs-ligands-vbox option-menu-ref-mol #t #f 2)
      (gtk-box-pack-start ncs-ligands-vbox ref-chain-hbox #f #f 2)
      (gtk-box-pack-start ncs-ligands-vbox lig-label #f #f 2)
      (gtk-box-pack-start ncs-ligands-vbox option-menu-lig-mol #t #f 2)
      (gtk-box-pack-start ncs-ligands-vbox specs-hbox #f #f 2)
      (gtk-box-pack-start ncs-ligands-vbox h-sep #f #f 2)
      (gtk-box-pack-start ncs-ligands-vbox buttons-hbox #f #f 6)

      (gtk-box-pack-start buttons-hbox     ok-button #t #f 4)
      (gtk-box-pack-start buttons-hbox cancel-button #t #f 4)

      (gtk-box-pack-start ref-chain-hbox chain-id-ref-label #f #f 2)
      (gtk-box-pack-start ref-chain-hbox chain-id-ref-entry #f #f 2)
      
      (gtk-box-pack-start specs-hbox chain-id-lig-label #f #f 2)
      (gtk-box-pack-start specs-hbox chain-id-lig-entry #f #f 2)
      (gtk-box-pack-start specs-hbox resno-start-label #f #f 2)
      (gtk-box-pack-start specs-hbox resno-start-entry #f #f 2)
      (gtk-box-pack-start specs-hbox to-label #f #f 2)
      (gtk-box-pack-start specs-hbox resno-end-entry #f #f 2)
      (gtk-box-pack-start specs-hbox (gtk-label-new " ") #f #f 2) ; neatness

      (gtk-widget-set-usize chain-id-lig-entry 32 -1)
      (gtk-widget-set-usize chain-id-ref-entry 32 -1)
      (gtk-widget-set-usize resno-start-entry  50 -1)
      (gtk-widget-set-usize resno-end-entry    50 -1)
      (gtk-entry-set-text chain-id-ref-entry "A")
      (gtk-entry-set-text chain-id-lig-entry "A")
      (gtk-entry-set-text resno-start-entry "1")

      (gtk-tooltips-set-tip (gtk-tooltips-new)
			    chain-id-ref-entry 
			    (string-append
			    "\"A\" is a reasonable guess at the NCS master chain id.  "
			    "If your ligand (specified below) is NOT bound to the protein's "
			    "\"A\" chain, then you will need to change this chain and also "
			    "make sure that the master molecule is specified appropriately "
			    "in the Draw->NCS Ghost Control window.")
			    "")
      (gtk-tooltips-set-tip (gtk-tooltips-new)
			    resno-end-entry 
			    "Leave blank for a single residue" "")
      
      (gtk-signal-connect ok-button "clicked"
			  (lambda ()
			    (format #t "ncs ligand function here~%")
			    (let* ((active-mol-no-ref
				    (get-option-menu-active-molecule 
				     option-menu-ref-mol
				     molecule-list-ref))
				   (active-mol-no-lig
				    (get-option-menu-active-molecule 
				     option-menu-lig-mol
				     molecule-list-lig))
				   (chain-id-lig (gtk-entry-get-text chain-id-lig-entry))
				   (chain-id-ref (gtk-entry-get-text chain-id-ref-entry))
				   (resno-start
				    (string->number (gtk-entry-get-text resno-start-entry)))
				   (resno-end-t (string->number (gtk-entry-get-text resno-end-entry)))
				   (resno-end (if (number? resno-end-t)
						  resno-end-t
						  resno-start)))

			      ;; (format #t "resnos: ~s ~s~%" resno-start resno-end)
					
			      (if (number? resno-start)
				  (if (number? resno-end)
				      (begin 
					(make-ncs-ghosts-maybe active-mol-no-ref)
					(format #t "ncs ligands with ~s %" 
						(list active-mol-no-ref
						      chain-id-ref
						      active-mol-no-lig
						      chain-id-lig
						      resno-start
						      resno-end))
					
					(ncs-ligand active-mol-no-ref
						    chain-id-ref
						    active-mol-no-lig
						    chain-id-lig
						    resno-start
						    resno-end)))))
			    (gtk-widget-destroy window)))

      (gtk-signal-connect cancel-button "clicked"
			  (lambda ()
			    (gtk-widget-destroy window)))

      (gtk-widget-show-all window))))


(define *ncs-jumping-time-step* 500)

(define (ncs-jumping-gui)

  ;; main body
  (let* ((window (gtk-window-new 'toplevel))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-hbox (gtk-hbox-new #f 2))
	 (cancel-hbox (gtk-vbox-new #f 2))
	 (h-sep (gtk-hseparator-new))
	 (jump-start-button (gtk-button-new-with-label "NCS Jump Start"))
	 (jump-stop-button  (gtk-button-new-with-label "Stop"))
	 (cancel-button     (gtk-button-new-with-label "Cancel"))
	 (ms-step *ncs-jumping-time-step*)
	 (timeout-function-token #f))

    (gtk-window-set-title window "Auto NCS Jumping")
    (gtk-container-add window outside-vbox)
    (gtk-box-pack-start outside-vbox inside-hbox #f #f 2)
    (gtk-box-pack-start outside-vbox h-sep #f #f 2)
    (gtk-box-pack-start outside-vbox cancel-hbox #f #f 2)
    (gtk-box-pack-start inside-hbox jump-start-button  #f #f 2)
    (gtk-box-pack-start inside-hbox jump-stop-button   #f #f 2)
    (gtk-box-pack-start cancel-hbox cancel-button #f #f 2)

    (gtk-signal-connect jump-start-button "clicked"
			(lambda ()
			  (if (not (number? timeout-function-token))
			      (set! timeout-function-token (gtk-timeout-add ms-step (lambda() (skip-to-next-ncs-chain 'forward)))))))

    (gtk-signal-connect jump-stop-button "clicked"
			(lambda ()
			  (if (number? timeout-function-token)
			      (gtk-timeout-remove timeout-function-token))
			  (set! timeout-function-token #f)))

    (gtk-signal-connect cancel-button "clicked"
			(lambda ()
			  (if (number? timeout-function-token)
			      (gtk-timeout-remove timeout-function-token))
			  (gtk-widget-destroy window)))

    (gtk-widget-show-all window)))


;; GUI for ligand superpositioning by graph matching
;;
(define (superpose-ligand-gui) 
  
  (let ((window (gtk-window-new 'toplevel))
	(ligands-vbox (gtk-vbox-new #f 2))
	(title (gtk-label-new "Superpose Ligands"))
	;; I (BL) would like to use generic-molecule-chooser but in 
	;; scheme we may run into packing problems, so use Paul/scheme style
	(ref-label (gtk-label-new "Model with reference ligand"))
	(ref-chain-hbox (gtk-hbox-new #f 2))
	(chain-id-ref-label (gtk-label-new "Ligand Chain ID: "))
	(chain-id-ref-entry (gtk-entry-new))
	(resno-ref-label (gtk-label-new " Residue Number "))
	(resno-ref-entry (gtk-entry-new))

	(mov-label (gtk-label-new "Model with moving ligand"))
	(mov-chain-hbox (gtk-hbox-new #f 2))
	(chain-id-mov-label (gtk-label-new "Ligand Chain ID: "))
	(chain-id-mov-entry (gtk-entry-new))
	(resno-mov-label (gtk-label-new " Residue Number "))
	(resno-mov-entry (gtk-entry-new))
	(h-sep (gtk-hseparator-new))
	(buttons-hbox (gtk-hbox-new #f 6))
	(ok-button (gtk-button-new-with-label "   Superpose Ligands  "))
	(cancel-button (gtk-button-new-with-label "    Cancel    "))
	(menu-ref (gtk-menu-new))
	(menu-mov (gtk-menu-new))
	(option-menu-ref-mol (gtk-option-menu-new))
	(option-menu-mov-mol (gtk-option-menu-new)))
	
    (let ((molecule-list-ref (fill-option-menu-with-coordinates-mol-options
				  menu-ref))
	  (molecule-list-mov (fill-option-menu-with-coordinates-mol-options
				  menu-mov)))

      (gtk-option-menu-set-menu option-menu-ref-mol menu-ref)
      (gtk-option-menu-set-menu option-menu-mov-mol menu-mov)

      (gtk-container-add window ligands-vbox)
      (gtk-box-pack-start ligands-vbox title #f #f 6)
      (gtk-box-pack-start ligands-vbox ref-label #f #f 2)
      (gtk-box-pack-start ligands-vbox option-menu-ref-mol #t #f 2)
      (gtk-box-pack-start ligands-vbox ref-chain-hbox #f #f 2)
      (gtk-box-pack-start ligands-vbox mov-label #f #f 2)
      (gtk-box-pack-start ligands-vbox option-menu-mov-mol #t #f 2)
      (gtk-box-pack-start ligands-vbox mov-chain-hbox #f #f 2)
      (gtk-box-pack-start ligands-vbox h-sep #f #f 2)
      (gtk-box-pack-start ligands-vbox buttons-hbox #f #f 6)

      (gtk-box-pack-start buttons-hbox     ok-button #t #f 4)
      (gtk-box-pack-start buttons-hbox cancel-button #t #f 4)

      (gtk-box-pack-start ref-chain-hbox chain-id-ref-label #f #f 2)
      (gtk-box-pack-start ref-chain-hbox chain-id-ref-entry #f #f 2)
      (gtk-box-pack-start ref-chain-hbox resno-ref-label #f #f 2)
      (gtk-box-pack-start ref-chain-hbox resno-ref-entry #f #f 2)
      
      (gtk-box-pack-start mov-chain-hbox chain-id-mov-label #f #f 2)
      (gtk-box-pack-start mov-chain-hbox chain-id-mov-entry #f #f 2)
      (gtk-box-pack-start mov-chain-hbox resno-mov-label #f #f 2)
      (gtk-box-pack-start mov-chain-hbox resno-mov-entry #f #f 2)

      ;;(gtk-widget-set-usize chain-id-lig-entry 32 -1)
      ;;(gtk-widget-set-usize chain-id-ref-entry 32 -1)
      ;;(gtk-widget-set-usize resno-start-entry  50 -1)
      ;;(gtk-widget-set-usize resno-end-entry    50 -1)
      ;;(gtk-entry-set-text chain-id-ref-entry "A")
      ;;(gtk-entry-set-text chain-id-lig-entry "A")
      ;;(gtk-entry-set-text resno-start-entry "1")

      ;;(gtk-tooltips-set-tip (gtk-tooltips-new)
	;;		    chain-id-ref-entry 
	;;		    (string-append
	;;		    "\"A\" is a reasonable guess at the NCS master chain id.  "
	;;		    "If your ligand (specified below) is NOT bound to the protein's "
	;;		    "\"A\" chain, then you will need to change this chain and also "
	;;		    "make sure that the master molecule is specified appropriately "
	;;		    "in the Draw->NCS Ghost Control window.")
	;;		    "")
      ;;(gtk-tooltips-set-tip (gtk-tooltips-new)
	;;		    resno-end-entry 
	;;		    "Leave blank for a single residue" "")
      
      (gtk-signal-connect ok-button "clicked"
			  (lambda ()
			    (let* ((active-mol-no-ref
				    (get-option-menu-active-molecule 
				     option-menu-ref-mol
				     molecule-list-ref))
				   (active-mol-no-mov
				    (get-option-menu-active-molecule 
				     option-menu-mov-mol
				     molecule-list-mov))
				   (chain-id-ref (gtk-entry-get-text chain-id-ref-entry))
				   (chain-id-mov (gtk-entry-get-text chain-id-mov-entry))
				   (resno-ref
				    (string->number (gtk-entry-get-text resno-ref-entry)))
				   (resno-mov
				    (string->number (gtk-entry-get-text resno-mov-entry))))

			      ;; (format #t "resnos: ~s ~s~%" resno-start resno-end)
					
			      (if (number? resno-ref)
				  (if (number? resno-mov)
				      (begin 
					(overlay-my-ligands active-mol-no-mov
							    chain-id-mov
							    resno-mov
							    active-mol-no-ref
							    chain-id-ref
							    resno-ref)))))
			    (gtk-widget-destroy window)))

      (gtk-signal-connect cancel-button "clicked"
			  (lambda ()
			    (gtk-widget-destroy window)))

      (gtk-widget-show-all window))))


(define (gui-overlap-ligands  imol-ligand imol-ref chain-id-ref res-no-ref)

    ;; we don't want to overlap-ligands if there is no dictionary
    ;; for the residue to be matched to.
    ;; 
    (let* ((res-name (residue-name imol-ref chain-id-ref res-no-ref ""))
	   (restraints (monomer-restraints res-name)))
      (if (not restraints)
	  #f
	  (begin 
	    (if (not (residue-has-hetatms? imol-ref chain-id-ref res-no-ref ""))
		#f
		(begin
		  (format #t "----------- overlap-ligands ~s ~s ~s ~s ------------ ~%"
			  imol-ligand imol-ref chain-id-ref res-no-ref)
		  ;; this can return the rtop operator or the #f (for fail of course).
		  (match-ligand-torsions imol-ligand imol-ref chain-id-ref res-no-ref)
		  (overlap-ligands imol-ligand imol-ref chain-id-ref res-no-ref)
		  ))))))





(define (key-bindings-gui)

  (define (dummy-func)
    (let ((s (string-append
	      "Cannot call given function with button,\n"
	      "probably a scheme function.\n"
	      "The shortcut should still work though.")))
      (info-dialog s)
      (format #t "INFO:: ~s ~%" s)))

  (define (box-for-binding item inside-vbox buttonize-flag)

    (let ((binding-hbox (gtk-hbox-new #f 2)))
      (let* ((txt (if (string? (car (cdr item)))
		      (car (cdr item))
		      (number->string (car (cdr item)))))
	     (key-label (gtk-label-new (string-append "   " txt "   ")))
	     (name-label (gtk-label-new (car (cdr (cdr item))))))

	(if buttonize-flag
	    (let* ((button-label (string-append "   " txt "   "
						(car (cdr (cdr item)))))
		   ;; (button (gtk-button-new-with-label button-label)))
		   (button (gtk-button-new))
		   (al (gtk-alignment-new 0 0 0 0)))

	      ;; (gtk-button-set-alignment button 0 0 0 0) not in guile-gtk?
	      
	      (let ((label (gtk-label-new button-label)))
		(gtk-container-add button al)
		(gtk-container-add al label))

	      (gtk-box-pack-start binding-hbox button #t #t 0)
	      (gtk-box-pack-start inside-vbox binding-hbox #f #f 0)
	      (gtk-signal-connect 
	       button "clicked" 
	       (if (procedure? (list-ref item 3))
		   (list-ref item 3)
		   (lambda() 
		     (format #t "WARNING:: placeholder for non-scheme function: ~s~%" )))))
	
	    (begin
	      (gtk-box-pack-start binding-hbox key-label  #f #f 2)
	      (gtk-box-pack-start binding-hbox name-label #f #f 2)
	      (gtk-box-pack-start inside-vbox binding-hbox #f #f 2))))))

  ;; main line
  ;; 
  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-vbox (gtk-vbox-new #f 0))
	 (dialog-name "Key Bindings")
	 (buttons-hbox (gtk-hbox-new #f 2))
	 (close-button (gtk-button-new-with-label "  Close  "))
	 (std-frame (gtk-frame-new "Standard Key Bindings:"))
	 (usr-frame (gtk-frame-new "User-defined Key Bindings:"))
	 (std-frame-vbox (gtk-vbox-new #f 2))
	 (usr-frame-vbox (gtk-vbox-new #f 2))
	 (std-key-bindings
	  (list 
	   '("^g" "keyboard-go-to-residue")
	   '("^s" "quick-save-as")
	   '("^i" "residue info")
	   '("^z" "undo")
	   '("^y" "redo")
	   '("a" "refine with auto-zone")
	   '("b" "toggle baton swivel")
	   '("c" "toggle cross-hairs")
	   '("d" "reduce depth of field")
	   '("f" "increase depth of field")
	   '("u" "undo last navigation")
	   '("i" "toggle spin mode")
	   '("l" "label closest atom")
	   '("m" "zoom out")
	   '("n" "zoom in")
	   '("o" "other NCS chain")
	   '("p" "update position to closest atom")
	   '("s" "update skeleton")
	   '("." "up in button list")
	   '("," "down in button list")
	   )))
	   
    (gtk-signal-connect close-button "clicked"
			(lambda ()
			  (gtk-widget-destroy window)))

    (gtk-window-set-default-size window 250 350)
    (gtk-window-set-title window dialog-name)
    (gtk-container-border-width inside-vbox 4)
    
    (gtk-container-add window outside-vbox)
    (gtk-container-add outside-vbox scrolled-win)
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

    (gtk-box-pack-start inside-vbox std-frame #f #f 2)
    (gtk-box-pack-start inside-vbox usr-frame #f #f 2)
    
    (gtk-container-add std-frame std-frame-vbox)
    (gtk-container-add usr-frame usr-frame-vbox)

    (let ((scm+py-keybindings 
	   (let ((py-key-bindings
		  (if (coot-has-python?)
		      (run-python-command "key_bindings")
		      '())))
	     ;; (format #t "==== We got this for py-key-bindings: ~s ~%" py-key-bindings)
	     (if (list? py-key-bindings)
	     (append *key-bindings* py-key-bindings)
	     *key-bindings*))))
      
      (let loop ((items scm+py-keybindings))
	(cond
	 ((null? items) 'done)
	 (else 
	  (if (string? (list-ref (car items) 3)) 
	      (begin
		; BL says:: maybe should check for doublicate keys and make info
		(list-set! (car items) 2 (string-append (list-ref (car items) 2) " (py)"))
		(list-set! (car items) 3 dummy-func)))
	  (box-for-binding (car items) usr-frame-vbox #t)
	  (loop (cdr items)))))
      
      (let loop ((items (map (lambda (x) (cons 'dum x)) std-key-bindings)))
	(cond
	 ((null? items) 'done)
	 (else 
	  (box-for-binding (car items) std-frame-vbox #f)
	  (loop (cdr items))))))
    
    (gtk-box-pack-end buttons-hbox close-button #f #f 6)
    (gtk-box-pack-start outside-vbox buttons-hbox #f #f 6)

    (gtk-widget-show-all window)))



(define coot-news-info
  (let ((status 'no-news)
	(text-1 #f)
	(text-2 #f)
	(news-string-1 #f)
	(news-string-2 #f)
	(thread #f)
	(url (string-append "http:"
			    "//www.biop.ox.ac.uk/coot/software"
			    "/binaries/pre-releases/PRE-RELEASE-NOTES")))
    (lambda args

      (define (test-string) 
	(sleep 2)
	(string-append
	 "assssssssssssssssssssssssssssssssssssssss\n\n"
	 "assssssssssssssssssssssssssssssssssssssss\n\n"
	 "assssssssssssssssssssssssssssssssssssssss\n\n"
	 "\n-----\n"
	 "bill asssssssssssssssssssssssssssssssssssss\n\n"
	 "percy asssssssssssssssssssssssssssssssssssss\n\n"
	 "fred asssssssssssssssssssssssssssssssssssss\n\n"
	 "george sssssssssssssssssssssssssssssssssssss\n\n"
	 "\n-----\n"))

      (define (stop)
	'nothing)

      ;; Return (cons pre-release-news-string std-release-news-string)
      ;;
      (define (trim-news s)
	
	(let* ((sm-pre (string-match "-----" s)))
	  (if (not sm-pre)
	      (cons "nothing" "nothing")
	      (let* ((pre-news (substring s 0 (car (vector-ref sm-pre 1))))
		     (post-pre (substring s (cdr (vector-ref sm-pre 1))))
		     (sm-std (string-match "-----" post-pre)))
					   
		(if (not sm-std)
		    (cons pre-news "nothing")
		    (cons pre-news
			  (substring post-pre 0 (car (vector-ref sm-std 1)))))))))

      (define (get-news-thread)
	(let ((s (coot-get-url-as-string url)))
	  (let ((both-news (trim-news s)))
	    ;; trim-news returns (cons pre-release-news-string
	    ;;                         std-release-news-string)
	    (set! news-string-2 (car both-news))
	    (set! news-string-1 (cdr both-news))
	    (set! status 'news-is-ready))))

      ;; Handle error on coot get news thread
      ;;
      (define (coot-news-error-handler key . args)
	(format #t "error: news: error in ~s with args ~s~%" key args))

      (define (get-news)
	;; you might not be able to do this gdk-threads-enter/leave
	;; thing with windows.
	;;
	; (gdk-threads-enter)
	(set! thread (call-with-new-thread 
		      get-news-thread
		      coot-news-error-handler))
	; (format #t "captured-thread: ~s~%" thread)
	; (gdk-threads-leave)
	)

      (define (insert-string s text)
	(let ((background-colour "#c0e6c0"))
	  (gtk-text-insert text #f "black" background-colour s -1)
	  #f)) ; stop the gtk timeout func.

      (define (insert-news)
	(insert-string news-string-1 text-1)
	(insert-string news-string-2 text-2))
      
      (define (insert-no-news)
	(insert-string "  No news\n" text-1)
	(insert-string "  Yep - there really is no news\n" text-2))
	
      ; (format #t "coot-news-info called with args: ~s~%" args)
      (cond
       ((= (length args) 1)
	  (cond
	   ((eq? (car args) 'stop) (stop))
	   ((eq? (car args) 'status)
	    status)
	   ((eq? (car args) 'insert-news) (insert-news))
	   ((eq? (car args) 'insert-no-news) (insert-no-news))
	   ((eq? (car args) 'get-news) (get-news))))
       ((= (length args) 3)
	(cond
	 ((eq? (car args) 'set-text) 
	  (set! text-1 (car (cdr args)))
	  (set! text-2 (car (cdr (cdr args)))))))))))
	  
	 


(define whats-new-dialog
  (let ((text-1 #f) ; the text widgets
	(text-2 #f)
	(timer-label #f)
	(ms-step 200))
    (lambda ()

      (define check-for-new-news
	(let ((count 0))
	  (lambda ()
	    (set! count (+ 1 count))
	    (let ((timer-string (string-append 
				 (number->string 
				  (exact->inexact
				   (/ (* count ms-step) 1000))) "s")))
	      (gtk-misc-set-alignment timer-label 0.96 0.5)
	      (gtk-label-set-text timer-label timer-string))
	    ;; (format #t "check-for-new-news: count ~s~%" count)
	    (if (> count 100)
		(begin
		  (gtk-label-set-text timer-label "Timeout")
		  (coot-news-info 'insert-no-news)
		  #f) ; turn off the gtk timeout function
		(if (eq? (coot-news-info 'status) 'news-is-ready)
		    (coot-news-info 'insert-news)
		    #t)))))

      (let ((window (gtk-window-new 'toplevel))
	    (vbox (gtk-vbox-new #f 2))
	    (inside-vbox (gtk-vbox-new #f 2))
	    (scrolled-win-1 (gtk-scrolled-window-new))
	    (scrolled-win-2 (gtk-scrolled-window-new))
	    (label (gtk-label-new "Lastest Coot Release Info"))
	    (text-1 (gtk-text-new #f #f))
	    (text-2 (gtk-text-new #f #f))
	    (h-sep (gtk-hseparator-new))
	    (close-button (gtk-button-new-with-label "   Close   "))
	    (notebook (gtk-notebook-new))
	    (notebook-label-pre (gtk-label-new "Pre-release"))
	    (notebook-label-std (gtk-label-new "Release")))
	
	(set! text-1 (gtk-text-new #f #f))
	(set! text-2 (gtk-text-new #f #f))
	(set! timer-label (gtk-label-new "0.0s"))

	(gtk-window-set-default-size window 540 400)
	(gtk-window-set-policy window #t #t #f)
	(gtk-box-pack-start vbox label  #f #f 10)
	(gtk-box-pack-start vbox timer-label  #f #f 2)
	(gtk-box-pack-start vbox notebook  #t #t 4)
	(gtk-notebook-append-page notebook scrolled-win-1 notebook-label-std)
	(gtk-notebook-append-page notebook scrolled-win-2 notebook-label-pre)
	(gtk-box-pack-start vbox h-sep  #f #f  4)
	(gtk-box-pack-start vbox close-button #f #f 2)
	(gtk-container-add window vbox)

	(coot-news-info 'set-text text-1 text-2)
	(coot-news-info 'get-news)
	
	;; not -add-with-viewport because a gtk-text include native
	;; scroll support - so we can just container-add to the
	;; scrolled-window
	;; 
	;; (gtk-scrolled-window-add-with-viewport scrolled-win text)
	(gtk-container-add scrolled-win-1 text-1)
	(gtk-container-add scrolled-win-2 text-2)
	(gtk-scrolled-window-set-policy scrolled-win-1 'automatic 'always)
	(gtk-scrolled-window-set-policy scrolled-win-2 'automatic 'always)

	(gtk-signal-connect close-button "clicked"
			    (lambda ()
			      (coot-news-info 'stop)
			      (gtk-widget-destroy window)))

	(gtk-timeout-add ms-step check-for-new-news)

	(gtk-widget-show-all window)))))



(define (map-sharpening-gui imol)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 2))
	 (hbox (gtk-hbox-new #f 2))
	 (adj (gtk-adjustment-new 0.0 -30 60 0.05 2 30.1))
	 (slider (gtk-hscale-new adj))
	 (label (gtk-label-new "\nSharpen Map:"))
	 (lab2  (gtk-label-new "Add B-factor: ")))

    (gtk-box-pack-start vbox label  #f #f 2)
    (gtk-box-pack-start vbox hbox   #f #f 2)
    (gtk-box-pack-start hbox lab2   #f #f 2)
    (gtk-box-pack-start hbox slider #t #t 2)
    (gtk-container-add window vbox)
    (gtk-window-set-default-size window 500 100)
    ;; (gtk-scale-add-mark slider -30 -30 "")  guile-gtk not up to it :-(

    (gtk-signal-connect adj "value_changed"
			(lambda ()
			  (sharpen imol (gtk-adjustment-value adj))))

    (gtk-widget-show-all window)))



;; Associate the contents of a PIR file with a molecule.  Select file from a GUI.
;; 
(define (associate-pir-with-molecule-gui do-alignment?)

  (format #t "in associate-pir-with-molecule-gui~%") 
  (generic-chooser-entry-and-file-selector 
   "Associate Sequence to Model: "
   valid-model-molecule?
   "Chain ID"
   ""
   "Select PIR Alignment file"
   (lambda (imol chain-id file-name)
     (associate-pir-file imol chain-id file-name)
     (if do-alignment?
        (alignment-mismatches-gui imol)))
    *generic-chooser-entry-and-file-selector-file-entry-default-text*))


;; Make a box-of-buttons GUI for the various modifications that need
;; to be made to match the model sequence to the assigned sequence(s).
;; 
;; Call this when the associated sequence(s) have been read in already.
;; 
(define (alignment-mismatches-gui imol)

  ;; Return CA if there is such an atom in the residue, else return
  ;; the first atom in the residue.
  ;; 
  (define (get-sensible-atom-name res-info)
    (let* ((chain-id (list-ref res-info 2))
	   (res-no   (list-ref res-info 3))
	   (ins-code (list-ref res-info 4)))
      (let ((residue-atoms (residue-info imol chain-id res-no ins-code)))
	(if (null? residue-atoms)
	    " CA " ;; won't work of course
	    (let loop ((atoms residue-atoms))
	      (cond 
	       ((null? atoms) (car (car (car residue-atoms))))
	       ((string=? (car (car (car atoms))) " CA ")
		       " CA ")
	       (else (loop (cdr atoms)))))))))


  (format #t "------------- alignment-mismatches-gui called with imol ~s~%" imol)

  ;; main line
  (let ((am (alignment-mismatches imol)))

    (cond
     ((not am) (info-dialog "Sequence not associated - no alignment"))
     ((null? am) (info-dialog "No sequence mismatches"))
     (else 

      ;; (format #t "mutations: ~s~%" (list-ref am 0))
      ;; (format #t "deletions: ~s~%" (list-ref am 1))
      ;; (format #t "insertions: ~s~%"(list-ref am 2))

      ;; a button is a (list button-label button-action)
      (let ((mutate-buttons 
	     (map (lambda (res-info)
		    (let* ((chain-id (list-ref res-info 2))
			   (res-no   (list-ref res-info 3))
			   (ins-code (list-ref res-info 4))
			   (button-1-label 
			    (string-append "Mutate "
					   (list-ref res-info 2)
					   " "
					   (number->string (list-ref res-info 3))
					   " " 
					   (residue-name imol chain-id res-no ins-code)
					   " to " 
					   (car res-info)))
			   (button-1-action
			    (lambda ()
			      (set-go-to-atom-molecule imol)
			      (set-go-to-atom-chain-residue-atom-name chain-id res-no 
								      (get-sensible-atom-name res-info)))))
		      (list button-1-label button-1-action)))
		  (list-ref am 0)))

	    (delete-buttons
	     (map (lambda (res-info)
		    (let* ((chain-id (list-ref res-info 2))
			   (res-no   (list-ref res-info 3))
			   (ins-code (list-ref res-info 4))
			   (button-1-label 
			    (string-append "Delete "
					   chain-id
					   " "
					   (number->string res-no)))
			   (button-1-action
			    (lambda ()
			      (let ((atom-name (get-sensible-atom-name res-info)))
				(set-go-to-atom-molecule imol)
				(set-go-to-atom-chain-residue-atom-name chain-id res-no atom-name)))))
		      (list button-1-label button-1-action)))
		  (list-ref am 1)))

	    (insert-buttons 
	     (map (lambda (res-info)
		    (let* ((chain-id (list-ref res-info 2))
			   (res-no   (list-ref res-info 3))
			   (ins-code (list-ref res-info 4))
			   (button-1-label 
			    (string-append "Insert "
					   (list-ref res-info 2)
					   " "
					   (number->string (list-ref res-info 3))))
			   (button-1-action
			    (lambda () 
			      (info-dialog button-1-label)
			      (let ((r (nearest-residue-by-sequence imol chain-id res-no ins-code)))
				(if r 
				    (begin
				      (set-go-to-atom-molecule imol)
				      (set-go-to-atom-chain-residue-atom-name 
				       chain-id (list-ref r 2) " CA ")))))))
		      (list button-1-label button-1-action)))
		  (list-ref am 2))))

	(let ((buttons (append delete-buttons mutate-buttons insert-buttons))
	      (alignments-as-text-list (list-ref am 3)))

	  (for-each (lambda (alignment-text)
		      (info-dialog-with-markup alignment-text))
		    alignments-as-text-list)

	  (dialog-box-of-buttons "Residue mismatches"
				 (cons 300 300)
				 buttons "  Close  ")))))))


;; Wrapper in that we test if there have been sequence(s) assigned to
;; imol before we look for the sequence mismatches
;; 
(define (wrapper-alignment-mismatches-gui imol)
  
  (let ((seq-info (sequence-info imol)))
    (format #t "DEBUG:: sequence-info: ~s~%" seq-info)
    (if seq-info
	(if (null? seq-info)
	    (associate-pir-with-molecule-gui #t)
	    (alignment-mismatches-gui imol))
	(associate-pir-with-molecule-gui #t))))


;; Use clustalw to do the alignment. Then mutate using that alignement.
;; 
(define (run-clustalw-alignment imol chain-id target-sequence-pir-file)

  (define (get-clustalw2-command)
    (let ((clustalw2-command "clustalw2"))
      (if (command-in-path-or-absolute? clustalw2-command)
          clustalw2-command
          (let ((s (getenv "CCP4")))
            (if (not (string? s))
                #f
                (let ((file-path (append-dir-file (append-dir-dir s "libexec") "clustalw2")))
                  (print-var file-path)
                  (if (command-in-path-or-absolute? file-path)
                      file-path
                      #f)))))))
  

  ;; write out the current sequence
  (let ((current-sequence-pir-file "current-sequence.pir")
        (aligned-sequence-pir-file "aligned-sequence.pir")
        (clustalw2-output-file-name "clustalw2-output-file.log")
        (clustalw2-command (get-clustalw2-command)))

    (if (not clustalw2-command)
        (begin
          (format #t "No clustalw2 command~%~!")
          #f)
        (begin

          ;; if these files are not deleted/renamed then the input to clustalw2
          ;; goes wonky.

          (if (file-exists? aligned-sequence-pir-file)
              (let ((new-file-name (string-append aligned-sequence-pir-file ".old")))
                (rename-file aligned-sequence-pir-file new-file-name)))

          (if (file-exists? "current-sequence.dnd")
              (let ((new-file-name (string-append "current-sequence.dnd" ".old")))
                (rename-file "current-sequence.dnd" new-file-name)))

          (if (file-exists? "current-sequence.aln")
              (let ((new-file-name (string-append "current-sequence.aln" ".old")))
                (rename-file "current-sequence.aln" new-file-name)))

          (print-sequence-chain-general imol chain-id 1 1 current-sequence-pir-file)
          (goosh-command
           clustalw2-command
           '()
           (list
            "3"
            "1"
            target-sequence-pir-file
            "2"
            current-sequence-pir-file
            "9"
            "2"
            ""
            "4"
            ""
            aligned-sequence-pir-file
            ""
            "x"
            ""
            "x")
           clustalw2-output-file-name
           #t)

          (associate-pir-alignment-from-file imol chain-id aligned-sequence-pir-file)
          (apply-pir-alignment imol chain-id)
          (simple-fill-partial-residues imol)
          (resolve-clashing-sidechains-by-deletion imol)))))


(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Calculate")))

      (add-simple-coot-menu-menuitem
       menu "Use Clustalw for Alignment, then Mutate"
       (lambda ()
         ;; gui with molecule number chain-id and file-name
         (generic-chooser-entry-and-file-selector "Target PIR file: " valid-model-molecule?
                                                  "Chain-ID" "A" "PIR file for target sequence"
                                                  (lambda (imol chain-id target-pir-file-name)
                                                    (run-clustalw-alignment imol chain-id target-pir-file-name)
                                                    (associate-pir-file imol chain-id target-pir-file-name)
                                                    (alignment-mismatches-gui imol)
                                                    ""))))))



;; Multiple residue ranges gui
;; 
;; Create a new top level window that contains a residue-range-vbox
;; which contains a set of hboxes that contain (or allow the user to
;; enter) a residue range (chain-id resno-start resno-end).
;; 
;; The '+' and '-' buttons on the right allow the addition of extra
;; residue ranges (and remove them of course).  The residue ranges are
;; added to residue-range-widgets in a somewhat ugly manner.  Note
;; also that fill-residue-range-widgets-previous-data generates
;; residue range widgets and adds them to residue-range-widgets.
;; 
;; The interesting function is make-residue-range-frame which returns
;; the outside vbox (that contains the frame and + and - buttons, so
;; it is not a great name for the function) and each of the entries -
;; so that they can be decoded (gtk-entry-get-text) when the "Go"
;; button is pressed.
;;
;; Using the variable saved-residue-ranges, the GUI can restore itself
;; from previous (saved) residue ranges.
;;
;; (Notice that we are not dealing with insertion codes).
;; 
(define residue-range-gui 
  (let ((residue-range-widgets '())
	(saved-residue-ranges '()))
    (lambda (func function-text go-button-label)

      ;; 
      (define (n-residue-range-vboxes residue-range-widgets-vbox)
	;; 20090718 WARNING to self, in new version of guile-gnome
	;; gtk-container-children is likely to be
	;; gtk-container-get-children.
	(let ((ls (gtk-container-children residue-range-widgets-vbox)))
	  (length ls)))

      ;; Remove widget from residue-range-widgets (on '-' button
      ;; pressed)
      ;; 
      (define (remove-from-residue-range-widget widget)
	(let loop ((ls residue-range-widgets)
		   (filtered-list '()))
	  (cond 
	   ((null? ls) (set! residue-range-widgets filtered-list))
	   ((equal? widget (car (car ls)))
	    (loop (cdr ls) filtered-list))
	   (else 
	    (loop (cdr ls)
		  (cons (car ls) filtered-list))))))

      ;; 
      (define (make-residue-range-frame residue-range-vbox)
	(let ((frame (gtk-frame-new ""))
	      (outside-hbox (gtk-hbox-new #f 2))
	      (hbox (gtk-hbox-new #f 2))
	      (text-1 (gtk-label-new "  Chain-ID:"))
	      (text-2 (gtk-label-new "  Resno Start:"))
	      (text-3 (gtk-label-new "  Resno End:"))
	      (entry-1 (gtk-entry-new))
	      (entry-2 (gtk-entry-new))
	      (entry-3 (gtk-entry-new))
	      ( plus-button (gtk-button-new-with-label "+"))
	      (minus-button (gtk-button-new-with-label " - ")))

	  (gtk-box-pack-start hbox text-1  #f #f 0)
	  (gtk-box-pack-start hbox entry-1 #f #f 0)
	  (gtk-box-pack-start hbox text-2  #f #f 0)
	  (gtk-box-pack-start hbox entry-2 #f #f 0)
	  (gtk-box-pack-start hbox text-3  #f #f 0)
	  (gtk-box-pack-start hbox entry-3 #f #f 0)

	  (gtk-box-pack-start outside-hbox frame #f #f 2)
	  (gtk-container-add frame hbox)
	  (gtk-box-pack-start outside-hbox  plus-button #f #f 2)
	  (gtk-box-pack-start outside-hbox minus-button #f #f 2)

	  (gtk-signal-connect  plus-button "clicked"
			       (lambda ()
				 ;; we need to add a new residue-range
				 ;; outside-hbox into the residue-range-widgets-vbox
				 (let ((rr-frame (make-residue-range-frame residue-range-vbox)))
				   (gtk-box-pack-start residue-range-vbox (car rr-frame) #f #f 2)
				   (gtk-widget-show (car rr-frame))
				   (set! residue-range-widgets
					 (cons rr-frame residue-range-widgets)))))
	  
	  (gtk-signal-connect minus-button "clicked"
			      (lambda ()
				(let ((n (n-residue-range-vboxes residue-range-vbox)))
				  (if (> n 1)
				      (begin 
					(remove-from-residue-range-widget outside-hbox)
					(gtk-widget-destroy outside-hbox))))))

	  (map gtk-widget-show (list frame outside-hbox hbox text-1 text-2 text-3
				     entry-1 entry-2 entry-3 plus-button minus-button))

	  ;; return the thing that we need to pack and the entries we
	  ;; need to read.
	  (list outside-hbox entry-1 entry-2 entry-3)))

      ;; 
      (define (make-residue-ranges residue-range-widgets)
	(reverse (map make-residue-range residue-range-widgets)))

      ;; Return a list (list "A" 2 3) or #f on failure to decode entries.
      ;; 
      (define (make-residue-range residue-range-widget)
	;; (format #t "make a residue range using ~s~%"  residue-range-widget)
	(let* ((entry-1 (list-ref residue-range-widget 1))
	       (entry-2 (list-ref residue-range-widget 2))
	       (entry-3 (list-ref residue-range-widget 3))
	       (chain-id (gtk-entry-get-text entry-1))
	       (res-no-1 (string->number (gtk-entry-get-text entry-2)))
	       (res-no-2 (string->number (gtk-entry-get-text entry-3))))
	  (if (and (number? res-no-1)
		   (number? res-no-2))
	      (list chain-id res-no-1 res-no-2)
	      (begin
		(format #t "did not understand ~s and ~s as numbers - fail resiude range~%"
			(gtk-entry-get-text entry-2)
			(gtk-entry-get-text entry-3))
		#f))))

      ;; 
      (define (save-ranges! residue-range-widgets)
	(let ((residue-ranges (reverse 
	       (map (lambda (residue-range-widget)
		      (make-residue-range residue-range-widget))
		    residue-range-widgets))))
	  (set! saved-residue-ranges residue-ranges)))

      
      ;; range-info is (list chain-id res-no-1 res-no-2)
      ;; 
      (define (fill-with-previous-range range-info vbox-info)
	(format #t "fill-with-previous-range using ~s~%" range-info)
	(let* ((entry-1 (list-ref vbox-info 1))
	       (entry-2 (list-ref vbox-info 2))
	       (entry-3 (list-ref vbox-info 3))
	       (chain-id (list-ref range-info 0))
	       (resno-1  (list-ref range-info 1))
	       (resno-2  (list-ref range-info 2)))

	  (gtk-entry-set-text entry-1 chain-id)
	  (gtk-entry-set-text entry-2 (number->string resno-1))
	  (gtk-entry-set-text entry-3 (number->string resno-2))))

      ;; 
      (define (fill-residue-range-widgets-previous-data previous-ranges
							first-vbox-info
							residue-range-vbox)
	(format #t "first one ~s~%" (car previous-ranges))
	(if (not (null? previous-ranges))
	    (if (not (eq? #f (car previous-ranges)))
		(fill-with-previous-range (car previous-ranges) first-vbox-info)))

	(if (not (null? (cdr previous-ranges)))
	    (map (lambda (range)
		   (if (not (eq? #f range))
		       (let ((vbox-info (make-residue-range-frame residue-range-vbox)))
			 (gtk-box-pack-start residue-range-vbox (car vbox-info) #f #f 2)
			 (format #t "next one ~s~%" range)
			 (set! residue-range-widgets
			       (cons vbox-info residue-range-widgets))
			 (fill-with-previous-range range vbox-info))))
		 (cdr previous-ranges))))
      
      
      ;; main line
      ;; 
      (let* ((window (gtk-window-new 'toplevel))
	     (vbox (gtk-vbox-new #f 0))
	     (residue-range-vbox (gtk-vbox-new #t 2))
	     (residue-range-widget-info (make-residue-range-frame residue-range-vbox))
	     (hbox-buttons (gtk-hbox-new #f 0))
	     (function-label (gtk-label-new function-text))
	     (cancel-button (gtk-button-new-with-label "  Cancel  "))
	     (go-button (gtk-button-new-with-label go-button-label))
	     (h-sep (gtk-hseparator-new))
	     ;; the first residue range
	     (outside-vbox-residue-range (car residue-range-widget-info)))
	
	(set! residue-range-widgets (list residue-range-widget-info))

	;; buttons
	(gtk-box-pack-end hbox-buttons cancel-button #f #f 6)
	(gtk-box-pack-end hbox-buttons     go-button #f #f 6)
	
	;; the vbox of residue ranges
	(gtk-box-pack-start residue-range-vbox outside-vbox-residue-range #f #f 0)

	(if (not (null? saved-residue-ranges))
	    (fill-residue-range-widgets-previous-data saved-residue-ranges
						      residue-range-widget-info
						      residue-range-vbox))

	;; main vbox
	(gtk-box-pack-start vbox function-label #f #f 0)
	(let ((mc-opt-menu+model-list (generic-molecule-chooser vbox "Molecule for Ranges:")))
	  (gtk-box-pack-start vbox residue-range-vbox #f #f 2)
	  (gtk-box-pack-start vbox h-sep #t #t 6)
	  (gtk-box-pack-start vbox hbox-buttons #f #f 0)
	  
	  (gtk-container-add window vbox)
	  (gtk-container-border-width vbox 6)
	  
	  (gtk-signal-connect cancel-button "clicked"
			      (lambda ()
				(save-ranges! residue-range-widgets)
				(gtk-widget-destroy window)))
	  (gtk-signal-connect go-button "clicked"
			      (lambda()
				(save-ranges! residue-range-widgets)
				(let ((residue-ranges (make-residue-ranges residue-range-widgets))
				      (imol (apply get-option-menu-active-molecule 
						   mc-opt-menu+model-list)))
				     (if (number? imol)
					 (func imol residue-ranges))
				  (gtk-widget-destroy window))))

	  (gtk-widget-show-all window))))))



(define *additional-solvent-ligands* '())

(define *solvent-ligand-list* 
  (append
   *additional-solvent-ligands*
   (list "EDO" "GOL" "DMS" "ACT" "MPD" "CIT" "SO4" "PO4" "TRS" "TAM" "PEG" "PG4" "PE8" "EBE" "BTB")))

(define *random-jiggle-n-trials* 50)

;; add solvent molecules 
;;
;; Change the translation jiggle-factor to 1.0, so the ligand doesn't
;; move so far and get sucked into protein density (this is just a
;; temporary hack, it would be better to mask the enviroment of the
;; ligand by the surrounding atoms of the molecule to which the ligand
;; is added - that is much harder).
;;
(define (solvent-ligands-gui)

  ;; 
  (define (add-ligand-func imol tlc)
    (format #t "Add a ~a to molecule ~s here ~%" tlc imol)
    (let ((imol-ligand (get-monomer tlc)))
      (if (valid-model-molecule? imol-ligand)
	  (begin
	    ;; delete hydrogens from the ligand if the master molecule
	    ;; does not have hydrogens.
	    (if (valid-model-molecule? imol)
		(if (not (molecule-has-hydrogens? imol))
		    (delete-residue-hydrogens imol-ligand "A" 1 "" "")))
	    (if (valid-map-molecule? (imol-refinement-map))
		(begin
		  (format #t "========  jiggling!  ======== ~%")
		  (fit-to-map-by-random-jiggle imol-ligand "A" 1 "" *random-jiggle-n-trials* 1.0)
		  (with-auto-accept
		   (refine-zone imol-ligand "A" 1 1 "")))
		(format #t "======== not jiggling - no map ======== ~%"))
	    (if (valid-model-molecule? imol)
		(begin
		  (merge-molecules (list imol-ligand) imol)
		  (set-mol-active    imol-ligand 0)
		  (set-mol-displayed imol-ligand 0)))))))
  
  ;; add a button for a 3-letter-code to the scrolled vbox that runs
  ;; add-ligand-func when clicked.
  ;; 
  (define (add-solvent-button comp-id button-label inside-vbox molecule-option-menu model-list)
    (let ((button (gtk-button-new-with-label button-label)))
      (gtk-box-pack-start inside-vbox button #f #f 1)
      (gtk-widget-show button)
      (gtk-signal-connect button "clicked"
			  (lambda ()
			    (let ((imol (get-option-menu-active-molecule
					 molecule-option-menu model-list)))
			      (add-ligand-func imol comp-id))))))

  (define (comp-id->button-label comp-id)
    (auto-load-dictionary comp-id)
    (let* ((comp-id-name (comp-id->name comp-id)))
      (if comp-id-name 
	  (string-append comp-id ": " comp-id-name)
	  comp-id)))

    
  ;; main 
  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-vbox  (gtk-vbox-new #f 2))
	 (label (gtk-label-new "\nSolvent molecules added to molecule: "))
	 (menu (gtk-menu-new))
	 (frame-for-option-menu (gtk-frame-new " Choose Molecule "))
	 (vbox-for-option-menu (gtk-vbox-new #f 6))
	 (molecule-option-menu (gtk-option-menu-new))
	 (model-list (fill-option-menu-with-coordinates-mol-options menu))
	 (add-new-button (gtk-button-new-with-label "  Add a new Residue Type..."))
	 (h-sep (gtk-hseparator-new))
	 (close-button (gtk-button-new-with-label "  Close  ")))
    
    (gtk-window-set-default-size window 450 500)
    (gtk-window-set-title window "Solvent Ligands")
    (gtk-container-border-width window 8)
    (gtk-container-add window outside-vbox)
    (gtk-box-pack-start outside-vbox label #f #f 2)
    (gtk-container-add frame-for-option-menu vbox-for-option-menu)
    (gtk-box-pack-start vbox-for-option-menu molecule-option-menu #f #f 8)
    (gtk-container-border-width frame-for-option-menu 6)
    (gtk-box-pack-start outside-vbox frame-for-option-menu #f #f 2)
    (gtk-box-pack-start outside-vbox scrolled-win #t #t 0)
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    (gtk-box-pack-start outside-vbox add-new-button #f #f 6)
    (gtk-box-pack-start outside-vbox h-sep #f #f 2)
    (gtk-box-pack-start outside-vbox close-button #f #f 2)
    (gtk-option-menu-set-menu molecule-option-menu menu)
    
    (map (lambda (comp-id)
	   (let ((button-label (comp-id->button-label comp-id)))
	     (add-solvent-button comp-id button-label inside-vbox molecule-option-menu model-list)))
	 (append *solvent-ligand-list* *additional-solvent-ligands*))

    (gtk-signal-connect add-new-button "clicked"
       (lambda ()
	 (generic-single-entry "Add new 3-letter-code/comp-id"
			       ""  "  Add  "
			       (lambda (txt)
				 (set! *additional-solvent-ligands*
				       (cons txt *additional-solvent-ligands*))
				 (add-solvent-button txt
						     (comp-id->button-label txt)
						     inside-vbox 
						     molecule-option-menu 
						     model-list)))))
	 
    (gtk-signal-connect close-button "clicked"
			(lambda () 
			  (gtk-widget-destroy window)))
		      
    (gtk-widget-show-all window)))



;; USER MODS gui
;; 
(define (user-mods-gui imol pdb-file-name)
  
  ;; no alt conf, no inscode
  (define (atom-spec->string atom-spec)
    (let ((chain-id  (list-ref atom-spec 1))
	  (res-no    (list-ref atom-spec 2))
	  (ins-code  (list-ref atom-spec 3))
	  (atom-name (list-ref atom-spec 4))
	  (alt-conf  (list-ref atom-spec 5)))
      (format #f "~a ~a ~a" chain-id res-no atom-name)))

  ;; 
  (define (make-flip-buttons flip-list)
    (map (lambda (flip)
	   (let 
	       ((atom-spec    (list-ref flip 0))
		(residue-type (list-ref flip 1))
		(info-string  (list-ref flip 2))
		(set-string   (list-ref flip 3))
		(score        (list-ref flip 4)))
	     (let ((chain-id  (list-ref atom-spec 1))
		   (res-no    (list-ref atom-spec 2))
		   (ins-code  (list-ref atom-spec 3))
		   (atom-name (list-ref atom-spec 4))
		   (alt-conf  (list-ref atom-spec 5)))
	       (let ((label (string-append 
			     set-string
			     " "
			     chain-id " " (number->string res-no)
			     atom-name " : "
			     info-string " " 
			     " score " 
			     (format #f "~2,2f" score)))
		     (func (lambda ()
			     (set-go-to-atom-molecule imol)
			     (set-go-to-atom-chain-residue-atom-name 
			      chain-id res-no atom-name))))
	       (list label func)))))
	   flip-list))

  ;; 
  (define (make-no-adj-buttons no-flip-list)
    (map (lambda (no-flip-item)
	   (let ((specs (car no-flip-item))
		 (info-string (car (cdr no-flip-item))))
	     (let ((label (string-append
			   "No Adjustment "
			   (string-append-with-spaces 
			    (map atom-spec->string specs))
			   " " 
			   info-string))
		   (func (lambda ()
			   (let ((atom-spec (car specs)))
		 	     (let ((chain-id  (list-ref atom-spec 1))
				   (res-no    (list-ref atom-spec 2))
				   (ins-code  (list-ref atom-spec 3))
				   (atom-name (list-ref atom-spec 4))
				   (alt-conf  (list-ref atom-spec 5)))
			       (set-go-to-atom-molecule imol)
			       (set-go-to-atom-chain-residue-atom-name 
				chain-id res-no atom-name))))))
	       (list label func))))
	 no-flip-list))

  ;; Return a list of buttons that are (in this case, (only) clashes,
  ;; unknown/problems and flips.
  ;; 
  (define (filter-buttons flip-list)
    (let loop ((ls flip-list))
      (cond
       ((null? ls) '())
       ((let* ((flip (car ls))
	       (atom-spec    (list-ref flip 0))
	       (residue-type (list-ref flip 1))
	       (info-string  (list-ref flip 2))
	       (set-string   (list-ref flip 3))
	       (score        (list-ref flip 4)))
	  (if (< (string-length info-string) 3)
	      #f
	      ;; keep-letter: K is for keep, 
	      ;; C is clash, X is "not sure"
	      ;; F is flip.
	      (let ((keep-letter (substring info-string 2 3)))
		(or (string=? keep-letter "X")
		    (string=? keep-letter "F")
		    (string=? keep-letter "C")))))
	(cons (car ls) (loop (cdr ls))))
       (else 
	(loop (cdr ls))))))

  ;; filter flag is #t or #f
  (define (clear-and-add-back vbox flip-list no-adj-list filter-flag)    
    ;; clear
    (let ((children (gtk-container-children vbox)))
      (map gtk-widget-destroy children))
    ;; add back
    (let ((buttons (if (not filter-flag)
		       (append (make-no-adj-buttons no-adj-list)
			       (make-flip-buttons flip-list))
		       ;; filter
		       (make-flip-buttons (filter-buttons flip-list)))))
      (map (lambda (button-info)
	     (add-button-info-to-box-of-buttons-vbox button-info vbox))
	   buttons)))


  ;; main line
  ;; 
  (if (using-gui?) 
      ;; user mods will return a pair of lists.
      (let* ((flips (user-mods pdb-file-name))
	     (flip-buttons (make-flip-buttons (car flips)))
	     (no-adj-buttons (make-no-adj-buttons (car (cdr flips))))
	     (all-buttons (append no-adj-buttons flip-buttons)))
	(dialog-box-of-buttons-with-check-button 
	 "  Flips " (cons 300 300) all-buttons "  Close  "
	 "Clashes, Problems and Flips only"
	 (lambda (check-button vbox)
	   (if (gtk-toggle-button-get-active check-button)
	       (clear-and-add-back vbox (car flips) (cadr flips) #t)
	       (clear-and-add-back vbox (car flips) (cadr flips) #f)))
	 #f))))




(define average-map-gui

  ;; mav-widgets (mav is "map average") is something like a
  ;; (local-to-the-function) "static" - it can be refered to in the
  ;; call-backs.  Presumably, if you refactor this as a class, it can
  ;; be a class data item.
  ;; 
  (let ((mav-widgets '()))
    (lambda ()

      ;; On pressing the - button, delete the widget from the
      ;; mav-widgets store.  The widget is actually gtk-widget-deleted
      ;; elsewhere.
      ;; 
      (define (remove-from-mav-widgets widget)
	(let loop ((ls mav-widgets)
		   (filtered-list '()))
	  (cond 
	   ((null? ls) (set! mav-widgets filtered-list))
	   ((equal? widget (car (car ls)))
	    (loop (cdr ls) filtered-list))
	   (else 
	    (loop (cdr ls)
		  (cons (car ls) filtered-list))))))

      (define (n-mav-vboxes maps-vbox)
	;; 20090718 WARNING to self, in new version of guile-gnome
	;; gtk-container-children is likely to be
	;; gtk-container-get-children.
	(let ((ls (gtk-container-children maps-vbox)))
	  (length ls)))
	   
      ;; Return a list of the hbox the optionmenu the model-mol-list
      ;; and the entry
      ;; 
      (define (add-average-molecule-widget maps-vbox)
	(let* ((frame (gtk-frame-new #f))
	       (hbox (gtk-hbox-new #f 2))
	       (label (gtk-label-new "Weight:  "))
	       (entry (gtk-entry-new))
	       ;; fill-option-menu-with-map-mol-options
	       (menu (gtk-menu-new))
	       (map-mol-list (fill-option-menu-with-map-mol-options menu))
	       (optionmenu (gtk-option-menu-new))
	       ( plus-button (gtk-button-new-with-label "+"))
	       (minus-button (gtk-button-new-with-label " - ")))
	  (gtk-box-pack-start hbox optionmenu #f #f 2)
	  (gtk-box-pack-start hbox label #f #f 2)
	  (gtk-box-pack-start hbox entry #f #f 2)
	  (gtk-box-pack-start hbox  plus-button #f #f 2)
	  (gtk-box-pack-start hbox minus-button #f #f 2)
	  (gtk-option-menu-set-menu optionmenu menu)
          (gtk-widget-set-usize entry 40 -1)
	  (gtk-entry-set-text entry "1.0")

	  (gtk-signal-connect 
	   plus-button "clicked"
	   (lambda ()
	     (let ((mav-bits (add-average-molecule-widget maps-vbox)))
	       (gtk-widget-show (car mav-bits))
	       (set! mav-widgets (cons mav-bits mav-widgets)))))

	  ;; when the - button is clicked, delete the hbox and
	  ;; everything in it.  Except if it was the only hbox/line,
	  ;; and in that case, don't do anything.
	  (gtk-signal-connect
	   minus-button "clicked"
	   (lambda ()
	     (let ((n (n-mav-vboxes maps-vbox)))
	       (if (> n 1)
		   (begin (remove-from-mav-widgets hbox)
			  (gtk-widget-destroy hbox))))))

	  (gtk-box-pack-start maps-vbox hbox #f #f 2)
	  ;; show everything we just created.
	  (map gtk-widget-show 
	       (list frame hbox label entry menu optionmenu 
		     plus-button minus-button))
	  
	  ;; (format #t "saving map-mol-list: ~s~%" map-mol-list)
	  (list hbox optionmenu map-mol-list entry)))

	  
      ;; main line,
      ;; 
      ;; create the usual outside vbox for everything, and the
      ;; inner-vbox which is for the map hboxes.
      ;; 
      (let* ((window (gtk-window-new 'toplevel))
	     (outer-vbox (gtk-vbox-new #f 0))
	     (inner-vbox (gtk-vbox-new #f 0))
	     (title (gtk-label-new "Average Maps"))
	     (h-sep (gtk-hseparator-new))
	     (buttons-hbox (gtk-hbox-new #f 2))
	     (mav-widget (add-average-molecule-widget inner-vbox))
	     (cancel-button (gtk-button-new-with-label "  Cancel  "))
	     (ok-button (gtk-button-new-with-label "  Average Maps  ")))
	
	(gtk-container-add window outer-vbox)
	(gtk-box-pack-start outer-vbox title #f #f 2)
	(gtk-box-pack-start outer-vbox inner-vbox #f #f 2)
	(gtk-box-pack-start outer-vbox h-sep #f #f 2)
	(gtk-box-pack-start outer-vbox buttons-hbox #f #f 6)
	(gtk-box-pack-end buttons-hbox     ok-button  #f #f 6)
	(gtk-box-pack-end buttons-hbox cancel-button  #f #f 6)

	;; reset mav-widget when we start a new dialog (otherwise, it
	;; retains the old mav-widgets, which is confusing/wrong.
	;; 
	(set! mav-widgets (list mav-widget))


	(gtk-signal-connect cancel-button "clicked"
			    (lambda ()
			      (gtk-widget-destroy window)))

	;; On clicking "Average Maps", we get a list of maps to
	;; average by looking at the widgets (the option menu and the
	;; entry) for each mav-bits (i.e. each widget list that is
	;; returned by add-average-molecule-widget).
	;; 
	(gtk-signal-connect 
	 ok-button "clicked"
	 (lambda ()

	   (let ((maps-to-average-list 
		  (reverse ;; so that the top map comes first and the gridding of the 
		           ;; resulting map is this map
		   (map 
		    (lambda (mav-bits)
		      (let* ((option-menu  (list-ref mav-bits 1))
			     (map-mol-list (list-ref mav-bits 2))
			     (entry        (list-ref mav-bits 3))
			     (nov (format #t "map-mol-list: ~s~%" 
					  map-mol-list))
			     (map-selected
			      (get-option-menu-active-molecule
			       option-menu map-mol-list))
			     (text (gtk-entry-get-text entry))
			     (weight (string->number text)))
			(list map-selected weight)))
		    mav-widgets))))
	     (format #t "maps to average: ~s~%" 
		     maps-to-average-list)
	     (average-map maps-to-average-list)
	     (gtk-widget-destroy window))))
	
	;; 
	(gtk-widget-show-all window)))))

;; simple rename residue GUI
(define (rename-residue-gui)
  (let ((active-atom (active-residue)))
    (if (not active-atom)
	(info-dialog "No Residue Here")
	(let* ((aa-imol      (list-ref active-atom 0))
	       (aa-chain-id  (list-ref active-atom 1))
	       (aa-res-no    (list-ref active-atom 2))
	       (aa-ins-code  (list-ref active-atom 3))
	       (label (string-append "Rename Residue [in molecule "
				     (number->string aa-imol)
				     "]: "
				     aa-chain-id
				     (number->string aa-res-no)
				     aa-ins-code
				     " to: ")))
	  (generic-single-entry label "ALA" "Rename Residue"
				(lambda (text)
				  (set-residue-name aa-imol aa-chain-id aa-res-no aa-ins-code text)))))))





;; molecule chooser
;; 
;; coordination number chooser
;; 
;; Max dist entry
;; 
;; results vbox
;;    containing buttons with atom-spec labels
;; 
;; h-sep
;;
;;      close-button
;; 
(define (water-coordination-gui)

  (define (atom-spec->text atom-spec)
    (apply string-append 
     (map (lambda (ele)
	    ;; (format #t "ele :~s: from spec: ~s~%" ele atom-spec)
	    (cond 
	     ((string? ele) (string-append ele " "))
	     ((number? ele) (string-append
			     (number->string ele) " "))
	     (else 
	      " "))) (cdr atom-spec)))) ;; remove model
  
  (define (get-ele imol atom-spec)
    (let ((atoms (residue-info imol 
			       (list-ref atom-spec 1)
			       (list-ref atom-spec 2)
			       (list-ref atom-spec 3)))
	  (input-atom-name (list-ref atom-spec 4)))
      (if (list? atoms)
	  (let loop ((atoms atoms))
	    (cond
	     ((null? atoms) #f)
	     ((string=? (car (car (car atoms))) input-atom-name)
	      (list-ref (list-ref (car atoms) 1) 2))
	     (else 
	      (loop (cdr atoms))))))))

  ;; add info about bumps (to non-H-bonding atoms or so).  given a
  ;; water info (a central atom spec and a list of its contacts).
  (define (make-bump-text imol water-info)
    (let ((central-atom (car water-info))
	  (contact-list (cadr water-info)))
      (let ((rv ""))
	(for-each (lambda (contact)
		    (let ((ele (get-ele imol contact)))
		      (if (string? ele)
			  (if (string=? ele " C")
			      (set! rv " Bump")))))
		  contact-list)
	rv)))

  ;; 
  (define (is-a-metal-site-too? atom-spec metal-results)
    (let ((metal-specs (map car metal-results)))
      (format #t "   checking if ~s is in ~s... result: ~s ~%" atom-spec metal-specs
	      (member? atom-spec metal-specs))
      (member? atom-spec metal-specs)))

    (let* ((window (gtk-window-new 'toplevel))
	   (vbox (gtk-vbox-new #f 0))
	   (results-vbox (gtk-vbox-new #f 0))
	   (water-results-label (gtk-label-new "Other Coordinated Waters"))
	   (metal-results-vbox (gtk-vbox-new #f 0))
	   (metal-results-label (gtk-label-new "Potential Metals: "))
	   (h-sep (gtk-hseparator-new))
	   (hint-text "Molecule: ")
	   (hbox-chooser (gtk-hbox-new #f 0))
	   (hbox-max-dist (gtk-hbox-new #f 0))
	   (hbox-number-chooser (gtk-hbox-new #f 0))
	   (number-text (gtk-label-new "Coordination Number: "))
	   (molecule-chooser-option-menu-and-model-list
	    (generic-molecule-chooser hbox-chooser hint-text))
	   (molecule-chooser-option-menu 
	    (car molecule-chooser-option-menu-and-model-list))
	   (model-list (cadr molecule-chooser-option-menu-and-model-list))
	   (scrolled-win (gtk-scrolled-window-new))
	   (metal-results-scrolled-win (gtk-scrolled-window-new))
	   (number-option-menu (gtk-option-menu-new))
	   (number-list (range 3 10))
	   (number-menu (gtk-menu-new))
	   (dist-label (gtk-label-new "Max Dist: "))
	   (dist-entry (gtk-entry-new))
	   (close-button (gtk-button-new-with-label "  Close  "))
	   (apply-button (gtk-button-new-with-label "  Apply  "))
	   (hbox-buttons (gtk-hbox-new #f 6))
	   (get-molecule
	    (lambda ()
	      (get-option-menu-active-molecule
	       molecule-chooser-option-menu
	       model-list)))
	   (get-number
	    (lambda ()
	      (get-option-menu-active-item
	       number-option-menu number-list)))
	   (get-distance (lambda ()
			   (let ((t (gtk-entry-get-text dist-entry)))
			     (string->number t))))
	   (clear-previous-results 
	    (lambda ()
	      (map (lambda (this-vbox)
		     (let ((children (gtk-container-children this-vbox)))
		       (map gtk-widget-destroy children)))
		   (list results-vbox metal-results-vbox))))
	   (update-water-results
	    (lambda (imol n d)
	      (let ((results (highly-coordinated-waters imol n d)))
		(clear-previous-results)
		(format #t "DEBUG:: got  ~s water coordinations ~%" (length results))
		(if results
		    (let ((coordination-results (cdr results))
			  (metal-results (car results)))
		      (for-each
		       (lambda (water-info)
			 ;; a water spec is a central-atom spec and a list
			 ;; of its neighbours.
			 (let* ((atom-spec (car water-info))
				(t (atom-spec->text atom-spec))
				(bump-text (make-bump-text imol water-info))
				(button (gtk-button-new-with-label
					 (string-append t bump-text))))
			   (if (not (is-a-metal-site-too? atom-spec metal-results))
			       (begin
				 (gtk-box-pack-start results-vbox button #f #f 1)
				 (gtk-widget-show button)
				 (gtk-signal-connect 
				  button "clicked"
				  (lambda ()
				    (let* ((water-spec (car water-info))
					   (chain-id  (list-ref water-spec 1))
					   (res-no    (list-ref water-spec 2))
					   (atom-name (list-ref water-spec 4)))
;			  (format #t "go to ~s~%" (car water-spec))
;			  (format #t "go to :~s: :~s: :~s:~%" 
;				  chain-id res-no atom-name)
				      (set-go-to-atom-molecule imol)
				      (set-go-to-atom-chain-residue-atom-name 
				       chain-id res-no atom-name))))))))
		       coordination-results)

		      ;; now handle metal results
		      (for-each 
		       (lambda (metal-site)
			 (let* ((metal-text (cadr metal-site))
				(t (atom-spec->text (car metal-site)))
;				(nov (format #t "metal-site: ~s    t: ~s   metal-text: ~s~%" 
;					     metal-site t metal-text))
				(button-text (string-append t " Potential "
							    metal-text))
				(button (gtk-button-new-with-label button-text)))
			   (gtk-box-pack-start metal-results-vbox button #f #f 1)
			   (gtk-widget-show button)
			   (gtk-signal-connect 
			    button "clicked"
			    (lambda ()
			      (let* ((metal-spec (car metal-site))
				     (chain-id  (list-ref metal-spec 1))
				     (res-no    (list-ref metal-spec 2))
				     (atom-name (list-ref metal-spec 4)))
				(set-go-to-atom-molecule imol)
				(set-go-to-atom-chain-residue-atom-name 
				 chain-id res-no atom-name))))))
		       metal-results)))))))
				    

      (gtk-container-add window vbox)

      (fill-option-menu-with-number-options number-menu number-list 5)
		
      (gtk-option-menu-set-menu number-option-menu number-menu)
      
      (gtk-scrolled-window-add-with-viewport scrolled-win results-vbox)
      (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

      (gtk-scrolled-window-add-with-viewport metal-results-scrolled-win metal-results-vbox)
      (gtk-scrolled-window-set-policy metal-results-scrolled-win 'automatic 'always)

      (gtk-box-pack-start hbox-max-dist dist-label #f #f 2)
      (gtk-box-pack-start hbox-max-dist dist-entry #f #f 2)
      
      (gtk-box-pack-start vbox hbox-chooser #f #f 6)

      (gtk-box-pack-start hbox-number-chooser number-text #f #f 2)
      (gtk-box-pack-start hbox-number-chooser number-option-menu #f #f 2)
			  
      (gtk-box-pack-start vbox hbox-number-chooser #f #f 6)

      (gtk-box-pack-start vbox hbox-max-dist #f #f 2)

      ;; metal sites
      (gtk-box-pack-start vbox metal-results-label #f #f 2)
      (gtk-box-pack-start vbox metal-results-scrolled-win #t #t 0)

      ;; interesting water sites
      (gtk-box-pack-start vbox water-results-label #f #f 2)
      (gtk-box-pack-start vbox scrolled-win #t #t 0) ; expand fill padding
      (gtk-box-pack-start vbox h-sep #f #f 2)
      (gtk-box-pack-end hbox-buttons close-button #f #f 2)
      (gtk-box-pack-end hbox-buttons apply-button #f #f 2)
      (gtk-box-pack-start vbox hbox-buttons #f #f 2)

      (gtk-entry-set-text dist-entry "2.7")

      (gtk-signal-connect close-button "clicked" 
			  (lambda ()
			    (gtk-widget-destroy window)))

      (gtk-signal-connect dist-entry "key-press-event"
			  (lambda (event)
			    (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
				(let ((n (get-number))
				      (imol (get-molecule))
				      (d (get-distance)))
				  (if (number? d)
				      (update-water-results imol n d))))
			    #f))
				
      ;; From the Nayal and Di Cera (1996) paper, it seems that 2.7
      ;; and at least 4 oxygens is a good test for Na+ or other
      ;; interesting things.
      ;; 

      
      (gtk-signal-connect 
       apply-button "clicked"
       (lambda ()
	 (let ((n (get-number))
	       (imol (get-molecule))
	       (d (get-distance)))
	   (if (number? d)
	       (update-water-results imol n d)))))

      (gtk-widget-show-all window)))


;; interface to the difference-map based pepflip finder
(define (pepflips-by-difference-map-gui)

  (define (fill-option-menu-with-difference-map-options menu)
    (fill-option-menu-with-mol-options menu is-difference-map?))

  (let* ((window (gtk-window-new 'toplevel))
         (chooser-label "Difference map")
	 (label (gtk-label-new chooser-label))
	 (vbox (gtk-vbox-new #f 6))
	 (hbox-buttons (gtk-hbox-new #f 5))
         (hbox-sigma (gtk-hbox-new #f 5))
	 (menu-map (gtk-menu-new))
	 (menu-coords (gtk-menu-new))
	 (option-menu-map    (gtk-option-menu-new))
	 (option-menu-coords (gtk-option-menu-new))
         (n-sigma-label (gtk-label-new "N-sigma cut-off"))
         (n-sigma-entry (gtk-entry-new))
	 (ok-button (gtk-button-new-with-label "  OK  "))
	 (cancel-button (gtk-button-new-with-label " Cancel "))
	 (h-sep (gtk-hseparator-new))
	 (model-mol-list (fill-option-menu-with-coordinates-mol-options menu-coords))
	 (map-mol-list (fill-option-menu-with-difference-map-options menu-map)))
    (gtk-window-set-default-size window 370 100)
    (gtk-container-add window vbox)
    (gtk-box-pack-start vbox label #f #f 5)
    (gtk-box-pack-start vbox option-menu-map    #t #t 6)
    (gtk-box-pack-start vbox option-menu-coords #t #t 6)
    (gtk-box-pack-start hbox-sigma n-sigma-label #f #f 6)
    (gtk-box-pack-start hbox-sigma n-sigma-entry #f #f 6)
    (gtk-box-pack-start vbox hbox-sigma  #t #f 6)
    (gtk-box-pack-start vbox h-sep #t #f 2)
    (gtk-box-pack-start vbox hbox-buttons #f #f 5)
    (gtk-box-pack-start hbox-buttons ok-button #t #f 5)
    (gtk-box-pack-start hbox-buttons cancel-button #t #f 5)
    
    (gtk-option-menu-set-menu option-menu-map    menu-map)
    (gtk-option-menu-set-menu option-menu-coords menu-coords)
    (gtk-entry-set-text n-sigma-entry "4.0")

    (gtk-signal-connect cancel-button "clicked" (lambda () (gtk-widget-destroy window)))
    (gtk-signal-connect ok-button "clicked"
                        (lambda ()
                          (let* ((imol-coords (get-option-menu-active-molecule
                                               option-menu-coords model-mol-list))
                                 (imol-map (get-option-menu-active-molecule
                                            option-menu-map map-mol-list))
                                 (n-sigma-str (gtk-entry-get-text n-sigma-entry))
                                 (n-sigma (string->number n-sigma-str)))
                            (let ((specs (pepflip-using-difference-map-scm imol-coords imol-map n-sigma)))
                              (format #t "dialog for specs\n")
                              (format #t "specs: ~s~%" specs)
                              (if (null? specs)
                                  (info-dialog "No pepflips found")
                                  (dialog-box-of-buttons-from-specs "Potential Pepflip"
                                                                    (cons 240 170)
                                                                    imol-coords specs))
                              (gtk-widget-destroy window)))))

    (gtk-widget-show-all window)))



;; return a list, or #f (e.g. if not in same chain and molecule)
;; 
(define (min-max-residues-from-atom-specs specs)

  (print-var specs)

  (let ((min-res-no #f)
	(max-res-no #f)
	(chain-id #f))

    (for-each (lambda (spec)
		(let ((spec-chain (list-ref spec 2))
		      (res-no     (list-ref spec 3))
		      (spec-model (list-ref spec 1)))

		  (if (string? chain-id)
		      (if (not (string=? chain-id spec-chain))
			  (begin
			    (set! chain-id #f))

			  (begin
			    (if (not (number? min-res-no))
				(set! min-res-no res-no)
				(if (< res-no min-res-no)
				    (set! min-res-no res-no)))
			    (if (not (number? max-res-no))
				(set! max-res-no res-no)
				(if (> res-no max-res-no)
				    (set! max-res-no res-no)))))
		      (begin
			(set! chain-id spec-chain)
			(if (not (number? min-res-no))
			    (set! min-res-no res-no)
			    (if (< res-no min-res-no)
				(set! min-res-no res-no)))
			(if (not (number? max-res-no))
			    (set! max-res-no res-no)
			    (if (> res-no max-res-no)
				(set! max-res-no res-no)))))))
	      specs)
			

    (if (or (eq? min-res-no #f)
	    (eq? max-res-no #f)
	    (eq? chain-id #f))
	#f
	(list min-res-no max-res-no chain-id))))
    

;; by default, rename loop residues to UNK.  If scheme true, then
;; leave them as the residue names found in the database.
(define *db-loop-preserve-residue-names* #f)

(define (click-protein-db-loop-gui)
  (generic-number-chooser (range 2 10) 4 
			  "Number of residues for basis"
			  "Pick Atoms..."
			  (lambda (n)
			    (user-defined-click 
			     n
			     (lambda (atom-specs)
			       (let ((residue-specs (map atom-spec->residue-spec atom-specs))
				     (imol (car (cdar atom-specs)))
				     (min-max-and-chain-id (min-max-residues-from-atom-specs atom-specs)))

				 (if (not (list? min-max-and-chain-id))
				     (info-dialog "WARNING:: Picked atoms not in same molecule and chain")
				     (let ((loop-mols
					    (protein-db-loops imol residue-specs 
							      (imol-refinement-map)
							      10 *db-loop-preserve-residue-names*)))

				       (let ((imol-loop-orig (car (car loop-mols)))
					     (imol-loops-consolodated (cadr (car loop-mols)))
					     (loop-mols (cadr loop-mols))
					     (min-resno (list-ref min-max-and-chain-id 0))
					     (max-resno (list-ref min-max-and-chain-id 1))
					     (chain-id  (list-ref min-max-and-chain-id 2)))

					 (if (valid-model-molecule? imol-loop-orig)
					     (if (> (length loop-mols) 0)
						 (let ((buttons
							(map (lambda (loop-mol)
							       (list 
								(string-append 
								 (number->string loop-mol)
								 " "
								 (molecule-name loop-mol))
								(lambda ()
								  (copy-residue-range imol chain-id
										      loop-mol chain-id
										      min-resno max-resno))))
							     loop-mols)))
						   (dialog-box-of-buttons
						    "Loop Candidates" 
						    (cons 360 200)
						    (append
						     (list 
						      (list "Original Loop"
							    (lambda () 
							      (copy-residue-range imol chain-id 
										  imol-loop-orig chain-id
										  min-resno max-resno)))
						      (list "Toggle All Candidate Loops"
							    (lambda () 
							      (toggle-display-mol imol-loops-consolodated)
							      (toggle-active-mol  imol-loops-consolodated)
							      ))
						      )
						     buttons)
						    " Close ")))))))))))))

(define (refmac-multi-sharpen-gui)
  (let ((window (gtk-window-new 'toplevel))
	;; boxes
	(vbox (gtk-vbox-new #f 4))
	(hbox-1 (gtk-hbox-new #f 4))
	(hbox-2 (gtk-hbox-new #f 4))
	(hbox-3 (gtk-hbox-new #f 4))
	;; menus
	(option-menu-map (gtk-option-menu-new))
	(option-menu-b-factor (gtk-option-menu-new))
	(option-menu-n-levels (gtk-option-menu-new))
	(menu-map (gtk-menu-new))
	(menu-b-factor (gtk-menu-new))
	(menu-n-levels (gtk-menu-new))
	;; labels
	(map-label (gtk-label-new "Map "))
	(sb-label  (gtk-label-new "Sharpen & Blur in "))
	(levels-label  (gtk-label-new " levels up to "))
	(A-label  (gtk-label-new " A*A"))
	;; separate
	(h-sep (gtk-hseparator-new))
	;; buttons
	(ok-button (gtk-button-new-with-label "   OK   "))
	(cancel-button (gtk-button-new-with-label " Cancel "))
	(n-levels-list (list 1 2 3 4 5 6 8 10 12 15))
	(b-factor-list (list 50 100 200 400 800 2000)))

    (let ((map-molecule-list (fill-option-menu-with-map-mol-options
				  menu-map)))

      ;; fill option menus
      (for-each (lambda (item)
		  (let* ((mlabel-str (number->string item))
			 (menuitem (gtk-menu-item-new-with-label mlabel-str)))
		    (gtk-menu-append menu-n-levels menuitem)))
		n-levels-list)
      (for-each (lambda (item)
		  (let* ((mlabel-str (number->string item))
			 (menuitem (gtk-menu-item-new-with-label mlabel-str)))
		    (gtk-menu-append menu-b-factor menuitem)))
		b-factor-list)
      (gtk-menu-set-active menu-n-levels 3)
      (gtk-menu-set-active menu-b-factor 2)
    
      (gtk-window-set-title window "Refmac for Sharpening & Blurring")
      (gtk-option-menu-set-menu option-menu-map menu-map)
      (gtk-option-menu-set-menu option-menu-b-factor menu-b-factor)
      (gtk-option-menu-set-menu option-menu-n-levels menu-n-levels)
      (gtk-box-pack-start hbox-1 map-label #f #f 2)
      (gtk-box-pack-start hbox-1 option-menu-map #f #f 2)
      (gtk-box-pack-start hbox-2 sb-label #f #f 2)
      (gtk-box-pack-start hbox-2 option-menu-n-levels #f #f 2)
      (gtk-box-pack-start hbox-2 levels-label #f #f 2)
      (gtk-box-pack-start hbox-2 option-menu-b-factor #f #f 2)
      (gtk-box-pack-end hbox-3 cancel-button #f #f 12)
      (gtk-box-pack-end hbox-3     ok-button #f #f 12)

      (gtk-box-pack-start vbox hbox-1 #f #f 6)
      (gtk-box-pack-start vbox hbox-2 #f #f 6)
      (gtk-box-pack-start vbox h-sep  #f #f 2)
      (gtk-box-pack-start vbox hbox-3 #f #f 6)

      (gtk-signal-connect cancel-button "clicked"
			  (lambda() (gtk-widget-destroy window)))

      (gtk-signal-connect ok-button "clicked"
			  (lambda()

			    ;; if starting-map.mtz exists, move it.

			    ;; get max-b and n-levels and map file name
			    (let* ((lb b-factor-list)
				   (ln n-levels-list)
				   (lmols map-molecule-list)
				   (max-b            (get-option-menu-active-item option-menu-b-factor lb))
				   (n-levels         (get-option-menu-active-item option-menu-n-levels ln))
				   (active-item-imol (get-option-menu-active-item option-menu-map lmols))
				   ;; There is no function to get a map file name from a molecule
				   ;; It is not stored.
				   (map-file-name (molecule-name active-item-imol))
				   (map-file-name-stub (strip-path (file-name-sans-extension map-file-name)))
				   (refmac-output-mtz-file-name (string-append "starting-map-"
									       map-file-name-stub
									       ".mtz"))
				   (log-file-name (string-append
						   "refmac-multisharp-"
						   map-file-name-stub
						   ".log")))

			      (if (not (file-exists? map-file-name))
				  (begin
				    (info-dialog (string-append "WARNING:: file not found " map-file-name)))

				  (begin

				    (format #t "active-item-imol: ~s~%" active-item-imol)

				    (let* ((step-size (/ max-b n-levels))
					   (numbers-string
					    (apply string-append (map (lambda(i)
									(let ((lev (* step-size (+ i 1))))
									  (string-append
									   (number->string (exact->inexact lev)) " ")))
								      (range n-levels))))
					   (blur-string  (string-append "SFCALC BLUR  " numbers-string))
					   (sharp-string (string-append "SFCALC SHARP " numbers-string)))

				      (let ((cmd-line-args (list "MAPIN" map-file-name))
					    (data-lines (list "MODE SFCALC"
							      blur-string
							      sharp-string
							      "END"))
					    (this-dir (getcwd)))
					(if (not (directory-is-modifiable? this-dir))
					    (info-dialog "WARNING:: Current directory is not writable")
					    (let ((s (goosh-command "refmac5"
								    cmd-line-args
								    data-lines
								    log-file-name
								    #f)))
					      (if (not (ok-goosh-status? s))

						  (begin
						    (info-dialog "WARNING:: refmac5 failed"))

						  ;; Happy path
						  (begin
						    (format #t "s: ~s~%" s)
						    (if (not (file-exists? "starting_map.mtz"))
							(begin
							  (format #t "WARNING:: starting_map.mtz does not exist~%"))
							(begin
							  (format #t "INFO renaming starting_map.mtz to ~s~%"
								  refmac-output-mtz-file-name)
							  (rename-file  "starting_map.mtz" refmac-output-mtz-file-name)
							  ;; offer a read-mtz dialog
							  (manage-column-selector refmac-output-mtz-file-name)
							  ))))))))))

			      (gtk-widget-destroy window))))

    (gtk-container-add window vbox)
    (gtk-widget-show-all window))))

(define (auto-assign-sequence-from-map)
  (using-active-atom
   (let ((imol aa-imol)
         (ch-id aa-chain-id))
     (let ((imol-map (imol-refinement-map)))
       (let ((fragment-residues (linked-residues-scm aa-res-spec aa-imol 1.7)))
         (let ((residue-number-list (map (lambda (item) (list-ref item 2)) fragment-residues)))
           (let ((resno-start (apply min residue-number-list))
                 (resno-end   (apply max residue-number-list)))
             (let ((new-sequence (sequence-from-map imol ch-id resno-start resno-end imol-map)))
               (set-rotamer-search-mode (ROTAMERSEARCHLOWRES))
               (format #t "INFO:: mutate-residue-range ~s ~s ~s ~s ~s~%" imol ch-id resno-start resno-end new-sequence)
               (format #t "INFO:: mutate-residue-range lengths: ~s ~s~%"
                       (+ (- resno-end resno-start) 1) (string-length new-sequence))
               (mutate-residue-range imol ch-id resno-start resno-end new-sequence)
               (backrub-rotamers-for-chain imol ch-id)
               (refine-residues imol fragment-residues)))))))))

;;
(define (add-module-cryo-em)
  (if (defined? 'coot-main-menubar)
      (add-module-cryo-em-gui)))

;;
(define (add-module-ccp4)
  (if (defined? 'coot-main-menubar)
      (add-module-ccp4-gui)))

(define (add-module-pdbe)
  (if (defined? 'coot-main-menubar)
      (add-module-pdbe-gui)))

;;
(define (add-module-cryo-em-gui)
  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "Cryo-EM")))

        (add-simple-coot-menu-menuitem
         menu "Go To Map Molecule Middle"
         (lambda ()
          (go-to-map-molecule-centre (imol-refinement-map))))

        (add-simple-coot-menu-menuitem
         menu "Go To Box Middle" go-to-box-middle)

        (add-simple-coot-menu-menuitem
         menu "Sharpen/Blur..."
         sharpen-blur-map-gui)

        (add-simple-coot-menu-menuitem
         menu "Multi-sharpen using Refmac..."
         refmac-multi-sharpen-gui)

        (add-simple-coot-menu-menuitem
         menu "Flip Hand of Map..."
         (lambda ()
           (map-molecule-chooser-gui "Select" (lambda (imol) (flip-hand imol)))))

        (add-simple-coot-menu-menuitem
         menu "Add molecular symmetry using MTRIX"
         (lambda ()
           (using-active-atom
            (set-show-symmetry-master 1)
            (add-molecular-symmetry-from-mtrix-from-self-file aa-imol))))

        (add-simple-coot-menu-menuitem
         menu "Align and Mutate using ClustalW2"
         (lambda ()
           'x
           (generic-chooser-entry-and-file-selector
            "Align Sequence to Model: "
            valid-model-molecule?
            "Chain ID"
            ""
            "Select PIR Alignment file"
            (lambda (imol chain-id target-sequence-pir-file)
              (run-clustalw-alignment imol chain-id target-sequence-pir-file)))))

        (add-simple-coot-menu-menuitem
         menu "Assign Sequence Based on Associated Sequence"
         (lambda ()
           (assign-sequence-to-active-fragment)))

        (add-simple-coot-menu-menuitem
         menu "Auto-assign Sequence Based on Map"
         (lambda ()
           (auto-assign-sequence-from-map)))

        (add-simple-coot-menu-menuitem
         menu "No Auto-Recontour Map Mode" (lambda () (set-auto-recontour-map 0)))

        (add-simple-coot-menu-menuitem
         menu "Enable Auto-Recontour Map Mode" (lambda () (set-auto-recontour-map 1)))

        (add-simple-coot-menu-menuitem
         menu "Interactive Nudge Residues..."
         (lambda ()
           (using-active-atom (nudge-residues-gui aa-imol aa-res-spec)))))))
;;
(define (add-module-ccp4-gui)
  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "CCP4")))

        (add-simple-coot-menu-menuitem
         menu "Make Link via Acedrg"
         (lambda ()
           (acedrg-link-generation-control-window))))))

(define (add-module-pdbe-gui)
  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "PDBe")))

        ;; ---------------------------------------------------------------------
        ;;     Recent structures from the PDBe
        ;; ---------------------------------------------------------------------
        ;;
        ;; 20110921 too crashy at the moment (something to do with lots of threads?)
        ;;
        (add-simple-coot-menu-menuitem
         menu "PDBe recent structures..."
         pdbe-latest-releases-gui)

        (add-simple-coot-menu-menuitem
         menu "Get from PDBe..."
         (lambda ()
           (let ((mess
                  (if (command-in-path? "refmac5")
                      "Get PDBe accession code"
                      (string-append
                       "\n  WARNING:: refmac5 not in the path - SF calculation will fail  \n\n"
                       "Get PDBe accession code"))))
             (generic-single-entry mess
                                   "" " Get it "
                                   (lambda (text)
                                     ;; fire off something that is controlled by a time-out -
                                     ;; doesn't return a useful value.
                                     (pdbe-get-pdb-and-sfs-cif 'include-sfs (string-downcase text))))))))))


;; let the c++ part of coot know that this file was loaded:
(set-found-coot-gui)

;;; Local Variables:
;;; mode: scheme
;;; End:
