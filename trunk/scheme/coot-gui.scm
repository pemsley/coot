#!/home/paule/build/bin/guile -s
!#

;;;; Copyright 2001 Neil Jerram
;;;; Copyright 2002, 2003, 2004, 2005, 2006, 2007 by Paul Emsley, 
;;;;   The University of York
 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 2 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

;;;; coot-gui is based on Neil Jerram's guile-gui, I modified the
;;;; widget somewhat to make it more pretty.
;;;; 
;;;; Wouldn't it be nice to have command completion?

(use-modules (gtk gdk)
             (gtk gtk)
             (gui event-loop)
             (gui entry-port)
             (gui text-output-port))

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
         (icolour (gdk-color-parse "red")))
    (gtk-window-set-policy window #t #t #f)
    (gtk-window-set-default-size window 400 250)
    ; (gtk-box-pack-start vbox text #t #t 5)
    (gtk-container-add window vbox)
    (gtk-window-set-title window "Coot Scripting")
    (gtk-container-border-width vbox 5)

    (gtk-box-pack-start hbox label #f #f 0) ; expand fill padding
    (gtk-box-pack-start hbox entry #t #t 0)
    (gtk-box-pack-start vbox hbox #f #f 5)
    (gtk-container-add scrolled-win text)
    (gtk-container-add vbox scrolled-win)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    (gtk-box-pack-end vbox close-button #f #f 5)

    (gtk-entry-set-text entry "(inactive)")
    ;(gtk-widget-set-sensitive entry #f)  ; testing 

    (add-hook! (entry-read-complete-hook entry)
               (lambda (str)
                 (with-fluids (; (text-output-font ifont)
                               (text-output-colour icolour))
                              (display str))
                 (newline)))

;    (gtk-signal-connect entry "key-press-event"
;			(lambda args
;			  (format #t "key press ~s~%" args)))

    (gtk-signal-connect text "key-press-event"
			(lambda args
			  ; (format #t "text key press ~s~%" args)))
			  (format #t "Don't type here.  Type in the white entry bar.~%")))

    (gtk-signal-connect close-button "clicked" 
			(lambda args 
			  (gtk-widget-destroy window)))

;    (gtk-signal-connect window "delete-event"
;			(lambda args
;			  (if (gtk-standalone?) ; it is
;			      (format #t "is standalone~%")  
;			      (format #t "is not standalone~%"))
;			  (gtk-widget-destroy window)))
			  ; (gtk-exit))) // kill the program

    (gtk-widget-show-all window)

    (set-current-input-port port)
    (set-current-output-port oport)
    (set-current-error-port oport)
    (top-repl)))

; let the c++ part of mapview know that this file was loaded:
(set-found-coot-gui)

(set-repl-prompt! "coot> ")

; (guile-gui)


;; The callback from pressing the Go button in the smiles widget, an
;; interface to run libcheck.
;; 
(define (handle-smiles-go tlc-entry smiles-entry)
  
  (let ((tlc-text (gtk-entry-get-text tlc-entry))
	(smiles-text (gtk-entry-get-text smiles-entry)))

    (if (> (string-length smiles-text) 0)

	(let ((three-letter-code 
	       (cond 
		((and (> (string-length tlc-text) 0)
		      (< (string-length tlc-text) 4))
		 tlc-text)
		((> (string-length tlc-text) 0)
		 (substring tlc-text 0 3))
		(else "DUM"))))
	  
	  ;; OK, let's run libcheck
	  (let* ((smiles-file (string-append "coot-" three-letter-code ".smi"))
		 (libcheck-data-lines
		  (list "N"
			(string-append "MON " three-letter-code)
			(string-append "FILE_SMILE " smiles-file)
			""))
		 (log-file-name (string-append "libcheck-" three-letter-code))
		 (pdb-file-name (string-append "libcheck_" three-letter-code ".pdb"))
		 (cif-file-name (string-append "libcheck_" three-letter-code ".cif")))
	    
					; write the smiles strings to a file
	    (call-with-output-file smiles-file
	      (lambda (port)
		(format port "~a~%" smiles-text)))
	    
	    (let ((status (goosh-command libcheck-exe '() libcheck-data-lines log-file-name #t)))
					; the output of libcheck goes to libcheck.lib, we want it in
					; (i.e. overwrite the minimal description in cif-file-name
	      (if (number? status)
		  (if (= status 0)
		      (begin
			(if (file-exists? "libcheck.lib")
			    (rename-file "libcheck.lib" cif-file-name))
			(let ((sc (rotation-centre))
			      (imol (handle-read-draw-molecule-with-recentre pdb-file-name 0)))
			  (if (valid-model-molecule? imol)
			      (let ((mc (molecule-centre imol)))
				(apply translate-molecule-by (cons imol (map - sc mc))))))
			(read-cif-dictionary cif-file-name)))
		  (format #t "OOPs.. libcheck returned exit status ~s~%" status))))))))


;; smiles GUI
(define (smiles-gui)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox1 (gtk-hbox-new #f 0))
	 (hbox2 (gtk-hbox-new #f 0))
	 (tlc-label (gtk-label-new "  3-letter code "))
	 (tlc-entry (gtk-entry-new))
	 (smiles-label (gtk-label-new "SMILES string "))
	 (smiles-entry (gtk-entry-new))
	 (go-button (gtk-button-new-with-label "  Go  ")))

    (gtk-box-pack-start vbox hbox1 #f #f 0)
    (gtk-box-pack-start vbox hbox2 #f #f 0)
    (gtk-box-pack-start vbox go-button #f #f 6)
    (gtk-box-pack-start hbox1 tlc-label #f 0)
    (gtk-box-pack-start hbox1 tlc-entry #f 0)
    (gtk-box-pack-start hbox2 smiles-label #f 0)
    (gtk-box-pack-start hbox2 smiles-entry #f 0)
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
				(gtk-widget-destroy window)))))
			  

    (gtk-widget-show-all window)))


;; Generic single entry widget
;; 
;; Pass the hint labels of the entries and a function that gets called
;; when user hits "Go".  The handle-go-function accepts one argument
;; that is the entry text when the go button is pressed.
;; 
(define generic-single-entry
  (lambda (function-label entry-1-default-text go-button-label handle-go-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox1 (gtk-hbox-new #f 0))
	 (hbox2 (gtk-hbox-new #f 0))
	 (function-label (gtk-label-new function-label))
	 (smiles-entry (gtk-entry-new))
	 (go-button (gtk-button-new-with-label go-button-label)))

    (gtk-box-pack-start vbox hbox1 #f #f 0)
    (gtk-box-pack-start vbox hbox2 #f #f 0)
    (gtk-box-pack-start vbox go-button #f #f 6)
    (gtk-box-pack-start hbox1 function-label #f 0)
    (gtk-box-pack-start hbox2 smiles-entry #f 0)
    (gtk-container-add window vbox)
    (gtk-container-border-width vbox 6)

    (if (string? entry-1-default-text)
	(gtk-entry-set-text smiles-entry entry-1-default-text))

    (gtk-signal-connect go-button "clicked"
			(lambda args
			  (handle-go-function (gtk-entry-get-text smiles-entry))
			  (gtk-widget-destroy window)))
    
    (gtk-signal-connect smiles-entry "key-press-event"
			(lambda (event)
			  (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
			      (begin
				(handle-go-function (gtk-entry-get-text smiles-entry))
				(gtk-widget-destroy window)))))

    (gtk-widget-show-all window))))

;; generic double entry widget
;; 
;; pass a the hint labels of the entries and a function that gets
;; called when user hits "Go" (which takes to string aguments).
;; 
(define (generic-double-entry label-1 label-2 entry-1-default-text entry-2-default-text go-button-label handle-go-function)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox1 (gtk-hbox-new #f 0))
	 (hbox2 (gtk-hbox-new #f 0))
	 (tlc-label (gtk-label-new label-1))
	 (tlc-entry (gtk-entry-new))
	 (smiles-label (gtk-label-new label-2))
	 (smiles-entry (gtk-entry-new))
	 (go-button (gtk-button-new-with-label go-button-label)))
    
    (gtk-box-pack-start vbox hbox1 #f #f 0)
    (gtk-box-pack-start vbox hbox2 #f #f 0)
    (gtk-box-pack-start vbox go-button #f #f 6)
    (gtk-box-pack-start hbox1 tlc-label #f 0)
    (gtk-box-pack-start hbox1 tlc-entry #f 0)
    (gtk-box-pack-start hbox2 smiles-label #f 0)
    (gtk-box-pack-start hbox2 smiles-entry #f 0)
    (gtk-container-add window vbox)
    (gtk-container-border-width vbox 6)
    
    (if (string? entry-1-default-text) 
	(gtk-entry-set-text tlc-entry entry-1-default-text))
    
    (if (string? entry-2-default-text )
	(gtk-entry-set-text smiles-entry entry-2-default-text))

    (gtk-signal-connect go-button "clicked"
			(lambda args
			  (handle-go-function (gtk-entry-get-text tlc-entry)
					      (gtk-entry-get-text smiles-entry))
			  (gtk-widget-destroy window)))
    
    (gtk-signal-connect smiles-entry "key-press-event"
			(lambda (event)
			  (if (= 65293 (gdk-event-keyval event)) ; GDK_Return
			      (begin
				(handle-go-function (gtk-entry-get-text tlc-entry)
						    (gtk-entry-get-text smiles-entry))
				(gtk-widget-destroy window)))))

    (gtk-widget-show-all window)))

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
(define interesting-things-with-fix-maybe
  (lambda (title baddie-list)

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

	    ; (gtk-box-pack-start inside-vbox button #f #f 2) no buttons hbox
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

      (gtk-widget-show-all window))))

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
		      (gtk-menu-append menu menuitem))
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


(define molecule-chooser-gui-generic
  (lambda (chooser-label callback-function option-menu-fill-function)

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
				    (format #t "INFO: operating on molecule number ~s~%" 
					    active-mol-no)
				    (callback-function active-mol-no))))
			      
			    (gtk-widget-destroy window)))

      (gtk-signal-connect cancel-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window)))
      
      (gtk-widget-show-all window))))

;; Fire up a molecule chooser dialog, with a given label and on OK we
;; call the call-back-fuction with an argument of the chosen molecule
;; number. 
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

(define map-molecule-chooser-gui
  (lambda (chooser-label callback-function)
    (molecule-chooser-gui-generic chooser-label callback-function 
				  fill-option-menu-with-map-mol-options)))

;; A pair of widgets, a molecule chooser and an entry.  The
;; callback-function is a function that takes a molecule number and a
;; text string.
;; 
(define (generic-chooser-and-entry chooser-label entry-hint-text defaut-entry-text callback-function)
  
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
				  (format #t "INFO: operating on molecule number ~s~%" 
					  active-mol-no)
				  (let ((text (gtk-entry-get-text entry)))
				    (callback-function active-mol-no text)))))
			  
			  (gtk-widget-destroy window)))
    
    (gtk-signal-connect cancel-button "clicked"
			(lambda args
			  (gtk-widget-destroy window)))
    
    (gtk-widget-show-all window)))


;; If a menu with label menu-label is not found in the coot main
;; menubar, then create it and return it. 
;; If it does exist, simply return it.
;;  
(define coot-menubar-menu
  (lambda (menu-label)
    
    (define menu-bar-label-list
      (lambda () 
	(map 
	 (lambda (menu-child)
	   (let* ((ac-lab-ls (gtk-container-children menu-child))
		  ;; ac-lab-ls is a GtkAccelLabel in a list
		  (ac-lab (car ac-lab-ls))
		  ;; ac-lab-ls is a simple GtkAccelLabel
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
	    menu)))))

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
(define add-simple-coot-menu-menuitem
  (lambda (menu menu-item-label activate-function)

    (let ((submenu (gtk-menu-new))
	  (sub-menuitem (gtk-menu-item-new-with-label menu-item-label)))
      
      (gtk-menu-append menu sub-menuitem)
      (gtk-widget-show sub-menuitem)
      
      (gtk-signal-connect sub-menuitem "activate"
			  activate-function))))




;;; Make an interesting things GUI for residues of molecule number
;;; imol that have alternate conformations.
;;;
(define (alt-confs-gui imol)

  (if (valid-model-molecule? imol)

      (let* ((alt-conf-residues (residues-with-alt-confs imol))
	     (centre-atoms (map (lambda (spec)
				  (if (list? spec)
				      (apply residue-spec->atom-for-centre (cons imol (cdr spec)))
				      #f))
				alt-conf-residues)))
	(interesting-things-gui
	 "Residues with Alt Confs"
	 (map (lambda (alt-conf-residue-cpmd centre-atom)
		(let* ((alt-conf-residue (cdr alt-conf-residue-cpmd))
		       (label (string-append 
			       (car alt-conf-residue) " "
			       (number->string (car (cdr alt-conf-residue))) " "
			       (car (cdr (cdr alt-conf-residue))) " "
			       (car centre-atom) " "
			       (car (cdr centre-atom)))))
		       (append (list label imol) alt-conf-residue centre-atom)))
	      alt-conf-residues centre-atoms)))))

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

;;  A gui that makes a generic number chooser
;; the go functionis a lambda function that takes the value of 
;; the active menu item - as a number.
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
				 (get-option-menu-active-molecule 
				  option-menu number-list)))
			    (go-function active-number)
			    (gtk-widget-destroy window))))
    (gtk-signal-connect cancel-button "clicked" 
			(lambda()
			  (gtk-widget-destroy window)))

    (gtk-option-menu-set-menu option-menu menu)

    (gtk-widget-show-all window)))

;; The gui for the strand placement function
;;
(define place-strand-here-gui
  (lambda ()

    (generic-number-chooser (number-list 4 12) 7
			    " Estimated number of residues in strand"
			    "  Go  "
			    (lambda (n)
			      (place-strand-here n)))))

	 
;;; Local Variables:
;;; mode: scheme
;;; End:
