;;;; Copyright 2007 by The University of York
;;;; Copyright 2007 by Paul Emsley
;;;; Copyright 2007, 2009 by The University of Oxford
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


(define (add-plugin-shelx)

  (if (defined? 'coot-main-menubar)

      (let ((menu (coot-menubar-menu "SHELX")))

	(add-simple-coot-menu-menuitem
	 menu "SHELXL Refine..."
	 (lambda ()

	   (let ((window (gtk-window-new 'toplevel))
		 (hbox (gtk-vbox-new #f 0))
		 (vbox (gtk-hbox-new #f 0))
		 (go-button (gtk-button-new-with-label "  Refine...  "))
		 (cancel-button (gtk-button-new-with-label "  Cancel  "))
		 (entry-hint-text "HKL data filename \n(leave blank for default)")
		 (chooser-hint-text " Choose molecule for SHELX refinement  ")
		 (h-sep (gtk-hseparator-new)))

	     (gtk-container-add window hbox)
	     (let ((option-menu-mol-list-pair (generic-molecule-chooser 
					       hbox chooser-hint-text))
		   (entry (file-selector-entry hbox entry-hint-text)))
	       (gtk-signal-connect go-button "clicked"
				   (lambda () 
				     (let ((txt (gtk-entry-get-text entry))
					   (imol (apply get-option-menu-active-molecule 
							option-menu-mol-list-pair)))
				       (if (number? imol)
					   (editable-shelx-gui imol txt))
				       (gtk-widget-destroy window))))
	       (gtk-signal-connect cancel-button "clicked"
				   (lambda ()
				     (gtk-widget-destroy window)))

	       (gtk-box-pack-start hbox h-sep #f #f 2)
	       (gtk-box-pack-start hbox vbox #f #f 2)
	       (gtk-box-pack-start vbox go-button #t #f 0)
	       (gtk-box-pack-start vbox cancel-button #t #f 0)
	       (gtk-widget-show-all window)))))

	(add-simple-coot-menu-menuitem
	 menu "Read SHELX Project..."
	 (lambda ()
	   (let* ((window (gtk-window-new 'toplevel))
		  (vbox (gtk-vbox-new #f 2))
		  (hbox (gtk-hbox-new #t 2))
		  (h-sep (gtk-hseparator-new))
		  (go-button (gtk-button-new-with-label "  Read Project   "))
		  (cancel-button (gtk-button-new-with-label "  Cancel   ")))

	     (let ((entry (file-selector-entry vbox " Project Name: ")))

	       (gtk-signal-connect cancel-button "clicked"
				   (lambda ()
				     (gtk-widget-destroy window)))
	       (gtk-signal-connect go-button "clicked"
				   (lambda ()
				     (let ((file-name (gtk-entry-get-text entry)))
				       (read-shelx-project file-name)
				       (gtk-widget-destroy window))))
	       
	       (gtk-container-add window vbox)
	       
	       (gtk-box-pack-start hbox go-button #f #f)
	       (gtk-box-pack-start hbox cancel-button #f #f)
	       (gtk-box-pack-start vbox h-sep  #f #f)
	       (gtk-box-pack-start vbox hbox  #f #f)
	       (gtk-widget-show-all window)))))


	(add-simple-coot-menu-menuitem
	 menu "Read LST file..."
	 (lambda ()
	   (generic-chooser-and-file-selector 
	    "Model Corresponding to LST file: "
	    valid-model-molecule?
	    "LST file"
	    ""
	    (lambda (imol lst-file-name)
	      (read-shelx-lst-file lst-file-name imol)))))


	(add-simple-coot-menu-menuitem
	 menu "Add SHELXL instruction..."
	 (lambda ()
	   (generic-chooser-and-entry
	    "Add new SHELXL command to model:"
	    "SHELX instruction:"
	    ""
	    (lambda (imol text)
	      (add-shelx-string-to-molecule imol text)))))
	)))


(define (shelx-ins-strings imol)      

  (let ((ins-tmp-file "coot-tmp.ins"))
    (write-shelx-ins-file imol ins-tmp-file)
    (call-with-input-file ins-tmp-file
      (lambda (port)
	(let f ((line (read-line port))
		(lines '()))
	  (cond 
	   ((eof-object? line) (reverse lines))
	   (else
	    (f (read-line port) 
	       (cons line lines)))))))))

(define (shelxl-refine-gui imol . hkl-file-name-maybe)

  (let ((window (gtk-window-new 'toplevel))
	(vbox (gtk-vbox-new #f 2))
	(buttons-hbox (gtk-hbox-new #f 2))
	(scrolled-win (gtk-scrolled-window-new))
	(text (gtk-text-new #f #f))
	(refine-button (gtk-button-new-with-label "   Refine   "))
	(cancel-button (gtk-button-new-with-label "   Cancel   ")))
	
    (gtk-container-add window vbox)
    (gtk-container-add scrolled-win text)
    (gtk-box-pack-start buttons-hbox refine-button #f #f 2)
    (gtk-box-pack-start buttons-hbox cancel-button #f #f 2)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    (gtk-box-pack-start vbox buttons-hbox #f #f 2)
    (gtk-box-pack-start vbox scrolled-win)

    (let ((text-strings (shelx-ins-strings imol))
	  (bg-col "#c0e6c0"))
      (for-each
       (lambda (str) 
	 (gtk-text-insert text #f "black" bg-col str -1)
	 (gtk-text-insert text #f "black" bg-col "\n" -1))
       text-strings))

    (gtk-signal-connect refine-button "clicked"
			(lambda ()
			  
			  (let* ((l (gtk-text-get-length text))
				 (text-text (gtk-editable-get-chars text 0 l)))

			    (shelxl-refine-primitive imol text-text hkl-file-name-maybe)
			    (gtk-widget-destroy window))))

    (gtk-signal-connect cancel-button "clicked"
			(lambda () 
			  (gtk-widget-destroy window)))

    (gtk-text-set-editable text #t)
    (gtk-window-set-default-size window 500 500)
    (gtk-widget-show-all window)))


(define (editable-shelx-gui imol hklin-file-name)

  (let ((window (gtk-window-new 'toplevel))
	(text (gtk-text-new #f #f))
	(scrolled-win (gtk-scrolled-window-new))
	(cancel-button (gtk-button-new-with-label "  Close  "))
	(run-button (gtk-button-new-with-label "  Run  "))
	(vbox (gtk-vbox-new #f 0))
	(buttons-vbox (gtk-hbox-new #f 0)))
    
    (gtk-window-set-default-size window 450 400)
    (gtk-container-add window vbox)
    (gtk-container-add scrolled-win text)
    (gtk-container-border-width vbox 5)
    (gtk-box-pack-start buttons-vbox    run-button #t #f 2)
    (gtk-box-pack-start buttons-vbox cancel-button #t #f 2)
    (gtk-box-pack-start vbox buttons-vbox #f #f 0)
    (gtk-box-pack-start vbox scrolled-win #t #t 2)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    
    (let ((shelx-ins-list (get-shelx-ins-list imol))
	  (background-colour "#c0e6c0"))
      
      ;; now fill the text with the shelx-ins strings 
      ;;
      (for-each
       (lambda (str)
	 (gtk-text-insert text #f "black" background-colour str -1)
	 (gtk-text-insert text #f "black" background-colour "\n" -1))
       shelx-ins-list))

    ;; button press actions
    ;; 
    (gtk-signal-connect cancel-button "clicked" 
			(lambda () (gtk-widget-destroy window)))
    (gtk-signal-connect run-button "clicked"
			(lambda () 
			  (let* ((l (gtk-text-get-length text))
				 (text-text (gtk-editable-get-chars text 0 l)))
			    (if (string? text-text)
				(let ((hklin-file-info 
				       (if (= (string-length hklin-file-name) 0)
					   '()
					   (list hklin-file-name))))
				(shelxl-refine-primitive imol text-text
							 hklin-file-info))))
			  (gtk-widget-destroy window)))
    

    (gtk-text-set-editable text #t)
    (gtk-widget-show-all window)))


(define (get-shelx-ins-list imol)
  (let ((tmp-ins "tmp.ins"))
    (write-shelx-ins-file imol tmp-ins)
    (if (file-exists? tmp-ins)
	(call-with-input-file tmp-ins
	  (lambda (port)
	    (let f ((obj (read-line port))
		    (lines '()))
	      (cond
	       ((eof-object? obj) (reverse lines))
	       (else 
		(f (read-line port)
		   (cons obj lines))))))))))

	


;;; This is a top-level dialog that add an HFIX instruction to a given .res 
;;; molecule.  It is a specific, not generally very useful, as it stands.
;;;
;(let ((buttons 
;        (list 
;	(list
; 	 "Rebuild Hydrogens: HFIX 137 C2BA "
;	 (lambda ()
;	   (molecule-chooser-gui "Add rebuild hydrogens instructions"
;				 (lambda (imol)
;				   (let ((s "HFIX 137 C2BA"))
;				     (add-shelx-string-to-molecule imol s)))))))))
;  (dialog-box-of-buttons "SHELX Demo" (cons 200 50) buttons "  Close  "))



