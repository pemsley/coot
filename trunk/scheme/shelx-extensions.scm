;;;; Copyright 2007 by The University of York
;;;; Copyright 2007 by Paul Emsley
;;;; Copyright 2007 by The University of Oxford
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

(if (defined? 'coot-main-menubar)

    (let ((menu (coot-menubar-menu "SHELX")))

      (add-simple-coot-menu-menuitem
       menu "SHELXL Refine..."
       (lambda ()

	 (let ((window (gtk-window-new 'toplevel))
	       (hbox (gtk-vbox-new #f 0))
	       (vbox (gtk-hbox-new #f 0))
	       (go-button (gtk-button-new-with-label "  Refine  "))
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
					 (if (= (string-length txt) 0)
					     (shelxl-refine imol)
					     (shelxl-refine imol txt)))
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
	  "LST file"
	  ""
	  (lambda (imol lst-file-name)
	    (read-shelx-lst-file lst-file-name imol)))))

))
      

;; This is a top-level dialog that add an HFIX instruction to a given .res 
;; molecule.  It is a specific, not generally very useful, as it stands.
;;
(let ((buttons 
       (list 
	(list
	 "Rebuild Hydrogens: HFIX 137 C2BA "
	 (lambda ()
	   (molecule-chooser-gui "Add rebuild hydrogens instructions"
				 (lambda (imol)
				   (let ((s "HFIX 137 C2BA"))
				     (add-shelx-string-to-molecule imol s)))))))))
  (dialog-box-of-buttons "SHELX Demo" (cons 200 50) buttons "  Close  "))



