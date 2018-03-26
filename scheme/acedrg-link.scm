
;; 20180309 verison 2.3

(define (acedrg-link-generation-control-window)

  ;; main body
  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 4))
	 (inside-hbox-1 (gtk-hbox-new #f 4))
	 (inside-hbox-2 (gtk-hbox-new #f 4))
	 (cancel-hbox (gtk-hbox-new #f 2))
	 (order-label  (gtk-label-new "  Order: "))
	 (delete-label (gtk-label-new "  Delete: "))
	 (other-label (gtk-label-new " from First Residue "))
	 (option-menu-order (gtk-option-menu-new))
	 (menu (gtk-menu-new))
	 (tte (gtk-tooltips-new))
	 (tto (gtk-tooltips-new))
	 (ttp (gtk-tooltips-new))
	 (entry (gtk-entry-new))
	 (h-sep (gtk-hseparator-new))
	 (cancel-button     (gtk-button-new-with-label "  Cancel  "))
	 (pick-button       (gtk-button-new-with-label "Start (Pick 2 Atoms)...")))

    (gtk-window-set-title window "Make a Link using Acedrg and Atom Click Click")
    (gtk-box-pack-start vbox inside-hbox-1 #f #f 2)
    (gtk-box-pack-start vbox inside-hbox-2 #f #f 2)
    (gtk-box-pack-start inside-hbox-1 order-label #f #f 2)
    (gtk-box-pack-start inside-hbox-1 option-menu-order #f #f 2)
    (gtk-box-pack-start inside-hbox-2 delete-label #f #f 2)
    (gtk-box-pack-start inside-hbox-2 entry #f #f 2)
    (gtk-box-pack-start inside-hbox-2 other-label #f #f 2)
    (gtk-box-pack-start vbox h-sep)
    (gtk-box-pack-start vbox cancel-hbox #f #f 6)
    (gtk-box-pack-start cancel-hbox   pick-button #f #f 6)
    (gtk-box-pack-start cancel-hbox cancel-button #f #f 6)
    (gtk-widget-set-usize entry 80 -1)

    (let* ((menu-item-1 (gtk-menu-item-new-with-label "Single"))
	   (menu-item-2 (gtk-menu-item-new-with-label "Double")))
      (gtk-menu-append menu menu-item-1)
      (gtk-menu-append menu menu-item-2)
      (gtk-menu-set-active menu 0))
    (gtk-option-menu-set-menu option-menu-order menu)

    (gtk-tooltips-set-tip tte entry "Type the atom name to be deleted" "")
    (gtk-tooltips-set-tip tte option-menu-order "Like a cheeseburger - you can only have single or double" "")
    (gtk-tooltips-set-tip ttp pick-button "Click on 2 atoms, Acedrg starts after the second click" "")
    (gtk-signal-connect cancel-button "clicked" (lambda () (gtk-widget-destroy window)))

    (gtk-signal-connect   pick-button "clicked" (lambda () (click-select-residues-for-acedrg window option-menu-order entry)))

    (gtk-container-add window vbox)
    (gtk-widget-show-all window)))
  

;; hack link cif file so that it works in coot 0.8.9 (:-/)
;;
;; return the new file name
(define (hack-link fn)
  (let ((stub (file-name-sans-extension fn)))
    (let ((new-file-name (string-append stub "-hack.cif")))
      (call-with-output-file new-file-name
	(lambda (port-out)
	  (call-with-input-file fn
	    (lambda (port-in)
	      (let loop ((line (read-line port-in)))
		(if (not (eof-object? line))
		    (begin
		      (if (string-match "L-PEPTIDE" line)
			  (let ((xx (coot-replace-string line "L-PEPTIDE" "L-peptide")))
			    (display xx port-out))
			  (display line port-out))
		      (newline port-out)
		      (loop (read-line port-in)))))))))
      new-file-name)))


(define (click-select-residues-for-acedrg window option-menu entry)

  (format #t "window: ~s~%" window)

  (user-defined-click 
   2
   (lambda (clicks)
     (format #t "we received these clicks: ~s~%" clicks)

     (let ((bond-order (get-option-menu-active-item option-menu (list 'single 'double))))
       (if (= (length clicks) 2)
	   (let ((click-1 (list-ref clicks 0))
		 (click-2 (list-ref clicks 1)))
	     (format #t "click-1: ~s~%" click-1)
	     (format #t "click-2: ~s~%" click-2)
	     (if (and (= (length click-1) 7)
		      (= (length click-2) 7))
		 (let ((resname-1 (residue-name 
				   (list-ref click-1 1)
				   (list-ref click-1 2)
				   (list-ref click-1 3)
				   (list-ref click-1 4)))
		       (resname-2 (residue-name
				   (list-ref click-2 1)
				   (list-ref click-2 2)
				   (list-ref click-2 3)
				   (list-ref click-2 4)))
		       (at-name-1 (list-ref click-1 5))
		       (at-name-2 (list-ref click-2 5))
		       (spec-1 (cddr click-1))
		       (spec-2 (cddr click-2))
		       (imol-click-1 (list-ref click-1 1))
		       (imol-click-2 (list-ref click-2 1))
		       (delete-atom-text (gtk-entry-get-text entry)))

		   (if (not (and (string? resname-1)
				 (string? resname-2)))
		       (begin
			 (format #t "Bad resnames: ~s and ~s~%"
				 resname-1 resname-2))
		       (begin
			 (if (not (= imol-click-1 imol-click-2))

			     #f

			     (let* ((imol imol-click-1)
				    (delete-stripped-1 (strip-spaces delete-atom-text))
				    (delete-atom-txt (if (> (string-length delete-stripped-1) 0) 
							 (string-append " DELETE ATOM " delete-stripped-1 " 1 ")
							 ""))
				    (s (string-append
					"LINK:"
					" RES-NAME-1 " resname-1 " ATOM-NAME-1 " at-name-1 
					" RES-NAME-2 " resname-2 " ATOM-NAME-2 " at-name-2)))

			       ;; I need to check here if resname-1 or resname-2 came from a file that 
			       ;; was read into Coot from somewhere other than the refmac monomer library
			       ;; (that acedrg knows about).

			       (let ((cif-fn-1 (cif-file-for-comp-id resname-1))
				     (cif-fn-2 (cif-file-for-comp-id resname-2))
				     (ns (non-standard-residue-names imol)))

				 ;; if the resnames are not non-standard-residue-names
				 ;; then we don't need to specify the file - if they
				 ;; are not, then use cif-fn-1 (or cif-fn-2).

				 (format #t "cif-fn-1: ~s~%" cif-fn-1)
				 (format #t "cif-fn-2: ~s~%" cif-fn-2)
				 (format #t "ns:       ~s~%" ns)

				 (if (eq? bond-order 'double)
				     (set! s (string-append s " BOND-TYPE DOUBLE")))

				 (if (string-member? resname-1 ns)
				     (set! s (string-append s " FILE-1 " cif-fn-1)))
				 (if (string-member? resname-2 ns)
				     (set! s (string-append s " FILE-2 " cif-fn-2)))

				 ;; and finally delete atom
				 (set! s (string-append s delete-atom-txt))

				 (format #t "LINK string: ~s~%" s)
				 (let* ((st-1 (string-append "acedrg-link-from-coot-"
							     resname-1 "-"
							     resname-2))
					(st (string-append st-1 "-link-instructions"))
					(log-file-name (string-append st ".log"))
					(ins-file-name (string-append st ".txt")))
				   (call-with-output-file ins-file-name
				     (lambda (port) (display s port) (newline port)))
				   (let ((goosh-status (goosh-command "acedrg"
								      (list "-L" ins-file-name "-o" st-1)
								      '()
								      log-file-name
								      #f)))
				     (if (not (ok-goosh-status? goosh-status))

					 (let ((m (string-append "WARNING:: acedrg failed.\nSee " log-file-name)))
					   (info-dialog m))

					 ;; happy path
					 (let ((link-file-name (string-append st-1 "_link.cif"))) ;; acedrg name

					   (let ((hack-link-file-name (hack-link link-file-name)))

					     (let ((dict-read-status (read-cif-dictionary hack-link-file-name)))
					       ;; dict-read-status is the number of bonds read
					       (format #t "dict-read-status: ~s~%" dict-read-status)
					       (if (> dict-read-status 0)
						   (make-link imol-click-1 spec-1 spec-2 "dummy-name" 1.0))))))))))))))))))

     (gtk-widget-destroy window))))




;;; beam in a new menu to the menu bar:
;(if (defined? 'coot-main-menubar)
;    (let* ((new-menu (coot-menubar-menu "Acedrg")))
      
;      (add-simple-coot-menu-menuitem 
;       new-menu "Make Link with Link Dictionary from Acedrg"
;       (lambda ()
;	 (acedrg-link-generation-control-window)))))

