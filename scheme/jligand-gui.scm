
(define *jligand-home-env* #f)

;; This happens when user clicks on the "Launch JLigand" button.
;; It starts a jligand and puts it in the background.
;; 
(define (launch-jligand-function)

  (start-jligand-listener)
  (if (not (command-in-path? "jligand"))

      ;; Boo.  Give us a warning dialog
      ;; 
      (let ((s (string-append "jligand not found in path"))
	    ;; make an extra message telling us that JLIGAND_HOME is
	    ;; not set if it is not set.
	    (env-message (if (string? *jligand-home-env*) 
			     ""
			     (string-append "Environment variable JLIGAND_HOME not set\n\n"))))
	(info-dialog (string-append env-message s)))
      
      ;; OK, it does exist - run it!
      ;;
      (let ((s "jligand -coot &"))
	(system s)
	;; beam in a new menu to the menu bar:
	(let* ((jligand-menu (coot-menubar-menu "JLigand")))

	  (add-simple-coot-menu-menuitem 
	   jligand-menu "Send Link to JLigand (click 2 monomers)"
	   (lambda ()
	     (click-select-residues-for-jligand)))))))




;; This happens when user clicks on the "Select Residues for JLigand"
;; (or some such) button.  It expects the user to click on atoms of
;; the two residues involved in the link.
;; 
(define (click-select-residues-for-jligand)
  (user-defined-click 
   2
   (lambda (clicks)
     (format #t "we received these clicks: ~s~%" clicks)
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
		     (imol-click-1 (list-ref click-1 1))
		     (imol-click-2 (list-ref click-2 1)))
		 (if (not (and (string? resname-1)
			       (string? resname-2)))
		     (begin 
		       (format #t "Bad resnames: ~s and ~s~%"
			       resname-1 resname-2))
		     (begin
		       (if (not (= imol-click-1 imol-click-2))
			   (begin
			     (set! *imol-jligand-link* #f))
			   (begin ;; happy path
			     (set! *imol-jligand-link* imol-click-1)
			     (write-file-for-jligand (click->res-spec click-1) resname-1 
						     (click->res-spec click-2) resname-2))))))))))))

