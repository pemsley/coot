
;;;; Copyright 2011 by The University of Oxford

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
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


;; (load "json-reader.scm") ;; loaded by coot.scm

(use-modules (json reader)
	     (ice-9 format)
	     (ice-9 q))


(define *coot-pdbe-image-cache-dir* "coot-pdbe-images") ;; can be shared (dir should
                                                        ;; be writable by sharers).

(define coot-thread-dispatcher 
  (let* ((q-locked #f)
	 (q (make-q))
	 (n-threads 0)
	 (n-threads-locked #f)
	 (get-q-lock (lambda () 
		       (while q-locked (begin 
					 (usleep (random 200)))) 
		       (set! q-locked #t)))
	 (release-q-lock (lambda () (set! q-locked #f)))
	 (get-n-threads-lock (lambda () 
			       (while n-threads-locked 
				      (begin 
					(usleep (random 200))))
			       (set! n-threads-locked #t)))
	 (release-n-threads-lock (lambda () (set! n-threads-locked #f)))
	 (in-queue? (lambda (file-name)
		      (if (q-empty? q)
			  (begin 
			    #f)

			  (let loop ((i 0)
				     (n (q-length q)))
			    (if (= i n) 
				(begin
				  #f) ;; got to end without finding string
				 
				(begin
				  (let* ((pv (car q))
					 (v (list-ref pv i)))
				    (if (string=? (car v) file-name)
					(begin 
					  #t)
					(begin
					  (loop (+ i 1) n))))))))))
		    
	 (timeout 
	  ;; pull off items from the queue and run them
	  ;; 
	  (gtk-timeout-add 400
			   (lambda ()
			     (if (= (q-length q) 0)
				 (begin
				   ;; (format #t "nothing in queue\n")
				   #t ;; nothing in queue
				   )
				 (begin
				   
				   (while (and (< n-threads 100) 
					       (not (= (q-length q) 0)))
					  
					  (begin
					    (get-q-lock)
					    (let ((q-image-name-and-func (q-pop! q))) ;; changes q
					      (release-q-lock)
					      
					      (let ((my-func
						     (lambda ()
						       (get-n-threads-lock)
						       (set! n-threads (+ n-threads 1))
						       (release-n-threads-lock)
						       ((cdr q-image-name-and-func))
						       (get-n-threads-lock)
						       (set! n-threads (- n-threads 1))
						       (release-n-threads-lock)
						       #t)))
						(call-with-new-thread my-func)))))))
			     #t)))) ;; lots of threads already, pass for this round.

    (lambda (image-name func)
      ;; add the func to the queue

      (let ((queue-status (in-queue? image-name)))
	(if (not queue-status)
	    (begin 
	      ;; (format #t ":::::::: ~s was not in queue already, adding it~%" image-name)
	      (get-q-lock)
	      (q-push! q (cons image-name func))
	      (release-q-lock)
	      'dummy)
	    (begin 
	      ;; (format #t ":::::::: ~s WAS in queue already, skipping...~%" image-name)
	      'previous-thread-getting ;; let caller know that a thread to get this image is in the queue
	      ))))))

	         

;;; 
(define (get-recent-json file-name)
  
  (if (not (file-exists? file-name))

      (begin
	(format #t "file not found: ~s~%" file-name)
	#f)
      
      (begin
	(call-with-input-file file-name
	  (lambda (port)
	    (let ((s (json:read-value port)))
	      s))))))

;; geometry is an improper list of ints.
;; 
;; return the h-box of the buttons.
;;
;; a button is a list of (label callback-thunk text-description)
;;
;; If check-button-label is #f, don't make one, otherwise create with with
;; the given label and "on" state.
;; 
(define (dialog-box-of-buttons-with-async-ligands window-name geometry buttons close-button-label)

  ;; the async function is evaluated here
  ;; 
  (define (add-button-info-to-box-of-buttons-vbox-for-ligand-images button-info vbox)
    (let* ((button-label (car button-info))
	   (callback (car (cdr button-info)))
	   (active-button-label-func
	    (if (= (length button-info) 2)
		#f ; it doesn't have one
		(list-ref button-info 2)))
	   (button (gtk-button-new))
	   (button-hbox         (gtk-hbox-new #f 0)) ;; pixmaps and labels get added here.
	   (protein-ribbon-hbox (gtk-hbox-new #f 0))
	   (ligands-hbox        (gtk-hbox-new #f 0)))
      
      (gtk-box-pack-start button-hbox (gtk-label-new button-label) #f 0)
      (gtk-box-pack-start button-hbox ligands-hbox #f 0)
      (gtk-box-pack-start button-hbox protein-ribbon-hbox #f 0)
      (gtk-container-add button button-hbox)
      (gtk-signal-connect button "clicked" callback)

      (if (procedure? active-button-label-func)
	  (active-button-label-func button-hbox ligands-hbox protein-ribbon-hbox))

      (gtk-box-pack-start vbox button #f #f 2)
      (gtk-widget-show button)))


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

    (gtk-box-pack-start outside-vbox scrolled-win #t #t 0) ; expand fill padding
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

    (map (lambda (button-info)
	   (add-button-info-to-box-of-buttons-vbox-for-ligand-images button-info inside-vbox))
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
      inside-vbox))


;; Get image-name (caller doesn't care how) and when it is in place run func.
;; 
;; This is a generally useful function, so it has been promoted outside of dig-table.
;; 
(define (cache-or-net-get-image image-url image-name func)

  (define show-image-when-ready
    (lambda ()
      (if (file-exists? image-name)
	  (begin
	    (func)
	    #f)
	  (begin 
	    #t))))

  
  (let ((curl-status 'start)
	(part-image-name (string-append image-name ".part")))
	  (if (and (file-exists? image-name)
		   (not (file-exists? part-image-name))
		   (> (stat:size (stat image-name)) 0))

	      ;; already here.
	      (begin
		(func))
	      
	      (begin
		;; OK, so the file was not already there... but we
		;; don't want to do this if there is already a
		;; download ongoing for this file/image (if the
		;; progress is non-empty then some other thread got
		;; there first).
		;; 
		(let ((image-curl-progress-info (curl-progress-info part-image-name)))

		  (if (not (eq? image-curl-progress-info #f))

		      ;; some other thread got there first.  Let's wait
		      ;; for that to put the image in place.
		      (begin
			(format #t "----------------- another thread already getting ~s~%" image-name)
			(gtk-timeout-add 600 show-image-when-ready))

		      ;; OK, this thread should get the image.
		      (begin
			;; (format #t "dispatch::: let's queue thread to get ~s ~s~%" image-url image-name)
			(let ((prev-thread-getting-status 
			       (coot-thread-dispatcher 
				part-image-name
				(lambda ()
				  (let* ((status (coot-get-url-and-activate-curl-hook 
						  image-url part-image-name 1)))
				    (if (not (= status 0))
					(begin
					  ;; (format #t "DEBUG:: curl status getting ~s is ~s~%" image-name status)
					  (set! curl-status 'fail))
					(begin
					  ;; (format #t "DEBUG:: HAPPY curl status getting ~s is ~s~%" image-name status)
					  (if (not (file-exists? image-name))
					      (if (file-exists? part-image-name)
						  (rename-file part-image-name image-name)))
					  (set! curl-status func))))))))

			  
			  (if (eq? prev-thread-getting-status 'previous-thread-getting)

			      ;; wait for other thread to get it
			      (gtk-timeout-add 600 show-image-when-ready)

			      ;; wait for this thread to get it
			      (gtk-timeout-add 600
					       (lambda ()
						 (cond 
						  ((procedure? curl-status)
						   (curl-status)
						   #f)
						  ((eq? curl-status 'fail)
						   #f) ; stop
						  (else 
						   #t ;; continue
						   )))))))))))))


;; return refmac-result or #f 
;; 
(define (refmac-calc-sfs-make-mtz  pdb-in-file-name mtz-file-name mtz-refmaced-file-name)
  
  (let* ((refmac-stub (append-dir-file "coot-refmac"
				       (strip-path 
					(file-name-sans-extension pdb-in-file-name))))
	 (pdb-out-file-name (string-append refmac-stub "-refmaced.pdb"))
	 (mtz-out-file-name mtz-refmaced-file-name)
	 (extra-cif-lib-filename "")
	 (imol-refmac-count 0)
	 (swap-map-colours-post-refmac? 0)
	 (imol-mtz-molecule #f)
	 (show-diff-map-flag 0)
	 (phase-combine-flag 0)
	 (phib-fom-pair #f)
	 (force-n-cycles 0)
	 (make-molecules-flag 0)) ;; don't make molecules, this may be a sub-thread.
        
    (let ((save-refmac-extra-params refmac-extra-params))
      
      (if (list? refmac-extra-params)
	  (set! refmac-extra-params (cons "MAKE NEWLIGAND CONTINUE") refmac-extra-params)
	  (set! refmac-extra-params (list "MAKE NEWLIGAND CONTINUE")))

      (let ((refmac-result
	     (run-refmac-by-filename pdb-in-file-name 
				     pdb-out-file-name 
				     mtz-file-name mtz-out-file-name
				     extra-cif-lib-filename 
				     imol-refmac-count
				     swap-map-colours-post-refmac? 
				     imol-mtz-molecule show-diff-map-flag
				     phase-combine-flag 
				     phib-fom-pair 
				     force-n-cycles
				     make-molecules-flag
				     "" "F.F_sigF.F" "F.F_sigF.sigF" "Rfree.Flag.flag")))
	                             ;; ccp4i-project-dir f-col sig-f-col . r-free-col

	;; restore refmac-extra-params to what it used to be
	;;
	(set! refmac-extra-params save-refmac-extra-params)
      
	(if (file-exists? mtz-refmaced-file-name)
	    refmac-result
	    #f)))))
    


(define (update-download-progress-bar-maybe curl-progress-info progress-bar)

  ;; (format #t "------------- curl-progress-info: ~s~%" curl-progress-info)
  (if curl-progress-info
      (let ((v1 (assoc 'content-length-download curl-progress-info))
	    (v2 (assoc 'size-download           curl-progress-info)))
	;; v1 and v2 are pairs but not lists.
	
	(if (pair? v1)
	    (if (pair? v2)
		(let ((f2 (cdr v2))
		      (f1 (cdr v1)))
		  (if (> f1 0) ;; content-length-download is -1 (or so).
		      (if (> f2 0)
			  (let ((f (/ f2 f1)))
			    (begin
			      ;; (format #t "---- update ~s to ~s/~s = ~s~%" progress-bar f2 f1 f)
			      (gtk-progress-bar-update progress-bar f))))
		      
		      (if (> f2 0)

			  ;; OK, the web server didn't tell us how big
			  ;; the file is (content-length-download) but
			  ;; we have downloaded something.  So, we
			  ;; should pulse.
			  ;; 
			  ;; But I don't have that in guile-gtk - boo.
			  ;; 
			  ;; So hack in something slightly better than nothing.
			  ;; 
			  (let ((f (/ f1 227085.0))) ;; some made up sensiblish number
			    (if (<= f 1)
				(begin
				  ;; (format #t "setting activity mode....~%")
				  (gtk-progress-set-activity-mode progress-bar f))))))))))))

;;; 
(define (dig-table t level)

  ;; is this a list of hash tables?
  ;; 
  (define (contains-pdb-entry-entities? ls)

    (define (non-null-list-of-hash-tables? ls)
      (if (not (list? ls))
	  #f
	  (if (null? ls)
	      #f
	      (let loop ((lsl ls))
		(cond
		 ((null? lsl) 
		  #t)
		 ((not (hash-table? (car lsl))) #f)
		 (else
		  (loop (cdr lsl))))))))

    (if (not (non-null-list-of-hash-tables? ls))
	#f ;; is null :)

	;; OK, a proper entry for a PDB had and EntryID
	(hash-ref (car ls) "EntryID")))

  ;;
  ;; 
  (define (files-exist? file-list)
    (cond
     ((null? file-list) #t)
     ((not (file-exists? (car file-list))) #f)
     (else
      (files-exist? (cdr file-list)))))


  ;; 
  ;; 
  (define (make-and-draw-map-local refmac-out-mtz-file-name)

    (make-and-draw-map refmac-out-mtz-file-name "FWT" "PHWT" "" 0 0)
    (make-and-draw-map refmac-out-mtz-file-name "DELFWT" "PHDELWT" "" 0 1))


	

  ;; Activate a progress bar here, I think.
  ;; 
  ;; include-get-sfs-flag is either 'no-sfs or 'include-sfs
  ;; 
  (define (pdbe-get-pdb-and-sfs-cif include-get-sfs-flag entry-id method-string)

    (let ((status (make-directory-maybe "coot-download")))
      (if (not (= status 0))
	  (info-dialog "Failed to make download directory")

	  ;; do it!
	  ;; 
	  (let ((curl-status 'start))

	    ;; we touch curl-status, that's why this is here.
	    ;; 
	    (define (get-sfs-run-refmac sfs-cif-url sfs-cif-file-name sfs-mtz-file-name pdb-file-name refmac-out-mtz-file-name)

	      ;; 
	      (define (convert-to-mtz-and-refmac sfs-cif-file-name sfs-mtz-file-name pdb-file-name)

		;; OK, let's run convert to mtz and run refmac
		;; 
		(set! curl-status 'converting-to-mtz)
		(let ((convert-status (mmcif-sfs-to-mtz sfs-cif-file-name sfs-mtz-file-name)))
		  (if (not (= convert-status 1))
		      ;; 
		      (begin
			;; we can't make a dialog of course
			(format #t "WARNING:: Failed to convert ~s to an mtz file~%" sfs-cif-file-name)
			(set! curl-status 'fail))
		      
		      (begin
			(refmac-inner pdb-file-name sfs-mtz-file-name refmac-out-mtz-file-name)))))

	      ;; 
	      ;; 
	      (define (refmac-inner pdb-file-name sfs-mtz-file-name refmac-out-mtz-file-name)

		(set! curl-status 'running-refmac-for-phases)
		(let ((refmac-result
		       (refmac-calc-sfs-make-mtz pdb-file-name 
						 sfs-mtz-file-name
						 refmac-out-mtz-file-name)))

;		  ;; if refmac-result is good?
		  (if (not (list? refmac-result))

		      (begin
			(set! curl-status 'fail-refmac))
		      
		      (begin
			(set! curl-status 
			      (lambda ()
				(let ((imol (read-pdb pdb-file-name)))
				  (if (not (valid-model-molecule? imol))
				      (let ((s  (string-append
						 "Oops - failed to correctly read "
						 pdb-file-name)))
					(info-dialog s))
				      
				      (make-and-draw-map-local refmac-out-mtz-file-name)))))))))



	      ;; main line of get-sfs-run-refmac
	      ;; 
;	      (format #t "in get-sfs-run-refmac ~s ~s ~s ~s~%"
;		      sfs-cif-file-name
;		      sfs-mtz-file-name
;		      pdb-file-name
;		      refmac-out-mtz-file-name)
	      
	      ;; check for cached results: Only run refmac if the
	      ;; output file does not exist or is empty.
	      ;; 
	      (if (not (or (not (file-exists? refmac-out-mtz-file-name))
			   (= (stat:size (stat refmac-out-mtz-file-name)) 0)))
		  
		  (begin
		    ;; using cached result
		    (lambda ()
		      (let ((imol (read-pdb pdb-file-name)))
			(if (not (valid-model-molecule? imol))
			    (let ((s  (string-append
				       "Oops - failed to correctly read "
				       pdb-file-name)))
			      (info-dialog s)))
			(make-and-draw-map refmac-out-mtz-file-name "FWT" "PHWT" "" 0 0)
			(make-and-draw-map refmac-out-mtz-file-name "DELFWT" "PHDELWT" "" 0 1))))
		      

		  ;; the coot refmac interface writes its
		  ;; output in coot-refmac directory.  If
		  ;; that doesn't exist and we can't make
		  ;; it, then give up.
		  ;; 
		  (if (not (= (make-directory-maybe "coot-refmac") 0))
		      (begin 
			(info-dialog "Can't make output directory coot-refmac"))
		      (begin
			(if (not (or (not (file-exists? sfs-cif-file-name))
				     (= (stat:size (stat sfs-cif-file-name)) 0)))
			    (begin
			      ;; OK we have sfs-cif-file-name already
			      (format #t ".... path 3~%")
			      (convert-to-mtz-and-refmac sfs-cif-file-name sfs-mtz-file-name pdb-file-name))

			    (begin
			      ;; need to get sfs-mtz-file-name 
			      ;; 
			      (set! curl-status 'downloading-sfs)
			      (let ((status-cif (coot-get-url-and-activate-curl-hook 
						 sfs-cif-url sfs-cif-file-name 1)))

				(if (not (= status-cif 0))
				    (begin
				      (format #t "failed to download cif ~s~%" sfs-cif-file-name)
				      (set! curl-status 'fail))

				    (begin
				      (convert-to-mtz-and-refmac sfs-cif-file-name sfs-mtz-file-name pdb-file-name))))))))))

	    (define (get-widget progress-widgets widget-symbol)
	      (let ((v (assoc-ref progress-widgets widget-symbol)))
		v))

	    ;; return a list of the progress bars and the window 
	    ;; 
	    ; (the pdb-file-name and sfs-cif-file-name are passed so
	    ; that the cancel button knows what transfers to cancel (if
	    ;; needed)).
	    ;; 
	    (define (progress-dialog pdb-file-name sfs-cif-file-name)
	      (let* ((window (gtk-window-new 'toplevel))
		     (dialog-name (string-append "Download and make SFS for " entry-id))
		     (main-vbox (gtk-vbox-new #f 6))
		     (cancel-button (gtk-button-new-with-label "  Cancel  "))
		     (buttons-hbox (gtk-hbox-new #f 2))
		     (pdb-hbox (gtk-hbox-new #f 6))
		     (cif-hbox (gtk-hbox-new #f 6))
		     (refmac-hbox (gtk-hbox-new #f 6))
		     (pdb-label (gtk-label-new "Download PDB:  "))
		     (cif-label (gtk-label-new "Download SFs cif:"))
		     (refmac-label (gtk-label-new "Running Refmac:"))
		     (refmac-fail-label (gtk-label-new "Running Refmac Failed"))
		     (fail-label        (gtk-label-new "Something Went Wrong"))

		     (   pdb-progress-bar (gtk-progress-bar-new))
		     (   cif-progress-bar (gtk-progress-bar-new))
		     (refmac-progress-bar (gtk-progress-bar-new))

		     (   pdb-execute-icon (gtk-image-new-from-stock "gtk-execute" 1))
		     (   cif-execute-icon (gtk-image-new-from-stock "gtk-execute" 1))
		     (refmac-execute-icon (gtk-image-new-from-stock "gtk-execute" 1))
		     (   pdb-good-icon (gtk-image-new-from-stock "gtk-ok" 1))
		     (   cif-good-icon (gtk-image-new-from-stock "gtk-ok" 1))
		     (refmac-good-icon (gtk-image-new-from-stock "gtk-ok" 1))
		     (   pdb-fail-icon (gtk-image-new-from-stock "gtk-no" 1))
		     (   cif-fail-icon (gtk-image-new-from-stock "gtk-no" 1))
		     (refmac-fail-icon (gtk-image-new-from-stock "gtk-no" 1))
		     (h-sep (gtk-hseparator-new)))
		     
		    
		(gtk-window-set-title window dialog-name)
		(gtk-box-pack-start buttons-hbox cancel-button #t #f 2)
		(gtk-box-pack-start pdb-hbox       pdb-label #t #f 2)
		(gtk-box-pack-start cif-hbox       cif-label #t #f 2)
		(gtk-box-pack-start refmac-hbox refmac-label #t #f 2)
		(gtk-box-pack-start pdb-hbox       pdb-progress-bar #t #f 3)
		(gtk-box-pack-start cif-hbox       cif-progress-bar #t #f 3)
		(gtk-box-pack-start refmac-hbox refmac-progress-bar #t #f 3)
		(gtk-box-pack-start pdb-hbox       pdb-execute-icon #f #f 2)
		(gtk-box-pack-start cif-hbox       cif-execute-icon #f #f 2)
		(gtk-box-pack-start refmac-hbox refmac-execute-icon #f #f 2)
		(gtk-box-pack-start pdb-hbox       pdb-good-icon #f #f 2)
		(gtk-box-pack-start cif-hbox       cif-good-icon #f #f 2)
		(gtk-box-pack-start refmac-hbox refmac-good-icon #f #f 2)
		(gtk-box-pack-start pdb-hbox       pdb-fail-icon #f #f 2)
		(gtk-box-pack-start cif-hbox       cif-fail-icon #f #f 2)
		(gtk-box-pack-start refmac-hbox refmac-fail-icon #f #f 2)
		
		(gtk-box-pack-start main-vbox pdb-hbox     #t #f 4)
		(gtk-box-pack-start main-vbox cif-hbox     #t #f 4)
		(gtk-box-pack-start main-vbox refmac-hbox  #t #f 4)
		(gtk-box-pack-start main-vbox refmac-fail-label #t #f 2)
		(gtk-box-pack-start main-vbox fail-label   #t #f 2)
		(gtk-box-pack-start main-vbox h-sep        #t #f 4)
		(gtk-box-pack-start main-vbox buttons-hbox #t #f 4)
		(gtk-container-border-width main-vbox 6)
		
		(gtk-container-add window main-vbox)
		(gtk-container-border-width window 4)
		(gtk-widget-show-all window)
		(gtk-widget-hide    pdb-good-icon)
		(gtk-widget-hide    cif-good-icon)
		(gtk-widget-hide refmac-good-icon)
		(gtk-widget-hide    pdb-fail-icon)
		(gtk-widget-hide    cif-fail-icon)
		(gtk-widget-hide refmac-fail-icon)
		(gtk-widget-hide refmac-fail-label)
		(gtk-widget-hide fail-label)

		(gtk-widget-set-sensitive    pdb-execute-icon #f)
		(gtk-widget-set-sensitive    cif-execute-icon #f)
		(gtk-widget-set-sensitive refmac-execute-icon #f)

		(gtk-progress-set-show-text    pdb-progress-bar #t)
		(gtk-progress-set-show-text    cif-progress-bar #t)
		(gtk-progress-set-show-text refmac-progress-bar #t)


		(gtk-signal-connect cancel-button "clicked"
				    (lambda ()
				      (stop-curl-download pdb-file-name)
				      (stop-curl-download sfs-cif-file-name)
				      (set! curl-status 'fail)
				      (gtk-widget-destroy window)))
		
		;; return these
		(list 
		 (cons 'pdb-progress-bar pdb-progress-bar)
		 (cons 'cif-progress-bar cif-progress-bar)
		 (cons 'refmac-progress-bar refmac-progress-bar)
		 (cons 'window window)
		 (cons    'refmac-fail-label   refmac-fail-label)
		 (cons    'fail-label          fail-label)
		 (cons    'pdb-execute-icon    pdb-execute-icon)
		 (cons    'cif-execute-icon    cif-execute-icon)
		 (cons 'refmac-execute-icon refmac-execute-icon)
		 (cons    'pdb-good-icon    pdb-good-icon)
		 (cons    'cif-good-icon    cif-good-icon)
		 (cons 'refmac-good-icon refmac-good-icon)
		 (cons    'pdb-fail-icon    pdb-fail-icon)
		 (cons    'cif-fail-icon    cif-fail-icon)
		 (cons 'refmac-fail-icon refmac-fail-icon))))

	    ;; ----------------------------------------
	    ;; 
	    (let* ((pdb-url (string-append 
			     "http://www.ebi.ac.uk/pdbe-srv/view/files/"
			     entry-id ".ent"))
		   (sfs-cif-url (string-append
				 "http://www.ebi.ac.uk/pdbe-srv/view/files/r"
				 entry-id "sf.ent"))
		   (pdb-file-name (append-dir-file "coot-download" (string-append entry-id ".ent")))
		   (sfs-cif-file-name (append-dir-file "coot-download" (string-append "r" entry-id "sf.cif")))
		   (sfs-mtz-file-name (append-dir-file "coot-download" (string-append "r" entry-id "sf.mtz")))
		   (refmac-out-mtz-file-name 
		    (append-dir-file "coot-download" (string-append "r" entry-id "-refmac.mtz")))
		   (refmac-log-file-name (string-append "refmac-from-coot-" 
							(number->string refmac-count)
							".log")) ;; set in run-refmac-by-filename
		   (progr-widgets (progress-dialog pdb-file-name sfs-cif-file-name)) 
		   (window              (get-widget progr-widgets    'window))
		   (pdb-progress-bar    (get-widget progr-widgets    'pdb-progress-bar))
		   (cif-progress-bar    (get-widget progr-widgets    'cif-progress-bar))
		   (refmac-progress-bar (get-widget progr-widgets 'refmac-progress-bar))
		   (fail-label          (get-widget progr-widgets 'fail-label))
		   (pdb-execute-icon    (get-widget progr-widgets    'pdb-execute-icon))
		   (cif-execute-icon    (get-widget progr-widgets    'cif-execute-icon))
		   (refmac-execute-icon (get-widget progr-widgets 'refmac-execute-icon))
		   (pdb-good-icon     (get-widget progr-widgets    'pdb-good-icon))
		   (cif-good-icon     (get-widget progr-widgets    'cif-good-icon))
		   (refmac-good-icon  (get-widget progr-widgets 'refmac-good-icon))
		   (pdb-fail-icon     (get-widget progr-widgets    'pdb-fail-icon))
		   (cif-fail-icon     (get-widget progr-widgets    'cif-fail-icon))
		   (refmac-fail-icon  (get-widget progr-widgets 'refmac-fail-icon))
		   (refmac-fail-label (get-widget progr-widgets 'refmac-fail-label)))
	      

	      (if (file-exists? refmac-log-file-name)
		  (delete-file refmac-log-file-name))

	      (gtk-timeout-add 200 
			       (lambda ()

				 (define (update-refmac-progress-bar refmac-progress-bar log-file-name)
				   ;; refmac puts out 100 lines of text before it starts running. 
				   ;; Let's not count those as progress of the computation (otherwise 
				   ;; we jump to 22% after a fraction of a second).
				   ;; 
				   (let ((max-lines 350) ;; thats 450 - 100
					 (n-lines (file-n-lines log-file-name)))
				     (if (number? n-lines)
					 (let ((n-lines-rest (- n-lines 100)))
					   (if (> n-lines-rest 0)
					       (let ((f (/ n-lines-rest max-lines)))
						 (format #t "refmac progress bar update to ~s/~s = ~s~%"
							 n-lines-rest max-lines f)
						 (if (< f 1)
						     (gtk-progress-bar-update refmac-progress-bar f))))))))
				     
				 
				 
				 ;; main line of timeout
				 ;; 
				 (let ((pdb-curl-progress-info (curl-progress-info pdb-file-name))
				       (cif-curl-progress-info (curl-progress-info sfs-cif-file-name)))
				   
;				   (format #t "curl-status: ~s~%" curl-status)
;				   (format #t "pdb info: ~s~%" pdb-curl-progress-info)
;				   (format #t "sfs info: ~s~%" cif-curl-progress-info)

				   (if (eq? curl-status 'downloading-pdb)
				       (begin
					 (gtk-widget-set-sensitive pdb-execute-icon #t)))

				   (if (eq? curl-status 'converting-to-mtz)
				       (begin
					 (gtk-progress-bar-update pdb-progress-bar 1)
					 (gtk-progress-bar-update cif-progress-bar 1)))

				   (if (eq? curl-status 'downloading-sfs)
				       (begin
					 (gtk-progress-bar-update pdb-progress-bar 1)
					 (gtk-widget-hide pdb-execute-icon)
					 (gtk-widget-show pdb-good-icon)
					 (gtk-widget-set-sensitive cif-execute-icon #t)))

				   (if (eq? curl-status 'running-refmac-for-phases)
				       (begin
					 (gtk-progress-bar-update pdb-progress-bar 1)
					 (gtk-progress-bar-update cif-progress-bar 1)
					 (gtk-widget-set-sensitive refmac-execute-icon #t)
					 (gtk-widget-hide cif-execute-icon)
					 (gtk-widget-show cif-good-icon)
					 (update-refmac-progress-bar refmac-progress-bar refmac-log-file-name)))
								  
				   (update-download-progress-bar-maybe pdb-curl-progress-info pdb-progress-bar)
				   (update-download-progress-bar-maybe cif-curl-progress-info cif-progress-bar))

				 (cond 

				  ((eq? curl-status 'fail-refmac)
				   (gtk-widget-show refmac-fail-label)
				   (gtk-widget-hide refmac-execute-icon)
				   (gtk-widget-show refmac-fail-icon)
				   #f ;; don't continue
				   )
				  
				  ;; generic fail (don't turn off execute icons because we 
				  ;; don't know *what* failed. (Not so useful).
				  ;; 
				  ((eq? curl-status 'fail)
				   (gtk-widget-show fail-label)
				   #f ;; don't continue
				   )
				  ((procedure? curl-status)
				   (curl-status) 
				   (gtk-widget-destroy window)
				   #f ;; we are all done
				   )
				  (else
				   #t ;; normal continue, downloading file(s) 
				      ;; (and computing sfs)
				   ))))

	      
	      (coot-thread-dispatcher
	       pdb-file-name
	       (lambda ()

		 ;; Get the PDB file if we don't have it already.
		 ;; 
		 (if (not (file-exists? pdb-file-name))
		     (begin
		       (set! curl-status 'downloading-pdb)
		       (let ((status (coot-get-url-and-activate-curl-hook pdb-url pdb-file-name 1)))
			 (if (not (= status 0))
			     
			     ;; OK failure, turn off the timeout function
			     (set! curl-status 'fail)

			     (if (not (eq? include-get-sfs-flag 'include-sfs))

				 ;; an NMR structure
				 ;; 
				 (begin
				   (set! curl-status 
					 (lambda () (read-pdb pdb-file-name)))))))))

		 
		 (if (eq? include-get-sfs-flag 'include-sfs)

		     ;; An X-ray structure
		     ;; 
		     ;; now read the sfs cif and if that is good then
		     ;; convert to mtz and run refmac, but for now, let's
		     ;; just show the PDB file.  But we can't do that in a
		     ;; subthread.  So we semaphore to the waiting timeout
		     ;; function that we are ready.  We do that by setting
		     ;; curl-status to the function that we want to
		     ;; timeout function to run when it's finishing up (on
		     ;; successful termination).
		     ;; 
		     (if (or (not (file-exists? sfs-cif-file-name))
			     (not (file-exists? refmac-out-mtz-file-name)))
			 (get-sfs-run-refmac sfs-cif-url
					     sfs-cif-file-name 
					     sfs-mtz-file-name 
					     pdb-file-name 
					     refmac-out-mtz-file-name)

			 ;; OK, the files exist already.
			 ;; 
			 (begin
			   (format #t "%%%%%%%%% files exist ~s ~s~%"
				   pdb-file-name refmac-out-mtz-file-name)
			   (set! curl-status
				 (lambda () 
				   (read-pdb pdb-file-name)
				   (make-and-draw-map-local refmac-out-mtz-file-name)))))))))))))
			 

  (define n-atoms-limit-small-ligands 6)
		
  ;; truncate (if needed) and newlineify string
  ;; 
  (define (pad-and-truncate-name name)
    (let ((sl (string-length name)))
      (if (< sl 60)
	  (string-append "\n     " name)
	  (string-append "\n     " (substring name 0 60) "..."))))

  (define (truncate-name name)
    (let ((sl (string-length name)))
      (if (< sl 60)
	  name
	  (substring name 0 60))))

  ;; 
  ;; 
  (define (make-ligands-string ligand-ht-list)

    (if (null? ligand-ht-list)
	""
	(let ((ligand-string-list
	       (map (lambda(lht)
		      (let ((name (hash-ref lht "Name"))
			    (n-atoms (hash-ref lht "NumberOfAtoms")))

			(if (not (string? name))
			    (begin
			      "")
			    (if (not (number? n-atoms))
				(begin
				  (pad-and-truncate-name name))

				(begin
				  (if (> n-atoms n-atoms-limit-small-ligands)
				      (begin
					(pad-and-truncate-name name))
				      (begin
					""))))))) ;; too small to be interesting
			   ligand-ht-list)))
	  
	  (let ((ligand-string (apply string-append ligand-string-list)))
	    (if (> (string-length ligand-string) 0)
		(string-append "\nLigands:" ligand-string)
		"")))))

  ;; return a list (possibly empty) of the three-letter-codes of the
  ;; ligands in the entry.  The tlc are strings -  (if some
  ;; problem ocurred on dehashing or it was a small/uninteresting
  ;; ligand), that gives a #f, which is filtered out before returning
  ;; 
  (define (make-ligand-tlc-list ligand-ht-list)

    (if (null? ligand-ht-list)
	'()
	(filter string? 
		(map (lambda (lht)
		       (let ((code (hash-ref lht "Code"))
			     (n-atoms (hash-ref lht "NumberOfAtoms")))
			 (if (not (string? code))
			     #f
			     (if (not (number? n-atoms))
				 code
				 (if (> n-atoms n-atoms-limit-small-ligands)
				     code
				     #f)))))
		     ligand-ht-list))))

  ;; The idea here is that we have a list of ligand tlcs.  From that,
  ;; we make a function, which, when passed a button-hbox (you can put
  ;; things into a button) will download images from pdbe and add them
  ;; to the button.  We download the images in the background
  ;; (sub-thread) and put them in coot-pdbe-ligand-images directory.
  ;; When/if the tranfer of the file is completed and good, set the
  ;; status to a function.
  ;;
  ;; a timeout function watches the value of the status and when it is
  ;; a function, it runs it.  (The function puts the image in the
  ;; button-hbox.)
  ;; 
  ;; Also we construct the url of a ribbon picture give the entry id
  ;; 
  (define (make-active-ligand-button-func entry-id ligand-tlc-list)

    (lambda (button-hbox ligands-hbox protein-ribbon-hbox) ;; returning a function

      (define (cached-or-net-get-image-func image-url image-name hbox)

	(let ((curl-status 'start)
	      (pack-image-func
	       (lambda ()
		 (let ((pixmap (gtk-pixmap-new-from-file image-name button-hbox)))
		   (gtk-box-pack-start hbox pixmap #f #f 1)
		   (gtk-widget-show pixmap)))))

	  (cache-or-net-get-image image-url image-name pack-image-func)))
	  

      ;; main line of make-active-ligand-button-func
      ;; 
      (for-each (lambda (tlc)
		  (let* ((image-size 100)
			 (image-url (string-append
				     "http://www.ebi.ac.uk/pdbe-srv/pdbechem/image/showNew?code=" 
				     tlc 
				     "&size=" (number->string image-size)))
			 (image-name (append-dir-file 
				      *coot-pdbe-image-cache-dir*
				      (string-append tlc "-" (number->string image-size) ".gif"))))

		    (cached-or-net-get-image-func image-url image-name ligands-hbox)))
		ligand-tlc-list)

      ;; now do the protein icon:
      (let* ((image-size 120)
	     (image-name-stub (string-append (string-append entry-id "_cbc" (number->string image-size) ".png")))
	     (image-url (string-append "http://www.ebi.ac.uk/pdbe-srv/view/images/entry/" image-name-stub))
	     (entry-image-file-name (append-dir-file *coot-pdbe-image-cache-dir* image-name-stub)))

	(cached-or-net-get-image-func image-url entry-image-file-name protein-ribbon-hbox))))

				   
	

  ;; hm is a hash map with keys: "Title", "Resolution", "EntryID" etc.
  ;; 
  (define (handle-pdb-entry-entity hm)

    ;; return a string.  Return "" on failure
    (define (make-authors-string citation-list)
      (if (list? citation-list)
	  (if (null? citation-list)
	      "No Authors listed"
	      (let* ((citation-hm (car citation-list)))
		(if (not (hash-table? citation-hm))
		    "No Citation table"
		    (let ((citation (hash-ref citation-hm "Authors")))
		      (string-append-with-spaces citation)))))))

    (make-directory-maybe *coot-pdbe-image-cache-dir*)

    ;; now make a button list (a label and what to do)
    ;; 
    (let* ((entry-id (hash-ref hm "EntryID"))
	   (entry-label (if (string? entry-id)
			    entry-id
			    "missing-entry-id"))) ;; should never happen
      (let ((method-item (hash-ref hm "Method"))
	    (resolution-item (hash-ref hm "Resolution"))
	    (title-item (hash-ref hm "Title"))
	    (authors-string (make-authors-string (hash-ref hm "Citation"))))
	(let* ((method-label 
		(if (not (string? method-item))
		    "Unknown method" ;; should never happen
		    method-item))
	       (title-label
		(if (not (string? title-item)) "" (truncate-name title-item)))
	       (authors-label (if (= (string-length authors-string) 0) "" (pad-and-truncate-name authors-string)))
	       (ligands (hash-ref hm "ContainsLigands"))
	       (ligands-string (make-ligands-string ligands))
	       (ligand-tlc-list (make-ligand-tlc-list ligands))
	       (resolution-string
		(if (not (string? resolution-item)) ;; (yes, really.)
		    ""
		    (string-append "     Resolution: " resolution-item)))
	       (label (string-append title-label "\n" entry-id ": " method-label resolution-string authors-label ligands-string)))

	  (if (string=? method-label "x-ray diffraction")
	      (list label
		    (lambda () 
		      (pdbe-get-pdb-and-sfs-cif 'include-sfs entry-id method-label))
		    (make-active-ligand-button-func entry-id ligand-tlc-list))
		    
	      ;; I am not interested in the ligands in NMR structures.  (Not yet).
	      (list label
		    (lambda () 
		      (pdbe-get-pdb-and-sfs-cif 'no-sfs      entry-id method-label))
		    (make-active-ligand-button-func entry-id ligand-tlc-list)))))))
		    
    
       
  ;; return a list of button
  ;; 
  (define (handle-pdb-entry-entities hash-table-list)
    (map handle-pdb-entry-entity hash-table-list))


  ;; This block is a hideous result of guess and fiddle - I don't know
  ;; anything about json, I just messed about until I got a something
  ;; interpretable.
  ;; 
  (if (not (hash-table? t))

      (begin
	(if (not (list? t))
	    t
	    (begin
	      (if (null? t)
		  '()
		  (if (null? (cdr t))
		      (map (lambda (x) (dig-table x (+ level 1))) t)
		      (begin

			;; not null? t, 
			;; not null? (cdr t)
			;; 
			(if (and (= level 3)
				 (contains-pdb-entry-entities? t))
			    (let ((button-list (handle-pdb-entry-entities t)))
			      (dialog-box-of-buttons-with-async-ligands "Recent Entries"
						     (cons 700 500)
						     button-list ;; cleverness here.  The third part of a 
						                 ;; button is a function (to add the 
                                                                 ;; ligands icon to the button).  The button
						                 ;; is passed as an argument.
						     "Close"))

			    (cons
			     (if (list? (car t))
				 (map (lambda(x) (dig-table x (+ level 1))) (car t))
				 (dig-table (car t) (+ level 1)))
			     (map (lambda(x) (dig-table x (+ level 1))) (cdr t))))))))))

      (begin 
	(let ((l (hash-map->list (lambda(x y) (list x y)) t)))
	  (if (null? (cdr l))
	      (map (lambda(x) (dig-table x (+ level 1))) (car l))
	      (cons (map (lambda(x) (dig-table x (+ level 1))) (car l))
		    (map (lambda(x) (dig-table x (+ level 1))) (cdr l))))))))


(define (recent-entries-progress-dialog)

  (let* ((window (gtk-window-new 'toplevel))
	 (dialog-name "Getting recent entries list...")
	 (main-vbox (gtk-vbox-new #f 4))
	 (label (gtk-label-new "Getting recent entries list..."))
	 (progress-bar (gtk-progress-bar-new)))

    (gtk-window-set-title window dialog-name)
    (gtk-box-pack-start main-vbox label #f #f 4)
    (gtk-box-pack-start main-vbox progress-bar #f #f 4)
    (gtk-container-add window main-vbox)
    (gtk-container-border-width window 4)
    (gtk-progress-set-show-text progress-bar #t)
    (gtk-widget-show-all window)

    (list progress-bar window)))



;; Sameer Velankar says to get this file for the latest releases
;; "http://www.ebi.ac.uk/pdbe-apps/jsonizer/latest/released/" (note the end "/").
;; 
(define (pdbe-latest-releases-gui)

  (let ((url "http://www.ebi.ac.uk/pdbe-apps/jsonizer/latest/released/")
	(json-file-name "latest-releases.json")
	(curl-status 'start))

    (add-status-bar-text "Retrieving list of latest releases...")
    (let ((progress-bars (recent-entries-progress-dialog)))

      (call-with-new-thread
       (lambda () 

	 (let ((status (coot-get-url-and-activate-curl-hook url json-file-name 1)))
	   (if (not (= status 0)) ;; i.e. not good
	       (set! curl-status 'fail-recent-entry-list)
	       (set! curl-status 
		     (lambda ()
		       (dig-table (get-recent-json json-file-name) 0)))))))

      (gtk-timeout-add 200 
		       (lambda ()
			 
			 (cond
			  ((eq? curl-status 'fail-recent-entry-list)
			   (info-dialog "Failed to get recent entry list from server.")
			   (gtk-widget-destroy (cadr progress-bars)) ;; the window
			   #f) ; done
			  ((procedure? curl-status) 
			   (gtk-widget-destroy (cadr progress-bars)) ;; the window
			   (curl-status) 
			   #f) ; all done
			  (else 
			   
			   (let ((pi (curl-progress-info json-file-name)))
			     (update-download-progress-bar-maybe pi (car progress-bars)))
			   
			   #t))))))) ;; keep going



;; put this on a menu item
; 
(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "PDB")))
      (add-simple-coot-menu-menuitem 
       menu
       "Latest Releases" pdbe-latest-releases-gui)))

