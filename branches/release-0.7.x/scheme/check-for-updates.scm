
(define my-format format)
(define my-format (lambda (args) #f)) ;; do nothing, because Coot will crash if user has the Scheme
				      ;; scripting window open (it captures the format output).

(defmacro with-working-directory args
  `(begin
     ;;     (format #t "args is ~s~%" ,args)
     (let ((current-dir (getcwd))
	   (funcs (cdr ,args)))
       ;; (format #t "funcs is ~s~%" funcs)
       (chdir (car ,args))
       @funcs
       (chdir current-dir))))


(defmacro test-macro args
  `(begin
     ;; (format #t "args is ~s~%" ',args)
     (let ((current-dir (getcwd))
	   (funcs (car (cdr ',args))))
       ;; (format #t "funcs is ~s~%" funcs)
       ,funcs
       )))




;; Is this true?
(define (get-stable-release-from-server-string str)
  str)

;; Needs testing.
(define (get-stable-release-from-coot-version)
  (let* ((s (coot-version))
	 (ls (separate-fields-discarding-char #\space s list)))
    (car ls)))


; Return #t or #f 
(define (new-version-on-server? str is-pre-release?)

  ;; (format #t "new-version-on-server? gets str: ~s~%" str)

  (if is-pre-release?

      ;; pre-releases are handle with svn revision numbers (as ints).
      ;; 
      (let ((server-rev (get-revision-from-string str)))
	;; for pre-release
	(if (number? server-rev)
	    (> server-rev (svn-revision))
	    #f))

      ;; For a stable release: 
      ;; 
      ;; The Coot release number (e.g. "0.6.4") is a string that can
      ;; be compared lexographically.
      ;; 
      (let ((server-version (get-stable-release-from-server-string str))
	    (this-build-version (get-stable-release-from-coot-version)))
	(if (string? server-version)
	    (string>? server-version this-build-version)))))



;; show the dialog
;; 
(define (notify-of-new-version str)
    (download-binary-dialog (coot-split-version-string str)))

;; version-string is something like: "coot-0.6-pre-1-revision-2060"
;; 
(define download-binary-dialog 
  (let ((pending-install-in-place #f)
	(file-name-for-progress-bar #f))
    (lambda (version-string)

      ;; 
      (define update-progress-bar
	(let ((count 0)
	      (active-count 0))
	(lambda (progress-bar)
	  (set! count (+ count 1))
	  (if (string? file-name-for-progress-bar)
	      (let ((curl-info (curl-progress-info file-name-for-progress-bar)))
		(set! active-count (+ active-count 1))
		;; (format #t "got curl-info: ~s for ~s ~%" curl-info file-name-for-progress-bar)
		(if curl-info
		    (let ((v1 (assoc 'content-length-download curl-info))
			  (v2 (assoc 'size-download           curl-info)))
		      ;; (format #t "v1: ~s v2: ~s~%" v1 v2)
		      (if (pair? v1)
			  (if (pair? v2)
			      (let ((f (/ (cdr v2) (cdr v1))))
				;; (format #t "count ~s, active-count ~s, f: ~s~%" count active-count f)
				(gtk-progress-bar-update progress-bar f)))))))))))

      ;; 
      (define (set-progress-bar-full progress-bar)
	(gtk-progress-bar-update progress-bar 1))
	

      (define (set-file-name-func file-name)
	(set! file-name-for-progress-bar file-name))

      (define (pending-install-in-place-func val)
	(set! pending-install-in-place val))

      (define (do-download-dialog)
	;; (format #t "running download-binary-dialog with version-string arg: ~s~%" version-string)
	;; 
	(let ((s (string-append "   New revision available " 
				"for this binary type:   \n"
				(coot-sys-build-type)
				"\n"
				"\n"
				version-string))
	      (revision (get-revision-from-string version-string)))

	  (let ((window (gtk-window-new 'toplevel))
		(dialog-name "Download binary")
		(main-vbox (gtk-vbox-new #f 6))
		(cancel-button (gtk-button-new-with-label "  Cancel  "))
		(ok-button (gtk-button-new-with-label "  Download and Pre-install "))
		(buttons-hbox (gtk-hbox-new #f 6))
		(progress-bar (gtk-progress-bar-new))
		(h-sep (gtk-hseparator-new))
		(info-string (gtk-label-new s)))

	    (gtk-window-set-title window dialog-name)
	    (gtk-box-pack-start buttons-hbox ok-button #t #f 6)
	    (gtk-box-pack-start buttons-hbox cancel-button #t #f 6)

	    (gtk-box-pack-start main-vbox info-string  #t #f 6) ;; not x padding, it is y padding
	    (gtk-box-pack-start main-vbox h-sep        #t #f 6)
	    (gtk-box-pack-start main-vbox buttons-hbox #t #f 6)
	    (gtk-box-pack-start main-vbox progress-bar #t #t 6)
	    (gtk-container-border-width main-vbox 6)
	    
	    (gtk-container-add window main-vbox)
	    (gtk-container-border-width window 6)

	    (gtk-signal-connect cancel-button "clicked"
				(lambda ()
				  ;; stop the download process(es) if
				  ;; they were running (they get set to
				  ;; #f when they are not running)
				  (set! pending-install-in-place 'cancelled) ;; not fail
				  (if (string? file-name-for-progress-bar)
				      (begin
					(stop-curl-download file-name-for-progress-bar))
				      (begin 
					(format #t "This shouldn't happen! (but it does)\n")))
				  (gtk-widget-destroy window)))

	    (gtk-signal-connect ok-button "clicked"
				(lambda ()
				  (if (not (number? revision))
				      (info-dialog "Failed to communicate with server")
				      
				      (begin
					(gtk-progress-set-show-text progress-bar #t)
					(call-with-new-thread
					 (lambda ()
					   (if (not (run-download-binary-curl revision version-string
									      pending-install-in-place-func
									      set-file-name-func))
					       (begin
						 (if (not (eq? pending-install-in-place 'cancelled))
						     (set! pending-install-in-place 'fail)))))
					 coot-updates-error-handler)
					
					;; and a timeout checking the progress of the download:
					(gtk-timeout-add 
					 500  (lambda ()
						 (cond
						  ((eq? pending-install-in-place 'fail)
						   (gtk-widget-destroy window)
						   (info-dialog "Failure to download and install binary")
						   #f)
						  ((eq? pending-install-in-place 'cancelled)
						   #f) ; do nothing and stop timeout
						  ((eq? pending-install-in-place 'full) ; don't go yet
						   (set-progress-bar-full progress-bar)
						   #t)
						  (pending-install-in-place
						   (gtk-widget-destroy window)
						   (restart-dialog)
						   #f)
						  (else
						   (update-progress-bar progress-bar)
						   #t))))))))
	    (gtk-widget-show-all window))))


      ;; main line
      ;; 
      (let ((coot-prefix (getenv "COOT_PREFIX")))

	(set! pending-install-in-place #f) ;; reset pending-install-in-place 
                                           ;; in case we do this a second (etc) time
	(if (not (directory-is-modifiable? coot-prefix))
	    (begin
	      (if (not (string? coot-prefix))
		  (info-dialog "COOT_PREFIX is not set.  Download not started.")
		  (info-dialog (string-append
				"Directory " coot-prefix " is not modifiable.\n"
				"Download/install not started."))))

	    (do-download-dialog))))))



		    
	
;; Test for prefix-dir 1) existing 2) being a directory 3) modifiable by user (ie. u+rwx)
;; 
;; Return #t or #f.
;; 
(define (directory-is-modifiable? prefix-dir)
  (if (not (file-exists? prefix-dir))
      #f 
      (let* ((s (stat prefix-dir))
	     (t (stat:type s)))
	(if (not (eq? t 'directory))
	    #f ; not a directory
	    (let ((p (stat:perms s)))
	      ;; test the file permissions for rwx for user using bitwise logical 
	      ;; operator on p (permissions). 448 is 256 + 128 + 64
	      (let ((b #b111000000))
		(= b (logand p b))))))))
  


(define (restart-dialog)

  (let ((window (gtk-window-new 'toplevel))
	(dialog-name "Restart Required to complete installation")
	(main-vbox (gtk-vbox-new #f 6))
	(cancel-button (gtk-button-new-with-label "  Later  "))
	(ok-button (gtk-button-new-with-label "  Restart Now "))
	(label (gtk-label-new "  Restart required to complete install "))
	(buttons-hbox (gtk-hbox-new #f 6))
	(h-sep (gtk-hseparator-new)))
    
    (gtk-window-set-title window dialog-name)
    (gtk-container-add window main-vbox)
    (gtk-container-border-width window 6)
    (gtk-box-pack-start main-vbox label #f #f 6)
    (gtk-box-pack-start main-vbox h-sep #f #f 5)
    (gtk-box-pack-start main-vbox buttons-hbox #f #f 6)
    (gtk-box-pack-start buttons-hbox ok-button #f #f 6)
    (gtk-box-pack-start buttons-hbox cancel-button #f #f 6)

    (gtk-signal-connect cancel-button "clicked"
			(lambda ()
			  (gtk-widget-destroy window)))

    (gtk-signal-connect ok-button "clicked"
			(lambda ()
			  (let ((coot-command "coot"))
			    ;; create a coot background subprocess
			    (run-concurrently coot-command)
			    (gtk-widget-destroy window)
			    (coot-real-exit 0))))
    (gtk-widget-show-all window)))
    




; 			(if success
; 			    (restart-dialog))))))))))
				  
	      
;; http://www.biop.ox.ac.uk/coot/software/binaries/pre-releases/coot-0.6-pre-1-revision-2535-binary-Linux-i386-centos-4-gtk2.tar.gz



(define check-for-updates-gui
  (let ((server-info-status #f))

    (define (handle-latest-version-server-response txt-from-server)
      ;; OK, so the server said something.
      ;; Set the status here, so that the
      ;; function that looks to see whether
      ;; or not the server responded is
      ;; notified.
      ;; 
      (if (string-match "The requested URL was not found on this server" txt-from-server)
	  (set! server-info-status 'file-not-found)
	  (set! server-info-status txt-from-server)))

    (define (get-server-info-status-thread)
      (call-with-new-thread
       (lambda ()
	 (let* ((url (make-latest-version-url))
		(latest-version-server-response (coot-get-url-as-string url)))
	   (handle-latest-version-server-response latest-version-server-response)))
       ;; the error handler
       coot-updates-error-handler))


    ;; main line
    ;; 
    (lambda ()

      (get-server-info-status-thread)
      (let ((is-pre-release? (pre-release?)))
	(let ((count 0))
	  (gtk-idle-add
	   (lambda ()
	     (cond
	      ((> count 2000) ;; try for 20 seconds, otherwise timeout.
	       ;; fail-with-timeout
	       ;; (format #t "final fail: server-info-status: ~s~%" server-info-status)
	       #f) ;; stop running this idle function
	      ((eq? server-info-status 'file-not-found)
	       (let ((s (string-append "No " 
				       (if is-pre-release? "pre-release" "release")
				       " binary for this system ("
				       (coot-sys-build-type) 
				       ") on the binary server")))
		 (info-dialog s)
		 #f)) ;; stop running this idle function
	      ((string? server-info-status)
	       (if (= (string-length server-info-status) 0)

		   ;; this happens when coot can't get past the proxy
		   ;; 
		   (info-dialog "Can't communicate with server")
		   
		   ;; otherwise we can get past the proxy, and it gave
		   ;; us a result, was the server version newer than
		   ;; this one?
		   ;; 
		   (let ((v (new-version-on-server? server-info-status is-pre-release?)))
		     ;; (format #t "new-version-on-server? returns ~s~%" v)
		     (if v 
			 (notify-of-new-version server-info-status)
			 (let ((s (string-append "No version newer than this revision ("
						 (number->string (svn-revision))
						 ").")))
			   (info-dialog s)))))
	       #f) ;; stop running idle function
	      (else 
	       (usleep 10000) 
	       ;; (format #t "server-info-status: ~s~%" server-info-status)
	       (set! count (+ count 1))
	       #t)))))))))


