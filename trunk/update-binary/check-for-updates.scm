
;; a thread handling function
(define (coot-updates-error-handler key . args)
  (format #t "error: finding updates: error in ~s with args ~s~%" key args))

;; The file name of the phone home script, executed with guile -s
;; FIXME - need to be installed in xxx/etc (or some such).
(define phone-home-cmd
  (string-append (getenv "HOME")
		 "/Projects/coot/update-binary/phone-home.scm"))


(define (get-revision-from-string str)
  ;; e.g. str is "coot-0.6-pre-1-revision-2060" (with a newline at the
  ;; end too).  We want to return 2060 (a number) from here (or #f).
  (let* ((s (sans-final-newline str))
	 (ls (separate-fields-discarding-char #\- s list)))
    (string->number (car (reverse ls)))))

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

  (if is-pre-release?

      ;; pre-releases are handle with svn revision numbers (as ints).
      ;; 
      (let ((server-rev (get-revision-from-string str)))
	;; for pre-release
	(if (number? server-rev)
	    (> server-rev (svn-revision))))

      ;; For a stable release: 
      ;; 
      ;; The Coot release number (e.g. "0.6.4") is a string that can
      ;; be compared lexographically.
      ;; 
      (let ((server-version (get-stable-release-from-server-string str))
	    (this-build-version (get-stable-release-from-coot-version)))
	(if (string? server-version)
	    (string>? server-version this-build-version)))))


(define (notify-of-new-version str)
  (format #t "notify-of-new-version given string str: ~s~%" str)
  (let* ((ls (split-before-char #\c str list))
	 (ls2 (split-before-char #\" (car (reverse ls)) list)))

    (format #t "notify-of-new-version str: ~s~%" str)
    (format #t "ls: ~s~%" ls)
    (format #t "ls2: ~s~%" ls2)
    (download-binary-dialog (car ls2))))

;; version-string is something like: "coot-0.6-pre-1-revision-2060"
(define (download-binary-dialog version-string)

  (format #t "running download-binary-dialog with version-string arg: ~s~%" version-string)

  (let ((s (string-append "   New revision available " 
			  "for this binary type:   \n"
			  (coot-sys-build-type)
			  "\n"
			  version-string))
	(revision (get-revision-from-string version-string)))

    (let ((window (gtk-window-new 'toplevel))
	  (dialog-name "Download binary")
	  (main-vbox (gtk-vbox-new #f 6))
	  (cancel-button (gtk-button-new-with-label "  Cancel  "))
	  (ok-button (gtk-button-new-with-label "  Download  "))
	  (buttons-hbox (gtk-hbox-new #f 6))
	  (h-sep (gtk-hseparator-new))
	  (info-string (gtk-label-new s)))

      (gtk-box-pack-start buttons-hbox ok-button #t #f 6)
      (gtk-box-pack-start buttons-hbox cancel-button #t #f 6)

      (gtk-box-pack-start main-vbox info-string  #t #f 6) ;; not x padding, it is y padding
      (gtk-box-pack-start main-vbox h-sep        #t #f 6)
      (gtk-box-pack-start main-vbox buttons-hbox #t #f 6)
      (gtk-container-border-width main-vbox 6)
      
      (gtk-container-add window main-vbox)

      (gtk-signal-connect cancel-button "clicked"
			  (lambda ()
			    (gtk-widget-destroy window)))

      (gtk-signal-connect ok-button "clicked"
			  (lambda () 
			    (run-download-binary-curl revision)))

      (gtk-widget-show-all window))))
  
;; get the binary
(define (run-download-binary-curl revision)
  (format #t "::::: run-download-binary-curl....~%")
  (let* ((install-prefix "something") ;; to get curl binary
	 (pre-release-flag (string-match "-pre-" (coot-version)))
	 (host-dir "www.biop.ox.ac.uk/coot/software/binaries/")
	 (binary-type (coot-sys-build-type))
	 (url (if pre-release-flag
		  (string-append
		   "http://" 
		   host-dir "/pre-releases/"
		   "coot-"
		   version
		   "-"
		   revision
		   "-binary-"
		   binary-type
		   ".tar.gz")

		  ;; stable
		  (string-append
		   "http://" 
		   host-dir "/releases/"
		   "coot-"
		   new-version
		   "-binary-"
		   binary-type
		   ".tar.gz"))))

    (format #t "url for curl: ~s~%" url)
	 
    (goosh-command (string-append install-prefix "/bin/curl")
		   (list url)
		   '() #t "tmp-download-coot.log")))
				  
	      
  
;; http://www.biop.ox.ac.uk/coot/software/binaries/pre-releases/coot-0.6-pre-1-revision-2535-binary-Linux-i386-centos-4-gtk2.tar.gz



(let ((menu (coot-menubar-menu "Updates")))
  (add-simple-coot-menu-menuitem
   menu "Check for updates..."
   (let ((server-info-status #f))


     (define (get-server-info-status-thread)
       (call-with-new-thread
	(lambda ()
	  (format #t "get updates info thread~%")
	  ;; here we construct args to goosh-command,
	  ;; adding in "pre-release" if this binary is a
	  ;; pre-release.
	  ;; args ends up as something like:
	  ;; ("-s" "xxx/phone-home.scm" "pre-release" 
	  ;;  "binary" "Linux-1386-fedora-10-python-gtk2"
	  ;;  "command-line" "/home/xx/coot/bin/coot")
	  (let* ((update-coot-log "tmp-update-coot.log")
		 (args-1 (list "binary" (coot-sys-build-type)
			       "command-line" (car (command-line))))
		 (st (coot-version))
		 (args-2 (if (string-match "-pre-" st)
			     (cons "pre-release" args-1)
			     args-1))
		 (args-3 (append 
			  (list "-s" phone-home-cmd) args-2)))
	    (if (file-exists? update-coot-log)
		(delete-file update-coot-log))
	    (format #t "about to: guile ~s~%" args-3)
	    (goosh-command "guile" args-3 '() update-coot-log #f)
	    (format #t "done: guile ~s~%" args-3)
	    (if (file-exists? update-coot-log)
		(begin
		  ;; OK, so the server said something.
		  ;; Set the status here, so that the
		  ;; function that looks to see whether
		  ;; or not the server responded is
		  ;; notified.
		  ;; 
		  (call-with-input-file update-coot-log
		    (lambda (port)
		      (let ((line (read-line port)))
			(format #t "got line ~s from ~s~%" line update-coot-log)
			;; recall that we can't do GUI
			;; things in a sub-thread.
			(format #t "setting server-info-status to ~s~%" line)
			(set! server-info-status line))))
		  ;; x(delete-file update-coot-log)
		  ))))
                               ;;; thread ends here.
	
	coot-updates-error-handler))


     (lambda ()

       (get-server-info-status-thread)

       (let ((is-pre-release? (string-match "-pre-" (coot-version))))
	 (let ((count 0))
	   (gtk-idle-add
	    (lambda ()
	      (cond
	       ((> count 2000) ;; try for 20 seconds, otherwise timeout.
		;; fail-with-timeout
		(format #t "final fail: server-info-status: ~s~%" server-info-status)
		#f) ;; stop running this idle function
	       ((string? server-info-status) 
		(if (new-version-on-server? server-info-status is-pre-release?)
		    (notify-of-new-version server-info-status)
		    (let ((s (string-append "No version newer than this revision ("
					    (number->string (svn-revision))
					    ").")))
		      (info-dialog s)))
		#f) ;; stop running idle function
	       (else 
		(usleep 10000) 
		;; (format #t "server-info-status: ~s~%" server-info-status)
		(set! count (+ count 1))
		#t))))))))))

	  
