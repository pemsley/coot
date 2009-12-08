
;; a thread handling function
(define (coot-updates-error-handler key . args)
  (format #t "error: finding updates: error in ~s with args ~s~%" key args))

;; The file name of the phone home script, executed with guile -s
;; FIXME - need to be installed in xxx/etc (or some such).
;; (define phone-home-cmd
;;   (string-append (getenv "HOME")
;; 		 "/Projects/coot/update-binary/phone-home.scm"))


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


;; return a string with no trailing newline.
;; 
(define (notify-of-new-version str)
  (format #t "notify-of-new-version given string str: ~s~%" str)
  (let* ((ls (split-before-char #\c str list))
	 (ls-2 (split-before-char #\" (car (reverse ls)) list))
	 (ls-3 (split-before-char #\newline (car ls-2) list)))

    (format #t "notify-of-new-version str: ~s~%" str)
    (format #t "ls: ~s~%" ls)
    (format #t "ls-2: ~s~%" ls-2)
    (format #t "ls-3: ~s~%" ls-3)
    (download-binary-dialog (car ls-3))))

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
      (gtk-container-border-width window 6)

      (gtk-signal-connect cancel-button "clicked"
			  (lambda ()
			    (gtk-widget-destroy window)))

      (gtk-signal-connect ok-button "clicked"
			  (lambda () 
			    (run-download-binary-curl revision version-string)))

      (gtk-widget-show-all window))))

;; return success status
(define (install-coot-tar-file tar-file-name)
  'xx)

;; 
(define (directory-is-modifiable? prefix-dir)
  #f)

;; return as a symbol, or #f
(define (get-target-md5-string file-name)
  (if (not (file-exists? file-name))
      #f
      (call-with-input-file file-name
	(lambda (port)
	  (read port)))))

(define (get-md5sum-string file-name)
  (if (not (file-exists? file-name))
      #f
      (let ((s (shell-command-to-string (string-append "md5sum " file-name))))
	(format #t "debug: s is ~s~%" s))))

  
	   
  

(define (match-md5sums tar-file-name target-md5sum-file-name)
  (if (not (file-exists? tar-file-name))
      #f
      (if (not (file-exists? target-md5sum-file-name))
	  #f
	  (let ((target-md5-string (get-target-md5-string target-md5sum-file-name))
		(md5-string (get-md5sum-string tar-file-name)))
	    (if (not (string? target-md5-string)
		#f
		(if (not (string? md5-string))
		    #f
		    (string=? target-md5-string md5-string))))))))

  
;; get the binary
;; 
;; 
(define (run-download-binary-curl revision version-string)
  (format #t "::::: run-download-binary-curl.... with revision ~s with version-string ~s~%" 
	  revision version-string)
  (let ((prefix (getenv "COOT_PREFIX")))
    (if (string? prefix)
	(let* ((curl-exe (string-append prefix "/bin/curl")) ;; to get curl binary
	       (pre-release-flag (string-match "-pre" (coot-version)))
	       (host-dir "www.biop.ox.ac.uk/coot/software/binaries/")
	       (binary-type (coot-sys-build-type))
	       (url (if pre-release-flag
			(string-append
			 "http://" 
			 host-dir "pre-releases/"
			 version-string
			 "-binary-"
			 binary-type
			 ".tar.gz")

			;; stable
			(string-append
			 "http://" 
			 host-dir "releases/"
			 "coot-"
			 new-version
			 "-binary-"
			 binary-type
			 ".tar.gz")))
	       (md5-url (string-append url ".md5sum")))

	  (format #t "url for curl: ~s~%" url)
	  (goosh-command curl-exe (list md5-url) '() #t "tmp-download-coot-md5sum.log")
	  (goosh-command curl-exe (list url) '() #t "tmp-download-coot.log")
	  (if (match-md5sums tar-file-name md5-tar-file-name)
	      (let ((success (install-coot-tar-file tar-file-name)))
		(if success
		    (restart-dialog))))))))
				  
	      
  
;; http://www.biop.ox.ac.uk/coot/software/binaries/pre-releases/coot-0.6-pre-1-revision-2535-binary-Linux-i386-centos-4-gtk2.tar.gz



(let ((menu (coot-menubar-menu "Updates")))
  (add-simple-coot-menu-menuitem
   menu "Check for updates..."
   (let ((server-info-status #f))

     ;; return a boolean
     (define (pre-release?)
       (string-match "-pre" (coot-version)))

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

     ;; here we construct args to goosh-command,
     ;; adding in "pre-release" if this binary is a
     ;; pre-release.
     ;; args ends up as something like:
     ;; ("-s" "xxx/phone-home.scm" "pre-release" 
     ;;  "binary" "Linux-1386-fedora-10-python-gtk2"
     ;;  "command-line" "/home/xx/coot/bin/coot")
     ;; 
     (define (make-latest-version-url)
       (let ((build-type (coot-sys-build-type)))
	 
	 ;; hack!
	 ;; 
	 (set! build-type "Linux-i386-centos-4-gtk2")

	 (string-append 
	  "http://www.biop.ox.ac.uk/coot/software/binaries/"
	  (if (pre-release?)
	      "pre-releases"
	      "releases")
	  "/"
	  "type-binary-"
	  build-type
	  "-latest.txt")))

     (define (get-server-info-status-thread)
       (call-with-new-thread
	(lambda ()
	  (let* ((url (make-latest-version-url))
		 (nov (format #t "DEBUG:: attempting to get URL ~s~%" url))
		 (latest-version-server-response (coot-get-url-as-string url)))
	    (handle-latest-version-server-response latest-version-server-response)))
	;; the error handler
	coot-updates-error-handler))


     (lambda ()

       (get-server-info-status-thread)

       (let ((is-pre-release? (pre-release?)))
	 (let ((count 0))
	   (gtk-idle-add
	    (lambda ()
	      (cond
	       ((> count 2000) ;; try for 20 seconds, otherwise timeout.
		;; fail-with-timeout
		(format #t "final fail: server-info-status: ~s~%" server-info-status)
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

	  
