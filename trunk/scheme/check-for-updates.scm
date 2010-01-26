
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




;; return an absolute file-name for file-name or #f
;;
(define (absolutify file-name)

  (if (not (string? file-name))
      #f
      (if (not (> (string-length file-name) 0))
	  "/"
	  (let ((first-char (substring file-name  0 1)))
	    (if (string=? first-char "/")
		file-name
		(append-dir-file (getcwd) file-name))))))


;; a thread handling function
(define (coot-updates-error-handler key . args)
  (my-format #t "error: finding updates: error in ~s with args ~s~%" key args))

;; The file name of the phone home script, executed with guile -s
;; FIXME - need to be installed in xxx/etc (or some such).
;; (define phone-home-cmd
;;   (string-append (getenv "HOME")
;; 		 "/Projects/coot/update-binary/phone-home.scm"))


(define (get-revision-from-string str)
  ;; e.g. str is "coot-0.6-pre-1-revision-2060" (with a newline at the
  ;; end too).  We want to return 2060 (a number) from here (or #f).
  (if (not (string? str))
      #f
      (if (= (string-length str) 0)
	  #f
	  (let* ((s (sans-final-newline str))
		 (ls (separate-fields-discarding-char #\- s list)))
	    (string->number (car (reverse ls)))))))

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


;; first generate a version string with no trailing newline.
;; 
(define (notify-of-new-version str)
  (let* ((ls (split-before-char #\c str list))
	 (ls-2 (split-before-char #\" (car (reverse ls)) list))
	 (ls-3 (split-before-char #\newline (car ls-2) list)))

;     (format #t "notify-of-new-version str: ~s~%" str)
;     (format #t "ls: ~s~%" ls)
;     (format #t "ls-2: ~s~%" ls-2)
;     (format #t "ls-3: ~s~%" ls-3)
    (download-binary-dialog (car ls-3))))

;; version-string is something like: "coot-0.6-pre-1-revision-2060"
(define download-binary-dialog 
  (let ((pending-install-in-place #f)
	(file-name-for-progress-bar #f))
    (lambda (version-string)

      ;; Get the binary (i.e. the action that happens when the download
      ;; button is pressed).  This is run in a thread, so it can't do
      ;; any graphics stuff. 
      ;; 
      ;; return #t if tar file was successfully downloaded and untared
      ;; and #f if not.
      ;; 
      (define (run-download-binary-curl revision version-string)
	(my-format "INFO:: run-download-binary-curl.... with revision ~s with version-string ~s~%" 
		   revision version-string)
	(let ((prefix (getenv "COOT_PREFIX")))
	  (if (not (string? prefix))
	      (begin
		(my-format "OOps! Can't find COOT_PREFIX~%")
		#f)
	      (let* ((pre-release-flag (string-match "-pre" (coot-version)))
		     (ys "www.ysbl.york.ac.uk/~emsley/software/binaries/")
		     (binary-type (coot-sys-build-type))
		     (host-dir (cond 
				((string=? binary-type "Linux-i386-fedora-3") ys)
				((string=? binary-type "Linux-i386-fedora-3-python") ys)
				((string=? binary-type "Linux-i386-fedora-8-python-gtk2") ys)
				((string=? binary-type "Linux-i386-fedora-8-gtk2") ys)
				((string=? binary-type "Linux-i386-fedora-10-python-gtk2") ys)
				((string=? binary-type "Linux-i686-ubuntu-8.04.3") ys)
				((string=? binary-type "Linux-i686-ubuntu-8.04.3-python") ys)
				(else 
				 "www.biop.ox.ac.uk/coot/software/binaries/")))
		     (tar-file-name (string-append version-string "-binary-" binary-type ".tar.gz"))
		     (url 
		      (if (string-match "ysbl.york.ac.uk" host-dir)
			  (if pre-release-flag 
			      (string-append "http://" host-dir "nightlies/pre-release/" tar-file-name)
			      (string-append "http://" host-dir "stable/"     tar-file-name))
			  (if pre-release-flag 
			      (string-append "http://" host-dir "pre-releases/" tar-file-name)
			      (string-append "http://" host-dir "releases/"     tar-file-name))))
		     (md5-url (string-append url ".md5sum"))
		     (md5-tar-file-name (string-append tar-file-name ".md5sum")))

		(set! file-name-for-progress-bar tar-file-name)
		(my-format "md5sum url for curl: ~s~%" md5-url)
		(my-format "url for curl: ~s~%" url)

		(coot-get-url-and-activate-curl-hook md5-url md5-tar-file-name 1)
		(coot-get-url-and-activate-curl-hook url tar-file-name 1)

		(if (not (file-exists? tar-file-name))
		    (begin
		      ;; (format #t "Ooops: ~s does not exist after attempted download~%" tar-file-name)
		      #f)
		    (if (not (file-exists? md5-tar-file-name))
			(begin
			  ;; (format #t "Ooops: ~s does not exist after attempted download~%" 
			  ;; md5-tar-file-name)
			  #f)
			(if (not (match-md5sums tar-file-name md5-tar-file-name))
			    #f 
			    (let ((success (install-coot-tar-file tar-file-name)))
			      (if success 
				  (begin
				    (set! pending-install-in-place #t)
				    #t)
				  (begin
				    ;; (format #t "Ooops: untar of ~s failed~%" tar-file-name)
				    #f))))))))))


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
		      (if (list v1)
			  (if (list v2)
			      (let ((f (/ (cdr v2) (cdr v1))))
				;; (format #t "count ~s, active-count ~s, f: ~s~%" count active-count f)
				(gtk-progress-bar-update progress-bar f)))))))))))


      ;; 
      ;; (format #t "running download-binary-dialog with version-string arg: ~s~%" version-string)
      ;; 
      (let ((s (string-append "   New revision available " 
			      "for this binary type:   \n"
			      (coot-sys-build-type)
			      "\n"
			      "\n"
			      version-string))
	    (revision (get-revision-from-string version-string)))

	;; (format #t "DEBUG:: revision ~s from version-string: ~s~%" revision version-string)

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
				(stop-curl-download file-name-for-progress-bar) 
				(set! pending-install-in-place 'cancelled) ;; not fail
				(gtk-widget-destroy window)))

	  (gtk-signal-connect ok-button "clicked"
			      (lambda ()
				(if (not (number? revision))
				    (info-dialog "Failed to communicate with server")
					 
				    (begin
				      (gtk-progress-set-show-text progress-bar #t)
				      (call-with-new-thread
				       (lambda ()
					 (if (not (run-download-binary-curl revision version-string))
					     (begin
					       (if (not (eq? pending-install-in-place 'cancelled))
						   (set! pending-install-in-place 'fail)))))
				       coot-updates-error-handler)
				      
				      ;; and a timeout checking the progress of the download:
				      (gtk-timeout-add 
				       1000  (lambda ()
					       (cond
						((eq? pending-install-in-place 'fail)
						 (gtk-widget-destroy window)
						 (info-dialog "Failure to download and install binary")
						 #f)
						((eq? pending-install-in-place 'cancelled)
						 #f) ; do nothing and stop timeout
						(pending-install-in-place
						 (gtk-widget-destroy window)
						 (restart-dialog)
						 #f)
						(else
						 (update-progress-bar progress-bar)
						 #t))))))))
				
	  (gtk-widget-show-all window))))))


;; return success status as a boolean
;;
(define (install-coot-tar-file tar-file-name)
  (let ((prefix-dir (getenv "COOT_PREFIX")))
    (if (not (string? prefix-dir))
	(begin
	  ;; (format #t "OOps could not get COOT_PREFIX~%")
	  #f)
	(if (not (directory-is-modifiable? prefix-dir))
	    (begin
	      ;; (format #t "OOps directory ~s is not modifiable~%" prefix-dir)
	      #f)
	    (let ((pending-dir (append-dir-file prefix-dir "pending-install")))
	      (if (not (file-exists? pending-dir))
		  (mkdir pending-dir))
	      (begin
		(if (not (file-exists? pending-dir))
		    (begin
		      ;; (format #t "OOps could not create ~s~%" pending-dir)
		      #f)
		    (let ((a-tar-file-name (absolutify tar-file-name)))
		      ;; with-working-directory 
		      (let ((current-dir (getcwd)))
			(chdir pending-dir)
			;; (format #t "now current-dir is ~s~%" (getcwd))
			(goosh-command "tar" (list "xzf" a-tar-file-name) '() "untar.log" #f)
			(chdir current-dir))
		      ;; (format #t "now current-dir is ~s~%" (getcwd))
		      ))))))))
		    
	
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

  

;; return as a string, or #f
(define (get-target-md5-string file-name)
  (if (not (file-exists? file-name))
      #f
      (call-with-input-file file-name
	(lambda (port)
	  (symbol->string (read port))))))

;; return a string
(define (get-md5sum-string file-name)
  (if (not (file-exists? file-name))
      #f
      (let* ((s (shell-command-to-string (string-append "md5sum " file-name)))
	     (s-bits (string->list-of-strings s)))
	(car s-bits))))


(define (match-md5sums tar-file-name target-md5sum-file-name)
  (if (not (file-exists? tar-file-name))
      #f
      (if (not (file-exists? target-md5sum-file-name))
	  (begin
	    ;; (format #t "OOps! ~s does not exist" target-md5sum-file-name)
	    #f)
	  (let ((target-md5-string (get-target-md5-string target-md5sum-file-name))
		(md5-string (get-md5sum-string tar-file-name)))
	    (if (not (string? target-md5-string))
		(begin
		  ;; (format #t "OOps ~s is not a string~%" target-md5-string)
		  #f)
		(if (not (string? md5-string))
		    (begin
		     ;; (format #t "OOps ~s is not a string~%" target-md5-string)
		      #f)
		    (if (not (string=? target-md5-string md5-string))
			(begin
			  ;;n(format #t "Oops: md5sums do not match ~s ~s.  Doing nothing~%"
			  ;; target-md5-string md5-string)
			  #f)
			#t)))))))

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

