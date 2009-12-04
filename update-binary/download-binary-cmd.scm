
;; download-binary-cmd.scm is a stand-alone executable script run from
;; coot's guile.  
;;
;; download-binary-cmd.scm constructs a string from its args and then
;; does something like wget on that file.  It needs to access the web,
;; so we currently use net http, but in future should use libcurl
;; (more tricky).


(use-modules (oop goops) 
	     (oop goops describe)
	     (net http)
	     (os process))

;; Return #t or #f:
(define (command-in-path? cmd)

  ;; test for command (see goosh-command-with-file-input description)
  ;; 
  (if (string? cmd) 
      (let ((have-command? (run "which" cmd)))
	
	(= have-command? 0)) ; return #t or #f
      #f))


;; Where cmd is e.g. "refmac" 
;;       args is (list "HKLIN" "thing.mtz")
;;       log-file-name is "refmac.log"      
;;       data-list is (list "HEAD" "END")
;; 
;; Return the exist status e.g. 0 or 1.
;; 
(define (goosh-command cmd args data-list log-file-name screen-output-also?)
    
  (if (not (command-in-path? cmd))

      (format #t "command ~s not found~%" cmd)

      (let* ((cmd-ports (apply run-with-pipe (append (list "r+" cmd) args)))
	     (pid (car cmd-ports))
	     (output-port (car (cdr cmd-ports)))
	     (input-port  (cdr (cdr cmd-ports))))
	
	(let loop ((data-list data-list))
	  (if (null? data-list)
	      (begin 
		(close input-port))
	      
	      (begin
		(format input-port "~a~%" (car data-list))
		(loop (cdr data-list)))))
	
	(call-with-output-file log-file-name
	  (lambda (log-file-port)
	    
	    (let f ((obj (read-line output-port)))
	      (if (eof-object? obj)
		  (begin 
		    (let* ((status-info (waitpid pid))
			   (status (status:exit-val (cdr status-info))))
		      (format #t "exit status: ~s~%" status)
		      status)) ; return status 
		  
		  (begin
		    (if (eq? screen-output-also? #t)
			(format #t ":~a~%" obj))
		    (format log-file-port "~a~%" obj)
		    (f (read-line output-port))))))))))

;; Return a boolean for whether or not pre-release was passed as a
;; command line arg.
(define (extract-pre-release-flag cla)
  (let loop ((ls cla))
    (cond 
     ((null? ls) #f)
     ((string=? (car ls) "pre-release"))
     (else 
      (loop (cdr ls))))))


;; Return the argument after interesting-string, or #f if no such
;; arg.
(define (extract-next-arg-from-command-line interesting-string cla)
  (let loop ((ls cla))
    (cond
     ((null? ls) #f)
     ((string=? (car ls) interesting-string)
      (if (null? (cdr ls))
	  #f
	  (car (cdr ls))))
     (else 
      (loop (cdr ls))))))
  
;; Return the binary type as a string (return #f for this part on
;; string-not-found).  The binary string is any string that is
;; proceeded with an arg of "binary"
(define (extract-binary-type cla)
  (extract-next-arg-from-command-line "binary" cla))

;; Return the binary type as a string (return #f for this part on
;; string-not-found).  The binary string is any string that is
;; proceeded with an arg of "binary"
(define (extract-version-type cla)
  (extract-next-arg-from-command-line "version" cla))


(format #t ":::::::::: Here: ~s~%" (command-line))


(let ((vers             (extract-version-type (command-line)))
      (binary-type      (extract-binary-type  (command-line)))
      (pre-release-flag (extract-pre-release-flag (command-line))))

  ;; debug
  (format #t ":::::::::: vers: ~s~%" vers) 
  (format #t ":::::::::: binary-type: ~s~%" binary-type) 
  (format #t ":::::::::: pre-release-flag: ~s~%" pre-release-flag)

  ;; kludge so that I can try to download a file that exists on the
  ;; server:
  ;; (set! binary-type "Linux-i386-fedora-10-python-gtk2")

  (if (and (string? binary-type)
	   (string? vers))

      (let ((binary-url 
	     (if pre-release-flag
		 (string-append
		  ;;"http://www.ysbl.york.ac.uk/"
		  ;; "~emsley/software/binaries/nightlies/pre-release/"
		  "http://www.biop.ox.ac.uk/coot/software/binaries/pre-releases/"
		  vers
		  "-binary-"
		  binary-type
		  ".tar.gz")

		 (string-append
		  "http://www.ysbl.york.ac.uk/"
		  "~emsley/software/binaries/release/"
		  vers
		  "-binary-"
		  binary-type
		  ".tar.gz"))))

	(format #t "download-bin-cmd: get url: ~s~%" binary-url)
	(goosh-command "wget" (list binary-url) '() #t "download-tmp.log"))))
