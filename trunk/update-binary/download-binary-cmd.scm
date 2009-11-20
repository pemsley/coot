

(use-modules (oop goops) (oop goops describe) (net http))

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




(let ((vers             (extract-version-type (command-line)))
      (binary-type      (extract-binary-type  (command-line)))
      (pre-release-flag (extract-pre-release-flag (command-line))))

  ;; debug
  (format #t "vers: ~s~%" vers) 
  (format #t "binary-type: ~s~%" binary-type) 
  (format #t "pre-release-flag: ~s~%" pre-release-flag)

  ;; kludge so that I can try to download a file that exists on the
  ;; server:
  ;; (set! binary-type "Linux-i386-fedora-10-python-gtk2")

  (if (and (string? binary-type)
	   (string? vers))

      (let ((binary-url 
	     (if pre-release-flag
		 (string-append
		  "http://www.ysbl.york.ac.uk/"
		  "~emsley/software/binaries/nightlies/pre-release/"
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

	(format #t "download-bin-cmd: get url: ~s~%" binary-url))))
	
