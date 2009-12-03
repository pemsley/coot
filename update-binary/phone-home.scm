
(use-modules (oop goops) (oop goops describe) (net http))

;;; We are given the binary sys type from the command line, and also
;;; whether or not this is a pre-release ("pre-release" is passed on
;;; the command line).
;;; 
;;; We need to construct the right url for the binary and go and check
;;; if it is there.


;; Return a boolean for whether or not pre-release was passed as a
;; command line arg.
(define (extract-pre-release-flag cla)
  (let loop ((ls cla))
    (cond 
     ((null? ls) #f)
     ((string=? (car ls) "pre-release"))
     (else 
      (loop (cdr ls))))))


;; Return the argument after interesting-string, or null if no such
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

;; Return what we presume to be the first argument in the command line
;; for the program that called this program.  We do this so that we
;; can find the bin dir and bin/.. dir and test whether we can write
;; into it.
;; 
(define (extract-coot-command-line cla)
  (extract-next-arg-from-command-line "command-line" cla))


(define get-url-as-string
  (lambda (my-url)
    
    (call-with-output-string
      (lambda (port)

	(let ((initial-output-port (current-output-port)))

	  (set-current-output-port port)
	  (http-get my-url)
	  (set-current-output-port initial-output-port))))))


(let ((pre-release-flag (extract-pre-release-flag (command-line)))
      (binary-type (extract-binary-type (command-line))))

(set! binary-type "Linux-i386-centos-4-python-gtk2")

(format #t "pre-release: ~s~%" pre-release-flag)
(format #t "binary-type: ~s~%" binary-type)

  (if (string? binary-type)
      
      (let ((check-url 
	     (if pre-release-flag
		 (string-append
		  "http://www.biop.ox.ac.uk/"
		  "coot/software/binaries/pre-releases/"
		  "type-binary-"
		  binary-type
		  "-latest.txt")
		 (string-append
		  "http://www.biop.ox.ac.uk/"
		  "coot/software/binaries/release/"
		  "type-binary-"
		  binary-type
		  "-latest.txt"))))

	(format #t "Now get: ~s~%" check-url)
	
	(let ((s (get-url-as-string check-url)))

	  (format #t "LATEST-BUILD-ON-SERVER: ~s~%" s)))))
