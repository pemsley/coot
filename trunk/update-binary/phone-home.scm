

;;; We are given the binary sys type from the command line, and also
;;; whether or not this is a pre-release ("pre-release" is passed on
;;; the command line).
;;; 
;;; We need to construct the right url for the binary and go and check
;;; if it is there.

(format #t "command line args: ~s~%" (command-line))


;; Return a boolean for whether or not pre-release was passed as a
;; command line arg.
(define (extract-pre-release-flag cla)
  (let loop ((ls cla))
    (cond 
     ((null? ls) #f)
     ((string=? (car ls) "pre-release"))
     (else 
      (loop (cdr ls))))))

;; Return the binary type as a string (return #f for this part on
;; string-not-found).  The binary string is any string that is
;; proceeded with an arg of "binary"
(define (extract-binary-type cla)
  (let loop ((ls cla))
    (cond 
     ((null? ls) #f)
     ((string=? (car ls) "binary")
      (if (null? (cdr ls))
	  #f
	  (car (cdr ls)))) ;; the string we want
     (else 
      (loop (cdr ls))))))
     


(let ((pre-release-flag (extract-pre-release-flag (command-line)))
      (binary-type (extract-binary-type (command-line))))


  (format #t "pre-release: ~s~%" pre-release-flag)
  (format #t "binary-type: ~s~%" binary-type)

  (if (string? binary-type)
      
      (let ((check-url 
	     (if pre-release-flag
		 (string-append
		  "http://www.ysbl.york.ac.uk/~emsley/software/binaries/pre-release/"
		  binary-type
		  "-latest.txt")
		 (string-append
		  "http://www.ysbl.york.ac.uk/~emsley/software/binaries/release/"
		  binary-type
		  "-latest.txt"))))

	(format #t "Now get: ~s~%" check-url))))
