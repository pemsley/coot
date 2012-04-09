
(use-modules (www main)
	     (www http)
	     (www url))

;;; Allow the user to set these variables in their .mapview file if
;;; they want some server other than the default choice.  (These
;;; values cribbed from pdb-mode.el, btw (thanks Charlie Bond)).
;;; 
;;; WARNING: Right now, this will never work as intended, because
;;; there is currently no way to set the the values in ~/.mapview that
;;; we can see here. We need to have called c++ code and stored the
;;; value there and do the define in c++ Hmmm... I have no idea if
;;; that will work.
;;;
;(if (not (defined? 'oca-server))
(define oca-server "http://oca.ebi.ac.uk")

;(if (not (defined? 'oca-request-stub))
(define oca-pdb-request-stub "oca-bin/send-pdb?id=")
(define oca-sfs-request-stub "oca-bin/send-sf?r2")
(define oca-sfs-request-tail "sf.ent.Z")

(define mapview-tmp-dir ".mapview")

(define is-directory? 
  (lambda (file-name)

    (eq? (stat:type (stat file-name)) 'directory)))


; e.g. (ebi-get-pdb "1crn")
; 
; no useful return value
; 
; Note that that for sf data, we need to construct something like the
; string: http://oca.ebi.ac.uk/oca-bin/send-sf?r2acesf.ent.Z
; and we don't need to strip any html (thank goodness). 
;
; data-type can be 'pdb or 'sfs (structure factors).  We might like to use
; coordinates rather than'pdb in the future.
; 
(define get-url-str
  (lambda (url-string 'data-type)

    (define output-pdb
      (lambda (id message)
		      
	(let ((file-name (string-append mapview-tmp-dir "/" id ".pdb")))
	  (let ((s (http:message-body message)))

	    (if (not (string? s))
		
		(format #t "s is not a string~%")
		
		(begin
		  ;; The first 2 lines are HTML junk, So we should
		  ;; remove those or else mmdb doesn't read the file.
		  ;; Notice that we get those 2 lines whether the
		  ;; http:get failed or succeeded (if it failed there
		  ;; are *only* those 2 lines in the message body, so
		  ;; if we count the number of lines after that, then
		  ;; if we wrote more than 1 of them then that file
		  ;; was OK (or so my testing tells me so far).
		  (call-with-output-file file-name
		    (lambda (port)
		      (call-with-input-string s
                        (lambda (str-port)
			  
			  (read-line str-port)
			  (read-line str-port)
			  
			  (let f ((str (read-line str-port))
				  (n-lines 0))
			    
			    (if (not (eof-object? str))
				
				(begin
				  (display str port)
				  (newline port)
				  (f (read-line str-port) (+ n-lines 1)))

				(if (> n-lines 1)
				    (handle-read-draw-molecule file-name)
				    (format #t "The EBI did not have this pdb~%"))))))))))))))

    (define output-sfs 
      (lambda (id message)
	
	(let ((file-name (string-append mapview-tmp-dir "/" id ".pdb"))
	      (s (http:message-body message)))
	  (if (not (string? s))
	      (format #t "s is not a string~%")
	      (begin 
		(call-with-output-file file-name
		  (lambda (port)
		    (call-with-input-string s
                      (lambda (str-port)

			(let f ((str (read-line str-port))
				(n-lines 0))
			  (if (not (eof-object? str))
			      (begin 
				(display str port)
				(newline port)
				(f (read-line str-port) (+ n-lines 1)))
			      (if (> n-lines 1)
				  (read-cif-data file-name)
				  (format #t "The EBI did not have sfs for the code~%")))))))))))))

    ; main body
    (let ((url (url:parse url-str))
	  (nov (format #t "getting pdb message...~%"))
	   (message (http:get url)))

      (format #t " got.~%")

      (if (file-exists? mapview-tmp-dir)
	  (if (is-directory? mapview-tmp-dir)
	      
	      (begin
		(if (eq? 'data-type 'pdb)
		    (output-pdb id message))
		(if (eq? 'data-type 'sfs)
		    (output-sfs id message)))

	      (format #t "FAIL: non-directory ~s.~%" mapview-tmp-dir))

	  (begin
	    (mkdir mapview-tmp-dir)
	    (if (is-directory?? mapview-tmp-dir)

		(begin
		  (if (eq? 'data-type 'pdb)
		      (output-pdb id message))
		  (if (eq? 'data-type 'sfs)
		    (output-sfs id message)))

		(format #t "FAIL: ~s directory creation failed.~%" 
			mapview-tmp-dir))))))) ; should never happen.

;
(define get-ebi-pdb-and-sfs 
  (lambda (id)
    
    (get-ebi-pdb (id))

    (let* ((up-id (string-upcase id))
	   (url-str (string-append 
		     oca-server
		     "/"
		     oca-sfs-request-stub
		     up-id
		     oca-sfs-request-tail)))

      (get-url-str url-str 'sfs))))

(define get-ebi-pdb
  (lambda (id)

    (let* ((up-id (string-upcase id))
	   (url-str (string-append 
		     oca-server
		     "/"
		     oca-pdb-request-stub
		     up-id)))

      (get-url-str url-str 'pdb))))



; return a string
(define my-get
  (lambda ()

    (with-output-to-string 
      (lambda ()
	(call-with-input-file ".mapview/1crn.pdb"
	  (lambda (port)
	    
	    ; 2 lines of html rubbish
	    (read-line port)
	    (read-line port)

	    (let f ((s (read-line port)))
	      (if (not (eof-object? s))
		  (begin
		    (display s)
		    (newline)
		    (f (read-line port)))))))))))
			      

;(ebi-get "1crn")
