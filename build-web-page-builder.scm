
(use-modules (srfi srfi-13)
             (srfi srfi-1)
             (sxml simple)
	     (ice-9 string-fun)
	     (ice-9 regex))



; (define build-prefix (add-dir-file (getenv "HOME") "autobuild"))
(define bin-tar-dir "/y/people/emsley/public_html/software/binaries/nightlies/pre-release")
(define source-tar-dir "/y/people/emsley/public_html/software/pre-release")

(define host-name (getenv "HOST"))

;; remove any trailing /s
;;
(define (directory-as-file-name dir)

  (if (= 0 (string-length dir))
      "."
      (let ((rchars (reverse (string->list dir))))
	(let f ((rls rchars))
	  (cond 
	   ((eq? #\/ (car rls)) (f (cdr rls)))
	   (else
	    (list->string (reverse rls))))))))



;;;  in a laughable attempt to minimise system dependence.
;;;
(define (append-dir-file dir-name file-name)

  (string-append (directory-as-file-name dir-name) "/" file-name))

;;; similarly attempting to minimise system dependence.
(define (append-dir-dir dir-name sub-dir-name)

  (string-append (directory-as-file-name dir-name) "/" sub-dir-name))

;; or define a utility function for this
(define (directory-files dir)

  (if (not (access? dir R_OK))
    '()
    (let ((p (opendir dir)))
      (do ((file (readdir p) (readdir p))
           (ls '()))
          ((eof-object? file) (closedir p) (reverse! ls))
        (set! ls (cons file ls))))))


;; where is the pre-release tar.gz file and how old is it?
;; binary-match-string is a string unique to that binary distribution
;; 
;; find the file with the biggest revision number in dir and return
;; its mtime
;; 
;; 
(define (latest-tar-gz dir binary-match-string)
  (let* ((files (directory-files dir))
	 (binary-match-string-tar-gz 
	  (string-append binary-match-string ".tar.gz")))

    ;; split on "-" and take the last field, convert it to a number
    (define (get-revision-number str)
      ; (format #t "==== trying to get revision number from ~s~%" str)
      (let ((bits (string-split str #\-)))
	(let loop ((bits bits))
	  (cond
	   ((null? bits) #f)
	   ((string=? "revision" (car bits))
	    (let ((n (string->number (car (cdr bits)))))
	      (if (number? n)
		  n
		  #f)))
	   (else 
	    (loop (cdr bits)))))))

		
    ; (format #t "files: ~s~%from ~s ~%" files dir)
    (let loop ((files files)
	       (latest-coot-tar (cons 0 #f)))
      (cond 
       ((null? files) 
	(if (not (cdr latest-coot-tar))
	    #f
	    (let ((mtime (stat:mtime (stat (cdr latest-coot-tar)))))
	      (list mtime 
		    (car latest-coot-tar)
		    (cdr latest-coot-tar)))))
       ((and (string-match "coot-" (car files))
	     (string-match ".tar.gz" (car files))
	     (string-match binary-match-string-tar-gz (car files))
	     (not (string-match ".md5sum" (car files))))
	 (let ((full-file (append-dir-file dir (car files))))
	   ;; what's the revision number?
	   (let ((bits (string-split (car files) #\.)))
	     ;; which bit has "revision"? use that to extract number
	     (let bits-loop ((bits bits))
	       (cond
		((null? bits) (loop (cdr files) latest-coot-tar))
		((string-match "revision" (car bits))
		 (let ((rn (get-revision-number (car bits))))
		   ; (format #t "==== revision number ~s~%" rn)
		   (if (not (number? rn))
		       (loop (cdr files) latest-coot-tar)
		       (if (< rn (car latest-coot-tar))
			   (loop (cdr files) latest-coot-tar)
			   (loop (cdr files) (cons rn full-file))))))
		(else 
		 (bits-loop (cdr bits))))))))
       (else 
	(loop (cdr files) latest-coot-tar))))))

;; time is time in seconds.
;; return a list of (sec minutes hours days months years)
;; 
(define (oldness time)

  (let* ((n-years   (floor (/ time (* 60 60 24 30 12))))
	 (n-months  (- (floor (/ time (* 60 60 24 30))) (* 12 n-years)))
	 (n-days    (- (floor (/ time (* 60 60 24))) (+ (* 30 n-months)
							(* 30 12 n-years))))
	 (n-hours   (- (floor (/ time (* 60 60))) (+ (* 24 n-days)
						     (* 24 30 n-months)
						     (* 24 30 12 n-years))))
	 (n-minutes (- (floor (/ time 60)) (+ (* 60 n-hours)
					      (* 60 24 n-days)
					      (* 60 24 30 n-months)
					      (* 60 24 30 12 n-years))))
	 (n-seconds (- time (+ (* 60 n-minutes)
			       (* 60 60 n-hours)
			       (* 60  60 24 n-days)
			       (* 60  60 24 30 n-months)
			       (* 60  60 24 30 12 n-years)))))


    (list (inexact->exact n-seconds)
	  (inexact->exact n-minutes)
	  (inexact->exact n-hours)
	  (inexact->exact n-days)
	  (inexact->exact n-months)
	  (inexact->exact n-years))))

;; return #f or the age of the file in seconds (relative to now)
;; 
(define (oldness-info file-info)
  ;; (format #t "oldness-info on ~s~%" file-info)
  (if (not file-info)
      #f
      (oldness (- (current-time) (car file-info)))))

;; 
(define (write-sxml block filename)
  
  (call-with-output-file filename
      (lambda (port)
	(sxml->xml block port))))

(define (time-text an-oldness)
  ;; 5: years
  ;; 4: months
  ;; 3: days
  ;; 2: hours
  ;; 1: minutes
  ;; 0; seconds
  (if (> (list-ref an-oldness 5) 0)
      (format #f "~s years ~s months ~s days ~s hours ~s minutes"
	      (list-ref an-oldness 5)
	      (list-ref an-oldness 4)
	      (list-ref an-oldness 3)
	      (list-ref an-oldness 2)
	      (list-ref an-oldness 1))
        (if (> (list-ref an-oldness 4) 0)
	    (format #f "~s month~a ~s days ~s hours ~s minutes"
		    (list-ref an-oldness 4)
		    (if (= (list-ref an-oldness 4) 1) "" "s")
		    (list-ref an-oldness 3)
		    (list-ref an-oldness 2)
		    (list-ref an-oldness 1))
	    (if (> (list-ref an-oldness 3) 0)
		(format #f "~s day~a ~s hour~a ~s minute~a"
			(list-ref an-oldness 3)
			(if (= (list-ref an-oldness 3) 1) "" "s")
			(list-ref an-oldness 2)
			(if (= (list-ref an-oldness 2) 1) "" "s")
			(list-ref an-oldness 1)
			(if (= (list-ref an-oldness 1) 1) "" "s"))
		(if (> (list-ref an-oldness 2) 0)
		    (format #f "~s hour~a ~s minute~a"
			    (list-ref an-oldness 2)
			    (if (= (list-ref an-oldness 2) 1) "" "s")
			    (list-ref an-oldness 1)
			    (if (= (list-ref an-oldness 1) 1) "" "s"))
		    (if (> (list-ref an-oldness 1) 0)
			(format #f "~s minute~a"
				(list-ref an-oldness 1)
				(if (= (list-ref an-oldness 1) 1) "" "s"))
			(format #f "~s seconds"
				(list-ref an-oldness 0))))))))

;; time-diff in seconds
(define (time-diff->n-spaces time-diff)
  ; (format #t "DEBUG:: time-diff: ~s~%" time-diff)
  (if time-diff
      (inexact->exact (/ time-diff (* 60 60 3)))
      0))

;; return a string of n spaces
(define (n-spaces n)

  (if (< n 1) 
      ""
      (let loop ((n n)
		 (str ""))
	(cond 
	 ((< n 1) str)
	 (else
	  (loop (- n 1) (string-append "|" str)))))))


(define (make-page bin-list file-name)

  (define (build-log-page file-info page-type)
    ;; (format #t "build-log-page on ~s~%" file-info)
    (string-append "http://www.ysbl.york.ac.uk/~emsley/build-logs/"
		   (car (cdr (list-ref file-info 0)))
		   (if (eq? page-type 'log)
		       ".log"
		       "-testing.log")))

  ;; return a time, or #f
  (define (make-time-diff file-info now-time)
    (if (= (length file-info) 2)
	#f
	(let* ((most-time 2000000)
	       (time-diff-pre (- now-time (car (list-ref file-info 3)))))
	  (if (< time-diff-pre most-time) time-diff-pre most-time))))


  (let ((latest-source-info (latest-tar-gz source-tar-dir ""))
	(binary-file-infos
	 (map
	  (lambda (bin)
	    (let ((l (latest-tar-gz bin-tar-dir (car bin))))
	      (if l 
		  (list bin (car (cdr l)) (oldness-info l) l)
		  (begin
		    (format #t "bad return from latest-tar-gz for ~s~%" bin)
		    (list bin #f)))))
	  bin-list)))

    (write-sxml 
     `(html (head (title "Coot Build Summary Page")
		  (meta (@ (http-equiv refresh) (content 600))))
	    (body 
	     (p ("Generated " ,(strftime "%a %d %b %H:%M:%S %G %Z" (localtime (current-time)))))
	     ;; source code
	     (p "Source code "
		,(basename (list-ref latest-source-info 2) ".tar.gz")
		" "
		,(list-ref latest-source-info 1)
		(br)
		,(time-text (oldness-info latest-source-info)))
	     ;; binary targets
	     ,(let ((now-time (current-time)))
		(map (lambda (file-info)
		       (format #t "file-info ~s~%" file-info)
		       (let* ((make-links? (car (cdr (car file-info))))
			      (time-diff (make-time-diff file-info now-time))
			      (ns (n-spaces (time-diff->n-spaces time-diff)))
			      (colour (cond
				       ((eq? time-diff #f) "#303030")
				       ((< time-diff (* 60 60 24)) "#80dd80")
				       (else 
					"#dd8080"))))

			 (if (= (length file-info) 2)
			     `(p 
			       ,(car (list-ref file-info 0)) '(b " not found")
			       (br)
			       (a (@ href ,(build-log-page file-info 'log)) build-log)
			       " "
			       (a (@ href ,(build-log-page file-info 'test)) test-log))
			     `(p 
			       ,(car (car file-info))
			       " "
			       ,(car (cdr file-info))
			       (br)
			       ,(time-text (list-ref file-info 2))
			       " "
			       ,(if make-links?
				    `(a (@ href ,(build-log-page file-info 'log)) build-log)
				    "")
			       " "
			       ,(if make-links?
				    `(a (@ href ,(build-log-page file-info 'test)) test-log)
				    "")
			       (table (@ (border 1) (bgcolor ,colour))
				      (tr (td ,ns)))
			       ))))
		     binary-file-infos))))
     file-name)))
      

(let ((bin-list 
       (list 
	(list "binary-Linux-i386-fedora-6" 
	      "Linux-bragg1.chem.york.ac.uk/bragg1.chem.york.ac.uk")
	(list "binary-Linux-i386-fedora-6-python"  
	      "Linux-bragg1.chem.york.ac.uk/bragg1.chem.york.ac.uk")
	(list "binary-Linux-i386-fedora-6-python-gtk2"
	      "Linux-bragg1.chem.york.ac.uk/gtk2-build")

	(list "binary-Linux-i386-fedora-5" 
	      "Linux-bragg3.chem.york.ac.uk/bragg3.chem.york.ac.uk")
	(list "binary-Linux-i386-fedora-5-python" 
	      "Linux-bragg3.chem.york.ac.uk/bragg3.chem.york.ac.uk")
	(list "binary-Linux-i386-fedora-5-python-gtk2"
	      "Linux-bragg3.chem.york.ac.uk/gtk2-build")

	(list "binary-Linux-i386-redhat-8.0" 
	      "Linux-bubbles/bubbles")
	(list "binary-Linux-i386-redhat-8.0-python" 
	      "Linux-bubbles/bubbles")

	(list "binary-Linux-i386-fedora-3-python-gtk2"
	      "Linux-bunyip.chem.york.ac.uk/gtk2-build")

	(list "binary-Linux-i386-fedora-3"
	      "Linux-bunyip.chem.york.ac.uk/bunyip.chem.york.ac.uk")

	(list "binary-Linux-i386-fedora-3-python"
	      "Linux-bunyip.chem.york.ac.uk/bunyip.chem.york.ac.uk")

	(list "binary-Linux-i386-fedora-8-python-gtk2"
	      "Linux-dragon.chem.york.ac.uk/gtk2-build")

	(list "binary-Linux-i386-fedora-8"
	      "Linux-dragon.chem.york.ac.uk/dragon.chem.york.ac.uk")

	(list "binary-Linux-i686-ubuntu-6.06.1-python-gtk2" #f)))

      (file-name "/y/people/emsley/public_html/coot/build-info.html"))
  (make-page bin-list file-name))


