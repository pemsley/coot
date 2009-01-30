
(use-modules (srfi srfi-13)
             (srfi srfi-1)
             (sxml simple)
	     (ice-9 string-fun)
	     (ice-9 popen)
	     (net http)
	     (goosh)
	     (ice-9 regex))

(use-syntax (ice-9 syncase))


; (define build-prefix (add-dir-file (getenv "HOME") "autobuild"))

(define source-tar-dir "/y/people/emsley/public_html/software/pre-release")

(define (bin-tar-dir person)
  (string-append "/y/people/"
		 person 
		 "/public_html/software/binaries"
		 "/nightlies/pre-release"))

(define host-name (getenv "HOST"))


(define-syntax my-ns
  (syntax-rules ()
    ((my-ns x n) (cons 
		  'td 
		  (let loop ((n1 (let ((max-n 80))
				   (if (< n max-n) n max-n)))
			     (x-list '()))
		    (cond 
		     ((= n1 0) x-list)
		     (else 
		      (loop (- n1 1) (cons x x-list)))))))))

;(define-syntax my-ns
;  (syntax-rules ()
;    ((my-ns x n) (cons 
;		  'td 
;		  (let loop ((n1 n))
;		    (cond 
;		     ((= n1 0) '())
;		     (else 
;		      (cons x (loop (- n1 1))))))))))


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


;;; range: works like the eponymous python function
;;; e.g. (range 3) -> '(0, 1, 2)
;;; e.g. (range 1 3) -> '(1, 2)
;;;
(define (range first . second)

  (if (number? first)
		
      ;; OK input
      (let ((n-low/n-high (if (= (length second) 0)
			      (cons 0 first)
			      (cons first (car second)))))
	
	(let loop ((count (car n-low/n-high))
		   (rng '()))
	  (cond 
	   ((>= count (cdr n-low/n-high)) (reverse rng))
	   (else 
	    (loop (+ count 1)
		  (cons count rng))))))))

;; where is the pre-release tar.gz file and how old is it?
;; binary-match-string is a string unique to that binary distribution
;; 
;; find the file with the biggest revision number in dir and return
;; its mtime
;; 
;; 
(define (latest-tar-gz dir binary-match-string suffix-string)
  (let* ((files (directory-files dir))
	 (binary-match-string-tar-gz (string-append binary-match-string suffix-string)))

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
    
    (format #t "binary-match-string: ~s  suffix: ~s~%" 
	    binary-match-string suffix-string)

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

		 
       ((or (and (string-match "coot-" (car files))
		 (string-match (string-append binary-match-string 
					      suffix-string)
			       (car files))
		 (not (string-match ".md5sum" (car files))))
	    
	    (and (string-match "WinCoot" (car files))
		 (not (string-match ".md5sum" (car files)))))
		 
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
(define (time-diff->n time-diff)
  ; (format #t "DEBUG:: time-diff: ~s~%" time-diff)
  (if time-diff
      (+ (inexact->exact (/ time-diff (* 60 60 1))) 1)
      0))


;; code from thi <ttn at mingle.glug.org>
;; 
;; run command and put the output into a string and return
;; it. (c.f. @code{run-command/strings})
(define (shell-command-to-string cmd)
  (with-output-to-string
    (lambda ()
      (let ((in-port (open-input-pipe cmd)))
	(let loop ((line (read-line in-port 'concat)))
	  (or (eof-object? line)
	      (begin
		(display line)
		(loop (read-line in-port 'concat)))))))))

;; return the revision number or #f
;; 
(define (get-svn-revision)
  (let ((coot-dir (string-append (getenv "HOME") "/Projects/coot"))
	(current-dir (getcwd)))
    ; (format #t "coot-dir: ~s~%" coot-dir)
    ; (format #t "curr-dir: ~s~%" current-dir)
    (chdir coot-dir)
    (let ((s (shell-command-to-string "svn -qu status")))
      (chdir current-dir)
      (string->number
       (car (reverse (string-split (car (cdr (reverse (string-split s #\newline)))) #\space)))))))


;; return a string of n spaces
(define (n-bars n)

  (if (= n 0) 
      "."
      (let loop ((n n)
		 (str ""))
	(cond 
	     ((< n 1) str)
	     (else
	      (loop (- n 1) (string-append "|" str)))))))
  

(define (make-page bin-list file-name)

  ;; link could be "Linux-host/gtk2" or '('absolute "http://x.ac.uk/build").
  ;; page-type is 'build, 'build-dir or 'test
  ;; 
  (define (build-log-page file-info page-type pythonic?)
    (let ((link (car (cdr (list-ref file-info 0)))))
      (format #t "build-log-page on file-info: ~s~%" file-info)
      (format #t "build-log-page on link: ~s~%" file-info)

      (if (eq? page-type 'build-dir)

	  (dirname
	   (if (list? link)
	       (car (cdr link))
	       (string-append "http://www.ysbl.york.ac.uk/~emsley/build-logs/"
			      link)))
	  
	  (string-append 
	   (if (list? link)
	       (car (cdr link))
	       (string-append "http://www.ysbl.york.ac.uk/~emsley/build-logs/"
			      link))
	   (cond
	    ((eq? page-type 'build) (string-append 
				     (if pythonic? 
					 "-python"
					 "")
				     "-build.log"))
	    ((eq? page-type 'test) "-test.log")
	    ((eq? page-type 'test-python) "-test-python.log")
	    (else 
	     "strange"))))))


  ;; return a time, or #f
  (define (make-time-diff file-info now-time)
    (if (= (length file-info) 2)
	#f
	(let* ((most-time 2000000)
	       (time-diff-pre (- now-time (car (list-ref file-info 3)))))
	  (if (< time-diff-pre most-time) time-diff-pre most-time))))

  ;; return "pass" or "fail" (actually, the contents of the link)
  ;; 
  (define (build-result-from-link file-info page-type is-python?)
    (let* ((link (car (cdr (list-ref file-info 0)))))
      (format #t "INFO:: getting build info from remote: ~s~%" link)
      (if (not (list? link))
	  #f
	  (let* ((url (string-append (car (cdr link))
				     "-" 
				     (symbol->string page-type)
				     (if is-python? "-python" "")
				     "-status"))
		 (res (call-with-output-string
		       (lambda (port)
			 (let ((initial-output-port (current-output-port)))
			   (set-current-output-port port)
			   (http-get url)
			   (set-current-output-port initial-output-port))))))
	    (if (string? res)
		(let ((l (string-length res)))
		  (if (not (> l 0))
		      ""
		      (if (< l 20)
			  (substring res 0 (- l 1))
			  "missing-file")))))))) ; we don't want to see
					         ; web server page
					         ; missing raw HTML


  ;; return some text
  (define (source-build-status-text)
    (let ((file (append-dir-file
		 (append-dir-dir "/y/people/emsley" 
				 "public_html/software/pre-release/build-status")
		 "status")))
      (call-with-input-file file
	(lambda (port)
	  (let* ((t1 (read-line port))
		 (t2 (read-line port))
		 (text (if (eof-object? t2)
			   t1
			   (string-append t1 " " t2))))
	    (markup text))))))
  
  (define (markup text)
    (if (not (string? text))
	""
	(cond 
	 ((string-match "pass" text)
	    `(b (font (@ color "#209920")
		      ,text)))
	 ((string-match "progress" text)
	    `(b (font (@ color "#606960")
		      ,text)))
	 (else 
	  `(b (font (@ color "#bb2020")
		    ,text))))))
	
  ;; page-type is 'build or 'test.
  ;; return "pass" "fail" or " ".
  (define (latest-build-result file-info page-type)
    
    (define (get-mtime file-name)
      (stat:mtime (stat file-name)))

    ;; is the file test-status-file-name more than 20 seconds older
    ;; than build-status-file-name (i.e. the test ran and the build
    ;; status is probably valid)
    ;;
    ;; return #t or #f
    ;; 
    (define (newer-file? test-status-file-name build-status-file-name)
      (if (not (file-exists? test-status-file-name))
	  #f
	  (if (not (file-exists? build-status-file-name))
	      #f
	      (let ((mtime-test  (get-mtime  test-status-file-name))
		    (mtime-build (get-mtime build-status-file-name)))
		(> mtime-test (+ ctime-build 20))))))

    (define (get-file link page-type)
      (if (list? link) 
	  'link
	  (string-append "/y/people/emsley/public_html/build-logs/"
			 link
			 (if (eq? page-type 'build)
			     "-build-status"
			     "-test-status"))))

    ;; main line of latest-build-result
    ;;
    ;; (format #t "file-info ~s~%" file-info)
    ;; 
    (let* ((link (car (cdr (list-ref file-info 0))))
	   (file-names-info (car file-info))
	   (is-python? (if (= (length file-names-info) 3)
			   (car (cdr (cdr (car file-info))))
			   #f))
	   (build-status-file (get-file link 'build))
	   ( test-status-file (get-file link 'test ))
	   (file (get-file link page-type)))

      (if (eq? file 'link)
	  (let ((br (build-result-from-link file-info page-type is-python?)))
	    (markup br))

	  (if (not file)
	      "-" 
	      (if (not (file-exists? file))
		  (begin 
		    (format #t "file does not exist ~s~%" file)
		    '(b  (font (@ color "#bb2020")) "not-found"))
		  (call-with-input-file file
		    (lambda (port)
		      (let ((text (read-line port)))
			(markup text)))))))))
			

  ;; string is /y/people/emsley/public_html... etc and make a URL
  ;; suffix (post server text)
  (define (url-from-file-name file-name person)
    (string-append "/~" person (substring file-name (+ 22 (string-length person)))))
    
  (define (get-binary-tar-file-url file-info person)
    ; (format #t "get-binary-tar-file-url from ~s~%" file-info)
    (let* ((bits (list-ref file-info 3))
	   (file-name (list-ref bits 2)))
      (url-from-file-name file-name person)))

  ;; 
  (define (format-binary-cell file-info now-time)
    ;;(format #t "formatting binary cell ~s~%" file-info)
    (if (not file-info)
	" "
	(let* (
	       (person (bin->person (car file-info))) ; fails on bad file-info?
	       (make-links? (car (cdr (car file-info))))
	       (time-diff (make-time-diff file-info now-time))
	       (n-for-spaces (time-diff->n time-diff))
; 	       (nov1 (format #t "...... n-for-spaces ~s~%" n-for-spaces))
	       (ns (n-bars n-for-spaces))
	       ; (nov2 (format #t ".......~s~%" (my-ns '(*ENTITY* "nbsp") n-for-spaces)))
	       (entities (my-ns '(*ENTITY* "nbsp") n-for-spaces))
	       (colour (cond
			((eq? time-diff #f) "#303030")
			((< time-diff (* 60 60 24)) "#80dd80")
			(else 
			 "#dd8080"))))

	  (if (= (length file-info) 2)
	      ;; a "not found" binary
	      `(p 
		,(car (list-ref file-info 0)) '(b " not found")
		(br)
		(table (@ (border 1) (bgcolor "#303030"))
		       (tr (td ,(n-bars 100))))
		;; this is pythonic, right?
		(a (@ href ,(build-log-page file-info 'build #t)) build-log)
		" "
		(a (@ href ,(build-log-page file-info 'test #t)) test-log))

	      ;; a regular binary cell
	      (let ((pythonic-build? (if (= (length (car file-info)) 2)
					 #t
					 (car (cdr (cdr (car file-info)))))))
		`(p 
		  (a (@ href ,(get-binary-tar-file-url file-info person)) 
		     ,(car (car file-info))) ; binary type string
		  " "
		  ,(car (cdr file-info)) ; revision number
		  ;; entities to pad the table
		  (*ENTITY* "nbsp")
		  (*ENTITY* "nbsp")
		  (*ENTITY* "nbsp")
		  (br)

		  ,(latest-build-result file-info 'build)
		  " " 
		  ,(latest-build-result file-info 'test)
		  (*ENTITY* "nbsp")
		  (*ENTITY* "nbsp")
		  ,(if make-links?
		       `(a (@ href ,(build-log-page file-info 'build pythonic-build?)) build-log)
		       "")
		  " "
		  ;; a quicker way to get the dir of the xx-yyy.txt files:
		  ,(if make-links? 
		       `(a (@ href ,(build-log-page file-info 'build-dir pythonic-build?)) build-dir)
		       "")
		  " "
		  ,(if make-links?
		       `(a (@ href ,(build-log-page file-info 'test pythonic-build?)) test-log)
		       "")
		  " "
		  ,(if make-links?
		       (let* ((bits (list-ref file-info 3))
;			    (nov-1 (format #t "==== pythonic?: ~s ~%" 
;					   (car file-info)))
;			    (nov-2 (format #t "==== pythonic - sample?:  ~s~%" 
;					   (if (= (length (car file-info)) 2)
;					       "nothing"
;					       (car (cdr (cdr (car file-info)))))))
			      (nov-3 (format #t "==== pythonic : ~s~%" 
					     pythonic-build?)))

			 (if pythonic-build?
			     `(a (@ href ,(build-log-page file-info 'test-python #t)) python-test-log)
			     ""))
		       "")
		  (table (@ (border 1) (bgcolor ,colour))
			 (tr ,entities))
		  ,(time-text (list-ref file-info 2))
		  ))))))


  ;; 
  (define (pair-up ls)
    
    (let loop ((ls ls)
	       (pairs '())
	       (running-item #f))
      (cond 
       ((null? ls) (reverse (if running-item
				(cons
				 (list running-item #f)
				 pairs)
				pairs)))
       (running-item
	(loop (cdr ls) 
	      (cons (list running-item (car ls))
		    pairs)
	      #f))
       (else (loop (cdr ls)
		   pairs
		   (car ls))))))

  (define (bin->person bin)
    (if (string=?  (car bin) "WinCoot")
	"lohkamp"
	"emsley"))


  ;; make-page main body
  ;; 
  (let ((svn-log-page "http://www.ysbl.york.ac.uk/~emsley/software/pre-release/svn-log")
	(latest-source-info (latest-tar-gz source-tar-dir "" ".tar.gz"))
	(binary-file-infos
	 (map
	  (lambda (bin)
	    (let* ((suffix-string (if (string=? (car bin) "WinCoot")
				      ".exe"
				      ".tar.gz"))
		   (person (bin->person bin))
		   (l (latest-tar-gz (bin-tar-dir person)
				     (car bin) 
				     suffix-string)))
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
	     (h2 "Coot SVN and Build Summary")
	     (p ("Generated " ,(strftime "%a %d %b %H:%M:%S %G %Z" (localtime (current-time)))))
	     "\n"
	     ;; repository version
	     (p ("SVN Repository Revision: " 
		 ,(get-svn-revision) 
		 " "
		 (a (@ href ,svn-log-page) "svn log")))
	     "\n"

	     ;; source code
	     (p "Source code "
		,(if (not latest-source-info)
		     "no latest prerelease source tar"
		     (let ((f (basename (list-ref latest-source-info 2) ".tar.gz"))
			   (l (url-from-file-name (list-ref latest-source-info 2) 
						  "emsley")))
		       `(a (@ href ,l) ,f)))
		" "
		,(list-ref latest-source-info 1)
		(br)
		,(time-text (oldness-info latest-source-info))
		" " " " " "
		,(source-build-status-text)
		(*ENTITY* "nbsp")
		(*ENTITY* "nbsp")
		(a (@ href "http://www.ysbl.york.ac.uk/~emsley/build-logs/bubbles-local-build.log")
		   "source build log"))

	     (p (a (@ href "burn-up.png") "Burn-up chart"))

	     ;; binary targets
	     (table 
	      (@ border 1)
	      ,(let ((now-time (current-time)))
		 (map (lambda (file-info-pair)
			(append 
			 `(tr
			   (td ,(format-binary-cell (car  file-info-pair) now-time))
			   (td ,(format-binary-cell (cadr file-info-pair) now-time)))
			 (list "\n"))) ; ease of reading
		      (pair-up binary-file-infos))))))
     file-name)))
      
     

(let ((file-name "/y/people/emsley/public_html/coot/build-info.html")
      (bin-list 
       (list 

; bragg1 is no longer a fedora-6 :-(
;	(list "binary-Linux-i386-fedora-6" 
;	      "Linux-bragg1.chem.york.ac.uk/bragg1.chem.york.ac.uk")
;	(list "binary-Linux-i386-fedora-6-python"  
;	      "Linux-bragg1.chem.york.ac.uk/bragg1.chem.york.ac.uk")
;	(list "binary-Linux-i386-fedora-6-python-gtk2"
;	      "Linux-bragg1.chem.york.ac.uk/gtk2")

	(list "binary-Linux-i386-fedora-3"
	      "Linux-bunyip.chem.york.ac.uk/gtk1" #f)
	(list "binary-Linux-i386-fedora-3-python"
	      "Linux-bunyip.chem.york.ac.uk/gtk1" #t)

	(list "binary-Linux-i386-fedora-4-python-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-cycle/gtk2" #t))

;	(list "binary-Linux-i386-fedora-6-python-gtk2"
;	      "Linux-setsuko.chem.york.ac.uk/gtk2" #t)

;	(list "binary-Linux-i386-fedora-6"
;	      "Linux-setsuko.chem.york.ac.uk/gtk1" #f)

;	(list "binary-Linux-i386-fedora-6-python"
;	      "Linux-setsuko.chem.york.ac.uk/gtk1" #t)

	(list "binary-Linux-i386-fedora-8-python-gtk2"
	      "fedora-8/gtk2" #t)

	(list "binary-Linux-i386-fedora-10-python-gtk2"
	      "fedora-10/gtk2" #t)

;	(list "binary-Linux-i386-fedora-8"
;	      "Linux-dragon.chem.york.ac.uk/gtk1" #f)

;	(list "binary-Linux-i386-fedora-8-python-gtk2"
;	      "Linux-dragon.chem.york.ac.uk/gtk2" #t)

	(list "binary-Linux-i386-redhat-8.0"
	      "Linux-bubbles/gtk1" #f)

	(list "binary-Linux-i386-redhat-8.0-python" 
	      "Linux-bubbles/gtk1" #t)

	(list "binary-Linux-x86_64-rhel-4-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-lemur/gtk2" #f))

	(list "binary-Linux-x86_64-rhel-4-python-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-lemur/gtk2" #t))

	(list "binary-Linux-i386-centos-4-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-jackal/gtk2" #f))

	(list "binary-Linux-i386-centos-4-python-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-jackal/gtk2" #t))

	(list "binary-Linux-x86_64-centos-5-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-koala.bioch/gtk2" #f))

	(list "binary-Linux-x86_64-centos-5-python-gtk2"
	      (list 'absolute "http://www.biop.ox.ac.uk/emsley/build-logs/Linux-koala.bioch/gtk2" #t))

	(list "binary-Linux-i686-ubuntu-6.06.1-python-gtk2" 
	      "ubuntu-6.06/gtk2" #t)

	(list "binary-Linux-i686-ubuntu-8.04.2" 
	      "Linux-bragg3/gtk1" #f)
	(list "binary-Linux-i686-ubuntu-8.04.2-python" 
	      "Linux-bragg3/gtk1" #t)
	(list "binary-Linux-i686-ubuntu-8.04.2-python-gtk2"
	      "Linux-bragg3/gtk2" #t)

	(list "WinCoot" 
	      (list 'absolute 
		    "http://www.ysbl.york.ac.uk/~lohkamp/build-logs/MINGW32_NT-5.1-sarabellum/gtk2" #t))

	; no guile-gtk for this one.
	; (list "binary-Linux-i386-fedora-3-python-gtk2"
        ;     "Linux-bunyip.chem.york.ac.uk/gtk2")

	)))

  (make-page bin-list file-name))


