
(use-modules (srfi srfi-13)
             (srfi srfi-1)
             (sxml simple)
	     (ice-9 string-fun)
	     (ice-9 popen)
	     (ice-9 rdelim)
	     (ice-9 receive)
	     (ice-9 format)
	     ;; (net http)
	     (www main)
	     (os process)
	     (ice-9 regex))

(use-syntax (ice-9 syncase))

(define devel-dir  "http://lmb.bioch.ox.ac.uk/coot/devel")
;;
(define svn-log-page  "http://lmb.bioch.ox.ac.uk/coot/devel/svn.log")
;; 
(define source-tar-dir (string-append
			(getenv "HOME")
			"/public_html"
			"/coot/software/source/pre-releases/"))
;; 
(define web-source-tar-dir "http://lmb.bioch.ox.ac.uk/coot/software/source/pre-releases/")


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


;; or define a utility function for this
(define (directory-files dir)

  (if (not (access? dir R_OK))
    '()
    (let ((p (opendir dir)))
      (do ((file (readdir p) (readdir p))
           (ls '()))
          ((eof-object? file) (closedir p) (reverse! ls))
        (set! ls (cons file ls))))))

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

;;;
(define (append-dir-file dir-name file-name)
  (string-append (directory-as-file-name dir-name) "/" file-name))

;; split input string into a number of strings - the split character
;; is a newline.
;;
;; a build info is (e.g. 
;; 	(list "binary-Linux-i686-ubuntu-8.04.2" 
;;	 "http://xx/build-logs/Linux-a/gtk1" 
;;	 "http://yy/binaries/nightlies/pre-release"
;;        #f #t)
;; 
;; Cache the result, so that we only do the http-get once for each web-dir
;; 
(define server-html-binary-page-lines 
  (let ((cached-pages '()))
    (lambda (build-info)
      
      (let ((web-dir (list-ref build-info 2)))
	(format #t "getting web page ~s for ~s" web-dir (car build-info))
	(let ((page (assoc web-dir cached-pages)))
	  (if page
	      (begin 
		(format #t " (cached)~%")
		page)
	      (begin
		(format #t "~%")
		(let* ((url (list-ref build-info 2))
                       ;; (nov (format #t "debug:: getting url ~s~%" url))
		       (s (www:get url))
			;;(call-with-output-string
			 ;;(lambda (port)
			  ;; (let ((initial-output-port (current-output-port)))
			   ;;  (set-current-output-port port)
			    ;; ;; (http-get url)
			     ;; (http-get "http://lmb.bioch.ox.ac.uk") 
			     ;;(www:get url)
			     ;;(set-current-output-port initial-output-port)
                             ;;))))
                       ;; (nov-2 (format #t "================== got http string : ~s~%" s))
		       (s-lines (string-split s #\newline)))
                  ;; (format #t "debug2:: done getting url~%")
		  (set! cached-pages (acons web-dir s-lines cached-pages))
		  s-lines))))))))

;; Return the revision (as a number) 
;; Return #f on failure.
;; 
(define (get-revision-from-tar binary-tar-file-name)
  (let ((bits (string-split binary-tar-file-name #\-)))
    ;; (format #t "bits in get-revision-from-binary-tar: ~s~%" bits)
    (let loop ((bits bits))
      (cond
       ((null? bits) (values #f #f))
       ((string-match "revision" (car bits))
	(if (> (length bits) 1)
	    (let ((revision+exe (car (cdr bits)))) ;; e.g. "1234.exe" for WinCoot
	      (let ((rev-bits (string-split revision+exe #\.)))
		(string->number (car rev-bits))))
	    #f))
       (else 
	(loop (cdr bits)))))))

;; Return the date the file was made from the tar file name
;; or #f
;;
(define (get-date-from-tar source-tar-file-name)
  (let ((full-path (string-append source-tar-dir source-tar-file-name)))

    (if (file-exists? full-path)
	(let ((mtime (stat:mtime (stat full-path))))
	  ;; (format #t "~s mtime ~s~%" full-path mtime)
	  (strftime "%a %d %b %H:%M:%S %G %Z" (localtime mtime)))
	(begin 
	  (format #t "source tar ~s not found" full-path)
	  #f))))




;; Return the url and the revision number (multiple values) of the
;; latest file in the http directory that matches build-info.  Return
;; #f #f any failure.
;; 
;; a build info is (e.g. 
;; 	(list "binary-Linux-i686-ubuntu-8.04.2" 
;;	 "http://xx/build-logs/Linux-a/gtk1" 
;;	 "http://yy/binaries/nightlies/pre-release"
;;        #f #t)
;;
;; 
;;
;;
(define (get-latest-binary-url-info build-info) 
  
  ;; In future, cache the result (the list of
  ;; coot-xxx-binary-xxx-xxx-tar.gz file names)
  ;;
  ;; python? is a flag to help out with the binary matching.  When
  ;; python? is #t "..python.." should be part of the file name (and
  ;; when it is #f it should not, of course).  This stops the
  ;; xx-python.tar.gz files matching the non-python build tar files.
  ;; 
  (define (get-binary-tar-file-names s-lines binary-match python? gtk2?)

    (let loop ((s-lines s-lines)
	       (binary-tar-file-names '()))
      
      (cond 
       ((null? s-lines) ;; (format #t "returning ~s~%" (reverse binary-tar-file-names))
                       	(reverse binary-tar-file-names))
       (else 
	(let ((line (car s-lines)))
	  (let ((line-bits (string-split line #\")))
	    ;; (format #t "line-bits: ~s~%" line-bits)
	    (if (> (length line-bits) 6)
		(let ((binary-file-name-bit (list-ref line-bits 6)))
		  ;; (format #t "binary-file-name-bit: ~s~%" binary-file-name-bit)
		  (if (or (string-match "coot" binary-file-name-bit)
			  (string-match "WinCoot" binary-file-name-bit))
		      (let ((bits-2 (string-split binary-file-name-bit #\>)))
			(if (> (length bits-2) 1)
			    (let ((bits-3 (string-split (list-ref bits-2 1) #\<)))
			      ;; (format #t "bit-3: ~s~%" bits-3)
			      (if (> (length bits-3) 1)
				  (let ((a-file-name (car bits-3)))
				    ;; (format #t "a-file-name: ~s~%" a-file-name)
				    ;; (format #t "::::: ~s~%" (string-match ".md5sum" a-file-name))
				    (if (not (string-match ".md5sum" a-file-name))
					(if (string=? "WinCoot" binary-match)

					    ;; WinCoot
					    (if (string-match ".exe" a-file-name)
						(loop (cdr s-lines) (cons a-file-name
									  binary-tar-file-names))
						(loop (cdr s-lines) binary-tar-file-names))
					    
					    ;; not WinCoot, e.g. Ubuntu
					    (if (string-match ".tar.gz" a-file-name)
						(if (string-match binary-match a-file-name)
						    (if (or (and (string-match "python" a-file-name)
								 python?)
							    (and (not (string-match "python" a-file-name))
								 (not python?)))
							(if (or (and (string-match "gtk2" a-file-name)
								     gtk2?)
								(and (not (string-match "gtk2" a-file-name))
								     (not gtk2?)))
							    (loop (cdr s-lines)
								  (cons a-file-name binary-tar-file-names))
							    (loop (cdr s-lines) binary-tar-file-names))
							(loop (cdr s-lines) binary-tar-file-names))
						    (loop (cdr s-lines) binary-tar-file-names))
						(loop (cdr s-lines) binary-tar-file-names)))
					(loop (cdr s-lines) binary-tar-file-names)))
				  (loop (cdr s-lines) binary-tar-file-names)))
			    (loop (cdr s-lines) binary-tar-file-names)))
		      (loop (cdr s-lines) binary-tar-file-names)))
		(loop (cdr s-lines) binary-tar-file-names))))))))


  ;; main line of get-latest-binary-url-info
  ;; 
  (let ((s-lines (server-html-binary-page-lines build-info))
	(binary-match (car build-info))
	(python? (list-ref build-info 3))
	(gtk2?   (list-ref build-info 4)))
    
    ;; (format #t "s-lines: ~s ~s lines~%" s-lines (length s-lines))
    (let ((file-ls (get-binary-tar-file-names s-lines binary-match python? gtk2?)))
      ;; (format #t "file-ls: ~s ~%" file-ls)
      
      (let loop ((file-ls file-ls)
		 (latest-revision #f)
		 (latest-tar-file #f))
	(cond
	 ((null? file-ls)
	  (if (number? latest-revision)
	      (values latest-tar-file latest-revision)
	      (values #f #f)))
	 (else 
	  (let ((revision-number (get-revision-from-tar (car file-ls))))
	    (if (not (number? revision-number))
		(loop (cdr file-ls) latest-revision latest-tar-file)
		(if (not (number? latest-revision))
		    (loop (cdr file-ls) revision-number (car file-ls))
		    (if (> revision-number latest-revision)
			(loop (cdr file-ls) revision-number (car file-ls))
			(loop (cdr file-ls) latest-revision latest-tar-file)))))))))))

		
  

(define (get-url url)
  (call-with-output-string
   (lambda (port)
     (let ((initial-output-port (current-output-port)))
       (set-current-output-port port)
       (www:get url)
       (set-current-output-port initial-output-port)))))

;; return the revision number or #f
;; 
(define (get-svn-revision)
  (let ((coot-dir (string-append (getenv "HOME") "/Projects/coot-svn/trunk"))
	(current-dir (getcwd)))
    ; (format #t "coot-dir: ~s~%" coot-dir)
    ; (format #t "curr-dir: ~s~%" current-dir)
    (chdir coot-dir)
    (let ((s (shell-command-to-string "svn -qu status")))
      (chdir current-dir)
      (string->number
       (car (reverse (string-split (car (cdr (reverse (string-split s #\newline)))) #\space)))))))





;;
;;
(define (print-build-info record port)
  (let ((rec-type (record-type-descriptor record)))
    (format port "\nbinary-type:      ~s~%" ((record-accessor rec-type 'binary-type)    record))
    (format port "binary-http-dir:  ~s~%" ((record-accessor rec-type 'binary-http-dir)  record))
    (format port "build-log-prefix: ~s~%" ((record-accessor rec-type 'build-log-prefix) record))
    (format port "binary-url:       ~s~%" ((record-accessor rec-type 'binary-url)       record))
    (format port "build-status:     ~s~%" ((record-accessor rec-type 'build-status)     record))
    (format port "test-status:      ~s~%" ((record-accessor rec-type 'test-status)      record))
    (format port "latest-bin-rev:   ~s~%" ((record-accessor rec-type 'latest-binary-revision) record))))

;; return a list (pair) of lists of whateverever is in ls.  Fill last
;; pair with #f if needed.
;; 
(define (pair-up ls)
  (let loop ((ls ls) (pairs '()) (running-item #f))
    (cond 
     ((null? ls) (reverse (if running-item
			      (cons (list running-item #f)
				    pairs) pairs)))
     (running-item
      (loop (cdr ls) (cons (list running-item (car ls)) pairs) #f))
     (else (loop (cdr ls) pairs (car ls))))))


;; return a table
(define (make-binary-table-from-records records source-code-revision-number)

  ;; 
  (define (colour-from-revision-number bin-rev source-rev)
    (if (= bin-rev source-rev)
	"#057705" ;; dark green
	(let* ((diff-1 (- source-rev bin-rev))
	       (max-diff 30)
	       (frac-diff (/ 1 max-diff))
	       (diff-2 (if (> diff-1 max-diff) max-diff diff-1)))
	  (if (< diff-2 0)
	      "#101010" ;; dark grey, for a strange situation
	      (let ((blue 0.05)
		    (red (+ 0.05 (* 0.9 diff-2 frac-diff)))
		    (green (- 0.9 (* 0.9 diff-2 frac-diff))))
		(apply string-append
		       (cons "#" (map (lambda(c) 
					(let ((s (format #f "~x" (inexact->exact (round (* 255 c))))))
					  (if (= (string-length s) 1)
					      (string-append "0" s)
					      s)))
				      (list red green blue)))))))))

  ;; 
  (define (format-binary-cell binary-record now-time)
    (if (not binary-record)
	" " ;; the last (unmatched) "fill-in" cell
	(let* ((rec-type (record-type-descriptor binary-record))
	       (binary-url ((record-accessor rec-type 'binary-url) binary-record)))
	  `(p 
	    ;; if the binary-tar-file is missing (typically target not
	    ;; built for a long time and the ones that did exist have
	    ;; been deleted), then we don't want to provide a link to
	    ;; the binary, just put the text of the binary type, not
	    ;; hyper-linked.
	    ,(if (string? binary-url)
		 `(a (@ href ,binary-url) ,((record-accessor rec-type 'binary-type) binary-record))
		 ((record-accessor rec-type 'binary-type) binary-record))
	      (*ENTITY* "nbsp")
	      (*ENTITY* "nbsp")
	      ,(let ((binary-revision ((record-accessor rec-type 'latest-binary-revision) binary-record)))
		 (if (not (number? source-code-revision-number))
		     (if (not (number? binary-revision))
			 "No-source-rev,no-bin-rev"
			 binary-revision)
		     (if (not (number? binary-revision))
			 "Missing-bin-rev"
			 (let ((font-colour (colour-from-revision-number 
					     binary-revision source-code-revision-number)))
			     `(b (font (@ color ,font-colour)) ,binary-revision)))))))))
			   
		     
  ;; main line of make-binary-table-from-records
  ;; 
  `(table 
    (@ border 1)
    ,(let ((now-time (current-time)))
       (map (lambda (record-pair)
	      (append
	       `(tr
		 (td ,(format-binary-cell (car  record-pair) now-time))
		 (td ,(format-binary-cell (cadr record-pair) now-time)))
	       (list "\n"))) ; ease of reading
	    (pair-up records)))))
  


;; colour up the text in a manner that depends on the text (and bold
;; the text).  This used to be a part of
;; make-builds-table-from-records, but it is useful to colour up the
;; results of the source code build too.

(define (colourize text)

  (cond
   ((not (string? text)) 	"status not found")
   ((string-match "404 Not Found" text)
    `(b (font (@ color "#991111") "Not Found")))
   ((string-match "Object not found!" text)
    `(b (font (@ color "#991111") "")))
   ((string-match "pass" text)
    `(b (font (@ color "#119911") ,text)))
   ((string-match "fail" text)
    `(b (font (@ color "#991111") ,text)))
   ((string-match "progress" text)
    `(b (font (@ color "#444444") ,text)))
   ((string-match "waiting" text)
    `(b (font (@ color "#444444") ,text)))
   (else 
    `(b ,text))))


(define (make-builds-table-from-records records source-code-revision-number)
	
  ;; 
  (define (format-binary-cell binary-record now-time)
    (if (not binary-record)
	" " ;; the last (unmatched) "fill-in" cell
	
	;; don't forget that we can come here with binary-record full of '#f's.
	(let* ((rec-type (record-type-descriptor binary-record))
	       (build-log-stub ((record-accessor rec-type 'build-log-prefix) binary-record))
	       (build-log-dir (if (string? build-log-stub)
				  (dirname build-log-stub)
				  "missing-build-log-prefix")))

	  `(p (a (@ href ,build-log-dir)
		 ,((record-accessor rec-type 'binary-type) binary-record))
	      " " 
	      ,(colourize ((record-accessor rec-type 'build-status) binary-record))
	      " "
	      ,(colourize ((record-accessor rec-type 'test-status) binary-record))))))
		 

  ;; main-line of make-builds-table-from-records
  `(table
    (@ border 1)
    ,(let ((now-time (current-time)))
       (map (lambda (record-pair)
	      (append
	       `(tr
		 (td ,(format-binary-cell (car  record-pair) now-time))
		 (td ,(format-binary-cell (cadr record-pair) now-time)))
	       (list "\n"))) ; ease of reading
	    (pair-up records)))))
    


;; 
(define (svn-details) 
  `(p ("SVN " (a (@ href "http://coot.googlecode.com/svn/") "Repository") " Revision: " 
       ,(get-svn-revision) 
       " "
       (a (@ href ,svn-log-page) "svn log"))))


;; Return values: the source code url, a file-name (for the link) and
;; a revision number (a number) and a date (string) for the tar file.
;; 
;; This is based on analysis of the file system.  
;; 
;; You can replace this function by something that looks at
;; web-source-tar-dir, you will have to write that.
;;
(define (source-code-url-info) 

  (let ((files (directory-files source-tar-dir)))
    (let loop ((files files)
	       (latest-file #f)
	       (latest-revision-number #f)
	       (tar-file-date #f))

      (cond
       ((null? files)
	;; (format #t "source-code-url-info returns using ~s ~s~%" latest-file latest-revision-number)
	(if (not latest-file)
	    (values #f #f #f)
	    (values (string-append web-source-tar-dir latest-file)
		    latest-file
		    latest-revision-number
		    tar-file-date)))
       ((string-match "md5sum" (car files))
	(loop (cdr files) latest-file latest-revision-number tar-file-date))
       ((string-match "coot" (car files))
	(if (not (string-match ".tar.gz" (car files)))
	    (loop (cdr files) latest-file latest-revision-number tar-file-date)
	    (let ((rev (get-revision-from-tar (car files))))
	      (if (not (number? rev))
		  (loop (cdr files) latest-file latest-revision-number tar-file-date)
		  (if (or (not (number? latest-revision-number))
			  (> rev latest-revision-number))
		      (let ((date (get-date-from-tar (car files))))
			(loop (cdr files) (car files) rev date))
		      (loop (cdr files) latest-file latest-revision-number tar-file-date))))))
       (else 
	(loop (cdr files) latest-file latest-revision-number tar-file-date))))))
  

;; return a sxml paragraph
;; 
(define (source-code-details)

  ;; 
  (define (get-running-status-contents)
    (let* ((running-status-file (append-dir-file devel-dir "source-build-running"))
	   (s (www:get running-status-file)))
      (if (string-match "Not Found" s)
	 #f
	 s)))
    
  
  ;; main line of source-code-details
  `(p 
    ,(receive (source-code-url source-code-file-name revision-number source-tar-date)
 	      (source-code-url-info)
	      `("Source code: "
		(a (@ href ,source-code-url) ,source-code-file-name)
		" "
		(a (@ href "latest-coot-build.log") "source build log")
		" "
		,(let ((s (get-running-status-contents)))
		   (if (string? s)
		       (colourize s)
		       ""))
		" "
		,revision-number
		(br)
		,(if (string? source-tar-date)
		     source-tar-date
		     "")
		))))

;; 
(define (burn-up-chart)
  '(p "Burn-up graph " (br) (a (@ href "burn-up.png") (image (@ src "burn-up-icon.png")))))


;; 
(define (write-sxml block filename)
  
  (call-with-output-file filename
      (lambda (port)
	(sxml->xml block port))))


;; a build info is, e.g. 
;; 	(list "binary-Linux-i686-ubuntu-8.04.2" 
;;	 "http://xx/build-logs/Linux-a/gtk1" 
;;	 "http://yy/binaries/nightlies/pre-release"
;;        #f #t)
;;
;; return a record (possibly full of #fs)
;; 
(define (fill-record build-info build-info-record-type make-build-record)

  ;; main line of fill-record
  ;; 
  (receive (tar-file revision-number)
	   (get-latest-binary-url-info build-info)

	   (let* ((binary-type      (car build-info))
		  (build-log-prefix (list-ref build-info 1))
		  (is-python?       (list-ref build-info 3))
		  (is-gtk2?         (list-ref build-info 4))
		  ;; python is part of the build/test status file name
		  ;; when we are looking at python builds (except for
		  ;; WinCoot, which does not use it).
		  (python-status-extra (if is-python? 
					   (if (string-match "WinCoot" binary-type)
					       ""
					       "-python")
					   ""))
		  (build-status-link (string-append build-log-prefix 
						    python-status-extra
						    "-build-status"))
		  (test-status-link (string-append build-log-prefix 
						   python-status-extra
						   "-test-status")))
  
	     (if (not (and (string? tar-file)
			   (number? revision-number)))
		 ;; make a dummy record then
		 (apply make-build-record (append build-info (list #f #f #f #f)))
		 (let ((binary-url (string-append 
				    (list-ref build-info 2) ;; "/" is included in pre-release dir, right?
				    tar-file))
		       (build-status (www:get build-status-link))
		       ( test-status (www:get  test-status-link)))

		   (format #t "============== test-status  for ~s is ~s~%" test-status-link test-status)
		   (format #t "============== build-status for ~s is ~s~%" build-status-link build-status)

		   ;; build-status is e.g. "passed build" or #f
		   ;;  test-status is e.g. "passed test" or #f
		   ;; 
		   (apply make-build-record (append build-info
						    (list revision-number
							  binary-url
							  build-status
							  test-status))))))))
      


;;
;; 
(define (make-page build-info-list html-file-name)

  ;; first thing: for each of the build-info-list, get the file name
  ;; of the latest tar file and its revision number.

  (let* ((record-type (make-record-type "build-info" 
					'(binary-type
					  build-log-prefix 
					  binary-http-dir 
					  is-python? 
					  is-gtk2?
					  latest-binary-revision
					  binary-url
					  build-status
					  test-status)
					print-build-info))
	 (make-build-record (record-constructor record-type))
	 (records (map (lambda(i) (fill-record i record-type make-build-record)) build-info-list)))
    
    (format #t "records: ~s~%" records)
    (format #t "getting source code info...~%")
    (receive (source-code-url source-code-file-name source-code-revision-number source-code-date)
	     (source-code-url-info)    
	     (write-sxml
	      `(html (head (title "Coot Build Summary Page")
			   (meta (@ (http-equiv refresh) (content 600))))
		     (link (@ (rel "icon")
			      (type "image/png")
			      (href "../coot-favicon.png")))
		     (body 
		      (h2 "Coot SVN and Build Summary")
		      (p ("Generated " ,(strftime "%a %d %b %H:%M:%S %G %Z" (localtime (current-time)))))
		      "\n"
		      ,(svn-details)
		      ,(source-code-details)
		      ,(burn-up-chart)
		      (h3 "Latest Binary Tars")
		      
		      ;;we want the souce code revision number so that we can
		      ;;colour the binary revision numbers in the binaries
		      ;;table.
		      ,(make-binary-table-from-records records source-code-revision-number))
		     
		     (h3 "Development Build Logs")
		     ,(make-builds-table-from-records records source-code-revision-number))
	      html-file-name))))
    
    



;; Do it.
;; 
;; element of the form
;; (list binary-tar-file-string-match
;;       web-log-stub
;;       binary-tar-file-dir
;;       is-a-python?
;;       is-a-gtk2-build?)
;;
(let* ((html-file-name (string-append (getenv "HOME") 
				      "/public_html/coot/devel/build-info.html"))
       (york-ubuntu-version "8.04.4")
       (build-list
	(list 
	 
	 (list "binary-Linux-i386-fedora-3"
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-bunyip.chem.york.ac.uk/gtk1"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #f #f)

	 (list "binary-Linux-i386-fedora-3-python"
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-bunyip.chem.york.ac.uk/gtk1"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #f)


;	 (list "binary-Linux-i386-fedora-4-python-gtk2"
;	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-cycle/gtk2" 
;	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
;	       #t #t)

;	 (list "binary-Linux-i386-fedora-4-gtk2"
;	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-cycle/gtk2" 
;	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
;	       #f #t)

	 (list "binary-Linux-i386-fedora-8-python-gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/fedora-8/gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #t)

	 (list "binary-Linux-i386-fedora-10-python-gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/fedora-10/gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #t)

	 ;; I make this occassionally.
	 (list "binary-Linux-i386-fedora-12-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-jackal-f12/gtk2"
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #f #t)

	;; kevin makes this occassionally.
 	(list "binary-Linux-i386-fedora-12-python-gtk2"
 	      "http://www.ysbl.york.ac.uk/~emsley/build-logs/fedora-12/gtk2"
  	      "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
 	      #t #t)

	 (list "binary-Linux-i386-centos-4-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-jackal/gtk2"
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #f #t)

	 (list "binary-Linux-i386-centos-4-python-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-jackal/gtk2"
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #t #t)

	 (list "binary-Linux-x86_64-centos-5-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-koala.bioch/gtk2" 
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #f #t)

	 (list "binary-Linux-x86_64-centos-5-python-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-koala.bioch/gtk2" 
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #t #t)

	 (list "binary-Linux-x86_64-rhel-4-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-lemur/gtk2" 
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #f #t)

	 (list "binary-Linux-x86_64-rhel-4-python-gtk2"
	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-lemur/gtk2" 
	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
	       #t #t)

	 (list "binary-Linux-i686-ubuntu-6.06.1-python-gtk2" 
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/ubuntu-6.06/gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #t)

	 (list (string-append "binary-Linux-i686-ubuntu-" york-ubuntu-version)
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-bragg3/gtk1" 
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #f #f)

	 (list (string-append "binary-Linux-i686-ubuntu-" york-ubuntu-version "-python")
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-bragg3/gtk1" 
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #f)

	 (list (string-append "binary-Linux-i686-ubuntu-" york-ubuntu-version "-python-gtk2")
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-bragg3/gtk2" 
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #t)

	 (list (string-append "binary-Linux-i686-ubuntu-" york-ubuntu-version "-gtk2")
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-bragg3/gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #f #t)

	 (list (string-append "binary-Linux-i386-rhel-5-python-gtk2")
	       "http://www.ysbl.york.ac.uk/~emsley/build-logs/Linux-krosp.chem.york.ac.uk/gtk2"
	       "http://www.ysbl.york.ac.uk/~emsley/software/binaries/nightlies/pre-release/"
	       #t #t)

; 	 (list "binary-Linux-x86_64-ubuntu-9.04-gtk2"
; 	       "http://lmb.bioch.ox.ac.uk/emsley/build-logs/Linux-scylla/gtk2" 
; 	       "http://lmb.bioch.ox.ac.uk/coot/software/binaries/pre-releases/" 
; 	       #f #t)

	 (list "WinCoot" 
	       "http://www.ysbl.york.ac.uk/~lohkamp/build-logs/MINGW32_NT-5.1-sarabellum/gtk2" 
	       "http://www.ysbl.york.ac.uk/~lohkamp/software/binaries/nightlies/pre-release/"
	       #t #t))))

  (make-page build-list html-file-name))
	
