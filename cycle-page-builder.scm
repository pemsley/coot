;;;; Copyright 2008 by The University of Oxford
;;;; Copyright 2013, 2015 by Medical Research Council
;;;;
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA

(use-modules (srfi srfi-13)
             (srfi srfi-1)
             (sxml simple)
	     (ice-9 string-fun)
	     (ice-9 popen)
	     (ice-9 rdelim)
	     (ice-9 receive)
	     (ice-9 format)
	     ;; (net http)

;; changed?
;	     (www main)
	     (www url)
	     (www http)

	     (os process)
	     (ice-9 regex))

(use-syntax (ice-9 syncase))

;; (define devel-dir  "http://lmb.bioch.ox.ac.uk/coot/devel")
(define devel-dir  "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/devel")


; no longer used.  Now we get source info from the web server, not the file system.
;
; ;; 
; (define source-tar-dir (string-append
; 			(getenv "HOME")
; 			"/public_html"
; 			"/coot/software/source/pre-releases/"))


;; this must be slashed
;; 
(define web-source-tar-dir "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/source/pre-releases/")


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




;; regex-split by Neil Jerram
;; 
(define (regex-split regex str . opts)
  (let* ((unique-char #\@)
         (unique-char-string (string unique-char)))
    (let ((splits (separate-fields-discarding-char
                   unique-char
                   (regexp-substitute/global #f
                                             regex
                                             str
                                             'pre
                                             unique-char-string
                                             0
                                             unique-char-string
                                             'post)
                   list)))
      (cond ((memq 'keep opts)
             splits)
            (else
             (let ((non-matches (map (lambda (i)
                                       (list-ref splits (* i 2)))
                                     (iota (floor (/ (1+ (length splits)) 2))))))
               (if (memq 'trim opts)
                   (filter (lambda (s)
                             (not (zero? (string-length s))))
                           non-matches)
                   non-matches)))))))



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
		       (url-slashed (string-append url "/"))
                       (nov (format #t "debug:: getting url ~s~%" url-slashed))

		       (s (get-url url-slashed)))

		  (if (not (= (vector-length s) 5))
		      #f
		      (let ((s-lines (string-split (vector-ref s 4) #\newline)))
			(set! cached-pages (acons web-dir s-lines cached-pages))
			s-lines))))))))))


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
        ;; (format #t "matched on revision: ~s~%" bits)
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

  ;; 
  (define (oxford-string->binary-tar-file binary-file-name-bit binary-match python? gtk2?)

    (if (or (string-match "coot" binary-file-name-bit)
	    (string-match "WinCoot" binary-file-name-bit))
	(let ((bits-2 (string-split binary-file-name-bit #\>)))

	  ;; (format #t "==== Here 2, bits-2 ~s ~%" bits-2)

	  (if (> (length bits-2) 1)
	      (let ((bits-3 (string-split (list-ref bits-2 1) #\<)))

		;; (format #t "==== Here 3 bit-3: ~s~%" bits-3)

		(if (> (length bits-3) 1)
		    (let ((a-file-name (car bits-3)))

		      ;; (format #t "**************** examining a-file-name: ~s~%" a-file-name)
		      
		      (cond 
		      ;; WinCoot
		      ((string-match ".exe" a-file-name)
		       a-file-name)
			  
		      ;; not WinCoot, e.g. Ubuntu
		      ((string-match ".tar.gz" a-file-name)
		       (if (not (string-match binary-match a-file-name))
			   (begin
			     ;; (format #t "+++++++++ fail 1 string-match ~s on ~s ~%" binary-match a-file-name)
			     #f)
			   (if (or (and (string-match "python" a-file-name)
					python?)
				   (and (not (string-match "python" a-file-name))
					(not python?)))
			       (if (or (and (string-match "gtk2" a-file-name)
					    gtk2?)
				       (and (not (string-match "gtk2" a-file-name))
					    (not gtk2?)))

				   (let ((dum 'dummy))
				     ;; (format #t "oxford-string->binary-tar-file returns ~s~%" a-file-name)
				     a-file-name))
			       (begin
				 (format #t "+++++++++ fail 2 ~%")
				 #f)
			       )))
		      (else
		       (begin
; 			 (format #t "+++++++++ oxford-string->binary-tar-file: fail 3 given: ~s~%~!" 
; 				 binary-file-name-bit)
			 #f)

		       )))

		    (begin
		      (format #t "+++++++++ fail 4 ~%")
		      #f)
		    ))
	      (begin
		(format #t "+++++++++ fail 5 ~%")
		#f)
	      ))
	#f))
  
  ;; return a string or #f
  ;; 
  (define (line->binary-tar-file line binary-match python? gtk2?)

    (let ((line-bits (string-split line #\")))

      ;; (format #t "    ----- binary-match: ~s line-bits: ~s~%" binary-match line-bits)

      (if (> (length line-bits) 6)
	  (let* ((tim-server (string-match "top" (list-ref line-bits 1)))
		 (tar-file-idx (if tim-server 7 6))
		 (binary-file-name-bit (list-ref line-bits tar-file-idx)))
	    
	    ;;(format #t "      ---- tim-server check binary-file-name-bit: idx: ~s string: ~s~%" 
            ;;       tar-file-idx binary-file-name-bit)
	    (if tim-server

		;; tim server
		(let ((smr-1 (string-match "debian-gnu-linux-" binary-file-name-bit)))

		  (if (not smr-1)
		      #f 
		      (begin
			;; (format #t "tim-server ===== pass 1! ~s ~s~%" 
			;; binary-file-name-bit smr-1)
			(let ((smr-2 (string-match ".tar.gz" binary-file-name-bit)))
			  (if (not smr-2)
			      #f 
			      binary-file-name-bit)))))

		;; oxford/york (-like) then
		;;
		(begin  
                   ;; (format #t "oxford/york-like~%")
		   (oxford-string->binary-tar-file binary-file-name-bit binary-match python? gtk2?))))

	  ;; it seems that the York server for WinCoot has changed.  Now file names are quoted with 's.
	  (if (not (string-match binary-match line))
	      #f
	      (let ((line-bits (string-split line #\')))
		;; (format #t "WinCoot line-bits: ~s~%" line-bits)
		(if (not (> (length line-bits) 2))
		    #f
		    (let ((file-name (list-ref line-bits 1)))
		      file-name)))))))

  
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
	(let ((binary-file-name (line->binary-tar-file (car s-lines) binary-match python? gtk2?)))
	  ;; (format #t "debug:: in get-binary-tar-file-names: binary-file-name: ~s~%" binary-file-name)
	  (if (string? binary-file-name)
	      (loop (cdr s-lines)
		    (cons binary-file-name binary-tar-file-names))
	      (loop (cdr s-lines) binary-tar-file-names)))))))
		    


  ;; main line of get-latest-binary-url-info
  ;; 
  (let ((s-lines (server-html-binary-page-lines build-info))
	(binary-match (car build-info))
	(python? (list-ref build-info 3))
	(gtk2?   (list-ref build-info 4)))
    
    ;; (format #t "  ----- in get-binary-tar-file-names s-lines: ~s ~s lines~%" s-lines (length s-lines))
    ;; (format #t " in get-binary-tar-file-names binary-match is ~s ~%" binary-match)

    (let ((file-ls (get-binary-tar-file-names s-lines binary-match python? gtk2?)))
      ;; (format #t "  ----- in get-binary-tar-file-names file-ls: ~s ~%" file-ls)
      
      (let loop ((file-ls file-ls)
		 (latest-revision #f)
		 (latest-tar-file #f))
	(cond
	 ((null? file-ls)
	  (if (number? latest-revision)
	      (values latest-tar-file latest-revision)
	      (values #f #f)))
	 (else 
	  (let ((revision-count (get-revision-from-tar (car file-ls))))
	    (if (not (number? revision-count))
		(loop (cdr file-ls) latest-revision latest-tar-file)
		(if (not (number? latest-revision))
		    (loop (cdr file-ls) revision-count (car file-ls))
		    (if (> revision-count latest-revision)
			(loop (cdr file-ls) revision-count (car file-ls))
			(loop (cdr file-ls) latest-revision latest-tar-file)))))))))))

		
  

(define (old-get-url url)
  (call-with-output-string
   (lambda (port)
     (let ((initial-output-port (current-output-port)))
       (set-current-output-port port)
       (www:get url)
       (set-current-output-port initial-output-port)))))

;; does this return a string?
;; No.
;; 
(define (get-url url)
  ;; (format #t "=== in get-url we get url: ~s~%" url)
  (http:request 'GET (url:parse url)  (list "User-Agent: Bespoke-0.0"
					    "Content-Type: text/plain")))


;; return the revision number or #f
; ;
(define (get-git-revision-count)
  (let ((coot-dir (string-append (getenv "HOME") "/Projects/for-source-tar/coot"))
	(current-dir (getcwd)))
    ; (format #t "coot-dir: ~s~%" coot-dir)
    ; (format #t "curr-dir: ~s~%" current-dir)
    (chdir coot-dir)
    (let ((s (shell-command-to-string "git rev-list --count HEAD")))
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
(define (make-binary-table-from-records records source-code-revision-count)

  ;; 
  (define (colour-from-revision-count bin-rev source-rev)
    
    (if (= bin-rev source-rev)
	"#057705" ;; dark green
	(let* ((diff-1 (- source-rev bin-rev))
	       (max-diff 30)
	       (frac-diff (/ 1 max-diff))
	       (diff-2 (if (> diff-1 max-diff) max-diff diff-1)))

	  (if (< diff-2 0)
	      "#101010" ;; dark grey, for a strange situation
	      (let ((blue 0.05)
		    (red   (+ 0.05 (* 0.9 diff-2 frac-diff)))
		    (green (- 0.9  (* 0.9 diff-2 frac-diff))))
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
 
                 ;; no colour means source-code-revision-count is #f
 
		 (if (not (number? source-code-revision-count))
		     (if (not (number? binary-revision))
			 "No-source-rev,no-bin-rev"
			 binary-revision)
		     (if (not (number? binary-revision))
			 "Missing-bin-rev"
			 (let ((font-colour (colour-from-revision-count 
					     binary-revision source-code-revision-count)))

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


(define (make-builds-table-from-records records source-code-revision-count)
	
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
	      ,(colourize ((record-accessor rec-type 'test-status) binary-record))
	      ))))
		 
  ;; --------------------------------
  ;; This for 2-column table
  ;; --------------------------------
  
;   ;; main-line of make-builds-table-from-records
;   `(table
;     (@ border 1)
;     ,(let ((now-time (current-time)))
;       (map (lambda (record-pair)
; 	      (append
; 	       `(tr
; 		 (td ,(format-binary-cell (car  record-pair) now-time))
; 		 (td ,(format-binary-cell (cadr record-pair) now-time)))
; 	       (list "\n"))) ; ease of reading
; 	    (pair-up records)))))

  ;; --------------------------------
  ;; 1-column table
  ;; --------------------------------
  `(table
    (@ border 1)
    ,(let ((now-time (current-time)))
        (map (lambda (record)
	       `(tr
		 (td ,(format-binary-cell record now-time))))
	     records))))

    


;; 
(define (git-revision-count-details) 
  `(p ((a (@ href "https://github.com/pemsley/coot") "Git Repository")
       " "
       " Revision Count: " ,(get-git-revision-count))))


;; Return values: the source code url, a file-name (for the link) and
;; a revision number (a number) and a date (string) for the tar file.
;; These values can be #f #f #f #f.
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
	       (latest-revision-count #f)
	       (tar-file-date #f))

      (cond
       ((null? files)
	(format #t "debug:: source-code-url-info returns using ~s ~s~%" latest-file latest-revision-count)
	(if (not latest-file)
	    (values #f #f #f #f)
	    (values (string-append web-source-tar-dir latest-file)
		    latest-file
		    latest-revision-count
		    tar-file-date)))
       ((string-match "md5sum" (car files))
	(loop (cdr files) latest-file latest-revision-count tar-file-date))
       ((string-match "coot" (car files))
	(if (not (string-match ".tar.gz" (car files)))
	    (loop (cdr files) latest-file latest-revision-count tar-file-date)
	    (let ((rev (get-revision-from-tar (car files))))
	      (if (not (number? rev))
		  (loop (cdr files) latest-file latest-revision-count tar-file-date)
		  (if (or (not (number? latest-revision-count))
			  (> rev latest-revision-count))
		      (let ((date (get-date-from-tar (car files))))
			(loop (cdr files) (car files) rev date))
		      (loop (cdr files) latest-file latest-revision-count tar-file-date))))))
       (else 
	(loop (cdr files) latest-file latest-revision-count tar-file-date))))))


;; Return values: the source code url, a file-name (for the link) and
;; a revision number (a number) and a date (string) for the tar file.
;; These values can be #f #f #f #f.
;;
;; use web-source-tar-dir
;; 
(define source-code-url-info
  (let ((cached #f))
    (lambda ()

      (if cached
	  cached
	  
	  (let* ((r (get-url web-source-tar-dir))
		 (s (vector-ref r 4)))

	    (let loop ((best-revision -1)
		       (best-date #f)
		       (best-file-name #f)
		       (lines (string-split s #\newline)))
	      (cond
	       ((null? lines)

		(set! cached
		      (if (> best-revision 0)
			  (values (string-append web-source-tar-dir "/" best-file-name) 
				  best-file-name best-revision best-date)
			  (values #f #f #f #f)))
		cached)
	       (else 
		(let* ((line (car lines))
		       (bits (regex-split  "<|>" line)))
		  (if (<= (length bits) 4)
		      (loop best-revision best-date best-file-name (cdr lines))

		      ;; 
		      (let ((potential-tar-string (list-ref bits 4)))
			(if (not (string-match ".tar.gz" potential-tar-string))
			    (loop best-revision best-date best-file-name (cdr lines))
			    (let ((tar-bits (string-split potential-tar-string #\-)))
			      ;; (format #t "bits: ~s~%" bits)
			      (if (not (> (length tar-bits) 4))
				  (loop best-revision best-date best-file-name (cdr lines))
				  (if (not (string-match ".tar.gz" (list-ref tar-bits 4)))
				      (loop best-revision best-date best-file-name (cdr lines))
				      (let* ((rev-bits (string-split (list-ref tar-bits 4) #\.))
					     (n (string->number (list-ref rev-bits 0))))
					(if (not (number? n))
					    (loop best-revision best-date best-file-name (cdr lines))
					    (if (not (> n best-revision))
						(loop best-revision best-date best-file-name (cdr lines))

						;; OK we can update best-revision.  What is the date file this file?
						;; 
						(let ((date-size-str (list-ref bits (- (length bits) 1))))
						  (let ((date-bits (string-split date-size-str #\space)))
						    (let ((this-date (list-ref date-bits 1)))
						      (loop n this-date potential-tar-string (cdr lines)))))))))))))))))))))))
  

  

;; return a sxml paragraph
;; 
(define (source-code-details)

  ;; 
  (define (get-running-status-contents)
    (let* ((running-status-file (append-dir-file devel-dir "source-build-running"))
	   (s (vector-ref (get-url running-status-file) 4)))
      (if (string-match "Not Found" s)
	 #f
	 s)))
    
  
  ;; main line of source-code-details
  `(p 
    ,(receive (source-code-url source-code-file-name revision-count source-tar-date)
 	      (source-code-url-info)
	      `("Source code: "
		(a (@ href ,source-code-url) ,source-code-file-name)
		" "
		;; (a (@ href "latest-coot-build.log") "source build log")
		" "
		,(let ((s (get-running-status-contents)))
		   (if (string? s)
		       (colourize s)
		       ""))
		" "
		,revision-count
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
  (receive (tar-file revision-count)
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
						    "build-status"))
		  (test-status-link (string-append build-log-prefix 
						   "test-status")))

	     (format #t "debug:: build-log-prefix: ~s~%"  build-log-prefix)
	     (format #t "debug:: test-status-link: ~s~%"  test-status-link)
	     (format #t "debug:: build-status-link: ~s~%" build-status-link)
	     ;; (format #t "---- in fill-record debug:: tar-file ~s revision-count: ~s~%" tar-file revision-count)

	     (let* ((build-status-pre (vector-ref (get-url build-status-link) 4))
		    ( test-status-pre (vector-ref (get-url  test-status-link) 4))
		    (build-status (if (< (string-length build-status-pre) 100) build-status-pre "Build-Status-Not-Found"))
		    ( test-status (if (< (string-length  test-status-pre) 100)  test-status-pre "Test-Status-Not-Found")))

  
	       ;; (format #t "============== test-status  for ~s is ~s~%" test-status-link test-status)
	       ;; (format #t "============== build-status for ~s is ~s~%" build-status-link build-status)
		     
	       (if (not (and (string? tar-file)
			     (number? revision-count)))

		   ;; make a dummy record then
		   (begin
		     (format #t "a dummy record created with build-info :~s: ~%" build-info)
		     (apply make-build-record (append build-info (list #f #f build-status test-status))))
		   
		   
		   (let ((binary-url (string-append 
				      (list-ref build-info 2) ;; "/" is included in pre-release dir, right?
				      "/" ;; not currently, it isn't.
				      tar-file)))
		     ;; build-status is e.g. "passed build" or #f
		     ;;  test-status is e.g. "passed test" or #f
		     ;; 
		     (apply make-build-record (append build-info
						      (list revision-count
							    binary-url
							    build-status
							    test-status)))))))))
      


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
    (receive (source-code-url source-code-file-name source-code-revision-count source-code-date)
	     (source-code-url-info)    
	     (write-sxml
	      `(html (head (title "Coot Build Summary Page")
			   (meta (@ (http-equiv refresh) (content 600))))
		     (link (@ (rel "icon")
			      (type "image/png")
			      (href "../web/coot-favicon.png")))
		     (body 
		      (h2 "Coot Pre-Release Build Summary")
		      (p ("Generated " ,(strftime "%a %d %b %H:%M:%S %G %Z" (localtime (current-time)))))
		      "\n"
		      ,(git-revision-count-details)
		      ,(source-code-details)
		      ;; ,(burn-up-chart)  ;; not at the moment.
		      (h3 "Latest Binary Tars")
		      
                      (p "Note: for Binay Tars, the colour merely denote age - not problems")
		      (br)

		      ;;we want the souce code revision number so that we can
		      ;;colour the binary revision numbers in the binaries
		      ;;table.
		      ,(make-binary-table-from-records records source-code-revision-count))
		     
		     (h3 "Development Build Logs")
		     ,(make-builds-table-from-records records source-code-revision-count))
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
(let* ( 
       ;; this puts the file on our local directory 
       ;;
       (html-file-name (string-append (getenv "HOME") 
				      "/public_html/coot/devel/build-info.html"))
       (h (getenv "HOST"))
       (this-host (if (string? h) h "pc"))
       (build-list
	(list 
	 
         ; Gone (it was an old thing)
	 ; (list "binary-Linux-x86_64-centos-5-python-gtk2"
	       ;"http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-pcterm37.lmb.internal/" 
	       ;"http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       ;#t #t)

	 (list "binary-Linux-x86_64-rhel-6-python-gtk2"
	       (string-append "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-" this-host "/")
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       #t #t)

	 (list "binary-Linux-x86_64-scientific-linux-7.5-python-gtk2"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-hal.lmb.internal/"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       #t #t)

; this VM can no longer build coot - make a new Ubuntu image
;	 (list "binary-Linux-x86_64-ubuntu-12.04.3-python-gtk2"
;	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-ubuntu-server/"
;	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
;	       #t #t)

	 (list "binary-Linux-x86_64-ubuntu-14.04-python-gtk2"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-emsley-vm-ubuntu1404/"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       #t #t)

	 (list "binary-Linux-x86_64-ubuntu-18.04-python-gtk2"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-gough-pc-2"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       #t #t)

	 (list "binary-Linux-x86_64-openSUSE-12.3-python-gtk2"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/Linux-emsley-vm-suse/"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       #t #t)

	 (list "binary-Linux-x86_64-debian-gnu-linux-8-python-gtk2"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/build-logs/build-logs/Linux-emsley-vm-debian-7-7/"
	       "http://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/pre-releases"
	       #t #t)

; this machine got upgraded.
;          (list "binary-Linux-x86_64-debian-gnu-linux-6.0-python-gtk2"
;                "http://shelx.uni-ac.gwdg.de/~tg/coot_deb/squeeze/build-logs/Linux-shelx10/gtk2"
;                "http://shelx.uni-ac.gwdg.de/~tg/coot_deb/squeeze/binaries/nightlies/pre-release/"
;                #t #t)

          (list "binary-Linux-x86_64-debian-gnu-linux-jessie-python-gtk2"
                "http://shelx.uni-ac.gwdg.de/~tg/coot_deb/stretch/build-logs/Linux-hilbert/"
                "http://shelx.uni-ac.gwdg.de/~tg/coot_deb/stretch/binaries/nightlies/pre-release"
                #t #t)


; gone - perhaps temporarily?
;           (list "binary-Linux-x86_64-debian-gnu-linux-7-python-gtk2"
;                 "http://shelx.uni-ac.gwdg.de/~tg/coot_deb/stretch/build-logs/Linux-shelx11/"
;                 "http://shelx.uni-ac.gwdg.de/~tg/coot_deb/stretch/binaries/nightlies/pre-release"
;                 #t #t)

	 (list "WinCoot" 
	       "http://www.ysbl.york.ac.uk/~lohkamp/build-logs/MINGW32_NT-6.1-bernie-pc/gtk2" 
	       "http://www.ysbl.york.ac.uk/~lohkamp/software/binaries/nightlies/pre-release"
	       #t #t)

	 )))

  (make-page build-list html-file-name))
	
