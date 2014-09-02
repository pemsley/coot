
(use-modules (ice-9 popen)
             (ice-9 string-fun)
	     (ice-9 rdelim)
	     (ice-9 regex))

(use-modules (os process))

(define is-valid-model-molecule 1) ;; hack to get past problem loading
				   ;; coot-utils on bubbles


(load "filter.scm")
;; (load "coot-utils.scm")

(define is-empty-file?
  (lambda (file-name)

    (let ((strs (run-command/strings "wc" (list "-l" file-name) '())))
      (if (list? strs)
	  (let ((ls (string->list-of-strings (car strs))))
	    (if (not (null? ls))
		(let ((n (string->number (car ls))))
		  (if (number? n)
		      (= 0 n)
		      #f))
		#f))
	  #f))))


;; return a list of file names that match pattern pat in directory dir.
(define (glob pat dir)
  (let ((rx (make-regexp (glob->regexp pat))))
    (filter (lambda (x) (regexp-exec rx x)) (directory-files dir))))


(define (glob->regexp pat)
  (let ((len (string-length pat))
        (ls '("^"))
        (in-brace? #f))
    (do ((i 0 (1+ i)))
        ((= i len))
      (let ((char (string-ref pat i)))
        (case char
;           ((#\*) (set! ls (cons "[^.]*" ls)))
          ((#\*) (set! ls (cons ".*" ls)))
          ((#\?) (set! ls (cons "[^.]" ls)))
          ((#\[) (set! ls (cons "[" ls)))
          ((#\]) (set! ls (cons "]" ls)))
          ((#\\)
           (set! i (1+ i))
           (set! ls (cons (make-string 1 (string-ref pat i)) ls))
           (set! ls (cons "\\" ls)))
          (else
           (set! ls (cons (regexp-quote (make-string 1 char)) ls))))))
    (string-concatenate (reverse (cons "$" ls)))))

;; "a.b.res" -> "a.b"
;; file-name-sans-extension
;; 
(define (strip-extension s) 
  (let ((ls (split-after-char-last #\. s list)))
    
    (cond 
     ((string=? (car ls) "") (car (cdr ls)))
     (else 
      (let ((dotted-s (car ls)))
	(car (split-discarding-char-last #\. (car ls) list)))))))


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
      
      (begin 
	(format #t "command ~s not found~%" cmd)
	255)

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


;; Return the strings screen output of cmd or #f if command was not found
;; 
(define (run-command/strings cmd args data-list)
    
  (if (not (command-in-path? cmd))

      (begin
	(format #t "command ~s not found~%" cmd)
	#f)

      (begin 
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
	  
	  (let f ((obj (read-line output-port))
		  (ls '()))
	    
	    (if (eof-object? obj)
		(begin 
		  (let* ((status-info (waitpid pid))
			 (status (status:exit-val (cdr status-info))))
		    (reverse ls))) ; return ls (in the right order)
		
		(begin
		  (f (read-line output-port) (cons obj ls)))))))))

;; Return a list if str is a string, else return '()
;; 
(define (string->list-of-strings str)

  (if (not (string? str))
      '()
      (let f ((chars (string->list str))
	      (word-list '())
	      (running-list '()))

	(cond 
	 ((null? chars) (reverse 
			 (if (null? running-list)
			     word-list
			     (cons (list->string (reverse running-list)) word-list))))
	 ((char=? (car chars) #\space) (f (cdr chars)
					  (if (null? running-list)
					      word-list
					      (cons (list->string 
						     (reverse running-list))
						    word-list))
					  '()))
	 (else 
	  (f (cdr chars)
	     word-list
	     (cons (car chars) running-list)))))))

;; The following functions from PLEAC (guile version thereof of course).
;; 
;; or define a utility function for this
(define (directory-files dir)

  (if (not (access? dir R_OK))
    '()
    (let ((p (opendir dir)))
      (do ((file (readdir p) (readdir p))
           (ls '()))
          ((eof-object? file) (closedir p) (reverse! ls))
        (set! ls (cons file ls))))))

;; texi2html on user-manual.texi with coot-scheme-functions.texi included gives:
;; 
;; *** Duplicate node found: mutate (in ./../../coot/scheme/coot-scheme-functions.texi l. 42)
;; *** Duplicate node found: povray (in ./../../coot/scheme/coot-scheme-functions.texi l. 341)
;;*** Duplicate node found: raster3d (in ./../../coot/scheme/coot-scheme-functions.texi l. 361)
;; 
(let* ((source-files (glob "*.scm" ".")))

  (map (lambda (source-file)
     
     (let* ((doc-file
	     (string-append 
	      (strip-extension source-file) ".texi"))
	    (tmp-file-1 (string-append doc-file ".tmp")))
       
       ;; problems snarfing coot docs?  check that the function-name
       ;; and args are define all on one line.
       (let ((args (list "doc-snarf"
			    "--texinfo"
			    "-o"
			    tmp-file-1
			    source-file)))
	 (format #t "guile-tools ~s > snarf.log~%" args)
	 (let ((status (goosh-command "guile-tools" args '() "snarf-log" #t)))
	   (format #t "guile-tools ~s returns status ~s~%" args status)
	   (if (not (= status 0))
	       (exit status)))
	 
	 (goosh-command "grep" (list "-v" "^" tmp-file-1) '() doc-file #f))

       (if (file-exists? tmp-file-1)
	   (begin
	     (delete-file tmp-file-1))
	   (begin
	     (format #t "OOPS: Problem snarfing ~s~%" source-file)
	     (format #t "OOPS: ~s does not exist - not deleting~%" tmp-file-1)))))

       source-files))
  
;; 
;; 
(define add-to-list-section-texis
  (lambda (texi-list)

    (let f ((ls texi-list)
	    (new-list '()))
      
      (cond 
       ((null? ls) (reverse new-list))
       (else 
	(let* ((stub (strip-extension (car ls)))
	       (section-texi-name (string-append stub "-section.texi")))
	  (f (cdr ls)
             ;; the section texi and the tex go in "backwards" because the list is 
             ;; reversed on return
	     (append (list (car ls) section-texi-name) new-list))))))))

;; 
;; 
(define delete-section-texi-files
  (lambda ()

    (let ((section-texi-files (glob "*-section.texi" ".")))
      (map delete-file section-texi-files))))
    

;; 
;; 
(delete-section-texi-files)
(let ((all-raw-texi-files (glob "*.texi" ".")))
  
  (let ((non-empty-texis
	 (let f ((ls all-raw-texi-files)
		 (non-null-list '()))
	   
	   (cond 
	    ((null? ls) non-null-list)
	    ;; we also want to reject coot-scheme.texi, the master
	    ;; texi file (we don't want to cat that) and the catted file
	    ((string=? (car ls) "coot-scheme.texi")
	     (f (cdr ls) non-null-list))
	    ((string=? (car ls) "coot-scheme-functions.texi")
	     (f (cdr ls) non-null-list))
	    ((is-empty-file? (car ls)) 
	     (delete-file (car ls))
	     (f (cdr ls) non-null-list))
	    (else 
	     (f (cdr ls) (cons (car ls) non-null-list)))))))

    ;; we want to add section/node info before each of the texi files
    ;; so let's create a section-texi file that has that info in.
    ;;
    ;; make the section texis here
    (call-with-output-file "menu-section.texi"
      (lambda (menu-port)
	(format menu-port "@menu~%")
	(for-each (lambda (file)
		    (format #t "menu section processing ~s~%" file)
		    (let* ((stub (strip-extension file))
			   (node-name (cond
				  ((string=? stub "mutate") "mutate-in-scheme")
				  ((string=? stub "povray") "mutate-from-scheme")
				  ((string=? stub "raster3d") "raster3d-from-scheme")
				  (else
				   stub)))
			   (section-name (string-append stub "-section.texi")))

		      (format menu-port "* ~a::~%" node-name)
		      
		      (call-with-output-file section-name
			(lambda (port)
			  
			  (format port "@node    ~a ~%" node-name)
			  (format port "@section ~a ~%" node-name)
			  ;; (format port "@cindex  ~a ~%" stub) not cindex
			  ))))
		  non-empty-texis)
	(format menu-port "@end menu~%")))
	
    
    (let ((all-useful-texis (cons "menu-section.texi"
				  (add-to-list-section-texis non-empty-texis))))

      (format #t "all-useful-texis: ~s~%" all-useful-texis)
      (goosh-command "cat" all-useful-texis '() "coot-scheme-functions.texi" #f))))


