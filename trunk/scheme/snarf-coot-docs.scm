
(use-modules (ice-9 popen)
             (ice-9 string-fun)
	     (ice-9 regex))

(use-modules (goosh))

(define is-valid-model-molecule 1) ;; hack to get past problem loading
				   ;; coot-utils on bubbles
(define  drag-intermediate-atom-scm #t)
(define  mark-atom-as-fixed-scm #t)
(define  ncs-chain-ids-scm #t)
(define  ncs-chain-differences-scm #t)
(define  refmac-parameters-scm #t)
(define  water-chain #f)
(define  water-chain-from-shelx-ins #f)


(load "filter.scm")
(load "coot-utils.scm")

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

(let* ((source-files (glob "*.scm" ".")))

  (map (lambda (source-file)
     
     (let* ((doc-file
	     (string-append 
	      (strip-extension source-file) ".texi"))
	    (tmp-file (string-append doc-file ".tmp")))
       
       ;; problems snarfing coot docs?  check that the function-name
       ;; and args are define all on one line.
       (format #t "running goosh-command guile-tools tmp-file: ~s source-file ~s ~%" 
	       tmp-file source-file)
       (goosh-command "guile-tools" 
		      (list "doc-snarf"
			    "--texinfo"
			    "-o"
			    tmp-file
			    source-file)
		      '() "snarf-log" #t)
       (format #t "   done goosh-command guile-tools~%")
       (goosh-command "grep" (list "-v" "^" tmp-file) '() doc-file #f)
       (if (file-exists? tmp-file)
	   (delete-file tmp-file)
	   (begin
	     (format #t "OOPS: Problem snarfing ~s~%" source-file)
	     (format #t "OOPS: ~s does not exist - not deleting~%" tmp-file)))))
   source-files))
  
(define add-to-list-section-texis
  (lambda (texi-list)

    (let f ((ls texi-list)
	    (new-list '()))
      
      (cond 
       ((null? ls) new-list)
       (else 
	(let* ((stub (strip-extension (car ls)))
	       (section-texi-name (string-append stub "-section.texi")))
	  (f (cdr ls)
	     (append (list section-texi-name (car ls)) new-list))))))))

(define delete-section-texi-files
  (lambda ()

    (let ((section-texi-files (glob "*-section.texi" ".")))
      (map delete-file section-texi-files))))
    

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
		    (format #t "processing ~s~%" file)
		    (let* ((stub (strip-extension file))
			   (section-name (string-append stub "-section.texi")))

		      (format menu-port "* ~a::~%" stub)
		      
		      (call-with-output-file section-name
			(lambda (port)
			  
			  (format port "@node    ~a ~%" stub)
			  (format port "@section ~a ~%" stub)
			  (format port "@cindex  ~a ~%" stub)))))
		  non-empty-texis)
	(format menu-port "@end menu~%")))
	
    
    (let ((all-useful-texis (cons "menu-section.texi"
				  (add-to-list-section-texis non-empty-texis))))

      (format #f "all-useful-texis: ~s~%" all-useful-texis)
      (goosh-command "cat" all-useful-texis '() "coot-scheme-functions.texi" #f))))

