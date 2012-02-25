

;; Define an environment variable for the place where jligand.jar
;; resides in a canonical distribution
;; 

(define *to-jligand-secret-file-name* ".coot-to-jligand-8lcs")
(define *from-jligand-secret-link-file-name* ".jligand-to-coot-link")

(define *jligand-home* (let ((s (getenv "JLIGAND_HOME")))
			 (if s s ".")))

(define *java-command* "java")
(define *jligand-jar* (append-dir-file *jligand-home* "jligand.jar"))
(define *jligand-args* (list "-jar" *jligand-jar* "-coot"))

;; jligand internal parameter
(define *imol-jligand-link* #f)

(define (write-file-for-jligand res-spec-1 resname-1 res-spec-2 resname-2)

  (call-with-output-file *to-jligand-secret-file-name*
    (let ((int-time (car (gettimeofday))))
      (format #t "res-spec-1: ~s ~%" res-spec-1)
      (format #t "res-spec-2: ~s ~%" res-spec-2)
      ;; (format #t "res-spec->resno res-spec-1: ~s~%" 
      (lambda (port)
	(let ((chain-id-1 (res-spec->chain-id res-spec-1))
	      (chain-id-2 (res-spec->chain-id res-spec-2)))
	  (display "CODE " port)
	  (display resname-1 port)
	  (display " " port)
	  (display chain-id-1 port)
	  (display " " port)
	  (display (res-spec->res-no res-spec-1) port)
	  (newline port)
	  (display "CODE " port)
	  (display resname-2 port)
	  (display " " port)
	  (display chain-id-2 port)
	  (display " " port)
	  (display (res-spec->res-no res-spec-2) port)
	  (newline port)
	  (display "TIME " port)
	  (write int-time port)
	  (newline port))))))

(define (click->res-spec click)
  (list (list-ref click 2)
	(list-ref click 3)
	(list-ref click 4)))

(if (defined? 'coot-main-menubar)

    (let* ((menu (coot-menubar-menu "JLigand")))

      (add-simple-coot-menu-menuitem
       menu "Launch JLigand"
       (lambda ()
	 (if (not (file-exists? *jligand-jar*))
	     (let ((s (string-append "jligand java jar file: " *jligand-jar* " not found")))
	       (info-dialog s)))
	 (let ((s (string-append-with-spaces (cons *java-command* (append *jligand-args* (list "&"))))))
	   (let ((status (system s)))
	     (print-var status)))))

      (add-simple-coot-menu-menuitem
       menu "Send Link to JLigand (click 2 monomers)"
       (lambda ()
	   (user-defined-click 2 (lambda (clicks)
				   (format #t "we received these clicks: ~s~%" clicks)
				   (if (= (length clicks) 2)
				       (let ((click-1 (list-ref clicks 0))
					     (click-2 (list-ref clicks 1)))
					 (format #t "click-1: ~s~%" click-1)
					 (format #t "click-2: ~s~%" click-2)
					 (if (and (= (length click-1) 7)
						  (= (length click-2) 7))
					     (let ((resname-1 (residue-name 
							       (list-ref click-1 1)
							       (list-ref click-1 2)
							       (list-ref click-1 3)
							       (list-ref click-1 4)))
						   (resname-2 (residue-name 
							       (list-ref click-2 1)
							       (list-ref click-2 2)
							       (list-ref click-2 3)
							       (list-ref click-2 4)))
						   (imol-click-1 (list-ref click-1 1))
						   (imol-click-2 (list-ref click-2 1)))
					       (if (not (and (string? resname-1)
							     (string? resname-2)))
						   (begin 
						     (format #t "Bad resnames: ~s and ~s~%"
							     resname-1 resname-2))
						   (begin
						     (if (not (= imol-click-1 imol-click-2))
							 (begin
							   (set! *imol-jligand-link* #f))
							 (begin ;; happy path
							   (set! *imol-jligand-link* imol-click-1)
							   (write-file-for-jligand (click->res-spec click-1) resname-1 
										   (click->res-spec click-2) resname-2)))))))))))))))

(define (handle-read-from-jligand-file)

  (define (bond-length pos-1 pos-2)
    (define (square x)
      (* x x))
    (sqrt (apply + (map square (map - pos-1 pos-2)))))

  (define (bond-length-from-atoms atom-1 atom-2)
    (if (not (list? atom-1))
	(format #t "   WARNING:: bond-length-from-atoms: atom-1 not a list: ~s~%" atom-1)
	(if (not (list? atom-2))
	    (format #t "   WARNING:: bond-length-from-atoms: atom-2 not a list: ~s~%" atom-2)
	    (bond-length (list-ref atom-1 2)
			 (list-ref atom-2 2)))))

  (define (get-dist atom-spec-1 atom-spec-2)
    (let ((atom-1 (apply get-atom (cons *imol-jligand-link* atom-spec-1)))
	  (atom-2 (apply get-atom (cons *imol-jligand-link* atom-spec-2))))
      (if (not (and (list? atom-1)
		    (list? atom-2)))
	  #f
	  (bond-length-from-atoms atom-1 atom-2))))
	       

  (if (file-exists? *from-jligand-secret-link-file-name*) 
      ;; this file consists of 2 lines:
      ;; the first is a file name that contains the cif dictionary for the cif
      ;; the second is a LINK link that should be added to the PDB file.
      ;; 
      (call-with-input-file *from-jligand-secret-link-file-name*
	(lambda (port)
	  (format #t "Here 2 ~%")
	  (let ((cif-dictionary (read-line port)))
	    (if (string? cif-dictionary)

		(begin
		  ;; was it the READY marker or a cif file name?
		  (if (string=? cif-dictionary "READY")
		      (begin 
			(format #t "JLigand is ready~%")
			;; maybe I need to set something here? (that
			;; from now a modification of .jligand-to-coot
			;; is means that we should read it?)
			))

		  ;; was it a link cif file?
		  (if (file-exists? cif-dictionary)
		      (begin
			(read-cif-dictionary cif-dictionary)
			(let ((link-line (read-line port)))
			  (format #t "Now handle this link line: ~s~%" link-line)
			  (if (> (string-length link-line) 72)
			      (let* ((atom-name-1 (substring link-line 12 16))
				     (atom-name-2 (substring link-line 42 46))
				     (chain-id-1  (substring link-line 21 22)) ;; 1-char, fixme in future
				     (chain-id-2  (substring link-line 51 52)) 
				     (res-no-1-str (substring link-line 22 26))
				     (res-no-2-str (substring link-line 52 56))
				     (res-no-1  (string->number (strip-leading-spaces res-no-1-str)))
				     (res-no-2  (string->number (strip-leading-spaces res-no-2-str)))
				     (ins-code-1 "") ;; too lazy to look these up.
				     (ins-code-2 "")
				     (alt-conf-1 "") ;; ditto
				     (alt-conf-2 "")
				     (link-type    (substring link-line 72)))
				
				(format #t "we parsed  these: \n")
				(format #t "       atom-name-1 ~s~%" atom-name-1)
				(format #t "       atom-name-2 ~s~%" atom-name-2)
				(format #t "        chain-id-1 ~s~%" chain-id-1)
				(format #t "        chain-id-2 ~s~%" chain-id-2)
				(format #t "      res-no-1-str ~s~%" res-no-1-str)
				(format #t "      res-no-1-str ~s~%" res-no-2-str)
				(format #t "        link-type  ~s~%" link-type)

				(if (number? res-no-1)
				    (if (number? res-no-2)
					(if (valid-model-molecule? *imol-jligand-link*)
					    (let* ((res-spec-1 (list chain-id-1 res-no-1 ins-code-1))
						   (res-spec-2 (list chain-id-2 res-no-2 ins-code-2))
						   (atom-spec-1 (append res-spec-1 (list atom-name-1 alt-conf-1)))
						   (atom-spec-2 (append res-spec-2 (list atom-name-2 alt-conf-2)))
						   (dist (get-dist atom-spec-1 atom-spec-2)))
					      (if (not (number? dist))
						  (begin 
						    (format #t "bad dist from ~s ~s~%" atom-spec-1 atom-spec-2))
						  (begin
						    (make-link *imol-jligand-link* atom-spec-1 atom-spec-2 link-type dist)))))))))))))))))))
			  
			  



(define (get-file-mtime file-name)
  (if (not (file-exists? file-name))
      #f
      (stat:mtime (stat file-name))))

(let ((startup-mtime (get-file-mtime *from-jligand-secret-link-file-name*)))
  (gtk-timeout-add 500 (lambda ()
			 (let ((now-time (get-file-mtime *from-jligand-secret-link-file-name*)))
			   ;; (format  #t "startup-mtime: ~s   now-time: ~s~%" startup-mtime now-time)
			   (if (number? now-time)
			       (if (not (number? startup-mtime))
				   (begin
				     (set! startup-mtime now-time)
				     (handle-read-from-jligand-file))
				   (if (> now-time startup-mtime) 
				       (begin
					 (set! startup-mtime now-time)
					 (format #t "just set startup-mtime to ~s~%" startup-mtime)
					 (handle-read-from-jligand-file)))))
			   #t))))

				   
