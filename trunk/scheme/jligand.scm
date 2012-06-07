

;; Coot and JLigand communicate via "secret" files.
;; 
;; Coot send to jligand the comp-ids it wants joined in *to-jligand-secret-file-name*.
;;
;; JLigand sends Coot a file that contains filenames for the link and
;; new dictionary/restraints.
;; 
(define *to-jligand-secret-file-name* ".coot-to-jligand-8lcs")
(define *from-jligand-secret-link-file-name* ".jligand-to-coot-link")

;; Define an environment variable for the place where jligand.jar
;; resides in a canonical distribution
;; 
(define *jligand-home-env* (getenv "JLIGAND_HOME"))
(define *jligand-home* (if *jligand-home-env* *jligand-home-env* "."))

(define *java-command* "java")  ;; what else would it be?

(define *jligand-jar* (append-dir-file *jligand-home* "jligand.jar"))
(define *jligand-args* (list "-jar" *jligand-jar* "-coot"))  ;; for coot-enable jligand

;; jligand internal parameter
(define *imol-jligand-link* #f)

(define (jligand-code-file-maybe comp-id port)

  (define (jligand-standard-amino-acid? comp-id)
    (or (string=? comp-id "ALA")
	(string=? comp-id "ARG")
	(string=? comp-id "ASN")
	(string=? comp-id "ASP")
	(string=? comp-id "CYS")
	(string=? comp-id "GLY")
	(string=? comp-id "GLU")
	(string=? comp-id "GLN")
	(string=? comp-id "PHE")
	(string=? comp-id "HIS")
	(string=? comp-id "ILE")
	(string=? comp-id "LEU")
	(string=? comp-id "LYS")
	(string=? comp-id "MET")
	(string=? comp-id "PRO")
	(string=? comp-id "SER")
	(string=? comp-id "TYR")
	(string=? comp-id "THR")
	(string=? comp-id "VAL")
	(string=? comp-id "TRP")))

  ;; if code was read from a cif file then write a line to port like
  ;; "CODE TLC FILE file.cif" (with a newline).
  ;; 
  (if (not (jligand-standard-amino-acid? comp-id))
      (let ((cif-file (cif-file-for-comp-id comp-id)))
	(if (> (string-length cif-file) 0)
	    (begin
	      (display "CODE " port)
	      (display comp-id port)
	      (display " " port)
	      (display "FILE " port)
	      (display cif-file port)
	      (newline port))))))


(define (write-file-for-jligand res-spec-1 resname-1 res-spec-2 resname-2)

  (call-with-output-file *to-jligand-secret-file-name*
    (let ((int-time (car (gettimeofday))))
      ;; (format #t "res-spec-1: ~s ~%" res-spec-1)
      ;; (format #t "res-spec-2: ~s ~%" res-spec-2)
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
	  (jligand-code-file-maybe resname-1 port)
	  (display "CODE " port)
	  (display resname-2 port)
	  (display " " port)
	  (display chain-id-2 port)
	  (display " " port)
	  (display (res-spec->res-no res-spec-2) port)
	  (newline port)
	  (jligand-code-file-maybe resname-2 port)
	  (display "TIME " port)
	  (write int-time port)
	  (newline port))))))


;; this could be in utils
(define (get-file-mtime file-name)
  (if (not (file-exists? file-name))
      #f
      (stat:mtime (stat file-name))))

(define (click->res-spec click)
  (list (list-ref click 2)
	(list-ref click 3)
	(list-ref click 4)))
;; every fraction of a second look to see if
;; *from-jligand-secret-link-file-name* has been updated.  If so, then
;; run the handle-read-from-jligand-file function.
;; 
(define (start-jligand-listener) 
  (let ((startup-mtime (get-file-mtime *from-jligand-secret-link-file-name*)))
    (gtk-timeout-add 700 
		     (lambda ()
		       (let ((now-time (get-file-mtime *from-jligand-secret-link-file-name*)))
			 (if (number? now-time)
			     (if (not (number? startup-mtime))
				 (begin
				   (set! startup-mtime now-time)
				   (handle-read-from-jligand-file))
				 (if (> now-time startup-mtime) 
				     (begin
				       (set! startup-mtime now-time)
				       (handle-read-from-jligand-file)))))
			 #t)))))




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
			  
			  



;;; stand-alone/development
;;; 
;(if (defined? 'stand-alone)n
;    (if (defined? 'coot-main-menubar)

;	(let* ((menu (coot-menubar-menu "JLigand")))

;	  (add-simple-coot-menu-menuitem
;	   menu "Launch JLigand"
;	   (lambda ()
;	     (launch-jligand-function))))))




