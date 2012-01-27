

;; Define an environment variable for the place where jligand.jar
;; resides in a canonical distribution
;; 

(define *to-jligand-secret-file-name* ".coot-to-jligand-8lcs")
(define *from-jligand-secret-file-name* ".jligand-to-coot")

(define *jligand-home* (let ((s (getenv "JLIGAND_HOME")))
			(if s s ".")))

(define *java-command* "java")
(define *jligand-jar* (append-dir-file *jligand-home* "jligand.jar"))
(define *jligand-args* (list "-jar" *jligand-jar* "-coot"))

(define (write-file-for-jligand resname-1 resname-2)

  (call-with-output-file *to-jligand-secret-file-name*
    (let ((int-time (car (gettimeofday))))
      (lambda (port)
	(display "CODE " port)
	(display resname-1 port)
	(newline port)
	(display "CODE " port)
	(display resname-2 port)
	(newline port)
	(display "TIME " port)
	(write int-time port)
	(newline port)))))


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
				   (format #t "we recieved these clicks: ~s~%" clicks)
				   (if (= (length clicks) 2)
				       (let ((click-1 (list-ref clicks 0))
					     (click-2 (list-ref clicks 1)))
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
							       (list-ref click-2 4))))
					       (if (not (and (string? resname-1)
							     (string? resname-2)))
						   (format #t "Bad resnames: ~s and ~s~%"
							   resname-1 resname-2)
						   (write-file-for-jligand resname-1 resname-2))))))))))))

(define (handle-read-from-jligand-file)
  (if (file-exists? *from-jligand-secret-file-name*)
      (read-cif-dictionary *from-jligand-secret-file-name*)))


(define (get-file-mtime file-name)
  (if (not (file-exists? file-name))
      #f
      (stat:mtime (stat file-name))))

(let ((startup-mtime (get-file-mtime *from-jligand-secret-file-name*)))
  (gtk-timeout-add 500 (lambda ()
			 (let ((now-time (get-file-mtime *from-jligand-secret-file-name*)))
			   (format  #t "startup-mtime: ~s   now-time: ~s~%" startup-mtime now-time)
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

				   
