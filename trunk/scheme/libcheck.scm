;;;; libcheck.scm
;;;;
;;;; Copyright 2004, 2006 The University of York
;;;; Copyright 2009 The University of Oxford
;;;; Author: Paul Emsley
;;;; 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3, or (at your option)
;;;; any later version.
;;;; 
;;;; This program is distributed in the hope that it will be useful,
;;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;;; GNU General Public License for more details.
;;;; 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this software; see the file COPYING.  If not, write to
;;;; the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
;;;; Boston, MA 02111-1307 USA
;;;;

;; this is override-able by the user in their .coot file (for example).
(define libcheck-exe "libcheck")

;;;; This files provides an interface to libcheck

;; Where @var{code} is @emph{e.g.} "3GP" and @var{ccp4i-project-dir}
;; is an optional arg which is the directory to which the libcheck log
;; files and pdb file will go.
;;
;; Return -2 on @var{code} is not a string
;; Return -3 on libcheck failure
;; Return @var{imol} on success
;; Return @code{handle-read-draw-molecule} error code on failure to 
;;       read resultant pdb file
;; 
;; Actually, it would be nice to know if this code represented a full
;; description - or a minimal one... perhaps we can parse the log
;; file and call a pop-up that will tell us.
;; 
;; @var{dict-cif-libin should be a string}.  If it is "" then it is
;; ignored.  If it is not "" then it is used to create input to
;; libcheck (not a command line argument) so that bespoke dictionary
;; libraries can produce coords using "Get Monomer".
;; 
(define (monomer-molecule-from-3-let-code code dict-cif-libin . ccp4i-project-dir)

  ;; return #t if log file has "has the minimal description" in it.
  ;; else return #f
  ;; 
  (define libcheck-minimal? 
    (lambda (filename)
      
      (call-with-input-file filename
	(lambda (port)
	  
	  (let loop ((line (read-line port)))
	    
	    (cond
	     ((eof-object? line) #f)
	     (else 
	      (let ((m (string-match "has the minimal description" line)))
		(if (eq? #f m)
		    (loop (read-line port))
		    #t)))))))))

  ;; move file-name to some file name with a date extension
  ;; 
  (define (move-aside file-name)

    (if (file-exists? file-name) 
	(let ((extension 
	       (let ((lt (localtime (current-time))))
		 
		 (string-append 
		  (number->string (+ 1900 (vector-ref lt 5)))
		  ":"
		  (number->string (vector-ref lt 4))
		  ":"
		  (number->string (vector-ref lt 3))
		  ":"
		  (number->string (vector-ref lt 2))
		  ":"
		  (number->string (vector-ref lt 1))
		  ":"
		  (number->string (vector-ref lt 0))))))

	  (let ((new-file-name (string-append 
				file-name
				"-"
				extension)))
	    (rename-file file-name new-file-name)))))

  ;; the exe and log-file-name (and command-lines-args perhaps)
  ;; relative to dir (not the calling dir)
  ;; 
  (define (run-command-in-dir dir exe-name command-lines-args data-lines log-file-name to-screen-flag)
    (let ((current (getcwd)))
      (chdir dir)
      (let ((status 
	     (goosh-command exe-name command-lines-args data-lines log-file-name to-screen-flag)))
	(chdir current)
	status)))

  ;; This gets called in non-graphics mode too.
  ;; 
  (define (libcheck-monomer-gui dir-prefix code-str cif-file-name 
				pdb-file-name post-refmac-pdb-file-name)

    ;; (format #t "================= debug: libcheck-monomer-gui: dir-prefix: ~s   code-str: ~s   cif-file-name: ~s   pdb-file-name: ~s   post-refmac-pdb-file-name: ~s~%"
    ;; dir-prefix code-str cif-file-name 
    ;; pdb-file-name post-refmac-pdb-file-name)

    (let* ((libcheck-input 
	    (if (= (string-length dict-cif-libin) 0)
		(list 
		 "N"
		 (string-append "MON " 
				;; This doesn't seem to work,
				;; but I'll leave it in anyway.
				;; Libcheck fails with both
				;; "GAL-b-D" and "GAL" fail.
				(if (> (string-length code-str) 3)
				    (substring code-str 0 3)
				    code-str))
		 "")
		(list 
		 "N"
		 (string-append "FILE_CIF " dict-cif-libin)
		 (string-append "MON " code-str)
		 "")))
	   (log-file-name (string-append "coot-libcheck-"
					 code-str ".log"))
	   (log-file-name-in-dir (append-dir-file dir-prefix log-file-name))
	   (refmac-input (list "MODE NEWENTRY" "END"))

	   (refmac-log-file-name  (string-append dir-prefix "coot-libcheck-refmac-"
						 code-str ".log"))
	   (refmac-command-line (list "LIBIN" cif-file-name "XYZIN" pdb-file-name 
				      "XYZOUT" post-refmac-pdb-file-name)))

      (move-aside (append-dir-file dir-prefix "libcheck.lib"))
      (let ((libstatus (run-command-in-dir dir-prefix libcheck-exe '() libcheck-input
					   log-file-name #t)))
	
	(format #t "INFO:: libcheck status: ~s~%" libstatus)
	
	(if (not (number? libstatus))
	    -3
	    (if (= libstatus 0)
		;; libcheck ran OK, 
		;; 
		;; But I now find that libcheck can run OK, but
		;; not produce an output file (using dict .cif
		;; file from PRODRG).
		;;
		;; So we first need to check that the output of
		;; libcheck exists.
		;; 
		(if (not (file-exists? cif-file-name))

		    (format #t "libcheck failed to write the output cif file.~%")

		    ;; OK, now let's run refmac:
		    ;; 
		    (let ((libcheck-minimal-desc-status (libcheck-minimal? log-file-name-in-dir))
			  (refmac-status (goosh-command refmac-exe
							refmac-command-line
							refmac-input
							refmac-log-file-name #t)))

		      (format #t "INFO:: libcheck-minimal? is ~s~%" libcheck-minimal-desc-status)
		      
		      (if (not (number? refmac-status))
			  -4 ; refmac fails 
			  (if (not (= 0 refmac-status))
			      -4 ; refmac fails elsewise
			      (begin 
				;; if there was a minimal description,
				;; we get the real cif file in
				;; libcheck.lib.
				(let ((libcheck-lib (append-dir-file dir-prefix "libcheck.lib")))
				  (if (file-exists? libcheck-lib)
				      (copy-file libcheck-lib cif-file-name)))
				(handle-libcheck-cif-and-pdb cif-file-name
							     pdb-file-name
							     post-refmac-pdb-file-name)))))))))))
  
  (define (handle-libcheck-cif-and-pdb cif-file-name pdb-file-name post-refmac-pdb-file-name)
    

      (if (and (file-exists? post-refmac-pdb-file-name)
	       (file-exists? cif-file-name))
	  (let ((pdb-status (handle-read-draw-molecule-with-recentre 
			     post-refmac-pdb-file-name 0)))
	    (if (valid-model-molecule? pdb-status) 
		(begin
		  (assign-hetatms pdb-status)
		  (move-molecule-here pdb-status)
		  (read-cif-dictionary cif-file-name)
		  pdb-status))))) ; return imol of the ligand

  

  ;; main body
  (if (not (or (string? code) (symbol? code)))
	  
	(begin
	  (format #t "WARNING:: Oops code was not a string~%")
	  -2)

	;; do the files exist already?  If so, just read them in.
	(let* ((dir-prefix 
	       (cond
		((null? ccp4i-project-dir)
		 (let ((dir-name "coot-ccp4"))
		   (make-directory-maybe dir-name)
		   (string-append dir-name "/")))
		(else 
		 (string-append 
		  (car ccp4i-project-dir) "/"))))
	       (code-str (if (string? code)
			     code
			     (symbol->string code)))

	       (pdb-file-name (string-append dir-prefix 
					     "libcheck_" 
					     code-str ".pdb"))
	       (cif-file-name (string-append dir-prefix
					     "libcheck_" 
					     code-str ".cif"))
	       (post-refmac-pdb-file-name (string-append dir-prefix 
							 "monomer-"
							 code-str ".pdb")))

	  (if (and (file-exists? post-refmac-pdb-file-name)
	       (file-exists? cif-file-name))
	      (handle-libcheck-cif-and-pdb cif-file-name pdb-file-name post-refmac-pdb-file-name)
	      
	      (if (not (command-in-path? libcheck-exe))
		  (begin
		    (info-dialog "You need to setup CCP4 (specifically LIBCHECK) first.")
		    -2)
		  (libcheck-monomer-gui dir-prefix code-str cif-file-name pdb-file-name post-refmac-pdb-file-name))))))
  
