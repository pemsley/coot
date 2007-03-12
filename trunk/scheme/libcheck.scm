;;;; libcheck.scm
;;;;
;;;; Copyright 2004, 2006 by Paul Emsley, The University of York
;;;; 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 2, or (at your option)
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
	
  ;; main body
  (if (not (or (string? code) (symbol? code)))
	  
	(begin
	  (format #t "WARNING:: Oops code was not a string~%")
	  -2)

	(begin
	  (let* ((dir-prefix 
		  (cond
		   ((null? ccp4i-project-dir) "")
		   (else 
		    (string-append 
		     (car ccp4i-project-dir) "/"))))
		 (code-str (if (string? code)
			       code
			       (symbol->string code)))
		 (libcheck-input 
		  (if (= (string-length dict-cif-libin) 0)
		      (list 
		       "N"
		       (string-append "MON " code-str)
		       "")
		      (list 
		       "N"
		       (string-append "FILE_CIF " dict-cif-libin)
		       (string-append "MON " code-str)
		       "")))
		 (pdb-file-name (string-append dir-prefix 
					       "libcheck_" 
					       code-str ".pdb"))
		 (cif-file-name (string-append dir-prefix
					       "libcheck_" 
					       code-str ".cif"))
		 (post-refmac-pdb-file-name (string-append dir-prefix 
							   "monomer-"
							   code-str ".pdb"))
		 (log-file-name (string-append dir-prefix "coot-libcheck-"
				 code-str ".log"))
		 (refmac-input (list "MODE NEWENTRY" "END"))

		 (refmac-log-file-name  (string-append dir-prefix "coot-libcheck-refmac-"
						       code-str ".log"))
		 (refmac-command-line (list "LIBIN" cif-file-name "XYZIN" pdb-file-name 
					    "XYZOUT" post-refmac-pdb-file-name)))


	    (let ((libstatus (goosh-command libcheck-exe '() libcheck-input log-file-name #t)))
	      
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
			  (let ((libcheck-minimal-desc-status (libcheck-minimal? log-file-name))
				(refmac-status (goosh-command refmac-exe
							      refmac-command-line
							      refmac-input
							      refmac-log-file-name #t)))
			    
			    (if (not (number? refmac-status))
				-4 ; refmac fails 
				(if (not (= 0 refmac-status))
				    -4 ; refmac fails elsewise
				    (let ((recentre-status (recentre-on-read-pdb))
					  (nov (set-recentre-on-read-pdb 0))
					  (pdb-status (read-pdb post-refmac-pdb-file-name)))
				      
				      ; move the coords to the centre of the screen 
				      ; (i.e. by rotation-centre - molecule-centre)
				      (if (>= pdb-status 0) 
					  (begin
					    (let ((rc (rotation-centre))
						  (mc (molecule-centre pdb-status)))
					      (apply translate-molecule-by (cons pdb-status (map - rc mc))))))
				      
				      (set-recentre-on-read-pdb recentre-status)
				      (if libcheck-minimal-desc-status
					  (read-cif-dictionary "libcheck.lib"))
				      pdb-status))))))))))))

