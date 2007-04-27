
;;;; Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 2 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


;;; Provide functions to handle shelx fcf file (which I believe to be old-style
;;; cif rather than mmCIF).
;;; 
;;; So, in Coot, if it sees the extention of .fcf then it calls the
;;; fcf handler and the handler does the conversion, writes a new
;;; mmCIF file and then reads it.
;;; 
(define handle-shelx-fcf-file
  (lambda (filename)

    ;; OK, so we have been passed a shelx fcf file.
    ;; We need to run this through the awk filter.
    ;; 
    (let ((output-mmCIF-file-name (string-append filename ".cif")))

      ; a goosh-command that creates output-mmCIF-file-name
      (let* ((awk-prog (string-append-with-string 
			convert-shelx-fcf-to-cif-awk-strings
			"\n"))
	     (args (list awk-prog filename)))

	; cmd args data-list log-file-name screen-output-also?
	(goosh-command "awk" args '() output-mmCIF-file-name #f)
      
	(auto-read-cif-data-with-phases output-mmCIF-file-name)))))


(define convert-shelx-fcf-to-cif-awk-strings
  (list 
"BEGIN { intable = 0; }"
""
"intable == 0 { "
"  if ($1 == \"_refln_F_squared_meas\") "
"    print \" _refln.F_meas\";"
"  else "
"    if ($1 == \"_refln_F_squared_sigma\") "
"      print \" _refln.F_meas_sigma\";"
"    else "
"      "
"      # get rid of shelx title anomolies"
"      if ($1 == \"_shelx_title\"           ||"
"	  $1 == \"_shelx_refln_list_code\" ||"
"	  $1 == \"_shelx_F_calc_maximum\"  ||"
"	  $1 == \"_exptl_crystal_F_000\"   ||"
"	  $1 == \"_reflns_d_resolution_high\") { "
"	# don't print it"
"      } else {"
"	# symm equivs:"
"	if ($1 == \"_symmetry_equiv_pos_as_xyz\") "
"	  print \"_symmetry_equiv.pos_as_xyz\""
"	else "
"	  if ($1 == \"_cell_length_a\")"
"	    print \"_cell.length_a     \", $2;"
"	  else "
"	    if ($1 == \"_cell_length_b\")"
"	      print \"_cell.length_b     \", $2;"
"	    else"
"	      if ($1 == \"_cell_length_c\")"
"		print \"_cell.length_c     \", $2;"
"	      else "
"		if ($1 == \"_cell_angle_alpha\")"
"		  print \"_cell.angle_alpha  \", $2;"
"		else "
"		  if ($1 == \"_cell_angle_beta\")"
"		    print \"_cell.angle_beta   \", $2;"
"		  else "
"		    if ($1 == \"_cell_angle_gamma\")"
"		      print \"_cell.angle_gamma  \", $2;"
"		    else "
"		      if ($1 == \"_refln_index_h\")"
"			print \" _refln.index_h\";"
"		      else "
"			if ($1 == \"_refln_index_k\")"
"			  print \" _refln.index_k\";"
"			else "
"			  if ($1 == \"_refln_index_l\")"
"			    print \" _refln.index_l\";"
"			  else "
"			    if ($1 == \"_refln_F_calc\")"
"			      print \" _refln.F_calc\";"
"			    else "
"			      if ($1 == \"_refln_phase_calc\")"
"				print \" _refln.phase_calc\";"
"			      else "
"				print $0;"
"      }"
"}"
"  "
"/^ _refln_phase_calc/ { "
""
"  if ($1 == \"_refln_phase_calc\") {"
"    intable = 1;"
"    next;"
"  }"
"}"
""
"intable == 1 {"
"  if (NF>5)"
"     if ($4+0 > 0)"
"        print $1, $2, $3, sqrt($4+0), 0.5*($5+0)/sqrt($4+0), $6, $7;"
"}"
""
))

(define *coot-shelxl-dir* "coot-shelxl")

(define remove-time-extension
  (lambda (str)
    
    ;; e.g. remove -2007-04-27_1326.24 from 03srv164-2007-04-27_1326.24

    ;; shelx and coot will retire before the 22nd century, won't they?
    (let ((pattern "-20[0-9][0-9]-[01][0-9]-[0-2][0-9]_[0-2][0-9][0-5][-0-9].[0-5][-0-9]"))
      (let ((result (string-match pattern str)))
	(if result
	    (substring str 0 (car (vector-ref result 1)))
	    str)))))


(define shelxl-refine 
  (lambda (imol . hkl-file-in)

    ; First write out the imol-th molecule in shelx format.
    ; What should the filename be? 
    ; Let's say that we read "abc.res"
    ; We want to create "abc-coot.ins".  Ah, but what about the data file?
    ; Urg.
    ; 
    ; This, as it stands, doesn't properly deal with incrementing the
    ; filenames.

    (if (not (command-in-path? "shelxl"))

	(format #t "coot warning: can't find executable shelxl~%")

	(let ((orig-hkl-file (cond 
			      ((eq? '() hkl-file-in) #f) 
			      (else (car hkl-file-in)))))

	  (if (not (file-exists? *coot-shelxl-dir*))
	      (mkdir *coot-shelxl-dir*))

	  (if (not (file-exists? *coot-shelxl-dir*))
	      ; We failed to make the directory for some reason.
	      ; Exit now.
	      (format #t "Failed to make shelxl directory ~s~% " 
		      *coot-shelxl-dir*)

	      (let* ((stub (string-append 
			    *coot-shelxl-dir* "/"
			    (strip-path (remove-time-extension (strip-extension (molecule-name imol))))
			    "-"
			    (unique-date/time-str)))
		     (ins-filename (string-append stub ".ins"))
		     (res-filename (string-append stub ".res"))
		     (fcf-filename (string-append stub ".fcf"))
		     (log-filename (string-append stub ".log"))
		     (hkl-filename (string-append stub ".hkl")))
		
		; make a link to the passed hkl-file-in, if it was passed.
		; (if not, we presume that this file (or the link) exists
		; already. (Let's check that, and stop early if it doesn't
		; exist)

		;; This is a bit of mind-mangling code.

		(if orig-hkl-file
		    (let ((symlink-target
			   (cond
			    ((= (string-length orig-hkl-file) 0) #f)
			    ((string=? (substring hkl-filename 0 1) "/") orig-hkl-file)
			    (else 
			     (if (slash-start? orig-hkl-file)
				 orig-hkl-file
				 (string-append "../" orig-hkl-file))))))
		      (format #t "make symlink ~s ~s~%" orig-hkl-file symlink-target)
		      (if symlink-target
			  (symlink symlink-target hkl-filename)))
		    
		    ;; hklin was not given, let's generate a filename
		    ;; from the filename of the coordinates molecule
		    ;; imol, trial names are derived from the name of
		    ;; the coordinates molecule
		    (let* ((trial-file-stub (strip-extension (molecule-name imol)))
			   (trial-hkl-file-name (string-append
						 trial-file-stub ".hkl")))
		      (format #t "debug:: looking for ~s~%" trial-hkl-file-name)
		      (if (file-exists? trial-hkl-file-name)
			  (begin 
			    (format #t "INFO:: hkl file ~s found~%" 
				    trial-hkl-file-name)
			    (format #t "need a link: ~s ~s~%" trial-hkl-file-name hkl-filename)
			    ;; if trial-hkl-file-name and hkl-filename
			    ;; are in the same directory (e.g
			    ;; shelx-coot), then (strange as it may
			    ;; seem) we *don't* put the directory name
			    ;; in the "file that exists" (the first
			    ;; argument - ie. the target).  Bizarre but true.
			    ;; (symlink target link)
			    (let ((hkl-filename-dir (dirname hkl-filename))
				  (trial-hkl-file-name-dir (dirname trial-hkl-file-name)))
			      (format #t "debug: comparing paths ~s ~s~%"
				      hkl-filename-dir
				      trial-hkl-file-name-dir)
				   
			      (if (string=? hkl-filename-dir
					    trial-hkl-file-name-dir)
				  ;; so lets strip off the dir of the target:
				  (let ((new-target (strip-path trial-hkl-file-name)))
				    (format #t "symlink: ~s ~s~%" new-target hkl-filename)
				    (symlink new-target hkl-filename))
				  ;; different dirs, keep the dirs:
				  (begin
				    (format #t "symlink2: ~s ~s~%" trial-hkl-file-name 
					    hkl-filename)
				    (symlink trial-hkl-file-name hkl-filename)))))
			  (begin 
			    (format #t "INFO:: hkl file ~s does not exist!~%" 
				    trial-hkl-file-name)))))
		

		(if (not (file-exists? hkl-filename))
		    (begin
		      (format #t "data (hkl) file ~s not found - not running shelxl~%"
			      hkl-filename))
		    
		    (begin
		      (write-shelx-ins-file imol ins-filename)
		      ; running shelxl creates stub.res
		      (goosh-command "shelxl" (list stub) '() log-filename #t)
		      ; it isn't a pdb file, but coot knows what to do.
		      (read-pdb res-filename)
		      (handle-shelx-fcf-file fcf-filename)))))))))

    
; (shelxl-refine 0 "data/shelx/problems/insulin.hkl")

;; do-shelx-lst-file
;;
;; ie. create a interesting-things GUI for split (and other things?)
;; in a shelx .lst file.
;; 
(define read-shelx-lst-file
  (lambda (file-name imol)

    ;; turn interesting-list into a GUI:
    ;; 
    (define gui-interesting-list 
      (lambda (interesting-list)
	
	(format #t "debug: interesting-list: ~s~%" interesting-list)
	(if (not (null? interesting-list))
	    (interesting-things-gui "Interesting Things from SHELX"
				    interesting-list))))

    ;; chop off last char of string
    (define (chop-end s)
      (let ((l (string-length s)))
	(if (= s 0)
	    ""
	    (substring s 0 (- l 1)))))

    ;; return last char of string:
    (define (last-char s)
      (let ((l (string-length s)))
	(substring s (- l 1) l)))

    ;; return a  chain resno inscode atom-name altconf list
    ;; 
    (define make-atom-parts
      (lambda (atom-id)

	(let ((ls (split-discarding-char-last #\_ atom-id list)))
	  (if (= (length ls) 2)

	      ;; chain resno inscode atom-name altconf
	      ;; 
	      ;; The resno part can contain an alt conf specifier
	      ;; e.g.: 34 or 34a. So let's try to make it a number
	      ;; simply.  If that fails, then strip off the last char
	      ;; and try again.  If *that* is successful, return the
	      ;; resno and the last char is the altconf specifier.
	      ;; 
	      (if (number? (string->number (car (cdr ls))))
		  (list "" (string->number (car (cdr ls))) "" (car ls) "")
		  ;; handle the alt conf
		  (let* ((p (chop-end s))
			 (alt-conf (last-char s))
			 (n (string->number p)))

		    (if (number? n)
			(list "" n "" (car ls) alt-conf)
			(list "" 1 "" "blank" "")))) ; failure
		    
	      (list "" 1 "" "blank" ""))))) ; failure
	  

    ;; main body
    (if (not (valid-model-molecule? imol))
	(format #t "WARNING:: Molecule number ~s not valid~%" imol)
	(if (not (file-exists? file-name))
	    (format #t "WARNING:: shelx lst file: ~s does not exist~%" 
		    file-name)
	    (call-with-input-file file-name
	      (lambda (port)
		
		(let loop ((line (read-line port))
			   (interesting-list '()))
		  
		  (cond 
		   ((eof-object? line) (gui-interesting-list interesting-list))
		   ((string-match "may be split into" line)
		    (let ((parts (string->list-of-strings line)))
		      (if (> (length parts) 6)
			  (let* ((atom-id (list-ref parts 3)))
			    ;; e.g; atom-id: CD1_1392
			    (let ((atom-parts (make-atom-parts atom-id)))
			      (let ((buton-label (string-append
						  "Atom "
						  (number->string (car (cdr atom-parts)))
						  " "
						  (list-ref atom-parts 3) ; atom name
						  " may be split?")))
				(loop (read-line port)
				      (cons 
				       (append
					(list buton-label imol) atom-parts)
				       interesting-list))))))))
		   (else 
		    (loop (read-line port)
			  interesting-list))))))))))


				    



	
