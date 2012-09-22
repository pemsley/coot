
;;;; Copyright 2004, 2005, 2006 by The University of York
;;;; Copyright 2007 by Paul Emsley
;;;; Copyright 2007, 2009 by The University of Oxford
 
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


;;; Provide functions to handle shelx fcf file (which I believe to be old-style
;;; cif rather than mmCIF).
;;; 
;;; So, in Coot, if it sees the extention of .fcf then it calls the
;;; fcf handler and the handler does the conversion, writes a new
;;; mmCIF file and then reads it.
;;; 
(define handle-shelx-fcf-file-old
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

(define handle-shelx-fcf-file handle-shelx-fcf-file-internal)


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

    ;; shelx and coot will retire before the 22nd century, won't they? (thank you Gabriele Balducci)
    (let ((pattern "-20[0-9][0-9]-[01][0-9]-[0-3][0-9]_[0-2][0-9][0-5][-0-9].[0-5][-0-9]"))
      (let ((result (string-match pattern str)))
	(if result
	    (substring str 0 (car (vector-ref result 1)))
	    str)))))

;; 
(define shelxl-refine 
  (lambda (imol . hkl-file-in)

    (let ((func (lambda (ins-file-name)
		  (write-shelx-ins-file imol ins-file-name))))
      ;; hkl-file-in can be null '() or (list astring).
      (shelxl-refine-inner imol hkl-file-in func))))


;; hkl-file-in can be null '() or (list astring).
(define shelxl-refine-primitive
  (lambda (imol ins-text hkl-file-in-maybe)

    (let ((func (lambda (ins-file-name)
		  (call-with-output-file ins-file-name
		    (lambda (port) 
		      (display ins-text port))))))

      (shelxl-refine-inner imol hkl-file-in-maybe func))))


(define shelxl-refine-inner
  (lambda (imol hkl-file-in-maybe func)

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
			      ((eq? '() hkl-file-in-maybe) #f) 
			      (else (car hkl-file-in-maybe)))))

	  (if (not (file-exists? *coot-shelxl-dir*))
	      (mkdir *coot-shelxl-dir*))

	  (if (not (file-exists? *coot-shelxl-dir*))
	      ; We failed to make the directory for some reason.
	      ; Exit now.
	      (format #t "Failed to make shelxl directory ~s~% " 
		      *coot-shelxl-dir*)

	      (let* ((stub (string-append 
			    *coot-shelxl-dir* "/"
			    (strip-path (remove-time-extension 
					 (strip-extension (molecule-name imol))))
			    "-"
			    (unique-date/time-str)))
		     (ins-filename (string-append stub ".ins"))
		     (res-filename (string-append stub ".res"))
		     (lst-filename (string-append stub ".lst"))
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

		      (func ins-filename)

		      ;; changed
		      ;; (write-shelx-ins-file imol ins-filename)
		      ;;
		      ;; need to write text to ins-filename instead.
		      ;; 

		      ; running shelxl creates stub.res
		      (goosh-command "shelxl" (list stub) '() log-filename #t)
		      ; it isn't a pdb file, but coot knows what to do.
		      (let ((imol-res (read-pdb res-filename)))
			(handle-shelx-fcf-file fcf-filename)
			(read-shelx-lst-file lst-filename imol-res))))))))))

    
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
	
	; (format #t "debug: interesting-list: ~s~%" interesting-list)

	(if (null? interesting-list)
	    (begin
	      (format #t "INFO:: Nothing Interesting from LST file~%")
	      (add-status-bar-text "Nothing Interesting from LST file"))
	    (interesting-things-gui "Interesting Things from SHELX"
				    interesting-list))))

    ;; chop off last char of string
    (define (chop-end s)
      (let ((l (string-length s)))
	(if (= l 0)
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
	      (let ((resno (string->number (car (cdr ls)))))
		(if (number? resno)
		    (let ((chain-id (chain-id-for-shelxl-residue-number imol resno)))
		      (if (eq? chain-id #f) 
			  (begin
			    (format #t "couldn't find chain id for resno ~s~%" resno)
			    (list "" resno  "" (car ls) ""))
			  (list chain-id resno "" (car ls) "")))

		    ;; handle the alt conf
		    (let* ((s (car (cdr ls)))
			   (p (chop-end s))
			   (alt-conf (last-char s))
			   (n (string->number p)))
		      
		      (if (number? n)
			  (let ((chain-id (chain-id-for-shelxl-residue-number imol n)))
			    (if (eq? #f chain-id) 
				(begin
				  (format #t "couldn't find chain id for resno ~s~%" resno)
				  (list "" n "" (car ls) alt-conf))
				(list chain-id n "" (car ls) alt-conf)))
			  (list "" 1 "" "blank" ""))))) ; failure
		    
		(list "" 1 "" "blank" ""))))) ; failure

    ;; So we have a line in a table something like these:
    ;; 
    ;; "    2.7499    2.5580    0.1919    0.0400    DANG CN_2011 C_2011"
    ;; "                        1.5156    0.5000    FLAT O_1012 CA_1012 N_1013 CA_1013"
    ;; 
    ;; It seems that there are either 4 leading numbers or 2 leading
    ;; numbers.  Let's parse up the line here and return (list
    ;; observed target error sigma restraint) where observed and
    ;; target are #f if there number is missing (e.g. a "FLAT"
    ;; restraint).
    ;; 
    ;; There are 2 regimes we understand, either 
    ;; number number number number restraint-descr atom-1-desc atom-2-desc 
    ;; or 
    ;; number number restraint-descr atom-stuff [.. atom-stuff]
    ;; 
    ;; Actually here we will just make a string list of the
    ;; atom-stuffs (after the restraint-descr) and parse that
    ;; elsewhere (in do-gui) to find the atoms.  The actual
    ;; meaning of the atom-stuffs depends on the restraint-descr.
    ;; 
    ;; return #f on something we don't understand 
    ;; 
    (define (parse-dr-line line)
      
      (let ((bits (string->list-of-strings line)))

	(if (<= (length bits) 3)
	    #f ; failed to understand
	    ;; 
	    (let ((n0 (string->number (list-ref bits 0)))
		  (n1 (string->number (list-ref bits 1)))
		  (n2 (string->number (list-ref bits 2)))
		  (n3 (string->number (list-ref bits 3))))

	      (if (not (and (number? n0)
			    (number? n1)))

		  #f ; failed to understand

		  (if (and (number? n2)
			   (number? n3))
		      ; a 4 number line
		      (if (<= (length bits) 6)
			  #f ; failed to understand
			  (list n0 n1 n2 n3 (list-ref bits 4)
				(cdr (cdr (cdr (cdr (cdr bits)))))))

		      ; a 2 number line
		      (if (<= (length bits) 3)
			  #f ; failed to understand
			  (list #f #f n0 n1 (list-ref bits 2)
				(cdr (cdr (cdr bits)))))))))))
				

    ;; 
    (define (do-gui disagreeable-restraints-list interesting-list)
      ; (format #t "DR: ~s~%" disagreeable-restraints-list)

      (let loop ((dr-list (reverse disagreeable-restraints-list))
		 (dis-res-button-list '()))
	
	(cond 
	 ((null? dr-list) 
	  (gui-interesting-list (append interesting-list dis-res-button-list)))
	 (else 
	  (let* ((restraint-type (list-ref (car dr-list) 4))
		 (drl-index
		  (cond
		   ((string=? restraint-type "BUMP") 0)
		   ((string=? restraint-type "DANG") 0)
		   ((string=? restraint-type "FLAT") 0)
		   ((string=? restraint-type "SIMU") 1)
		   ((string=? restraint-type "ISOR") 1)
		   (else 
		    0)))
		 (atom-parts (make-atom-parts 
			      (list-ref (list-ref (car dr-list) 5) drl-index)))
		 (stats-string
		  (let ((n2 (list-ref (car dr-list) 2))
			(n3 (list-ref (car dr-list) 3)))
		    (if (not (and (number? n2)
				  (number? n3)))
			""
			(let ((z (/ (abs n2) n3)))
			  (if (not (number? z))
			      ""
			      (string-append 
			       " "
			       (number->string n2)
			       " [Z="
			       (number->string z)
			       "]"))))))
		 (button-label (string-append "Disagreeable Restraint " 
					      (if (= drl-index 0) 
						  restraint-type
						  (string-append
						   restraint-type " "
						   (car (list-ref (car dr-list) 5))))
					      "  "
					      (list-ref atom-parts 0)
					      "  "
					      (number->string (list-ref atom-parts 1))
					      "  "
					      (list-ref atom-parts 3)
					      stats-string))
		 (interesting-thing (cons button-label (cons imol atom-parts))))

	    (loop (cdr dr-list)
		  (cons interesting-thing dis-res-button-list)))))))
	    

    ;; main body
    (if (not (valid-model-molecule? imol))
	(format #t "WARNING:: Molecule number ~s not valid~%" imol)
	(if (not (file-exists? file-name))
	    (format #t "WARNING:: shelx lst file: ~s does not exist~%" 
		    file-name)
	    (call-with-input-file file-name
	      (lambda (port)
		
		(let loop ((line (read-line port))
			   (split-list '())
			   (disagreeable-restraints-list '())
			   (dr-count 0))
		  
		  (cond 
		   ((eof-object? line) (do-gui (reverse disagreeable-restraints-list)
					       (reverse split-list)))
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
				       split-list)
				      disagreeable-restraints-list
				      0)))))))

		   ((string-match 
		     "   Observed   Target    Error     Sigma     Restraint" line)
		    (loop (read-line port)
			  split-list
			  '()
			  1))

		   ((= dr-count 1)
		    (loop (read-line port)
			  split-list
			  '() ; reset the disagreeable-restraints-list
			  2))

		   ((= dr-count 2) ; OK, in the disagreeable table
		    ;; (format #t "DR: ~s~%" line)
		    (let ((dr-bits (parse-dr-line line)))
		      ;; (format #t "dr-bits: ~s~%" dr-bits)
		      (loop (read-line port)
			    split-list
			    (if (list? dr-bits)
				(cons dr-bits disagreeable-restraints-list)
				disagreeable-restraints-list)
			    (if (> (string-length line) 0) 2 0))))

		   (else 
		    (loop (read-line port)
			  split-list
			  disagreeable-restraints-list dr-count))))))))))

;;; Read a shelx project (i.e. the .res file, the .lst file and the
;;; .fcf.cif file (if it exists (else read the .fcf file (if it
;;; exists)))
;;;
;;; If the file-name has an extension of .lst .res .ins .fcf .log .hkl
;;; then strip that off before adding new extensions.
;;;
(define (read-shelx-project file-name)

  (let ((extension (file-name-extension file-name)))
    (let ((file-stub 
	   (cond
	    ((string=? extension "lst")
	     (file-name-sans-extension file-name))
	    ((string=? extension "ins")
	     (file-name-sans-extension file-name))
	    ((string=? extension "pdb")
	     (file-name-sans-extension file-name))
	    ((string=? extension "log")
	     (file-name-sans-extension file-name))
	    ((string=? extension "hkl")
	     (file-name-sans-extension file-name))
	    ((string=? extension "res")
	     (file-name-sans-extension file-name))
	    ((string=? extension "fcf")
	     (file-name-sans-extension file-name))
	    (else
	     file-name))))

      (let ((res-file-name (string-append file-stub ".res"))
	    (lst-file-name (string-append file-stub ".lst"))
	    (fcf-file-name (string-append file-stub ".fcf"))		       
	    (fcf-cif-file-name (string-append file-stub ".fcf.cif")))

	(format #t "  file-name: ~s~%" file-name)
	(format #t "  file-stub: ~s~%" file-stub)
	(format #t "  extension: ~s~%" extension)

	(let ((imol-res
	       (if (file-exists? res-file-name)
		   (begin
		     (format #t "Read res file ~s~%" res-file-name)
		     (handle-read-draw-molecule-with-recentre res-file-name 0))
		   (begin
		     (format #t "   No res file ~s~%" res-file-name)
		     -1))))

	  (if (not (valid-model-molecule? imol-res))
	      (format #t "WARNING:: Bad molecule from res file read.~%")
	      (begin
		(if (file-exists? fcf-cif-file-name)
		    (if (not (file-exists? fcf-file-name))
			(begin
			  (format #t "   Read fcf-cif file ~s~%" fcf-cif-file-name)
			  (auto-read-cif-data-with-phases fcf-cif-file-name))

			; OK both xxx.fcf and xxx.fcf.cif exist, we
			; only want to read the xxx.fcf.cif if it is
			; more recent than xxx.fcf (if it is not, we
			; want to handle-shelx-fcf-file.
			(let ((    fcf-date (stat:mtime (stat     fcf-file-name)))
			      (fcf-cif-date (stat:mtime (stat fcf-cif-file-name))))
;			  (format #t "    fcf-date ~s~%fcf-cif-data ~s~%" 
;				  fcf-date fcf-cif-date)
			  (if (< fcf-date fcf-cif-date)
			      (auto-read-cif-data-with-phases fcf-cif-file-name)
			      (handle-shelx-fcf-file fcf-file-name))))
			  
		    ;; xxx.fcf.cif did not exist:
		    (if (file-exists? fcf-file-name)
			(begin
			  (format #t "   Read fcf file ~s~%" fcf-file-name)
			  (handle-shelx-fcf-file fcf-file-name))))
		
		(if (file-exists? lst-file-name)
		    (begin 
		      (format #t "   ::Read lst file ~s~%" lst-file-name)
		      (read-shelx-lst-file lst-file-name imol-res))
		    (format #t "   ::No lst file ~s~%" lst-file-name)))))))))

