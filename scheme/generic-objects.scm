;;;; Copyright 2006 by The University of York
;;;; Author: Paul Emsley
;;;; Copyright 2008 by The University of York
;;;; Copyright 2009 by the University of Oxford
;;;; Author: Paul Emsley
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
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


;; map to scheme names:
(define generic-object-is-displayed? generic-object-is-displayed-p)

;; map to scheme names:
(define is-closed-generic-object? is-closed-generic-object-p)

;; return a new generic object number for the given object obj-name.
;; If there is no such object with name obj-name, then create a new
;; one.  (Maybe you want the generic object to be cleaned if it exists
;; before returning, this function does not do that).
;; 
(define generic-object-with-name 
  (lambda (obj-name)
      
    (let ((t (generic-object-index obj-name)))
      (if (= t -1)
	  (new-generic-object-number obj-name)
	  t))))



;; display a GUI for generic objects
;; 
(define generic-objects-gui
  (lambda ()
    
    (if (using-gui?) 
	(let ((n-objects (number-of-generic-objects)))

	  (if (> n-objects 0)
	      (begin

		(let* ((window (gtk-window-new 'toplevel))
		       (vbox (gtk-vbox-new #f 10))
		       (label (gtk-label-new " Generic Objects ")))

		  (gtk-container-add vbox label)
		  (gtk-widget-show label)

		  (for-each 
		   (lambda (generic-object-number)

		     (format #t "~s ~s ~s~%" generic-object-number (generic-object-name generic-object-number)
			     (is-closed-generic-object? generic-object-number))

		     (if (= (is-closed-generic-object? generic-object-number) 0)
			 (let ((name (generic-object-name generic-object-number)))
			   
			   (if name 
			       (let* ((label (string-append 
					      (number->string generic-object-number)
					      "  "
					      name))
				      (frame (gtk-frame-new #f))
				      (check-button (gtk-check-button-new-with-label label)))
				 
					; this callback gets called by the
					; gtk-toggle-button-set-active just below,
					; which is why we need the state and active
					; test.
				 (gtk-signal-connect check-button "toggled" 
						     (lambda () 
						       (let ((button-state (gtk-toggle-button-active check-button))
							     (object-state (generic-object-is-displayed? 
									    generic-object-number)))
							 (if (and (eq? button-state #t)
								  (= object-state 0))
							     (set-display-generic-object generic-object-number 1))
							 (if (and (eq? button-state #f)
								  (= object-state 1))
							     (set-display-generic-object generic-object-number 0)))))
				 
				 (let ((current-state 
					(generic-object-is-displayed? 
					 generic-object-number)))
				   (if (= current-state 1)
				       (gtk-toggle-button-set-active check-button #t)))
				 
				 (gtk-container-add vbox frame)
				 (gtk-container-add frame check-button)
				 (gtk-widget-show frame)
				 (gtk-widget-show check-button))))))
		   
		   (number-list 0 (- n-objects 1)))

		  (gtk-widget-show vbox)
		  (gtk-container-add window vbox)
		  (gtk-container-border-width window 10)
		  (gtk-widget-show window))))))))


;; return status.
;; 
(define (reduce-on-pdb-file imol pdb-in pdb-out)
    (reduce-on-pdb-file-generic imol 'build pdb-in pdb-out))


;; return status.
;; 
(define (reduce-on-pdb-file-no-flip imol pdb-in pdb-out)
    (reduce-on-pdb-file-generic imol 'no-flip pdb-in pdb-out))

;; return status.
;; 
(define (reduce-on-pdb-file-generic imol no-flip-or-build pdb-in pdb-out)

  (format #t "running reduce on ~s~%" pdb-in)
  (if (not (command-in-path-or-absolute? *reduce-command*))
      (format #t "command for reduce (~s) is not in path~%" *reduce-command*)
      (begin
	
	(let* ((nshl (non-standard-residue-names imol))
	       (ext (apply string-append (map (lambda (item)
						(string-append item "-"))
					      nshl)))
               (mode "-build")
	       (reduce-het-dict-file-name (string-append 
					   "coot-molprobity/reduce-het-dict" ext ".txt")))

	  (write-reduce-het-dict imol reduce-het-dict-file-name)

	  ;; As Bernie, let's try to set REDUCE_HET_DICT if
	  ;; we can, assuming that the het-dict is in the
	  ;; same directory as is reduce-command.  If it
	  ;; isn't then do nothing.
	  ;; 
	  (let ((h (getenv "REDUCE_HET_DICT")))
	    (if (not (string? h))
		(let* ((d (file-name-directory *reduce-command*)))
		  (if (string? d)
		      (let ((f (append-dir-file d "reduce_wwPDB_het_dict.txt")))
			(if (file-exists? f)
			    (let ((env-string (string-append "REDUCE_HET_DICT=" f)))
			      (format #t "puting env ~s~%" env-string)
			      (putenv env-string))))))))

	  (if (eq? no-flip-or-build 'no-flip)
	      (set! mode "-NOFLIP"))

	  (format #t "======= reduce-on-pdb-file: command ~s args ~s with pdb-out: ~s~%"
		  *reduce-command*
		  (list mode pdb-in "-DB" reduce-het-dict-file-name)
		  pdb-out)
	  
	  (let ((status (goosh-command *reduce-command* 
				       (list mode pdb-in "-DB" reduce-het-dict-file-name)
				       '() pdb-out #f)))
	    
	    (format #t "======== status ~s~%" status)
	    status)))))



(define *old-pdb-style* #f)

(define reduce-molecule-updates-current #f)

(define *probe-command* "probe")
(define *reduce-command* "reduce")

;; run molprobity (well reduce and probe) to make generic objects (and
;; display the generic objects gui)
;; 
(define probe
  (lambda (imol)

    (if (valid-model-molecule? imol)
	(if (not (command-in-path-or-absolute? *probe-command*))
	    (format #t "command for probe (~s) is not in path~%" *probe-command*)
	    (if (not (command-in-path-or-absolute? *reduce-command*))
		(format #t "command for reduce (~s) is not in path~%" *reduce-command*)
		(begin
		  (make-directory-maybe "coot-molprobity")
		  (let ((mol-pdb-file  "coot-molprobity/for-reduce.pdb")
			(reduce-out-pdb-file "coot-molprobity/reduced.pdb")
			(reduce-het-dict-file-name "coot-molprobity/reduce-het-dict.txt"))
		    (write-pdb-file imol mol-pdb-file)
		    (write-reduce-het-dict imol reduce-het-dict-file-name)

		    ;; As Bernie, let's try to set REDUCE_HET_DICT if
		    ;; we can, assuming that the het-dict is in the
		    ;; same directory as is reduce-command.  If it
		    ;; isn't then do nothing.
		    ;; 
		    (let ((h (getenv "REDUCE_HET_DICT")))
		      (if (not (string? h))
			  (let* ((d (file-name-directory *reduce-command*)))
			    (if (string? d)
				(let ((f (append-dir-file d "reduce_wwPDB_het_dict.txt")))
				  (if (file-exists? f)
				      (let ((env-string (string-append "REDUCE_HET_DICT=" f)))
					(putenv env-string))))))))


		    (let ((arg-list 
			   (if *old-pdb-style*
			       ;; old 
			       (list "-build" "-oldpdb" mol-pdb-file "-DB" reduce-het-dict-file-name)
			       ;; modern
			       (list "-build" mol-pdb-file "-DB" reduce-het-dict-file-name))))

		      (format #t "============= running reduce: ~s ~s and output to: ~s~%"
			      *reduce-command* arg-list reduce-out-pdb-file)
		      (format #t "============= running reduce: REDUCE_HET_DICT env var: ~s~%" 
			      (getenv "REDUCE_HET_DICT"))
		      (goosh-command *reduce-command* arg-list '() reduce-out-pdb-file #f)
		      (let* ((probe-name-stub (strip-extension (strip-path (molecule-name imol))))
			     (probe-pdb-in  (string-append 
					     "coot-molprobity/" probe-name-stub "-with-H.pdb"))
			     (probe-out "coot-molprobity/probe-dots.out"))
			
			(goosh-command "grep" (list "-v" "^USER" reduce-out-pdb-file)
				       '() probe-pdb-in #f)
			(goosh-command *probe-command* (list "-u"  "-mc"
					; "'(chainA,chainZ) alta'" 
							     "ALL"
							     probe-pdb-in)
				       '() probe-out #f)
					; by default, we don't want to click on the
					; imol-probe molecule (I think :-)
			(let* ((recentre-status (recentre-on-read-pdb))
			       (novalue (set-recentre-on-read-pdb 0))
			       (imol-probe 
				(if reduce-molecule-updates-current
				    (begin
				    (format #t "======= update molecule =======~%")
				    (clear-and-update-model-molecule-from-file imol probe-pdb-in))
				    (begin
				      (format #t "======= read new pdb file =======~%")
				      (read-pdb probe-pdb-in)
				      imol))))
			  (if (= 1 recentre-status)
			      (set-recentre-on-read-pdb 1))
			  
			  ;; show the GUI for the USER MODS
			  (if (using-gui?)
			      (user-mods-gui imol-probe reduce-out-pdb-file))
			  
					; (toggle-active-mol imol-probe) let's not do
					; that actually.  I no longer think that the
					; new probe molecule should not be clickable
					; when it is initally displayed (that plus
					; there is some active/displayed logic problem
					; for the molecules, which means that after
					; several probes, the wrong molecule is
					; getting refined).
			  (handle-read-draw-probe-dots-unformatted probe-out imol-probe 1)
			  (generic-objects-gui)
			  (graphics-draw)))))))))))

;; Write the connectivity for the non-standard (non-water) residues in
;; the given molecule for which we have the dictionary.
;; 
;; Don't return anything interesting.  
;; 
(define (write-reduce-het-dict imol reduce-het-dict-file-name)
  
  (let ((con-file-names 
	 (filter string? 
		 (map (lambda (res-name)
			(let* ((f-name (string-append 
					"coot-molprobity/conn-"
					res-name
					".txt"))
			       (status (write-connectivity res-name f-name)))
			  (if (= status 1)
			      f-name
			      #f)))
		      (non-standard-residue-names imol)))))

    (if (not (null? con-file-names))
	(call-with-output-file reduce-het-dict-file-name
	  (lambda (port-out)
	    (for-each (lambda (conn-file)
			(call-with-input-file conn-file
			  (lambda (port)
			    (let f ((line (read-line port)))
			      (if (not (eof-object? line))
				  (begin
				    (display line port-out)
				    (newline port-out)
				    (f (read-line port))))))))
		      con-file-names))))))
			

		    
;; gets set the first time we run interactive-probe.  Takes values
;; unset (initial value) 'yes and 'no)
;; 			
(define *interactive-probe-is-OK?* 'unset)

;; run "probe" interactively, 
;; which in the current implementation, means that this function 
;; can run during a edit-chi angles manipulation, or after
;; a real space refine zone.
;; 
;; Thus function presumes that there are 2 pdb files in the current 
;; directory, one of which is the reference pdb file and the other
;; is a pdb file containing the tmp/moving atom set.
;; 
;; The function takes arguments for the centre of the probe dots
;; and the radius of the probe dots sphere.  The chain id and 
;; residue number are also needed to pass as arguments to probe.
;; 
(define interactive-probe
  (lambda (x-cen y-cen z-cen radius chain-id res-no)

;    (make-directory-maybe "coot-molprobity")
;     (if (not (command-in-path? *probe-command*))
; 	(format #t "command for probe (~s) is not in path~%" *probe-command*)

    (let* ((probe-pdb-in-1 "molprobity-tmp-reference-file.pdb")
	   (probe-pdb-in-2 "molprobity-tmp-moving-file.pdb")
	   (probe-out "coot-molprobity/molprobity-tmp-probe-dots.out")
	   (atom-sel (string-append 
		      "(file1 within "
		      (number->string radius)
		      " of "
		      (number->string x-cen) ", "
		      (number->string y-cen) ", "
		      (number->string z-cen) ", "
		      "not water not ("
		      (if (string=? chain-id "") 
			  ""
			  "chain")
		      (string-downcase chain-id) " "
		      (number->string res-no)
		      ")),file2")))

      (format #t "probe command ~s ~s~%" *probe-command* 
	      ;; no longer use "--std-bonds"
	      (list "-mc" "-u" "-quiet" "-drop" "-both" 
		    atom-sel "file2" probe-pdb-in-1 probe-pdb-in-2))

      ;; if unset, then set it.
      (if (eq? *interactive-probe-is-OK?* 'unset)
	  (if (command-in-path-or-absolute? *probe-command*)
	      (set! *interactive-probe-is-OK?* 'yes)
	      (set! *interactive-probe-is-OK?* 'no)))

      (if (eq? *interactive-probe-is-OK?* 'yes)
	  (begin
	    (goosh-command *probe-command* 
			   (list "-mc" "-u" "-quiet" "-drop" "-stdbonds" 
				 "-both" atom-sel "file2"
				 probe-pdb-in-1 probe-pdb-in-2)
			   '() probe-out #f)
	    ;; don't show the gui, so the imol is not needed/dummy.
	    (handle-read-draw-probe-dots-unformatted probe-out 0 0)
	    (graphics-draw))))))

;; 
;;
(define (get-probe-dots-from pdb-file-name point radius)

  ;; if unset, then set it, try to make dir too
  (if (eq? *interactive-probe-is-OK?* 'unset)
      (if (not (command-in-path-or-absolute? *probe-command*))
	  (set! *interactive-probe-is-OK?* 'no)
	  (let ((dir-name (get-directory "coot-molprobity")))
	    (if (string? dir-name) ;; OK, we had it or made it
		(set! *interactive-probe-is-OK?* 'yes)
		(set! *interactive-probe-is-OK?* 'no)))))
  
  (if (eq? *interactive-probe-is-OK?* 'yes)
      (let* ((probe-out "coot-molprobity/molprobity-tmp-probe-dots.out")
	     (within-str (string-append "(within "
					(number->string radius)
					" of "
					(number->string (list-ref point 0))
					", "
					(number->string (list-ref point 1))
					", "
					(number->string (list-ref point 2))
					")"))
	     (args (list "-mc" "-u" "-quiet" "-drop" "-stdbonds" 
			 "ALL" ;; whole residues from sphere selection
			       ;; was needed to make this work
			 ;; within-str ;; problems with atom selection
			 pdb-file-name)))
	(format #t "goosh-command on ~s ~s~%" *probe-command* args)
	(goosh-command *probe-command* args '() probe-out #f)
	(handle-read-draw-probe-dots-unformatted probe-out 0 0)
	(graphics-draw))))



;; Update the generic objects probe dots from residues within radius
;; of the screen centre.
;; 
;; Return nothing interesting.
;; 
(define (probe-local-sphere imol radius)

  ;; We need to select more atoms than we probe because if the atom
  ;; selection radius and the probe radius are the same, then
  ;; sometimes the middle atom of a bonded angle set is missing
  ;; (outside the sphere) and that leads to bad clashes.
  ;; There are also bad clashed at the edge when alt-confed atoms 
  ;; are not selected but non-alt-confs are which leads to missing atoms
  ;; in a bond angle and therefore clashes.

  (let* ((pt (rotation-centre))
	 (imol-new (apply new-molecule-by-sphere-selection 
		    imol (append pt (list radius 0)))))

    (set-mol-displayed imol-new 0)
    (set-mol-active    imol-new 0)
    (let ((pdb-name "molprobity-tmp-reference-file.pdb"))
      (make-directory-maybe "coot-molprobity")
      (write-pdb-file imol-new pdb-name)
      (get-probe-dots-from pdb-name pt radius)
      (close-molecule imol-new))))

