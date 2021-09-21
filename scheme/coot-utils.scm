
;;;; Copyright 2000 by Paul Emsley
;;;; Copyright 2004, 2005, 2006, 2007 by The University of York
;;;; Copyright 2008 by The University of York
;;;; Copyright 2008, 2009, 2010, 2011, 2012 by The University of Oxford
;;;; Copyright 2013, 2014, 2015, 2016 by Medical Research Council

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

(use-modules (ice-9 popen)
	     (ice-9 string-fun)
	     (ice-9 format)
	     (ice-9 rdelim)
	     (srfi srfi-1))
(use-modules (os process)) 

;; 3D annotations - a bit of a hack currently
(define *annotations* '())

;; used in Extensions -> Representation -> Ball & Stick
(define *default-ball-and-stick-selection* "//A/1-2")

(load-by-search "redefine-functions.scm")


;; documented functions


;; rigid body refine using residue ranges.  Takes 2 arguments, the
;; first is the molecule number, the second is a list of
;; residue-ranges.  A residue range is 
;; @code{(list @emph{chain-id resno-start resno-end})}.
;; 
(define rigid-body-refine-by-residue-ranges rigid-body-refine-by-residue-ranges-scm)

;; add terminal residue is the normal thing we do with an aligned
;; sequence, but also we can try ton find the residue type of a
;; residue in the middle of the chain that is modelled as an ALA, say.
;; 
(define find-aligned-residue-type find-terminal-residue-type)

(define (using-gui?)
  (defined? 'coot-main-menubar))


;; Macro to tidy up a a setup of functions to be run with no backup
;; for a particular molecule.
;;
;; func is a thunk.
;;  
(defmacro with-no-backups (imol func)
  
  `(begin
     (let ((b-state (backup-state ,imol)))
       (turn-off-backup ,imol)
       (let ((v ,func))
	 (if (= 1 b-state)
	     (turn-on-backup ,imol))
	 v)))) ;; return the value generated from running the function

;; Macro to tidy up a set of functions to be run with automatic
;; accepting of the refinement
;; 
;; funcs is a normal set of functions (not a thunk)
;; 
(defmacro with-auto-accept funcs

  `(begin
     (let ((replace-state (refinement-immediate-replacement-state)))
       (set-refinement-immediate-replacement 1)
       (let ((val ,@funcs))
	 (accept-regularizement)
	 (if (= replace-state 0)
	     (set-refinement-immediate-replacement 0))
	 val))))

;; 
(defmacro using-active-atom funcs

  `(begin
     (begin
       (let ((active-atom (closest-atom-simple-scm)))
	 (if (not active-atom)
	     (begin
	       (add-status-bar-text "No residue found"))

	     (let* ((aa-imol      (list-ref active-atom 0))
		    (aa-chain-id  (list-ref active-atom 1))
		    (aa-res-no    (list-ref active-atom 2))
		    (aa-ins-code  (list-ref active-atom 3))
		    (aa-atom-name (list-ref active-atom 4))
		    (aa-alt-conf  (list-ref active-atom 5))
                    (aa-res-spec (list aa-chain-id aa-res-no aa-ins-code)))

	       ,@funcs
	       
	       ))))))

; e.g.:
; (with-auto-accept 
;   (format #t "tum tee tumm...~%")
;   (format #t "tra la la...~%"))

(defmacro print-var (var)
  (let ((s (format #f "~s" var)))
    `(format #t "DEBUG:: ~a is ~s~%" ,s ,var)))

;; not tail-recursive
;; 
(define (all-true? ls)
  (cond 
   ((null? ls) #t)
   ((car ls) (all-true? (cdr ls)))
   (else #f)))


;; make the directory and return the directory name.  If you cant make
;; it directly in this directory try to make it in $HOME.  Return #f if
;; complete failure.  e.g. coot-ccp4 or coot-backup
;; 
(define (get-directory dir-name)

  (if (file-exists? dir-name)

      (if (is-directory? dir-name)
	  dir-name
	  #f)

      (let ((status (make-directory-maybe dir-name)))

	(if (= status 0)
	    dir-name

	    (let ((h (getenv "HOME")))
	      (if (not (string? h))
		  #f ;; couldnt find home.  Baah.
		  (begin
		    (let ((new-dir (append-dir-dir h dir-name)))
		      (let ((status (make-directory-maybe new-dir)))
			(if (= status 0)
			    new-dir
			    #f))))))))))
      

;; schemify function
;; 
(define (molecule-has-hydrogens? imol)
  (= (molecule-has-hydrogens-raw imol) 1))


(define (add-hydrogens-using-refmac imol)
  (let ((out-file-name (append-dir-file "coot-refmac"
					(string-append (molecule-name-stub imol 0)
						       "-needs-H.pdb")))
	(in-file-name (append-dir-file "coot-refmac"
				       (string-append (molecule-name-stub imol 0)
						      "-with-H.pdb"))))
    (make-directory-maybe "coot-refmac")
    (write-pdb-file imol out-file-name)
    (add-hydrogens-using-refmac-inner imol in-file-name out-file-name)))


(define (add-hydrogens-to-chain-using-refmac imol chain-id)
  (let ((out-file-name (append-dir-file "coot-refmac"
					(string-append (molecule-name-stub imol 0)
						       "-chain-"
						       chain-id
						       "-needs-H.pdb")))
	(in-file-name (append-dir-file "coot-refmac"
				       (string-append (molecule-name-stub imol 0)
						      "-chain-"
						      chain-id
						      "-with-H.pdb"))))
    (make-directory-maybe "coot-refmac")
    (write-chain-to-pdb-file imol chain-id out-file-name)
    (add-hydrogens-using-refmac-inner imol in-file-name out-file-name)))


;; out-file-name has been written, in-file-name needs to be read.
;; 
(define (add-hydrogens-using-refmac-inner imol in-file-name out-file-name)

    (let ((goosh-status (goosh-command
			 "refmac5"
			 (list "XYZIN" out-file-name  ;; yes (coot out)
			       "XYZOUT" in-file-name)
			 (list "MAKE HOUT YES" "NCYCLE 0" "END")
			 "refmac-H-addition.log"
			 #f)))
      (if (not (ok-goosh-status? goosh-status))
	  (begin 
	    (format #f "WARNING:: problem running refmac5 on ~s~%" out-file-name)
	    #f)
	  (add-hydrogens-from-file imol in-file-name))))



;; set this to a function accepting two argument (the molecule number
;; and the manipulation mode) and it will be run after a model
;; manipulation.
;;
;; The manipulation mode will be one of (MOVINGATOMS), (DELETED) or
;; (MUTATED) and these can be tested with "=".
;;
;; e.g.
;;
;; (if (= mode (DELETED))
;      (display "Something was deleted"))
;; 
(define post-manipulation-hook #f)

;; do something based on the active residue (presumably)
(define post-set-rotation-centre-hook #f)

;; a function that takes one argument (the model molecule number)
(define post-read-model-hook #f)

;; Return a boolean
;; 
(define (pre-release?)
  (string-match "-pre" (coot-version)))


;; Return a list of molecule numbers (closed and open) The elements of
;; the returned list need to be tested against
;; is-valid-model-molecule?
;; 
(define molecule-number-list
  (lambda ()
    (let loop ((ls (range (graphics-n-molecules)))
	       (acc '()))
      (cond
       ((null? ls) (reverse acc))
       ((or (valid-map-molecule? (car ls))
	    (valid-model-molecule? (car ls)))
	(loop (cdr ls) (cons (car ls) acc)))
       (else 
	(loop (cdr ls) acc))))))

(define (model-molecule-number-list)
  (filter valid-model-molecule? (molecule-number-list)))

(define (display-all-maps)
  (let ((map-list (map-molecule-list)))
    (for-each (lambda (imol) (set-map-displayed imol 1)) map-list)))


;; c.f. graphics_info_t::undisplay_all_model_molecules_except(int imol)
(define (undisplay-all-maps-except imol-map)

  (format #t "undisplay-all-maps-except imol-map: ~s~%" imol-map)

  (let ((map-list (map-molecule-list)))
    (for-each (lambda (imol)
		(if (not (= imol imol-map))
		    (set-map-displayed imol 0)))
	      map-list)
    (set-map-displayed imol-map 1)))

(define (display-cycle-through-maps)

  ;; return an index or #f, lst must be a list
  (define (find-in-list item lst)
    (let loop ((idx 0))
      (cond
       ((= (length lst) idx) #f)
       ((eq? (list-ref lst idx) item) idx)
       (else
	(loop (+ idx 1))))))

  (define (next-map current-map-number map-number-list)
    (let ((current-idx (find-in-list current-map-number map-number-list))
	  (l (length map-number-list)))
      (format #t "current-idx: ~s from list map-number-list ~s~%" current-idx map-number-list)
      (if (number? current-idx)
	  (let ((next-index (if (= (+ current-idx 1) l)
				0
				(+ current-idx 1))))
	    (list-ref map-number-list next-index))
	  (list-ref map-number-list 0))))

  (let ((map-list (map-molecule-list)))
    (let ((current-displayed-maps
	   (filter (lambda(imol)
		     (= (map-is-displayed imol) 1))
		   map-list)))
      (let ((n-displayed (length current-displayed-maps)))

      ;; if nothing is displayed, display the first map in map-list
      ;; if one map is displayed, display the next map in map-list
      ;; if more than one map is displayed, display only the last map
      ;;    in the current-displayed-maps

      (cond
       ((= n-displayed 0) (if (> (length map-list) 0)
			      (undisplay-all-maps-except (car map-list))))
       ((= n-displayed 1) (if (> (length map-list) 1)
                              (let ((nm (next-map (car current-displayed-maps) map-list))
                                    (currently-displayed-map (car current-displayed-maps)))
                                (if (> nm currently-displayed-map)
                                    (undisplay-all-maps-except nm)
                                    (display-all-maps)))))

       (else (undisplay-all-maps-except (car current-displayed-maps))))))))

;; is isn't quite because one of the options is "all"
(define (just-one-or-next-map)
  (display-cycle-through-maps))


;; first n fields of ls. if length ls is less than n, return ls.
;; if ls is not a list, return ls.  If n is negative, return ls.
;; 
(define (first-n n ls)
  (if (not (list? ls))
      ls
      (if (< n 0)
	  ls 
	  (if (<= (length ls) n)
	      ls
	      (let loop ((r (reverse ls))
			 (count (- (length ls) n)))
		(cond
		 ((= 0 count) (reverse r))
		 (else 
		  (loop (cdr r) (- count 1)))))))))

		    
	
;; Test for prefix-dir (1) being a string (2) existing (3) being a
;; directory (4) modifiable by user (ie. u+rwx).  prefix-dir must be a
;; string.
;; 
;; Return #t or #f.
;; 
(define (directory-is-modifiable? prefix-dir)
  (if (not (string? prefix-dir))
      #f
      (if (not (file-exists? prefix-dir))
	  #f 
	  (let* ((s (stat prefix-dir))
		 (t (stat:type s)))
	    (if (not (eq? t 'directory))
		#f ; not a directory
		(let ((p (stat:perms s)))
		  ;; test the file permissions for rwx for user using bitwise logical 
		  ;; operator on p (permissions). 448 is 256 + 128 + 64
		  (let ((b #b111000000))
		    (= b (logand p b)))))))))



;; return an absolute file-name for file-name or #f
;; 
(define (absolutify file-name)

  (if (not (string? file-name))
      #f
      (if (not (> (string-length file-name) 0))
	  "/"
	  (let ((first-char (substring file-name  0 1)))
	    (if (string=? first-char "/")
		file-name
		(append-dir-file (getcwd) file-name))))))

;; return the directory component of file-name, leave "/" on, if it's
;; there. Note "x", "", "/" -> ""
;; 
(define (file-name-directory file-name)
  ;; (format #t "in file-name-directory: trying to split ~s~%" file-name)
  (if (not (string? file-name))
      #f
      (let ((ls (split-after-char-last #\/ file-name list)))
	(let ((r (car ls)))
	  (if (string=? r "/")
	      ""
	      r)))))

;; return the file component of file-name
;; "x/y" -> "y" 
;; "x" -> "x"
(define (file-name-file file-name)

  (if (not (string? file-name))
      #f
      (let ((ls (split-after-char-last #\/ file-name list)))
	(car (reverse ls)))))



;; Find the most recently created file from the given glob and dir
;;
;; return #f on no-such-file
;; 
(define most-recently-created-file
  (lambda (glob-str dir)

    (define add-dir-prefix
      (lambda (file)
	(if (string? file)
	    (append-dir-file dir file)
	    file)))

    (let ((files (glob glob-str dir)))

      (let loop ((files files)
                 (latest-file #f)
                 (latest-mtime 0))
        (cond
         ((null? files) (add-dir-prefix latest-file))
         (else 
          (let ((this-mtime (stat:mtime (stat (add-dir-prefix (car files))))))
            (if (> this-mtime latest-mtime)
                (loop (cdr files)
                      (car files)
                      this-mtime)
                (loop (cdr files)
                      latest-file
                      latest-mtime)))))))))

		

;; Convert a residue-spec to an mmdb atom selection string.
;; 
(define (residue-spec->atom-selection-string centre-residue-spec)
  (string-append "//" (car centre-residue-spec)
		 "/" (number->string 
		      (car (cdr centre-residue-spec)))))

;; residue-info atom
(define (residue-atom->atom-name ra)
  (if (not (list? ra))
      #f
      (car (car ra))))

;; residue-info atom
(define (residue-atom->alt-conf ra)
  (if (not (list? ra))
      #f
      (cadr (car ra))))

;; residue-info atom
(define (residue-atom->occupancy ra)
  (if (not (list? ra))
      #f
      (car (cadr ra))))

;; residue-info atom
(define (residue-atom->position ra)
  (if (not (list? ra))
      #f
      (list-ref ra 2)))

;; residue-info atom needs other parameters to make a spec for an atom
(define (residue-atom->atom-spec ra chain-id res-no ins-code)
  (list chain-id res-no ins-code (residue-atom->atom-name ra) (residue-atom->alt-conf ra)))

;; residue spec (e.g. from residue-near-residue)
(define (residue-spec->chain-id rs)
  (if (not (list? rs))
      #f
      (if (= (length rs) 3)
	  (car rs)
	  (if (= (length rs) 4)
	      (list-ref rs 1)
	      #f))))

(define (residue-spec->res-no rs)
  (if (not (list? rs))
      #f
      (if (= (length rs) 3)
	  (list-ref rs 1)
	  (if (= (length rs) 4)
	      (list-ref rs 2)
	      #f))))

(define (residue-spec->ins-code rs)
  (if (not (list? rs))
      #f
      (if (= (length rs) 3)
	  (list-ref rs 2)
      (if (= (length rs) 4)
	  (list-ref rs 3)
	  #f))))

(define (residue-specs-match? spec-1 spec-2)
  (if (string=? (residue-spec->chain-id spec-1)
		(residue-spec->chain-id spec-2))
      (if (= (residue-spec->res-no spec-1)
	     (residue-spec->res-no spec-2))
	  (if (string=? (residue-spec->ins-code spec-1)
			(residue-spec->ins-code spec-2))
	      #t)
	  #f)
      #f))

(define (atom-spec->imol atom-spec)
  (if (not (list? atom-spec))
      #f
      (if (= (length atom-spec) 6)
	  (car atom-spec)
	  (if (= (length atom-spec) 7)
	      (cadr atom-spec)))))

(define (residue-spec->residue-name imol spec)
  (if (= (length spec) 4)
      (residue-name imol
		    (list-ref spec 1)
		    (list-ref spec 2)
		    (list-ref spec 3))
      (residue-name imol
		    (list-ref spec 0)
		    (list-ref spec 1)
		    (list-ref spec 2))))

;; for sorting residue specs
(define (residue-spec-less-than spec-1 spec-2)
  (let ((chain-id-1 (residue-spec->chain-id spec-1))
	(chain-id-2 (residue-spec->chain-id spec-2)))
    (if (string<? chain-id-2 chain-id-1)
	#t
	(let ((rn-1 (residue-spec->res-no spec-1))
	      (rn-2 (residue-spec->res-no spec-2)))
	  (if (< rn-2 rn-1)
	      #t
	      (let ((ins-code-1 (residue-spec->ins-code spec-1))
		    (ins-code-2 (residue-spec->ins-code spec-2)))
		(string<? ins-code-2 ins-code-1)))))))

(define (residue-spec->string spec)
  (string-append
   (residue-spec->chain-id spec)
   " "
   (number->string (residue-spec->res-no spec))
   (residue-spec->ins-code spec)))


;; Return a list of molecules that are maps
;; 
(define (map-molecule-list)
  
  (let loop ((molecule-list (number-list 0 (- (graphics-n-molecules) 1)))
	     (map-list '()))
    (cond
     ((null? molecule-list) (reverse map-list))
     ((valid-map-molecule? (car molecule-list))
      (loop (cdr molecule-list)
	    (cons (car molecule-list) map-list)))
     (else 
      (loop (cdr molecule-list)
	    map-list)))))

    
;; Return a list of molecules that are maps
;; 
(define (model-molecule-list)
  
  (let loop ((molecule-list (number-list 0 (- (graphics-n-molecules) 1)))
	     (model-list '()))
    (cond
     ((null? molecule-list) (reverse model-list))
     ((valid-model-molecule? (car molecule-list))
      (loop (cdr molecule-list)
	    (cons (car molecule-list) model-list)))
     (else 
      (loop (cdr molecule-list)
	    model-list)))))

;; Return #t (#f) if @var{imol} is (isn't) a shelx molecule.
;; 
(define (shelx-molecule? imol)
  (= (is-shelx-molecule imol) 1))

;; return an int. 0 means no, 1 means yes, -1 on error
;;
(define is-protein-chain? is-protein-chain-p)

;; Is a nucleotide chain?
;;
;; Now return a boolean.
;;
(define (is-nucleotide-chain? imol chain-id)

  (let ((result (is-nucleotide-chain-p imol chain-id)))
    (= result 1)))


(define (delete-residue-by-spec imol spec)
  (delete-residue imol
		  (residue-spec->chain-id spec)
		  (residue-spec->res-no   spec)
		  (residue-spec->ins-code spec)))

;;; No! don't define this.  It is misleading.  It can return 0, which
;;; is true!  use instead valid-model-molecule?
;;; 
;;; define is-valid-model-molecule? is-valid-model-molecule

    
;; Set the virtual trackball behaviour.
;; 
;; trackball @var{type} is a symbol: either 'flat or 'spherical-surface.
;; 
(define (set-virtual-trackball-type type)
  
  (cond 
   ((eq? type 'flat) (vt-surface 1))
   ((eq? type 'spherical-surface) (vt-surface 0))
   (else
					; usually not output anywhere
    (format #t "virtual trackball type ~s not understood~%" type))))
    
;; Is @var{ls} a list of strings?  Return #t or #f
;; 
(define (list-of-strings? ls)

  (if (not (list? ls))
      #f
      (let f ((ls ls))
	(cond 
	 ((null? ls) #t)
	 ((string? (car ls))
	  (f (cdr ls)))
	 (else #f)))))


;; (coot-replace-string "one two three" " " "_") -> "one_two_three"
;; 
(define (coot-replace-string string-in target replacing-str)

  (let ((sm (string-match target string-in)))
    (if (not sm)
	string-in
	(let ((new-string (string-replace string-in
					  replacing-str 
					  (car (vector-ref sm 1))
					  (cdr (vector-ref sm 1)))))
	  (coot-replace-string new-string target replacing-str)))))

		    

;; string concat with spaces, @var{ls} must be a list of strings.
;; 
(define (string-append-with-spaces ls)

    (if (null? ls) 
	""
	(string-append (car ls) " " (string-append-with-spaces (cdr ls)))))

;; The screen centre.
;; 
;; return the rotation centre as a 3 membered list of numbers
;; 
(define rotation-centre
  (lambda ()

    (map rotation-centre-position (list 0 1 2))))

(define rotation-center rotation-centre) ; maybe there is a better
					 ; place for US spellings.  (There is now.)



;; Make list of integers, @var{a} to @var{b}: eg (number-list 2 5) ->
;; (2 3 4 5)
;; 
(define (number-list a b)

  (let loop ((count a)
	     (acc '()))
    (cond
     ((> count b) (reverse acc))
     (else
      (loop (+ count 1) (cons count acc))))))



;; ls must be a list of strings, atom must be a string.
;; 
;; return either #t or #f.
;; 
(define (string-member? atom ls)
  
  (cond 
   ((null? ls) #f)
   ((string=? atom (car ls)) #t)
   (else
    (string-member? atom (cdr ls)))))

;; 
(define (member? atom ls)
  (cond 
   ((null? ls) #f)
   ((equal? atom (car ls)) #t)
   (else 
    (member? atom (cdr ls)))))
    

;; range: works like the eponymous python function
;; e.g. (range 3) -> '(0, 1, 2)
;; e.g. (range 1 3) -> '(1, 2)
;; 
(define (range first . second)

  (if (number? first)
		
      ;; OK input
      (let ((n-low/n-high (if (= (length second) 0)
			      (cons 0 first)
			      (cons first (car second)))))
	
	(let loop ((count (car n-low/n-high))
		   (rng '()))
	  (cond 
	   ((>= count (cdr n-low/n-high)) (reverse rng))
	   (else 
	    (loop (+ count 1)
		  (cons count rng))))))))

(define (file-n-lines file-name)
  (if (not (file-exists? file-name))
      #f
      (call-with-input-file file-name
	(lambda (port)
	  (let loop ((line (read-line port))
		     (n 0))
	    (cond
	     ((eof-object? line) n)
	     (else (loop (read-line port) (+ n 1)))))))))


;; code from thi <ttn at mingle.glug.org>
;; 
;; run command and put the output into a string and return
;; it. (c.f. @code{run-command/strings})
;; 
(define (shell-command-to-string cmd)
  (with-output-to-string
    (lambda ()
      (let ((in-port (open-input-pipe cmd)))
	(let loop ((line (read-line in-port 'concat)))
	  (or (eof-object? line)
	      (begin
		(display line)
		(loop (read-line in-port 'concat)))))))))

;; run @var{cmd} putting output to @var{file-name} and reading
;; commands data from the list of strings @var{data-list}.
;; 
(define (shell-command-to-file-with-data cmd file-name data-list)

;    (call-with-output-file file-name
;	(lambda (port)

  (let ((write-port (open-pipe cmd OPEN_WRITE)))
    
    (let f ((data-list data-list))
      
      (if (null? data-list)
	  
	  (close-port write-port)
	  
	  (begin
	    (format #t "writing ~s~%" (car data-list))
	    (write (car data-list) write-port)
	    (newline write-port)
	    (f (cdr data-list)))))))


;; Is the command @var{f} in the path?
;; return #t or #f.
;; 
;; Thank you for this, rixed at happyleptic.org
(define (command-in-path? f)
  (let ((path (string-split (getenv "PATH") #\:))
	(make-absolute (lambda (path) (string-append path "/" f))))
    (any file-exists? (map make-absolute path))))


;; Return #t or #f
;; 
(define (command-in-path-or-absolute? cmd)

  (if (command-in-path? cmd)
      #t
      (begin
	(if (not (string? cmd))
	    #f
	    (let ((l (string-length cmd)))
	      (if (= l 0)
		  #f
		  (file-exists? cmd)))))))

      
;; Where cmd is e.g. "refmac" 
;;       args is (list "HKLIN" "thing.mtz")
;;       log-file-name is "refmac.log"      
;;       data-list is (list "HEAD" "END")
;; 
;; Return the exist status e.g. 0 or 1.
;; 
(define (goosh-command cmd args data-list log-file-name screen-output-also?)
   
  (if (not (command-in-path-or-absolute? cmd))
      
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
                      ;; (format #t "exit status: ~s~%" status) silence
                      status)) ; return status

                  (begin
                    (if (eq? screen-output-also? #t)
                        (format #t ":~a~%" obj))
                    (format log-file-port "~a~%" obj)
                    (f (read-line output-port))))))))))

(define (ok-goosh-status? status)

  (if (not (number? status))
      #f
      (= 0 status)))

      

; example usage:
;(goosh-command "mtzdump" (list "HKLIN" "a.mtz") (list "HEAD" "END") "test.log" #t)

;; run commands from an input file.
;; 
(define (goosh-command-with-file-input cmd args input-file log-file-name)

  ;; we had a problem here (may 2004).  If command was not in path,
  ;; this function either hung or caused coot to immediately stop.
  ;; So now we test to see if cmd exists in the path first, by
  ;; running "which" on it.
  ;; 
  (if (string? cmd) 
      (let ((have-command? (command-in-path-or-absolute? cmd)))

	(if have-command? ; we *do* have the command

	    (with-output-to-file log-file-name
	      (lambda ()
		(with-input-from-file input-file
		  (lambda() 
		    (apply run (cons cmd args))))))
	    '()))
      '()))

;; Return the strings screen output of cmd or #f if command was not found
;; 
(define (run-command/strings cmd args data-list)
    
  (if (not (command-in-path-or-absolute? cmd))

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


;(define goosh-command-sans-log
;  (lambda (cmd args data-list)

;    (if (string? cmd) 
;	(let ((have-command? (run "which" cmd)))

;	  (if (= have-command? 0) ; we *do* have the command

;	      )))))

;; Old stack-overflowing version	      
;;; concatonate st-ls with intermediate tag-str (tag-str is often " ").
;;;
;(define (string-append-with-string str-ls tag-str)
    
;  (string-concatenate

;   (let f ((ls str-ls))
;     (cond
;      ((null? ls) '())
;      (else 
;       (cons (car ls)
;	     (cons tag-str
;		   (f (cdr ls)))))))))

;; Crude test to see of 2 floats are the same (more or less).
;; Used in a greg test after setting the atom position.
;; 
(define (close-float? x1 x2)
  
  (< (abs (- x1 x2)) 0.001))


;; if passed a string, return a string with no spaces,
;; else return #f.
;; 
(define (strip-spaces str)

  (if (not (string? str))
      #f
      (let f ((str str))
	(apply string-append (string-split str #\space)))))

;; "  53" -> "53", " " -> ""
(define (strip-leading-spaces str)
  (let ((l (string-length str)))
    (if (= l 0)
	str
	(if (not (string=? (substring str 0 1) " "))
	    str
	    (strip-leading-spaces (substring str 1))))))


;; Append strings with tag-str between them
;; 
(define (string-append-with-string str-ls tag-str)

  (let ((r-string ""))

    (for-each (lambda (s)
		(set! r-string (string-append r-string 
					      s tag-str)))
	      str-ls)
    r-string))


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

;; What is the extension of file-name?
;; 
;; "a.pdb" -> "pdb"
;; "" -> ""
;; 
(define file-name-extension
  (lambda (file-name)

    (let ((ls (split-after-char-last #\. file-name list)))
      
      (cond 
       ((string=? (car ls) "") "") ; no extension
       (else 
	(car (cdr ls)))))))

;; e.g. "a.pdb" -> "a-tmp.pdb"
;; 
(define add-tmp-extension-to
  (lambda (file-name)

    (if (not (string? file-name))
	"tmp"
	(let ((f (strip-extension file-name)))
	  (string-append f "-tmp." (file-name-extension file-name))))))
	  

;; Same function as strip-extension, different name, as per scsh, in fact.
;; 
(define file-name-sans-extension
  (lambda (s)
    (strip-extension s)))


;; /a/b.t -> b.t   d/e.ext -> e.ext
;; file-name-sans-path
;; 
(define (strip-path s)

  (let ((ls (split-after-char-last #\/ s list)))
    (car (cdr ls))))

;; does s start with a "/" ?
;; return #t or #f
;; 
(define (slash-start? s)

  (if (string? s)
      (if (> (string-length s) 0)
	  (string=? (substring s 0 1) "/")
	  #f)
      #f))

;; simple scheme functions to concat the strings in ls (ls must contain
;; only strings)
;; 
(define (string-concatenate ls)
  (apply string-append ls))


;; return a string that contains the date/time
;; e.g.	"2006-01-02_2216.03"
;; 
(define unique-date/time-str
  (lambda ()

    (let ((lt (localtime (current-time))))

      (strftime "%Y-%m-%d_%H%M.%S" lt))))
    


;; return a list that has only every-nth members; 
;; e.g. @code{(every-nth '(0 1 2 3 4 5 6 7 8) 2)} -> '(0 2 3 6 8)
;;      @code{(every-nth '(0 1 2 3 4 5 6 7 8) 3)} -> '(0 3 6)
;; 
;; @var{n} must be positive
;; 
(define (every-nth ls n)

    (reverse
     (let f ((res '())
	     (ls ls)	    
	     (count 0))
       (cond
	((null? ls) res)
	((= count 0)
	 (f (cons (car ls) res) (cdr ls) (- n 1)))
	(else 
	 (f res (cdr ls) (- count 1)))))))


;; now useful in testing
;; 
;; residue-atoms must be a list
(define (get-atom-from-residue atom-name residue-atoms alt-conf)
  (let loop ((residue-atoms residue-atoms))
    (cond
     ((null? residue-atoms) #f)
     ((and (string=? atom-name (car (car (car residue-atoms))))
	   (string=? alt-conf (car (cdr (car (car residue-atoms))))))
      (car residue-atoms))
     (else 
      (loop (cdr residue-atoms))))))

;;
(define (get-atom-from-spec imol atom-spec)
  (get-atom imol
	    (atom-spec->chain-id  atom-spec)
	    (atom-spec->res-no    atom-spec)
	    (atom-spec->ins-code  atom-spec)
	    (atom-spec->atom-name atom-spec)
	    (atom-spec->alt-loc   atom-spec)))

;; Return an atom info or #f (if atom not found).
;; 
(define (get-atom imol chain-id resno ins-code atom-name . alt-conf)
  
  ;; main line
  ;; 
  ;; (format #t "DEBUG::  get-atom is passed imol ~s: chain-id: ~s  resno: ~s  atom-name: ~s alt-conf: ~s ~%" imol chain-id resno atom-name alt-conf)

  (let ((res-info (residue-info imol chain-id resno ""))
	(alt-conf-internal (if (null? alt-conf)
			       ""
			       (car alt-conf))))
    (if (not res-info)
	#f
	(get-atom-from-residue atom-name res-info alt-conf-internal))))

;;
(define (residue-info-dialog-displayed?)
  (= (residue-info-dialog-is-displayed) 1))

;; multi-read-pdb reads all the files matching
;; @code{@emph{glob-pattern}} in
;; directory @code{@emph{dir}}.  Typical usage of this might be:
;; @code{(multi-read-pdb "a*.pdb" ".")}
;;  
(define (multi-read-pdb glob-pattern dir)

  (let ((mol-list
	 (map (lambda (file)
		(format #t "Reading ~s in ~s~%" file dir)
		(let ((full-path (append-dir-file dir file)))
		  (handle-read-draw-molecule-with-recentre full-path 0)))
       (glob glob-pattern dir))))
    (if (not (null? mol-list))
	(let ((last-model (car (reverse mol-list))))
	  (if (valid-model-molecule? last-model)
	      (apply set-rotation-centre (molecule-centre last-model)))))))


;; read-pdb-all reads all the "*.pdb" files in the current directory.
;; 
(define read-pdb-all
  (lambda ()
    (let ((recentre-status (recentre-on-read-pdb)))
      (set-recentre-on-read-pdb 0)
      (map read-pdb (glob "*.pdb" "."))
      (set-recentre-on-read-pdb recentre-status))))

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
	 ((or (char=? (car chars) #\space) 
	      (char=? (car chars) #\tab)
	      (char=? (car chars) #\newline))
	  (f (cdr chars)
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



;;; These 2 functions from Chart and are copyrighted by Paul Emsley.

;; In a laughable attempt to minimise system dependence.
;; 
(define (append-dir-file dir-name file-name)

  (if (not (string? dir-name))
      #f
      (if (not (string? file-name))
	  #f
	  (if (> (string-length dir-name) 0)
	      (string-append (directory-as-file-name dir-name) "/" file-name)
	      file-name))))

;; Similarly attempting to minimise system dependence.
;; 
(define (append-dir-dir dir-name sub-dir-name)

  (if (not (string? dir-name))
      #f
      (if (not (string? sub-dir-name))
	  #f
	  (string-append (directory-as-file-name dir-name) "/" sub-dir-name))))

;; the first elements are directories, the last is the file name
;; c.f. os.path.join
;; 
(define (join-dir-file path-elements)
  (if (null? path-elements)
      ""
      (let ((fn (car (reverse path-elements)))
	    (dirs (reverse (cdr (reverse path-elements)))))
	(let ((running ""))
	  (for-each (lambda (d)
		      (set! running (append-dir-dir running d)))
		    dirs)
	  ;; (format #t "dirs: ~s~%" dirs)
	  ;; (format #t "fn: ~s~%" fn)
	  (append-dir-file running fn)))))


;; remove any trailing /s
;; 
(define (directory-as-file-name dir)

  (if (= 0 (string-length dir))
      "."
      (let ((rchars (reverse (string->list dir))))
	(let f ((rls rchars))
	  (cond 
	   ((eq? #\/ (car rls)) (f (cdr rls)))
	   (else
	    (list->string (reverse rls))))))))


;; return #t or #f depending on if file-name (which must be a string)
;; is a directory.
;; 
(define (is-directory? file-name)

  (eq? (stat:type (stat file-name)) 'directory))


;; return the dir-name on success.
;; 
;; return #f if dir-name is a file or we can't do the mkdir.
;; 
(define (coot-mkdir dir-name)

  (if (file-exists? dir-name)
      (if (is-directory? dir-name)
	  dir-name ; often we will be returning this
	  #f)
      (mkdir dir-name)))
	    

;; The following functions from PLEAC (guile version thereof of course).
;; 
;; Return: a list of files in the given directory
;; 
(define (directory-files dir)

  (if (not (access? dir R_OK))
    '()
    (let ((p (opendir dir)))
      (do ((file (readdir p) (readdir p))
           (ls '()))
          ((eof-object? file) (closedir p) (reverse! ls))
        (set! ls (cons file ls))))))

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

;; return a list of file names that match pattern pat in directory dir.
;; 
(define (glob pat dir)
  (let ((rx (make-regexp (glob->regexp pat))))
    (filter (lambda (x) (regexp-exec rx x)) (directory-files dir))))

;; prepend the files with the directory dir
;; 
;; (note that the args are reversed to "natural" order (which is the
;; reverse of glob args)
;; 
(define (useful-glob dir pat)
   (map (lambda (x) (append-dir-file dir x)) (glob pat dir)))


;; return the view matrix (useful for molscript, perhaps).
;; 
(define view-matrix
  (lambda ()

    (map (lambda (row-number)
	   (map (lambda (column-number)
		  (get-view-matrix-element row-number column-number))
		(list 0 1 2)))
	 (list 0 1 2))))

;; return the view quaternion
;; 
(define view-quaternion
  (lambda ()
    
    (map (lambda(ele)
	   (get-view-quaternion-internal ele))
	 (list 0 1 2 3))))

;; Return the view number 
;; 
(define (add-view position quaternion zoom view-name)

  (apply add-view-raw 
	 (append position
		 quaternion
		 (list zoom view-name))))
	 

;; Convert a view matrix to a view quaternion to set Coot view internals.
;; 
(define matrix->quaternion
  (lambda (m00 m10 m20
	   m01 m11 m21
	   m02 m12 m22)

    ;; From an idea by "Christian" at euclidianspace.com.  The
    ;; rotation matrix is special orthogonal, so (1 + trace) > 0. So
    ;; we can simply do a sqrt on the 4 trace variants.  Most of the
    ;; code here is simply recovering the sign.
    
    ;; return x with the sign of y
    (define convert-sign 
      (lambda (x y)
	
	(cond
	 ((and (> x 0)
	       (> y 0)) x)
	 ((and (< x 0)
	       (> y 0)) (- x))
	 ((and (> x 0)
	       (< y 0)) (- x))
	 (else x))))

    (let ((pw (+ 1 m00 m11 m22))
	  (px (+ 1 m00 (- m11) (- m22)))
	  (py (+ 1 (- m00) m11 (- m22)))
	  (pz (+ 1 (- m00) (- m11) m22)))

      (let ((pr (map (lambda (v)
		       (let ((v1 (if (< v 0) 0 v)))
			 (/ (sqrt v1) 2)))
		     (list pw px py pz))))

	(append (map (lambda (v signed)
		       (convert-sign v signed))
		     (cdr pr) (list (- m21 m12)
				    (- m02 m20)
				    (- m10 m01)))
	      (list (car pr)))))))
  

; e.g (quat-convert 0.0347695872187614 0.773433089256287 0.632923781871796
;                   0.774806916713715 0.379149734973907 -0.505885183811188
;                  -0.631241261959076 0.507983148097992 -0.586078405380249)
; -> 
; (-0.55715757608 -0.694704711 -7.549694273e-4 0.45492890477)

;(matrix->quaternion 0.0347695 0.7734330  0.632923
;		    0.7748069 0.3791497 -0.50588
;		    -0.631241 0.5079831 -0.58607)

; (matrix->quaternion 0.034 0.77 0.63 0.774  0.379 -0.505 -0.631 0.507 -0.586)
; (matrix->quaternion 1 0 0 0 1 0 0 0 1)

;; Set the view matrix using matrix->quaternion. 
;; 
;; Useful for using a view matrix from another program, perhaps.
;; 
(define (set-view-matrix m00 m10 m20 m01 m11 m21 m02 m12 m22)

  (apply set-view-quaternion (matrix->quaternion m00 m10 m20
						 m01 m11 m21
						 m02 m12 m22)))


;; Miguel's molecular orientation axes
;;
(define miguels-axes
  (lambda ()

    (apply set-axis-orientation-matrix (apply append (view-matrix)))
    (set-axis-orientation-matrix-usage 1)))


;; Return the molecule centre as a list of 3 numbers.
;; 
;; Note: mol-cen could contain values less than -9999.  
;; 
(define (molecule-centre imol) 

  (map (lambda(iaxis) 
	 (molecule-centre-internal imol iaxis))
       (list 0 1 2)))
      
;; Move the centre of molecule number imol to the current screen centre
;; 
(define (move-molecule-to-screen-centre imol)

  (if (valid-model-molecule? imol)
      (apply translate-molecule-by 
	     (cons imol (map - (rotation-centre) (molecule-centre imol))))))

;; This is a short name for the above.
;; 
(define move-molecule-here move-molecule-to-screen-centre)

;; this is an americanism
;; 
(define move-molecule-to-screen-center move-molecule-to-screen-centre)

;; Return a nine-membered list of numbers.
;; 
(define identity-matrix
  (lambda ()

    '(1 0 0 0 1 0 0 0 1)))

;; e.g. (translation 'x 2)
;;  -> '(2 0 0)
;; Return: "scheme false" on error
;; 
(define (translation axis length)

  (cond 
   ((not (symbol? axis)) (format #t "incomprehensible axis argument: ~s~%" axis))
   ((not (number? length)) (format #t "incomprehensible length argument: ~s~%" length))
   ((eq? axis 'x) (list length 0 0))
   ((eq? axis 'y) (list 0 length 0))
   ((eq? axis 'z) (list 0 0 length))
   (else 
    (format #t "symbol axis: ~s incomprehensible~%" axis)
    #f)))

;; Rotate degrees about screen axis, where axis is either 'x, 'y or 'z.
;; 
(define (rotate-about-screen-axis axis degrees)

  ;; 
  (define deg-to-rad
    (lambda (degs)
      (/ (* 3.1415926 degs) 180.0)))
  
  ;; 
  (define simple-rotation-x
    (lambda (alpha) 
      (let ((cos-alpha (cos alpha))
	    (sin-alpha (sin alpha)))
	(list (list 1 0 0)
	      (list 0 cos-alpha (- sin-alpha))
	      (list 0 sin-alpha cos-alpha)))))

  ;; 
  (define simple-rotation-y
    (lambda (alpha) 
      (let ((cos-alpha (cos alpha))
	    (sin-alpha (sin alpha)))
	(list (list cos-alpha 0 sin-alpha)
	      (list 0 1 0)
	      (list (- sin-alpha) 0 cos-alpha)))))

  ;; 
  (define simple-rotation-z
    (lambda (alpha) ;; radians
      (let ((cos-alpha (cos alpha))
	    (sin-alpha (sin alpha)))
	(list (list cos-alpha (- sin-alpha) 0)
	      (list sin-alpha cos-alpha 0)
	      (0 0 1)))))

  (define vm (view-matrix))
  
  (define mult 
    (lambda (mat1 mat2)
      mat2))
  
  (cond
   ((not (symbol? axis)) (format #t "incomprehensible axis argument: ~s~%" axis))
   ((not (number? degrees)) (format #t "incomprehensible length argument: ~s~%" degrees))
   ((eq? axis 'x) (mult view-matrix (simple-rotation-x (deg-to-rad degrees))))
   ((eq? axis 'y) (mult view-matrix (simple-rotation-y (deg-to-rad degrees))))
   ((eq? axis 'z) (mult view-matrix (simple-rotation-z (deg-to-rad degrees))))
   (else
    (format #t "symbol axis: ~s incomprehensible~%" axis)
    #f)))

;; Support for old toggle functions.  (consider instead the raw
;; functions use the direct set_displayed functions).
;; 
(define (toggle-display-map imol idummy)
  (if (= (map-is-displayed imol) 0)
      (set-map-displayed imol 1)
      (set-map-displayed imol 0) ))

;; toggle the display of imol
;; 
(define (toggle-display-mol imol)
  (if (= (mol-is-displayed imol) 0)
      (set-mol-displayed imol 1)
      (set-mol-displayed imol 0)))

;; toggle the active state (clickability) of imol
;; 
(define (toggle-active-mol imol)
  (if (= (mol-is-active imol) 0)
      (set-mol-active imol 1)
      (set-mol-active imol 0)))


;; return a scheme representation of molecule imol, or #f if we can't
;; do it (imol is a map, say).
;; 
(define (scheme-representation imol)

  (if (not (valid-model-molecule? imol))
      #f
      (list 
       (map 
	(lambda (chain-id)
	  (list chain-id 
		(map (lambda (serial-number)
		       (let ((res-name (resname-from-serial-number imol chain-id serial-number))
			     (res-no   (seqnum-from-serial-number  imol chain-id serial-number))
			     (ins-code (insertion-code-from-serial-number imol chain-id serial-number)))
			 (let ((r-info (residue-info imol chain-id res-no ins-code)))
			   (list res-no ins-code res-name r-info))))
		     (range 0 (chain-n-residues chain-id imol)))))
	(chain-ids imol)))))

;; reorder chains
;; 
(define (reorder-chains imol)
  
  ;; reorder elements of chain-list: e.g.
  ;; 
  ;; chain-list: (list '("C" 'xx) '("A" 'xx) '("B" 'xx))
  ;; 
  (define (reorder-chains-in-model chain-list)
    (sort-list! chain-list (lambda (ele-1 ele-2) (string<? (car ele-1) (car ele-2)))))

  ;; main line
  (let* ((s-rep (scheme-representation imol)))
    (if (list? s-rep)
	(let ((sorted-rep (map reorder-chains-in-model s-rep)))
	  (clear-and-update-molecule imol sorted-rep)))))


;; transform a coordinates molecule by a coot-rtop (which is a SCM
;; expression of a clipper::RTop), i.e. a list of a 9-element list and
;; a 3 element list. e.g. (list (list 1 0 0 0 1 0 0 0 1) (list 4.5 0.4
;; 1.2)).
;; 
(define (transform-coords-molecule imol rtop)
  (apply transform-molecule-by imol
	 (apply append rtop)))



;; @code{(transform-map imol mat trans about-pt radius space-group cell)}
;; 
;; where space-group is a HM-symbol and cell is a list of 6
;; parameters, where the cell angles are in degrees.
;; 
;; or @code{(transform-map imol trans about-pt radius)} for a simple translation
;; 
;; or @code{(transform-map imol trans radius)} when using the default
;; rotation-centre as the about-pt
;; 
(define transform-map
  (lambda args

    (define tf
      (lambda (imol mat trans about-pt radius space-group cell)

	(format #t "DEBUG:: tf was passed imol: ~s, mat: ~s trans: ~s, about-pt: ~s, radius: ~s, space-group: ~s, cell: ~s~%"
		imol mat trans about-pt radius space-group cell)

	(transform-map-raw imol 
			   (list-ref mat 0) 
			   (list-ref mat 1) 
			   (list-ref mat 2) 
			   (list-ref mat 3) 
			   (list-ref mat 4) 
			   (list-ref mat 5) 
			   (list-ref mat 6) 
			   (list-ref mat 7) 
			   (list-ref mat 8) 
			   (list-ref trans 0)
			   (list-ref trans 1)
			   (list-ref trans 2)
			   (list-ref about-pt 0)
			   (list-ref about-pt 1)
			   (list-ref about-pt 2)
			   radius
			   space-group
			   (list-ref cell 0) 
			   (list-ref cell 1) 
			   (list-ref cell 2) 
			   (list-ref cell 3) 
			   (list-ref cell 4) 
			   (list-ref cell 5))))

    (format #t  "debug:: transform-map was passed ~s~%" args)

    ;; main line
    (cond 
     ((= (length args) 7)
      (tf (list-ref args 0)
	  (list-ref args 1)
	  (list-ref args 2)
	  (list-ref args 3)
	  (list-ref args 4)
	  (list-ref args 5)
	  (list-ref args 6)))
     ((= (length args) 5) ;; imol-map mat trans about-pt radius
      (let ((imol (car args)))
	(tf imol
	    (list-ref args 1)
	    (list-ref args 2)
	    (list-ref args 3)
	    (list-ref args 4)
	    (space-group imol)
	    (cell imol))))
     ((= (length args) 4) ; no matrix specified
      (let ((imol (car args)))
	(tf imol
	    (identity-matrix)
            (list
             (list-ref args 1)
             (list-ref args 2)
             (list-ref args 3))
            (rotation-centre)
            (list-ref (cell imol) 0)
	    (space-group imol)
	    (cell imol))))
     ((= (length args) 3) ; no matrix or about point specified
      (let ((imol (car args)))
	(tf imol
	    (identity-matrix)
	    (list-ref args 1)
	    (rotation-centre)
	    (list-ref args 2)
	    (space-group imol)
	    (cell imol))))
     (else
      (format #t "arguments to transform-map incomprehensible: args: ~s~%" args)
      #f))))

;; return the NCS master of the first molecule that has NCS.
;; 
;; return "" on fail to find an ncs chain
;; 
(define (get-first-ncs-master-chain) 

  (let ((r ""))
    (let loop ((mols (model-molecule-list)))
      (cond
       ((null? mols) "")
       (else
	(let ((ncs-masters (ncs-master-chains (car mols))))
	  (if (not (null? ncs-masters))
	      (if (list? ncs-masters)
		  (car ncs-masters)
		  (loop (cdr mols)))
	      (loop (cdr mols)))))))))
		

	 
     
;; Define a map transformation function that obeys Lapthorn's Law of
;; NCS Handling Programs
;; 
;; typical usage: (transform-map-using-lsq-matrix 1 "A" 10 30 0 "A" 10 30 2 (rotation-centre) 6)
;;
;; Remember, that now the about-pt is the "to" point, i.e. the maps are brought from 
;; somewhere else and generated about the about-pt.
;; 
(define (transform-map-using-lsq-matrix imol-ref ref-chain ref-resno-start ref-resno-end imol-mov mov-chain mov-resno-start mov-resno-end imol-map about-pt radius)

  (clear-lsq-matches)
  (add-lsq-match ref-resno-start ref-resno-end ref-chain 
		 mov-resno-start mov-resno-end mov-chain 1)
  (let ((space-group (symmetry-operators->xHM 
		      (symmetry-operators imol-ref)))
	(cell-params (cell imol-ref)))

    (format #t "debug:: transform-map-using-lsq-matrix: imol-ref ~s~%" imol-ref)
    (format #t "debug:: transform-map-using-lsq-matrix: imol-mov ~s~%" imol-mov)
    (format #t "debug:: transform-map-using-lsq-matrix: imol-map ~s~%" imol-map)
    (format #t "debug:: transform-map-using-lsq-matrix: symmetry-operators imol-ref ~s~%" (symmetry-operators imol-ref))
    (format #t "debug:: transform-map-using-lsq-matrix: space-group ~s~%" space-group)
    (format #t "debug:: transform-map-using-lsq-matrix: cell-params ~s~%" cell-params)
    
    (if (not (and space-group cell-params))
	(let ((message (format #f "Bad cell or symmetry ~s ~s for molecule ~s~%"
			       cell-params space-group imol-ref)))
	  message) ;; fix syntax
	(let ((rtop (apply-lsq-matches imol-ref imol-mov)))
	  (transform-map imol-map (car rtop) (car (cdr rtop)) about-pt radius space-group cell-params)))))


;; Make the imol-th map brighter.
;; 
(define (brighten-map imol scale-factor)

  (if (valid-map-molecule? imol)
      (let ((current-colour (map-colour-components imol)))
	(if (list? current-colour)
	    (apply set-map-colour imol
		   (map (lambda (v)
			  (let ((new-v (* scale-factor v)))
			    (cond 
			     ((< new-v 0.05) 0.05)
			     ((> new-v 1.0) 1.0)
			     (else new-v))))
			current-colour))
	    (format #t "bad non-list current-colour ~s~%"
		    current-colour))
	(graphics-draw))))


;; Make all maps brighter
;; 
(define (brighten-maps)
  (map (lambda (imap)
	 (brighten-map imap 1.25))
       (map-molecule-list)))

;; Make all maps darker.
;; 
(define (darken-maps)
  (map (lambda (imap)
	 (brighten-map imap 0.8))
       (map-molecule-list)))
    

;; Return a list of chain ids for given molecule number @var{imol}.
;; return empty list on error
;; 
(define (chain-ids imol)

  (let ((number-of-chains (n-chains imol)))
    
    (map (lambda (chain-no)
	   (chain-id imol chain-no))
	 (number-list 0 (- number-of-chains 1)))))


;; convert from interface name to schemish name
;; 
;; Return #t or #f. 
;; 
(define (is-solvent-chain? imol chain-id)

  (= (is-solvent-chain-p imol chain-id) 1))

;; schemey interface to eponymous scripting interface function.
;; Return scheme true or false
;; 
(define (valid-model-molecule? imol)

  (if (not (number? imol))
      #f
      (= (is-valid-model-molecule imol) 1)))

;; schemey interface to eponymous scripting interface function.
;; Return scheme true or false.
;; 
(define (valid-map-molecule? imol)
  (= (is-valid-map-molecule imol) 1))


;; convenience function (slightly less typing).  
;;
;; Return #t or #f
;; 
(define (valid-refinement-map?)
  (valid-map-molecule? (imol-refinement-map)))


;; schemey interface to shelx molecule test
;; 
;; Return #t or #f.
;; 
(define (shelx-molecule? imol)
  (= (is-shelx-molecule imol) 1))

;; schemey interface to the function that returns whether or not a map
;; is a difference map.  
;;
;; Return #t or #f.
;; 
(define (is-difference-map? imol-map)
  (if (not (valid-map-molecule? imol-map))
      #f
      (= (map-is-difference-map imol-map) 1)))
      

;; Does residue resno with insertion code ins-code of chain chain-id
;; and in molecule number imol exist?  
;; 
;; Return #t or #f.
;; 
(define (residue-exists? imol chain-id resno ins-code)

  (= 1 (does-residue-exist-p imol chain-id resno ins-code)))

;; Does the residue contain hetatoms? 
;; Return #t or #f.
;; 
(define (residue-has-hetatms? imol chain-id res-no ins-code)
  (= 1 (residue-has-hetatms imol chain-id res-no ins-code)))


;; Return a list of 3 float for the centre of mas of molecule number imol.
;;
;; on faiure return #f.
;;
(define (centre-of-mass imol)

  (let ((centre (centre-of-mass-string imol)))
    
    (if (eq? #f centre)
	(format #t "molecule number ~s is not valid~%" imol)
	(call-with-input-string centre (lambda (port) (read port))))))

;; Return as a list the occupancy temperature-factor element x y z coordinates
;; of the given atom.
;; (the x,y,z are in Cartesian Angstroms).
;; 
;; on error (e.g. atom not found) return #f
;; 
(define (atom-specs imol chain-id resno ins-code atom-name alt-conf)
  (atom-info-string imol chain-id resno ins-code atom-name alt-conf))

(define (atom-spec->string spec)

  (string-append
   (atom-spec->chain-id spec)
   " "
   (number->string (atom-spec->res-no spec))
   (atom-spec->ins-code spec)
   " "
   (atom-spec->atom-name spec)
   (let ((al (atom-spec->alt-loc spec)))
     (if (= (string-length al) 0)
	 ""
	 (string-append " " al)))))

(define (atom-spec->residue-spec atom-spec)
  (if (= (length atom-spec) 5)
      (list-head atom-spec 3)
      (list-head (cddr atom-spec) 3)))


;; return a guess at the map to be refined (usually called after
;; imol-refinement-map returns -1)
;; 
(define (guess-refinement-map)
  
  (let loop ((map-list (map-molecule-list)))
    (cond
     ((null? map-list) -1) ; failed to find a map
     ((= (map-is-difference-map (car map-list)) 0) (car map-list))
     (else 
      (loop (cdr map-list))))))


;; Ian Tickle says (as far as I can understand) that the target rmsd
;; should be 0.25 or thereabouts.  You can over-ride it now.
;; 
(define target-auto-weighting-value 1.0) ;; set to 1.0 so that it
					 ;; converges.  More thought
					 ;; needed.


;; Set the refinement weight (matrix) by iterating the refinement and
;; varying the weight until the chi squares (not including the
;; non-bonded terms) reach 1.0 =/- 10%.  It uses sphere refinement.
;; The refinement map must be set.  At the end show the new weight in
;; the status bar.  Seems to take about 5 rounds.
;; 
(define (auto-weight-for-refinement)

  ;; return a pair of the imol and a list of residue specs.
  ;; or #f if that is not possible
  ;; 
  (define (sphere-residues radius)
    (let ((active-atom (active-residue)))
        (if (not (list? active-atom))
	    (begin
	      (format #t "No active atom~%")
	      #f)
            (let* ((centred-residue (list-head (cdr active-atom) 3))
                   (imol (car active-atom))
                   (other-residues (residues-near-residue imol centred-residue radius))
                   (all-residues (if (list? other-residues)
                                     (cons centred-residue other-residues)
                                     (list centred-residue))))
	      (list imol all-residues)))))

  ;; the refinent function that is run and returns nice refinement
  ;; results.
  ;; 
  (define (refinement-func)
    (let ((sr (sphere-residues 3.5)))
      (if sr
	  (with-auto-accept
	   (apply refine-residues sr)))))

  ;; get rid of non-bonded chi-squared results from the input list ls.
  ;; 
  (define (no-non-bonded ls)
    (cond 
     ((null? ls) '())
     ((string=? (car (car ls)) "Non-bonded")
      (no-non-bonded (cdr ls)))
     (else 
      (cons (car ls) (no-non-bonded (cdr ls))))))
  

  ;; return #f or a number, which is the current overweighting of the
  ;; density terms.  (of course, when not overweighted, the geometric
  ;; chi squareds will be about 1.0).
  ;; 
  (define (weight-scale-from-refinement-results rr)

    (if (not (list? rr))
	#f
        (let ((results-inner (list-ref rr 2)))
          (if (null? results-inner)
              #f
              (let* ((nnb-list (no-non-bonded results-inner))
                     (chi-squares (map (lambda (x) (list-ref x 2)) nnb-list))
                     (n (length chi-squares))
                     (sum (apply + chi-squares)))
                (if (= n 0)
                    #f
                    (/ sum n)))))))

  
  ;; main body
  ;; 
  (let loop ((results (refinement-func))
	     (n-trials 0))
    (let ((av-rms-d (weight-scale-from-refinement-results results)))
      (format #t " av-rms-d: ~s~%" av-rms-d)

      (cond 
       ;; something bad happened?
       ((not (number? av-rms-d))
	#f)

       ;; failed to converge?
       ((> n-trials 20)
	(begin 
	  (format #t "Unconverged~%")
	  #f))

       ;; done?
       ((and  (< av-rms-d (* target-auto-weighting-value 1.1))
	      (> av-rms-d (* target-auto-weighting-value 0.9)))
		  ;; done
	(let ((s (string-append 
		  "Success: Set weight matrix to "
		  (number->string (matrix-state)))))
		    (add-status-bar-text s)))

       (else 
	;; more refinement required
	
	;; squared causes ringing, 
	;; as does 1.5.
	;; Simple is overdamped.
	;; 
	(let* ((current-weight (matrix-state))
	       (new-weight (/ (* target-auto-weighting-value current-weight) av-rms-d)))

	  (format #t "INFO:: ==== setting refinement weight to ~s from * ~s / ~s ~%" 
		   new-weight 
		   current-weight 
		   av-rms-d)

	  (if (< new-weight 2)
	      (begin 
		;; weight refinement not converging
		#f)
	  
	      (begin
		(set-matrix new-weight)
		(loop (refinement-func) (+ n-trials 1))))))))))



;; Print the sequence of molecule number @var{imol}
;; 
;; This is not really a util, perhaps it should be somewhere else?
;; 
(define (print-sequence imol)

  (map (lambda (chain)
	 (print-sequence-chain imol chain))
       (chain-ids imol))
  'done) ; return value


;; simple utility function to return the contents of a file as a string.
;; 
(define (pir-file-name->pir-sequence pir-file-name)
  (if (not (file-exists? pir-file-name))
      #f
      (call-with-input-file pir-file-name
	(lambda (port)
	  (let loop  ((lines '())
		      (line (read-line port)))
	    (cond
	     ((eof-object? line) 
	      (string-append-with-string (reverse lines) "\n"))
	     (else
	      (loop (cons line lines) (read-line port)))))))))

;; Associate the contents of a PIR file with a molecule.
;; 
(define (associate-pir-file imol chain-id pir-file-name)
  (let ((seq-text (pir-file-name->pir-sequence pir-file-name)))
    (if seq-text
	(assign-pir-sequence imol chain-id seq-text) ;; it isn't assigning sequence!
	(format #t "WARNING:: associate-pir-file: bad text for ~s~%" pir-file-name))))



;; comma key hook
(define graphics-comma-key-pressed-hook
  (lambda ()
    'empty))


;; dot key hook
(define graphics-dot-key-pressed-hook
  (lambda ()
    'empty))

;; a list of (code key name thunk) 
;; e.g. '(103 "g" "Goto Blob" (blob-under-pointer-to-screen-centre))
(define *key-bindings* (list))

;; Add a key binding 
;; 
;; with a given name, key (e.g. "x" or "S") and the function to run
;; (a thunk) when that key is pressed.
;; 
(define (add-key-binding name key thunk)
  (if (number? key)
      (set! *key-bindings* (cons (list key key name thunk) *key-bindings*)))
  (if (string? key)
      (let ((code (key-sym-code key)))
	(if (not (= code -1))
	    (set! *key-bindings* (cons (list code key name thunk) *key-bindings*))
	    (format #t "WARNING:: could not find key-sym-code for ~s~%" key)))))


;; general key press hook, not for public use.
;; 
(define (graphics-general-key-press-hook key)
  ;; (format #t "Key ~s was pressed~%" key)
  (let ((field (assoc key *key-bindings*)))
    ;; (format #t "Field: ~s~%" field)
    (if (not field)
        (begin
          ;; try python first
          (if (coot-has-python?)
              (run-python-command (string-append 
                                   "graphics_general_key_press_hook("
                                   (number->string key)
                                   ")")))
          (format #t "Key ~s not found in (scheme) key bindings~%" key))
	(begin
	  ((car (cdr (cdr (cdr field)))))))))

;; Function requested by Mark White.
;; 
;; read XtalView (and maybe other) .vu files and convert them into generic 
;; objects.  
;; 
;; Pass the filename and an object name e.g.
;; (read-vu-file "axes.vu" "axes")
;; 
;; Returns: nothing interesting.
;; 
(define read-vu-file
  (lambda (filename obj-name)


    (define colour-from-number
      (lambda (obj)
	
	(if (string? obj)
	    obj
	    (if (number? obj)
		(cond
		 ((= 0 obj) "white")
		 (else 
		  "red"))
		"white"))))
	    
    ; main body
    (call-with-input-file filename
      (lambda (port)

	(let ((n (new-generic-object-number obj-name)))
	
	(let loop ((o (read port))
		   (current-line '()))

	  (cond 
	   ((eof-object? o) 
	    (set-display-generic-object n 1))
	   ((= (length current-line) 6)
	    (let ((colour (colour-from-number o)))
	      (apply to-generic-object-add-line 
		     (append (list n colour 2) (reverse current-line))))
	    (loop (read port) '()))
	   (else 
	    (loop (read port) (cons o current-line))))))))))

;; residue-test-function is a function that takes 4 arguments, the
;; chain-id, resno, inscode and residue-serial-number (should it be
;; needed) and returns either #f or return something interesting
;; (e.g. text for a button label).
;;
;; Return a list of residues, each of which has a return value at the
;; start, ie. (list return-value chain-id res-no ins-code)
;; 
(define (residues-matching-criteria imol residue-test-func)

  ;; these specs are prefixed by the serial number
  (let loop ((molecule-residue-specs (all-residues-with-serial-numbers imol))
	     (matchers '()))

      (cond
       ((null? molecule-residue-specs)
	(reverse matchers))
       (else
	(let ((rs (cdr (car molecule-residue-specs))))
	  (if (residue-test-func (residue-spec->chain-id rs)
				 (residue-spec->res-no rs)
				 (residue-spec->ins-code rs)
				 (car (car molecule-residue-specs)))
	      (loop (cdr molecule-residue-specs) (cons rs matchers))
	      (loop (cdr molecule-residue-specs) matchers)))))))

;; now this is in the API
;;
;; Return residue specs for all residues in imol (each spec is preceeded by #t)
;; 
;; (define (all-residues imol)
;;  (residues-matching-criteria imol (lambda (chain-id resno ins-code serial) #t)))
;;
(define (all-residues imol)
  (let ((l (all-residues-with-serial-numbers imol)))
    (if (not (list? l))
	l
	(map cdr l))))

(define (all-residues-sans-water imol)
  (residues-matching-criteria imol (lambda (chain-id res-no ins-code serial)
				     (let ((rn (residue-name imol chain-id res-no ins-code)))
				       (not (string=? rn "HOH"))))))

;; Return a list of all the residues in the chain
;; 
(define (residues-in-chain imol chain-id-in)
  (residues-matching-criteria imol (lambda (chain-id resno ins-code serial) (string=? chain-id chain-id-in))))

  
;; Return a list of all residues that have alt confs: where a residue
;; is specified thusly: (list chain-id resno ins-code)
;; 
(define (residues-with-alt-confs imol)
  
  (residues-matching-criteria
   imol 
   (lambda (chain-id res-no ins-code res-serial-no)

     (let ((atom-ls (residue-info imol chain-id res-no ins-code)))
       
       ;; return #f if there are no atoms with alt-confs, else return
       ;; a list of the residue's spec (chain-id resno ins-code)
       ;; 
       (let g ((atom-ls atom-ls)
	       (r #f))
	 
	 (cond 
	  ((null? atom-ls) #f)
	  (else 
	   (let* ((atom (car atom-ls))
		  (compound-name (car atom))
		  (alt-conf-str (car (cdr compound-name))))
	     (if (string=? alt-conf-str "")
		 (g (cdr atom-ls) #f)
		 #t)))))))))

;; Return a list of all the altconfs in the residue. 
;; Typically this will return (list "") or (list "A" "B")
;; 
(define (residue-alt-confs imol chain-id res-no ins-code) 
  
  (reverse 
   (let f ((atom-ls (residue-info imol chain-id res-no ins-code))
	   (alt-confs '()))

     (cond 
      ((null? atom-ls) alt-confs)
      ((not (list? atom-ls)) alt-confs) ;; might be #f
      (else
       (let* ((atom (car atom-ls))
	      (compound-name (car atom))
	      (alt-conf-str (car (cdr compound-name))))
	 (if (string-member? alt-conf-str alt-confs)
	     (f (cdr atom-ls) alt-confs)
	     (f (cdr atom-ls) (cons alt-conf-str alt-confs)))))))))

;; not to be confused with residue-atom->atom-name (which uses the output of residue-info
;; (which is not atom specs)
;;
;; extraction function
(define (atom-spec->chain-id atom-spec)
  ;; atom-spec example (list "A" 7 "" " SG " "")
  (cond
   ((null? atom-spec) #f)
   ((= (length atom-spec) 5) (list-ref atom-spec 0))
   ((= (length atom-spec) 6) (list-ref atom-spec 1))
   (else
    #f)))

;; extraction function
(define (atom-spec->res-no atom-spec)
  ;; atom-spec example (list "A" 7 "" " SG " "")
  (cond
   ((null? atom-spec) #f)
   ((= (length atom-spec) 5) (list-ref atom-spec 1))
   ((= (length atom-spec) 6) (list-ref atom-spec 2))
   (else
    #f)))

;; extraction function
(define (atom-spec->ins-code atom-spec)
  ;; atom-spec example (list "A" 7 "" " SG " "")
  (cond
   ((null? atom-spec) #f)
   ((= (length atom-spec) 5) (list-ref atom-spec 2))
   ((= (length atom-spec) 6) (list-ref atom-spec 3))
   (else
    #f)))

;; extraction function
(define (atom-spec->atom-name atom-spec)
  ;; atom-spec example (list "A" 7 "" " SG " "")
  (cond
   ((null? atom-spec) #f)
   ((= (length atom-spec) 5) (list-ref atom-spec 3))
   ((= (length atom-spec) 6) (list-ref atom-spec 4))
   (else
    #f)))

;; extraction function
(define (atom-spec->alt-loc atom-spec)
  ;; atom-spec example (list "A" 7 "" " SG " "")
  (cond
   ((null? atom-spec) #f)
   ((= (length atom-spec) 5) (list-ref atom-spec 4))
   ((= (length atom-spec) 6) (list-ref atom-spec 5))
   (else
    #f)))


;; simple extraction function
(define (res-spec->chain-id res-spec)
  (cond 
   ((null? res-spec) #f)
   ((= (length res-spec) 4)
    (list-ref res-spec 1))
   ((= (length res-spec) 3)
    (list-ref res-spec 0))
   (else 
    #f)))
    
;; simple extraction function 
(define (res-spec->res-no res-spec)
  (cond 
   ((null? res-spec) #f)
   ((= (length res-spec) 4)
    (list-ref res-spec 2))
   ((= (length res-spec) 3)
    (list-ref res-spec 1))
   (else 
    #f)))

;; simple extraction function 
(define (res-spec->ins-code res-spec)
  (cond 
   ((null? res-spec) #f)
   ((= (length res-spec) 4)
    (list-ref res-spec 3))
   ((= (length res-spec) 3)
    (list-ref res-spec 2))
   (else 
    #f)))
    

;; Return #f if no atom can be found given the spec else return a list
;; consisting of the atom name and alt-conf specifier.  
;; 
;; Choose an atom that is called " CA ".  Failing that choose the
;; first atom.
;; 
(define (residue-spec->atom-for-centre imol chain-id res-no ins-code)

  ;; residue-info can return #f
  (let ((atom-ls (residue-info imol chain-id res-no ins-code)))

    (if (list? atom-ls)
	(let f ((atom-ls atom-ls)
		(centre-atom-name-alt-conf #f))

	  (cond 
	   ((null? atom-ls) centre-atom-name-alt-conf)
	   ((let* ((atom (car atom-ls))
		   (compound-name (car atom))
		   (atom-name (car compound-name))
		   (alt-conf-str (car (cdr compound-name))))
	      (if (eq? centre-atom-name-alt-conf #f)
		  (f (cdr atom-ls) (list atom-name alt-conf-str))
		  (if (string=? atom-name " CA ")
		      (f (cdr atom-ls) (list " CA " alt-conf-str))
		      (f (cdr atom-ls) centre-atom-name-alt-conf))))))))))

(define (set-go-to-atom res-spec)
  (set-go-to-atom-chain-residue-atom-name 
   (res-spec->chain-id res-spec)
   (res-spec->res-no   res-spec)
   " CA "))


;; 
;; 
(define (update-go-to-atom-from-current-atom)

  (let ((active-atom (active-residue)))
    (if active-atom
	(let ((imol      (list-ref active-atom 0))
	      (chain-id  (list-ref active-atom 1))
	      (resno     (list-ref active-atom 2))
	      (ins-code  (list-ref active-atom 3))
	      (atom-name (list-ref active-atom 4))
	      (alt-conf  (list-ref active-atom 5))
	      (go-to-atom-imol-current (go-to-atom-molecule-number)))
	  (set-go-to-atom-molecule imol)
	  ; (if (not (= imol go-to-atom-imol-current))
	  (update-go-to-atom-window-on-other-molecule-chosen imol)
	  (set-go-to-atom-chain-residue-atom-name chain-id
						  resno
						  atom-name)))))
;; 
;; 
(define (flip-active-ligand)
  (let ((active-atom (active-residue)))
    (if active-atom
	(let ((imol      (list-ref active-atom 0))
	      (chain-id  (list-ref active-atom 1))
	      (resno     (list-ref active-atom 2))
	      (ins-code  (list-ref active-atom 3))
	      (atom-name (list-ref active-atom 4))
	      (alt-conf  (list-ref active-atom 5)))
	  (flip-ligand imol chain-id resno)))))


;; Typically one might want to use this on a water, but it deletes the
;; nearest CA currently...  Needs a re-think.  Should active-atom just
;; return the nearest atom and not be clever about returning a CA.
;; 
(define (delete-atom-by-active-residue)

  (let ((active-atom (active-residue)))
    (if active-atom
	(apply delete-atom active-atom))))

;; 
(define (merge-solvent-chains imol)

  (if (valid-model-molecule? imol)
      (let ((solvent-chains (filter (lambda (c) (is-solvent-chain? imol c)) 
				    (chain-ids imol))))
	(if (> (length solvent-chains) 1)
	    (let* ((master-solvent-chain (car solvent-chains)))
	      (renumber-waters imol)
	      (for-each (lambda (c)
			  (let* ((master-n-solvent-residues (chain-n-residues master-solvent-chain imol))
				 (last-prev-water
				  (seqnum-from-serial-number imol master-solvent-chain
							     (- master-n-solvent-residues 1)))
				 (n-residues (chain-n-residues master-solvent-chain imol)))
			    (renumber-residue-range imol c 1 n-residues last-prev-water)
			    (change-chain-id imol c master-solvent-chain 1
					     (+ last-prev-water 1)
					     (+ last-prev-water n-residues))))
			(cdr solvent-chains)))))))
						 
				


;; general mutate
;; 
;; typically:
;; 
;; overlay PTY onto given TYR
;
;; delete speced TYR
;; 
;; merge molecules PTY molecule int molecule number imol
;; 
;; change resno of PTY to that of the speced TYR
;; 
;; change chain id of PTY to that of speced TYR
;;
;; change chain ids with residue range for the PTY
;; 
(define (mutate-by-overlap imol chain-id-in resno tlc)

  ;; residue is standard residues or phosphorylated version?
  ;;
  (define (is-amino-acid? imol chain-id res-no)

    (let ((rn (residue-name imol chain-id res-no "")))
      (if (not (string? rn))
	  #f
	  (or (string=? rn "ALA")
	      (string=? rn "ARG")
	      (string=? rn "ASN")
	      (string=? rn "ASP")
	      (string=? rn "CYS")
	      (string=? rn "GLY")
	      (string=? rn "GLU")
	      (string=? rn "GLN")
	      (string=? rn "PHE")
	      (string=? rn "HIS")
	      (string=? rn "ILE")
	      (string=? rn "LEU")
	      (string=? rn "LYS")
	      (string=? rn "MET")
	      (string=? rn "PRO")
	      (string=? rn "SER")
	      (string=? rn "TYR")
	      (string=? rn "THR")
	      (string=? rn "VAL")
	      (string=? rn "TRP")
	      (string=? rn "SEP")
	      (string=? rn "PTR")
	      (string=? rn "TPO")))))

  ;;
  (define (is-nucleotide? imol chain-id res-no)
    (let ((rn (residue-name imol chain-id res-no "")))
      (if (not (string? rn))
	  #f
	  (or (string=? rn "G")
              (string=? rn "A")
              (string=? rn "T")
              (string=? rn "U")
              (string=? rn "C")
              (string=? rn "DG")
              (string=? rn "DA")
              (string=? rn "DT")
              (string=? rn "DC")))))


  ;;
  (define (overlap-by-main-chain imol-mov chain-id-mov res-no-mov ins-code-mov
				 imol-ref chain-id-ref res-no-ref ins-code-ref)

    (format #t "debug:: in overlap-by-main-chain:: ---------------- imol-mov: ~s imol-ref: ~s~%" imol-mov imol-ref)
    (clear-lsq-matches)
    (for-each (lambda (atom-name)
		(add-lsq-atom-pair (list chain-id-ref res-no-ref ins-code-ref atom-name "")
				   (list chain-id-mov res-no-mov ins-code-mov atom-name "")))
	      (list " CA " " N  " " C  "))  ;; PDBv3 FIXME
    (apply-lsq-matches imol-ref imol-mov))

  (define (is-purine? res-name)
    (or (string=? res-name "G")
        (string=? res-name "A")
        (string=? res-name "DA")
        (string=? res-name "DG")))

  (define (is-pyrimidine? res-name)
    (or (string=? res-name "C")
        (string=? res-name "T")
        (string=? res-name "U")
        (string=? res-name "DC")
        (string=? res-name "DT")))

  (define (overlap-by-base imol-mov chain-id-mov res-no-mov ins-code-mov
                           imol-ref chain-id-ref res-no-ref ins-code-ref)

    (format #t "debug:: in overlap-by-base:: ---------------- imol-mov: ~s imol-ref: ~s~%" imol-mov imol-ref)
    (clear-lsq-matches)
    ;; G and A

    (let ((rn-1 (residue-name imol-mov chain-id-mov res-no-mov ins-code-mov))
          (rn-2 (residue-name imol-ref chain-id-ref res-no-ref ins-code-ref))
          (purine-set (list " N9 " " N7 " " C5 " " N1 " " N3 "))
          (pyrimidine-set (list " N1 " " C5 " " N3 "))
          (purine->pyrimidine-set (list " N1 " " C2 " " N3 "))
          (pyrimidine->purine-set (list " N9 " " C4 " " N5 ")))

      (let ((atom-list-1 '())
            (atom-list-2 '()))

        (if (is-purine? rn-1)
            (if (is-purine? rn-2)
                (begin
                  (set! atom-list-1 purine-set)
                  (set! atom-list-2 purine-set))))
    
        (if (is-pyrimidine? rn-1)
            (if (is-pyrimidine? rn-2)
                (begin
                  (set! atom-list-1 pyrimidine-set)
                  (set! atom-list-2 pyrimidine-set))))

        (if (is-pyrimidine? rn-1)
            (if (is-purine? rn-2)
                (begin
                  (set! atom-list-1 pyrimidine->purine-set)
                  (set! atom-list-2 purine->pyrimidine-set))))
        
        (if (is-purine? rn-1)
            (if (is-pyrimidine? rn-2)
                (begin
                  (set! atom-list-1 purine->pyrimidine-set)
                  (set! atom-list-2 pyrimidine->purine-set))))
        
        (for-each (lambda (atom-name-1 atom-name-2)
                    (add-lsq-atom-pair (list chain-id-ref res-no-ref ins-code-ref atom-name-1 "")
                                       (list chain-id-mov res-no-mov ins-code-mov atom-name-2 "")))
                  atom-list-1
                  atom-list-2)

        (format #t "applying matches~%~!")
        (apply-lsq-matches imol-ref imol-mov)
        (format #t "done matches~%~!"))))

  ;; get-monomer-and-dictionary, now we check to see if we have a
  ;; molecule already loaded that matches this residue, if we have,
  ;; then use it.
  ;;
  (define (get-monomer-and-dictionary tlc)
    
    ;; (get-monomer tlc))
    (let ((have-tlc-molecule #f))
      (for-each (lambda (imol)
		  (let ((nc (n-chains imol)))
		    (if (= nc 1)
			(let* ((ch-id (chain-id imol 0))
			       (nr (chain-n-residues ch-id imol)))
			  (if (= nr 1)
			      (let ((rn (resname-from-serial-number imol ch-id 0)))
				(if (string? rn)
				    (if (string=? rn tlc)
					(set! have-tlc-molecule imol)))))))))
		(model-molecule-list))
      (let ((have-dict-for-tlc (monomer-restraints tlc)))
	(if (or (not have-tlc-molecule) (not have-dict-for-tlc))
	    (begin
	      (get-monomer tlc))
	    (begin
	      (format #t "we have dict and model for tlc already~%")
	      have-tlc-molecule)))))

  ;;
  (define (mutate-it)
    (let ((imol-ligand (get-monomer-and-dictionary tlc)))
      (if (not (valid-model-molecule? imol-ligand))
	  (let ((s (string-append " Oops.  Failed to get monomer " tlc)))
	    (add-status-bar-text s))
	  (begin
	    (delete-residue-hydrogens imol-ligand "A" 1 "" "")
	    (delete-atom imol-ligand "A" 1 "" " OXT" "")
	    (if (and (is-amino-acid? imol-ligand "A" 1)
		     (is-amino-acid? imol chain-id-in resno))
		(overlap-by-main-chain imol-ligand "A" 1 "" imol chain-id-in resno "")
		(overlap-ligands imol-ligand imol chain-id-in resno))

            (if (is-nucleotide? imol chain-id-in resno)
                (if (residue-exists? imol chain-id-in (- resno 1) "")
                    (delete-atom imol-ligand "A" 1 "" " OP3" "")))

            (if (and (is-nucleotide? imol-ligand "A" 1)
                     (is-nucleotide? imol chain-id-in resno))
		(overlap-by-base imol-ligand "A" 1 "" imol chain-id-in resno "")
		(overlap-ligands imol-ligand imol chain-id-in resno))

            (if (not (is-nucleotide? imol-ligand "A" 1))
                (match-ligand-torsions imol-ligand imol chain-id-in resno))
	    (delete-residue imol chain-id-in resno "")
	    (let* ((new-chain-id-info (merge-molecules (list imol-ligand) imol))
		   (nov (format #t "DEBUG:: ------ new-chain-id-info: ~s~%" new-chain-id-info)))
	      (let ((merge-status (car new-chain-id-info)))
                ;; merge-status is sometimes a spec, sometimes a chain-id pair
		(format #t "DEBUG:: ------ merge-status: ~s~%" merge-status)
		(if (= merge-status 1)
		    (let* ((new-res-spec (car (car (cdr new-chain-id-info))))
			   (new-chain-id (residue-spec->chain-id new-res-spec)))
		      (format #t "debug:: ------ new-chain-id-info: ~s~%" new-chain-id-info)
		      (format #t "debug:: ------ new-chain-id: ~s~%" new-chain-id)
		      (format #t "debug:: ------ new-res-spec: ~s~%" new-res-spec)
		      (format #t "debug:: ------ change-residue-number to ~s~%" resno)

 		      (change-residue-number imol
                                             (residue-spec->chain-id new-res-spec)
                                             (residue-spec->res-no   new-res-spec)
                                             (residue-spec->ins-code new-res-spec)
                                             resno "")

                      ;; (change-residue-number imol new-chain-id 1 "" resno "")

		      (format #t "debug:: ------ chain ids : ~s ~s~%" new-chain-id chain-id-in)
                      
		      (if (not (string=? new-chain-id chain-id-in))
                          (change-chain-id imol new-chain-id chain-id-in 1 resno resno)) ;; 1 means "use range"

		      (let ((replacement-state (refinement-immediate-replacement-state))
			    (imol-map (imol-refinement-map)))
			(set-refinement-immediate-replacement 1)
			(if (= imol-map -1)
			    (regularize-zone imol chain-id-in resno resno "")
			    (let ((spin-atoms (list " P  " " O1P" " O2P" " O3P"))
				  (dir-atoms (cond 
					      ((string=? tlc "PTR") (list " CZ " " OH "))
					      ((string=? tlc "SEP") (list " CB " " OG "))
					      ((string=? tlc "TPO") (list " CB " " OG1"))
					      (else 
					       #f))))
			      ;; (refine-zone imol chain-id-in resno resno "")
			      (if dir-atoms
				  (spin-search imol-map imol chain-id-in resno "" 
					       dir-atoms spin-atoms))
			      ;; (refine-zone imol chain-id-in resno resno "")
			      ))
			(accept-regularizement)
			(set-refinement-immediate-replacement replacement-state))

		      (set-mol-displayed imol-ligand 0)
		      (set-mol-active imol-ligand 0)))))))))

  ;; main line
  ;; 
  ;; First, if there are multiple maps, force the user to choose one,
  ;; rather than continuing.
  ;; 
  (let ((imol-map (imol-refinement-map)))
    (if (= imol-map -1)
	(let ((map-mols (map-molecule-list)))
	  (if (> (length map-mols) 1)
	      (show-select-map-dialog)
	      (mutate-it)))
	(mutate-it))))

  
;; A bit of fun 
;; 
(define (phosphorylate-active-residue)

  (let ((active-atom (active-residue)))
	
    (if (list? active-atom)
	(let* ((imol     (list-ref active-atom 0))
	       (chain-id (list-ref active-atom 1))
	       (resno    (list-ref active-atom 2))
	       (inscode  (list-ref active-atom 3))
	       (res-name (residue-name imol chain-id resno inscode)))
	  (cond
	   ((string=? res-name "TYR") (mutate-by-overlap imol chain-id resno "PTR"))
	   ((string=? res-name "SER") (mutate-by-overlap imol chain-id resno "SEP"))
	   ((string=? res-name "THR") (mutate-by-overlap imol chain-id resno "TPO"))
	   (else 
	    (let ((s (string-append "Can't Phosphorylate residue of type" res-name)))
	      (info-dialog s))))))))


;; A function for Overlaying ligands.  The transformation is applied
;; to all the atoms of the molecule that contains the moving ligand.
;; 
(define (overlay-my-ligands imol-mov chain-id-mov resno-mov imol-ref chain-id-ref resno-ref)

  (let* ((imol-frag (new-molecule-by-atom-selection 
		     imol-mov (string-append "//" chain-id-mov "/" (number->string resno-mov))))
	 (rtop+i (overlap-ligands imol-frag imol-ref chain-id-ref resno-ref)))
    (set-mol-displayed imol-frag 0)
    (set-mol-active    imol-frag 0)
    (if (not rtop+i)
	(format #t "WARNING:: Failure to get RToperator in overlap-ligands~%")
	(transform-coords-molecule imol-mov (car rtop+i)))))



;; 
(define (label-all-CAs imol)

  (map (lambda (chain-id)
	 (if (not (is-solvent-chain? imol chain-id))
	     (let ((n-residues (chain-n-residues chain-id imol)))
	       
	       (for-each 
		(lambda (serial-number)
		  
		  (let ((res-name (resname-from-serial-number imol chain-id serial-number))
			(res-no   (seqnum-from-serial-number  imol chain-id serial-number))
			(ins-code (insertion-code-from-serial-number imol chain-id serial-number)))

		    (add-atom-label imol chain-id res-no " CA ")))
		(number-list 0 (- n-residues 1))))))
       (chain-ids imol))
  (graphics-draw))



;; 
(define (label-all-atoms-in-residue imol chain-id resno inscode)

  (let ((atom-list (residue-info imol chain-id resno inscode)))
    (if (list? atom-list)
	(begin
	  (for-each (lambda (atom-info)
		      (add-atom-label imol chain-id resno (car (car atom-info))))
		    atom-list)
	  (graphics-draw)))))


  

;; 
(define (label-all-active-residue-atoms)

  (let* ((active-atom (active-residue)))
    (if (list? active-atom)
	(let ((imol     (list-ref active-atom 0))
	      (chain-id (list-ref active-atom 1))
	      (resno    (list-ref active-atom 2))
	      (inscode  (list-ref active-atom 3)))

	  (label-all-atoms-in-residue imol chain-id resno inscode))))
  (graphics-draw))



;; Resets alt confs and occupancies of atoms in residue that have
;; orphan alt-loc attributes.
;; 
(define (sanitise-alt-confs atom-info atom-ls)

  ;; return a matching atom (name match) if it exists.  Else return #f
  (define (name-match? atom-1 atom-ls)
    (let* ((compound-name-1 (car atom-1))
	   (atom-name-1 (car compound-name-1)))
      (let loop ((atom-ls atom-ls)
		 (matchers '()))
	(cond 
	 ((null? atom-ls) matchers) ; return the matching atoms
	 (else 
	  (let* ((compound-name-2 (car (car atom-ls)))
		 (atom-name-2 (car compound-name-2)))

	    (if (string=? atom-name-1 atom-name-2)
		(loop (cdr atom-ls) (cons (car atom-ls) matchers))
		(loop (cdr atom-ls) matchers))))))))



  ;; main body
  (let ((imol     (list-ref atom-info 0))
	(chain-id (list-ref atom-info 1))
	(resno    (list-ref atom-info 2))
	(inscode  (list-ref atom-info 3))
	(atom-attribute-settings '())) ; add to this with set!
    (for-each
     (lambda (atom)
       
       (let* ((compound-name (car atom))
	      (atom-name (car compound-name))
	      (alt-conf (car (cdr compound-name))))
	 (if (not (string=? alt-conf ""))
	     (let ((matchers (name-match? atom atom-ls)))
	       (if (= (length matchers) 1)
		   (set! atom-attribute-settings 
			 (cons (list imol chain-id resno inscode atom-name alt-conf
				     "alt-conf" "")
			       atom-attribute-settings))
		   (set! atom-attribute-settings
			 (cons (list imol chain-id resno inscode atom-name alt-conf
				     "occ" (if (shelx-molecule? imol) 11.0 1.0))
			       atom-attribute-settings)))))))
     atom-ls)
    ; (format #t "DEBUG:: atom-attribute-settings: ~s~%" atom-attribute-settings)
    (if (not (null? atom-attribute-settings))
	(begin
	  (set-atom-attributes atom-attribute-settings)
	  (if (residue-info-dialog-displayed?)
	      (residue-info-dialog imol chain-id resno inscode))))))

;; 
(define (sanitise-alt-confs-in-residue imol chain-id resno inscode)

  (let ((atom-info (list imol chain-id resno inscode 'dummy 'dummy))
	(atom-ls (residue-info imol chain-id resno inscode)))
    (sanitise-alt-confs atom-info atom-ls)))

  

;; Resets alt confs and occupancies of atoms in residue that have
;; orphan alt-loc attributes.  Use the active-residue.
;; 
(define (sanitise-alt-confs-active-residue)
  (let* ((active-atom (active-residue)))
    (if (list? active-atom)
	(let ((imol     (list-ref active-atom 0))
	      (chain-id (list-ref active-atom 1))
	      (resno    (list-ref active-atom 2))
	      (inscode  (list-ref active-atom 3)))
	  
	  (let ((atom-ls (residue-info imol chain-id resno inscode)))
	    
	    (if (list? atom-ls) 
		(sanitise-alt-confs active-atom atom-ls)))))))

;; 
(define (print-molecule-names)

  (map (lambda (molecule-number)
	 (format #t "    ~s    ~s~%" molecule-number (molecule-name molecule-number)))
       (molecule-number-list)))

;;   
(define (save-dialog-positions-to-init-file)

  ;; return #f on failure to find ~/.coot
  (define (dump-positions-to-init-file positions)
    (let ((init-file (string-append (getenv "HOME") "/.coot")))
      (if (not (file-exists? init-file))
	  #f
	  (let ((port (open-file init-file "a")))
	    (format port "; ----------------------------------------------------------~%")
	    (format port "; the following were written by extraction from a state file~%")
	    (format port "; ----------------------------------------------------------~%")
	    (let loop ((positions positions))
	      (cond
	       ((null? positions) 'done)
	       (else 
		(format port "~a~%" (car positions))
		(loop (cdr positions)))))
	    (format port "; -------------------------~%")
	    (close port)))))

  ;; main line
  (save-state)

  ;; FYI, the graphics window is set using
  ;;
  ;; (set-graphics-window-size 643 500)
  ;; (set-graphics-window-position 0 1)
  ;; They are not dialogs

  (let ((state-file  "0-coot.state.scm"))
    (if (not (file-exists? state-file))
	(format #t "Ooops ~s does not exist~%" state-file)
	(call-with-input-file "0-coot.state.scm"
	  (lambda (port)
	    (let loop ((line (read-line port))
		       (positions '()))
	      (cond
	       ((eof-object? line) (dump-positions-to-init-file positions))
	       ((string-match "-dialog-position" line)
		(format #t " Adding dialog position: ~a~%" line)
		(loop (read-line port)
		      (cons line positions)))
	       ((string-match "set-go-to-atom-window-position" line)
		(format #t " Adding dialog position: ~a~%" line)
		(loop (read-line port)
		      (cons line positions)))
	       (else 
		(loop (read-line port) positions)))))))))


;; multiple maps of varying colour from a given map.
;; 
(define (multi-chicken imol . n-colours) 

  (define (rotate-colour-map col degress)
    (list (* (/ degress 360)) (cadr col) (- (caddr col) (* (/ degress 360)))))
  
  (let ((start 1.0)
	(stop 6.0)
	(initial-colour '(0.2 0.2 0.8))
	(colour-range 360))
    
    (if (valid-map-molecule? imol)
	(let ((n-col (if (null? n-colours) 10 (car n-colours)))
	    (sigma (map-sigma imol)))

	  (format #t "range n-col returns: ~s ~%" (range n-col))
	  
	  (for-each 
	   (lambda (icol)
	     (let* ((new-map (copy-molecule imol))
		    (frac (/ icol n-col))
		    (contour-level-sigma (+ start (* (- stop start) frac))))
	       (set-last-map-contour-level (* sigma contour-level-sigma))
	       (apply set-last-map-colour (rotate-colour-map initial-colour (* colour-range frac)))))
	   (range n-col))))))


;; simple enumeration
;; 
(define BALL_AND_STICK 2)


;; hilight-colour is specified in degrees (round the colour wheel -
;; starting at yellow (e.g. 230 is purple))
;; 
(define (hilight-binding-site imol centre-residue-spec hilight-colour radius)

  (if (valid-model-molecule? imol)
      
      (let ((other-residues (residues-near-residue imol centre-residue-spec radius))
	    (atom-sel-str (residue-spec->atom-selection-string centre-residue-spec)))
	
	(let ((imol-new (new-molecule-by-atom-selection imol atom-sel-str))
	      (bb-type 1)
	      (draw-hydrogens-flag (draw-hydrogens-state imol)))

	  (set-mol-active imol-new 0)
	  (set-show-environment-distances 1)
	  (set-molecule-bonds-colour-map-rotation imol-new hilight-colour)
	  (additional-representation-by-attributes imol-new 
						   (car centre-residue-spec)
						   (cadr centre-residue-spec)
						   (cadr centre-residue-spec)
						   (caddr centre-residue-spec)
						   BALL_AND_STICK
						   bb-type 0.14
						   draw-hydrogens-flag)
	  
	  (map (lambda (spec)
		 (additional-representation-by-attributes imol
							  (car spec)
							  (cadr spec)
							  (cadr spec)
							  (caddr spec)
							  BALL_AND_STICK
							  bb-type 0.14
							  draw-hydrogens-flag))
		 
	       other-residues)))))

;; Function based on Davis et al. (2007) Molprobity: all atom contacts
;; and structure validation for proteins and nucleic acids, Nucleic
;; Acids Research 35, W375-W383.  
;; 
;;    "RNA sugar puckers (C3'endo or C2'endo) is strongly correlated
;;    to the perpendicular distance between the following (3')
;;    phosphate and either the plane of the base or the C1'-N1/9
;;    glycosidic bond vector. [] .. a sugar pucker is very difficult
;;    to determine directly from the electron density at resolutions
;;    typical for RNAs."
;;
;; To paraphrase:
;; The distance of the plane of the base to the following phosphate
;; is highly correlated to the pucker of the ribose. 
;; 
;; An analysis of the structures in RNADB2005 shows that a critical
;; distance of 3.3A provides a partition function to separate C2' from
;; C3' endo puckering.  Not all ribose follow this rule.  There may be
;; some errors in the models comprising RNADB2005. So we check the
;; distance of the following phosphate to the plane of the ribose and
;; record the riboses that are inconsitent.  We also report puckers
;; that are not C2' or C3'.  The puckers are determined by the most
;; out-of-plane atom of the ribose (the rms deviation of the 4 atoms
;; in the plane is calculated, but not used to determine the
;; puckering atom).
;; 
(define (pukka-puckers? imol)

  (let ((residue-list '())
	(crit-d 3.0)) ;; Richardson's grup value to partition C2'-endo from C3'-endo

    ;; 
    (define (add-questionable r)
      (set! residue-list (cons r residue-list)))

    ;;
    (define (get-ribose-residue-atom-name imol residue-spec pucker-atom)
      (let ((r-info (apply residue-info (cons imol residue-spec)))
	    (t-pucker-atom (string-append (substring pucker-atom 0 3) "*")))
	(if 
	 (string-member? pucker-atom (map (lambda (at) (car (car at))) r-info))
	 pucker-atom
	 t-pucker-atom)))

    ;; main line
    (map (lambda (chain-id)
           (if (not (is-solvent-chain? imol chain-id))
               (let ((n-residues (chain-n-residues chain-id imol)))
                 
                 (for-each 
                  (lambda (serial-number)
                    
                    (let ((res-name (resname-from-serial-number imol chain-id serial-number))
                          (res-no   (seqnum-from-serial-number  imol chain-id serial-number))
                          (ins-code (insertion-code-from-serial-number imol chain-id 
                                                                       serial-number)))
                      (if (not (string=? res-name "HOH"))
                          
                          (let* ((residue-spec (list chain-id res-no ins-code))
                                 (pi (pucker-info imol residue-spec 1)))
                            (if pi
                                (if (list? pi)
                                    (if (= 4 (length pi))
					(let ((pucker-atom (car (cdr pi))))
					  ; (format #t " ~s ~s~%" residue-spec pi)
					  (if (and (> (abs (car pi)) crit-d)
						   (string=? pucker-atom " C2'"))
					      (add-questionable
					       (list pucker-atom residue-spec
						     "Inconsistent phosphate distance for C2' pucker")))
					  (if (and (< (abs (car pi)) crit-d)
						       (string=? pucker-atom" C3'"))
					      (add-questionable
					       (list pucker-atom residue-spec
						     "Inconsistent phosphate distance for C3' pucker")))
					  (if (not (or (string=? pucker-atom " C2'")
						       (string=? pucker-atom " C3'")))
					      (add-questionable
					       (list pucker-atom residue-spec
						     (string-append "Puckered atom: " pucker-atom))))))))))))

                  (range n-residues)))))
         (chain-ids imol))

    (if (= (length residue-list) 0)
	(info-dialog "No bad puckers.")
	(let ((buttons (map (lambda (residue-info)
			      (let* ((residue-spec (car (cdr residue-info)))
				     (pucker-atom (car residue-info))
				     (at-name (get-ribose-residue-atom-name imol residue-spec pucker-atom)))
				(list 
				 (string-append 
				  (list-ref residue-spec 0)
				  " "
				  (number->string (list-ref residue-spec 1))
				  (list-ref residue-spec 2)
				  ": "
				  (list-ref residue-info 2))
				 (lambda ()
				   (set-go-to-atom-molecule imol)
				   (set-go-to-atom-chain-residue-atom-name (car residue-spec) 
									   (car (cdr residue-spec)) 
									   at-name)))))
			    (reverse residue-list))))
	  (dialog-box-of-buttons "Non-pukka puckers"
				 (cons 370 250)
				 buttons
				 "  Close  ")))))

;; if you don't have mogul, set this to #f
(define *use-mogul* #t)
(if (not (command-in-path? "mogul"))
    (set! *use-mogul* #f))


;; Generate restraints from the residue at the centre of the screen
;; using PRODRG. Delete hydrogens from the residue because PRODRG has
;; anomalous hydrogens.
;; 
(define (prodrg-ify imol chain-id res-no ins-code)
  
  (let ((new-mol (new-molecule-by-atom-selection imol
		  (string-append "//" chain-id "/" (number->string res-no)))))

    (set-mol-active    new-mol 0)
    (set-mol-displayed new-mol 0)

    (let ((prodrg-dir "coot-ccp4")
	  (res-name (residue-name imol chain-id res-no ins-code)))

      (if res-name 
	  (begin
	    
	    (make-directory-maybe prodrg-dir)
	    (let ((prodrg-xyzin  (append-dir-file prodrg-dir "prodrg-in.pdb"))
		  (prodrg-xyzout (append-dir-file prodrg-dir
						  (string-append "prodrg-" res-name ".pdb")))
		  (prodrg-cif    (append-dir-file prodrg-dir 
						  (string-append "prodrg-" res-name ".cif")))
		  (prodrg-log    (append-dir-file prodrg-dir "prodrg.log")))
	      
	      (delete-residue-hydrogens new-mol chain-id res-no ins-code "")
	      (delete-residue-hydrogens imol    chain-id res-no ins-code "") ; otherwise they fly
	      (write-pdb-file new-mol prodrg-xyzin)
	      (close-molecule new-mol)
	      (let ((status 
		     (goosh-command "cprodrg" 
				    (list "XYZIN"  prodrg-xyzin
					  "XYZOUT" prodrg-xyzout
					  "LIBOUT" prodrg-cif)
				    (list "MINI PREP" "END")
				    prodrg-log #t)))
		;; (format #t "cprodrg gave status ~s~%" status)
		(if (number? status)
		    (if (= status 0)
			(begin
			  (read-cif-dictionary prodrg-cif)
			  (let ((imol-new (handle-read-draw-molecule-with-recentre prodrg-xyzout 0))
				(rn (residue-name imol chain-id res-no ins-code)))
			    ;; (move-molecule-here imol-new)
			    (with-auto-accept
			     (regularize-zone imol-new "" 1 1 ""))
			    (overlap-ligands imol-new imol chain-id res-no)
			    (match-ligand-torsions imol-new imol chain-id res-no) ;; broken
			    (overlap-ligands imol-new imol chain-id res-no)
			    (set-residue-name imol-new "" 1 "" rn)
			    (change-chain-id imol-new "" chain-id 1 1 1)
			    (renumber-residue-range imol-new chain-id 1 1 (- res-no 1))
			    (set-mol-displayed imol-new 0)
			    (set-mol-active    imol-new 0)
			    ; (set-mol-displayed imol 0)
			    ; (set-mol-active    imol 0)

			    ;; I don't think that replace-fragment is the right
			    ;; function because that does not copy across the hydrogen
			    ;; atoms - and we want those, probably

;			    (replace-fragment imol imol-new 
;					      (string-append "//"
;							     chain-id
;							     "/"
;							     (number->string res-no))))))))))))))
								     
			    (let* ((imol-replacing (add-ligand-delete-residue-copy-molecule
						    imol-new chain-id res-no imol chain-id res-no))
				   (col (get-molecule-bonds-colour-map-rotation imol))
				   (new-col (+ col 5)))
			      (set-molecule-bonds-colour-map-rotation imol-replacing new-col)
			      (set-mol-active imol-replacing 0)
			      (set-mol-displayed imol-replacing 0)
			      (graphics-draw)))))))))))))

			    







;; ---------- annotations ---------------------

;; 
(define (add-annotation-here text)
  (let ((rc (rotation-centre))
	(ann (cons text (rotation-centre))))

    (set! *annotations* (cons ann *annotations*))
    (apply place-text text (append (rotation-centre) (list 0)))
    (graphics-draw)))

;; 
(define (save-annotations file-name)
  (call-with-output-file file-name
    (lambda (port)
      (format port "~s~%" *annotations*))))

;; 
(define (load-annotations file-name)
  (if (file-exists? file-name)
      (begin
	(call-with-input-file file-name
	  (lambda (port)
	    (let ((ls (read port)))
	      (if (list? ls)
		  (begin
		    (set! *annotations* ls)
		    (for-each (lambda (ann)
				(apply place-text (append ann (list 0))))
			      *annotations*)))
	      (graphics-draw)))))))

;; Here we construct the url that contains the latest (pre)
;; release info adding in "pre-release" if this binary is a
;; pre-release.  args ends up as something like: ("-s"
;; "xxx/phone-home.scm" "pre-release" "binary"
;; "Linux-1386-fedora-10-python-gtk2" "command-line"
;; "/home/xx/coot/bin/coot")
;; 
(define (make-latest-version-url)
  (let ((build-type (coot-sys-build-type)))
    
    (string-append 
     "http://www.biop.ox.ac.uk/coot/software/binaries/"
     (if (pre-release?)
	 "pre-releases"
	 "releases")
     "/"
     "type-binary-"
     build-type
     "-latest.txt")))



;; Get the binary (i.e. the action that happens when the download
;; button is pressed).  This is run in a thread, so it can't do
;; any graphics stuff. 
;; 
;; return #t if tar file was successfully downloaded and untared
;; and #f if not.
;; 
(define (run-download-binary-curl revision version-string pending-install-in-place-func set-file-name-func)

  (define my-format 
    (lambda args
      (apply format (cons #t args))))

  (define (match-md5sums tar-file-name target-md5sum-file-name)
    (if (not (file-exists? tar-file-name))
	#f
	(if (not (file-exists? target-md5sum-file-name))
	    (begin
	      ;; (format #t "OOps! ~s does not exist" target-md5sum-file-name)
	      #f)
	    (let ((target-md5-string (get-target-md5-string target-md5sum-file-name))
		  (md5-string (get-md5sum-string tar-file-name)))
	      (if (not (string? target-md5-string))
		  (begin
		    ;; (format #t "OOps ~s is not a string~%" target-md5-string)
		    #f)
		  (if (not (string? md5-string))
		      (begin
			;; (format #t "OOps ~s is not a string~%" target-md5-string)
			#f)
		      (if (not (string=? target-md5-string md5-string))
			  (begin
			    ;;n(format #t "Oops: md5sums do not match ~s ~s.  Doing nothing~%"
			    ;; target-md5-string md5-string)
			    #f)
			  #t)))))))

  ;; return success status as a boolean
  ;;
  (define (install-coot-tar-file tar-file-name)
    (let ((prefix-dir (getenv "COOT_PREFIX")))
      (if (not (string? prefix-dir))
	  (begin
	    ;; (format #t "OOps could not get COOT_PREFIX~%")
	    #f)
	  (if (not (directory-is-modifiable? prefix-dir))
	      (begin
		;; (format #t "OOps directory ~s is not modifiable~%" prefix-dir)
		#f)
	      (let ((pending-dir (append-dir-file prefix-dir "pending-install")))
		(if (not (file-exists? pending-dir))
		    (mkdir pending-dir))
		(begin
		  (if (not (file-exists? pending-dir))
		      (begin
			;; (format #t "OOps could not create ~s~%" pending-dir)
			#f)
		      (let ((a-tar-file-name (absolutify tar-file-name)))
			;; with-working-directory 
			(let ((current-dir (getcwd)))
			  (chdir pending-dir)
			  (goosh-command "tar" (list "xzf" a-tar-file-name) '() "untar.log" #f)
			  (chdir current-dir))
			))))))))


  ;; return as a string, or #f
  (define (get-target-md5-string file-name)
    (if (not (file-exists? file-name))
	#f
	(call-with-input-file file-name
	  (lambda (port)
	    (symbol->string (read port))))))

  ;; return a string
  (define (get-md5sum-string file-name)
    (if (not (file-exists? file-name))
	#f
	(let* ((s (shell-command-to-string (string-append "md5sum " file-name)))
	       (s-bits (string->list-of-strings s)))
	  (format #t "get md5sum on ~s returns ~s~%" file-name s)
	  (car s-bits))))


  ;; main line
  ;; 
;  (format #t "INFO:: run-download-binary-curl.... with revision ~s with version-string ~s~%" 
;	     revision version-string)
  (let ((prefix (getenv "COOT_PREFIX")))
    (if (not (string? prefix))
	(begin
	  (my-format "OOps! Can't find COOT_PREFIX~%")
	  #f)
	(let* ((pre-release-flag (string-match "-pre" (coot-version)))
	       (ys "www.ysbl.york.ac.uk/~emsley/software/binaries/")
	       (binary-type (coot-sys-build-type))
	       (host-dir (cond 
			  ((string=? binary-type "Linux-i386-fedora-3") ys)
			  ((string=? binary-type "Linux-i386-fedora-3-python") ys)
			  ((string=? binary-type "Linux-i386-fedora-8-python-gtk2") ys)
			  ((string=? binary-type "Linux-i386-fedora-8-gtk2") ys)
			  ((string=? binary-type "Linux-i386-fedora-10-python-gtk2") ys)
			  ((string=? binary-type "Linux-i686-ubuntu-8.04.3") ys)
			  ((string=? binary-type "Linux-i686-ubuntu-8.04.3-python") ys)
			  (else 
			   "www.biop.ox.ac.uk/coot/software/binaries/")))
	       (tar-file-name (string-append version-string "-binary-" binary-type ".tar.gz"))
	       (url 
		(if (string-match "ysbl.york.ac.uk" host-dir)
		    (if pre-release-flag 
			(string-append "http://" host-dir "nightlies/pre-release/" tar-file-name)
			(string-append "http://" host-dir "stable/"     tar-file-name))
		    (if pre-release-flag 
			(string-append "http://" host-dir "pre-releases/" tar-file-name)
			(string-append "http://" host-dir "releases/"     tar-file-name))))
	       (md5-url (string-append url ".md5sum"))
	       (md5-tar-file-name (string-append tar-file-name ".md5sum")))
	  
	  (if (procedure? set-file-name-func)
	      (begin
		(set-file-name-func tar-file-name)))

	  ;; (my-format "md5sum url for curl: ~s~%" md5-url)

	  ;; (format #t  "DEBUG:: getting (md5sum) URL: ~s~%" md5-tar-file-name)
	  (coot-get-url-and-activate-curl-hook md5-url md5-tar-file-name 1)
	  (format #t  "INFO:: getting URL: ~s~%" url)
	  (let ((curl-get-url-status (coot-get-url-and-activate-curl-hook url tar-file-name 1)))
	    ;; (format #t  "DEBUG:: done getting of URL: ~s~%" url)
	    ;; (format #t  "DEBUG:: now setting pending-install-in-place-func to full~%")

	    ;; OK, was the downloading of the binary successful?  If
	    ;; so, go on to match md5sums etc
	    ;; 
	    (if (= 0 curl-get-url-status)
		(begin
		  (pending-install-in-place-func 'full))
		

		(if (not (file-exists? tar-file-name))
		    (begin
		      ;; (format #t "Ooops: ~s does not exist after attempted download~%" tar-file-name)
		      #f)
		    (if (not (file-exists? md5-tar-file-name))
			(begin
			  ;; (format #t "Ooops: ~s does not exist after attempted download~%" 
			  ;; md5-tar-file-name)
			  #f)
			(if (not (match-md5sums tar-file-name md5-tar-file-name))
			    #f 
			    (let ((success (install-coot-tar-file tar-file-name)))
			      (if success 
				  (begin
				    (pending-install-in-place-func #t)
				    #t)
				  (begin
				    ;; (format #t "Ooops: untar of ~s failed~%" tar-file-name)
				    #f))))))))))))

;; get revision number from string
;;
;; (used in downloading new version)
;; 
(define (get-revision-from-string str)
  ;; e.g. str is "coot-0.6-pre-1-revision-2060" (with a newline at the
  ;; end too).  We want to return 2060 (a number) from here (or #f).
  (if (not (string? str))
      #f
      (if (= (string-length str) 0)
	  #f
	  (let* ((s (sans-final-newline str))
		 (ls (separate-fields-discarding-char #\- s list)))
	    (string->number (car (reverse ls)))))))

;; first generate a version string with no trailing newline.
;;
;; e.g. input:  "coot-0.6.2-pre-1-revision-2765\n"
;;      output: "coot-0.6.2-pre-1-revision-2765"
;;  
(define (coot-split-version-string str)
  (let* ((ls (split-before-char #\c str list))
	 (ls-2 (split-before-char #\" (car (reverse ls)) list))
	 (ls-3 (split-before-char #\newline (car ls-2) list)))

;     (format #t "notify-of-new-version str: ~s~%" str)
;     (format #t "ls: ~s~%" ls)
;     (format #t "ls-2: ~s~%" ls-2)
;     (format #t "ls-3: ~s~%" ls-3)
     (car ls-3)))


(define (file->string file-name)

  (if (not (file-exists? file-name))
      #f
      (call-with-output-string 
       (lambda (sport)
	 (call-with-input-file file-name
	   (lambda (port)
	     (let loop ((obj (read-line port)))
	       (if (not (eof-object? obj))
		   (begin
		     (display obj sport)
		     (newline sport)
		     (loop (read-line port)))))))))))

      
;; If "default.seq" (a simple text file with the sequence (not PIR or
;; FASTA)) exists in the current directory, then try to assign it to
;; each chain of each molecule.

;; In the first case the sequence is assigned to the closest match
;; (model sequence to target sequence), subsequently only chains
;; without a sequence associated with them are candidates for
;; matching.  The protein sequence has to have at least 95% sequence
;; identity with the target sequence in "default.seq"
;; 
(define (load-default-sequence)

  (let ((default-seq "default.seq"))
    (if (file-exists? default-seq)
	(let ((s (file->string default-seq)))
	  (align-to-closest-chain s 0.95)))))




;; a thread handling function
;;
;; not really for public manipulation.
;; 
(define (coot-updates-error-handler key . args)
  (format #t "error: finding updates: error in ~s with args ~s~%" key args))


;; update self 
;;
;; keep a copy of the old directories around in a directory named
;; after expiration time.
;; 
(define (update-self)
  (let* ((url (make-latest-version-url))
	 (version-string (coot-split-version-string (coot-get-url-as-string url)))
	 (revision (get-revision-from-string version-string))
	 (file-name-for-progress-bar #f)
	 (pending-install-in-place #f))

    (define (set-file-name-func file-name)
      (set! file-name-for-progress-bar file-name))

    (define (pending-install-in-place-func val)
      (set! pending-install-in-place val))

    (let ((continue-status #t))
      (call-with-new-thread
       (lambda()
	 (run-download-binary-curl revision version-string
				   pending-install-in-place-func
				   set-file-name-func)
	 (set! continue-status #f))
       coot-updates-error-handler)

      (while continue-status
	     (if (string? file-name-for-progress-bar)
		 (let ((curl-info (curl-progress-info file-name-for-progress-bar)))
		   ;; (format #t "curl-info: ~s~%" curl-info)
		   (if curl-info
		       (let ((v1 (assoc 'content-length-download curl-info))
			     (v2 (assoc 'size-download           curl-info)))
			 (if (list v1)  ;; Do I mean "list?"?
			     (if (list v2)
				 (let ((f (/ (cdr v2) (cdr v1))))
				   (format #t "~3,2f% " (* f 100))
				   (if (> f 0.999) 
				       (set! continue-status #f)))))))))
	     (sleep 2))))
  (coot-real-exit 0))

;; Invert the chiral centre of the atom we are centred on.  
;; If not centred on a chiral atom, then give a dialog.
;; 
;; The restraints for this monomer type are copied and renamed
;; (changing comp-id, 3-letter-code and name too).  The monomer is
;; regularized.  Chiral Hydrogen (if needed) is enabled now in the
;; minimizer.
;;
;; This should almost all be c++ code so that Bernie doesn't have to
;; redo it.  This is temporary then.
;; 
(define (chiral-centre-inverter-scm-version)
  
  ;; return "LI2" given "LIG", 
  ;;        "LI3" given "LI2"
  ;;        "L10" given "L11"  (etc.)
  ;; 
  ;; return #f if current-res-name is not a string of 1, 2 or 3
  ;; characters.
  ;; 
  ;; 20140321 This used to return #f when current-res-name was 1XX.
  ;; Now it returned LIG, which is not great, but better than it was.
  ;; 
  (define (new-res-name current-res-name)

    (define (generate-name current-res-name)
      (format #t "----- generate new name from current-res-name: ~s~%" current-res-name)
      (if (not (string? current-res-name))
	  #f
	  (cond
	   ((= (string-length current-res-name) 3)
	    (let ((one   (substring current-res-name 0 1))
		  (two   (substring current-res-name 1 2))
		  (three (substring current-res-name 2 3)))
	      (let ((n-1 (string->number one))
		    (n-2 (string->number two))
		    (n-3 (string->number three)))

		;; LLL LLN LNN
		;;
		;; others LNL NNN NNL NLN NLL 

		(cond
		 ;; "LIG" -> "LI2", "LI9" -> "L10"
		 ;; 
		 ((and (not (number? n-1))
		       (not (number? n-2))
		       (number? n-3))
		  (if (< n-3 9)
		      (string-append one two (number->string (+ n-3 1)))
		      (string-append one "10")))
		 
		 ;; "L10" -> "L11",  "L19" -> "L20"
		 ((and (not (number? n-1))
		       (number? n-2)
		       (number? n-3))
		  (if (< n-3 9)
		      (string-append one two (number->string (+ n-3 1)))
		      (if (< n-2 9)
			  (string-append one (number->string (+ n-2 1)) "0")
			  ;; "L99" -> #f
			  #f)))

		 ((and (not (number? n-1))
		       (number? n-2)
		       (not (number? n-3)))
		  "LIG")

		 ((and (number? n-1)
		       (number? n-2)
		       (number? n-3))
		  "LIG")

		 ((and (number? n-1)
		       (number? n-2)
		       (not (number? n-3)))
		  "LIG")

		 ((and (number? n-1)
		       (not (number? n-2))
		       (number? n-3))
		  "LIG")

		 ((and (number? n-1)
		       (not (number? n-2))
		       (not (number? n-3)))
		  "LIG")

		 ;; LIG -> LI2
		 ((and (not (number? n-1))
		       (not (number? n-2))
		       (not (number? n-3)))
		  
		  (string-append one two "2"))

		 (else #f)))))

	   (else 
	    (string-append current-res-name "2")))))

    ;; check for name collision, if yes, then try next.
    ;; 
    (let f ((s (generate-name current-res-name)))
      (if (eq? s #f) ;; e.g. current-res-name was "L99"
	  #f
	  (let ((r (monomer-restraints s)))
	    (if (not r)
		s
		(f (generate-name s)))))))
  

  ;; 
  (define (is-chiral-centre-atom? atom-name chirals)

    (let loop ((chiral-restraints (cdr chirals)))
      (cond
       ((null? chiral-restraints) #f)
       ((string=? (list-ref (car chiral-restraints) 1) atom-name) #t)
       (else
	(loop (cdr chiral-restraints))))))


  ;; main line
  ;; 
  (using-active-atom 
   (let* ((res-name (residue-name aa-imol aa-chain-id aa-res-no aa-ins-code))
	  (restraints (monomer-restraints res-name)))

     (if (list? restraints)
	 (begin
	   (let ((chirals (assoc "_chem_comp_chir" restraints)))
	     (if (not chirals)
		 (info-dialog "Not chiral restraints for residue")
		 (if (not (is-chiral-centre-atom? aa-atom-name chirals))
		     (info-dialog "Not a chiral centre")
		     (begin

		       ;; switch the chiral centre
		       (for-each
			(lambda (chiral-index)
			  (let ((chiral-atom (list-ref (list-ref (cdr chirals) chiral-index) 1)))
			    (if (string=? chiral-atom aa-atom-name)

				(let* ((current-sign (list-ref (list-ref (cdr chirals) chiral-index) 5))
				       (new-sign
					(cond 
					 ((=  1 current-sign) -1)
					 ((= -1 current-sign)  1)
					 (else 
					  'unset))))
				  (format #t "current sign: ~s~%" current-sign)
				  (format #t "chiral atom  pre manip ~s~%"
					  (list-ref (cdr chirals) chiral-index))
				  
				  (if (not (eq? new-sign 'unset))
				      (list-set! (list-ref (cdr chirals) chiral-index) 5 new-sign))
				  
				  (format #t "chiral atom post manip ~s~%"
					  (list-ref (cdr chirals) chiral-index))
				  ))))
			(range (length (cdr chirals))))
		       
		       (format #t "post manip chirals ~s~%" 
			       (assoc "_chem_comp_chir" restraints))

		       
		       (let ((new-residue-name (new-res-name res-name)))

			 ;; Now new-res-name should not return #f, so
			 ;; this is belt and braces protection.
			 ;; 
			 (if (not (string? new-residue-name))
			     (begin
			       (format #t "ERROR:: Ooops new-residue-name: ~s not a string ~%" 
				       new-residue-name)
			       #f)

			     ;; Happy Path
			     ;; 
			     (let ((new-mol (new-molecule-by-atom-selection 
					     aa-imol (string-append "//"
								    aa-chain-id
								    "/"
								    (number->string aa-res-no)

								    aa-ins-code))))
			       ;; switch the comp-id three-letter-code and name
			       ;; 
			       (let ((chem-comp (assoc "_chem_comp" restraints)))
				 (if chem-comp
				     (begin
				       (list-set! (cdr chem-comp) 0 new-residue-name) ;; comp_id
				       (list-set! (cdr chem-comp) 1 new-residue-name) ;; 3-letter-code
				       (list-set! (cdr chem-comp) 2 new-residue-name) ;; name
				       )))

			       (format #t "restraints: ~s~%" restraints)
			       (format #t "new-residue-name: ~s~%" new-residue-name)
			       (format #t " ------------ pre set-monomer-restraints ~%")
			       (set-monomer-restraints new-residue-name restraints)
			       (format #t " ------------ post set-monomer-restraints ~%")
			       (set-residue-name new-mol aa-chain-id aa-res-no aa-ins-code new-residue-name)
			       (regularize-zone new-mol aa-chain-id aa-res-no aa-res-no aa-alt-conf)
			       ))))))))))))

(define (residue-is-close-to-screen-centre? imol chain-id res-no ins-code)
  (define (square x) (* x x))
  (let ((rc (residue-centre imol chain-id res-no ins-code)))
    (if (not (list? rc))
	#f
	(let ((sc (rotation-centre)))
	  (let ((dist-sum-sq (apply + (map square (map - rc sc)))))
	    (< dist-sum-sq 25))))))


 
;; return a string "DB00xxxx.mol" or some such - this is the file name
;; of the mdl mol file from drugbank. Or #f/undefined on fail.  
;; Test result with string?.
;; 
;; If you are debugging this function, you probably want to be looking
;; at handle-rev-string or drugbox->drugitem.
;; 
(define (get-drug-via-wikipedia drug-name-in)

  ;; To Bernie: I imagine that this would be confusing.  I don't think
  ;; that this is the way a clever scheme xml programmer would do
  ;; it... My advice is to ignore the code here for the most part and
  ;; make a function that (only) externally acts as this does.
  ;; i.e. goes to wikipedia, downloads the xml file - following the
  ;; redirect if needed and then construct the url for the mol file
  ;; from drugbank and get that.


  (define (drugbox->drugitem key s) 
    ;; (format #t "===== here in drugbox->drugitem : ~s ~s ================== ~%" key s)
    ;; (format #t "==========================================================~%" )
    (let ((n (string-length s)))
      (let ((ls (string-split s #\newline)))
	(let loop ((ls ls))
	  (cond
	   ((null? ls) #f) ;; no drugbank reference found
	   ((string-match key (car ls))
	    (let ((sls (string-split (car ls) #\space)))
	      (let ((n-sls (length sls)))
		(let ((r (list-ref sls (- n-sls 1))))
		  (if (string=? r "")
		      #f
		      r)))))
	   (else 
	    (loop (cdr ls))))))))

  (define (drugbox->drugbank s)
    (drugbox->drugitem "DrugBank *=" s))

  (define (drugbox->pubchem s)
    (let ((i (drugbox->drugitem "PubChem *=" s)))
      (format #t "in drugbox->pubchem i: ~s~%" i)
      i))

  (define (drugbox->chemspider s)
    (let ((i (drugbox->drugitem "ChemSpiderID *=" s)))
      (format #t "in drugbox->chemspider i: ~s~%" i)
      i))

  (define (last-element list-of-strings)
    (car (reverse list-of-strings)))

  ;; With some clever coding, these handle-***-value functions could
  ;; be consolidated.  There is likely something clever in Python to
  ;; process XML files. 
  ;; 
  ;; In fact, now that I read the returned XML file more carefully,
  ;; you probably don't need to process XML at all! Just split to
  ;; strings and look for "DrugBank = " and take the last of the words
  ;; on that line.

  (define (handle-rev-string-old rev-string)

    (format #t "in handle-revi-string rev-string: ~s~%~!" rev-string)

    (let ((open-match (string-match "#[Rr][Ee][Dd][Ii][Rr][Ee][Cc][Tt] \\[\\[" rev-string))
	  (close-match (string-match "\\]\\]" rev-string))
	  (re (make-regexp "[{][{][Dd]rugbox|[{][{][Cc]hembox")))

      (cond 
       ((and open-match close-match)
	(let ((s (substring rev-string 12 (car (vector-ref close-match 1)))))
	  (get-drug-via-wikipedia s)))
      
       ((let ((re-result (regexp-exec re rev-string)))
	  (if re-result
	      (begin
		;; (format #t "INFO:: matched a Drugbox~%")
		(let ((dbi (drugbox->drugbank rev-string)))

		  (if (not (string? dbi))
		      (begin 

			(format #t "DEBUG:: dbi (drugbank result) not a string: ~s~%" dbi)
			;; try pubchem as a fallback
			(let ((pc (drugbox->pubchem rev-string)))
			  (if (not (string? pc))

			      ;; OK, try chemspider extraction
			      (let ((cs (drugbox->chemspider rev-string)))
				(if (not (string? cs))
				    (begin 
				      (format #t "not a string ~s ~s~%~!" cs)
				      #f)

				    ;; chemspider extraction worked
				    (let ((cs-mol-url
					   (string-append "http://www.chemspider.com/"
							  "FilesHandler.ashx?type=str&striph=yes&id="
							  cs))
					  (file-name (string-append "cs-" cs ".mol")))
				      (coot-get-url cs-mol-url file-name)
				      file-name)))

			      ;; pubchem extraction worked
			      (let ((pc-mol-url

				     ;; Old style
				     ;; (string-append "http://pubchem.ncbi.nlm.nih.gov"
				     ;;                 "/summary/summary.cgi?cid=" 
				     ;;                 pc 
				     ;;                 "&disopt=DisplaySDF")

				     ;; chasing new style:
				     (string-append "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
						    pc 
						    "/record/SDF/?record_type=2d&response_type=display")
				     )
				    (file-name (string-append "pc-" pc ".mol")))
				(format #t "========== pubchem pc-mol-url: ~s~%" pc-mol-url)
				(coot-get-url pc-mol-url file-name)
				file-name))))

		      ;; normal path hopefully
		      ;; 
		      (let ((db-mol-uri (string-append 
					 "http://www.drugbank.ca/structures/structures/small_molecule_drugs/" 
					 dbi ".mol"))
			    (file-name (string-append dbi ".mol")))
			(format #t "getting drugbank url: ~s ~s~%" db-mol-uri file-name)
			(coot-get-url db-mol-uri file-name)
			file-name)))))))

       (else 
	(format #t "handle-rev-string none of the above~%")
	#f))))

  ;;
  (define (file-seems-good? file-name)
    (if (not (file-exists? file-name))
	#f
	(let* ((stat-result (stat file-name))
	       (stat-size (stat:size stat-result)))
	  (> stat-size 20))))

  ;; return a mol file name
  (define (handle-rev-string-2016 rev-string)

    (let ((lines (string-split rev-string #\newline)))

      (cond 
       ((= (length lines) 1)
	(let* ((line (car lines))
	       (open-match (string-match "#[Rr][Ee][Dd][Ii][Rr][Ee][Cc][Tt] \\[\\[" line)))
	  (if open-match
	      (let ((close-match (string-match "\\]\\]" line)))
		(if close-match
		    (begin
; 				  (format #t "-----------:  open-match ~%" open-match)
; 				  (format #t "-----------: close-match ~%" close-match)
; 				  (format #t "-----------: line ~%" line)
		      (let ((s (substring line 12 (car (vector-ref close-match 1)))))
			(get-drug-via-wikipedia s)))))))) ;; returns a file anme

       (else
	(let ((db-id-list '())) ;; drugbank-id e.g. DB01098
	  (for-each (lambda (line)
		      (format #t "debug:: line: ~s~%" line)

		      ;; we don't want to hit xxx_Ref - hence the trailing space

		      (if (string-match "DrugBank[ \t]" line)
                (let ((parts (string->list-of-strings line)))
                  (format #t "   debug:: drugbank parts: ~s~%" parts)
                  (let ((id-string (last-element parts)))
                    (if (number? (string->number id-string))
                        (set! db-id-list (cons (cons "DrugBank" id-string) db-id-list))))))

		      (if (string-match "ChemSpiderID[ \t]" line)
                (let ((parts (string->list-of-strings line)))
                  (format #t "   debug:: ChemSpiderID parts: ~s~%" parts)
                  (let ((id-string (last-element parts)))
                    (if (not (string-match "correct" id-string))
                        (if (number? (string->number id-string))
                            (set! db-id-list (cons (cons "ChemSpider" id-string) db-id-list)))))))

		      (if (string-match "PubChem[ \t]" line)
                (let ((parts (string->list-of-strings line)))
                  (format #t "   debug:: PubChem parts: ~s~%" parts)
                  (let ((id-string (last-element parts)))
                    (if (not (string-match "correct" id-string))
                        (set! db-id-list (cons (cons "PubChem" id-string) db-id-list))))))

		      (if (string-match "ChEMBL[ \t]" line)
                (let ((parts (string->list-of-strings line)))
                  (format #t "   debug:: ChEMBL parts: ~s~%" parts)
                  (let ((id-string (last-element parts)))
                    (if (number? (string->number id-string))
                        (set! db-id-list (cons (cons "ChEMBL" id-string) db-id-list))))))
		      )
		    lines)

	  (format #t "DEBUG:: handle-rev-string-2016: db-id-list: ~s~%" db-id-list)

	  ;; now db-id-list is something like (("DrugBank" . 12234) ("ChEMBL" . 6789))
	  ;; can we find one of them that works?

	  (let loop ((db-id-list db-id-list))

	    (cond
	     ((null? db-id-list) "Failed-to-find-a-molecule-file-name")
	     ((string=? (car (car db-id-list)) "DrugBank")

	       (let ((db-id (car db-id-list)))
		 (if (pair? db-id)

		     (let ((DBWebsite (car db-id))
			   (id (cdr db-id)))

		       (let ((db-mol-uri (string-append
					  "https://www.drugbank.ca/structures/small_molecule_drugs/" id ".mol"))
			     (file-name (string-append "drugbank-" id ".mol")))
			 (format #t "DEBUG:: DrugBank path: getting url: ~s to file ~s~%"
				 db-mol-uri file-name)
			 (coot-get-url db-mol-uri file-name)
			 ;; check that file-name is good here
			 (if (file-seems-good? file-name)
			     (format #t "DEBUG:: yes db file-name: ~s seems good ~%" file-name))
			 (if (file-seems-good? file-name)
			     file-name
			     (loop (cdr db-id-list))))))))

	      ((string=? (car (car db-id-list))  "ChemSpider")
	       (let ((db-id (car db-id-list)))
		 (if (pair? db-id)

		     (let ((DBWebsite (car db-id))
			   (id (cdr db-id)))

		       (let ((cs-mol-url
			      (string-append "http://www.chemspider.com/"
					     "FilesHandler.ashx?type=str&striph=yes&id="
					     id))
			     (file-name (string-append "cs-" id ".mol")))
			 (format #t "Get ---ChemSpider--- ~s ~s~%" cs-mol-url id)
			 (coot-get-url cs-mol-url file-name)
			 (if (file-seems-good? file-name)
			     (format #t "DEBUG:: yes cs file-name: ~s seems good ~%" file-name))
			 (if (file-seems-good? file-name)
			     file-name
			     (loop (cdr db-id-list))))))))

	      ((string=? (car (car db-id-list))  "ChEMBL")
	       (let ((db-id (car db-id-list)))
		 (if (pair? db-id)

		     (let ((DBWebsite (car db-id))
			   (id (cdr db-id)))

		       (let ((mol-url
			      (string-append "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL"
					     id ".sdf"))
			     (file-name (string-append "chembl-" id ".sdf")))
			 (coot-get-url mol-url file-name)
			 (if (file-seems-good? file-name)
			     (format #t "DEBUG:: yes chembl file-name: ~s seems good ~%" file-name))
			 (if (file-seems-good? file-name)
			     file-name
			     (loop (cdr db-id-list))))))))

	      ((string=? (car (car db-id-list))  "PubChem")
	       (let ((db-id (car db-id-list)))
		 (if (pair? db-id)

		     (let ((DBWebsite (car db-id))
			   (id (cdr db-id)))

		       (let ((pc-mol-url (string-append "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
							id "/record/SDF/?record_type=2d&response_type=display"))
			     (file-name (string-append "pc-" id ".mol")))
			 (format #t "========== pubchem pc-mol-url: ~s~%" pc-mol-url)
			 (coot-get-url pc-mol-url file-name)
			 (if (file-seems-good? file-name)
			     (format #t "DEBUG:: yes pc file-name: ~s seems good ~%" file-name))
			 (if (file-seems-good? file-name)
			     file-name
			     (loop (cdr db-id-list))))))))

	      (else
	       (loop (cdr db-id-list))))))))))

  (define (handle-sxml-rev-value sxml)
    ;; (format #t "handle-sxml-rev-value sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond 
       ((null? ls) #f)
       ((not (list? ls)) #f) ;; shouldn't happen
       (else 
	(let ((entity (car ls)))
	  (if (string? entity)
	      (let ((aaa (handle-rev-string-2016 entity)))
		(format #t "---------- aaa: ~s~%" aaa)
		aaa)
	      (loop (cdr ls))))))))

  (define (handle-sxml-revisions-value sxml)
    ;; (format #t "handle-sxml-revisions-value sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond
       ((null? ls) #f)
       ((let ((entity (car ls)))
	  (if (not (list? entity))
	      #f
	      (let ((l (length entity)))
		(if (= l 0)
		    #f
		    (let ((first (car entity)))
		      (if (not (eq? first 'rev))
			  #f
			  (let ((rev (cdr entity)))
			    (handle-sxml-rev-value rev))))))))))))

  (define (handle-sxml-page-value sxml)
    ;; (format #t "handle-sxml-page-value sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond
       ((null? ls) #f)
       ((let ((entity (car ls)))
	  (if (not (list? entity))
	      (loop (cdr ls))
	      (let ((l (length entity)))
		(if (not (> l 0))
		    (loop (cdr ls))
		    (let ((first (car entity)))
		      (if (not (eq? first 'revisions))
			  (loop (cdr ls))
			  (let ((revisions (cdr entity)))
			    (handle-sxml-revisions-value revisions))))))))))))
	      
  (define (handle-sxml-pages-value sxml)
    ;; (format #t "handle-sxml-pages-value sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond
       ((null? ls) #f)
       ((let ((entity (car ls)))
	  (if (not (list? entity))
	      (loop (cdr ls))
	      (let ((l (length entity)))
		(if (not (> l 0))
		    (loop (cdr ls))
		    (let ((first (car entity)))
		      (if (not (eq? first 'page))
			  (loop (cdr ls))
			  (let ((page (cdr entity)))
			    (handle-sxml-page-value page))))))))))))


  (define (handle-sxml-query-value sxml)
    ;; (format #t "handle-sxml-query-value sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond
       ((null? ls) #f)
       ((let ((entity (car ls)))
	  (if (not (list? entity))
	      (loop (cdr ls))
	      (let ((l (length entity)))
		(if (not (> l 0))
		    (loop (cdr ls))
		    (let ((first (car entity)))
		      (if (not (eq? first 'pages))
			  (loop (cdr ls))
			  (let ((pages (cdr entity)))
			    (handle-sxml-pages-value pages))))))))))))
			      

  (define (handle-sxml-api-value sxml)
    ;; (format #t "handle-sxml-api-value sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond
       ((null? ls) #f)
       ((let ((entity (car ls)))
	  (if (not (list? entity))
	      (loop (cdr ls))
	      (let ((l (length entity)))
		(if (not (> l 0))
		    (loop (cdr ls))
		    (let ((first (car entity)))
		      (if (not (eq? first 'query))
			  (loop (cdr ls))
			  (let ((query-value (cdr entity)))
			    (handle-sxml-query-value query-value))))))))))))

  (define (handle-sxml sxml)
    ;; (format #t "handle-sxml sxml: ~s~%" sxml)
    (let loop ((ls sxml))
      (cond
       ((null? ls) #f)
       ((not (list? ls)) #f)
       ((let ((entity (car ls)))
	  (if (not (list? entity))
	      (loop (cdr ls))
	      (let ((l (length entity)))
		(if (not (> l 0))
		    (loop (cdr ls))
		    (let ((first (car entity)))
		      (if (not (eq? first 'api))
			  (loop (cdr ls))
			  (let ((value (cdr entity)))
			    (handle-sxml-api-value value))))))))))))


  ;; main line
  ;;
  (if (string? drug-name-in)
      (if (> (string-length drug-name-in) 0)
	  ;; we need to downcase the drug name for wikipedia
	  (let* ((drug-name (string-downcase drug-name-in))
		 (url (string-append 
		       ;; "http://en.wikipedia.org/wiki/" 
		       "https://en.wikipedia.org/w/api.php?format=xml&action=query&titles="
		       drug-name
		       "&prop=revisions&rvprop=content"
		       ))
		 (xml (coot-get-url-as-string url)))

	    (format #t "INFO:: get-drug-via-wikipedia: url: ~s~%" url)
	    (let ((l (string-length xml)))

	      (if (= l 0)
		  ;; something bad happened - wikipedia replied with an empty string
		  (begin
		    (format #t "INFO:: Badness - wikipedia replied with an empty string~%"))

		  (begin
		    (call-with-output-file (string-append drug-name ".xml")
		      (lambda (port)
			(display xml port)))

		    ;;
		    ;; (format #t "about to get captured ~%")
		    (let ((captured (call-with-input-string
				     xml (lambda (string-port)
					   ;; (format #t "about to get sxml ~%") paracetamol.xml fails here
					   (let ((sxml (xml->sxml string-port)))
					     (handle-sxml sxml))))))
		      (format #t "------------- captured: ~s~%" captured)
		      captured))))))))


;; --------------------------------------------------------


(define (get-SMILES-for-comp-id-from-pdbe comp-id)

  (if (not (string? comp-id))
      #f
      (let ((s (SMILES-for-comp-id comp-id)))
	(if (string? s)
	    s
	    (let ((file-name (get-pdbe-cif-for-comp-id comp-id)))
	      (if (string? file-name)
		  (begin
		    (read-cif-dictionary file-name)
		    (let ((s (SMILES-for-comp-id comp-id)))
		      s))))))))


;; return #f or a file-name
;;
(define (get-pdbe-cif-for-comp-id comp-id)

  (define (wrapped-fetch url cif-file-name)
	;; use the network then 
	(let ((nov (format #t "getting url: ~s~%" url))
	      (state (coot-get-url url cif-file-name)))
	  (if (not (= state 0))
	      (begin
		(let ((s (apply string-append (list "Problem downloading\n"
						    url "\n to file \n"
						    cif-file-name
						    "."
						    ))))
		  (info-dialog s)
		  #f))
	      
	      ;; it worked
	      cif-file-name)))
		

  ;; main line of get-pdbe-cif-for-comp-id
  ;;
  (let ((download-dir (get-directory "coot-download")))
    (let ((cif-file-name (append-dir-file download-dir
					(string-append "PDBe-" comp-id ".cif")))
	(url (string-append
	      "ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/files/mmcif/"
	      comp-id
	      ".cif")))

    (if (file-exists? cif-file-name)
	;; try the file system cache
	(begin
	  (let* ((stat-data (stat cif-file-name))
		 (l (stat:size stat-data)))
	    (if (> l 0)
		cif-file-name
		;; overwrite the empty file
		(wrapped-fetch url cif-file-name))))

	;; go get it then
	(wrapped-fetch url cif-file-name)))))



(define (template-keybindings-to-preferences)

  ;; (let* ((pkg-data-dir
  ;;         (if (file-exists? (pkgdatadir))
  ;;             (pkgdatadir)
  ;;             (append-dir-dir (append-dir-dir (getenv "COOT_PREFIX") "share") "coot"))))

  (let ((pkg-data-dir
         (let ((coot-prefix-dir (getenv "COOT_PREFIX"))) ;; try this first
           (if (string? coot-prefix-dir)
               (append-dir-dir (append-dir-dir coot-prefix-dir "share") "coot")
               (pkgdatadir))))) ;; self-install directory

    (let* ((bindings-file-name "template-key-bindings.scm")
	 (scm-dir (append-dir-dir pkg-data-dir "scheme"))
	 (ref-scm (append-dir-file scm-dir bindings-file-name)))
    (if (not (string? ref-scm))
	(begin
	  (add-status-bar-text "bad get-pkgdata-dir."))
	(begin
	  (if (not (file-exists? ref-scm))
	      (begin
		(add-status-bar-text "Missing reference template key bindings."))
	      (begin
		;; happy path
		(let ((home (getenv "HOME")))
		  (if (string? home)
		      (let* ((pref-dir  (append-dir-dir home ".coot-preferences"))
			     (pref-file (append-dir-file pref-dir bindings-file-name)))
			;; don't install it if it is already in place.
			(if (file-exists? pref-file)
			    (begin 
			      (let ((s (string-append "keybinding file "
						      pref-file
						      " already exists.  Not overwritten.")))
				(add-status-bar-text s)))
			    (begin
			      ;; check the directory first
			      (if (not (file-exists? pref-dir))
				  (make-directory-maybe pref-dir))
			      (copy-file ref-scm pref-file)
			      (if (file-exists? pref-file)
				  (load pref-file))))))))))))))



;; something like this for intermediate atoms also?
;;
;; n-neighbs is either 0 or 1
;;
(define (rebuild-residues-using-db-loop imol middle-residue-spec n-neighbs)

  ;; utility function
  ;;
  (define (remove-any-GLY-CBs imol-db-loop saved-residue-names ch-id resno-low)

    (for-each (lambda (residue-idx res-name)
		(if (string? res-name) ;; might be #f for missing residues
		    (if (string=? res-name "GLY")
			(let ((res-no-gly (+ resno-low residue-idx)))
			  (delete-atom imol-db-loop ch-id res-no-gly "" " CB " "")))))
	      (range (length saved-residue-names)) saved-residue-names))

   ;; main line

   (let* ((resno-mid (residue-spec->res-no middle-residue-spec))
	  (resno-low  (- resno-mid n-neighbs))
	  (resno-high (+ resno-mid n-neighbs))
          (r (append (range (- resno-mid 4 n-neighbs) (+ resno-mid n-neighbs 0))
		     (range (+ resno-mid 1 n-neighbs) (+ resno-mid 3 n-neighbs))))
          (nov (print-var r))
          (ch-id (residue-spec->chain-id middle-residue-spec))
          (residue-specs (map (lambda(res-no)
                                (list ch-id res-no ""))
                              r)))

     (let ((loop-mols (protein-db-loops
		       imol residue-specs (imol-refinement-map) 1 *db-loop-preserve-residue-names*)))
        (let ((residue-spec-of-residues-to-be-replaced (if (= n-neighbs 0)
                                                           (list middle-residue-spec)
							   (map (lambda (r)
								  (list ch-id r ""))
								(range resno-low
								       (+ resno-high 1))))))

          (let ((saved-residue-names (map (lambda (r) (residue-spec->residue-name imol r))
					  residue-spec-of-residues-to-be-replaced)))
	    (if (> (length loop-mols) 0)

		;; loop-mols have the correct residue numbering but the wrong chain-id and
		;; residue type

		(let ((imol-db-loop (car (list-ref loop-mols 1)))
		      (tmp-loop-mols (car loop-mols)))
		  (let ((chain-id-db-loop (chain-id imol-db-loop 0)))
		    (if (not (string=? ch-id chain-id-db-loop))
			(change-chain-id imol-db-loop chain-id-db-loop ch-id 0 0 0))
		    (let ((selection (string-append "//" ch-id "/"
						    (number->string resno-low)
						    "-"
						    (number->string resno-high)
						    )))

		      ;; if the original residue was a GLY, the imol-db-loop might (probably
		      ;; will) have CB. If that is the case, then we should remove the CBs now
		      ;;
		      (remove-any-GLY-CBs imol-db-loop saved-residue-names ch-id resno-low)

		      ;; this moves atoms, adds atoms if needed, doesn't change the residue name
		      ;;
		      (replace-fragment imol imol-db-loop selection)
		      ;; tidy up
		      (for-each (lambda (i) (close-molecule i)) tmp-loop-mols))))))))))


	
;; Americans...
;; (define chiral-center-inverter chiral-centre-inverter)

;; to determine if we have pygtk
;; 
(define (coot-has-pygtk?)
  (if (coot-has-python?)
      ;; sometimes coot can start up (in development environment with
      ;; broken embeddded python), in which case coot_has_pygtk()
      ;; fails to return a value.  So try to account for that here.
      (let ((cap-var (run-python-command "coot_has_pygtk()")))
	(if (and (not (eq? cap-var #t))
		 (not (eq? cap-var #f)))
	    #f
	    cap-var))
      #f))




(define (pdb-validate accession-code imol)

  (let* ((down-id (string-downcase accession-code))
	 (python-string (string-append
			 "pdb_validate(\""
			 down-id "\"," (number->string imol)
			 ")")))
    (format #t "--------- run this python-string: ~s~%" python-string)
    (run-python-command python-string)))



(define (go-to-box-middle)

  (let ((ls (map-molecule-list)))
    (if (not (null? ls))
       (let ((imol-map (car ls)))
         (let ((c (cell imol-map)))
           (format #t "c ~s~%" c)
           (apply set-rotation-centre (map (lambda (a) (* a 0.5)) (list-head c 3))))))))




(define (write-current-sequence-as-pir imol chain-id file-name)
  (print-sequence-chain-general imol chain-id 1 1 file-name))

(define (run-clustalw-alignment imol chain-id target-sequence-pir-file)

  ;; write out the current sequence
  (let ((current-sequence-pir-file "current-sequence.pir")
        (aligned-sequence-pir-file "aligned-sequence.pir")
        (clustalw2-output-file-name "clustalw2-output-file.log"))

    (if (file-exists? aligned-sequence-pir-file)
        (delete-file aligned-sequence-pir-file))
    (if (file-exists? "aligned-sequence.dnd")
        (delete-file "aligned-sequence.dnd"))
    (if (file-exists? "current-sequence.dnd")
        (delete-file "current-sequence.dnd"))

    (write-current-sequence-as-pir imol chain-id current-sequence-pir-file)
    (goosh-command
     "clustalw2"
     '()
     (list
     "3"
     "1"
     target-sequence-pir-file
     "2"
     current-sequence-pir-file
     "9"
     "2"
     ""
     "4"
     ""
     aligned-sequence-pir-file
     ""
     "x"
     ""
     "x")
     clustalw2-output-file-name
     #t)

    (associate-pir-alignment-from-file imol chain-id aligned-sequence-pir-file)
    (apply-pir-alignment imol chain-id)
    (simple-fill-partial-residues imol)
    (resolve-clashing-sidechains-by-deletion imol)
    ))
