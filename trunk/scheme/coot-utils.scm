
;;;; Copyright 2000 by Paul Emsley
;;;; Copyright 2004, 2005, 2006, 2007 by The University of York
;;;; Copyright 2008 by The University of York
;;;; Copyright 2008, 2009 by The University of Oxford

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
	     (ice-9 rdelim))
(use-modules (os process)) 

;; 3D annotations - a bit of a hack currently
(define *annotations* '())

;; scm aliases.  Don't forget to add them to snarf-coot-docs too.
;; 
(define drag-intermediate-atom drag-intermediate-atom-scm)
(define mark-atom-as-fixed mark-atom-as-fixed-scm)
(define ncs-chain-ids ncs-chain-ids-scm)
(define ncs-chain-differences ncs-chain-differences-scm)
(define refmac-parameters refmac-parameters-scm)
(define water-chain water-chain-scm)
(define water-chain-from-shelx-ins water-chain-from-shelx-ins-scm)
(define generic-object-name generic-object-name-scm)
(define additional-representation-info additional-representation-info-scm)
(define missing-atom-info missing-atom-info-scm)
(define key-sym-code key-sym-code-scm)
(define ncs-ghosts ncs-ghosts-scm)
(define inverse-rtop inverse-rtop-scm)
(define coot-has-python? coot-has-python-p)
(define map-sigma map-sigma-scm)
(define pucker-info pucker-info-scm)
(define map-parameters map-parameters-scm)
(define cell cell-scm)
(define map-cell cell)
(define ccp4i-projects ccp4i-projects-scm)
(define map-sigma map-sigma-scm)
(define get-rotamer-name get-rotamer-name-scm)
(define test-internal test-internal-scm)
(define atom-info-string atom-info-string-scm)
(define get-refmac-sad-atom-info get-refmac-sad-atom-info-scm)
;; fix typo
(define set-find-hydrogen-torsion set-find-hydrogen-torsions)
(define residues-near-position residues-near-position-scm)
(define non-standard-residue-names non-standard-residue-names-scm)
(define refine-residues refine-residues-scm)
(define map-peaks map-peaks-scm)
(define map-peaks-near-point map-peaks-near-point-scm)
(define add-dipole add-dipole-scm)
(define add-dipole-for-residues add-dipole-for-residues-scm)
(define get-torsion get-torsion-scm)
(define test-internal-single test-internal-single-scm)
(define user-defined-click user-defined-click-scm)
(define add-lsq-atom-pair add-lsq-atom-pair-scm)
(define coot-sys-build-type coot-sys-build-type-scm)
(define add-alt-conf add-alt-conf-scm)
(define origin-pre-shift origin-pre-shift-scm)
(define alignment-mismatches alignment-mismatches-scm)
(define rigid-body-refine-by-residue-ranges rigid-body-refine-by-residue-ranges-scm)
(define average-map average-map-scm)
(define symmetry-operators symmetry-operators-scm)
(define symmetry-operators->xHM symmetry-operators-to-xHM-scm)
(define user-mods user-mods-scm)

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
       ,@funcs
       (accept-regularizement)
       (if (= replace-state 0)
	   (set-refinement-immediate-replacement 0)))))

;; 
(defmacro using-active-atom funcs

  `(begin
     (let ((active-atom (active-residue)))
       (if (not active-residue)
	   (add-status-bar-text "No residue found")
	   (let ((aa-imol      (list-ref active-atom 0))
		 (aa-chain-id  (list-ref active-atom 1))
		 (aa-res-no    (list-ref active-atom 2))
		 (aa-ins-code  (list-ref active-atom 3))
		 (aa-atom-name (list-ref active-atom 4))
		 (aa-alt-conf  (list-ref active-atom 5)))
	     
		 ,@funcs
		 
		 )))))

; e.g.:
; (with-auto-accept 
;   (format #t "tum tee tumm...~%")
;   (format #t "tra la la...~%"))

(defmacro print-var (var)
  (let ((s (format #f "~s" var)))
    `(format #t "DEBUG:: ~a is ~s~%" ,s ,var)))

;; schemify function
;;
(define (molecule-has-hydrogens? imol)
  (= (molecule-has-hydrogens-raw imol) 1))


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
    (format #t "virtual trackball type ~s not understood~%"))))
    
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



;;; Make list of integers, @var{a} to @var{b}: eg (number-list 2 5) ->
;;; (2 3 4 5)
;; 
(define (number-list a b)

  (let loop ((count a)
	     (acc '()))
    (cond
     ((> count b) (reverse acc))
     (else
      (loop (+ count 1) (cons count acc))))))



;;; ls must be a list of strings, atom must be a string.
;;; 
;;; return either #t or #f.
;;;
(define (string-member? atom ls)
  
  (cond 
   ((null? ls) #f)
   ((string=? atom (car ls)) #t)
   (else
    (string-member? atom (cdr ls)))))

;;; range: works like the eponymous python function
;;; e.g. (range 3) -> '(0, 1, 2)
;;; e.g. (range 1 3) -> '(1, 2)
;;;
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

;; code from thi <ttn at mingle.glug.org>
;; 
;; run command and put the output into a string and return
;; it. (c.f. @code{run-command/strings})
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

;; Return #t or #f:
(define (command-in-path? cmd)

  ;; test for command (see goosh-command-with-file-input description)
  ;; 
  (if (string? cmd) 
      (let ((have-command? (run "which" cmd)))
	
	(= have-command? 0)) ; return #t or #f
      #f))

      
;; Where cmd is e.g. "refmac" 
;;       args is (list "HKLIN" "thing.mtz")
;;       log-file-name is "refmac.log"      
;;       data-list is (list "HEAD" "END")
;; 
;; Return the exist status e.g. 0 or 1.
;; 
(define (goosh-command cmd args data-list log-file-name screen-output-also?)
    
  (if (not (command-in-path? cmd))

      (format #t "command ~s not found~%" cmd)

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
		      (format #t "exit status: ~s~%" status)
		      status)) ; return status 
		  
		  (begin
		    (if (eq? screen-output-also? #t)
			(format #t ":~a~%" obj))
		    (format log-file-port "~a~%" obj)
		    (f (read-line output-port))))))))))

; example usage:
;(goosh-command "mtzdump" (list "HKLIN" "a.mtz") (list "HEAD" "END") "test.log" #t)

;; 
;;
(define (goosh-command-with-file-input cmd args input-file log-file-name)

  ;; we had a problem here (may 2004).  If command was not in path,
  ;; this function either hung or caused coot to immediately stop.
  ;; So now we test to see if cmd exists in the path first, by
  ;; running "which" on it.
  ;; 
  (if (string? cmd) 
      (let ((have-command? (run "which" cmd)))

	(if (= have-command? 0) ; we *do* have the command

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
    
  (if (not (command-in-path? cmd))

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


;; Return an atom info or #f (if atom not found).
;;
(define (get-atom imol chain-id resno ins-code atom-name . alt-conf)
  

  (define (get-atom-from-res atom-name residue-atoms alt-conf)
    (let loop ((residue-atoms residue-atoms))
      (cond 
       ((null? residue-atoms) #f)
       ((and (string=? atom-name (car (car (car residue-atoms))))
	     (string=? alt-conf (car (cdr (car (car residue-atoms))))))
	(car residue-atoms))
       (else 
	(loop (cdr residue-atoms))))))

  ;; main line
  ;; 
  ;; (format #t "DEBUG::  get-atom is passed imol ~s: chain-id: ~s  resno: ~s  atom-name: ~s alt-conf: ~s ~%" imol chain-id resno atom-name alt-conf)

  (let ((res-info (residue-info imol chain-id resno ""))
	(alt-conf-internal (if (null? alt-conf)
			       ""
			       (car alt-conf))))

    (if (not res-info)
	#f
	(get-atom-from-res atom-name res-info alt-conf-internal))))
	


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
	 ((char=? (car chars) #\space) (f (cdr chars)
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

;;;  in a laughable attempt to minimise system dependence.
;;;
(define (append-dir-file dir-name file-name)

  (if (> (string-length dir-name) 0)
      (string-append (directory-as-file-name dir-name) "/" file-name)
      file-name))

;;; similarly attempting to minimise system dependence.
(define (append-dir-dir dir-name sub-dir-name)

  (string-append (directory-as-file-name dir-name) "/" sub-dir-name))

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


;; return #f if dir-name is a file or we can't do the mkdir
;; 
(define (coot-mkdir dir-name)

  (if (file-exists? dir-name)
      (if (is-directory? dir-name)
	  dir-name ; often we will be returning this
	  #f)
      (mkdir dir-name)))
	    

;; The following functions from PLEAC (guile version thereof of course).
;; 
;; or define a utility function for this
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
(define (glob pat dir)
  (let ((rx (make-regexp (glob->regexp pat))))
    (filter (lambda (x) (regexp-exec rx x)) (directory-files dir))))


;;; return the view matrix (useful for molscript, perhaps).
(define view-matrix
  (lambda ()

    (map (lambda (row-number)
	   (map (lambda (column-number)
		  (get-view-matrix-element row-number column-number))
		(list 0 1 2)))
	 (list 0 1 2))))

;;; return the view quaternion
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


;;; Miguel's molecular orientation axes
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
;; return #f on error
(define (translation axis length)

  (cond 
   ((not (symbol? axis)) (format #t "incomprehensible axis argument: ~s~%" axis))
   ((not (number? length)) (format #t "incomprehensible length argument: ~s~%" length))
   ((eq? axis 'x) (list length 0 0))
   ((eq? axis 'y) (list 0 length 0))
   ((eq? axis 'z) (list 0 0 length))
   (else 
    (format #t "symbol axis: ~s incomprehensible~%")
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
    (format #t "symbol axis: ~s incomprehensible~%")
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

    (cond 
     ((= (length args) 7)
      (tf (list-ref args 0)
	  (list-ref args 1)
	  (list-ref args 2)
	  (list-ref args 3)
	  (list-ref args 4)
	  (list-ref args 5)
	  (list-ref args 6)))
     ((= (length args) 4) ; no matrix specified
      (tf (list-ref args 0)
	  (identity-matrix)
	  (list-ref args 1)
	  (list-ref args 2)
	  (list-ref args 3)
	  (symmetry-operators->xHM (symmetry-operators imol))
	  (cell imol)))
     ((= (length args) 3) ; no matrix or about point specified
      (tf (list-ref args 0)
	  (identity-matrix)
	  (list-ref args 1)
	  (rotation-centre)
	  (list-ref args 2)
	  (symmetry-operators->xHM (symmetry-operators imol))
	  (cell imol)))
     (else
      (format #t "arguments to transform-map incomprehensible: args: ~s~%" args)
      #f))))
	 
     
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
    
    (if (not (and space-group cell-params))
	(let ((message (format #f "Bad cell or symmetry for molecule ~s~%"
			       cell space-group imol-ref))))

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
;; Return #t or #f. 
(define (valid-model-molecule? imol)

  (if (not (number? imol))
      #f
      (= (is-valid-model-molecule imol) 1)))

;; schemey interface to eponymous scripting interface function.
;; Return #t or #f. 
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
  (if (not (is-valid-map-molecule? imol-map))
      #f
      (= (map-is-difference-map imol-map) 1)))
      

;; Does residue resno with insertion code ins-code of chain chain-id
;; and in molecule number imol exist?  
;; 
;; Return #t or #f.
;; 
(define (residue-exists? imol chain-id resno ins-code)

  (= 1 (does-residue-exist-p imol chain-id resno ins-code)))


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
	(assign-pir-sequence imol chain-id seq-text)
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

(define (add-key-binding name key thunk)
  (if (number? key)
      (set! *key-bindings* (cons (list key key name thunk) *key-bindings*)))
  (if (string? key)
      (let ((code (key-sym-code key)))
	(if (not (= code -1))
	    (set! *key-bindings* (cons (list code key name thunk) *key-bindings*))))))


;; general key press hook
;; 
(define (graphics-general-key-press-hook key)
  ;; (format #t "Key ~s was pressed~%" key)
  (let ((field (assoc key *key-bindings*)))
    ;; (format #t "Field: ~s~%" field)
    (if (not field)
	(format #t "Key ~s not found in bindings~%" key)
	(begin
	  ((car (cdr (cdr (cdr field)))))))))


;;; Function requested by Mark White.
;;; 
;;; read XtalView (and maybe other) .vu files and convert them into generic 
;;; objects.  
;;;
;;; Pass the filename and an object name e.g.
;;; (read-vu-file "axes.vu" "axes")
;;;
;;; Returns: nothing interesting.
;;;
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
  
  (reverse 
   (let f ((chain-list (chain-ids imol))
	   (seq-residue-list #f)
	   (alt-conf-residue-list '()))

     ;; "double-shuffle" for starting a new chain...
     (cond
      ((null? chain-list) alt-conf-residue-list)
      ((null? seq-residue-list) (f (cdr chain-list)
				   #f
				   alt-conf-residue-list))
      ((eq? #f seq-residue-list) (f chain-list
				    (number-list 0 (- (chain-n-residues 
						       (car chain-list) imol) 1))
				    alt-conf-residue-list))
      (else 
       (let* ((chain-id (car chain-list))
	      (serial-number (car seq-residue-list))
	      (res-no (seqnum-from-serial-number imol chain-id (car seq-residue-list)))
	      (ins-code (insertion-code-from-serial-number imol chain-id serial-number)))

	 (let ((r (residue-test-func chain-id
				     res-no
				     ins-code
				     serial-number)))
	   (cond
	    ((eq? r #f) (f chain-list (cdr seq-residue-list) alt-conf-residue-list))
	    (else
	     (f chain-list
		(cdr seq-residue-list)
		(cons (list r chain-id res-no ins-code) 
		      alt-conf-residue-list)))))))))))

;; Return residue specs for all residues in imol (each spec is preceeded by #t)
;; 
(define (all-residues imol)
  (residues-matching-criteria imol (lambda (chain-id resno ins-code serial) #t)))

  
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
(define (mutate-by-overlap imol chain-id resno tlc)

  (define (mutate-it)
    (let ((imol-ligand (get-monomer tlc)))
      (if (not (valid-model-molecule? imol-ligand))
	  (let ((s (string-append " Oops.  Failed to get monomer " tlc)))
	    (add-status-bar-text s))
	  (begin
	    (delete-residue-hydrogens imol-ligand "A" 1 "" "")
	    (delete-atom imol-ligand "A" 1 "" " OXT" "")
	    (overlap-ligands imol-ligand imol chain-id resno)
	    (delete-residue imol chain-id resno "")
	    (let* ((new-chain-id-info (merge-molecules (list imol-ligand) imol))
		   (nov (format #t "new-chain-id-info: ~s~%" new-chain-id-info)))
	      (let ((merge-status (car new-chain-id-info)))
		(if (= merge-status 1)
		    (let ((new-chain-id (car (car (cdr new-chain-id-info)))))
		      (change-residue-number imol new-chain-id 1 "" resno "")
		      (change-chain-id imol new-chain-id chain-id 1 resno resno)

		      (let ((replacement-state (refinement-immediate-replacement-state))
			    (imol-map (imol-refinement-map)))
			(set-refinement-immediate-replacement 1)
			(if (= imol-map -1)
			    (regularize-zone imol chain-id resno resno "")
			    (let ((spin-atoms (list " P  " " O1P" " O2P" " O3P"))
				  (dir-atoms (cond 
					      ((string=? tlc "PTR") (list " CZ " " OH "))
					      ((string=? tlc "SEP") (list " CB " " OG "))
					      ((string=? tlc "TPO") (list " CB " " OG1"))
					      (else 
					       #f))))
			      (if dir-atoms
				  (spin-search imol-map imol chain-id resno "" 
					       dir-atoms spin-atoms))
			      (refine-zone imol chain-id resno resno "")))
			(accept-regularizement)
			(set-refinement-immediate-replacement replacement-state))
		      
		      (set-mol-displayed imol-ligand 0)
		      (set-mol-active imol-ligand 0)))))))))

  ;; First, if there are multiple maps, force the user to choose one,
  ;; rather than continuing.
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
				     "occ" (shelx-molecule? imol) 11.0 1.0))))))))
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


(define (print-molecule-names)

  (map (lambda (molecule-number)
	 (format #t "    ~s    ~s~%" molecule-number (molecule-name molecule-number)))
       (molecule-number-list)))

  
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


(define BALL_AND_STICK 2)

;; hilight-colour is specified in degrees (round the colour wheel -
;; starting at yellow (e.g. 230 is purple))
;; 
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
						   bb-type 6
						   draw-hydrogens-flag)
	  
	  (map (lambda (spec)
		 (additional-representation-by-attributes imol
							  (car spec)
							  (cadr spec)
							  (cadr spec)
							  (caddr spec)
							  BALL_AND_STICK
							  bb-type 6
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


;; ---------- annotations ---------------------

(define (add-annotation-here text)
  (let ((rc (rotation-centre))
	(ann (cons text (rotation-centre))))

    (set! *annotations* (cons ann *annotations*))
    (apply place-text text (append (rotation-centre) (list 0)))
    (graphics-draw)))

(define (save-annotations file-name)
  (call-with-output-file file-name
    (lambda (port)
      (format port "~s~%" *annotations*))))

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

	
; to determine if we have pygtk
;
(define (coot-has-pygtk?)
  (if (coot-has-python?)
      (run-python-command "coot_has_pygtk()")
      #f))


