
(define *mapman-exe* (string-append (getenv "HOME") "/usf/lx_mapman"))

(define (f n) 
  (- (if (even? n)
	 (/ n 2) ; un-fortran indexing
	 (/ (- n 1) 2))
     1))


(define (generic-object-from-bones bones-file)
  (let ((atom-xyz-list '())
	(bone-list '())
	(conn-list '()))

    (define add-ordinate
      (let ((running-list '()))
	(lambda (x)
	  (cond
	   ((< (length running-list) 2)
	    (set! running-list (cons x running-list)))
	   (else 
	    (set! atom-xyz-list (cons (reverse (cons x running-list))
				      atom-xyz-list))
	    (set! running-list '()))))))

    (define add-bone
      (lambda (x)
	(set! bone-list (cons x bone-list))))

    (define add-connection-in-pairs
      (let ((running-list '()))
      (lambda (con)
	(cond 
	 ((null? running-list)
	  (set! running-list (list con)))
	 (else 
	  (set! conn-list (cons (reverse (cons con running-list))
				conn-list))
	  (set! running-list '()))))))

    (define add-connection
      (lambda (con)
	(cond 
	 ((null? conn-list)
	  (set! conn-list (list con)))
	 (else 
	  (set! conn-list (cons con conn-list))))))

    (define (my-indexer n)
      (let ((num (f n)))
	(if (integer? num)
	    num
	    (inexact->exact (+ (/ 1 2) num)))))
    

    ;; main body
    (if (string? bones-file)
	(if (file-exists? bones-file)

	    (call-with-input-file bones-file
	      (lambda (port)
		
		(let ((obj #t)
		      (reading-atoms-flag #f)
		      (reading-atom-bone-flag #f)
		      (reading-connectivity-flag #f)
		      (reading-atom-colour-flag #f)
		      (reading-residue-name-flag #f)
		      (reading-residue-type-flag #f)
		      (reading-residue-pointers-flag #f)
		      (continue-flag #t))
		  
		  (while continue-flag
		       
			 (set! obj (read port))
			 
			 (cond
			  ((eof-object? obj)
			   (begin
			     (format #t ":::: reached end-of-file~%")
			     (set! continue-flag #f)))

			  ;; atom colour
			  ((eq? obj 'SKL_ATOM_COLOUR)
			   (format #t "found atom-colour: ~s~%" obj)
			   (set! reading-connectivity-flag #f))
			   
			  ;; connectivity
			  ((eq? reading-connectivity-flag #t)
			   (add-connection obj))
			  ((number? reading-connectivity-flag)
			   (format #t "number reading-connectivity-flag: ~s~%" 
				   reading-connectivity-flag)
			   (if (< reading-connectivity-flag 3)
			       (set! reading-connectivity-flag (+ reading-connectivity-flag 1))
			       (set! reading-connectivity-flag #t)))
			  ((eq? obj 'SKL_CONNECTIVITY)
			   (format #t "found atom-connectivity: ~s~%" obj)
			   (set! reading-connectivity-flag 1)
			   (set! reading-atom-bone-flag #f))

			  ;; bones
			  ((eq? reading-atom-bone-flag #t)
			   (add-bone obj))
			  ((number? reading-atom-bone-flag)
			   (format #t "number reading-atom-bone-flag: ~s~%" 
				   reading-atom-bone-flag)
			   (if (< reading-atom-bone-flag 3)
			       (set! reading-atom-bone-flag (+ reading-atom-bone-flag 1))
			       (set! reading-atom-bone-flag #t)))
			  ((eq? obj 'SKL_ATOM_BONE)
			   (format #t "found atom-bone: ~s~%" obj)
			   (set! reading-atoms-flag #f)
			   (set! reading-atom-bone-flag 1))

			  ;; atoms 
			  ((eq? #t reading-atoms-flag)
			   (add-ordinate obj))
			  ((number? reading-atoms-flag)
			   (format #t "number reading-atoms-flag: ~s~%" reading-atoms-flag)
			   (if (< reading-atoms-flag 3)
			       (set! reading-atoms-flag (+ reading-atoms-flag 1))
			       (set! reading-atoms-flag #t)))
			  ((eq? obj 'SKL_ATOM_XYZ)
			   (format #t "found atom-xyz: ~s~%" obj)
			   (set! reading-atoms-flag 1))

			  ;; something else..
			   (else 'something))))))))

    
    (format #t "length atom-xyz-list: ~s~% length bone-list: ~s~% length conn-list ~s~%" 
	    (length atom-xyz-list) (length bone-list) (length conn-list))

    (let ((ordered-atoms (reverse atom-xyz-list))
	  (lines-obj (new-generic-object-number "Lines")))

;       (format #t "rcon: ~s~%" (reverse conn-list))
;       (format #t "ordered-atoms: ~s~%" ordered-atoms)

      (if #t
	  (let loop ((connections (reverse conn-list))
		     (current-position #f))
	    (cond
	     ((null? connections)
	      (set-display-generic-object lines-obj 1))
	     (else 
	      (let ((next-index (car connections)))
		(if (not current-position)
		    (begin
;		      (format #t "(from ~s) init move to : ~s ~%" (car connections) 
;			      (my-indexer (car connections)))
		      (loop (cdr connections) (my-indexer (car connections))))
		    (let ((xyz-index (my-indexer (car connections))))
		      (if (even? (car connections))
			  (begin
			    ;; a move to: 
;			    (format #t "(from ~s) move to : ~s ~%" (car connections)
;				    xyz-index)
			    (loop (cdr connections) xyz-index))
			  (begin
			    ;; a draw to:
;			    (format #t "(from ~s) line: ~s ~s ~%" (car connections)
;				    current-position xyz-index)
			    (apply to-generic-object-add-line lines-obj "green" 3
				   (append (list-ref ordered-atoms current-position)
					   (list-ref ordered-atoms xyz-index)))
			    (loop (cdr connections) xyz-index))))))))))

; debugging: output the points too.					   
;      (let ((points-obj (new-generic-object-number "Points")))
;	(let loop ((xyz-list ordered-atoms)
;		   (point-number 0))
;	  (cond
;	   ((null? xyz-list)
;	    (set-display-generic-object points-obj 1))
;	   (else 
;	    (apply to-generic-object-add-point points-obj "green" 10 (car xyz-list))
;	    ; (apply place-text (number->string point-number) (append (car xyz-list) (list 2)))
;	    (loop (cdr xyz-list) (+ point-number 1))))))

      )))



;;;
(define (bones-it map-file-name)
  
  (let* ((bones-file "my.bones")
	 (data-lines (list (string-append "read m1 " map-file-name " ccp4")
			   "bo sk m1 0.5 0.15 1"
			   "bones connect"
			   bones-file
			   "skl"
			   "5"
			   "quit")))
				   
    (goosh-command *mapman-exe* '() data-lines "coot-mapman.log" #t)
    ;; now read in bones-file
    (generic-object-from-bones bones-file)))

								    
;;; Add a menu item:
(let ((menu (coot-menubar-menu "Mapman")))
  (add-simple-coot-menu-menuitem 
   menu "Mapman Bones..."
   (lambda ()
     (map-molecule-chooser-gui "Map to Bonesify:"
			       (lambda (imol)
				 (format #t "bonesing ~s~%" imol)
				 (export-map imol "tmp.map")
				 (bones-it "tmp.map"))))))
				 


