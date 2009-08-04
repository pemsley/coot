
(use-modules (sxml simple))


(define (pisa-xml imol file-name)
  (if (file-exists? file-name)
      (begin
	(format #t "opened ~s~%" file-name)
	(call-with-input-file file-name
	  (lambda (port)
	    (let ((sm (xml->sxml port)))
	      ;; (format #t "sm: ~s~%" sm)
	      (parse-pisa imol sm)))))))


;; pisa_results
;;    name
;;    status
;;    total_asm
;;    asm_set
;;       ser_no
;;       assembly
;;          id
;;          size
;;          mmsize
;;          diss_energy
;;          asa
;;          bas
;;          entropy
;;          diss_area
;;          int_energy
;;          n_uc
;;          n_diss
;;          symNumber
;;          molecule
;;             chain_id
;;             rxx
;;             rxy
;;             rxz
;;             tx
;;             ryx
;;             ryy
;;             ryz
;;             ty
;;             rzx
;;             rzy
;;             rzz
;;             tz
;;             rxx-f
;;             rxy-f
;;             rxz-f
;;             tx-f
;;             ryx-f
;;             ryy-f
;;             ryz-f
;;             ty-f
;;             rzx-f
;;             rzy-f
;;             rzz-f
;;             tz-f


  
(define parse-pisa
  (lambda (imol entity)

    ;; was it a chain? (or residue selection disguised as
    ;; a chain?)
    ;; 
    (define (make-atom-selection-string chain-id-raw)

      ;; first try to split the chain-id-raw on a "]".  If there was
      ;; no "]" then we had a simple chain-id.  If there was a "]"
      ;; then we have something like "[CL]D:32", from which we need to
      ;; extract a residue number and a chain id.  Split on ":" and
      ;; construct the left and right hand sides of s.  Then use those
      ;; together to make an mmdb selection string.
      ;; 
      (let ((match-info (string-match "]" chain-id-raw)))
			 
	(if match-info
	    (begin
	      ;; (format #t "found a residue selection chain_id ~s ~s ~%" chain-id-raw match-info)
	      (let* ((l (string-length chain-id-raw))
		     (s (substring chain-id-raw (match:end match-info) l)) ; e.g. "D:32"
		     (ls (string-length s))
		     (s-match-info (string-match ":" s)))
		(if s-match-info
		    (let ((residue-string (substring s (match:end s-match-info) ls))
			  (chain-string (substring s 0 (match:start s-match-info))))
		      ;; (print-var residue-string)
		      ;; (print-var chain-string)
		      (string-append "//" chain-string "/" residue-string))
		    s)))
	    (begin
	      ;; (format #t "found a simple chain_id ~s~%" chain-id-raw)
	      (string-append "//" chain-id-raw)))))

    
    ;; Return a molecule number or #f
    ;; 
    (define (handle-molecule molecule)
      (let ((symbols (list 'rxx 'rxy 'rxz 'ryx 'ryy 'ryz 'rzx 'rzy 'rzz 'tx 'ty 'tz))) ;; orthogonal 
                                                                                       ;; matrix elements.

	(let ((ass-symbols '())
	      (atom-selection-string "//") ;; default, everything
	      (symm-name-part "")) ;; the atom selection info without
				   ;; "/" because that would chop the
				   ;; name in the display manager.
	  (for-each
	   (lambda (mol-ele)
	     ;; (format #t "mol-ele: ~s~%" mol-ele)
	     (if (list? mol-ele)
		 (begin
		   
		   (if (eq? (car mol-ele) 'chain_id)
		       (let* ((chain-id-raw (cadr mol-ele)))
			 (set! atom-selection-string 
			       (make-atom-selection-string chain-id-raw))
			 (if (string=? chain-id-raw atom-selection-string)
			     (set! symm-name-part (string-append "chain " chain-id-raw))
			     (set! symm-name-part chain-id-raw))))
		 
		   ;; was it one of the rotation or translation symbols?
		   (for-each 
		    (lambda (symbol)
		      (if (eq? (car mol-ele) symbol)
			  (let ((n (string->number (cadr mol-ele))))
			    ;; (format #t "~s: ~s~%" symbol n)
			    (set! ass-symbols 
				  (cons (list symbol n)
					ass-symbols)))))
		   
		    symbols))))
	   molecule)
	
	  (if (= (length ass-symbols) 12)
	      (let ((mat 
		     (map (lambda (sym)
			    (cadr (assoc sym ass-symbols)))
			  symbols)))
		;; (format #t "mat:::: ~s~%" mat)
		(format "currently ~s molecules~%" (graphics-n-molecules))
		(let* ((new-molecule-name (string-append "Symmetry copy of "
							 (number->string imol)
							 " using " symm-name-part))
		       (new-mol-no (apply new-molecule-by-symmetry-with-atom-selection
					 imol
					 new-molecule-name
					 atom-selection-string
					 (append mat (list 0 0 0)))))
		  ;; (format #t "created molecule number ~s~%" new-mol-no)
		  new-mol-no))
	      #f ;; ass-symbols were not all set
	      ))))

    ;; Return the model number of the new assembly molecule
    ;; 
    (define (create-assembly-set-molecule assembly-set-molecule-numbers)
      (if (null? assembly-set-molecule-numbers)
	  #f
	  (let ((first-copy (copy-molecule (car assembly-set-molecule-numbers))))
	    (if (not (valid-model-molecule? first-copy))
		#f
		(let ((rest (cdr assembly-set-molecule-numbers)))
		  (merge-molecules rest first-copy)
		  (set-molecule-name first-copy "An Assembly Set")
		  first-copy)))))

    ;; Return a list of model numbers
    ;; 
    (define (create-assembly-molecule assembly-molecule-numbers)
      assembly-molecule-numbers)

    ;;
    (define (handle-assembly assembly)
      
      (let ((assembly-molecule-numbers '())
	    (assembly-size #f)
	    (assembly-mmsize #f) 
	    (assembly-id #f)
	    (assembly-symm-number #f)
	    (assembly-n-uc #f)
	    (assembly-n-diss #f)
	    (assembly-asa #f)
	    (assembly-bsa #f)
	    (assembly-diss-energy #f)
	    (assembly-diss-area #f)
	    (assembly-entropy #f)
	    (assembly-int-energy #f))
	
	(for-each
	 (lambda (ele)
	   ;; (format #t "ele: ~s~%" ele)
	   (if (not (list? ele))
	       (format #t "=== ele: ~s not list~%" ele)
	       (cond 
		((eq? (car ele) 'molecule)
		 (let ((mol-no (handle-molecule ele)))
		   (set! assembly-molecule-numbers 
			 (append assembly-molecule-numbers (list mol-no)))))
		((eq? (car ele) 'symNumber)
		 (set! assembly-symm-number (cadr ele)))
		((eq? (car ele) 'size)
		 (set! assembly-size (string->number (cadr ele))))
		((eq? (car ele) 'mmsize)
		 (set! assembly-mmsize (string->number (cadr ele))))
		((eq? (car ele) 'id)
		 (set! assembly-id (cadr ele)))
		((eq? (car ele) 'asa)
		 (set! assembly-asa (string->number (cadr ele))))
		((eq? (car ele) 'bsa)
		 (set! assembly-bsa (string->number (cadr ele))))
		((eq? (car ele) 'diss_energy)
		 (set! assembly-diss-energy (string->number (cadr ele))))
		((eq? (car ele) 'diss_area)
		 (set! assembly-diss-area (string->number (cadr ele))))
		((eq? (car ele) 'entropy)
		 (set! assembly-entropy (string->number (cadr ele))))
		((eq? (car ele) 'int_energy)
		 (set! assembly-int-energy (string->number (cadr ele))))
		((eq? (car ele) 'n_diss)
		 (set! assembly-n-diss (string->number (cadr ele))))
		((eq? (car ele) 'n_uc)
		 (set! assembly-n-diss (string->number (cadr ele))))
		(else 
		 (format #t "unhandled ele: ~s~%" ele)))))
	 
	 assembly)

	;; use an assembly record here?
	(map (lambda (mol-no)
	       (if (valid-model-molecule? mol-no)
		   (set-mol-displayed mol-no 0)))
	     assembly-molecule-numbers)
	assembly-molecule-numbers))

    
    ;; handle assembly-set
    (define (handle-assembly-set assembly-set)
      
      (let ((assembly-set-molecules '()))
	(for-each (lambda (assembly-set-part)
		    (if (list? assembly-set-part)
			(if (eq? 'assembly (car assembly-set-part))
			    (let ((assembly (handle-assembly assembly-set-part)))
			      (set! assembly-set-molecules
				    (append assembly-set-molecules (list assembly)))))))
		  assembly-set)

	(format #t "do something with this assembly-set molecule list: ~s~%" 
		assembly-set-molecules)

	(create-assembly-set-molecule (apply append assembly-set-molecules))))

      

    


    ;; main line
    (let loop ((entity entity))
      (cond
       ((list? entity) (if (not (null? entity))
			   (begin
			     (if (eq? 'asm_set (car entity))
				 (begin
				   (handle-assembly-set entity))
				 (map loop entity)))))
; 			     (if (eq? 'assembly (car entity))
; 				 (begin 
; 				   (handle-assembly entity))
; 				 (map loop entity)))))
       (else 
	#f)))))
     


;;(pisa-xml "pisa.xml")

(define (t)
  (pisa-xml (read-pdb "coot-download/3FCS.pdb") "pisa.xml"))



(let ((menu (coot-menubar-menu "PISA")))

    (add-simple-coot-menu-menuitem
     menu "PISA..."
     (lambda ()
       (pisa-xml (read-pdb "coot-download/pdb3lz2.ent") "pisa.xml")))

    (add-simple-coot-menu-menuitem
     menu "entry code..."
     (lambda ()
       (generic-chooser-and-entry "PISA xml for" " XML file-name:" "" 
			     (lambda (imol text)
			       (pisa-xml imol text))))))


