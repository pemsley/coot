
(use-modules (sxml simple)
	     (ice-9 receive))


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

    ;; 
    (define (print-assembly record port)
      (let ((rec-type (record-type-descriptor record)))
	(format port "assembly id: ~s~%"          ((record-accessor rec-type 'id)          record))
	(format port "assembly size: ~s~%"        ((record-accessor rec-type 'size)        record))
	(format port "assembly symm-number: ~s~%" ((record-accessor rec-type 'symm-number) record))
	(format port "assembly asa: ~s~%"         ((record-accessor rec-type 'asa)         record))
	(format port "assembly bsa: ~s~%"         ((record-accessor rec-type 'bsa)         record))
	(format port "assembly diss-energy: ~s~%" ((record-accessor rec-type 'diss-energy) record))
	(format port "assembly entropy: ~s~%"     ((record-accessor rec-type 'entropy)     record))
	(format port "assembly diss-area: ~s~%"   ((record-accessor rec-type 'diss-area)   record))
	(format port "assembly int-energy: ~s~%"  ((record-accessor rec-type 'int-energy)  record))
	(format port "assembly components: ~s~%"  ((record-accessor rec-type 'molecule-numbers) record))))
	

    ;; Return a list of model numbers
    ;; 
    (define (create-assembly-molecule assembly-molecule-numbers)
      assembly-molecule-numbers)


    ;; If the molecule does not "complexate in solution" then we get
    ;; "0" for total_asm element.
    ;;
    (define (check-there-are-assemblies? entity)
      (if (string=? (list-ref entity 1) "0")
	  (info-dialog "This molecule does not form a complex in solution\n(most probably)")))



    ;; Return an assembly record
    ;; 
    (define (handle-assembly assembly make-assembly-record)

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
	       (format #t "ele not list: ~s~%" ele) ; caution when deleting
	       (cond 
		((eq? (car ele) 'molecule)
		 (let ((mol-no (handle-molecule ele)))
		   (set! assembly-molecule-numbers 
			 (append assembly-molecule-numbers (list mol-no)))))
		((eq? (car ele) 'symNumber)
		 (set! assembly-symm-number (string->number (cadr ele))))
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

	;; assembly-molecule-numbers ;; return this
	(make-assembly-record assembly-id assembly-size assembly-symm-number assembly-asa assembly-bsa
			      assembly-diss-energy assembly-entropy assembly-diss-area 
			      assembly-int-energy assembly-molecule-numbers)))

    
    ;; handle assembly-set.
    ;; 
    ;; return values [#f or the molecule number list of the
    ;; assembly-set] and the assembly-record-set
    ;; 
    (define (handle-assembly-set assembly-set make-assembly-record record-type)
      
      (let ((assembly-set-molecules '())
	    (assembly-records '()))
	(for-each (lambda (assembly-set-part)
		    (if (list? assembly-set-part)
			(if (eq? 'assembly (car assembly-set-part))
			    (let* ((assembly-record (handle-assembly assembly-set-part make-assembly-record))
				   (assembly-molecules-numbers (lambda (assembly-record)
								 ((record-accessor record-type 'molecule-numbers)
								  assembly-record)))
				   (assembly-molecules (assembly-molecules-numbers assembly-record)))
			      (set! assembly-records 
				    (append assembly-records (list assembly-record)))
			      (set! assembly-set-molecules
				    (append assembly-set-molecules (list assembly-molecules)))))))
		  assembly-set)

	(let ((new-mol 
	       (create-assembly-set-molecule (apply append assembly-set-molecules))))
	  (values new-mol assembly-records))))

      

    ;; main line
    ;; 
    ;; first set up the assembly record:
    ;; 
    (let* ((rt (make-record-type "assembly" '(id size symm-number asa bsa diss-energy entropy diss-area
						 int-energy molecule-numbers) print-assembly))
	   (make-assembly-record (record-constructor rt)))
      ;; 
      (let ((first-assembly-set-is-displayed-already? #f)
	    (top-assembly-set #f))
	(let loop ((entity entity))
	  (cond
	   ((list? entity) (if (not (null? entity))
			       (begin
				 (if (eq? 'total_asm (car entity))
				     (check-there-are-assemblies? entity))
				 (if (eq? 'asm_set (car entity))
				     (receive (molecule-number assembly-record-set)
					      (handle-assembly-set entity make-assembly-record rt)
					      (if (not top-assembly-set)
						  (set! top-assembly-set assembly-record-set))
					      (if first-assembly-set-is-displayed-already?
						  (set-mol-displayed molecule-number 0)
						  (set! first-assembly-set-is-displayed-already? #t)))
				     (map loop entity)))))
	   (else 
	    #f)))

	(format #t "top assembly-set:~% ~s~%" top-assembly-set)
	;; it is ugly if we get a dialog saying that the
	;; top-assembly-set is #f.  So inhibit the dialog in that case
	;; (we already check above and give a dialog if total_asm is
	;; 0).
	(if top-assembly-set
	    (let ((s (format #f "top assembly-set:~% ~s~%" top-assembly-set)))
	      (info-dialog s)))))))
     


;;(pisa-xml "pisa.xml")


(define (pisa-assemblies imol)

  ;; 
  (define (make-stubbed-name imol)
    (strip-extension (basename (molecule-name imol))))

  ;; 
  (define (make-pisa-config pisa-coot-dir config-file-name)
    (call-with-output-file config-file-name
      (lambda (port)
	(let ((s (getenv "CCP4")))
	  (if (file-exists? s)
	      (for-each (lambda (item-pair)
			  (display (car item-pair) port)
			  (newline port)
			  (display (cdr item-pair) port)
			  (newline port))
			(list (cons "DATA_ROOT"   (append-dir-file s "share/pisa"))
			      (cons "PISTORE_DIR" (append-dir-file s "share/pisa/pisastore"))
			      (cons "PISA_WORK_ROOT"  pisa-coot-dir)
			      (cons "SBASE_DIR"  (append-dir-file s "share/sbase"))
			      ;; that we need these next 3 is ridiculous
			      (cons "MOLREF_DIR" (append-dir-file s "share/pisa/molref"))
			      (cons "RASMOL_COM" (append-dir-file s "bin/rasmol"))
			      (cons "CCP4MG_COM" (append-dir-file s "bin/ccp4mg"))
			      (cons "SESSION_PREFIX" "pisa_"))))))))

	
  
  ;; main line
  (if (valid-model-molecule? imol) 
      (let* ((pisa-coot-dir "coot-pisa")
	     (stubbed-name (make-stubbed-name imol))
	     (pdb-file-name (append-dir-file pisa-coot-dir (string-append stubbed-name ".pdb")))
	     (pisa-config (append-dir-file pisa-coot-dir (string-append stubbed-name "-pisa.cfg")))
	     (pisa-xml-file-name (append-dir-file pisa-coot-dir (string-append stubbed-name ".xml"))))
	
	(make-directory-maybe pisa-coot-dir)
	(make-pisa-config pisa-coot-dir pisa-config)
	(write-pdb-file imol pdb-file-name)
	(if (file-exists? pdb-file-name)

	    (begin
	      (format #t "pisa args ~s~%" (list "pisa" "-analyse" pdb-file-name pisa-config))
	      (let ((status (goosh-command "pisa" (list "pisa" "-analyse" pdb-file-name pisa-config) 
					   '() "pisa.log" #f)))
		(if (not (number? status))
		    (info-dialog "Ooops PISA failed to deliver the goods!\n\n(Go for a curry instead?)")
		    (if (= status 0) ;; good
			(let ((status-2 (goosh-command "pisa" (list "pisa" "-xml" pisa-config) '() 
						       pisa-xml-file-name #f)))
			  (if (= status 0)
			      (pisa-xml imol pisa-xml-file-name)))))))))))
  


