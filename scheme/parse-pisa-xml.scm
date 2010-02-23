
(use-modules (sxml simple)
	     (ice-9 receive))



;; (define *pisa-command* "/home/paule/ccp4/ccp4-6.1.2/bin/pisa")
(define *pisa-command* "pisa")	


(define (pisa-assemblies imol)

  ;; 
  ;; main line
  (receive (pdb-file-name pisa-config pisa-xml-file-name)
	   (prep-for-pisa 'assemblies imol)

	   (begin
	     (format #t "DEBUG:: ========= pisa with args ~s~%" 
		     (list *pisa-command* "-analyse" pdb-file-name pisa-config))
	     (let ((status (goosh-command "pisa" (list "pisa" "-analyse" pdb-file-name pisa-config) 
					  '() "pisa.log" #f)))
	       (if (not (number? status))
		   (info-dialog "Ooops PISA failed to deliver the goods!\n\n(Go for a curry instead?)")
		   (if (= status 0) ;; good
		       (begin
			 (let ((status-2 (goosh-command *pisa-command* (list "pisa" "-xml" 
								     "assemblies" pisa-config) '() 
							pisa-xml-file-name #f)))
			   (if (= status 0)
			       (pisa-assemblies-xml imol pisa-xml-file-name))))))))))


;; called by pisa-assemblies, which is the entry point from the main
;; program gui.
;;
;;  it calls parse-pisa-assemblies which does the work
;; 
(define (pisa-assemblies-xml imol file-name)
  (if (file-exists? file-name)
      (begin
	(format #t "opened ~s~%" file-name)
	(call-with-input-file file-name
	  (lambda (port)
	    (let ((sm (xml->sxml port)))
	      (parse-pisa-assemblies imol sm)))))))



;; Exported to the main level.  A molecule is common to assemblies and interfaces.
;; 
;; If pisa-result-type is 'assemblies, return a molecule number or #f
;; 
;; If pisa-result-type is 'interfaces, return a interface-molecule record or #f
;;
;; If pisa-result-type is neither of the above, return #f
;;
;; a interface-molecule record contains information about pvalue and residues.
;; 
(define (pisa-handle-sxml-molecule imol molecule pisa-results-type)

  ;; was it a chain? (or residue selection disguised as
  ;; a chain?)
  ;; 
  (define (pisa-make-atom-selection-string chain-id-raw)

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

  ;; Return a list of residue records (they can be #f, we will use map)
  ;; 
  ;; indeed they will be #f if bsa is not > 0.1 (or so)
  ;; 
  (define (handle-residues residues make-residue-record)

    (map (lambda (residue)

	   (let ((ser-no #f)
		 (name #f)
		 (seq-num #f)
		 (ins-code #f)
		 (asa #f)
		 (bsa #f)
		 (bonds #f)
		 (solv-en #f))
	     (for-each (lambda (residue-part)
			 (cond
			  ((eq? (car residue-part) 'ser_no)
			   (set! ser-no (cadr residue-part)))
			  ((eq? (car residue-part) 'name)
			   (set! name (cadr residue-part)))
			  ((eq? (car residue-part) 'seq_num)
			   (set! seq-num (cadr residue-part)))
			  ((eq? (car residue-part) 'ins_code)
			   (if (> (length residue-part) 1)
			       (set! ins-code (cadr residue-part))))
			  ((eq? (car residue-part) 'bonds)
			   (if (> (length residue-part) 1)
			       (set! bonds (cadr residue-part))))
			  ((eq? (car residue-part) 'asa)
			   (if (> (length residue-part) 1)
			       (set! asa (cadr residue-part))))
			  ((eq? (car residue-part) 'bsa)
			   (if (> (length residue-part) 1)
			       (set! bsa (cadr residue-part))))
			  ((eq? (car residue-part) 'solv_en)
			   (if (> (length residue-part) 1)
			       (set! solv-en (cadr residue-part))))))
		       (cdr residue)) ;; strip of the initial 'residue
	     (if (not (and ser-no seq-num))
		 #f ;; badness!
		 (if (not (string? bsa))
		     #f
		     (let ((bsa-number (string->number bsa)))
		       (if (< bsa-number 0.1)
			   #f
			   (make-residue-record ser-no name seq-num (if ins-code ins-code "") asa bsa solv-en)))))))
	 (cdr residues))) ;; just the residue list (strip off initial 'residues)
  
  ;; remove #fs from the residue list
  (define (filter-residues residue-list)
    (let ((result '()))
      (cond
       ((null? residue-list) '())
       ((car residue-list) (cons (car residue-list) (filter-residues (cdr residue-list))))
       (else 
	(filter-residues (cdr residue-list))))))

  (define (print-molecule record port)
    (let ((rec-type (record-type-descriptor record)))
      (format port "[molecule: id ~s ~%"       ((record-accessor rec-type 'id)       record))
      (format port " molecule: chain-id ~s ~%" ((record-accessor rec-type 'chain-id) record))
      (format port " molecule: class ~s ~%"    ((record-accessor rec-type 'class)    record))
      (format port " molecule: symop ~s ~%"    ((record-accessor rec-type 'symop-no) record))
      (format port " molecule: symop-no ~s ~%" ((record-accessor rec-type 'symop-no) record))
      (format port " molecule: natoms ~s ~%"   ((record-accessor rec-type 'natoms)   record))
      (format port " molecule: nres ~s ~%"     ((record-accessor rec-type 'nres)     record))
      (format port " molecule: area ~s ~%"     ((record-accessor rec-type 'area)     record))
      (format port " molecule: solv-en ~s ~%"  ((record-accessor rec-type 'solv-en)  record))
      (format port " molecule: pvalue ~s~%"   ((record-accessor rec-type 'pvalue)   record))
      (format port " molecule: with ~s residues]~%"  (length ((record-accessor rec-type 'residues) record)))
      ))


  (define (print-residue record port)
    (let ((rec-type (record-type-descriptor record)))
      (format port "[residue: ser-no ~s,"   ((record-accessor rec-type 'ser-no)   record))
      (format port " name ~s,"     ((record-accessor rec-type 'name)     record))
      (format port " seq-num ~s,"  ((record-accessor rec-type 'seq-num)  record))
      (format port " ins-code ~s," ((record-accessor rec-type 'ins-code) record))
      (format port " asa ~s,"      ((record-accessor rec-type 'asa)      record))
      (format port " bsa ~s,"      ((record-accessor rec-type 'bsa)      record))
      (format port " solv-en ~s]"  ((record-accessor rec-type 'solv-en)  record))))

  (define (make-pisa-molecule-record-type print-molecule)
    (make-record-type "molecule" '(id chain-id class symop symop-no natoms nres area 
				      solv-en pvalue residues)
		      print-molecule))
  ;; 
  (define (make-pisa-residue-record-type print-residue)
    (make-record-type "residue" '(ser-no name seq-num ins-code asa bsa solv-en) print-residue))


  ;; 
  (define (process-molecule-internal)

    ;;
    ;; rtop-symbols are common to interfaces and assemblies
    (let* ((rtop-symbols (list 'rxx 'rxy 'rxz 'ryx 'ryy 'ryz 'rzx 'rzy 'rzz 'tx 'ty 'tz)) ;; orthogonal 
	   ;; matrix elements.
	   ;; these symbols only interfaces have
	   (extra-symbols (list 'pvalue 'residues))
	   ;; 
	   ;; somewhere to store the (interface) pvalue and residues:
	   ;; 
	   (id #f)
	   (chain-id #f)
	   (class #f)
	   (symop #f)
	   (symop-no #f)
	   (natoms #f)
	   (nres #f)
	   (area #f)
	   (solv-en #f)
	   (pvalue #f)
	   (residues '())
	   (rt-mol (make-pisa-molecule-record-type print-molecule))
	   (rt-res (make-pisa-residue-record-type print-residue))
	   (make-molecule-record (record-constructor rt-mol))
	   (make-residue-record  (record-constructor rt-res)))



      (let ((ass-rtop-symbols '()) ;; association
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
			     (pisa-make-atom-selection-string chain-id-raw))
		       (if (string=? chain-id-raw atom-selection-string)
			   (set! symm-name-part (string-append "chain " chain-id-raw))
			   (set! symm-name-part chain-id-raw))))
		 
		 ;; was it one of the rotation or translation symbols?
		 (for-each 
		  (lambda (symbol)
		    (if (eq? (car mol-ele) symbol)
			(let ((n (string->number (cadr mol-ele))))
			  ;; (format #t "~s: ~s~%" symbol n)
			  (set! ass-rtop-symbols 
				(cons (list symbol n)
				      ass-rtop-symbols)))))
		  rtop-symbols)

		 ;; was it an interface symbol?
		 ;; 
		 (for-each
		  (lambda (symbol)
		    (cond 
		     ((eq? (car mol-ele) 'pvalue) (set! pvalue (cadr mol-ele)))
		     ((eq? (car mol-ele) 'id)          (if (> (length mol-ele) 1) (set! id       (cadr mol-ele))))
		     ((eq? (car mol-ele) 'chain_id)    (if (> (length mol-ele) 1) (set! chain-id (cadr mol-ele))))
		     ((eq? (car mol-ele) 'class)       (if (> (length mol-ele) 1) (set! class    (cadr mol-ele))))
		     ((eq? (car mol-ele) 'symop)       (if (> (length mol-ele) 1) (set! symop    (cadr mol-ele))))
		     ((eq? (car mol-ele) 'symop_no)    (if (> (length mol-ele) 1) (set! symop-no (cadr mol-ele))))
		     ((eq? (car mol-ele) 'int_natoms)  (if (> (length mol-ele) 1) (set! natoms   (cadr mol-ele))))
		     ((eq? (car mol-ele) 'int_nres)    (if (> (length mol-ele) 1) (set! nres     (cadr mol-ele))))
		     ((eq? (car mol-ele) 'int_area)    (if (> (length mol-ele) 1) (set! area     (cadr mol-ele))))
		     ((eq? (car mol-ele) 'int_solv_en) (if (> (length mol-ele) 1) (set! solv-en  (cadr mol-ele))))
		     ((eq? (car mol-ele) 'residues) 
		      (let ((filtered-residues (filter-residues (handle-residues mol-ele make-residue-record))))
			(set! residues filtered-residues)))))
		  extra-symbols))))
	 molecule) ;; passed to (at the main level) to pisa-handle-sxml-molecule.
	
	(if (not (= (length ass-rtop-symbols) 12))
	    #f  ;; ass-symbols were not all set
	    (let ((mat 
		   (map (lambda (sym)
			  (cadr (assoc sym ass-rtop-symbols)))
			rtop-symbols)))
;	       (format #t "==================== atom-selection string ~s  mat:::: ~s~%" 
;		       atom-selection-string mat)
	      ;; (format "currently ~s molecules~%" (graphics-n-molecules))
	      (let* ((new-molecule-name (string-append "Symmetry copy of "
						       (number->string imol)
						       " using " symm-name-part))
		     ;; new-molecule-by-symop-with-atom-selection,
		     ;; perhaps? (20100222 doesn't seem needed because
		     ;; the transformation is contained in mat.)
		     (new-mol-no (apply new-molecule-by-symmetry-with-atom-selection
					imol
					new-molecule-name
					atom-selection-string
					(append mat (list 0 0 0)))))

		;; the return value depends on pisa-results-type
		(if (eq? pisa-results-type 'assemblies)
		    ;; assemblies:
		    new-mol-no
		    ;; interfaces:
		    (list (list new-mol-no symop)
			  (make-molecule-record id chain-id class symop symop-no natoms nres area 
						solv-en pvalue residues)))))))))

  ;; main line
  ;; 
  (if (eq? pisa-results-type 'assemblies)
      (process-molecule-internal)
      (if (eq? pisa-results-type 'interfaces)
	  (process-molecule-internal)
	  #f))) ;; neither of the above.
	  

;; ----------------------------------------------------
;;                      pisa assemblies:
;; ----------------------------------------------------
;; 
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
;; 
(define parse-pisa-assemblies
  (lambda (imol entity)


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
		 (let ((mol-no (pisa-handle-sxml-molecule imol ele 'assemblies #f #f)))
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


;; 20100213 prep-for-pisa needs to make directory, config file, write
;; the pdb file and thre return value should be #f if there was a
;; problem or some value where we can check that pisa -analyse ran
;; (probably a directory).  It is not clear right now where the output
;; is going.  config files has PISA_WORK_ROOT coot-pisa but things
;; seems to be written to DATA_ROOT
;; /home/emsley/ccp4/ccp4-6.1.3/share/pisa which seems like a bug (or
;; something like it) in pisa.  Needs rechecking.
;; 
(define (prep-for-pisa mode imol)

  (define (make-stubbed-name imol)
    (strip-extension (basename (molecule-name imol))))
  
  (if (valid-model-molecule? imol) 
      (let* ((pisa-coot-dir "coot-pisa")
	     (stubbed-name (make-stubbed-name imol))
	     (pdb-file-name (append-dir-file pisa-coot-dir (string-append stubbed-name ".pdb")))
	     (pisa-config (append-dir-file pisa-coot-dir (string-append stubbed-name "-pisa.cfg")))
	     (pisa-xml-file-name (append-dir-file pisa-coot-dir (string-append stubbed-name 
									       "-"
									       (symbol->string mode)
									       ".xml"))))
	
	(make-directory-maybe pisa-coot-dir)
	(make-pisa-config pisa-coot-dir pisa-config)
	(write-pdb-file imol pdb-file-name)
	(if (file-exists? pdb-file-name)
	    (values pdb-file-name pisa-config pisa-xml-file-name)
	    (values #f #f #f)))))

	    

;; needs fleshing out (see notes for prep-for-pisa).
;; 
(define (cached-pisa-analysis dir)
  #f)



;;;; -----------------------------------------------------------------------------------
;;;; -----------------------------------------------------------------------------------
;;;;                                   interfaces
;;;; -----------------------------------------------------------------------------------
;;;; -----------------------------------------------------------------------------------

(define (pisa-interfaces imol)

  (receive (pdb-file-name pisa-config pisa-xml-file-name)
	   (prep-for-pisa 'interfaces imol)

	   (if pisa-config
	       (begin
		 (if (not (cached-pisa-analysis pisa-config))
		     ;; pisa analysis
		     (let ((status (goosh-command *pisa-command*
						  (list "pisa"  "-analyse" pdb-file-name pisa-config)
						  '()
						  "pisa-analysis.log" #f)))
					; check status?
		       status))

		 ;; OK, let's do interfaces
		 ;; 
		 (let ((status (goosh-command *pisa-command*
					      (list "pisa" "-xml" "interfaces" pisa-config) '()
					      pisa-xml-file-name #f)))
		   (if (= status 0) ; good
		       (pisa-interfaces-xml imol pisa-xml-file-name)))))))
		     
			 
(define (pisa-interfaces-xml imol file-name)
  (if (not (file-exists? file-name))
      (begin
	(format #t "WARNING: in pisa-interfaces-xml: ~s does not exist~%" 
		file-name))
      (begin
	(call-with-input-file file-name
	  (lambda (port)
	    (let ((sm (xml->sxml port)))
	      (parse-pisa-interfaces imol sm)))))))



;; pdb_entry
;;    pdb_code
;;    status
;;    n_interfaces
;;    interface
;;       id
;;       type
;;       n_occ
;;       int_area
;;       int_solv_en
;;       pvalue
;;       stab_en
;;       css
;;       overlap
;;       x-rel
;;       fixed
;;       h-bonds
;;          n_bonds
;;          bond
;;             chain-1
;;             res-1
;;             seqnum-1
;;             inscode-1
;;             atname-1
;;             chain-2
;;             res-2
;;             seqnum-2
;;             inscode-2
;;             atname-2
;;             dist
;;       salt-bridges
;;          n-bonds
;;          bond
;;             chain-1
;;             res-1
;;             seqnum-1
;;             inscode-1
;;             atname-1
;;             chain-2
;;             res-2
;;             seqnum-2
;;             inscode-2
;;             atname-2
;;             dist
;;       ss-bonds
;;          n-bonds
;;          bond
;;             chain-1
;;             res-1
;;             seqnum-1
;;             inscode-1
;;             atname-1
;;             chain-2
;;             res-2
;;             seqnum-2
;;             inscode-2
;;             atname-2
;;             dist
;;       cov-bonds
;;          n-bonds
;;          bond
;;             chain-1
;;             res-1
;;             seqnum-1
;;             inscode-1
;;             atname-1
;;             chain-2
;;             res-2
;;             seqnum-2
;;             inscode-2
;;             atname-2
;;             dist
;;       molecule
;;          id
;;          chain_id
;;          class
;;          symop
;;          symop_no
;;          cell_i
;;          cell_j
;;          cell_k
;;          rxx
;;          rxy
;;          rxz
;;          tx
;;          ryx
;;          ryy
;;          ryz
;;          ty
;;          rzx
;;          rzy
;;          rzz
;;          tz
;;          int_natoms
;;          int_nres
;;          int_area
;;          int_solv_en
;;          pvalue
;;          residues
;;             residue
;;                ser_no
;;                name
;;                seq_num
;;                ins_code
;;                bonds
;;                asa
;;                bsa
;;                solv_en
;;                
(define (parse-pisa-interfaces imol sxml-entity)

  ;; return a list of bonds:
  ;;
  ;; (bond-type atom-spec-1 atom-spec-2)
  ;; 
  ;; where bond-type is 'hbond, 'salt-bridge, 'ss-bond, or 'cov-bond.
  ;; 
  (define (molecule-bonds entity)

    ;; return a list of 2 atom specs given something like: 
    ;; 
    ;; (bond (chain-1 "A") (res-1 "HIS") (seqnum-1 "268") (inscode-1)
    ;; (atname-1 " NE2") (chain-2 "A") (res-2 "GLU") (seqnum-2 "270")
    ;; (inscode-2) (atname-2 " OE2") (dist "2.2738189697"))
    ;; 
    ;; return #f if the atom specs are not fully set
    ;; 
    (define (parse-bond bond-description)
      (let ((chain-id-1  #f)
	    (chain-id-2  #f)
	    (res-no-1    #f)
	    (res-no-2    #f)
	    (ins-code-1  #f)
	    (ins-code-2  #f)
	    (atom-name-1 #f)
	    (atom-name-2 #f))
	(for-each 
	 (lambda (ele)
	   (cond 
	    ((not (pair? ele)) 'pass)
	    ((not (= (length ele) 2)) 'pass)
	    ((eq? 'chain-1   (car ele)) (set! chain-id-1  (cadr ele)))
	    ((eq? 'chain-2   (car ele)) (set! chain-id-2  (cadr ele)))
	    ((eq? 'seqnum-1  (car ele)) (set! res-no-1    (cadr ele)))
	    ((eq? 'seqnum-2  (car ele)) (set! res-no-2    (cadr ele)))
	    ((eq? 'inscode-1 (car ele)) (set! ins-code-1  (cadr ele)))
	    ((eq? 'inscode-2 (car ele)) (set! ins-code-2  (cadr ele)))
	    ((eq? 'atname-1  (car ele)) (set! atom-name-1 (cadr ele)))
	    ((eq? 'atname-2  (car ele)) (set! atom-name-2 (cadr ele)))))
	 bond-description)
	(if (not (and chain-id-1 chain-id-2 res-no-1 res-no-2 atom-name-1 atom-name-2))
	    #f
	    (list
	     (list chain-id-1 (string->number res-no-1) (if ins-code-1 ins-code-1 "") atom-name-1 "")
	     (list chain-id-2 (string->number res-no-2) (if ins-code-2 ins-code-2 "") atom-name-2 "")))))
      
      

    ;; (format #t "molecule-bonds handle entity: ~s~%" entity)
    (let loop ((elements entity)
	       (bond-type #f)
	       (bonds '()))
      (cond
       ((null? elements) bonds)
       ((not (pair? (car elements)))
	(if (not (member? (car elements) (list 'h-bonds 'salt-bridges 'ss-bonds 'cov-bonds)))
	    (loop (cdr elements) bond-type bonds)
	    (begin
	      (loop (cdr elements) (car elements) bonds))))
       ((eq? (car (car elements)) 'bond)
	(let ((atom-specs (parse-bond (cdr (car elements)))))
	  ;; (format #t "=== got atom specs: ~s from ~s~%" atom-specs (car elements))
	  (if (not atom-specs)
	      (loop (cdr elements) bond-type bonds)
	      (loop (cdr elements)
		bond-type
		(cons (cons bond-type atom-specs) bonds)))))
       (else
	(loop (cdr elements) bond-type bonds)))))
		  
	    

  (define (pisa-handle-sxml-interface interface-entity)
    (let ((molecules '())
	  (bonds '())
	  (area #f)
	  (solv-en #f)
	  (pvalue #f)
	  (stab-en #f))
      (let loop ((interface-entity interface-entity))
	(cond
	 ((null? interface-entity) 
	  (list (reverse molecules) (reverse bonds) area solv-en pvalue stab-en)) ; return this
	 ((not (pair? (car interface-entity))) (loop (cdr interface-entity)))
	 ((eq? (car (car interface-entity)) 'molecule)
	  ;;
	  ;; molecule info is either #f or a pair of a new
	  ;; molecule number and a pisa molecule record (which contains
	  ;; things like pvalue, residues, id, class).
	  ;; 
	  (let ((molecule-info (pisa-handle-sxml-molecule imol (car interface-entity) 'interfaces)))
	    (set! molecules (cons molecule-info molecules))
	    (loop (cdr interface-entity))))
	 ;; or was it a description of the interface bonds?
	 ((member? (car (car interface-entity)) (list 'h-bonds 'salt-bridges 'ss-bonds 'cov-bonds))
	  (set! bonds (append (molecule-bonds (car interface-entity)) bonds))
	  (loop (cdr interface-entity)))
	 ((eq? 'pvalue (car (car interface-entity)))
	  (set! pvalue (string->number (cadr (car interface-entity))))
	  (loop (cdr interface-entity)))
	 ((eq? 'int_area (car (car interface-entity)))
	  (set! area (string->number (cadr (car interface-entity))))
	  (loop (cdr interface-entity)))
	 ((eq? 'int_solv_en (car (car interface-entity)))
	  (set! solv-en (string->number (cadr (car interface-entity))))
	  (loop (cdr interface-entity)))
	 ((eq? 'stab_en (car (car interface-entity)))
	  (set! stab-en (string->number (cadr (car interface-entity))))
	  (loop (cdr interface-entity)))
	 (else 
	  (loop (cdr interface-entity)))))))

  (define (pisa-handle-pdb-entry pdb-entry-entity)
    (let ((interfaces '()))
      (let loop ((entity pdb-entry-entity))
	(cond
	 ((null? entity) 
	  (reverse interfaces)) ;; return this
	 ((not (pair? (car entity))) (loop (cdr entity)))
	 ((eq? (car (car entity)) 'interface)
	  (set! interfaces (cons (pisa-handle-sxml-interface (car entity)) interfaces))
	  (loop (cdr entity)))
	  ;; (format #t "interfaces now: ~s~%" interfaces)
	 (else 
	  (loop (cdr entity)))))))

  ;; 

  ;; main line
  ;; 
  (let ((pdb-entry #f))
    (let loop ((sxml-entity sxml-entity))
      (cond
       ((null? sxml-entity) pdb-entry)
       ((not (pair? sxml-entity)) 'pass)
       ((eq? (car sxml-entity) 'pdb_entry)
	(set! pdb-entry (pisa-handle-pdb-entry (cdr sxml-entity)))
	;; (format #t "DEBUG:: pdb-entry: ~s~%" pdb-entry)
	(handle-pisa-interfaces pdb-entry))
       (else 
	(map loop sxml-entity))))))
      
