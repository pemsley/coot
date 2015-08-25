

(define (mogul-file-name-stub prefix prefix-2 rn chain-id res-no ins-code)
  (string-append "mogul-"
                 (if (not (string? prefix-2))
	             ""
                     (string-append prefix-2 "-"))
		 "-"
		 prefix
		 "-"
		 rn 
		 "-" 
		 (if (string=? chain-id " ") "" chain-id)
		 "-"
		 (number->string res-no)))

(define (write-mogul-ins mogul-ins-file-name mogul-out-file-name sdf-out-file-name mode)

  (call-with-output-file mogul-ins-file-name
    (lambda (port)
      (display (string-append "mogul molecule file " sdf-out-file-name)   port)
      (newline port)
      (display (string-append "mogul output   file " mogul-out-file-name) port)
      (newline port)
      (display "config output format csv" port)
      (newline port)
      (display "config output items fragment_type  atom_indices query_value nhits mean median sd z-score dmin " port)
      (newline port)
      (display "mogul output distribution all on" port)
      (newline port)
      (display "config distribution bond  bin_width 0.01" port)
      (newline port)
      (display "config distribution angle bin_width 0.4" port)
      (newline port)
      (display "config search all filter exclude_solvents" port)
      (newline port)
      (display "config output invalid_fragments exclude" port)
      (newline port)
      (display "bond all" port)
      (newline port)
      (display "angle all" port)
      (newline port)
      (display "torsion all" port)
      (newline port)
      (if (eq? mode 'all)
	  (begin
	    (display "torsion all" port)
	    (newline port)
	    (display "ring all" port)
	    (newline port))))))


;; prefix-str can be blank
;; 
;; return:
;;   a string that is the mogul-out file-name
;;   a symbol: something went wrong
;; 
(define (run-mogul mode imol chain-id res-no ins-code prefix-str use-cache?)

  (let* ((rn (residue-name imol chain-id res-no ins-code)))

    (if (not (string? rn))
	
	'no-residue-ligand-name

	(let* ((stub (string-append (mogul-file-name-stub  (molecule-name-stub imol 0)
							   prefix-str rn chain-id res-no ins-code)))
	       (mogul-ins-file-name (string-append stub ".ins"))
	       (mogul-out-file-name (string-append stub ".out"))
	       (mogul-log-file-name (string-append stub ".log"))
	       (sdf-out-file-name   (string-append stub ".sdf")))
	  (write-mogul-ins mogul-ins-file-name mogul-out-file-name sdf-out-file-name mode)
	  ;; residue-to-sdf-file only works if all of the non-hydrogen atoms
	  ;; in the dictionary are present in the coordinates. So, this will
	  ;; often fail, for example, with covalently-attached ligands (that
	  ;; lose an atom in the attachment).
	  (if (and use-cache? (file-exists? mogul-out-file-name))
	      ;; do this for testing
	      mogul-out-file-name

	      ;; usual happy path
	      (if (not (residue-to-mdl-file-for-mogul imol chain-id res-no ins-code sdf-out-file-name))
		  'bad-sdf-file
		  (let ((args (list "-ins" mogul-ins-file-name)))
		    (let ((status (goosh-command "mogul" args '() mogul-log-file-name #f)))
		      (if (not (ok-goosh-status? status))
			  'bad-run-mogul-status
			  mogul-out-file-name)))))))))
	

;; return a string that is the mogul results file name that gets
;; parsed to make results (return #f on failure)
;; mode is 'all or 'bonds-and-angles
;; 
(define (run-mogul-show-results imol code mode chain-id res-no ins-code)
  (let* ((rn (residue-name imol chain-id res-no ins-code))
	 (stub (string-append (mogul-file-name-stub  (molecule-name-stub imol 0)
						     code rn chain-id res-no ins-code)))
	 (mogul-ins-file-name (string-append stub ".ins"))
	 (mogul-out-file-name (string-append stub ".out"))
	 (sdf-out-file-name   (string-append stub ".sdf")))

    (if (file-exists? mogul-out-file-name)

	;; cached mogul out
	(begin
	  (mogul-markup imol chain-id res-no ins-code mogul-out-file-name)
	  mogul-out-file-name)

; CONFIG OUTPUT FORMAT  CSV
; CONFIG OUTPUT ITEMS fragment_type  atom_labels query_value nhits mean
; median sd z-score dmin
; MOGUL OUTPUT DISTRIBUTION ALL ON
; CONFIG SEARCH ALL FILTER EXCLUDE_ORGANOMETALLICS
; CONFIG SEARCH ALL FILTER EXCLUDE_SOLVENTS
; CONFIG OUTPUT INVALID_FRAGMENTS EXCLUDE


	;; create a new by running mogul.
	(begin 
	  (write-mogul-ins mogul-ins-file-name mogul-out-file-name sdf-out-file-name mode)

	  ;; residue-to-sdf-file only works if all of the non-hydrogen atoms
	  ;; in the dictionary are present in the coordinates. So, this will
	  ;; often fail, for example, with covalently-attached ligands (that
	  ;; lose an atom in the attachment).
	  ;; 
	  (if (not (residue-to-mdl-file-for-mogul imol chain-id res-no ins-code sdf-out-file-name))
	   
	   'bad-sdf-file
	   
	   ;; happy path
	   ;; 
	   ;; sub-thread/semaphore/timeout for this?
	   ;; 
	   (let ((args (list "-ins" mogul-ins-file-name)))
	     (format #t "------ running mogul with ~s ~%" args)
	     (let ((status (goosh-command "mogul" args '() "run-mogul.log" #f)))
	       (if (not (ok-goosh-status? status))
		   'bad-run-mogul-status
		   (begin
		     (mogul-markup imol chain-id res-no ins-code mogul-out-file-name)
		     mogul-out-file-name)))))))))
				 

