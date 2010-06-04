

(define *cprodrg* "cprodrg")
(define *cprodrg* "/home/paule/ccp4/ccp4-6.1.2/bin/cprodrg")

(define (import-from-prodrg minimize-mode)

  (let ((prodrg-dir "coot-ccp4")
	(res-name "DRG"))
    
    (make-directory-maybe prodrg-dir)
    (let ((prodrg-xyzin  "../../coot/lbg/prodrg-in.mdl")
	  (prodrg-xyzout (append-dir-file prodrg-dir
					  (string-append "prodrg-" res-name ".pdb")))
	  (prodrg-cif    (append-dir-file prodrg-dir 
					  (string-append "prodrg-out.cif")))
	  (prodrg-log    (append-dir-file prodrg-dir "prodrg.log")))
      (let* ((mini-mode (if (eq? minimize-mode 'mini-no) "NO" "PREP"))
	     (status 
	      (goosh-command *cprodrg*
			     (list "XYZIN"  prodrg-xyzin
				   "XYZOUT" prodrg-xyzout
				   "LIBOUT" prodrg-cif)
			     (list (string-append "MINI " mini-mode) "END")
			     prodrg-log #t)))
	(if (number? status)
	    (if (= status 0)
		(begin
		  (read-cif-dictionary prodrg-cif)
		  (let ((imol (handle-read-draw-molecule-and-move-molecule-here prodrg-xyzout)))

		    (using-active-atom 
		     (match-ligand-torsions imol aa-imol aa-chain-id aa-res-no))

		    ;; question, can we overlap this new molecule on
		    ;; the residue that we are sitting on?
		    ;; 
		    (using-active-atom
		     (overlap-ligands imol aa-imol aa-chain-id aa-res-no))
		    
		    (with-auto-accept
		     (regularize-residues imol (list (list "" 1 ""))))
		    
		    (additional-representation-by-attributes imol "" 1 1 "" 2 2 0.2 1)

		    ))))))))




(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "PRODRG")))
      (add-simple-coot-menu-menuitem
       menu "Import (using MINI PREP)" 
       (lambda () 
	 ;; run prodrg, read its output files, and run regularisation
	 ;; on the imported PDB file.
	 (import-from-prodrg 'mini-prep)))

      (add-simple-coot-menu-menuitem
       menu "Import (no pre-minimisation)" 
       (lambda () 
	 ;; run prodrg, read its output files, and run regularisation
	 ;; on the imported PDB file.
	 (import-from-prodrg 'mini-no)))

      (add-simple-coot-menu-menuitem
       menu "Export to lbg"
       (lambda ()
	 (using-active-atom 
	  (prodrg-flat aa-imol aa-chain-id aa-res-no))))

      (add-simple-coot-menu-menuitem
       menu "FLE-View"
       (lambda ()
	 (using-active-atom (fle-view aa-imol aa-chain-id aa-res-no aa-ins-code))))

      (add-simple-coot-menu-menuitem
       menu "Load SBase monomer..."
       (lambda ()
	 (generic-single-entry "Load SBase Monomer from three-letter-code: " ""
			       " Load "
			       (lambda (tlc)
				 (get-sbase-monomer tlc)))))))

		  
			   


;; if there is a prodrg-xyzin set the current-time to its mtime, else #f
;; 
(define prodrg-xyzin       "../../coot/lbg/prodrg-in.mdl")
(define sbase-to-coot-tlc  "../../coot/lbg/.sbase-to-coot-comp-id")

(define (get-mdl-latest-time file-name)
  (if (not (file-exists? file-name))
      #f
      (stat:mtime (stat file-name))))

(let ((mdl-latest-time (get-mdl-latest-time prodrg-xyzin))
      (sbase-transfer-latest-time (get-mdl-latest-time sbase-to-coot-tlc)))
  (let ((func (lambda ()
		(let ((mdl-now-time (get-mdl-latest-time prodrg-xyzin))
		      (sbase-now-time (get-mdl-latest-time sbase-to-coot-tlc)))

;		  (format #t "sbase-now-time ~s   sbase-transfer-latest-time ~s~%" 
;			  sbase-now-time sbase-transfer-latest-time)

		  (if (number? mdl-now-time)
		      (if (> mdl-now-time mdl-latest-time)
			  (begin
			    (set! mdl-latest-time mdl-now-time)
			    (import-from-prodrg 'mini-prep))))

		  (if (number? sbase-now-time)
		      (if (> sbase-now-time sbase-transfer-latest-time)
			  (begin
			    (set! sbase-transfer-latest-time sbase-now-time)
			    (let ((tlc-symbol 
				   (call-with-input-file sbase-to-coot-tlc
				     (lambda (port)
				       (read-line port)))))
			      (let ((imol (get-sbase-monomer tlc-symbol)))
				(if (not (valid-model-molecule? imol))
				    (format #t "failed to get SBase molecule for ~s~%"
					    tlc-symbol))))))))
		
		#t))) ;; return value, keep running
    (gtk-timeout-add 500 func)))


;; return #f (if fail) or a list of: the molecule number of the
;; selected residue, the prodrg output mol file-name, the prodrg
;; output pdb file-name
;; 
(define (prodrg-flat imol-in chain-id-in res-no-in)

  (let* ((selection-string (string-append "//" chain-id-in "/" (number->string res-no-in)))
	 (imol (new-molecule-by-atom-selection imol-in selection-string))
	 (prodrg-input-file-name (append-dir-file "coot-ccp4"
						  "tmp-residue-for-prodrg.pdb"))
	 (prodrg-output-mol-file (append-dir-file "coot-ccp4" ".coot-to-lbg-mol"))
	 (prodrg-output-pdb-file (append-dir-file "coot-ccp4" ".coot-to-lbg-pdb"))
	 (prodrg-output-lib-file (append-dir-file "coot-ccp4" ".coot-to-lbg-lib"))
	 (prodrg-log (append-dir-file "coot-ccp4" "tmp-prodrg-flat.log")))
    (set-mol-displayed imol 0)
    (set-mol-active    imol 0)
    (write-pdb-file imol prodrg-input-file-name)
    (let* ((arg-list (list "XYZIN"  prodrg-input-file-name
			   "MOLOUT" prodrg-output-mol-file
			   "XYZOUT" prodrg-output-pdb-file
			   "LIBOUT" prodrg-output-lib-file))
	   (nov (format #t "arg-list: ~s~%" arg-list))
	   (status 
	    (goosh-command *cprodrg*
			   arg-list
			   (list "COORDS BOTH" "MINI FLAT" "END")
			    prodrg-log #t)))
      (if (not (number? status))
	  (begin
	    (info-dialog "Ooops: cprodrg not found?")
	    #f)
	  (if (not (= status 0))
	      (let ((mess (string-append "Something went wrong running cprodrg\n"
					 (if (< (random 100) 10)
					     "(quelle surprise)" ""))))
		(info-dialog mess)
		#f)

	      ;; normal return value (hopefully)
	      (list imol 
		    prodrg-output-mol-file 
		    prodrg-output-pdb-file 
		    prodrg-output-lib-file))))))



(define (fle-view imol chain-id res-no ins-code)

  (using-active-atom
   (let ((imol aa-imol)
	 (chain-id aa-chain-id)
	 (res-no aa-res-no))
     (let ((r (prodrg-flat imol chain-id res-no)))
       (if r
	   (let ((imol-ligand-fragment (car r))
		 (prodrg-output-flat-mol-file-name (list-ref r 1))
		 (prodrg-output-flat-pdb-file-name (list-ref r 2))
		 (prodrg-output-cif-file-name      (list-ref r 3)))
	     (fle-view-internal imol chain-id res-no "" ;; from active atom
	      imol-ligand-fragment
	      prodrg-output-flat-mol-file-name
	      prodrg-output-flat-pdb-file-name
	      prodrg-output-cif-file-name)
	     (goosh-command "touch" 
			    (list
			     (append-dir-file "coot-ccp4" ".coot-to-lbg-mol-ready"))
			    '() "/dev/null" #f)))))))

