

(define *cprodrg* "cprodrg")
;; (define *cprodrg* "/home/paule/ccp4/ccp4-6.1.2/bin/cprodrg")

;; if there is a prodrg-xyzin set the current-time to its mtime, else #f
;; 
(define prodrg-xyzin       "prodrg-in.mdl")
(define sbase-to-coot-tlc  ".sbase-to-coot-comp-id")


;; what is this rubbish?
;; (define prodrg-xyzin       "/home/emsley/build-coot-f12/coot/lbg/prodrg-in.mdl")
;; (define sbase-to-coot-tlc  "/home/emsley/build-coot-f12/coot/lbg/.sbase-to-coot-comp-id")


;; (define prodrg-xyzin       "/lmb/wear/emsley/Projects/coot/lbg/prodrg-in.mdl")
;; (define sbase-to-coot-tlc  "/lmb/wear/emsley/Projects/coot/lbg/.sbase-to-coot-comp-id")



(define (import-from-prodrg minimize-mode)

  ;; return #t or #f
  ;; 
  (define (have-restraints-for? res-name)
    (let ((restraints (monomer-restraints res-name)))
      (if restraints #t #f)))


  ;; overlap the imol-ligand residue if there are restraints for the
  ;; reference residue/ligand.
  ;; 
  (define (overlap-ligands-maybe imol-ligand imol-ref chain-id-ref res-no-ref)

    ;; we don't want to overlap-ligands if there is no dictionary
    ;; for the residue to be matched to
    (let* ((res-name (residue-name imol-ref chain-id-ref res-no-ref ""))
	   (restraints (monomer-restraints res-name)))
      (if (not restraints)
	  #f
	  (begin 
	    (overlap-ligands imol-ligand imol-ref chain-id-ref res-no-ref)
	    #t))))
		     

  ;; return the new molecule number
  ;; 
  (define (read-and-regularize prodrg-xyzout)
    (let ((imol (handle-read-draw-molecule-and-move-molecule-here prodrg-xyzout)))
      (with-auto-accept
       ;; speed up the minisation (and then restore setting).
       (let ((s (dragged-refinement-steps-per-frame)))
	 (set-dragged-refinement-steps-per-frame 500)
	 (regularize-residues imol (list (list "" 1 "")))
	 (set-dragged-refinement-steps-per-frame s)))
      imol)) 

  ;; return the new molecule number
  ;; (only works with aa-ins-code of "")
  (define (read-regularize-and-match-torsions prodrg-xyzout aa-imol aa-chain-id aa-res-no)
    (let ((imol (handle-read-draw-molecule-and-move-molecule-here prodrg-xyzout)))

      (if (have-restraints-for? (residue-name aa-imol aa-chain-id aa-res-no ""))
	  (overlap-ligands-maybe imol aa-imol aa-chain-id aa-res-no))

      (with-auto-accept
       ;; speed up the minisation (and then restore setting).
       (let ((s (dragged-refinement-steps-per-frame)))
	 (set-dragged-refinement-steps-per-frame 500)
	 (regularize-residues imol (list (list "" 1 "")))
	 (if (have-restraints-for? (residue-name aa-imol aa-chain-id aa-res-no ""))
	     (match-ligand-torsions imol aa-imol aa-chain-id aa-res-no))
	 (regularize-residues imol (list (list "" 1 "")))
	 (set-dragged-refinement-steps-per-frame s)))
      imol))


  ;; main line
  ;; 
  (let ((prodrg-dir "coot-ccp4")
	(res-name "DRG"))
    
    (make-directory-maybe prodrg-dir)
    (let ((prodrg-xyzout (append-dir-file prodrg-dir
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

		  ;; OK, so here we read the PRODRG files and
		  ;; manipulate them.  We presumen that the active
		  ;; residue is quite like the input ligand from
		  ;; prodrg.
		  ;;
		  ;; Read in the lib and coord output of PRODRG.  Then
		  ;; overalay the new ligand onto the active residue
		  ;; (just so that we can see it approximately
		  ;; oriented). Then match the torsions from the new
		  ;; ligand to the those of the active residue.  Then
		  ;; overlay again so that we have the best match.
		  ;;
		  ;; We want to see just one molecule with the protein
		  ;; and the new ligand.
		  ;; add-ligand-delete-residue-copy-molecule provides
		  ;; that for us.  We just colour it and undisplay the
		  ;; other molecules.

		  (read-cif-dictionary prodrg-cif)

		  ;; we do different things depending on whether or
		  ;; not there is an active residue.  We need to test
		  ;; for having an active residue here (currently we
		  ;; presume that there is).
		  ;; 
		  ;; Similarly, if the aa-ins-code is non-null, let's
		  ;; presume that we do not have an active residue.

		  (let ((active-atom (active-residue)))

		    (if (or (not active-atom)
			    (not (string=? (list-ref active-atom 3) ""))) ;; aa-ins-code
			
			;; then there is no active residue to match to 
			(begin 
			  (read-and-regularize prodrg-xyzout))
			
			;; we have an active residue to match to 
			(using-active-atom
			 (let ((imol (read-regularize-and-match-torsions 
				      prodrg-xyzout aa-imol aa-chain-id aa-res-no )))
			   
			   (let ((overlapped-flag
				  (overlap-ligands-maybe imol aa-imol aa-chain-id aa-res-no)))
			   
			     ;; (additional-representation-by-attributes imol "" 1 1 "" 2 2 0.2 1)
			     
			     (if overlapped-flag
				 (set-mol-displayed aa-imol 0))
			     (set-mol-displayed imol 0)
			     (let* ((col (get-molecule-bonds-colour-map-rotation aa-imol))
				    (new-col (+ col 5)) ;; a tiny amount.
				    (imol-replaced 
				     ;; new ligand specs, then "reference" ligand (to be deleted)
				     (add-ligand-delete-residue-copy-molecule imol "" 1 
									      aa-imol aa-chain-id aa-res-no)))
			       (set-molecule-bonds-colour-map-rotation imol-replaced new-col)
			       (graphics-draw))))
			 #t
			 ))))))))))




(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Lidia")))
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

;      (add-simple-coot-menu-menuitem
;       menu "Export to LIDIA"
;       (lambda ()
;	 (using-active-atom 
;	  (prodrg-flat aa-imol aa-chain-id aa-res-no))))

      (add-simple-coot-menu-menuitem
       menu "View in LIDIA"
       (lambda ()
	 (using-active-atom (fle-view aa-imol aa-chain-id aa-res-no aa-ins-code))))

      (add-simple-coot-menu-menuitem
       menu "Load SBase monomer..."
       (lambda ()
	 (generic-single-entry "Load SBase Monomer from three-letter-code: " ""
			       " Load "
			       (lambda (tlc)
				 (get-sbase-monomer tlc)))))))

		  
			   



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
		      (if (number? mdl-latest-time)
			  (if (> mdl-now-time mdl-latest-time)
			      (begin
				(set! mdl-latest-time mdl-now-time)
				(import-from-prodrg 'mini-prep)))))

		  (if (number? sbase-transfer-latest-time)
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
						tlc-symbol)
					
					;; it was read OK, do an overlap:
					(using-active-atom
					 (overlap-ligands imol aa-imol aa-chain-id aa-res-no))
					
					))))))))
		
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


(define (prodrg-plain mode imol-in chain-id-in res-no-in)
  
  (let* ((selection-string (string-append "//" chain-id-in "/" 
					  (number->string res-no-in)))
	 (imol (new-molecule-by-atom-selection imol-in selection-string)))
    (let* ((stub (append-dir-file "coot-ccp4" (string-append 
					       "prodrg-tmp-" 
					       (number->string (getpid)))))
	   (prodrg-xyzin (string-append stub "-xyzin.pdb"))
	   (prodrg-xyzout (string-append stub "-xyzout.pdb"))
	   (prodrg-cif  (string-append stub "-dict.cif"))
	   (prodrg-log  (string-append stub ".log")))
      
      (write-pdb-file imol prodrg-xyzin)
      (let ((result (goosh-command "cprodrg" 
				   (list "XYZIN"  prodrg-xyzin
					 "XYZOUT" prodrg-xyzout
					 "LIBOUT" prodrg-cif)
				   (list "MINI PREP" "END")
				   prodrg-log #t)))
	(close-molecule imol)
	(list result prodrg-xyzout prodrg-cif)))))
    


(define (fle-view imol chain-id res-no ins-code)

  (using-active-atom
   (let ((imol aa-imol)
	 (chain-id aa-chain-id)
	 (res-no aa-res-no))
     (let ((r-flat  (prodrg-flat  imol chain-id res-no))
	   (r-plain (prodrg-plain 'mini-no  imol chain-id res-no)))
           ;; (r-plain (list 0 0))) ;; dummy value
       (if (and r-flat
		(and (number? (car r-plain))
		     (= (car r-plain) 0)))
	   (let ((imol-ligand-fragment (car r-flat))
		 (prodrg-output-flat-mol-file-name (list-ref r-flat  1))
		 (prodrg-output-flat-pdb-file-name (list-ref r-flat  2))
		 (prodrg-output-cif-file-name      (list-ref r-flat  3))
		 (prodrg-output-3d-pdb-file-name   (list-ref r-plain 1)))
	     (fle-view-internal imol chain-id res-no "" ;; from active atom
	      imol-ligand-fragment
	      prodrg-output-flat-mol-file-name
	      prodrg-output-flat-pdb-file-name
	      prodrg-output-3d-pdb-file-name
	      prodrg-output-cif-file-name)
	     (goosh-command "touch" 
			    (list
			     (append-dir-file "coot-ccp4" ".coot-to-lbg-mol-ready"))
			    '() "/dev/null" #f)))))))

