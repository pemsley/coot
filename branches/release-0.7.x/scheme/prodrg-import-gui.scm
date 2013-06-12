(if (defined? 'coot-main-menubar)
    (let ((menu (coot-menubar-menu "Lidia")))

;   Not for public use.
;
;       (add-simple-coot-menu-menuitem
;        menu "Import (using MINI PREP)" 
;        (lambda () 
; 	 ;; run prodrg, read its output files, and run regularisation
; 	 ;; on the imported PDB file.
; 	 (import-from-prodrg 'mini-prep)))

;   Not for public use.
; 
;       (add-simple-coot-menu-menuitem
;        menu "Import (no pre-minimisation)" 
;        (lambda () 
; 	 ;; run prodrg, read its output files, and run regularisation
; 	 ;; on the imported PDB file.
; 	 (import-from-prodrg 'mini-no)))

      (add-simple-coot-menu-menuitem
       menu "Hydrogenate region"
       (lambda () 
	 (hydrogenate-region 6)))
      
;   This doesn't work at the moment - let's activate it later...
;
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
				 (get-sbase-monomer tlc)))))

      (add-simple-coot-menu-menuitem
       menu "Activate prodrg flat mode"
       (lambda ()
	 (using-active-atom 
	  (prodrg-flat aa-imol aa-chain-id aa-res-no))))))

