
(define (add-module-prosmart) 

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "ProSMART")))
	
	(add-simple-coot-menu-menuitem 
	 menu "Cut to 6"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -6 6))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 4"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -4 4))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 2.5"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -2.5 2.5))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 2.0"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -2.0 2.0))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 1.5"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -1.5 1.5))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 1.0"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -1.0 1.0))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 0.5"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -0.5 0.5))))

	(add-simple-coot-menu-menuitem 
	 menu "Cut to 0"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol 0 0 ))))

	(add-simple-coot-menu-menuitem
	 menu "Restraint Representation To CA"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-representation-for-bonds-go-to-CA aa-imol 1))))
        
	(add-simple-coot-menu-menuitem
	 menu "Restraint Representation To Home Atom"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-representation-for-bonds-go-to-CA aa-imol 0))))

	(add-simple-coot-menu-menuitem
	 menu "Generate Self Restraints 4.3"
	 (lambda ()
	   (using-active-atom
	    (generate-local-self-restraints aa-imol aa-chain-id 4.3))))

	(add-simple-coot-menu-menuitem
	 menu "Generate Self Restraints 6"
	 (lambda ()
	   (using-active-atom
	    (generate-local-self-restraints aa-imol aa-chain-id 6))))
        
	)))

