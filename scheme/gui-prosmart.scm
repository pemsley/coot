;;;; Copyright 2014 by Medical Research Council

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


;; target is my molecule, ref is the homologous (high-res) model
;; 
(define (run-prosmart imol-target imol-ref include-side-chains?)
  (let ((dir-stub "coot-ccp4"))
    (make-directory-maybe dir-stub)
    (let ((target-pdb-file-name (append-dir-file dir-stub
						 (string-append (molecule-name-stub imol-target 0)
								"-prosmart.pdb")))
	  (reference-pdb-file-name (append-dir-file dir-stub
						    (string-append (molecule-name-stub imol-ref 0)
								   "-prosmart-ref.pdb")))
	  (prosmart-out (append-dir-file "ProSMART_Output"
					 (string-append
					  (coot-replace-string (molecule-name-stub imol-target 0)
							       " " "_" )
					  "-prosmart.txt")))
	  (prosmart-rmax 6.0))

      (write-pdb-file imol-target target-pdb-file-name)
      (write-pdb-file imol-ref reference-pdb-file-name)
      (goosh-command "prosmart"
		     (let ((l (list "-p1" target-pdb-file-name
				    "-p2" reference-pdb-file-name
				    "-restrain_seqid" "30" "-rmax" (number->string prosmart-rmax))))
		       (if include-side-chains?
			   (append l (list "-side"))
			   l))
		     '()
		     (append-dir-file dir-stub "prosmart.log")
		     #f)
      (if (not (file-exists? prosmart-out))
	  (begin
	    (format #t "file not found ~s~%" prosmart-out))
	  (begin
	    (format #t "Reading ProSMART restraints from ~s~%" prosmart-out)
	    (add-refmac-extra-restraints imol-target prosmart-out))))))

(define (add-prosmart-secondard-structure-restraints imol)

  (let ((dir-stub (get-directory "coot-ccp4"))
	(stub-name (molecule-name-stub imol 0)))

    (let* ((helix-pdb-file-name-rwd  (string-append stub-name "-helix.pdb"))
	   (helix-pdb-file-name (append-dir-file dir-stub helix-pdb-file-name-rwd))
	   (strand-pdb-file-name-rwd (string-append stub-name "-strand.pdb"))
	   (strand-pdb-file-name (append-dir-file dir-stub strand-pdb-file-name-rwd))
	   (helix-out (join-dir-file
		       (list dir-stub
			     "ProSMART_Output"
			     (string-append "LIB_" stub-name "-helix" ".txt"))))
	   (strand-out (join-dir-file
			(list
			 dir-stub
			 "ProSMART_Output"
			 (string-append "LIB_" stub-name "-strand" ".txt")))))

      (write-pdb-file imol  helix-pdb-file-name)
      (write-pdb-file imol strand-pdb-file-name)

      ;; Prosmart writes results in ProSMART_Output, so we change to the coot-ccp4 directory
      ;; so that it puts it in the place we expect it

      (let ((current-dir (getcwd)))

	(format #t "dir-stub: ~s~%" dir-stub)
	(chdir dir-stub)
	(format #t "(getcwd): ~s~%" (getcwd))

	(for-each (lambda (p)
		    (goosh-command "prosmart"
				   (list "-p1" (car p) (cadr p))
				   '()
				   (string-append "prosmart-" stub-name "-" (caddr p) ".log")
				   #f))
		  (list (list helix-pdb-file-name-rwd   "-helix" "helix")
			(list strand-pdb-file-name-rwd "-strand" "strand")))

	;; (format #t "helix-out:  ~s~%" helix-out)
	;; (format #t "strand-out: ~s~%" strand-out)

	(chdir current-dir))

      (for-each (lambda (fn)
		  (format #t "INFO:: reading ProSMART output file fn: ~s~%" fn)
		  (if (file-exists? fn)
		      (add-refmac-extra-restraints imol fn)
		      (let ((s (string-append "Missing file: " fn)))
			(info-dialog s))))
		(list helix-out strand-out)))))

(define (add-module-restraints)

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "Restraints")))
	
	(add-simple-coot-menu-menuitem
	 menu "Generate Self Restraints 3.7 for Chain"
	 (lambda ()
	   (using-active-atom
	    (generate-local-self-restraints aa-imol aa-chain-id 3.7))))

	(add-simple-coot-menu-menuitem
	 menu "Generate Self Restraints 4.3 for Chain"
	 (lambda ()
	   (using-active-atom
	    (generate-local-self-restraints aa-imol aa-chain-id 4.3))))

        ;; note to self: make this a loop next time
	(add-simple-coot-menu-menuitem
	 menu "Generate Self Restraints 6 for Chain"
	 (lambda ()
	   (using-active-atom
	    (generate-local-self-restraints aa-imol aa-chain-id 6))))

	(add-simple-coot-menu-menuitem
	 menu "Generate All-Molecule Self Restraints 4.3"
	 (lambda ()
	   (using-active-atom
	    (generate-self-restraints aa-imol 4.3))))

	(add-simple-coot-menu-menuitem
	 menu "Generate All-Molecule Self Restraints 5.0"
	 (lambda ()
	   (using-active-atom
	    (generate-self-restraints aa-imol 5.0))))

	(add-simple-coot-menu-menuitem
	 menu "Generate Local Self Restraints 6"
	 (lambda ()
	   (using-active-atom

            (let* ((centred-residue (list-head (cdr active-atom) 3))
		   (radius 10)
		   (local-dist-max 4.2)
                   (imol (car active-atom))
                   (other-residues (residues-near-residue imol centred-residue radius))
                   (residue-specs (if (list? other-residues)
				      (cons centred-residue other-residues)
				      (list centred-residue))))

	      (generate-local-self-restraints-by-residues-scm aa-imol residue-specs local-dist-max)))))


	(add-simple-coot-menu-menuitem
	 menu "Undisplay Extra Restraints"
	 (lambda ()
	   (using-active-atom
	    (set-show-extra-restraints aa-imol 0))))

	(add-simple-coot-menu-menuitem
	 menu "Display Extra Restraints"
	 (lambda ()
	   (using-active-atom
	    (set-show-extra-restraints aa-imol 1))))

	(add-simple-coot-menu-menuitem
	 menu "Show Only Deviant Distances Beyond 6"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -6 6))))

	(add-simple-coot-menu-menuitem
	 menu "Show Only Deviant Distances Beyond 4"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -4 4))))

	(add-simple-coot-menu-menuitem
	 menu "Show Only Deviant Distances Beyond 2.0"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -2.0 2.0))))

	(add-simple-coot-menu-menuitem
	 menu "Show Only Deviant Distances Beyond 1.0"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol -1.0 1.0))))

	(add-simple-coot-menu-menuitem
	 menu "Undisplay All Extra Distance Restraints"
	 (lambda ()
	   (using-active-atom
	    (set-extra-restraints-prosmart-sigma-limits aa-imol 0 0 ))))

;	(add-simple-coot-menu-menuitem
;	 menu "Restraint Representation To CA"
;	 (lambda ()
;	   (using-active-atom
;	    (set-extra-restraints-representation-for-bonds-go-to-CA aa-imol 1))))

;	(add-simple-coot-menu-menuitem
;	 menu "Restraint Representation To Home Atom"
;	 (lambda ()
;	   (using-active-atom
;	    (set-extra-restraints-representation-for-bonds-go-to-CA aa-imol 0))))

	(load-by-search "user-define-restraints.scm")

	(add-simple-coot-menu-menuitem
	 menu "Delete All Extra Restraints"
	 (lambda ()
	   (using-active-atom
	    (delete-all-extra-restraints aa-imol)))))))



(define (add-module-prosmart)

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "ProSMART")))

	(add-simple-coot-menu-menuitem
	 menu "Add ProSMART Secondary Structure Restraints"
	 (lambda ()
	   (using-active-atom
	    (add-prosmart-secondard-structure-restraints aa-imol))))

	(add-simple-coot-menu-menuitem
	 menu "ProSMART..."
	 (lambda ()
	   (let ((window (gtk-window-new 'toplevel))
		 (hbox (gtk-hbox-new #f 0))
		 (vbox (gtk-vbox-new #f 0))
		 (h-sep (gtk-hseparator-new))
		 (chooser-hint-text-1 " Target molecule ")
		 (chooser-hint-text-2 " Reference (high-res) molecule ")
		 (go-button (gtk-button-new-with-label " ProSMART "))
		 (cancel-button (gtk-button-new-with-label " Cancel "))
		 (check-button (gtk-check-button-new-with-label "Include Side-chains")))

	     (let ((option-menu-mol-list-pair-tar (generic-molecule-chooser
						   vbox chooser-hint-text-1))
		   (option-menu-mol-list-pair-ref (generic-molecule-chooser
						   vbox chooser-hint-text-2)))

	       (gtk-box-pack-start vbox check-button  #f #f 2)
	       (gtk-box-pack-start vbox h-sep         #f #f 2)
	       (gtk-box-pack-start vbox hbox          #f #f 2)
	       (gtk-box-pack-start hbox go-button     #f #f 6)
	       (gtk-box-pack-start hbox cancel-button #f #f 6)
	       (gtk-container-add window vbox)

	       (gtk-signal-connect cancel-button "clicked"
				   (lambda ()
				     (gtk-widget-destroy window)))

	       (gtk-signal-connect go-button "clicked"
				   (lambda ()
				     (let ((imol-tar
					    (apply get-option-menu-active-molecule
						   option-menu-mol-list-pair-tar))
					   (imol-ref
					    (apply get-option-menu-active-molecule
						   option-menu-mol-list-pair-ref))
					   (do-side-chains? (gtk-toggle-button-get-active check-button)))
				       (run-prosmart imol-tar imol-ref do-side-chains?)
				       (gtk-widget-destroy window))))
	       (gtk-widget-show-all window)))))

	)))

