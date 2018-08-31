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


(define (add-module-prosmart) 

  (if (defined? 'coot-main-menubar)
      (let ((menu (coot-menubar-menu "ProSMART")))
	
	(add-simple-coot-menu-menuitem
	 menu "Generate Self Restraints 4.3 for Chain"
	 (lambda ()
	   (using-active-atom
	    (generate-local-self-restraints aa-imol aa-chain-id 4.3))))

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

	)))

