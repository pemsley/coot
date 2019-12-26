;;;; Copyright 2013 by the University of Oxford
;;;; Copyright 2013 by Medical Research Council

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

(define (res-name-from-atom-spec atom-spec)

  (let ((imol (cadr atom-spec)))
    (residue-name imol 
		  (list-ref atom-spec 2)
		  (list-ref atom-spec 3)
		  (list-ref atom-spec 4))))


(define (user-defined-add-single-bond-restraint)

  (add-status-bar-text "Click on 2 atoms to define the additional bond restraint")
  (user-defined-click 2 
		      (lambda (atom-specs)
			(if (= (cadr (list-ref atom-specs 0))
			       (cadr (list-ref atom-specs 1)))
			    (let ((imol (cadr (list-ref atom-specs 0))))
			      (format #t "imol: ~s   spec-1: ~s   spec-2: ~s~% "
				      imol
				      (list-ref atom-specs 0)
				      (list-ref atom-specs 1))
			      
			      (add-extra-bond-restraint imol
							(list-ref (list-ref atom-specs 0) 2)
							(list-ref (list-ref atom-specs 0) 3)
							(list-ref (list-ref atom-specs 0) 4)
							(list-ref (list-ref atom-specs 0) 5)
							(list-ref (list-ref atom-specs 0) 6)
							(list-ref (list-ref atom-specs 1) 2)
							(list-ref (list-ref atom-specs 1) 3)
							(list-ref (list-ref atom-specs 1) 4)
							(list-ref (list-ref atom-specs 1) 5)
							(list-ref (list-ref atom-specs 1) 6)
							1.54 0.05))))))

(define (user-defined-add-arbitrary-length-bond-restraint)
  (generic-single-entry "Add a User-defined extra distance restraint"
			"2.8"
			"OK..."
			(lambda (text) 
			  (let ((s  "Now click on 2 atoms to define the additional bond restraint"))
			    (add-status-bar-text s))
			  (let ((bl (string->number text)))
			    (if (not (number? bl))
				(add-status-bar-text "Must define a number for the bond length")

				(user-defined-click 
				 2 
				 (lambda (atom-specs)
				   (if (= (cadr (list-ref atom-specs 0))
					  (cadr (list-ref atom-specs 1)))
				       (let ((imol (cadr (list-ref atom-specs 0))))
					 (format #t "imol: ~s   spec-1: ~s   spec-2: ~s~% "
						 imol
						 (list-ref atom-specs 0)
						 (list-ref atom-specs 1))
					 
					 (add-extra-bond-restraint
					  imol
					  (list-ref (list-ref atom-specs 0) 2)
					  (list-ref (list-ref atom-specs 0) 3)
					  (list-ref (list-ref atom-specs 0) 4)
					  (list-ref (list-ref atom-specs 0) 5)
					  (list-ref (list-ref atom-specs 0) 6)
					  (list-ref (list-ref atom-specs 1) 2)
					  (list-ref (list-ref atom-specs 1) 3)
					  (list-ref (list-ref atom-specs 1) 4)
					  (list-ref (list-ref atom-specs 1) 5)
					  (list-ref (list-ref atom-specs 1) 6)
					  bl 0.01))))))))))

;; spec-1 and spec-2 are 7-element atom-specs
;; 
(define (add-base-restraint imol spec-1 spec-2 atom-name-1 atom-name-2 dist)
  (format #t "add-base-restraint ~s ~s ~s ~s ~s ~s~%" imol spec-1 spec-2 atom-name-1 atom-name-2 dist)
  (add-extra-bond-restraint imol
			    (list-ref spec-1 2)
			    (list-ref spec-1 3)
			    (list-ref spec-1 4)
			    atom-name-1
			    (list-ref spec-1 6)
			    (list-ref spec-2 2)
			    (list-ref spec-2 3)
			    (list-ref spec-2 4)
			    atom-name-2
			    (list-ref spec-2 6)
			    dist 0.035))

(define (a-u-restraints spec-1 spec-2)

  (let ((imol (cadr spec-1)))
    (add-base-restraint imol spec-1 spec-2 " N6 " " O4 " 3.12)
    (add-base-restraint imol spec-1 spec-2 " N1 " " N3 " 3.05)
    (add-base-restraint imol spec-1 spec-2 " C2 " " O2 " 3.90)
    (add-base-restraint imol spec-1 spec-2 " C2 " " O2 " 3.90)
    (add-base-restraint imol spec-1 spec-2 " N3 " " O2 " 5.12)
    (add-base-restraint imol spec-1 spec-2 " C6 " " O4 " 3.92)
    (add-base-restraint imol spec-1 spec-2 " C4 " " C6 " 8.38)))

(define (g-c-restraints spec-1 spec-2)
    
  (let ((imol (cadr spec-1)))
    (add-base-restraint imol spec-1 spec-2 " O6 " " N4 " 3.08)
    (add-base-restraint imol spec-1 spec-2 " N1 " " N3 " 3.04)
    (add-base-restraint imol spec-1 spec-2 " N2 " " O2 " 3.14)
    (add-base-restraint imol spec-1 spec-2 " C4 " " N1 " 7.73)
    (add-base-restraint imol spec-1 spec-2 " C5 " " C5 " 7.21)))

(define (dna-a-t-restraints spec-1 spec-2)

  (let ((imol (cadr spec-1)))
    (add-base-restraint imol spec-1 spec-2 " C2 " " O2 " 3.49)
    (add-base-restraint imol spec-1 spec-2 " N1 " " N3 " 2.85)
    (add-base-restraint imol spec-1 spec-2 " N6 " " O4 " 3.23)
    (add-base-restraint imol spec-1 spec-2 " C6 " " C8 " 9.94)
    (add-base-restraint imol spec-1 spec-2 " N1 " " O2 " 3.52)))

(define (dna-g-c-restraints spec-1 spec-2)
    
  (let ((imol (cadr spec-1)))
    (add-base-restraint imol spec-1 spec-2 " O6 " " N4 " 2.72)
    (add-base-restraint imol spec-1 spec-2 " N1 " " N3 " 2.81)
    (add-base-restraint imol spec-1 spec-2 " N2 " " O2 " 2.83)
    (add-base-restraint imol spec-1 spec-2 " N9 " " N1 " 8.83)))



(define (user-defined-RNA-A-form)
  (user-defined-click 2
		      (lambda (atom-specs)
			(let ((spec-1 (list-ref atom-specs 0))
			      (spec-2 (list-ref atom-specs 1)))
			  (let ((res-name-1 (res-name-from-atom-spec spec-1))
				(res-name-2 (res-name-from-atom-spec spec-2)))

			    (if (and (string=? res-name-1 "G")
				     (string=? res-name-2 "C"))
				(g-c-restraints spec-1 spec-2))

			    (if (and (string=? res-name-1 "C")
				     (string=? res-name-2 "G"))
				(g-c-restraints spec-2 spec-1))

			    (if (and (string=? res-name-1 "A")
				     (string=? res-name-2 "U"))
				(a-u-restraints spec-1 spec-2))

			    (if (and (string=? res-name-1 "U")
				     (string=? res-name-2 "A"))
				(a-u-restraints spec-2 spec-1)))))))

(define (user-defined-DNA-B-form)
  (user-defined-click 2
		      (lambda (atom-specs)
			(let ((spec-1 (list-ref atom-specs 0))
			      (spec-2 (list-ref atom-specs 1)))
			  (let ((res-name-1 (res-name-from-atom-spec spec-1))
				(res-name-2 (res-name-from-atom-spec spec-2)))

			    (if (and (string=? res-name-1 "DG")
				     (string=? res-name-2 "DC"))
				(dna-g-c-restraints spec-1 spec-2))

			    (if (and (string=? res-name-1 "DC")
				     (string=? res-name-2 "DG"))
				(dna-g-c-restraints spec-2 spec-1))

			    (if (and (string=? res-name-1 "DA")
				     (string=? res-name-2 "DT"))
				(dna-a-t-restraints spec-1 spec-2))

			    (if (and (string=? res-name-1 "DT")
				     (string=? res-name-2 "DA"))
				(dna-a-t-restraints spec-2 spec-1)))))))

(define (user-defined-add-helix-restraints)
    (user-defined-click 
     2 (lambda (atom-specs)
	 (let* ((spec-1 (list-ref atom-specs 0))
		(spec-2 (list-ref atom-specs 1))
		(chain-id-1 (list-ref spec-1 2))
		(chain-id-2 (list-ref spec-2 2))
		(resno-1 (list-ref spec-1 3))
		(resno-2 (list-ref spec-2 3))
		(imol (cadr spec-1)))


	   (if (string=? chain-id-1 chain-id-2)
	       (begin
		 ;; if backwards, swap them
		 (if (< resno-2 resno-1)
		     (let ((tmp resno-1))
		       (set! resno-1 resno-2)
		       (set! resno-2 tmp)))
		 (let loop ((rn resno-1))
		   (cond
		    ((> (+ rn 3) resno-2) 'done)
		    ((<= (+ rn 4) resno-2)
		     (add-extra-bond-restraint imol 
					       chain-id-1 rn       "" " O  " ""
					       chain-id-1 (+ rn 3) "" " N  " ""
					       3.18 0.035)
		     (add-extra-bond-restraint imol 
					       chain-id-1 rn       "" " O  " ""
					       chain-id-1 (+ rn 4) "" " N  " ""
					       2.91 0.035)
		     (loop (+ rn 1)))
		     (else 
		      (add-extra-bond-restraint imol 
						chain-id-1 rn       "" " O  " ""
						chain-id-1 (+ rn 3) "" " N  " ""
						3.18 0.035)
		      (loop  (+ rn 1)))))))))))



(define (user-defined-delete-restraint)
  (user-defined-click 2 
		      (lambda (atom-specs)
			(let* ((spec-1 (cddr (list-ref atom-specs 0)))
			       (spec-2 (cddr (list-ref atom-specs 1)))
			       (imol (cadr (list-ref atom-specs 0))))
			  (delete-extra-restraint imol (list 'bond spec-1 spec-2))))))

;; exte dist first chain A resi 19 ins . atom  N   second chain A resi 19 ins . atom  OG  value 2.70618 sigma 0.4
;;
(define (extra-restraints->refmac-restraints-file imol file-name)
  (let ((restraints (list-extra-restraints imol)))
    (if (list? restraints)
	(call-with-output-file file-name
	  (lambda (port)
	    (for-each (lambda (restraint)
			(if (eq? (car restraint) 'bond)
			    (begin
			      (let ((chain-id-1 (list-ref (cadr restraint) 1))
				    (resno-1    (list-ref (cadr restraint) 2))
				    (inscode-1  (list-ref (cadr restraint) 3))
				    (atom-1     (list-ref (cadr restraint) 4))
				    (chain-id-2 (list-ref (caddr restraint) 1))
				    (resno-2    (list-ref (caddr restraint) 2))
				    (inscode-2  (list-ref (caddr restraint) 3))
				    (atom-2     (list-ref (caddr restraint) 4))
				    (value      (list-ref restraint 3))
				    (esd        (list-ref restraint 4)))
				
			      (format port "EXTE DIST FIRST CHAIN ~a RESI ~s INS ~a ATOM ~a "
				      (if (or (string=? chain-id-1 " ")
					      (string=? chain-id-1 ""))
					  "."
					  chain-id-1)
				      resno-1
				      (if (or (string=? inscode-1 " ")
					      (string=? inscode-1 ""))
					  "."
					  inscode-1)
				      atom-1)
			      (format port " SECOND CHAIN ~a RESI ~s INS ~a ATOM ~a "
				      (if (or (string=? chain-id-2 " ")
					      (string=? chain-id-2 ""))
					  "."
					  chain-id-2)
				      resno-2 ;; fix Richard Bunker bug
				      (if (or (string=? inscode-2 " ")
					      (string=? inscode-2 ""))
					  "."
					  inscode-2)
				      atom-2)
			      (format port "VALUE ~s SIGMA ~s~%" value esd)))))
		      restraints))))))


(define (res-name->plane-atom-name-list res-name)

  (cond 
   ((not (string? res-name)) #f)
   ((string=? res-name  "DG")
    (list "N1" "C6" "O6" "C2" "N2" "N3" "C5" "C4" "N9" "N7" "C8" "C1'"))
   ((string=? res-name  "DA")
    (list "N1" "C6" "C2" "N6" "N3" "C5" "C4" "N9" "N7" "C8" "C1'"))
   ((string=? res-name  "DT")
    (list "N3" "C2" "O2" "C4" "O4" "C5" "C7" "C6" "N1" "C1'"))
   ((string=? res-name  "DC")
    (list "N3" "C2" "O2" "C4" "N4" "C5"      "C6" "N1" "C1'"))

   ((string=? res-name  "G")
    (list "N1" "C6" "O6" "C2" "N2" "N3" "C5" "C4" "N9" "N7" "C8" "C1'"))
   ((string=? res-name  "A")
    (list "N1" "C6" "C2" "N6" "N3" "C5" "C4" "N9" "N7" "C8" "C1'"))
   ((string=? res-name  "T")
    (list "N3" "C2" "O2" "C4" "O4" "C5" "C7" "C6" "N1" "C1'"))
   ((string=? res-name  "U")
    (list "N3" "C2" "O2" "C4" "O4" "C5"      "C6" "N1" "C1'"))
   ((string=? res-name  "C")
    (list "N3" "C2" "O2" "C4" "N4" "C5"      "C6" "N1" "C1'"))
   (else (list))))


;; Example of a stack restraint: "exte stac plan 1  firs resi 99 chai A atoms substitute-open-brace CB CG CD1 CD2 CE1 CE2 CZ OH substitute-close-brace  plan 2  firs resi 61 chai B atoms substitute-open-brace CB CG CD1 CD2 CE1 CE2 CZ substitute-close-brace  dist 3.4 sddi 0.2  sdan 6.0 type 1"
;; 
(define (write-refmac-parallel-plane-restraint file-name res-spec-0 res-spec-1 atom-list-0 atom-list-1)

  (call-with-output-file file-name
    (lambda (port)
      (display "EXTE STACK PLAN 1 FIRST RESIDUE " port)
      (display (residue-spec->res-no res-spec-0) port)
      (display " CHAIN " port)
      (display (residue-spec->chain-id res-spec-0) port)
      (display " ATOMS { " port)
      (for-each (lambda (atom-name)
		  (display " " port)
		  (display atom-name port)
		  (display " " port))
		atom-list-0)
      (display " } PLAN 2 FIRST RESIDUE " port)
      (display (residue-spec->res-no res-spec-1) port)
      (display " CHAIN " port)
      (display (residue-spec->chain-id res-spec-1) port)
      (display " ATOMS { " port)
      (for-each (lambda (atom-name)
		  (display " " port)
		  (display atom-name port)
		  (display " " port))
		atom-list-1)
      (display " } DIST 3.4 SDDI 0.2 SDAN 6.0 TYPE 1 " port))))


(define (add-parallel-planes-restraint imol rs-0 rs-1)

  (format #t "in add-parallel-planes-restraint: rs-0: ~s rs-1 ~s ~%" rs-0 rs-1)

  (let* ((rn-0 (residue-spec->residue-name imol rs-0))
	 (rn-1 (residue-spec->residue-name imol rs-1))
	 (atom-ls-0 (res-name->plane-atom-name-list rn-0))
	 (atom-ls-1 (res-name->plane-atom-name-list rn-1)))

  (write-refmac-parallel-plane-restraint "tmp.rst" rs-0 rs-1 atom-ls-0 atom-ls-1)
  (add-refmac-extra-restraints imol "tmp.rst")))



(if (defined? 'coot-main-menubar)

    (let* ((menu (coot-menubar-menu "Restraints")))
      
      (add-simple-coot-menu-menuitem 
       menu "Add Simple C-C Single Bond Restraint..."
       user-defined-add-single-bond-restraint)

      (add-simple-coot-menu-menuitem 
       menu "Add Distance Restraint..."
       user-defined-add-arbitrary-length-bond-restraint)

      (add-simple-coot-menu-menuitem 
       menu "Add Helix Restraints..."
       user-defined-add-helix-restraints)

      (add-simple-coot-menu-menuitem
       menu "RNA A form bond restraints..."
       user-defined-RNA-A-form)

      (add-simple-coot-menu-menuitem
       menu "DNA B form bond restraints..."
       user-defined-DNA-B-form)

      (add-simple-coot-menu-menuitem
       menu "Read Refmac Extra Restraints..."
       (lambda ()
         (generic-chooser-and-file-selector 
          "Apply restraints to molecule"
          valid-model-molecule?  "File:" ""
          (lambda (imol file-name)
            (add-refmac-extra-restraints imol file-name)))))

      (add-simple-coot-menu-menuitem
       menu "Delete Extra Restraints..."
       (lambda ()
	 (molecule-chooser-gui "Delete Extra Restraints for molecule:"
			       (lambda (imol) 
				 (delete-all-extra-restraints imol)))))

      (add-simple-coot-menu-menuitem 
       menu "Interesting limit to 0.5"
       (lambda ()
	 (using-active-atom
	  (set-extra-restraints-prosmart-sigma-limits aa-imol -0.5 0.5))))

      (add-simple-coot-menu-menuitem 
       menu "Interesting limit to 2.5"
       (lambda ()
	 (using-active-atom
	  (set-extra-restraints-prosmart-sigma-limits aa-imol -2.5 2.5))))

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
       menu "Extra Restraints to CA"
       (lambda ()
         (using-active-atom
          (set-extra-restraints-representation-for-bonds-go-to-CA aa-imol 1))))

      (add-simple-coot-menu-menuitem
       menu "Extra Restraints Standard Representation"
       (lambda ()
         (using-active-atom
          (set-extra-restraints-representation-for-bonds-go-to-CA aa-imol 0))))

      (add-simple-coot-menu-menuitem
       menu "Delete an Extra Restraint..."
       user-defined-delete-restraint)

      (add-simple-coot-menu-menuitem
       menu "Delete Restraints for this residue"
       (lambda()
	 (using-active-atom
	  (delete-extra-restraints-for-residue aa-imol aa-chain-id aa-res-no aa-ins-code))))

      (add-simple-coot-menu-menuitem
       menu "Delete Deviant Extra Restraints..."
       (lambda ()
	 (generic-single-entry "Delete Restraints worse than " "4.0" " Delete Outlying Restraints "
			       (lambda (text)
				 (let ((n (string->number text)))
				   (if (number? n)
				       (using-active-atom
					(delete-extra-restraints-worse-than aa-imol n))))))))

    (let ((menu (coot-menubar-menu "Restraints")))
      (add-simple-coot-menu-menuitem 
       menu "Add Parallel Planes restraint..."
       (lambda ()

	 (add-status-bar-text "Click on 2 residues to define the additional parallel planes restraint")
	 (user-defined-click 2 (lambda (atom-specs)

				 (let* ((atom-0 (list-ref atom-specs 0))
					(atom-1 (list-ref atom-specs 1))
					(rs-0 (atom-spec->residue-spec atom-0))
					(rs-1 (atom-spec->residue-spec atom-1))
					(imol (atom-spec->imol atom-0)))

				   (let ((rn-0 (residue-name imol
							     (residue-spec->chain-id (atom-spec->residue-spec atom-0))
							     (residue-spec->res-no   (atom-spec->residue-spec atom-0))
							     (residue-spec->ins-code (atom-spec->residue-spec atom-0))))
					 (rn-1 (residue-name imol
							     (residue-spec->chain-id (atom-spec->residue-spec atom-1))
							     (residue-spec->res-no   (atom-spec->residue-spec atom-1))
							     (residue-spec->ins-code (atom-spec->residue-spec atom-1)))))
					 
				     (format #t "got resname 0 ~s ~%" rn-0)
				     (format #t "got resname 1 ~s ~%" rn-1)

				     (let ((atom-ls-0 (res-name->plane-atom-name-list rn-0))
					   (atom-ls-1 (res-name->plane-atom-name-list rn-1)))

				       (write-refmac-parallel-plane-restraint "tmp.rst" rs-0 rs-1 atom-ls-0 atom-ls-1)
				       (add-refmac-extra-restraints imol "tmp.rst")

				       ))))))))))
