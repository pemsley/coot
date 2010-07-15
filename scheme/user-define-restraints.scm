

(define (res-name-from-atom-spec atom-spec)

  (let ((imol (cadr atom-spec)))
    (residue-name imol 
		  (list-ref atom-spec 2)
		  (list-ref atom-spec 3)
		  (list-ref atom-spec 4))))

  

(define (user-defined-add-single-bond)

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
							1.54 0.02))))))

(define (add-base-restraint imol spec-1 spec-2 atom-name-1 atom-name-2 dist)
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
			    dist 0.05))

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

(define (user-defined-RNA-A-form)
  (user-defined-click 2
		      (lambda (atom-specs)
			(let ((spec-1 (list-ref atom-specs 0))
			      (spec-2 (list-ref atom-specs 1)))
			  (let ((res-name-1 (res-name-from-atom-spec spec-1))
				(res-name-2 (res-name-from-atom-spec spec-2)))

			    (if (and (string=? res-name-1 "Gr")
				     (string=? res-name-2 "Cr"))
				(g-c-restraints spec-1 spec-2))

			    (if (and (string=? res-name-1 "Cr")
				     (string=? res-name-2 "Gr"))
				(g-c-restraints spec-2 spec-1))

			    (if (and (string=? res-name-1 "Ar")
				     (string=? res-name-2 "Ur"))
				(a-u-restraints spec-1 spec-2))

			    (if (and (string=? res-name-1 "Ur")
				     (string=? res-name-2 "Ar"))
				(a-u-restraints spec-2 spec-1)))))))

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
					       3.18 0.05)
		     (add-extra-bond-restraint imol 
					       chain-id-1 rn       "" " O  " ""
					       chain-id-1 (+ rn 4) "" " N  " ""
					       2.91 0.05)
		     (loop (+ rn 1)))
		     (else 
		      (add-extra-bond-restraint imol 
						chain-id-1 rn       "" " O  " ""
						chain-id-1 (+ rn 3) "" " N  " ""
						3.18 0.05)
		      (loop  (+ rn 1)))))))))))

(define (user-defined-delete-restraint)
  (user-defined-click 2 
		      (lambda (atom-specs)
			(let* ((spec-1 (cddr (list-ref atom-specs 0)))
			       (spec-2 (cddr (list-ref atom-specs 1)))
			       (imol (cadr (list-ref atom-specs 0))))
			  (delete-extra-restraint imol (list 'bond spec-1 spec-2))))))


(if (defined? 'coot-main-menubar)

    (let* ((menu (coot-menubar-menu "Extras")))
      
      (add-simple-coot-menu-menuitem 
       menu "Add Single Bond..."
       user-defined-add-single-bond)

      (add-simple-coot-menu-menuitem 
       menu "Add Helix Restraints..."
       user-defined-add-helix-restraints)

      (add-simple-coot-menu-menuitem
       menu "RNA A form bond restraints..."
       user-defined-RNA-A-form)

      (add-simple-coot-menu-menuitem
       menu "Delete an Extra Restraint..."
       user-defined-delete-restraint)))

;; Now make helix restraints, demonstrate with Bill Weiss, GPCR bent
;; helix

