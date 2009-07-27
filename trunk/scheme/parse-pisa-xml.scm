
(use-modules (sxml simple))


(define (pisa-xml imol file-name)
  (if (file-exists? file-name)
      (begin
	(format #t "opened ~s~%" file-name)
	(call-with-input-file file-name
	  (lambda (port)
	    (let ((sm (xml->sxml port)))
	      (format #t "sm: ~s~%" sm)
	      (parse-pisa imol sm)))))))


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


  
(define parse-pisa
  (lambda (imol entity)
    
    (define (handle-molecule molecule)
      (let ((symbols (list 'rxx 'rxy 'rxz 'ryx 'ryy 'ryz 'rzx 'rzy 'rzz 'tx 'ty 'tz)))

	(let ((ass-symbols '()))
	  (for-each
	   (lambda (mol-ele)
	     ;; (format #t "mol-ele: ~s~%" mol-ele)
	     (if (list? mol-ele) 
		 (for-each 
		  (lambda (symbol)
		    (if (eq? (car mol-ele) symbol)
			(let ((n (string->number (cadr mol-ele))))
			  ;; (format #t "~s: ~s~%" symbol n)
			  (set! ass-symbols 
				(cons (list symbol n)
				      ass-symbols)))))
		  symbols)))
	   molecule)
	
	  (if (= (length ass-symbols) 12)
	      (let ((mat 
		     (map (lambda (sym)
			    (cadr (assoc sym ass-symbols)))
			  symbols)))
		(format #t "mat:::: ~s~%" mat)
		(format "currently ~s molecules~%" (graphics-n-molecules))
		(let ((new-mol-no (apply new-molecule-by-symmetry imol
					 ""
					 (append mat (list 0 0 0)))))
		  (format #t "created molecule number ~s~%" new-mol-no)))))))


    (define (handle-assembly assembly)
      
      (for-each
       (lambda (ele)
	 ;; (format #t "ele: ~s~%" ele)
	 (if (list? ele) 
	     (if (eq? (car ele) 'molecule)
		 (begin
		   (format #t "molecule: ~s~%" ele)
		   (handle-molecule ele)))))
       assembly))

    
    ;; main line
    (let loop ((entity entity))
      (cond
       ((list? entity) (if (not (null? entity))
			   (begin
			     (if (eq? 'assembly (car entity))
				 (begin 
				   (handle-assembly entity))
				 (map loop entity)))))
       (else 
	#f)))))
     


;;(pisa-xml "pisa.xml")

(define (t)
  (pisa-xml (read-pdb "coot-download/3FCS.pdb") "pisa.xml"))



(let ((menu (coot-menubar-menu "PISA")))

    (add-simple-coot-menu-menuitem
     menu "PISA..."
     (lambda ()
       (pisa-xml (read-pdb "coot-download/pdb3lz2.ent") "pisa.xml")))

    (add-simple-coot-menu-menuitem
     menu "entry code..."
     (lambda ()
       (generic-chooser-and-entry "PISA xml for" " XML file-name:" "" 
			     (lambda (imol text)
			       (pisa-xml imol text))))))


