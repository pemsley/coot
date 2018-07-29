
(define (get-glyco-model-dir)
  (append-dir-dir (append-dir-dir (pkgdatadir) "data") "cho-models"))

(define (get-glyco-tree-acedrg-pyranose-dictionaries-dir) 
  (append-dir-dir (append-dir-dir (pkgdatadir) "data") "cho-acedrg"))


(define (remove-duplicates ls)
   (cond
      ((null? ls) '())
      (else
         (cons (car ls) (remove-duplicates (filter (lambda (item) (not (equal? item (car ls)))) (cdr ls)))))))
;; (let ((ls (list 1 2 3 3 4 4))) (let ((result (remove-duplicates ls))) (format #t "result: ~s~%" result)))

(define (read-acedrg-pyranose-dictionaries)
   (let ((cho-list (list "NAG" "MAN" "BMA" "SIA" "GLC" "FUC" "XYP")))
      (for-each (lambda (comp-id)
          (let ((fn (append-dir-file (get-glyco-tree-acedrg-pyranose-dictionaries-dir) (string-append comp-id "-acedrg.cif"))))
	      (if (file-exists? fn)
                 (read-cif-dictionary fn))))
          cho-list)))

(define (add-cho-restraints-for-residue imol residue-spec)
   (if (list? residue-spec)
       (let ((id (glyco-tree-residue-id imol residue-spec)))
	 ;; (format #t "debug:: --------------------- add-cho-restraints-for-residue: glyco-tree-residue-id: ~s~%" id)
	 ;; (format #t "debug:: --------------------- add-cho-restraints-for-residue: glyco-tree-residues: ~s~%" (glyco-tree-residues-scm imol residue-spec))
	 (add-cho-restraints-for-residue-with-id imol residue-spec id))))

;; the smaller the sigmas the more weight the restraints have
;; This pushes up the weight a bit (c.f. 1.0).
;;
(define *cho-geman-mcclure-sigma-scale* 0.5)

(define (add-cho-restraints-for-residue-with-id imol residue-spec glyco-id)
 
   ;; hacky
   (define (pad-name str)
      (let ((l (string-length str)))
        (if (> l 3)
           str
           (if (= l 3)
              (string-append " " str)
              (if (= l 2)
                 (string-append " " str " ")
                 (string-append " " str "  "))))))

   ;; list (atom-1 atom-2 d esd)
   ;; or #f
   (define (line->extra-bond-restraint-spec parent-residue-spec line)
       (let ((parts (string->list-of-strings line)))
           (if (= (length parts) 7)
               (let ((at-name-1 (pad-name (list-ref parts 0)))
                     (at-name-2 (pad-name (list-ref parts 1)))
                     (mean-str  (list-ref parts 2))
                     (sd-str    (list-ref parts 3))
                     (n-str     (list-ref parts 4))
                     (d-str     (list-ref parts 5))
		     (mod-sarle-str (list-ref parts 6)))
		 (let ((mean (string->number mean-str))
		       (sd   (string->number sd-str))
		       (n    (string->number n-str))
		       (d    (string->number d-str))
		       (mod-sarle (string->number mod-sarle-str)))

		   (if (not (>= n 20))
		       #f
		       (if (not (< d 0.42))
			   #f
			   (list (list (residue-spec->chain-id residue-spec)
				       (residue-spec->res-no   residue-spec)
				       (residue-spec->ins-code residue-spec)
				       at-name-1 "")
				 (list (residue-spec->chain-id parent-residue-spec)
				       (residue-spec->res-no   parent-residue-spec)
				       (residue-spec->ins-code parent-residue-spec)
				       at-name-2 "")
				 mean (* sd *cho-geman-mcclure-sigma-scale*)))))))))

  (define (glyco-id->level-number glyco-id)
    (list-ref glyco-id 0))

  (define (glyco-id->prime-arm-sym glyco-id)
    (list-ref glyco-id 1))

  (define (glyco-id->residue-type glyco-id)
    (list-ref glyco-id 2))

  (define (glyco-id->link-type glyco-id)
    (list-ref glyco-id 3))

  (define (glyco-id->parent-residue-type glyco-id)
    (list-ref glyco-id 4))

  (define (glyco-id->residue-spec glyco-id)
    (list-ref glyco-id 5))

   ;; main line
   ;;
   (if (list? glyco-id)
       (let* ((dir (get-glyco-model-dir))
	      (level               (glyco-id->level-number        glyco-id))
	      (prime-arm-sym       (glyco-id->prime-arm-sym       glyco-id))
	      (res-name            (glyco-id->residue-type        glyco-id))
	      (link-type           (glyco-id->link-type           glyco-id))
	      (parent-res-name     (glyco-id->parent-residue-type glyco-id))
	      (parent-res-spec     (glyco-id->residue-spec        glyco-id))
	      (fn-file (string-append "model-level-" (number->string level) "-" res-name "-" link-type "-" parent-res-name ".tab")))
         (let ((model-fn (append-dir-file dir fn-file)))
	   (if (not (file-exists? model-fn))
	      (begin
                (format #t "model file does not exist ~s~%~!" model-fn))
              (let ((lines (reverse
	            (call-with-input-file model-fn
                       (lambda (port)
                          (let loop ((lines '()) (line (read-line port)))
	                     (cond
                                ((eof-object? line) lines)
                                (else
                                   (loop (cons line lines) (read-line port))))))))))
                 (format #t "INFO:: read ~s lines from file ~s~%" (length lines) model-fn)
                 (let ((new-restraints (filter list? (map (lambda(line) (line->extra-bond-restraint-spec parent-res-spec line)) lines))))
		    ;; (for-each (lambda (nr) (format #t "debug ~s~%" nr)) new-restraints)
		    (add-extra-bond-restraints-scm imol new-restraints))
   ))))))


(define (test-get-cho-restraints imol)
   (let ((raw-carbo-tree-list '()))

      (for-each (lambda (chain-id)
         (for-each (lambda (res-serial)
            (let ((res-no (seqnum-from-serial-number imol chain-id res-serial)))
               (let ((rn (residue-name imol chain-id res-no "")))
                  (if (string? rn)
                     (if (string=? rn "NAG")
                        (let* ((residue-spec (list chain-id res-no ""))
                               (rl (glyco-tree-residues-scm imol residue-spec)))
                           (format #t "rl: ~s~%" rl)
                           (print-glyco-tree imol chain-id res-no "")
                           (if (not (list? rl))
                              (begin
                                 (format #t "bad glyco-tree for residue ~s~%~!" residue-spec))
                              (if (> (length rl) 3)
                                 (let ((id (glyco-tree-residue-id imol residue-spec)))
                                    (format #t "### residue ~s residue id: ~s~%" residue-spec id)
                                    (add-cho-restraints-for-residue-with-id imol residue-spec id))))))))))
            (range (chain-n-residues chain-id imol))))
         (chain-ids imol))))

(define (correlation-coefficient-of-this-tree)
  (using-active-atom
   (let ((residues (glyco-tree-residues aa-imol aa-res-spec)))
     (let ((cc (map-to-model-correlation aa-imol residues '() 0 (imol-refinement-map))))
       (format #t "residues ~s\ncc: ~5,3f~%" residues cc)
       (info-dialog (format #f "residues ~s~%cc ~5,3f" residues cc))))))


