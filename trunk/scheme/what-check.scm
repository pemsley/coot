
;; "  53" -> "53", " " -> ""
(define (strip-leading-spaces str)
  (let ((l (string-length str)))
    (if (= l 0)
	str
	(if (not (string=? (substring str 0 1) " "))
	    str
	    (strip-leading-spaces (substring str 1))))))


;; consider for coot-utils
(define (go-to-residue-by-spec imol spec)

  (define (atom-spec? spec)
    (= (length spec) 5))
  
  ;; main line
  (set-go-to-atom-molecule imol)
  (let ((atom-name (if (atom-spec? spec)
		       (list-ref spec 3) " CA ")))

;    (format #t "debug::=== go to atom ~s ~s ~s atom-spec? ~s from spec ~s ~%"
;	    (car spec) (cadr spec) atom-name (atom-spec? spec) spec)
    (set-go-to-atom-chain-residue-atom-name (car spec)
					    (cadr spec)
					    atom-name)))


;; consider for coot-utils
(define (residue-spec->string spec)
  (string-append (list-ref spec 0)
		 (number->string (list-ref spec 1))
		 (list-ref spec 2)))

;; make this internal 
;; 
;; where line is something like  "  53 HIS  (  53 ) A12"
;;                                   ignore   xxxxI Cyy
;; xxxx is 4char resno, I insertion code, C chain id, yy is model number
;; 
;; return a residue spec, or #f
;; 
(define (get-flip-residue line)
  (let ((l (string-length line)))
    (if (< l 19)
	#f
	(let* ((resno-string (substring line 11 15))
	       (resno (string->number (strip-leading-spaces resno-string)))
	       (inscode-pre (substring line 15 16))
	       (inscode (if (string=? inscode-pre " ") "" inscode-pre))
	       (chain-id (substring line 18 19)))
	  (format #t "found ~s ~s ~s~%" resno inscode chain-id)
	  (if (not (number? resno))
	      #f
	      (list chain-id resno inscode))))))

(define convert-to-button-info
  (lambda (problemed-res-list)
    (format #t "convert-to-button-info got ~s~%" problemed-res-list)
    (let ((button-label      (list-ref problemed-res-list 1))
	  (problem-string    (list-ref problemed-res-list 0))
	  (action-func       (list-ref problemed-res-list 2))
	  (residue-spec-list (list-ref problemed-res-list 3)))

      (let ((o (map (lambda (x)
		      (list problem-string button-label action-func x))
		    residue-spec-list)))
	(format #t "convert-to-button-info returns ~s~%" o)
	o))))


;; A problemed-res-list is improper pair where the car is a pair of
;; string, one describing the problem, and the other being the "fix
;; it" button label.  The cdr is a list of residue specs.
;; 
;; A problemed-flip-list-list is a list of those things.
;; 
(define (problem-residues->dialog imol problemed-res-list-list)

  (if (not (null? problemed-res-list-list))
      (dialog-box-of-pairs-of-buttons imol
       "What-Check Report" (cons 220 250)
       ;; buttons is a list of: (list (list button-1-label button-1-action
       ;;                                   button-2-label button-2-action))
       ;; The button-1-action function takes as an argument the imol.
       ;; The button-2-action function takes as an argument the imol.
       (map 
	(lambda (flip-res-info)
	  (let ((flip-res-spec  (list-ref flip-res-info 3))
		(button-label   (list-ref flip-res-info 1))
		(action-func    (list-ref flip-res-info 2))
		(problem-string (list-ref flip-res-info 0)))
	    (list 
	     (string-append problem-string "  " (residue-spec->string flip-res-spec))
	     (lambda (imol-local)
	       (go-to-residue-by-spec imol-local flip-res-spec))
	     button-label
	     (lambda (imol-local)
	       (action-func imol-local flip-res-spec)))))
	
	;; Convert a list of: ("Needs Flip" "Flip" '( '("A" 54 "")
	;;                                            '("B" 12 "")))
	;;        to list of: '("Needs Flip" "Flip" '("A" 54 ""))
	;;                    '("Needs Flip" "Flip" '("B" 12 ""))
	;; 
	(apply append (map convert-to-button-info problemed-res-list-list)))
       "  Close  ")))



;; action is either 'gui (then we get a gui) or 'apply-actions, then the model
;; modifications are automatically applied.
;; 
(define (parse-check-db imol file-name action)

  ;; 
  (define my-flip
    (lambda (imol-local res-spec)
      (do-180-degree-side-chain-flip imol-local
				     (car res-spec)
				     (cadr res-spec)
				     (caddr res-spec)
				     "")))

  (define my-rsr
    (lambda (imol-local res-spec)
      (with-auto-accept
       (refine-zone imol-local (car res-spec) (cadr res-spec) (cadr res-spec) ""))))

  ;; 
  (define my-delete-atom
    (lambda (imol atom-spec)
      (let ((ls (cons imol atom-spec)))
	(format #t "calling delete-atom with args: ~s~%" ls)
	(apply delete-atom ls))))

  ;;
  (define my-rotamer-search
    (lambda (imol res-spec) 

      (let ((chain-id (list-ref res-spec 0))
	    (resno    (list-ref res-spec 1))
	    (ins-code (list-ref res-spec 2))
	    (altloc   ""))

	(let ((imol-map (imol-refinement-map)))
	  (if (valid-map-molecule? imol-map)
	      (begin
		(format #t "rotamer search with imol: ~s chain: ~s resno: ~s inscode: ~s~% " 
			imol chain-id resno ins-code)
		(auto-fit-best-rotamer resno altloc ins-code chain-id imol imol-map 1 0.01)))))))


  ;; return #f or a residue spec:
  (define (line->residue-spec line) 
    (if (> (string-length line) 21)
	(let* ((resno-string (substring line 26 30))
	       (resno (string->number (strip-leading-spaces resno-string)))
	       (inscode-pre (substring line 40 41))
	       (inscode (if (string=? inscode-pre "_") "" inscode-pre))
	       (chain-id-pre (substring line 19 20))
	       (chain-id (if (string=? chain-id-pre "_")
			     " " 
			     chain-id-pre)))
;	  (format #t "found ~s (from ~s) inscode ~s ~s chain-id ~s in ~s ~%" 
;		  resno resno-string inscode-pre inscode chain-id line)
	  (if (not (number? resno))
	      #f
	      (list chain-id resno inscode)))))

  ;; Return #f or an atom spec.
  ;; an atom-spec is (list chain-id resno ins-code atom-name alt-conf)
  ;; 
  (define (line->atom-spec line) 
    ;; (format #t "===== line->atom-spec ====== ~s ~s ~%" line (string-length line))
    (if (> (string-length line) 53)
	(let* ((resno-string (substring line 26 30))
	       (resno (string->number (strip-leading-spaces resno-string)))
	       (inscode-pre (substring line 40 41))
	       (inscode (if (string=? inscode-pre "_") "" inscode-pre))
	       (chain-id-pre (substring line 19 20))
	       (chain-id (if (string=? chain-id-pre "_")
			     ""
			     chain-id-pre))
	       (atom-name-pre (substring line 47 51))
	       (atom-name atom-name-pre)
	       (alt-conf ""))
	  (format #t "found atom name ~s resno ~s (from ~s) inscode ~s (~s) chain-id ~s ~%" 
		  atom-name resno resno-string inscode-pre inscode chain-id)
	  (format #t "line: ~s~%" line)
	  (if (not (number? resno))
	      #f
	      (list chain-id resno inscode atom-name alt-conf)))))

  (define (apply-actions imol action-list)
    
    (for-each 
     (lambda (action)

       (let* ((func (car action))
	      (baddie-list (car (cdr action)))
	      (imol-list (map (lambda(x) imol) baddie-list)))
	 (format #t "applying ~s to ~s~%" func baddie-list)
	 (map func imol-list baddie-list)))
     
     action-list))

  ;; main line
  (format #t "parse-whatcheck-report: parsing output of what-if...~%")
  (let ((found-flip-table #f)
	(found-h2ohbo-table #f)
	(h2o-h-bond-less-list '())
	(found-BMPCHK-table #f)
	(found-NAMCHK-table #f)
	(found-PLNCHK-table #f)
	(plane-list '())
	(flip-list '())
	(namchk-list '())
	(bump-list '()))

    (if (not (file-exists? file-name))
	(format #t "file not found ~s~%" file-name)
	
	(call-with-input-file file-name
	  (lambda (port)

	    (let loop ((line (read-line port)))

	      ;; (format #t "checking: ~s~%" line)

	      (cond 
	       ((eof-object? line) '()) ; done

	       ((string=? "//" line)
		(set! found-flip-table   #f)
		(set! found-h2ohbo-table #f)
		(set! found-BMPCHK-table #f)
		(set! found-NAMCHK-table #f)
		(set! found-PLNCHK-table #f)
		(loop (read-line port)))

	       ((string=? "CheckID   : HNQCHK" line)
		(set! found-flip-table #t)
		; (format #t "found flip-table from ~s~%" line)
		(loop (read-line port)))
	       ((eq? found-flip-table #t)
		(if (> (string-length line) 12)
		    (let ((attr-type (substring line 0 11)))
		      (if (string=? attr-type "   Name   :")
			  (begin
			    ; (format #t "found flip name ~s~%" line)
			    (let ((res (line->residue-spec line)))
			      (if res
				  (set! flip-list (cons res flip-list))))))))
		(loop (read-line port)))

	       ((string=? "CheckID   : H2OHBO" line)
		(set! found-h2ohbo-table #t)
		(loop (read-line port)))
	       ((eq? found-h2ohbo-table #t)
		(if (> (string-length line) 12)
		    (let ((attr-type (substring line 0 11)))
		      (if (string=? attr-type "   Name   :")
			  (begin
			    (format #t "found h2obmo name ~s~%" line)
			    (let ((atom (line->atom-spec line)))
			      (if atom
				  (set! h2o-h-bond-less-list 
					(cons atom h2o-h-bond-less-list))))))))
		(loop (read-line port)))

	       ((string=? "CheckID   : BMPCHK" line)
		(set! found-BMPCHK-table #t)
		(loop (read-line port)))
	       ((eq? found-BMPCHK-table #t)
		(if (> (string-length line) 12)
		    (let ((attr-type (substring line 0 11)))
		      (if (string=? attr-type "   Name   :")
			  (begin
			    (let ((res (line->residue-spec line)))
			      (if res
				  (set! bump-list (cons res bump-list))))))))
		(loop (read-line port)))

	       ((string=? "CheckID   : NAMCHK" line)
		(set! found-NAMCHK-table #t)
		(loop (read-line port)))
	       ((eq? found-NAMCHK-table #t)
		(if (> (string-length line) 12)
		    (let ((attr-type (substring line 0 11)))
		      (if (string=? attr-type "   Name   :")
			  (begin
			    (let ((res (line->residue-spec line)))
			      (if res
				  (set! namchk-list (cons res namchk-list))))))))
		(loop (read-line port)))

	       ((string=? "CheckID   : PLNCHK" line)
		(set! found-PLNCHK-table #t)
		(loop (read-line port)))
	       ((eq? found-PLNCHK-table #t)
		(if (> (string-length line) 12)
		    (let ((attr-type (substring line 0 11)))
		      (if (string=? attr-type "   Name   :")
			  (begin
			    (let ((res (line->residue-spec line)))
			      (if res
				  (set! plane-list (cons res plane-list))))))))
		(loop (read-line port)))

		
	       (else
		(loop (read-line port)))))

	    ;; (format #t "Here flip-list is ~s~%" flip-list)
	    ;; (format #t "Here h2o-h-bond-less-list is ~s~%" h2o-h-bond-less-list)
	    ;; (format #t "Here bump-list is ~s~%" bump-list)

	    (if (eq? action 'gui)
		(problem-residues->dialog 
		 imol (list 
		       (list "Needs flipping" "Flip it" my-flip flip-list)
		       (list "Water with no H-bonds" "Delete it" my-delete-atom h2o-h-bond-less-list)
		       (list "Has Bumps" "Rotamer Search" my-rotamer-search bump-list)
		       (list "Nomenclature error" "180 degree flip" my-flip namchk-list)
		       )))

	    (if (eq? action 'apply-actions)
		(apply-actions imol 
			       (list
				;; (list my-flip flip-list)
				;; (list my-delete-atom h2o-h-bond-less-list)
				;; (list my-rotamer-search bump-list)
				(list my-rsr plane-list)
				;; (list my-flip namchk-list)))))))))
				))))))))

