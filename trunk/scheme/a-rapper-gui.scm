
(use-modules (gtk gdk)
             (gtk gtk)
             (gui event-loop)
             (gui entry-port)
             (gui text-output-port))


;; something that the user sets:
(define rapper-dir (append-dir-dir (getenv "HOME") "rappermc"))
; (define *rapper-command* (append-dir-file rapper-dir "rapper"))
(define *rapper-command* "rapper")

(define *arp/warp-loop-command* "unset")

(define rapper-process
  (let ((pid #f))
    (lambda args

      (if (= (length args 1))
	  (cond
	   ((eq? (car args) 'get)
	    pid)
	   ((eq? (car args) 'kill)
	    (if (number? pid)
		(begin
		  (run "kill" (number->string pid))
		  (set! pid #f))))
	   (else 
	    (format #t "ignored: ~s~%" args))))

      (if (= (length args 2))
	  (cond 
	   ((eq? (car args) 'store)
	    (set! pid (cdr args)))
	   (else 
	    (format #t "ignored: ~s~%" args)))))))


(define (rename-dir-by-date dir-name)
  (let* ((date-str (localtime (current-time)))
	 (new-dir-name (string-append
			dir-name "-"  (strftime "%Y-%h-%d--%H-%M--%S" date-str))))
    (rename-file dir-name new-dir-name)
    new-dir-name))

		 
;; return a string or #f
;;
(define (sequence-string imol chain-id resno-start resno-end)

  (define (all-chars? ls)
    (cond
     ((null? ls) #t)
     ((not (char? (car ls)))
      #f)
     (else
      (all-chars? (cdr ls)))))

  (if (not (valid-model-molecule? imol))
      #f
      (let ((single-letter-code-list
	     (map (lambda (resno)
		    (let ((res-name (residue-name imol chain-id resno "")))
		      (3-letter-code->single-letter res-name)))
		  (range resno-start (+ resno-end 1)))))
	(if (not (all-chars? single-letter-code-list))
	    (begin
	      (format #t "bad sequence chars ~s~%" single-letter-code-list)
	      #f)
	    (list->string single-letter-code-list)))))





(define (rapper-it imol chain-id start-resno end-resno sequence number-of-models)

  (let ((imol-map (imol-refinement-map)))
    (if (not (valid-map-molecule? imol-map))
	(format #t "No valid map molecule given (possibly ambiguous)~%")
	(let* ((str (string-append "//" chain-id "/" (number->string start-resno)
				   "-" (number->string end-resno)))
	       (frag-mol (new-molecule-by-atom-selection imol str))
	       (fragment-pdb "coot-rapper-fragment-in.pdb")
	       (rapper-out-pdb "rapper-out.pdb")
	       (rapper-mode "model-loops-benchmark") ; maybe "ca-trace" perhaps.
	       ;; (rapper-mode "ca-trace")
	       (length-for-rapper 3) ; test
	       (sequence-string (if (string? sequence)
				    sequence
				    "AAA"))
	       (map-file "coot-rapper.map")
	       (whole-pdb-file-name "rapper-all-atoms.pdb"))
	  (write-pdb-file frag-mol fragment-pdb)
	  (write-pdb-file imol whole-pdb-file-name)
	  (close-molecule frag-mol)
	  (export-map imol-map map-file)
	  (format #t "running rapper: ~s ~s ~s ~s ~s ~s~%"
		  imol chain-id start-resno end-resno sequence number-of-models)
;	  (let ((rapper-pid (run-concurrently "rapper" "params.xml"
;					      "--pdb" fragment-pdb 
;					      "--map" map-file 
;					      "--start" (number->string start-resno)
;					      "--stop"  (number->string   end-resno)
;					      "--mainchain-restraint-threshold" "2.0"
;					      "--sidechain-centroid-restraint-threshold"  "2.0"
;					      "--sidechain-mode" "smart"
;					      "--sidechain-radius-reduction" "0.75"
;					      "--enforce-mainchain-restraints" "true"
;					      "--enforce-sidechain-centroid-restraints" "true"
;					      "--edm-fit" "true"
;					      "--rapper-dir" rapper-dir)))
	    ; (rapper-process 'store rapper-pid))))))
	  (if (file-exists? "TESTRUNS")
	      (rename-dir-by-date "TESTRUNS"))
	  (let* ((rapper-command-line-args
		  (list ; "params.xml"
		   rapper-mode
		   "--pdb" fragment-pdb
		   ;; "--framework" whole-pdb-file-name
		   "--pdb-out" rapper-out-pdb
		   "--map" map-file 
		   "--models" (number->string number-of-models)
		   "--start" (number->string start-resno)
		   "--stop"  (number->string   end-resno)
		   ;; "--length" (number->string length-for-rapper)
		   "--seq" sequence-string
		   "--mainchain-restraint-threshold" "2.0"
		   "--sidechain-centroid-restraint-threshold"  "2.0"
		   "--sidechain-mode" "smart"
		   "--sidechain-radius-reduction" "0.75"
		   "--enforce-mainchain-restraints" "true"
		   "--enforce-sidechain-centroid-restraints" "true"
		   "--edm-fit" "true"
		   "--rapper-dir" rapper-dir))
		 (rapper-status (goosh-command *rapper-command*
					       rapper-command-line-args
					      '()
					      "rapper.log"
					      #t)))
	    (format #t "rapper-status: ~s~%" rapper-status)
	    (if (file-exists? "TESTRUNS")
		(let* ((new-dir-name (rename-dir-by-date "TESTRUNS"))
		       (result-pdb-file-name (append-dir-file new-dir-name "looptest-best.pdb")))
		  (if (file-exists? result-pdb-file-name)
		      (read-pdb result-pdb-file-name)
		      (begin
			(info-dialog "RAPPER failed - no results"))))))))))

    

(define (stop-rapper)
  (format #t "stopping rapper process...~%"))


(define (cancel-dialog-func window)
  
  ;; do something more clever if the rapper process is still running.
  (gtk-widget-destroy window))

;; loop-building-tool is either 'rapper or 'ARP/wARP 
(define (a-rapper-gui loop-building-tool)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 2))
	 (pdb-hbox (gtk-hbox-new #f 2))
	 (chain-hbox (gtk-hbox-new #f 2))
	 (start-resno-hbox (gtk-hbox-new #f 2))
	 (end-resno-hbox (gtk-hbox-new #f 2))
	 (sequence-vbox (gtk-vbox-new #f 2))
	 (models-hbox (gtk-hbox-new #f 2))
	 (buttons-hbox (gtk-hbox-new #f 2))
	 (h-sep (gtk-hseparator-new))
	 (stop-button (gtk-button-new-with-label "  Stop  "))
	 (go-button (gtk-button-new-with-label "  Go  "))
	 (cancel-button (gtk-button-new-with-label "  Cancel  "))
	 
	 ; pdb
	 (label-pdb (gtk-label-new "Rebuild fragment of Model: "))
	 (option-menu-pdb (gtk-option-menu-new))
	 (menu-pdb (gtk-menu-new))

	 ; chain 
	 (label-chain (gtk-label-new "Chain"))
	 (entry-chain (gtk-entry-new))

	 ; start resno
	 (label-start-resno (gtk-label-new "Starting Res Number: "))
	 (entry-start-resno (gtk-entry-new))
	 
	 ; end resno
	 (label-end-resno (gtk-label-new "        End Res Number: "))
	 (entry-end-resno (gtk-entry-new))
	 
	 ; :: sequence

	 ; sequence as it currently is:
	 (sequence-as-is-check-button 
	  (gtk-check-button-new-with-label "As is"))

	 ; sequence as text
	 (label-sequence (gtk-label-new "Sequence"))
	 (text-sequence (gtk-text-new #f #f))

	 ; number of rapper models
	 (label-models (gtk-label-new "Number of Models: "))
	 (entry-models (gtk-entry-new))

	 ; 
	 (model-mol-list (fill-option-menu-with-coordinates-mol-options 
			  menu-pdb)))

    (gtk-text-set-editable text-sequence #t)

    ;; use instead gtk-widget-set-size-request? (not supported on
    ;; penelope yet.
    ;; 
    (gtk-widget-set-usize text-sequence -1 80)
    (gtk-option-menu-set-menu option-menu-pdb menu-pdb)
    (gtk-box-pack-start pdb-hbox label-pdb #f #f 2)
    (gtk-box-pack-start pdb-hbox option-menu-pdb #f #f 2)
    (gtk-box-pack-start chain-hbox label-chain #f #f 2)
    (gtk-box-pack-start chain-hbox entry-chain #f #f 2)
    (gtk-box-pack-start start-resno-hbox label-start-resno #f #f 2)
    (gtk-box-pack-start start-resno-hbox entry-start-resno #f #f 2)
    (gtk-box-pack-start end-resno-hbox   label-end-resno   #f #f 2)
    (gtk-box-pack-start end-resno-hbox   entry-end-resno   #f #f 2)

    (gtk-box-pack-start sequence-vbox label-sequence #f #f 2)
    (gtk-box-pack-start sequence-vbox sequence-as-is-check-button #f #f 2)
    (gtk-box-pack-start sequence-vbox text-sequence  #f #f 2)

    (gtk-box-pack-start models-hbox label-models)
    (gtk-box-pack-start models-hbox entry-models)
    
    (gtk-window-set-title window "A Rapper GUI")
    (gtk-box-pack-start vbox pdb-hbox  #f #f 2)
    (gtk-box-pack-start vbox chain-hbox  #f #f 2)
    (gtk-box-pack-start vbox start-resno-hbox  #f #f 2)
    (gtk-box-pack-start vbox   end-resno-hbox  #f #f 2)
    (gtk-box-pack-start vbox  sequence-vbox  #f #f 2)
    (gtk-box-pack-start vbox  models-hbox  #f #f 2)
    (gtk-box-pack-start vbox  h-sep  #f #f 2)

    (gtk-box-pack-start buttons-hbox go-button   #t #f 6)
    (gtk-box-pack-start buttons-hbox stop-button #t #f 6)
    (gtk-box-pack-start buttons-hbox cancel-button #t #f 6)

    (gtk-box-pack-start vbox  buttons-hbox  #t #f 2)
    (gtk-container-add window vbox)

    (gtk-entry-set-text entry-models "2")

    (gtk-toggle-button-set-active sequence-as-is-check-button #t)

    (gtk-signal-connect stop-button "clicked" (lambda () (stop-rapper)))

    (gtk-signal-connect cancel-button "clicked" (lambda () (cancel-dialog-func window)))

    (gtk-signal-connect go-button "clicked"
			(lambda ()
			  (let ((imol (get-option-menu-active-molecule
				       option-menu-pdb
				       model-mol-list)))
			    (if (not (number? imol))
				(format #t "bad active model~%")
				(let* ((chain-text (gtk-entry-get-text entry-chain))
				       (start-resno-text (gtk-entry-get-text entry-start-resno))
				       (end-resno-text (gtk-entry-get-text entry-end-resno))
				       (nov (format #t "gtk-text? for ~s returns ~s~%"
						    text-sequence (gtk-text? text-sequence)))
				       ; (text-buffer (gtk-text-get-buffer text-sequence))
				       (l (gtk-text-get-length text-sequence))
				       (text-sequence-text (gtk-editable-get-chars text-sequence 0 l))
				       (sequence "")
				       (number-of-models-text (gtk-entry-get-text entry-models))
				       (start-resno (string->number start-resno-text))
				       (end-resno (string->number end-resno-text))
				       (number-of-models (string->number number-of-models-text)))

				  (format #t "text-sequence-text: ~s~%" text-sequence-text)
				  
				  (if (not (and (number? start-resno)
						(number? end-resno)
						(number? number-of-models)))
				      (format #t "Something incomprehensible: ~s ~s ~s"
					      start-resno end-resno number-of-models)
				      (begin
					(if (eq? loop-building-tool 'rapper)
					    (let ((seq
						   (sequence-string 
						    imol chain-text start-resno end-resno)))
					      (rapper-it imol chain-text start-resno end-resno
							 seq number-of-models)))
					(if (eq? loop-building-tool 'ARP/wARP)
					    (let* ((new-start (- start-resno 1))
						   (new-end   (+ end-resno 1))
						   (seq (sequence-string imol chain-text
									 new-start new-end)))
					      (arp/warp-it imol chain-text 
							   (cons start-resno end-resno)
							   (cons new-start new-end)
							   seq number-of-models)))))))
			    (gtk-widget-destroy window))))

    (gtk-widget-show-all window)))


(let ((menu (coot-menubar-menu "Loop")))

  (add-simple-coot-menu-menuitem 
   menu "RAPPER..."
   (lambda()
     (a-rapper-gui 'rapper)))

  (add-simple-coot-menu-menuitem 
   menu "ARP/wARP Loopy..."
   (lambda()
     (a-rapper-gui 'ARP/wARP))))


; I have no idea why I wanted a timer here...
; 
; Oh... maybe to read in results as they are generated?  I need rapper
; to test this script now.
;
;(let ((timeout-handle (gtk-timeout-add 
;		       1000 ; ms
;		       (lambda () (display "timed out\n") 0))))
;  (format #t "timeout-handle: ~s~%" timeout-handle))
;(gtk-main)

