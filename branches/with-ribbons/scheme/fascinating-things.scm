; -*-scheme-*-

;;(molprobity-fascinating-clusters-things-gui
;;   dialog-name
;;   sorting-options
;;   list-of-clusters)
;; 
;; where a cluster is:
;;    (list 
;;     cluster-name-string
;;     cluster-center-go-button-label-string
;;     ccgb-x ccgb-y ccgb-z
;;     ; a list of specific items:
;;     (list 
;;         (list specific-button-label-string button-red button-green button-blue
;;               specific-x specific-y specific-z))))
;;
;; 
;;(molprobity-fascinating-clusters-things-gui
;; gui-name-string
;; (list
;;   (list "Active Site" (list 0 1 2 3 4))
;;   (list "Worst First" (list 3 4 2 1 0)))
;; ; now a list of clusters:
;; (list
;;    (list cluster-name-string
;;       cluster-center-go-button-label-string
;;       ccgb-x ccgb-y ccgb-z
;;       ; now a list of specific items
;;       (list
;;         (list specific-button-label-string button-red button-green button-blue
;;               specific-x specific-y specific-z)
;;         (list specific-button-label-string button-red button-green button-blue
;;               specific-x specific-y specific-z)))

;;    (list cluster-name-string
;;       cluster-center-go-button-label-string
;;       ccgb-x ccgb-y ccgb-z
;;       ; now a list of specific items
;;       (list
;;         (list specific-button-label-string button-red button-green button-blue
;;               specific-x specific-y specific-z)
;;         (list specific-button-label-string button-red button-green button-blue
;;               specific-x specific-y specific-z)))))
;; 
(define (fascinating-clusters-gui window-name sorting-options cluster-list)

  ;; utility function
  (define (add-feature-buttons feature-list cluster-vbox)
    (let ((frame (gtk-frame-new "Cluster Features"))
	  (vbox (gtk-vbox-new #f 0)))
      (gtk-box-pack-start cluster-vbox frame #f #f 2)
      (gtk-container-add frame vbox)

      ;; add buttons to vbox for each feature
      ;; 
      (map (lambda (feature)
	    ; (format #t "feature: ~s~%" feature)
	     (let ((button (gtk-button-new-with-label (car feature))))
	       (gtk-signal-connect button "clicked"
				   (lambda ()
				     (set-rotation-centre 
				      (list-ref feature 4)
				      (list-ref feature 5)
				      (list-ref feature 6))))
	       (gtk-box-pack-start vbox button #f #f 1)))
	     feature-list)))
	     
  ;; main body
  (let* ((window (gtk-window-new 'toplevel))
	 (scrolled-win (gtk-scrolled-window-new))
	 (outside-vbox (gtk-vbox-new #f 2))
	 (inside-vbox (gtk-vbox-new #f 0)))
    
    (gtk-window-set-default-size window  300 200)
    (gtk-window-set-title window window-name)
    (gtk-container-border-width inside-vbox 2)
    (gtk-container-add window outside-vbox)
    (gtk-box-pack-start outside-vbox scrolled-win #t #t 0) ; expand fill padding
    (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
    (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)
    
    (let loop ((cluster-list cluster-list)
	       (count 0))
      
      (cond
       ((null? cluster-list) 'done)
       ((= 180 count) 'done)
       (else 

	(let ((cluster-info (car cluster-list)))
	   (let* ((frame (gtk-frame-new #f))
		  (vbox (gtk-vbox-new #f 2)))

	     (gtk-container-border-width frame 6)
	     (gtk-container-add frame vbox)
	     (gtk-box-pack-start inside-vbox frame #f #f 10)
	     (let ((go-to-cluster-button (gtk-button-new-with-label
					  (car cluster-info))))
	       (gtk-signal-connect go-to-cluster-button "clicked"
				   (lambda ()
				     (set-rotation-centre
				      (list-ref cluster-info 1)
				      (list-ref cluster-info 2)
				      (list-ref cluster-info 3))))
	       (gtk-box-pack-start vbox go-to-cluster-button #f #f 2)
	       
	       ;; now we have a list of individual features:
	       (let ((features (list-ref cluster-info 4)))
		 (if (> (length features) 0)
		     (add-feature-buttons features vbox)))
	       (loop (cdr cluster-list) (+ count 1))))))))
		   

    (gtk-container-border-width outside-vbox 2)
    (let ((ok-button (gtk-button-new-with-label "  Close  ")))
      (gtk-box-pack-end outside-vbox ok-button #f #f 0)
      (gtk-signal-connect ok-button "clicked"
			  (lambda args
			    (gtk-widget-destroy window))))
    (gtk-widget-show-all window)))


;;; example:
;;; 
;(molprobity-fascinating-clusters-gui
; "Testing the GUI" 
; (list 
;  (list "Active Site" (list 0 1 2 3 4))
;  (list "Worst First" (list 3 4 1 2 0)))
; (list 
;  (list "The first cluster"
;	11 12 15
;	(list 
;	 (list "A bad thing" 0.4 0.6 0.7 10 13 16)	
;	 (list "Another bad thing" 0.4 0.6 0.7 12 15 16)))
;  (list "Another cluster of baddies"
;	-11 12 15
;	(list
;	 (list "A quite bad thing" 0.4 0.6 0.7 -10 -13 16)	
;	 (list "A not so bad thing" 0.4 0.6 0.7 -12 -15 16)))
;  (list "A third cluster of baddies"
;	11 12 -15
;	(list
;	 (list "A quite bad rotamer" 0.4 0.6 0.7 10 13 -16)	
;	 (list "A hydrogen clash" 0.4 0.6 0.7 12 15 -16)
;	 (list "A not so bad H-H clash" 0.4 0.6 0.7 12 15 -16)))))
