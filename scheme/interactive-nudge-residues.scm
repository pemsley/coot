
(define (new-molecule-with-nudged-residues imol residue-spec residue-delta nudge-by)

  (let ((imol-new (copy-molecule imol))
	(chain-id (residue-spec->chain-id residue-spec))
	(resno-start (- (residue-spec->res-no residue-spec) residue-delta))
	(resno-end   (+ (residue-spec->res-no residue-spec) residue-delta)))

    (if #f ;; debugging
	(begin
	  (format #t "imol: ~s~%" imol)
	  (format #t "residue-spec: ~s~%" residue-spec)
	  (format #t "residue-delta: ~s~%" residue-delta)
	  (format #t "resno-start: ~s~%" resno-start)
	  (format #t "resno-end: ~s~%"   resno-end)
	  (format #t "nudge-by: ~s~%" nudge-by)))

    (let ((status (nudge-residue-sequence imol-new chain-id resno-start resno-end nudge-by 1)))

      ;;(format #t "nudge-residue-sequence status: ~s~%" status)

      (if (= status 0) ;; fail
	  (let ((s (string-append "Failed to nudge around " chain-id " "
				  (number->string (residue-spec->res-no residue-spec )))))
	    (add-status-bar-text s)
	    (close-molecule imol-new)
	    -1) ;; return a bad new molecule id.
	  (begin
	    imol-new)))))

(define (nudge-residues-gui imol residue-spec)

  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 4))
	 (hbox-0 (gtk-hbox-new #f 4))
	 (hbox-1 (gtk-hbox-new #f 4))
	 (hbox-2 (gtk-hbox-new #f 4))
	 (label-1 (gtk-label-new " Nudge by "))
	 (label-2 (gtk-label-new " residues"))
	 (label-n (gtk-label-new " Nudge "))
	 (res-lab (string-append " residues up and down from "
				 (residue-spec->chain-id residue-spec) " "
				 (number->string (residue-spec->res-no residue-spec))))
	 (label-a (gtk-label-new res-lab))
	 (m-lab (string-append " Nudging residues from Molecule:\n   "
			       (number->string imol) ": "
			       (molecule-name imol)))
	 (label-m (gtk-label-new m-lab))
	 (entry   (gtk-entry-new))
	 (h-sep (gtk-hseparator-new))
	 (adj (gtk-adjustment-new 0 -30 59 0.01 1 29)) ; value lower upper step-inc page-inc page-size
	 (slider (gtk-hscale-new adj))
	 (buttons-hbox (gtk-hbox-new #f 4))
	 (cancel-button (gtk-button-new-with-label " Close "))
	 (residue-delta 5))

    (gtk-window-set-title window "Nudge Residues")
    (gtk-box-pack-start vbox hbox-0 #f #f 4)
    (gtk-box-pack-start vbox hbox-1 #f #f 4)
    (gtk-box-pack-start vbox hbox-2 #f #f 4)
    (gtk-box-pack-start vbox h-sep #f #f 4)
    (gtk-box-pack-start vbox buttons-hbox #f #f 4)
    (gtk-box-pack-start hbox-0 label-m #f #f 4)
    (gtk-box-pack-start hbox-1 label-n #f #f 2)
    (gtk-box-pack-start hbox-1 entry   #f #f 2)
    (gtk-box-pack-start hbox-1 label-a #f #f 2)
    (gtk-box-pack-start hbox-2 label-1 #f #f 4)
    (gtk-box-pack-start hbox-2 slider  #t #t 4)
    (gtk-box-pack-start hbox-2 label-2 #f #f 4)
    (gtk-box-pack-end   buttons-hbox cancel-button #f #f 4)
    (gtk-container-add window vbox)
    (gtk-window-set-default-size window 350 100)
    (gtk-entry-set-text entry "5")
    (gtk-widget-set-usize entry 30 -1)
    (gtk-scale-set-digits slider 0)

    (gtk-signal-connect cancel-button "clicked" (lambda () (gtk-widget-destroy window)))

    (gtk-signal-connect adj "value_changed"
			(lambda ()
			  (let* ((v (gtk-adjustment-value adj))
				 (v-int (inexact->exact v))
				 (rdt (gtk-entry-get-text entry))
				 (rd (string->number rdt)))
			    (if (number? rd)
				(set! residue-delta rd))

			    ;; (format #t "v: ~s~%" v)
			    ;; (format #t "rdt: ~s~%" rdt)
			    ;; (format #t "rd: ~s~%" rd)
			    ;; (format #t "residue-delta: ~s~%" residue-delta)

			    (let ((imol-new (new-molecule-with-nudged-residues imol residue-spec residue-delta v-int)))
			      imol-new)))) ;; not very interesting

    (gtk-widget-show-all window)))
