;;;; Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 2 of the License, or (at
;;;; your option) any later version.
 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA

(use-modules (ice-9 regex))


;; For each residue in the protein (molecule number @var{imol}), do a
;; rotamer fit and real-space refinement.  Update the graphics and
;; rotate the scene at each residue for eye candy goodness.
;;
;; Note that residue with alt confs do not undergo auto-fit-rotamer.
;; This is because that autofit-rotamer then refine will tend to put
;; both rotamers into the same place.  Not good.  It seems a
;; reasonable expectation that residues with an alternate conformation
;; are already resonably well-fitted.  So residues with alternate
;; conformations undergo only real space refinement.
;; 
(define (fit-protein imol)

  (set-go-to-atom-molecule imol)
  (make-backup imol) ; do a backup first
  (let ((backup-mode (backup-state imol))
	(imol-map (imol-refinement-map))
	(replacement-state (refinement-immediate-replacement-state)))

    (if (= imol-map -1)
	(info-dialog "Oops.  Must set a map to fit")

	(begin
	  (turn-off-backup imol)
	  (set-refinement-immediate-replacement 1)
	  
	  (map (lambda (chain-id)
		 (if (not (is-solvent-chain? imol chain-id))
		     (let ((n-residues (chain-n-residues chain-id imol)))
		       (format #t "There are ~s residues in chain ~s~%" n-residues chain-id)
		       
		       (for-each 
			(lambda (serial-number)
			  
			  (let ((res-name (resname-from-serial-number imol chain-id serial-number))
				(res-no   (seqnum-from-serial-number  imol chain-id serial-number))
				(ins-code (insertion-code-from-serial-number imol chain-id serial-number)))
			    (if (not (string=? res-name "HOH"))
				(map (lambda (alt-conf) 
				       (format #t "centering on ~s ~s ~s~%" chain-id res-no "CA")
				       (set-go-to-atom-chain-residue-atom-name chain-id res-no "CA")
				       (rotate-y-scene 30 0.3) ; n-frames frame-interval(degrees)
				       (if (string=? alt-conf "")
					   (auto-fit-best-rotamer res-no alt-conf ins-code chain-id imol 
								  imol-map 1 0.1))
				       (if (>= imol-map 0)
					   (begin
					     ;; (refine-auto-range imol chain-id res-no "")
					     (refine-zone imol chain-id res-no res-no alt-conf)
					     (accept-regularizement)))
				       (rotate-y-scene 30 0.3))
				     (residue-alt-confs imol chain-id res-no ins-code)))))
			(number-list 0 (- n-residues 1))))))
	       (chain-ids imol))
	  
	  (if (= replacement-state 0)
	      (set-refinement-immediate-replacement 0))
	  (if (= backup-mode 1)
	      (turn-on-backup imol))))))

;; For each residue in chain chain-id of molecule number imol, do a
;; rotamer fit and real space refinement of each residue.  Don't
;; update the graphics while this is happening (which makes it faster
;; than fit-protein, but much less interesting to look at).
;; 
(define (fit-chain imol chain-id)

  (make-backup imol)
  (let ((backup-mode (backup-state imol))
	(imol-map (imol-refinement-map))
	(alt-conf "")
	(replacement-state (refinement-immediate-replacement-state)))

    (turn-off-backup imol)
    (set-refinement-immediate-replacement 1)

    (if (= imol-map -1)
	(format #t "WARNING:: fit-chain undefined imol-map. Skipping~%")
	(let ((n-residues (chain-n-residues chain-id imol)))
	  (for-each 
	   (lambda (serial-number)
	     
	     (let ((res-name (resname-from-serial-number imol chain-id serial-number))
		   (res-no   (seqnum-from-serial-number  imol chain-id serial-number))
		   (ins-code (insertion-code-from-serial-number imol chain-id serial-number)))
	       (format #t "centering on ~s ~s ~s~%" chain-id res-no "CA")
	       (set-go-to-atom-chain-residue-atom-name chain-id res-no "CA")
	       (auto-fit-best-rotamer res-no alt-conf ins-code chain-id imol imol-map 1 0.1)
	       (if (>= imol-map 0)
		   (begin
		     (refine-zone imol chain-id res-no res-no alt-conf)
		     (accept-regularizement)))))
	   (number-list 0 (- n-residues 1)))))
    
    (if (= replacement-state 0)
	(set-refinement-immediate-replacement 0))
    (if (= backup-mode 1)
	(turn-on-backup imol))))

(define (fit-residue-range imol chain-id resno-start resno-end)
  
  (make-backup imol)
  (let ((backup-mode (backup-state imol))
	(imol-map (imol-refinement-map))
	(alt-conf "")
	(replacement-state (refinement-immediate-replacement-state)))

    (turn-off-backup imol)
    (set-refinement-immediate-replacement 1)

    (if (= imol-map -1)
	(format #t "WARNING:: fit-chain undefined imol-map. Skipping~%")
	(let ((n-residues (chain-n-residues chain-id imol))
	      (ins-code ""))
	  (for-each 
	   (lambda (res-no)
	     (format #t "centering on ~s ~s ~s~%" chain-id res-no "CA")
	     (set-go-to-atom-chain-residue-atom-name chain-id res-no "CA")
	     (auto-fit-best-rotamer res-no alt-conf ins-code chain-id imol imol-map 1 0.1)
	     (if (>= imol-map 0)
		 (begin
		   (refine-zone imol chain-id res-no res-no alt-conf)
		   (accept-regularizement))))
	   (number-list resno-start resno-end))))
	     
    (if (= replacement-state 0)
	(set-refinement-immediate-replacement 0))
    (if (= backup-mode 1)
	(turn-on-backup imol))))
    

;; For each residue in the solvent chains of molecule number
;; @var{imol}, do a rigid body fit of the water to the density.
;; 
(define (fit-waters imol . animate?)

  (format #t "animate?: ~s~%" animate?)
  (let ((imol-map (imol-refinement-map))
	(do-animate? (if (null? animate?) #f #t)))

    (format #t "do-animate?: ~s~%" do-animate?)
    (if (= imol-map -1)
	(add-status-bar-text "You need to define a map to fit the waters")
	(let ((replacement-state (refinement-immediate-replacement-state))
	      (backup-mode (backup-state imol))
	      (alt-conf ""))

	  (turn-off-backup imol)
	  (set-refinement-immediate-replacement 1)
	  (set-go-to-atom-molecule imol)
	  
	  ;; refine waters
	  (let ((chain-identifiers (chain-ids imol)))
	    
	    (for-each 
	     (lambda (chain-id)
	       (if (is-solvent-chain? imol chain-id)
		   
		   (let ((n-residues (chain-n-residues chain-id imol)))
		     
		     (format #t "There are ~s residues in chain ~s~%" n-residues chain-id)
		     (for-each 
		      (lambda (serial-number)
			
			(let ((res-no (seqnum-from-serial-number imol chain-id serial-number)))
			  
			  (if do-animate?
			      (let ((res-info (residue-info imol chain-id res-no "")))
				(if (not (null? res-info))
				    (let ((atom (car res-info)))
				      (let ((atom-name (car (car atom))))
					
					(set-go-to-atom-chain-residue-atom-name chain-id res-no atom-name)
					(refine-zone imol chain-id res-no res-no alt-conf)
					(rotate-y-scene 30 0.6))) ; n-frames frame-interval(degrees)
				    (refine-zone imol chain-id res-no res-no alt-conf)))
			      (begin
				(refine-zone imol chain-id res-no res-no alt-conf)))

			  (accept-regularizement)))
		      (number-list 0 (- n-residues 1))))))
	     chain-identifiers))
	  
	  (if (= replacement-state 0)
	      (set-refinement-immediate-replacement 0))
	  (if (= backup-mode 1)
	      (turn-on-backup imol))))))


;; Step through the residues of molecule number imol and at each step
;; do a residue range refinement (unlike fit-protein for example,
;; which does real-space refinement for every residue).
;; 
;; The step is set internally to 2.
;; 
(define (stepped-refine-protein imol . res-step)

  (set-go-to-atom-molecule imol)
  (make-backup imol)
  (let ((res-step (if (and (list res-step)
			   (not (null? res-step))
			   (number? (car res-step)))
		      (car res-step)
		      2))
	(backup-mode (backup-state imol))
	(imol-map (imol-refinement-map))
	(alt-conf "")
	(replacement-state (refinement-immediate-replacement-state)))

    (if (= imol-map -1)
	(add-status-bar-text "Oops.  Must set a map to fit")

	;; we jump through this hoop with range-step because
	;; (inexact->exact (/ 1 2) now returns 1/2.  In guile 1.6.x it
	;; returned 1
	(let* ((step-inter (inexact->exact (/ (- res-step 1) 2)))
	       (range-step (if (integer? step-inter) 
			       step-inter
			       (+ step-inter (/ 1 2)))))
	  (turn-off-backup imol)
	  (set-refinement-immediate-replacement 1)
	  (set-refine-auto-range-step range-step)
	  
	  (map (lambda (chain-id)
		 (let ((n-residues (chain-n-residues chain-id imol)))
		   (format #t "There are ~s residues in chain ~s~%" n-residues chain-id)

		   (for-each
		    (lambda (serial-number)

		      (let ((res-name (resname-from-serial-number imol chain-id serial-number))
			    (res-no   (seqnum-from-serial-number  imol chain-id serial-number))
			    (ins-code (insertion-code-from-serial-number imol chain-id serial-number)))
			(format #t "centering on ~s ~s ~s~%" chain-id res-no "CA")
			(set-go-to-atom-chain-residue-atom-name chain-id res-no "CA")
			(rotate-y-scene 30 0.3) ; n-frames frame-interval(degrees)
			(if (>= imol-map 0)
			    (begin
			      (refine-auto-range imol chain-id res-no "")
			      (accept-regularizement)))
			(rotate-y-scene 30 0.3)))
		    (every-nth (number-list 0 (- n-residues 1)) res-step))))
	       (chain-ids imol))

	  (if (= replacement-state 0)
	      (set-refinement-immediate-replacement 0))
	  (if (= backup-mode 1)
	      (turn-on-backup imol))))))

;; The GUI that you see after ligand finding. 
;; 
(define post-ligand-fit-gui 
  (lambda ()

    (molecules-matching-criteria 
     (lambda (imol)
       (if (not (valid-model-molecule? imol))
	   #f
	   (let ((name (molecule-name imol)))
	     (if (string-match "Fitted ligand" name)
		 (cons name (molecule-centre imol))
		 #f)))))))
    
;; test-func is a function given one argument (a molecule number) that
;; returns either #f if the condition is not satisfied or something
;; else if it is.  And that "something else" can be a list like
;; (list label x y z)
;; or 
;; (list "Bad Chiral" 0 "A" 23 "" "CA" "A")
;; 
;; It is used in the create a button label and "what to do when the
;; button is pressed".
;; 
(define (molecules-matching-criteria test-func)
    
  ;; first count the number of fitted ligands, and post this if is
  ;; is greater than 0.
  
  (let ((passed-molecules 
	 (let loop ((molecule-numbers (molecule-number-list))
		    (passed-molecules '()))
	   
	   (cond
	    ((null? molecule-numbers) passed-molecules)
	    (else
	     (let ((var (test-func (car molecule-numbers))))
	       (if (eq? var #f)
		   (loop (cdr molecule-numbers)
			 passed-molecules)
		   (loop (cdr molecule-numbers)
			 (cons (list var (car molecule-numbers))
			       passed-molecules)))))))))

    (if (null? passed-molecules)
	
	;; no matching molecules
	(add-status-bar-text "No matching molecules!")
	
	;; OK, proceed.
	(let* ((window (gtk-window-new 'toplevel))
	       (scrolled-win (gtk-scrolled-window-new))
	       (outside-vbox (gtk-vbox-new #f 2))
	       (inside-vbox (gtk-vbox-new #f 0)))
	  
	  (gtk-window-set-default-size window 200 140)
	  (gtk-window-set-title window "Fitted Ligands")
	  (gtk-container-border-width inside-vbox 2)
	  
	  (gtk-container-add window outside-vbox)
	  (gtk-box-pack-start outside-vbox scrolled-win #t #t 0) ; expand fill padding
	  (gtk-scrolled-window-add-with-viewport scrolled-win inside-vbox)
	  (gtk-scrolled-window-set-policy scrolled-win 'automatic 'always)

	  ; (format #t "debug:: passed-molecules ~s~%" passed-molecules)
	  ;; ((("Fitted ligand #9" 68.4 11.9 4.6) 21)
	  ;;  (("Fitted ligand #8" 68.3 12.8 8.1) 20))
	  
	  (let loop ((passed-molecules passed-molecules))
	    
	    (cond 
	     ((null? passed-molecules) 'done)
	     (else 
	      (let* ((imol (car (cdr (car passed-molecules))))
		     (var (car (car passed-molecules)))
		     (name (molecule-name imol)))
		
		(let ((button (gtk-button-new-with-label name)))
		  (gtk-box-pack-start inside-vbox button #f #f 1)
		  (gtk-signal-connect button "clicked"
				      (lambda args 
					(let ((s (format #f "Centred on ~a" name)))
					  (add-status-bar-text s)
					  (format #t "debug:: var: ~s~%" var)
					  (apply set-rotation-centre (cdr var))))))
		(loop (cdr passed-molecules))))))
	  
	  (gtk-container-border-width outside-vbox 2)
	  (let ((ok-button (gtk-button-new-with-label "OK")))
	    (gtk-box-pack-end outside-vbox ok-button #f #f 0)
	    (gtk-signal-connect ok-button "clicked"
			      (lambda args
				(gtk-widget-destroy window))))
	  
	  (gtk-widget-show-all window)))))

;; This totally ignores insertion codes.  A clever algorithm would
;; need a re-write, I think.  Well, we'd have at this end a function
;; that took a chain-id res-no-1 ins-code-1 res-no-2 ins-code-2 
;; 
;; And refine-zone would need to be re-written too, of course.  So
;; let's save that for a rainy day (days... (weeks)).
;; 
(define (refine-active-residue-generic side-residue-offset)

  (let ((active-atom (active-residue)))
    
    (if (not active-atom)
	(format #t "No active atom~%")
	(let ((imol      (list-ref active-atom 0))
	      (chain-id  (list-ref active-atom 1))
	      (res-no    (list-ref active-atom 2))
	      (ins-code  (list-ref active-atom 3))
	      (atom-name (list-ref active-atom 4))
	      (alt-conf  (list-ref active-atom 5)))

	  (format #t "active-atom: ~s~%" active-atom)
	  (let ((imol-map (imol-refinement-map))
		(replacement-state (refinement-immediate-replacement-state)))
	    
	    (if (= imol-map -1)
		(info-dialog "Oops.  Must Select Map to fit to!")
		
		(begin
		  (set-refinement-immediate-replacement 1)
		  (refine-zone imol chain-id 
			       (- res-no side-residue-offset)
			       (+ res-no side-residue-offset)
			       alt-conf)
		  (accept-regularizement)))
		  
	    (if (= replacement-state 0)
		(set-refinement-immediate-replacement 0)))))))
		  

(define (refine-active-residue)
  (refine-active-residue-generic 0))

(define (refine-active-residue-triple)
  (refine-active-residue-generic 1))

;; Another cool function that needs a key binding.
;; 
(define (auto-fit-rotamer-active-residue)

  (let ((active-atom (active-residue)))
    (if (not active-atom)
	(format #t "No active atom~%")
	(let ((imol      (list-ref active-atom 0))
	      (chain-id  (list-ref active-atom 1))
	      (res-no    (list-ref active-atom 2))
	      (ins-code  (list-ref active-atom 3))
	      (atom-name (list-ref active-atom 4))
	      (alt-conf  (list-ref active-atom 5)))

	  (format #t "active-atom: ~s~%" active-atom)
	  (let ((imol-map (imol-refinement-map)))
	    
	    (if (= imol-map -1)
		(info-dialog "Oops.  Must Select Map to fit to!")
		
		(auto-fit-best-rotamer res-no alt-conf ins-code chain-id imol imol-map 1 0.1)))))))

		  
