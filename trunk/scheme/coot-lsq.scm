;;;; Copyright 2005, 2006 by The University of York
 
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
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


;; Internal type conversion for LSQ fitting.  Return a number
;; according to the symbol match-type-in
;; 
(define (lsq-match-type-symbol->number match-type-in)

  (cond
   ((number? match-type-in) match-type-in)
   ((symbol? match-type-in)
    (cond 
     ((eq? 'CA match-type-in) 2)
     ((eq? 'ca match-type-in) 2)
     ((eq? 'Ca match-type-in) 2)
     ((eq? 'Main match-type-in) 1)
     ((eq? 'main match-type-in) 1)
     ((eq? 'mainchain match-type-in) 1)
     ((eq? 'Mainchain match-type-in) 1)
     ((eq? 'ALL match-type-in) 0)
     ((eq? 'all match-type-in) 0)
     ((eq? 'All match-type-in) 0)))
   (else -1))) ; unknown

;; Create matchers, 7 elements: 
;;   (list ref-start-resno ref-end-resno ref-chain-id imol-ref
;; 	   mov-start-resno mov-end-resno mov-chain-id imol-mov
;; 	   match-type)
;; 
(define set-match-element
  (lambda (m)
    (if (= 7 (length m))
	(let ((match-type-in (list-ref m 6)))
	  (let ((match-type (lsq-match-type-symbol->number match-type-in)))
	    (add-lsq-match (list-ref m 0)
			   (list-ref m 1)
			   (list-ref m 2)
			   (list-ref m 3)
			   (list-ref m 4)
			   (list-ref m 5)
			   (list-ref m 6)
			   match-type)))
	(format #t "Wrong number of elements in match (was ~s should be 7)~%"
		(length m)))))

;; The scripting interface to LSQ matching.  Pass molecule numbers for
;; the reference (imol-ref) and moving (imol-moving) molecules and a
;; match list.  The match list format is described in the manual.
;; 
(define (lsq-match imol-ref imol-moving match-list)
    
  (clear-lsq-matches)
  (map set-match-element match-list)
  
  (apply-lsq-matches imol-ref imol-moving))

;; Simple interface to LSQ fitting.  More often than not this is what
;; you will want, I imagine, e.g. (simple-lsq-match 940 950 "A" 0 940
;; 950 "A" 1 'main)
;; 
(define (simple-lsq-match ref-start-resno ref-end-resno ref-chain-id imol-ref mov-start-resno mov-end-resno mov-chain-id imol-mov match-type)

  (let ((internal-match-type (lsq-match-type-symbol->number match-type)))
    (clear-lsq-matches)
    (add-lsq-match ref-start-resno ref-end-resno ref-chain-id
		   mov-start-resno mov-end-resno mov-chain-id
		   internal-match-type)
    (apply-lsq-matches imol-ref imol-mov)))

		    
; (simple-lsq-match 940 950 "A" 0 940 950 "A" 1 'main)
