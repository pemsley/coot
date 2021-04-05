
;;;; Copyright 2004, 2005, 2006 by The University of York
 
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


;; enable debugging.  We don't care about the speed of scheme code,
;; but we do care about the ease of debugging.
(debug-enable 'debug)
(debug-enable 'backtrace)

;; globals 

;; This is full pathname of molprobity's probe program
(define *probe-command* "probe")

;; This is full pathname of molprobity's reduce program
(define *reduce-command* "reduce")

;; This has be be here (in a general place) because tips-gui (where it
;; used to be is conditionally loaded).
;; (default to tips-gui displayed is true).
(define *do-coot-tips-flag* #f)

(define load-by-search 
  (lambda (file)

    (let ((f (%search-load-path file)))
      (if f
	  (begin
	    ;; (format #t "load ~s~%" file) ;; remove verbose
	    (load f))
	  (format #t "Error finding ~s~%" file)))))

;; set the prompt to be coot, not guile.
;; (set-repl-prompt! "coot> ")

;; Note the position of coot-gui is important.  It seem that if there
;; are too many files in front of it (even blank ones!) coot barfs
;; when it gets to coot-gui.scm.
;; 
;; 20060203 I have now enabled coot in scripting mode (no graphics
;; (--no-graphics command line option)).  In that case, we need to not
;; load up scheme files which load up gtk (currently coot-gui and
;; tips-gui).
;; 
(define load-all-scheme
  (lambda (use-gui?)

    (let ((pre-list (list "filter.scm" 
			  "matrices.scm"
			  "coot-utils.scm"
			  "json-reader.scm"))
	  (post-list (list "coot-lsq.scm"
			   "shelx.scm"
			   "get-ebi.scm" 
			   "local-code.scm" 
			   "hello.scm"  ; this gives some people problems.  If you 
					; get problems too, just comment out 
					; "hello.scm", nothing bad will happen.
			   "mutate.scm" "refmac.scm"
			   "brute-lsqman.scm" "libcheck.scm" "gap.scm"
			   "fitting.scm"
			   "raster3d.scm" "povray.scm"
			   "remote-control.scm"
			   "generic-objects.scm"
			   "fascinating-things.scm"
			   "ncs.scm"
			   "parse-pisa-xml.scm"
			   "cns2coot.scm"
			   "clear-backup.scm"
			   "tips.scm"
			   "3d-generator-import.scm"
			   "dictionary-generators.scm"
			   "cho-restraints-from-models.scm"
			   "add-linked-cho.scm"
			   "jligand.scm"
			   "americanisms.scm"
			   "group-settings.scm")))

      (let ((scheme-list 
	     (cond 
	      ((eq? #t use-gui?) 
	       (append pre-list (list "coot-gui.scm") post-list (list "tips-gui.scm"
								      ;; "check-for-updates.scm"
								      "gui-hole.scm"
								      "gui-prosmart.scm"
								      "gui-add-linked-cho.scm"
								      "jligand-gui.scm"
								      "get-recent-pdbe.scm"
								      "extensions.scm"
								      "shelx-extensions.scm"
								      "enhanced-ligand.scm"
								      "ligand-check.scm"
								      "acedrg-link.scm"
								      "sharpen-blur.scm"
								      "interactive-nudge-residues.scm"
								      "gui-ligand-sliders.scm"
								      "find-baddies.scm"
								      )))
	      (else 
	       (append pre-list post-list)))))
	      
	(map load-by-search scheme-list)))))

