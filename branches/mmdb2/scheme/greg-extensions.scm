;;;; Copyright 2007 by The University of York
;;;; Copyright 2007 by Paul Emsley
;;;; Copyright 2007 by The University of Oxford
;;;; 
;;;; This program is free software; you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation; either version 3 of the License, or (at
;;;; your option) any later version.
;;;; 
;;;; This program is distributed in the hope that it will be useful, but
;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;;;; General Public License for more details.
 
;;;; You should have received a copy of the GNU General Public License
;;;; along with this program; if not, write to the Free Software
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;;; 02110-1301, USA

(use-modules (ice-9 greg))

(set-run-state-file-status 0) 

(if (defined? 'coot-main-menubar)

    (let ((menu (coot-menubar-menu "Greg")))

      (add-simple-coot-menu-menuitem
       menu "Run All Greg Tests"
       (lambda ()

	 ;; use set! not define 
	 (set! greg-tools (list "greg-tests"))
	 (set! greg-debug #t)
	 (set! greg-verbose 5)
	 (greg-test-run)
	 (force-output)))

      (add-simple-coot-menu-menuitem
       menu "Run One Greg Test"
       (lambda ()

	 ;; use set! not define 
	 (set! greg-tools (list "greg-one-test"))
	 (set! greg-debug #t)
	 (set! greg-verbose 5)
	 (greg-test-run)
	 (force-output)))))



