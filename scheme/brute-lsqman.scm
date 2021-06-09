;;;; Copyright 2003 by Charlie Bond, The University of Dundee
 
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
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
;;;; 02110-1301, USA

;; brute-lsqman - run lsqman on chains A of two pdb-files and read in
;; the result to coot. Charlie Bond 2003.
;; Can keep a count of the number of successful runs if necessary
(define lsqman-count 0)
(define brute-lsqman
  (lambda (pdb1-imol
           pdb2-imol)
    (let* (
           (pdbout-filename (string-append "coot-lsqman-" (number->string lsqman-count) ".pdb"))
           (lsqman-command "lsqman")
           (command-line-args (append (list " -b ")))
           (data-lines (list "re m1 coot-tmp1.pdb"
                             "re m2 coot-tmp2.pdb"
                             "brute m1 a m2 a 50 25 100"
                             "imp m1 * m2 *"
                             "app m1 m2"
                             "wr m2 " pdbout-filename
                             "quit"))
           (lsqman-log (string-append "coot-lsqman" (number->string lsqman-count) ".log")))

      (write-pdb-file pdb1-imol "coot-tmp1.pdb")
      (write-pdb-file pdb2-imol "coot-tmp2.pdb")
      (let ((status (goosh-command lsqman-command command-line-args data-lines 
				   lsqman-log #t)))
        (if (and (number? status) (= status 0)) ; lsqman ran OK
                                                ; (need to check for number as error returns #f)
            (begin
            (handle-read-draw-molecule pdbout-filename)
            (set! lsqman-count (+ 1 lsqman-count)))
            (begin
            (format #t "lsqman failed - sorry. I don't know what to say.")))))))
