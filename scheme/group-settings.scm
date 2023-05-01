;;;; Copyright 2004, 2006 The University of York
;;;; Author: Paul Emsley
 
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

;;; Add here settings for your group - rather than expecting each user
;;; to add them to their own .coot files.
;;; 

(let ((os-type (vector-ref (uname) 0)))

  (if (string=? os-type "Linux")
      (set-display-lists-for-maps 1))

  (if (string=? os-type "Linux")
      (set-browser-interface "firefox"))

  (if (string=? os-type "Darwin")
      (begin
	(set-display-lists-for-maps 0)
	(set-browser-interface "open"))))

;; (set! *probe-command* "probe")
;; (set! *reduce-command* "reduce")
