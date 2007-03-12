;;;; Copyright 2004, 2006 by Paul Emsley, The University of York
 
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

;;; Add here settings for your group - rather than expecting each user
;;; to add them to their own .coot files.
;;; 

(let ((os-type (vector-ref (uname) 0)))

  (if (string=? os-type "Linux")
      (set-display-lists-for-maps 0))

  (if (string=? os-type "Linux")
      (set-browser-interface "firefox"))

  (if (string=? os-type "Darwin")
      (set-browser-interface "open")))

;; York setting for molprobity. 
;; 
(set! *probe-command* "/y/people/emsley/coot/Linux/bin/probe.2.11.050121.linux.RH9")
(set! *reduce-command* "/y/people/emsley/coot/Linux/bin/reduce.2.21.030604")

    
