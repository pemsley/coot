;;;; Copyright 2000 by Paul Emsley 
;;;; Copyright 2006 by Paul Emsley, The University of York
 
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
;;;; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
;;;; 02110-1301, USA

;; Basic scheme function, filter the objects in list ls by function
;; fn.  e.g. (filter even? (list 0 1 2 3) -> '(0 2)
;; 
(define (filter fn ls)

  (let f ((ls ls)
	  (acc '()))
    (cond
     ((null? ls) (reverse acc))
     ((fn (car ls))  (f (cdr ls) (cons (car ls) acc)))
     (else 
      (f (cdr ls) acc)))))

