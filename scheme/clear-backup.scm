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

(define (delete-coot-backup-files action-type)
  (let* ((files (glob "*.pdb*" "coot-backup"))
	 (now (current-time))
	 (last-week (- now (* 60 60 24 7))) ;; clear out more than 7 days old.
	 (dir "coot-backup"))
    
    (define add-dir-prefix
      (lambda (file)
	(if (string? file)
	    (append-dir-file dir file)
	    file)))

    (define (old-files-list)
      (let loop ((files files)
		 (old-files '()))
	(cond
	 ((null? files) (reverse old-files))
	 (else
	  (let ((this-mtime (stat:mtime (stat (add-dir-prefix (car files))))))
	    ; (format #t "comparing ~s ~s~%" last-week this-mtime)
	    (if (< this-mtime last-week)
		(loop (cdr files) (cons (car files) old-files))
		(loop (cdr files) old-files)))))))

    (define (file-size file)
      (stat:size (stat (add-dir-prefix file))))

    ; main 
    (let ((olds (old-files-list)))
      (if (eq? action-type 'count)
	  (let* ((total-size (apply + (map file-size olds)))
		 (mb-size-str-1 (number->string (/ total-size (* 1024 1024))))
		 (mb-size-str-2 (if (> (string-length mb-size-str-1) 5)
				    (substring mb-size-str-1 0 5)
				    mb-size-str-1)))
		
	    (format #t "There are ~s old files (~aMb) in ~s~%"
		    (length olds) 
		    mb-size-str-2
		    dir)
	    (list (length olds) mb-size-str-2))

	  (if (eq? action-type 'delete)
	      (begin
		(map (lambda (file)
		       (format #t "deleting old file ~s from ~s~%" file dir)
		       (let ((full-name (add-dir-prefix file)))
			 (delete-file full-name)))
		     olds)
		;; now create a last-backup file with a time stamp:
		(let ((last-cleaned-file (append-dir-file
					  "coot-backup" "last-cleaned")))
		  (call-with-output-file last-cleaned-file
		    (lambda (port)
		      (display (current-time) port)
		      (newline port))))
		#t))))))


;; Make a GUI:
(define (clear-backup-gui)

  (let* ((file-stats (delete-coot-backup-files 'count)))

    ;; more than 1 file to possibly delete?
    (if (> (car file-stats) 1) 

	(let* ((window (gtk-window-new 'toplevel))
	       (frame (gtk-frame-new "Old Backups"))
	       (vbox (gtk-vbox-new #f 3))
	       (hbox (gtk-hbox-new #f 10))
	       (ok-button (gtk-button-new-with-label " Clear up "))
	       (cancel-button (gtk-button-new-with-label " Stay messy "))
	       (label (gtk-label-new (string-append "  There are "
						    (number->string (car file-stats))
						    " old backup files ("
						    (car (cdr file-stats))
						    
						    "Mb)  \n"
						    "Delete Them?"))))
	  (gtk-signal-connect ok-button "clicked"
			      (lambda ()
				(delete-coot-backup-files 'delete)
				(gtk-widget-destroy window)))

	  (gtk-signal-connect cancel-button "clicked"
			      (lambda ()
				;; (coot-real-exit 0)
				(gtk-widget-destroy window)))

	  (gtk-tooltips-set-tip (gtk-tooltips-new) ok-button 
				" Consider yourself patted on the back " "")
	  (gtk-tooltips-set-tip (gtk-tooltips-new) 
				cancel-button
				(string-append 
				 "A less pejorative label here might be \"Keep\" or \"Cancel\" "
				 "but seeing as I like my intestines where they are and not used "
				 "as hosiery fastenings for Systems Adminstrators then we get this "
				 "rather nannying label...")
				"")

	  (gtk-container-add window frame)
	  (gtk-container-border-width frame 6)
	  (gtk-container-add frame vbox)
	  (gtk-box-pack-start vbox label #f #f 6)
	  (gtk-box-pack-start vbox hbox)
	  (gtk-container-border-width hbox 10)
	  (gtk-container-border-width vbox 6)
	  (gtk-box-pack-start hbox ok-button #t #t 6)
	  (gtk-box-pack-start hbox cancel-button #t #t 6)

	  (gtk-widget-set-flags ok-button     '(can-default))
	  (gtk-widget-set-flags cancel-button '(can-default))
	  (gtk-widget-grab-default ok-button)
	  
	  (gtk-widget-show-all window)
	  (gtk-standalone-main window)))))


(define (clear-backups-maybe)

  (let* ((now (current-time))
	 (last-week (- now (* 60 60 24 7)))) ;; Clear out every 7 days... Hmmm.
    (let ((last-cleaned-file (append-dir-file "coot-backup"
					      "last-cleaned")))
      (if (not (file-exists? last-cleaned-file))
	  (clear-backup-gui)
	  (begin
	    (call-with-input-file last-cleaned-file
	      (lambda (port)
		(let ((val (read port)))
		  (if (number? val)
		      (if (< val last-week)
			  (clear-backup-gui)
			  (format #t "INFO:: backup clearout done ~s days ago~%"
				  (/ (- now val) (* 60 60 24)))))))))))))

; (gtk-main)


