;;;; Copyright 2010 by the University of Oxford
;;;; Copyright 2013, 2014, 2015 by Medical Research Council

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
;;;; Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA


(use-modules (oop goops) 
	     (oop goops describe))

;;; Allow the user to set these variables in their .mapview file if
;;; they want some server other than the default choice.  (These
;;; values cribbed from pdb-mode.el, btw (thanks Charlie Bond)).
;;; 
;;;
;;; Sill the OCA Server is bugged so that the send-x-pdb (compresed
;;; file) returns a short file.  So we have to use the uncompressed
;;; form which returns a file that had META tags.  Argh.

(define pdbe-server "https://www.ebi.ac.uk")
(define pdbe-pdb-file-dir "pdbe-srv/view/files")

(define pdbe-file-name-tail "ent")

;; sf exmaple http://www.ebi.ac.uk/pdbe-srv/view/files/r4hrhsf.ent


;; 20151126-PE No, we can't have coot-download created on coot-startup, it must be
;;             made only when we need it.
;; (define coot-tmp-dir (get-directory "coot-download"))

; e.g. (ebi-get-pdb "1crn")
; 
; no useful return value
;
; data-type can be 'pdb or 'sfs (structure factors).  We might like to use
; 'coordinates rather than 'pdb in the future.
; 
; The optional argument imol-coords-arg-list is necessary for
; ouptutting sfs, because we need coordinates from which we can
; calculate phases.
;


(define (net-get-url my-url file-name)
  ;; (format #t "URL:: ~s~%" my-url)
  (coot-get-url my-url file-name))


;; check the directory and get url url-string.
;; 
(define check-dir-and-get-url
  (lambda (dir file-name url-string)

      (if (file-exists? dir)
	  (begin
	    (if (is-directory? dir)
		(net-get-url url-string file-name)
		(format #t "ERROR:: Oops - Can't write to ~s directory!" dir)))
	  (begin
	    (mkdir dir)
	    (if (is-directory? dir)
		(net-get-url url-string file-name)
		(format #t "ERROR:: Oops - create-directory ~s failed!" dir))))))


;; get url-string for data type 'pdb or 'sfs
;; 
(define get-url-str
  (lambda (id url-string data-type . imol-coords-arg-list)

    (format #t "DEBUG:: in get-url-str: ~s ~s ~s~%" id url-string data-type)
    (let ((coot-tmp-dir (get-directory "coot-download")))

      (cond 
       
       ((eq? data-type 'pdb)
	(let ((pdb-file-name (string-append coot-tmp-dir "/" id ".pdb" "."
					    pdbe-file-name-tail)))
	  (check-dir-and-get-url coot-tmp-dir pdb-file-name url-string)
	  (handle-read-draw-molecule pdb-file-name)))

       ((eq? data-type 'sfs)
	(let ((sfs-file-name (string-append coot-tmp-dir "/" id ".cif"))
	      (imol (get-ebi-pdb id)))
	  
	  (if (and (number? imol)
		   (not (= imol -1)))
	      (begin
		(check-dir-and-get-url coot-tmp-dir sfs-file-name url-string)
		(read-cif-data sfs-file-name (car imol-coords-arg-list))))))

       (else 
	"unknown")))))

	    

(define get-ebi-pdb-and-sfs 
  (lambda (id)
    
    (let ((imol-coords (get-ebi-pdb id)))
      (if (not (number? imol-coords))
	  (format #t "failed at reading coordinates. imol-coords was ~s~%"
		  imol-coords)

	  (if (< imol-coords 0) ; -1 is coot code for failed read.

	      (format #t "failed to read coordinates. ~%")
	      
	      (let* ((down-id (string-downcase id))
		     (url-str (string-append 
			       pdbe-server
			       "/"
			       pdbe-pdb-file-dir
			       down-id
			       "." pdbe-file-name-tail)))

		(get-url-str id url-str 'sfs imol-coords)))))))

;; Return a molecule number on success
;; or not a number (#f) or -1 on error.
;; 
(define (get-ebi-pdb id)

  ;; (format #t "======= id: ~s~%" id)

  (let* ((down-id (string-downcase id))
	 (url-str (string-append
		   pdbe-server
		   "/"
		   pdbe-pdb-file-dir
		   "/"
		   down-id ".ent")))
    
    (let ((url-status (get-url-str id url-str 'pdb)))

      ;; e.g. http://ftp.ebi.ac.uk/pub/databases/pdb + 
      ;;      /validation_reports/cb/1cbs/1cbs_validation.xml.gz

      (if (valid-model-molecule? url-status)
	  (pdb-validate down-id url-status)))))


;(ebi-get-pdb "1crn")


;; Return a list of molecules (i.e. the model molecule and the 2 maps).
;; 
;; or, if it didn't work then return #f
;; 
(define (get-eds-pdb-and-mtz id)

  ;; Gerard DVD Kleywegt says we can find the coords/mtz thusly:
  ;;
  ;; - model = http://eds.bmc.uu.se/eds/sfd/1cbs/pdb1cbs.ent
  ;; - mtz   = http://eds.bmc.uu.se/eds/sfd/1cbs/1cbs_sigmaa.mtz
  ;;
  ;; 20091222
  ;; - newprefix: http://eds.bmc.uu.se/eds/dfs/cb/1cbs/
  ;; 
  ;; URL::       "http://eds.bmc.uu.se/eds/sfd/sa/2sar/pdb2sar.ent"
  ;; URL:: "http://eds.bmc.uu.se/eds/sfd/sa/2sar/2sar_sigmaa.mtz"
  ;;
  ;; 20161010 new prefix
  ;; http://www.ebi.ac.uk/pdbe/coordinates/
  ;; http://www.ebi.ac.uk/pdbe/coordinates/files/1cbs_map.mtz
  ;; 

  ;; Return a list of 3 molecule numbers or #f.
  ;; 
  (define (get-cached-eds-files accession-code)

    (let* ((coot-tmp-dir (get-directory "coot-download"))
	   (down-code (string-downcase accession-code))
	   (pdb-file-name (append-dir-file coot-tmp-dir
					   (string-append "pdb" down-code ".ent")))
	   (mtz-file-name (append-dir-file coot-tmp-dir
					   (string-append down-code "_map.mtz"))))

      (format #t "::::::::: pdb-file-name: ~s~%" pdb-file-name)
      (format #t "::::::::: mtz-file-name: ~s~%" mtz-file-name)

      (if (not (file-exists? pdb-file-name))
	  #f
	  (if (not (file-exists? mtz-file-name))
	      #f
	      (let ((imol (read-pdb pdb-file-name))
		    (imol-map (make-and-draw-map mtz-file-name "FWT" "PHWT" "" 0 0))
		    (imol-map-d (make-and-draw-map mtz-file-name "DELFWT"  "PHDELWT" "" 0 1)))
		(if (not (and (valid-model-molecule? imol)
			      (valid-map-molecule? imol-map)
			      (valid-map-molecule? imol-map-d)))
		    (begin
		      (close-molecule imol)
		      (close-molecule imol-map)
		      (close-molecule imol-map-d)
		      #f)

		    (begin
		      (list imol imol-map imol-map-d))))))))


  (define eds-site "https://www.ebi.ac.uk/pdbe/coordinates")
  (define eds-core "some://thing") ;; for web pages
  (define eds-coords-site "https://www.ebi.ac.uk/pdbe/entry-files/download")

  ;; "1cbds" -> "cb/"
  ;; 
  (define (mid-chars id-code)
    (if (not (string? id-code))
	"//fail//"
	(if (not (= (string-length id-code) 4))
	    "/FAIL/"
	    (string-append (substring id-code 1 3) "/"))))

  ;; main line
  ;; 
  (let ((cached-status (get-cached-eds-files id)))

    (if (list? cached-status)

	cached-status ;; a list of molecules

	(let ((coot-tmp-dir (get-directory "coot-download")))

	  (if (eq? #f coot-tmp-dir)
	      (format #t "Can't make coot-download directory~%")
	      
	      (let* ((down-id (string-downcase id))
		     (target-pdb-file (string-append "pdb" down-id ".ent"))
		     (dir-target-pdb-file (string-append coot-tmp-dir "/" target-pdb-file))
		     (model-url (string-append eds-coords-site "/" target-pdb-file))
		     (target-mtz-file (string-append down-id "_map.mtz"))
		     (dir-target-mtz-file (string-append coot-tmp-dir "/" target-mtz-file))
		     (mtz-url (string-append eds-site "/files/" target-mtz-file))
		     (eds-info-page (string-append eds-core "/cgi-bin/eds/uusfs?pdbCode=" down-id)))

		(print-var model-url)
		(print-var mtz-url)
		
		(let* ((pre-download-info (coot-get-url-as-string eds-info-page))
		       ;; (pre-download-info-sxml (xml->sxml pre-download-info))
		       (s1 (net-get-url model-url dir-target-pdb-file))
		       (s2 (net-get-url mtz-url   dir-target-mtz-file)))

		  ;; (format #t "INFO:: --------------- pre-download-info-sxml: ~s~%" pre-download-info-sxml)
		  ;; (format #t "INFO:: --------------- pre-download-info: ~s~%" pre-download-info)

		  (let ((bad-map-status (string-match "No reliable map available" pre-download-info)))
		    (if bad-map-status
			(let ((s  (string-append 
				   "This map (" down-id ") is marked by the EDS as \"not a reliable map\"")))
			  (info-dialog s))))

		  (format #t "INFO:: read model status: ~s~%" s1)
		  (format #t "INFO:: read mtz   status: ~s~%" s2)
		  
		  (let ((r-imol (handle-read-draw-molecule dir-target-pdb-file))
			(map-1 (make-and-draw-map dir-target-mtz-file "FWT" "PHWT" "" 0 0))
			(map-2 (make-and-draw-map dir-target-mtz-file  "DELFWT"  "PHDELWT" "" 0 1)))
		    (set-scrollable-map map-1)
		    (if (valid-model-molecule? r-imol)
			(list r-imol map-1 map-2)
			#f)))))))))




(define (get-pdb-redo text) 

  (define mid-string (lambda (text-in) (substring text-in 1 3)))

  (if (not (string? text))
      "Pass an accession code"

      ;; 
      (if (not (= (string-length text) 4))

	  "Give an acession code"
	  
	  (let* ((stub (string-append 
			"https://pdb-redo.eu/db/"
			text "/" text "_final"))
		 (pdb-file-name (string-append text "_final.pdb"))
		 (mtz-file-name (string-append text "_final.mtz"))
		 (scm-file-name (string-append text "_final.scm"))
		 (url-pdb (string-append stub ".pdb"))
		 (url-mtz (string-append stub ".mtz"))
		 (url-scm (string-append stub ".scm")))

	    (format #t "getting ~s~%" url-pdb)
	    (let ((status (net-get-url url-pdb pdb-file-name)))
	      (if (not (= status 0))
		  (format #t "Failed to get ~s ~s status ~s ~%" url-pdb pdb-file-name status)))
	    (format #t "getting ~s~%" url-mtz)
	    (let ((status (net-get-url url-mtz mtz-file-name)))
	      (if (not (= status 0))
		  (format #t "Failed to get ~s ~s status ~s ~%" url-mtz mtz-file-name status)))
	    (format #t "getting ~s~%" url-scm)
	    (let ((status (net-get-url url-scm scm-file-name)))
	      (if (not (= status 0))
		  (format #t "Failed to get ~s ~s status ~s ~%" url-scm scm-file-name status)))

	    (read-pdb pdb-file-name)
	    (format #t "make-and-draw-map with ~s~%" mtz-file-name)
	    (make-and-draw-map mtz-file-name "FWT" "PHWT" "" 0 0)
	    (make-and-draw-map mtz-file-name "DELFWT" "PHDELWT" "" 0 1)
	    (let ((anom-map (make-and-draw-map mtz-file-name "FAN" "PHAN" "" 0 1)))
	      (set-map-colour anom-map 0.5 0.5 0))
	    (run-script scm-file-name)))))
