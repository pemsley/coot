
(use-modules (oop goops) (oop goops describe))

;;; Allow the user to set these variables in their .mapview file if
;;; they want some server other than the default choice.  (These
;;; values cribbed from pdb-mode.el, btw (thanks Charlie Bond)).
;;; 
;;;
;;; Sill the OCA Server is bugged so that the send-x-pdb (compresed
;;; file) returns a short file.  So we have to use the uncompressed
;;; form which returns a file that had META tags.  Argh.

;(if (not (defined? 'oca-server))
;(define oca-server "http://oca.ebi.ac.uk")
(define oca-server "http://bip.weizmann.ac.il")
;(define oca-server "http://structure.embl-hamburg.de")
(define oca-server "http://www.ebi.ac.uk/msd-srv/oca")

;(if (not (defined? 'oca-request-stub))
(define oca-pdb-request-stub "oca-bin/save-pdb?id=")
;(define oca-pdb-request-stub "oca-bin/send-pdb?id=")
;(define oca-pdb-request-stub "oca-bin/send-x-pdb?id=") ; get a compressed file

(define pdb-file-name-tail "")
; (define pdb-file-name-tail ".gz") ; if the server gives us a compressed
				  ; file, we need to add .Z to the
				  ; output filename so that mmdb
				  ; realizes it's compressed rather
				  ; than gobbledegook.

(define oca-sfs-request-stub "oca-bin/send-sf?r")
(define oca-sfs-request-tail "sf.ent.Z")

(define coot-tmp-dir "coot-download")


; e.g. (ebi-get-pdb "1crn")
; 
; no useful return value
; 
; Note that that for sf data, we need to construct something like the
; string: http://oca.ebi.ac.uk/oca-bin/send-sf?r2acesf.ent.Z and we
; don't need to strip any html (thank goodness). Also not that the
; accession code now is lower case.
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
    
    (if (eq? data-type 'pdb)
	(begin
	  (let ((pdb-file-name (string-append coot-tmp-dir "/" id ".pdb"
					      pdb-file-name-tail)))
	    (check-dir-and-get-url coot-tmp-dir pdb-file-name url-string)
	    (handle-read-draw-molecule pdb-file-name))))

    (if (eq? data-type 'sfs)
	(begin
	  (let ((sfs-file-name (string-append coot-tmp-dir "/" id ".cif"))
		(imol (get-ebi-pdb id)))
	    
	    (if (and (number? imol)
		     (not (= imol -1)))
		(begin
		  (check-dir-and-get-url coot-tmp-dir sfs-file-name url-string)
		  (read-cif-data sfs-file-name (car imol-coords-arg-list)))))))))
	    

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
			       oca-server
			       "/"
			       oca-sfs-request-stub
			       down-id
			       oca-sfs-request-tail)))

		(get-url-str id url-str 'sfs imol-coords)))))))

;; Return a molecule number on success
;; or not a number (#f) or -1 on error.
;; 
(define (get-ebi-pdb id)

  (let* ((up-id (string-upcase id))
	 (url-str (string-append 
		   oca-server
		   "/"
		   oca-pdb-request-stub
		   up-id)))
    
    (get-url-str id url-str 'pdb)))

;(ebi-get-pdb "1crn")


;;
;; Get data and pdb for accession code id from the Electron Density
;; Server.
;; 
;; @var{id} is the accession code.
;; 
;; 20050725 EDS code
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

  ;; Return a list of 3 molecule numbers or #f.
  ;; 
  (define (get-cached-eds-files accession-code)
    
    (let* ((down-code (string-downcase accession-code))
	   (pdb-file-name (append-dir-file "coot-download"
					   (string-append "pdb" down-code ".ent")))
	   (mtz-file-name (append-dir-file "coot-download"
					   (string-append down-code "_sigmaa.mtz"))))
      
      (if (not (file-exists? pdb-file-name))
	  #f
	  (if (not (file-exists? mtz-file-name))
	      #f
	      (let ((imol (read-pdb pdb-file-name))
		    (imol-map (make-and-draw-map mtz-file-name "2FOFCWT" "PH2FOFCWT" "" 0 0))
		    (imol-map-d (make-and-draw-map mtz-file-name "FOFCWT"  "PHFOFCWT" "" 0 1)))
		(if (not (and (valid-model-molecule? imol)
			      (valid-map-molecule? imol-map)
			      (valid-map-molecule? imol-map-d)))
		    (begin
		      (close-molecule imol)
		      (close-molecule imol-map)
		      (close-molecule imol-map-d)
		      #f)
		    (list imol imol-map imol-map-d)))))))


  (define eds-site "http://eds.bmc.uu.se/eds")

  ;; "1cbds" -> "cb/"
  (define (mid-chars id-code)
    (if (not (string? id-code))
	"//fail//"
	(if (not (= (string-length id-code) 4))
	    "/FAIL/"
	    (string-append (substring id-code 1 3) "/"))))

  ;; main line
  ;; 
  (let ((cached-status (get-cached-eds-files id)))
    (if (not cached-status)
	(let ((r (coot-mkdir coot-tmp-dir)))

	  (if (eq? #f r)
	      (format #t "Can't make directory ~s~%" coot-tmp-dir)
	      
	      (let* ((down-id (string-downcase id))
		     (eds-url (string-append eds-site "/dfs/"))
		     (target-pdb-file (string-append "pdb" down-id ".ent"))
		     (dir-target-pdb-file (string-append coot-tmp-dir "/" target-pdb-file))
		     (mc (mid-chars down-id))
		     (model-url (string-append eds-url mc down-id "/" target-pdb-file))
		     (target-mtz-file (string-append down-id "_sigmaa.mtz"))
		     (dir-target-mtz-file (string-append coot-tmp-dir "/" target-mtz-file))
		     (mtz-url (string-append eds-url mc down-id "/" target-mtz-file)))
		
		(let ((s1 (net-get-url model-url dir-target-pdb-file))
		      (s2 (net-get-url mtz-url   dir-target-mtz-file)))
		  
		  (format #t "INFO:: read model status: ~s~%" s1)
		  (format #t "INFO:: read mtz   status: ~s~%" s2)
		  
		  (let ((r-imol (handle-read-draw-molecule dir-target-pdb-file))
			(sc-map (make-and-draw-map dir-target-mtz-file "2FOFCWT" "PH2FOFCWT" "" 0 0)))
		    (make-and-draw-map dir-target-mtz-file  "FOFCWT"  "PHFOFCWT" "" 0 1)
		    (set-scrollable-map sc-map)
		    (if (valid-model-molecule? r-imol)
			r-imol
			#f)))))))))



