
(define coot-povray-file-name "coot.pov")
(define coot-povray-png-file-name "coot-povray.png") ; see coot-povray in callbacks.c
(define coot-povray-log-file-name "coot-povray.log")

; (define coot-png-display-program "display") ; defined in raster3d.scm

;; args not including the output filename 
;; 
(define povray-args
  (lambda ()
    
    ; this is directory that contains colors.inc
    (list "+L/sw/share/povray-3.5/include"
	  "+FN16")))

;; Run provray using current displayed image and write .pov file to
;; default filename
;; 
(define povray-image
  (lambda ()
    (povray coot-povray-file-name)
    (format #t "calling povray with args: ~s~%" (povray-args))
    (let ((extra-args (list "+H600" "+W600"))) ; just a guess
    ; goosh-command cmd args data-list log-file-name screen-output-also?
      (goosh-command "povray" (append (povray-args) extra-args) 
		     '() coot-povray-log-file-name #t)
      (format #t "INFO:: displaying...~%")
      (run-concurrently coot-png-display-program coot-povray-png-file-name))))

