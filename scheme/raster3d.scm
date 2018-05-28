
(define coot-r3d-file-name "coot.r3d")
(define coot-png-file-name "coot.png")
(define coot-png-display-program "display")
(let ((os-type (vector-ref (uname) 0)))
  (if (string=? os-type "Darwin")
      (set! coot-png-display-program "open")))

;; run raster3d
(define render-image
  (lambda ()
    (raster3d coot-r3d-file-name)
    (format #t "calling render...~%")
;     (goosh-command-with-file-input "render" coot-r3d-file-name 
; 				   coot-png-file-name)
    (goosh-command-with-file-input "render" (list "-png" "-labels") 
				   coot-r3d-file-name coot-png-file-name)
		   
    (format #t "calling display...~%")
    (run-concurrently coot-png-display-program coot-png-file-name)
    ))


;; Run either raster3d or povray
;; 
(define raytrace
  (lambda (image-type source-file-name image-file-name x-size y-size)

    (cond
     ((eq? image-type 'raster3d) 
      (begin
	(goosh-command-with-file-input "render" (list "-png" "-labels")
				       source-file-name image-file-name)
	(run-concurrently coot-png-display-program image-file-name)))

     ((eq? image-type 'povray) 
      (let ((args (append (povray-args)
			  (list
			   "-UV"
			   (string-append "+W" (number->string x-size))
			   (string-append "+H" (number->string y-size))
			   source-file-name))))
	(format #t "INFO:: povray with args: ~s~%" args)
	(goosh-command "povray" args '() coot-povray-log-file-name #t)
	(run-concurrently coot-png-display-program coot-povray-png-file-name)))
     (else 
      (format #t "Image type ~s unknown~%" image-type)))))


