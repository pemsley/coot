
(use-modules (ice-9 threads)) ; must come after coot-gui or else coot-gui barfs.

;; a place-holder for the coot listener socket
(define %coot-listener-socket #f)

;; Open a coot listener socket, return nothing much, but do set! 
;; %coot-listener-socket.
;; 
;; Hmmm... why make a side-effect like that?  Why not return
;; %coot-listener-socket so that the caller can set it?  There may be
;; a reason...
;; 
;; And the reason is that I can then call
;; coot-listener-idle-function-proc without having to use a c++
;; variable.
;; 
(define (open-coot-listener-socket port-number host-name)

    (format #t "in open-coot-listener-socket port: ~s host: ~s~%"
	     port-number host-name)

    (let ((soc (socket AF_INET SOCK_STREAM 0))
	  (host-address "127.0.0.1"))
      
      (connect soc AF_INET (inet-aton host-address) port-number)
      
      (display "Coot listener socket ready!\n" soc)
      (set! %coot-listener-socket soc)
;        (if soc 
; 	   (set-coot-listener-socket-state-internal 1))))
      ))

;      (call-with-new-thread coot-listener-idle-function-proc
;			    coot-listener-error-handler)
      


;; Do this thing when idle
;; 
;; currently set to listen to the %coot-listener-socket
;; 
(define coot-listener-idle-function-proc
  (lambda ()
    
    (if %coot-listener-socket
	(begin
	  (listen-coot-listener-socket-v2 %coot-listener-socket)
	  (usleep 2000))
	(begin
	  (format #t "coot-listener-idle-function-proc bad sock: ~s~%" 
		  %coot-listener-socket)))))

;; must return 
(define (listen-coot-listener-socket-v2 soc)

  (define evaluate-char-list
    (lambda (read-bits)
      
      ;; 
      ; (format #t "so far we have: ~s~%" (reverse read-bits))

      (let ((s (list->string (reverse read-bits))))
	;; is it python?
	(if (and (> (string-length s) 8) ; (checked in order)
		 (string-match "# PYTHON" (substring s 0 8)))
	    (safe-python-command s)
	    ;; not python
	    (begin 
	      (format #t "this string will be evaluated: ~s~%" s)
	      (with-input-from-string s
		(lambda () 
		  
		  (let loop ((thing (read)))
		    (if (not (eof-object? thing))
			(begin 
			  (format #t "this code be evaluated: ~s~%" thing)
			  (eval thing (interaction-environment))
			  (loop (read))))))))))))


  (if (not (char-ready? soc))
      (begin
;   	(format #t "nothing on the line...~%")
	(usleep 10000)
	#t)
      (let f ((c (read-char soc))
	      (read-bits '()))
	(cond
	 ((eof-object? c) (format #t "server gone\n"))
	 (else 
	  (let ((ready-char-flag (char-ready? soc)))
	   ; (format #t "DEBUG: ready-char-flag: ~s~%" ready-char-flag)
	    (if ready-char-flag
		(begin
		  (f (read-char soc) (cons c read-bits)))
		(begin
		  ;; was the end? Let's wait a bit to see if the
		  ;; buffer gets refilled - and if it still is not
		  ;; ready, then let's evaluate what we have...
		  (usleep 500000)
		  (format #t "waiting for new char...\n")
		  (let ((ready-char-2-flag (char-ready? soc)))
		    (if ready-char-2-flag 
			(begin
			  (f (read-char soc) (cons c read-bits)))
			(begin (evaluate-char-list read-bits))))))))))))
			  
