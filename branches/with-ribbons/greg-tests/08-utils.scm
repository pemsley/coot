
(greg-testcase "Test key symbols" #t
   (lambda ()

     (add-key-binding "name" "missing key" (lambda () 'flunk))

     (test-list (list 
		 (cons (key-sym-code 'a-symbol) -1)
		 (cons (key-sym-code (list 'a 'b)) -1)
		 (cons (key-sym-code "A") 65)
		 (cons (key-sym-code "a") 97)
		 (cons (key-sym-code "9") 57)
		 (cons (key-sym-code "cent") 162)
		 (cons (key-sym-code ">") 62)
		 (cons (key-sym-code ";") 59)
		 (cons (key-sym-code "|") 124)))))
		 


(greg-testcase "Test running a Python function" #t 
   (lambda ()

     (if (not (coot-has-python?))
	 #t ; pass without testing
	 (let ((tot (run-python-command "2 + 4")))
	   (if (not (= 6 tot))
	       #f
	       (begin
		 (run-python-command "test_val = 2")
		 (let ((rv (run-python-command "test_val")))
		   (if (not (= rv 2))
		       #f
		       ;; coot does not convert a tuple.
		       (let ((rv2 (run-python-command "test_val_2 = (1,2,3)")))
			 (run-python-command "test_val_2") ; doesn't crash?
			 (format #t "   test_val_2 passed~%")
			 ;; test c.f. running cootaneer_gui_bl() when we don't 
			 ;; have pygtk loaded.  Shouldn't crash. (Roni Gordon bug).
			 (let ((test_val_3 (run-python-command "rr_not_found_thing()")))
			   ;; test_val_3 should be unspecified.
			   (format #t "   test_val_3 passed~%") ; no crash
			   #t))))))))))
			   
		     
(greg-testcase "Internal/External Molecule Numbers match" #t 
   (lambda ()

     (let ((m (molecule-number-list)))
       ;; (format #t "   m: ~s~%" m)
       ;; (format #t " own: ~s~%" (map own-molecule-number m))
       (equal? m (map own-molecule-number m)))))

