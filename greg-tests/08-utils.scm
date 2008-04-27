
(greg-testcase "Test key symbols" #t
   (lambda ()

     (add-key-binding "name" "missing key" (lambda () 'flunk))

     (test-list (list 
		 (cons (key-sym-code 'a-symbol) -1)
		 (cons (key-sym-code (list 'a 'b)) -1)
		 (cons (key-sym-code "A") 65)
		 (cons (key-sym-code "a") 97)
		 (cons (key-sym-code "9") 57)
		 (cons (key-sym-code "cent") 162)))))


