
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


(greg-testcase "spacegroup operators to space group conversion" #t 
   (lambda ()

     (let ((test-groups (list
 "P 1" "P -1" "P 1 2 1" "P 1 1 2" "P 2 1 1" "P 1 21 1" "P 1 1 21"
 "P 21 1 1" "C 1 2 1" "A 1 2 1" "I 1 2 1" "A 1 1 2" "B 1 1 2" "I 1 1 2"
 "B 2 1 1" "C 2 1 1" "I 2 1 1" "P 1 m 1" "P 1 1 m" "P m 1 1" "P 1 c 1"
 "P 1 n 1" "P 1 a 1" "P 1 1 a" "P 1 1 n" "P 1 1 b" "P b 1 1" "P n 1 1"
 "P c 1 1" "C 1 m 1" "A 1 m 1" "I 1 m 1" "A 1 1 m" "B 1 1 m" "I 1 1 m"
 "B m 1 1" "C m 1 1" "I m 1 1" "C 1 c 1" "A 1 n 1" "I 1 a 1" "A 1 a 1"
 "C 1 n 1" "I 1 c 1" "A 1 1 a" "B 1 1 n" "I 1 1 b" "B 1 1 b" "A 1 1 n" "I 1 1 a"
 "B b 1 1" "C n 1 1" "I c 1 1" "C c 1 1" "B n 1 1" "I b 1 1" "P 1 2/m 1" "P 1 1 2/m" 
 "P 2/m 1 1" "P 1 21/m 1" "P 1 1 21/m" "P 21/m 1 1" "C 1 2/m 1" "A 1 2/m 1" "I 1 2/m 1"
 "A 1 1 2/m" "B 1 1 2/m" "I 1 1 2/m" "B 2/m 1 1" "C 2/m 1 1" "I 2/m 1 1" "P 1 2/c 1" 
 "P 1 2/n 1" "P 1 2/a 1" "P 1 1 2/a" "P 1 1 2/n" "P 1 1 2/b" "P 2/b 1 1" "P 2/n 1 1"
 "P 2/c 1 1" "P 1 21/c 1" "P 1 21/n 1" "P 1 21/a 1" "P 1 1 21/a" "P 1 1 21/n" "P 1 1 21/b"
 "P 21/b 1 1" "P 21/n 1 1" "P 21/c 1 1" "C 1 2/c 1" "A 1 2/n 1" "I 1 2/a 1" "A 1 2/a 1"
 "C 1 2/n 1" "I 1 2/c 1" "A 1 1 2/a" "B 1 1 2/n" "I 1 1 2/b" "B 1 1 2/b" "A 1 1 2/n"
 "I 1 1 2/a" "B 2/b 1 1" "C 2/n 1 1" "I 2/c 1 1" "C 2/c 1 1" "B 2/n 1 1" "I 2/b 1 1"
 "P 2 2 2" "P 2 2 21" "P 21 2 2" "P 2 21 2" "P 21 21 2" "P 2 21 21" "P 21 2 21"
 "P 21 21 21" "C 2 2 21" "A 21 2 2" "B 2 21 2" "C 2 2 2" "A 2 2 2" "B 2 2 2"
 "F 2 2 2" "I 2 2 2" "I 21 21 21" "P m m 2" "P 2 m m" "P m 2 m" "P m c 21" "P c m 21"
 "P 21 m a" "P 21 a m" "P b 21 m" "P m 21 b" "P c c 2" "P 2 a a" "P b 2 b" "P m a 2"
 "P b m 2" "P 2 m b" "P 2 c m" "P c 2 m" "P m 2 a" "P c a 21" "P b c 21" "P 21 a b"
 "P 21 c a" "P c 21 b" "P b 21 a" "P n c 2" "P c n 2" "P 2 n a" "P 2 a n" "P b 2 n"
 "P n 2 b" "P m n 21" "P n m 21" "P 21 m n" "P 21 n m" "P n 21 m" "P m 21 n"
 "P b a 2" "P 2 c b" "P c 2 a" "P n a 21" "P b n 21" "P 21 n b" "P 21 c n"
 "P c 21 n" "P n 21 a" "P n n 2" "P 2 n n" "P n 2 n" "C m m 2" "A 2 m m" "B m 2 m"
 "C m c 21" "C c m 21" "A 21 m a" "A 21 a m" "B b 21 m" "B m 21 b" "C c c 2" "A 2 a a"
 "B b 2 b" "A m m 2" "B m m 2" "B 2 m m" "C 2 m m" "C m 2 m" "A m 2 m" "A b m 2"
 "B m a 2" "B 2 c m" "C 2 m b" "C m 2 a" "A c 2 m" "A m a 2" "B b m 2" "B 2 m b" "C 2 c m"
 "C c 2 m" "A m 2 a" "A b a 2" "B b a 2" "B 2 c b" "C 2 c b" "C c 2 a" "A c 2 a" "F m m 2"
 "F 2 m m" "F m 2 m" "F d d 2" "F 2 d d" "F d 2 d" "I m m 2" "I 2 m m" "I m 2 m" "I b a 2" 
 "I 2 c b" "I c 2 a" "I m a 2" "I b m 2" "I 2 m b" "I 2 c m" "I c 2 m" "I m 2 a"
 "P m m m"   "P c c m" "P m a a" "P b m b"
"P m m a" "P m m b"
 "P b m m" "P c m m" "P m c m" "P m a m" "P n n a" "P n n b" "P b n n" "P c n n"
 "P n c n" "P n a n" "P m n a" "P n m b" "P b m n" "P c n m" "P n c m" "P m a n"
 "P c c a" "P c c b" "P b a a" "P c a a" "P b c b" "P b a b" "P b a m" "P m c b"
 "P c m a" "P c c n" "P n a a" "P b n b" "P b c m" "P c a m" "P m c a" "P m a b"
 "P b m a" "P c m b" "P n n m" "P m n n" "P n m n"
"P b c n" "P c a n" "P n c a" "P n a b" "P b n a"
 "P c n b" "P b c a" "P c a b" "P n m a" "P m n b" "P b n m" "P c m n" "P m c n"
 "P n a m" "C m c m" "C c m m" "A m m a" "A m a m" "B b m m" "B m m b" "C m c a"
 "C c m b" "A b m a" "A c a m" "B b c m" "B m a b" "C m m m" "A m m m" "B m m m"
 "C c c m" "A m a a" "B b m b" "C m m a" "C m m b" "A b m m" "A c m m" "B m c m"
 "B m a m" 
"F m m m"  "I m m m" "I b a m" "I m c b"
 "I c m a" "I b c a" "I c a b" "I m m a" "I m m b" "I b m m" "I c m m" "I m c m" 
 "I m a m" "P 4" "P 41" "P 42" "P 43" "I 4" "I 41" "P -4" "I -4" "P 4/m" "P 42/m"
 "I 4/m" 
 "P 4 2 2" "P 4 21 2" "P 41 2 2" "P 41 21 2" "P 42 2 2" "P 42 21 2" "P 43 2 2"
 "P 43 21 2" "I 4 2 2" "I 41 2 2" "P 4 m m" "P 4 b m" "P 42 c m" "P 42 n m"
 "P 4 c c" "P 4 n c" "P 42 m c" "P 42 b c" "I 4 m m" "I 4 c m" "I 41 m d" "I 41 c d"
 "P -4 2 m" "P -4 2 c" "P -4 21 m" "P -4 21 c" "P -4 m 2" "P -4 c 2" "P -4 b 2"
 "P -4 n 2" "I -4 m 2" "I -4 c 2" "I -4 2 m" "I -4 2 d" "P 4/m m m" "P 4/m c c"
 "P 4/m b m" "P 4/m n c"
"P 42/m m c"
 "P 42/m c m"
 "P 42/m b c" "P 42/m n m" 
 "I 4/m m m" "I 4/m c m" 
 "P 3" "P 31" "P 32" "P -3"
 "P 3 1 2" "P 3 2 1" "P 31 1 2" "P 31 2 1" "P 32 1 2" "P 32 2 1"
"P 3 m 1" "P 3 1 m" "P 3 c 1" "P 3 1 c" 
 "P -3 1 m" "P -3 1 c" "P -3 m 1" "P -3 c 1"
 "P 6" "P 61" "P 65" "P 62" "P 64" "P 63" "P -6" "P 6/m"
 "P 63/m" "P 6 2 2" "P 61 2 2" "P 65 2 2" "P 62 2 2" "P 64 2 2" "P 63 2 2"
 "P 6 m m" "P 6 c c" "P 63 c m" "P 63 m c" "P -6 m 2" "P -6 c 2" "P -6 2 m"
 "P -6 2 c" "P 6/m m m" "P 6/m c c" "P 63/m c m" "P 63/m m c" "P 2 3" "F 2 3"
 "I 2 3" "P 21 3" "I 21 3" "P m -3" "F m -3"
 "I m -3" "P a -3" "I a -3" "P 4 3 2" "P 42 3 2" "F 4 3 2" "F 41 3 2"
 "I 4 3 2" "P 43 3 2" "P 41 3 2" "I 41 3 2" "P -4 3 m" "F -4 3 m" "I -4 3 m" "P -4 3 n"
 "F -4 3 c" "I -4 3 d" "P m -3 m" "P m -3 n"
 "F m -3 m" "F m -3 c" 
 "I m -3 m" "I a -3 d")))

       ;; failures: "P n n n :1" "P n n n :2"  "P b a n :1"  "P b a n :2" 
       ;; "P n c b :1" "P n c b :2" "P c n a :1" "P c n a :2"  "P m m n :1"
       ;; "P m m n :2"  "P n m m :1"  "P n m m :2" "P m n m :1" "P m n m :2" 
       ;; "R 3 :H"  "R 3 :R" "R -3 :H" "R -3 :R" "R 3 2 :H"  "R 3 2 :R" 
       ;; 
       ;; probable failures (ran out of patience to test the indiviually)
       ;; "C c c a :1" "C c c a :2" "C c c b :1" "C c c b :2" "A b a a :1"
       ;; "A b a a :2" "A c a a :1" "A c a a :2" "B b c b :1" "B b c b :2" "B b a b :1"
       ;; "B b a b :2" "F d d d :1" "F d d d :2"
       ;; "P 4/n :1" "P 4/n :2" "P 42/n :1" "P 42/n :2"
       ;; "I 41/a :1" "I 41/a :2" 
       ;; "P 4/n b m :1" "P 4/n b m :2" "P 4/n n c :1" "P 4/n n c :2" 
       ;;  "P 4/n m m :1" "P 4/n m m :2" "P 4/n c c :1" "P 4/n c c :2" 
       ;;  "P 42/n b c :1" "P 42/n b c :2" "P 42/n n m :1" "P 42/n n m :2" 
       ;; "P 42/n m c :1" "P 42/n m c :2" "P 42/n c m :1" "P 42/n c m :2" 
       ;; "I 41/a m d :1" "I 41/a m d :2" "I 41/a c d :1" "I 41/a c d :2"
       ;;  "P n -3 :1" "P n -3 :2"
       ;; "F d -3 :1"
       ;; "F d -3 m :1" "F d -3 m :2" "F d -3 c :1"
       ;; "P n -3 m :1" "F d -3 :2"
       ;; "R 3 m :H" "R 3 m :R" "R 3 c :H"  "R 3 c :R"  "R -3 m :H" "R -3 m :R"
       ;;  "R -3 c :H" "R -3 c :R" "P n -3 n :1" "P n -3 n :2" 
       ;; "P n -3 m :2"
       ;; "F d -3 c :2" 
       
       ;; clipper doesn't know about these space groups (mmdb returns
       ;; the symop strings correctly (below called symops) (in that
       ;; they correspond to syminfo.lib):
       ;;
       ;; "I 1 21 1" "C 1 21 1"  "B 1 1 m"

       (let ((imol (greg-pdb "monomer-ACT.pdb")))
	 
	 (map (lambda (space-group)

		(let ((set-success (set-space-group imol space-group))
		      (symops (symmetry-operators imol)))

		  (if (not (= set-success 1))
		      (begin 
			(format #t "   bad status on setting space group ~s~%" space-group)
			(throw 'fail)))
		      
		  (if (not (list? symops))
		      (begin 
			(format #t "   bad symops for ~s: ~s~%" space-group symops)
			(throw 'fail)))
			
		  (let* ((derived-HM (symmetry-operators->xHM symops)))
		    ;; (format #t "comparing ~s and ~s~%" space-group derived-HM)

		    (if (not (string? derived-HM))
			(begin 
			  (format #t "   bad derived-HM ~s~%" derived-HM)
			  (throw 'fail)))

		    (if (not (string=? space-group derived-HM))
			(begin
			  (format #t "   No match ~s and ~s~%" space-group derived-HM)
			  (throw 'fail))))))
	      test-groups)
	 #t)))) ;; every passed to get here.
