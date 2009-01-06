
(use-modules (ice-9 greg))

;; use set! not define 
(set! greg-tools (list "greg-tests"))
(set! greg-debug #t)
(set! greg-verbose 5)


(greg-test-run)


