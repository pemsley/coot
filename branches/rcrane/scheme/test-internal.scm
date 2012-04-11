
(let ((r (test-internal)))
   (format #t "test-internal returned status: ~s ~%" r)
   if r
   (coot-real-exit 0)
   (coot-real-exit 1))


   
