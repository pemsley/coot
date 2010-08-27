
(let loop ((height-scale 0.7)
	   (var-scale 0.5))
  
  (cond 
;    ((> height-scale 0.8) 'done)
   ((> var-scale 1.3) 'done)
   (else 
    (let ((r (fit-rama height-scale var-scale)))
      
      (format #t "Final Remainder: at height-scale var-scale ~s ~s ~s~%" 
	      height-scale var-scale r)
      (loop height-scale
	    (+ 0.05 var-scale))))))
    

