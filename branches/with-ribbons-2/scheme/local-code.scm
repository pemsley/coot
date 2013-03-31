
(define do-code 
  (lambda (id-code)

    (let* ((pdb-code-str (string-append id-code ".pdb"))
	   (cif-code-str (string-append id-code ".cif"))
	   (imol (handle-read-draw-molecule pdb-code-str)))
      (read-cif-data cif-code-str imol))))

