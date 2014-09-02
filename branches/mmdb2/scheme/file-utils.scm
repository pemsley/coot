
(define (my-file-info fd/port/fname . maybe-chase?)
  (let ((chase? (:optional maybe-chase? #t)))
    (let ((info (if (or chase?
			(not (string? fd/port/fname)))
		    (stat fd/port/fname)
		    (lstat fd/port/fname))))
      (make-file-info (stat:type info)
		      (stat:dev info)
		      (stat:ino info)
		      (stat:mode info)
		      (stat:nlink info)
		      (stat:uid info)
		      (stat:gid info)
		      (stat:size info)
		      (stat:atime info)
		      (stat:mtime info)
		      (stat:ctime info)))))

(define (file-not-exists? fd/port/fname . maybe-chase?)
  (with-errno-handler
      ((err data)
       ((errno/acces) 'search-denied)
       ((errno/noent errno/notdir) #t))
    (apply my-file-info fd/port/fname maybe-chase?)
    #f))


