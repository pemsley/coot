(use-modules (gtk gtk))

(read-set! keywords 'prefix)

(let* ((window (gtk-widget-new 'GtkWindow
			       :type         'toplevel
			       :title        "hello world"
			       :allow_grow   #f
			       :allow_shrink #f
			       :GtkContainer::border_width 10))
       (label  (gtk-widget-new 'GtkLabel
			       :label        "hello world"
			       :visible      #t))
       (button (gtk-widget-new 'GtkButton
			       :child        label
			       :parent       window
			       :visible      #t)))

  (gtk-signal-connect button "clicked" 
		      (lambda ()
			(display (gtk-object-get label :label))
			(newline)
			(gtk-widget-set label :label "yo!")))
  (gtk-widget-show window)
  (gtk-standalone-main window))
