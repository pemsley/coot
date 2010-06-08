
(if (defined? 'coot-main-menubar)

    (let* ((menu (coot-menubar-menu "PISA"))
	   (submenu-pisa (gtk-menu-new))
	   (menuitem-pisa (gtk-menu-item-new-with-label "PISA Assemblies...")))

      (add-simple-coot-menu-menuitem
       menu "PISA assemblies..." 
       (lambda ()
	 (molecule-chooser-gui "Choose molecule for PISA assembly analysis"
			       (lambda (imol)
				 (pisa-assemblies imol)))))

      (add-simple-coot-menu-menuitem
       menu "PISA interfaces..." 
       (lambda ()
	 (molecule-chooser-gui "Choose molecule for PISA interfaces analysis"
			       (lambda (imol)
				 (pisa-interfaces imol)))))))


