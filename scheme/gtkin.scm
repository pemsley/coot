;; A Simple relativistic two body decay calculator and display.

;;    Copyright (C) 2000  Brett Viren
;;
;;    This program is free software; you can redistribute it and/or modify
;;    it under the terms of the GNU General Public License as published by
;;    the Free Software Foundation; either version 2 of the License, or
;;    (at your option) any later version.
;;
;;    This program is distributed in the hope that it will be useful,
;;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;    GNU General Public License for more details.
;;
;;    You should have received a copy of the GNU General Public License
;;    along with this program; if not, write to the Free Software
;;    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


(use-modules (gnome gnome)
	     (gtk gtk)
	     (gtk gdk))

(read-set! keywords 'prefix)
(gnome-init-hack "canvas" #f '())

;(load "kin.scm")

(define (make-line parent col)
;  (display "make-line\n")
  (gnome-canvas-item-new parent 'GnomeCanvasLine
                         :points #(10 10 20 20)
                         :fill_color col))

(define (set-line item pts)
;  (display pts) (newline)
  (gnome-canvas-item-set item :points (list->array 1 pts)))



(define (make-disp)
  (let* ((vbox (gtk-vbox-new #f 0))
	 (zoom-a (gtk-adjustment-new 0.05 0.001 1 0.01 0.05 0.0))
	 (zoom-s (gtk-hscale-new zoom-a))
	 (canvas (gnome-canvas-new))
	 (root (gnome-canvas-root canvas))
	 (group (gnome-canvas-item-new root 'GnomeCanvasGroup
				       :x 150 :y 150))
	 (v (make-line group "black"))
	 (v1 (make-line group "red"))
	 (v2 (make-line group "blue")))

    (gtk-widget-set-usize canvas 300 300)
    (gnome-canvas-set-pixels-per-unit canvas 0.05)

    (gtk-container-add vbox canvas)
    (gtk-container-add vbox zoom-s)
    
    (gtk-signal-connect zoom-a "value-changed"
			(lambda ()
			  (let ((val (gtk-adjustment-value zoom-a)))
			    (gnome-canvas-set-pixels-per-unit 
			     canvas val))))
    
    (list vbox v v1 v2)))


(define (update-disp disp vmag vec1 vec2)
;  (display "vmag=") (display vmag) (newline)
  (set-line (cadr disp)  (list (- vmag) 0.0 0.0 0.0))
;  (display "vec1=") (display vec1) (newline)
  (set-line (caddr disp) (list 0.0 0.0 (car vec1) (cadr vec1)))
;  (display "vec2=") (display vec2) (newline)
  (set-line (cadddr disp) (list 0.0 0.0 (car vec2) (cadr vec2))))


(define (make-window)
  (let* ((window (gtk-window-new 'toplevel))
	 (vbox (gtk-vbox-new #f 0))
	 (hbox #f)
	 	 
	 ;; user changeable
	 (angle-adjustment (gtk-adjustment-new 90 0 179.9 1 10 0))
	 (angle-scale (gtk-hscale-new angle-adjustment))
	 (momentum-adjustment (gtk-adjustment-new 0 0 10000 100 1000 1000))
	 (momentum-scale (gtk-hscale-new momentum-adjustment))
	 (m-adjustment  (gtk-adjustment-new 139.6 0 1000 1 10 10))
	 (m-spin        (gtk-spin-button-new m-adjustment 1 1)) 
	 (m1-adjustment (gtk-adjustment-new 105.7 0 1000 1 10 10))
	 (m1-spin       (gtk-spin-button-new m1-adjustment 1 1)) 
	 (m2-adjustment (gtk-adjustment-new 0 0 1000 1 10 10))
	 (m2-spin       (gtk-spin-button-new m2-adjustment 1 1)) 
	 
	 ;; user readonly
	 (momentum1-adjustment (gtk-adjustment-new 0 0 10000 100 1000 1000))
	 (momentum1-scale (gtk-hscale-new momentum1-adjustment))
	 (momentum2-adjustment (gtk-adjustment-new 0 0 10000 100 1000 1000))
	 (momentum2-scale (gtk-hscale-new momentum2-adjustment))
	 (separation-adjustment (gtk-adjustment-new 0 0 180 1 10 10))
	 (separation-scale (gtk-hscale-new separation-adjustment))
	 
	 (disp (make-disp))

	 (button (gtk-button-new-with-label "Quit")))

    ;;Called whenever something changes - updates outputs
    (define (register-change)
      (let* ((P (gtk-adjustment-value momentum-adjustment))
	     (a (gtk-adjustment-value angle-adjustment))
	     (m (gtk-adjustment-value m-adjustment))
	     (m1 (gtk-adjustment-value m1-adjustment))
	     (m2 (gtk-adjustment-value m2-adjustment))
	     (res (lab-2body-decay P m m1 m2 (deg-to-rad a))))
;	(display "r=") (display res) (newline)
	(gtk-adjustment-set-value momentum1-adjustment (magnitude (car res)))
	(gtk-adjustment-set-value momentum2-adjustment (magnitude (cadr res)))
	(gtk-adjustment-set-value separation-adjustment 
				  (rad-to-deg (angle-between
					       (car res) 
					       (cadr res))))
	(update-disp disp P (car res) (cadr res))
))

    (gtk-signal-connect momentum-adjustment "value_changed" register-change)
    (gtk-signal-connect angle-adjustment "value_changed" register-change)
    (gtk-signal-connect m-adjustment "value_changed" register-change)
    (gtk-signal-connect m1-adjustment "value_changed"  register-change)
    (gtk-signal-connect m2-adjustment "value_changed"  register-change)
    
    
    (gtk-widget-set-usize separation-scale 150 30)
    (gtk-widget-set-usize momentum-scale 150 30)
    (gtk-widget-set-usize momentum1-scale 150 30)
    (gtk-widget-set-usize momentum2-scale 150 30)

    (gtk-container-add window vbox)

    (gtk-container-add vbox (car disp))

    (set! hbox (gtk-hbox-new #f 0))
    (gtk-container-add vbox (gtk-label-new "Center of Mass values"))
    (gtk-container-add vbox hbox)
    (gtk-container-add hbox (gtk-label-new "Angle (deg):"))
    (gtk-container-add hbox angle-scale)

    (set! hbox (gtk-hbox-new #f 0))
    (gtk-container-add vbox hbox)
    (gtk-container-add hbox (gtk-label-new "M (MeV/c^2):"))
    (gtk-container-add hbox m-spin)
    (gtk-container-add hbox (gtk-label-new "P (MeV/c):"))
    (gtk-container-add hbox momentum-scale)

    (set! hbox (gtk-hbox-new #f 0))
    (gtk-container-add vbox hbox)
    (gtk-container-add hbox (gtk-label-new "m1 (MeV/c^2):"))
    (gtk-container-add hbox m1-spin)
    (gtk-container-add hbox (gtk-label-new "m2 (MeV/c^2):"))
    (gtk-container-add hbox m2-spin)


    (gtk-container-add vbox (gtk-hseparator-new))
    (gtk-container-add vbox (gtk-label-new "Lab coordinate values"))

    (set! hbox (gtk-hbox-new #f 0))
    (gtk-container-add vbox hbox)
    (gtk-container-add hbox (gtk-label-new "p1 (MeV/c):"))
    (gtk-container-add hbox momentum1-scale)
    (gtk-container-add hbox (gtk-label-new "p2 (MeV/c):"))
    (gtk-container-add hbox momentum2-scale)

    (set! hbox (gtk-hbox-new #f 0))
    (gtk-container-add vbox hbox)
    (gtk-container-add hbox (gtk-label-new "separation angle(1,2) (deg):"))
    (gtk-container-add hbox separation-scale)

    (gtk-container-add vbox button)
    (gtk-signal-connect button "clicked"
			(lambda ()
			  (gtk-widget-destroy window)))
    (gtk-widget-show-all window)
    (gtk-standalone-main window)))

(define (lorentz vec beta)
  (let* ((gamma (/ 1.0 (sqrt (- 1 (* beta beta)))))
	 (gammabeta (* gamma beta))
	 (vt (+ (* gamma (car vec)) (* gammabeta (cadr vec))))
	 (vx (+ (* gamma (cadr vec)) (* gammabeta (car vec))))
	 (vy (caddr vec)))
    (list vt vx vy)))

(define (energy-com m1 m2 m)
  (let ((e (/ (+ (* m m)
		 (- (* m2 m2))
		 (* m1 m1))
	      (* 2 m))))
;    (display (list "m1=" m1 "m2=" m2 "m=" m "e=" e)) (newline)
    e))

(define (decay-2body m m1 m2 theta)
  (let* ((e1 (energy-com m1 m2 m))
	 (p1 (sqrt (- (* e1 e1) (* m1 m1)))))
;    (display "e1= ") (display e1) (display " p1=") (display p1) (newline)
    (cons (list e1 (* p1 (cos theta)) (* p1 (sin theta)))
	  (list (- m e1) (* -1.0 p1 (cos theta)) (* -1.0 p1 (sin theta))))))

(define (decay-2body-in-flight P m m1 m2 theta)
  (let* ((beta (/ P (sqrt (+ (* P P) (* m m)))))
	 (prod (decay-2body m m1 m2 theta))
	 (v1 (car prod))
	 (v2 (cdr prod)))
;    (display "prod=") (display prod) (newline)
    (cons (lorentz v1 beta) (lorentz v2 beta))))

(define (dot r1 r2)
  (+ (* (car r1) (car r2)) (* (cadr r1) (cadr r2))))

(define (angle-between r1 r2)
;  (display r1) (newline)
;  (display r2) (newline)
  (acos (/ (dot r1 r2) (sqrt (* (dot r1 r1) (dot r2 r2))))))

(define (magnitude r1)
  (sqrt (dot r1 r1)))

(define (lab-2body-angle P m m1 m2 theta)
  (let* ((prod (decay-2body-in-flight P m m1 m2 theta))
	 (r1 (cdr (car prod)))
	 (r2 (cdr (cdr prod))))
;    (display prod) (newline)
;    (display r1) (newline)
;    (display r2) (newline)
    (angle-between r1 r2)))

(define (lab-2body-decay P m m1 m2 theta)
  (let* ((prod (decay-2body-in-flight P m m1 m2 theta))
    (r1 (cdar prod))
    (r2 (cddr prod)))
  (list r1 r2)))

(define (rad-to-deg r)
  (/ (* r 180.0) 3.141592682))

(define (deg-to-rad r)
  (* (/ r 180.0) 3.141592682))


(make-window)

  