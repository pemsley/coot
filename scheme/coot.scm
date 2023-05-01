(use-modules (ice-9 popen)
             (ice-9 string-fun)
             (ice-9 format)
             (ice-9 rdelim)
             (srfi srfi-1))
(use-modules (os process))



(let ((pkgdd (get-pkgdatadir-scm)))
  (let ((l 2))

    (display pkgdd)
    (newline)

    (let ((ls (split-after-char-last #\/ pkgdd list)))
      (let ((dir (car ls))) ;; ends with a /
        (let ((coot-scheme-dir (string-append dir "coot/scheme")))
          (format #t "coot-scheme-dir: ~s~%" coot-scheme-dir)
          (set! %load-path (cons coot-scheme-dir %load-path))

          (let ((pre-list (list "redefine-functions.scm"
                                "filter.scm"
                                "matrices.scm"
                                "coot-utils.scm"
                                "json-reader.scm"))
                (post-list (list "coot-lsq.scm"
                                 "shelx.scm"
                                 "get-ebi.scm"
                                 "local-code.scm"
                                 "hello.scm"  ; this gives some people problems.  If you
                                              ; get problems too, just comment out
                                              ; "hello.scm", nothing bad will happen.
                                 "mutate.scm" "refmac.scm"
                                 "brute-lsqman.scm" "gap.scm"
                                 "fitting.scm"
                                 "raster3d.scm" "povray.scm"
                                 "remote-control.scm"
                                 "generic-objects.scm"
                                 "fascinating-things.scm"
                                 "ncs.scm"
                                 "parse-pisa-xml.scm"
                                 "cns2coot.scm"
                                 "clear-backup.scm"
                                 ;; "tips.scm"
                                 "3d-generator-import.scm"
                                 "dictionary-generators.scm"
                                 "cho-restraints-from-models.scm"
                                 "add-linked-cho.scm"
                                 "jligand.scm"
                                 "americanisms.scm"
                                 ;;; "group-settings.scm"
                             ))

                (gui-list (list ;; "tips-gui.scm"
                                "gui-hole.scm"
                                "gui-prosmart.scm"
                                "gui-add-linked-cho.scm"
                                "jligand-gui.scm"
                                "get-recent-pdbe.scm"
                                "extensions.scm"
                                "shelx-extensions.scm"
                                "enhanced-ligand.scm"
                                "ligand-check.scm"
                                "acedrg-link.scm"
                                "sharpen-blur.scm"
                                "interactive-nudge-residues.scm"
                                "gui-ligand-sliders.scm"
                                "find-baddies.scm"
                                )))

            (let ((list-for-loading (append pre-list post-list)))
              (for-each (lambda (f)
                          (let ((ff (%search-load-path f)))
                            (if ff
                                (begin
                                  ;; (format #t "load ~s~%" ff) ;; remove verbose
                                  (load ff)))))
                        list-for-loading))))))))
