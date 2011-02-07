
EXTRA_DIST = coot.tex derivatives.tex tutorial.tex coot-plain-2.ps \
	function-filter.scsh functions.tex derivatives-part.tex \
	reference.tex crib-sheet.tex user-manual.texi \
	coot-reference-manual.texi 3-amino-benzamide-2.eps  3-amino-benzamide-2.png  


user-manual.dvi: $(top_srcdir)/doc/user-manual.texi
	if [ ! -e coot-plain-2.eps ] ; then cp $(top_srcdir)/doc/coot-plain-2.ps coot-plain-2.eps ; fi
	texi2dvi $^

user-manual.ps: user-manual.dvi
	dvips -o $@ $^

user-manual.pdf: user-manual.ps
	ps2pdf $^

user-manual.html: $(top_srcdir)/doc/user-manual.texi
	texi2html $^
	splitstr='-split_chapter' ; \
	texi2html --version ; \
	if [ $$? = 0 ] ; then splitstr='-split=chapter' ; fi ; \
	texi2html $$splitstr $^ ; \
	cd ..

coot.html: $(top_srcdir)/doc/user-manual.texi
	makeinfo --html --no-split --force $^
	makeinfo --html            --force $^

coot.dvi: $(top_srcdir)/doc/coot.tex
	if [ ! -e coot-plain-2.ps ] ; then cp $(top_srcdir)/doc/coot-plain-2.ps . ; fi
	-latex -interaction=nonstopmode $^
	-makeindex coot
	-latex -interaction=nonstopmode $^
	latex -interaction=nonstopmode $^

booklet.ps: user-manual.dvi
	dvips $^ -f | psbook | psnup -2 | psset -t > $@

coot-manual.ps: user-manual.dvi
	dvips $^ -o $@

functions.tex: $(top_srcdir)/src/c-interface.h
# scsh -s $(top_srcdir)/doc/function-filter.scsh $^ > $@
	echo passing on making function.tex

functions.texi: $(top_srcdir)/src/c-interface.h
# scsh -s $(top_srcdir)/doc/function-filter.scsh $^ > $@
	echo passing on making functions.texi

reference.dvi: $(top_srcdir)/doc/reference.tex functions.tex
	-latex -interaction=nonstopmode $^
	latex -interaction=nonstopmode $^

coot-reference-manual.dvi: $(top_srcdir)/doc/coot-reference-manual.texi functions.texi
	texi2dvi $<

coot-reference-manual.html: $(top_srcdir)/doc/coot-reference-manual.texi functions.texi
	texi2html $<

reference-booklet.ps: reference.dvi
	dvips $^ -f | psbook | psnup -2 | psset -t > $@

crib-sheet.ps: crib-sheet.dvi
	dvips -o $@ $^

crib-sheet.pdf: crib-sheet.ps
	ps2pdf $^

coot-keys-and-buttons.pdf: crib-sheet.pdf
	cp $^ $@

crib-sheet.dvi: $(top_srcdir)/doc/crib-sheet.tex
	latex $^

tutorial.dvi: $(top_srcdir)/doc/tutorial.tex
	-latex -interaction=nonstopmode $^
	latex -interaction=nonstopmode $^

# this relies on the image directory being in place with ps images in it:
tutorial.ps: tutorial.dvi
	dvips -o $@ $^

tutorial.pdf: tutorial.ps
	ps2pdf $^


tutorial-2.dvi: $(top_srcdir)/doc/tutorial-2.tex
	-latex -interaction=nonstopmode $^
	latex -interaction=nonstopmode $^

# this relies on the image directory being in place with ps images in it:
tutorial-2.ps: tutorial-2.dvi
	dvips -o $@ $^

tutorial-2.pdf: tutorial-2.ps
	ps2pdf $^

