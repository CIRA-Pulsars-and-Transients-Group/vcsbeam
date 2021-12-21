TARGETS = vcsbeam_documentation.pdf

all: $(TARGETS) doxygen

%.pdf: %.tex
	pdflatex $(basename $<)
	makeindex $(basename $<)
	bibtex $(basename $<)
	pdflatex $(basename $<)
	pdflatex $(basename $<)

doxygen:
	doxygen && make -C latex
