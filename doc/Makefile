TARGETS = vcsbeam_documentation.pdf

all: $(TARGETS)

%.pdf: %.tex
	pdflatex $(basename $<)
	makeindex $(basename $<)
	bibtex $(basename $<)
	pdflatex $(basename $<)
	pdflatex $(basename $<)
