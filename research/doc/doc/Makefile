PDFLATEX = pdflatex -interaction=batchmode -synctex=1
SH 		   = /bin/bash
ASCRIPT  = /usr/bin/osascript
SOURCE   = researchplan.tex
BASE     = "$(basename $(SOURCE))"

default: pdf clean

.PHONY: pdf graphics
pdf: $(SOURCE)
	# run pdflatex and bibtex multiple times to get references right
	$(PDFLATEX) $(BASE).tex && bibtex $(BASE) && $(PDFLATEX) $(BASE).tex && $(PDFLATEX) $(BASE).tex

.PHONY: view
view: 
	# open the document in default PDF viewer
	open $(BASE).pdf

.PHONY: clean
clean:
	# remove all TeX-generated files in your local directory
	$(RM) -f -- *.aux *.bak *.bbl *.blg *.log *.out *.toc *.tdo _region.*
