all: pdf

pdf:
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main

clean:
	rm *.aux *.log *.bbl *.blg
