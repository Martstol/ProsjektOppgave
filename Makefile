all: pdf

pdf:
	pdflatex main
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main

clean:
	rm *.aux *.log *.bbl *.blg *.lof *.lot *.toc
