HOCKING-svm-compare.pdf: HOCKING-svm-compare.tex refs.bib figure-norm-data.tex figure-hard-margin.tex aistats2014.sty figure-simulation-samples.tex figure-norm-level-curves.tex
	rm -f *.aux *.bbl
	pdflatex HOCKING-svm-compare
	bibtex HOCKING-svm-compare
	pdflatex HOCKING-svm-compare
	pdflatex HOCKING-svm-compare
figure-norm-data.tex: figure-norm-data.R tikz.R
	R --no-save < $<
figure-hard-margin.tex: figure-hard-margin.R tikz.R
	R --no-save < $<
simulation.samples.RData: simulation.samples.R svmlight.R
	R --no-save < $<
figure-simulation-samples.tex: figure-simulation-samples.R simulation.samples.RData tikz.R
	R --no-save < $<
figure-norm-level-curves.tex: figure-norm-level-curves.R tikz.R simulation.RData
	R --no-save < $<
simulation.RData: simulation.R svmlight.R
	R --no-save < $<
