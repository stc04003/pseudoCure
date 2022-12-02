#! /bin/bash

Rscript -e "knitr::knit('pseudo.Rnw')"

pdflatex pseudo.tex
bibtex pseudo.aux
pdflatex pseudo.tex
pdflatex pseudo.tex

rm *log
rm *out
rm *aux
rm *bbl
rm *blg
rm *cb
rm *cb2
rm *fff
rm *synctex.gz
rm *ttt
