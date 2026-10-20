# Build a chapter:   make ch01        (runs the R code, then compiles the PDF)
# Build everything:  make all
# Remove LaTeX junk: make clean

CHAPTERS := $(sort $(wildcard ch[0-9][0-9]))

.PHONY: all clean $(CHAPTERS)

all: $(CHAPTERS)

$(CHAPTERS):
	Rscript $@/code/$@.R
	cd $@ && latexmk -pdf -interaction=nonstopmode -halt-on-error $@.tex

clean:
	for ch in $(CHAPTERS); do (cd $$ch && latexmk -c); done
