# gnuplot Makefile

FILE_GP = *.gp

all: plots pdf clean

plots: 
	gnuplot $(FILE_GP)

pdf:
	pdflatex -interaction=batchmode trajectory_raffinate.tex
	pdflatex -interaction=batchmode trajectory_extract.tex
	pdflatex -interaction=batchmode chromatogram.tex

clean:
	rm *.aux *.log *-inc.eps *-inc-eps-converted-to.pdf *.tex
