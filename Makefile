# MATLAB Makefile

FILE = simulatedMovingBed

all: run clean

run:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); $(FILE); quit"

clean: 
	rm *.m~ *.txt~
	@echo "all cleaned up"
