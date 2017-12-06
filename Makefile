# MATLAB Makefile

FILE = simulatedMovingBed

all: run clean

run:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); simulatedMovingBed(1,'raffinate'); simulatedMovingBed(2,'raffinate'); quit"

clean: 
	rm *.m~ *.txt~ Makefile~
	@echo "all cleaned up"
