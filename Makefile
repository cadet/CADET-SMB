# MATLAB Makefile

FILE = simulatedMovingBed

all: run clean

run:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); $(FILE) quit"

clean: 
	rm *.m~
	@echo "all cleaned up"
