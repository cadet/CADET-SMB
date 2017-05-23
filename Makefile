# MATLAB Makefile

FILE = simulatedMovingBed

all: run clean

run:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); cascade; quit"

clean: 
	rm *.m~
	@echo "all cleaned up"
