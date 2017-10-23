# MATLAB Makefile

FILE = simulatedMovingBed

all: run clean

run:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); $(FILE); quit"

smb:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); smbOperatOpt; quit"

clean: 
	rm *.m~ *.txt~ Makefile~
	rm examples/Forward/*.m~ examples/Optimization/*.m~
	@echo "all cleaned up"
