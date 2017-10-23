# MATLAB Makefile

FILE = simulatedMovingBed

all: run clean

run:
	matlab -nodesktop -nosplash -r "cd ..; installCADET; cd $(FILE); $(FILE); quit"

smb:
	matlab -nodesktop -nosplash -nodisplay -nojvm -r "cd ..; installCADET; cd $(FILE); smbOperatOpt; quit"

clean: 
	rm *.m~ *.txt~ Makefile~
	@echo "all cleaned up"
