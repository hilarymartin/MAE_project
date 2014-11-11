LDFLAGS =  -lm

CC = gcc

#Can set POPULATION
# by adding items to OPTFLAGS using -D.
# E.g., change the following line (with no hash sign) to
# OPTFLAGS = -DPOPULATION=60000 
OPTFLAGS = 

####################################################################
####  Shouldn't need to change anything beyond this point       ####
####################################################################

ALL_EXE = asp

ASP_OBJ = asp.o dbaccess.o list_functions.o table_download.o printing.o

all: $(ALL_EXE)

asp: $(ASP_OBJ)
	$(CC) $(ASP_OBJ) $(DBLIBS) $(LDFLAGS) -o asp

clean:
	rm -f $(ASP_OBJ) asp

# DO NOT DELETE THIS LINE -- make depend depends on it.

