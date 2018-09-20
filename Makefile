# A list of program
#PROGS = bruss_main robotarm_main vdp_main twelve_main fvdp_main vdp_main_new ms_main ms_main_old

PROGS = mainCLRT  

############## GiNac Package ######################
GINACLIB = -lcln -lginac

############## Eigen Package ######################
EIGENLIB = -I/usr/local/include/eigen3/

############## NLOPTLIB Package ######################
NLOPTLIB = -lnlopt -lm

############ CAPD related #################
# a list of all your units to be linked with your programs
#OTHERS = utils output 
# directory where capd scripts are (e.g. capd-config)
CAPDBINDIR =/Users/mdarifui/Work/capd/bin

# setting compiler and linker flags
CAPDFLAGS = `${CAPDBINDIR}/capd-config --cflags`
CAPDLIBS = `${CAPDBINDIR}/capd-config --libs`
CXXFLAGS += ${CAPDFLAGS} -O2 


# directory where object and dependency files will be created
OBJDIR = .obj/

# Defining Object files
#OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
#OBJ_FILES = ${OTHERS_OBJ} ${PROGS:%=${OBJDIR}%.o}
OBJ_FILES = ${PROGS:%=${OBJDIR}%.o}


#
.PHONY: all
all: ${PROGS}

#rule to link executables
#${PROGS}: % : ${OBJDIR}%.o ${OTHERS_OBJ}
${PROGS}: % : ${OBJDIR}%.o	
	#$(CXX) -o $@ $< ${OTHERS_OBJ} ${CAPDLIBS} ${GINACLIB} ${NLOPTLIB}
	$(CXX) -o $@ $< ${CAPDLIBS} ${GINACLIB} ${NLOPTLIB}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}

#rule to compile cpp files
${OBJ_FILES}: ${OBJDIR}%.o: %.cpp
	@mkdir -p ${OBJDIR}
	$(CXX) $(EIGENLIB) $(CXXFLAGS) -std=c++11 -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<

# rule to clean all object files, dependencies and executables
.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o ${OBJDIR}*.o.d ${PROGS} cd
