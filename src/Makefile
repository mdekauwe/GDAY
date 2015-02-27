##############################################################################
HOME     = /Users/mdekauwe
#CFLAGS   = -g -Wall -Wformat -fsanitize=bounds -fsanitize-undefined-trap-on-error -O0#-Wextra
CFLAGS   = -O3
ARCH     =  x86_64
INCLS    = -I./include #-I/opt/local/include
LIBS     = -lm #-L/opt/local/lib -lgsl -lgslcblas
CC       =  gcc
PROGRAM  =  gday


SOURCES  =  \
$(PROGRAM).c version.c read_param_file.c read_met_file.c \
litter_production.c utilities.c plant_growth.c mate.c water_balance.c \
simple_moving_average.c soils.c optimal_root_model.c initialise_model.c \
write_output_file.c phenology.c

OBJECTS = $(SOURCES:.c=.o)
RM       =  rm -f
##############################################################################

# top level create the program...
all: 		$(PROGRAM)


version.c:
		git rev-parse HEAD | awk ' BEGIN {print "#include \"version.h\""} \
		                    {print "const char *build_git_sha = \"" $$0"\";"}\
		                     END {}' > version.c
		date | awk 'BEGIN {} {print "const char *build_git_time = \""$$0"\";"}\
		                     END {} ' >> version.c

# Compile the src file...
$(OBJECTS):	$(SOURCES)
		$(CC) ${INCLS} $(CFLAGS) -c $(SOURCES)

# Linking the program...
$(PROGRAM):	$(OBJECTS)
		$(CC) $(OBJECTS) $(LIBS) ${INCLS} $(CFLAGS) -o $(PROGRAM)

clean:
		$(RM) $(OBJECTS) $(PROGRAM) version.c

install:
		cp $(PROGRAM) $(HOME)/bin/$(ARCH)/.
		$(RM) $(OBJECTS)
##############################################################################
