##############################################################################
USER    := $(shell whoami)
HOME     = /Users/$(USER)
#CFLAGS   = -g
#CFLAGS   = -g -Wall -Wextra -Wpedantic \
           -Wformat=2 -Wno-unused-parameter -Wshadow \
           -Wwrite-strings -Wstrict-prototypes -Wold-style-definition \
           -Wredundant-decls -Wnested-externs -Wmissing-include-dirs \
		   -fsanitize=bounds -fsanitize-undefined-trap-on-error
CFLAGS   = -g #-O3
ARCH     =  x86_64
INCLS    = -I./include #-I/opt/local/include
LIBS     = -lm #-L/opt/local/lib -lgsl -lgslcblas
CC       =  gcc
PROGRAM  =  gday


SOURCES  =  \
$(PROGRAM).c version.c read_param_file.c read_met_file.c litter_production.c \
utilities.c plant_growth.c photosynthesis.c water_balance.c \
water_balance_sub_daily.c simple_moving_average.c soils.c optimal_root_model.c \
initialise_model.c write_output_file.c phenology.c disturbance.c canopy.c \
radiation.c zbrent.c odeint.c nrutil.c rkqs.c rkck.c

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
