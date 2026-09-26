##############################################################################
CC       = gcc
CFLAGS   = -O3 -Wall -Wno-unused-parameter
#CFLAGS   = -g -O0 -Wall -Wextra -Wno-unused-parameter -Wshadow \
#           -fsanitize=address,undefined
INCLS    = -I./include
LIBS     = -lm
PROGRAM  = gday
ARCH     = x86_64
PREFIX   = $(HOME)/bin/$(ARCH)

SOURCES  =  \
$(PROGRAM).c version.c read_param_file.c read_met_file.c litter_production.c \
utilities.c plant_growth.c photosynthesis.c water_balance.c \
water_balance_sub_daily.c simple_moving_average.c soils.c optimal_root_model.c \
initialise_model.c write_output_file.c phenology.c disturbance.c canopy.c \
radiation.c zbrent.c odeint.c nrutil.c rkqs.c rkck.c

OBJECTS  = $(SOURCES:.c=.o)
DEPS     = $(OBJECTS:.o=.d)
RM       = rm -f
##############################################################################

all: $(PROGRAM)

# Regenerated every build so the git hash written to the outputs is current
version.c: FORCE
	@(git rev-parse HEAD 2>/dev/null || echo unknown) | \
	    awk 'BEGIN {print "#include \"version.h\""} \
	         {print "const char *build_git_sha = \"" $$0 "\";"}' > version.c.tmp
	@date | awk '{print "const char *build_git_time = \"" $$0 "\";"}' \
	    >> version.c.tmp
	@if cmp -s version.c.tmp version.c; then rm version.c.tmp; \
	 else mv version.c.tmp version.c; fi

%.o: %.c
	$(CC) $(INCLS) $(CFLAGS) -MMD -MP -c $< -o $@

$(PROGRAM): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(PROGRAM)

clean:
	$(RM) $(OBJECTS) $(DEPS) $(PROGRAM) version.c version.c.tmp

install: $(PROGRAM)
	cp $(PROGRAM) $(PREFIX)/.

.PHONY: all clean install FORCE
FORCE:

-include $(DEPS)
##############################################################################
