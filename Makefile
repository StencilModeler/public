# Makefile for model
#

ifdef DEBUG
  OPTS+=-DDEBUG
endif

ifdef PLOT
	OPTS+=-DPLOT
endif

ifdef NEHALEM
  OPTS+=-DNEHALEM
endif

ifdef OPTERON
  OPTS+=-DOPTERON
endif

ifdef SHANGHAI
  OPTS+=-DSHANGHAI
endif

ifdef POWER6
  OPTS+=-DPOWER6
endif

ifdef BGP
  OPTS+=-DBGP
endif

ifdef PPC970
  OPTS+=-DPPC970
endif

ifdef SEMI
  OPTS+=-DSEMI
endif

ifdef SANDY
  OPTS+=-DSANDY
endif

ifdef IVY
  OPTS+=-DIVY
endif

ifdef KNC
  OPTS+=-DKNC
endif

ifdef THREADS
  OPTS+=-DTHREADS=$(THREADS)
endif

ifdef THRCORE
  OPTS+=-DTHRCORE=$(THRCORE)
endif

ifdef SEMI
  OPTS+=-DSEMI
  BNAME:=_semi
endif

ifdef INTERP
#  OPTS+=-DINTERP=$(INTERP)
else
#  INTERP=NO_INTERP
endif

ifdef NO_SCALE
#  OPTS+=-DNO_SCALE
#  SCALE:=.NO_SCALE
endif


CC = gcc
CFLAGS = -mtune=core2 -march=core2 -O -g -std=c99 -I./includes $(OPTS)
#CFLAGS += -Wall -Wextra -pedantic 

#PAPI = /usr/local/lib/libpapi.a
LDFLAGS = $(PAPI) -lm

# the line below defines timers.  if not defined, will attempt to automatically
# detect available timers.  See cycle.h.
# should be set to -DHAVE_PAPI or -DHAVE_GETTIMEOFDAY or unset.
#TIMER = -DHAVE_PAPI
TIMER = -DHAVE_GETTIMEOFDAY

BINS = model$(BNAME).$(THREADS)thr.$(THRCORE)xCore

all: $(BINS)

model$(BNAME).$(THREADS)thr.$(THRCORE)xCore: source/model.c source/util.c includes/util.h includes/cycle.h
	$(CC) $(CFLAGS) -DINTERP=NO_INTERP $(TIMER) $(filter %.c,$^) $(LDFLAGS) -o bin/$@.NO_INTERP
	$(CC) $(CFLAGS) -DINTERP=NO_INTERP -DNO_SCALE $(TIMER) $(filter %.c,$^) $(LDFLAGS) -o bin/$@.NO_INTERP.NO_SCALE
	$(CC) $(CFLAGS) -DINTERP=LIN_INTERP $(TIMER) $(filter %.c,$^) $(LDFLAGS) -o bin/$@.LIN_INTERP
	$(CC) $(CFLAGS) -DINTERP=EXP_INTERP $(TIMER) $(filter %.c,$^) $(LDFLAGS) -o bin/$@.EXP_INTERP
	$(CC) $(CFLAGS) -DINTERP=LOG_INTERP $(TIMER) $(filter %.c,$^) $(LDFLAGS) -o bin/$@.LOG_INTERP

clean:
	rm -f *.o bin/$(BINS)

