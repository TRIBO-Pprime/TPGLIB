include ../makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = prg

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ       = $(SRC:.$(EXT)=.o)
ALL_OBJS  = $(addprefix $(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(utils_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(stats_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(qsort_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(digis_dir)/$(OBJDIR)/, mod_*.o)

ALL_MODS  = -I$(MODDIR)
ALL_MODS += -I$(utils_dir)/$(MODDIR)
ALL_MODS += -I$(stats_dir)/$(MODDIR)
ALL_MODS += -I$(qsort_dir)/$(MODDIR)
ALL_MODS += -I$(digis_dir)/$(MODDIR)

CFLAGS  = $(ALL_MODS) -fopenmp -fPIC
CFLAGS += -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -finit-real=nan -ffpe-summary=none

LFLAGS  = -lm -lgomp -lpthread

ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
	CFLAGS += -$(OPTC)
endif

ifneq ('$(GPROF)','')
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

%.o:	%.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	@find . -maxdepth 1 -name '*.mod*' -type f -print0 | xargs -0r mv -t ./$(MODDIR)

$(EXE):	$(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)
	@rm $(OBJDIR)/prg.o

mod_morpho.o :
prg.o : mod_morpho.o

#--------------------------------------------------
.PHONY: clean debug gprof all last

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f $(EXE)

debug:
	make "DEBUG=TRUE"

gprof:
	make "GPROF=TRUE"

all:
	(make clean ; make)

