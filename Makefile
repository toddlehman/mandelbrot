#------------------------------------------------------------------------------
#  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
#------------------------------------------------------------------------------

COMPILE =	cc -I/opt/local/include -Wall -std=c11 -O4
LINK =		cc -L/opt/local/lib -lmpfr

#------------------------------------------------------------------------------

PROGS =		mset_image

MODULES =	Common \
                MPReal \
		Memory \
		DeepReal \
		Mandelbrot \
		RGB \
		Palette \
		Pixel \
		Image \
		MainControl

H_FILES =	$(MODULES:%=%.h)

C_FILES =	$(MODULES:%=%.c)

O_FILES =	$(MODULES:%=%.o)

LIBS =

ALL_SRC =	$(H_FILES) \
		$(C_FILES) \
		mset_image \
		mset_zoom_make \
		mset_queue_runner \
		mset_blend_frames \
		Makefile

SNAPSHOT =	snapshot.tar.gz

#------------------------------------------------------------------------------

default:	all

all:		$(PROGS) $(SNAPSHOT)

distclean:	clean

clean:
	@echo Removing $(O_FILES) $(PROGS) $(SNAPSHOT)
	@rm -f $(O_FILES) $(PROGS) $(SNAPSHOT)

foo:
	@echo MODULES = $(MODULES)
	@echo H_FILES = $(H_FILES)
	@echo C_FILES = $(C_FILES)
	@echo O_FILES = $(O_FILES)

#------------------------------------------------------------------------------

snapshot.tar.gz:	$(ALL_SRC)
	@echo Making $@
	@tar --options gzip:compression-level=9 -zcf $@ $(ALL_SRC)


#------------------------------------------------------------------------------

mset_image:	$(O_FILES)
	@echo Linking $@
	@$(LINK) -o $@ $(O_FILES)

.c.o:
	@echo Compiling $<
	@$(COMPILE) -c $< -o $@

#------------------------------------------------------------------------------

MPReal.o:       MPReal.h       MPReal.c

Common.o:	Common.h       Common.c
Memory.o:	Memory.h       Memory.c        Common.o

DeepReal.o:	DeepReal.h     DeepReal.c      Common.o
Mandelbrot.o:	Mandelbrot.h   Mandelbrot.c

RGB.o:		RGB.h          RGB.c           Common.o
Palette.o:	Palette.h      Palette.c       RGB.o
Pixel.o:	Pixel.h        Pixel.c         RGB.o
Image.o:	Image.h        Image.c         Mandelbrot.o Pixel.o Palette.o

MainControl.o:	MainControl.h  MainControl.c   Image.o

