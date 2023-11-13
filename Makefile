#------------------------------------------------------------------------------
#  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
#------------------------------------------------------------------------------

# FIXME: IMPORTANT:
# Something is wrong with parsing of command-line parameters when MPFR is
# enabled.  For example, y=.1 and y=.1000 result in different y values.
USE_MPFR =	0
#               ^ Beware of setting this to 1.

COMPILE =	cc -Wall -std=c11 -O3
LINK =		cc

ifeq ($(USE_MPFR), 1)
	COMPILE += -I/opt/local/include -DUSE_MPFR
	LINK +=	-L/opt/local/lib -lmpfr
endif

#------------------------------------------------------------------------------

PROGS =		mset_image

MODULES =	Common \
		Vector3 \
                MPReal \
		Memory \
		Mandelbrot \
		RGB \
		Palette \
		Pixel \
                Camera \
		Image \
		MainControl

H_FILES =	$(MODULES:%=%.h)

C_FILES =	$(MODULES:%=%.c)

O_FILES =	$(MODULES:%=%.o)

LIBS =

ALL_SRC =	$(H_FILES) \
		$(C_FILES) \
		mset_make_map_tile \
		mset_make_zoom \
		mset_make_flyover \
		mset_blend_frames \
                mset_collate_stats \
                mset_anim \
                mset_anim_Math.pm \
                mset_anim_Vector2D.pm \
                mset_anim_CubicBezier2D.pm \
                mset_anim_Frame.pm \
                mset_anim_Animation.pm \
		Makefile

SNAPSHOT =	snapshot.tar.gz

#------------------------------------------------------------------------------

default:	all

all:		$(PROGS) $(SNAPSHOT)

distclean:	clean

clean:
	@echo Removing $(O_FILES) $(PROGS) $(SNAPSHOT)
	@rm -f $(O_FILES) $(PROGS) $(SNAPSHOT)

info:
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
	@sleep 1  # So that target has most recent timestamp.
	@$(LINK) -o $@ $(O_FILES)

.c.o:
	@echo Compiling $<
	@$(COMPILE) -c $< -o $@

#------------------------------------------------------------------------------

MPReal.o:       MPReal.h       MPReal.c

Common.o:	Common.h       Common.c
Vector3.o:	Vector3.h      Vector3.c       Common.o
Memory.o:	Memory.h       Memory.c        Common.o

Mandelbrot.o:	Mandelbrot.h   Mandelbrot.c

RGB.o:		RGB.h          RGB.c           Common.o
Palette.o:	Palette.h      Palette.c       RGB.o
Pixel.o:	Pixel.h        Pixel.c         RGB.o
Camera.o:       Camera.h       Camera.c        Common.o
Image.o:	Image.h        Image.c         Mandelbrot.o Pixel.o Palette.o

MainControl.o:	MainControl.h  MainControl.c   Image.o

