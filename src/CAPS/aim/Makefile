#
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
#

export AIM_UTILSDIR = $(realpath ./utils)

all:
	(cd utils; make)
	(cd avl; make)
	(cd delaundo; make)
	(cd egadsTess; make)
	(cd friction; make)
	(cd fun3d; make)
	(cd mystran; make)
	(cd masstran; make)
	(cd nastran; make)
	(cd astros; make)
	(cd hsm; make)
	(cd skeleton; make)
	(cd su2; make)
	(cd tsfoil; make)
	(cd xfoil; make)
ifdef AFLR
	(cd aflr2; make)
	(cd aflr3; make)
	(cd aflr4; make)
endif
ifdef AWAVE
	(cd awave; make)
endif
ifdef TETGEN
	(cd tetgen; make)
endif
ifneq ("$(ESP_ARCH)","CYGWIN64")
	(cd cart3d; make)
endif
	(cd pointwise; make)
	(cd meshWriter; make)
	(cd custom; make)

lint:
	(cd avl; make lint)
	(cd delaundo; make lint)
	(cd egadsTess; make lint)
	(cd friction; make lint)
	(cd fun3d; make lint)
	(cd mystran; make lint)
	(cd masstran; make lint)
	(cd nastran; make lint)
	(cd astros; make lint)
	(cd hsm; make lint)
	(cd skeleton; make lint)
	(cd su2; make lint)
	(cd tsfoil; make lint)
	(cd xfoil; make lint)
ifdef AFLR
	(cd aflr2; make lint)
	(cd aflr3; make lint)
	(cd aflr4; make lint)
endif
ifdef AWAVE
	(cd awave; make lint)
endif
ifdef TETGEN
	(cd tetgen; make lint)
endif
ifneq ("$(ESP_ARCH)","CYGWIN64")
	(cd cart3d; make lint)
endif
	(cd pointwise; make lint)
	(cd meshWriter; make lint)
	(cd utils; make lint)

clean:
	(cd utils; make clean)
	(cd avl; make clean)
	(cd delaundo; make clean)
	(cd egadsTess; make clean)
	(cd friction; make clean)
	(cd fun3d; make clean)
	(cd mystran; make clean)
	(cd masstran; make clean)
	(cd nastran; make clean)
	(cd astros; make clean)
	(cd hsm; make clean)
	(cd skeleton; make clean)
	(cd su2; make clean)
	(cd tsfoil; make clean)
	(cd xfoil; make clean)
ifdef AFLR
	(cd aflr2; make clean)
	(cd aflr3; make clean)
	(cd aflr4; make clean)
endif
ifdef AWAVE
	(cd awave; make clean)
endif
ifdef TETGEN
	(cd tetgen; make clean)
endif
ifneq ("$(ESP_ARCH)","CYGWIN64")
	(cd cart3d; make clean)
endif
	(cd pointwise; make clean)
	(cd meshWriter; make clean)
	(cd custom; make clean)

cleanall:
	(cd utils; make cleanall)
	(cd avl; make cleanall)
	(cd delaundo; make cleanall)
	(cd egadsTess; make cleanall)
	(cd friction; make cleanall)
	(cd fun3d; make cleanall)
	(cd skeleton; make cleanall)
	(cd su2; make cleanall)
	(cd mystran; make cleanall)
	(cd masstran; make cleanall)
	(cd nastran; make cleanall)
	(cd astros; make cleanall)
	(cd hsm; make cleanall)
	(cd tsfoil; make cleanall)
	(cd xfoil; make cleanall)
ifdef AFLR
	(cd aflr2; make cleanall)
	(cd aflr3; make cleanall)
	(cd aflr4; make cleanall)
endif
ifdef AWAVE
	(cd awave; make cleanall)
endif
ifdef TETGEN
	(cd tetgen; make cleanall)
endif
ifneq ("$(ESP_ARCH)","CYGWIN64")
	(cd cart3d; make cleanall)
endif
	(cd pointwise; make cleanall)
	(cd meshWriter; make cleanall)
	(cd custom; make cleanall)

dox:
	(cd avl; make dox)
	(cd awave; make dox)
	(cd friction; make dox)
	(cd delaundo; make dox)
	(cd egadsTess; make dox)
	(cd fun3d; make dox)
	(cd mystran; make dox)
	(cd masstran; make dox)
	(cd nastran; make dox)
	(cd astros; make dox)
	(cd hsm; make dox)
	(cd su2; make dox)
	(cd aflr2; make dox)
	(cd aflr3; make dox)
	(cd aflr4; make dox)
	(cd tetgen; make dox)
	(cd tsfoil; make dox)
	(cd xfoil; make dox)
	(cd cart3d; make dox)
	(cd pointwise; make dox)

doxclean:
	(cd avl; make doxclean)
	(cd awave; make doxclean)
	(cd friction; make doxclean)
	(cd delaundo; make doxclean)
	(cd egadsTess; make doxclean)
	(cd fun3d; make doxclean)
	(cd mystran; make doxclean)
	(cd masstran; make doxclean)
	(cd nastran; make doxclean)
	(cd astros; make doxclean)
	(cd hsm; make doxclean)
	(cd su2; make doxclean)
	(cd aflr2; make doxclean)
	(cd aflr3; make doxclean)
	(cd aflr4; make doxclean)
	(cd tetgen; make doxclean)
	(cd tsfoil; make doxclean)
	(cd xfoil; make doxclean)
	(cd cart3d; make doxclean)
	(cd pointwise; make doxclean)
