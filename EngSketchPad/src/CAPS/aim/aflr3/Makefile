#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
ifndef ESP_ARCH
$(error ESP_ARCH must be set -- Please fix the environment...)
endif
#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

ifndef AFLR
lint default:
	@echo "Not compiling aflr3AIM: AFLR must be set -- Please fix the environment...."
else
ifndef AFLR_ARCH
lint default:
	@echo "Not compiling aflr3AIM: AFLR_ARCH must be set -- Please fix the environment...."
else
default:	$(LDIR)/aflr3AIM.so

lint: aflrlint
endif
endif

AFLR3_LIB  = -laflr3 -laflr43 -laflr4 -lsurf_auto_cpp -lsurf_auto -lrec3
AFLR3_LIB += -lug_mp -lanbl3 -laflr2c -lrec2 -lbl1
AFLR3_LIB += -lotb -legen -lice2 -lice3  -ldftr2 -ldftr3 -ldgeom
AFLR3_LIB += -lug_io -lug_gq -lug_cpp -lcgns -lqtb -lug3 -lug -lug2

ifeq ("$(ESP_ARCH)","LINUX64")
AFLR3_LIB += -lrt
endif

# CAPS compiles egads_aflr4 library to use the latest egads.h file
VPATH             = $(AFLR)/src/egads_aflr4
EGADS_INC_FILES   = $(wildcard $(IDIR)/egads*.h)
# EGADS_AFLR4_FILES = $(wildcard $(AFLR)/src/egads_aflr4/*.c)
# EGADS_AFLR4_OBJS  = $(subst $(AFLR)/src/egads_aflr4,$(ODIR),$(EGADS_AFLR4_FILES:.c=.o))

EGADS_AFLR4_FILES =               \
EG_getEffBodyTopos.c              \
egads_aflr4_check_tess_coord.c    \
egads_aflr4_check_tess_ordering.c \
egads_aflr4_get_bedge_coord.c     \
egads_aflr4_get_edge_spacing.c    \
egads_aflr4_get_faceOrder.c       \
egads_aflr4_get_face_spacing.c    \
egads_aflr4_get_surf_grid.c       \
egads_aflr4_create_tess.c         \
egads_aflr4_tess_to_dgeom.c       \
egads_aflr4_version.c             \
egads_auto_cad_geom_setup.c       \
egads_cad_geom_add_attr.c         \
egads_cad_geom_add_ff.c           \
egads_cad_geom_add_ret_attr.c     \
egads_cad_geom_bodyFlag.c         \
egads_cad_geom_check_bc.c         \
egads_cad_geom_data_cleanup.c     \
egads_cad_geom_edgeData.c         \
egads_cad_geom_edgeFace.c         \
egads_cad_geom_edit.c             \
egads_cad_geom_faceEdge.c         \
egads_cad_geom_faceFace.c         \
egads_cad_geom_faceFlag.c         \
egads_cad_geom_faceLoop.c         \
egads_cad_geom_file_read.c        \
egads_cad_geom_file_write.c       \
egads_cad_geom_idef2idef.c        \
egads_cad_geom_loopEdge.c         \
egads_cad_geom_nodeData.c         \
egads_cad_geom_reset_attr.c       \
egads_cad_geom_setup.c            \
egads_cad_geom_write_topo_data.c  \
egads_eval_curv_at_uv.c           \
egads_eval_edge_arclen.c          \
egads_eval_norm_at_uv.c           \
egads_eval_uv_at_t.c              \
egads_eval_uv_at_xyz.c            \
egads_eval_uv_bounds.c            \
egads_eval_xyz_at_u.c             \
egads_eval_xyz_at_uv.c            \
egads_face_node_map.c             \
egads_cad_geom_faceMatch.c        \
egads_data.c                      \
egads_mzone.c

EGADS_AFLR4_OBJS  = $(addprefix $(ODIR)/,$(EGADS_AFLR4_FILES:.c=.o))

OBJS=$(EGADS_AFLR4_OBJS)

###############################################################################
# Used to add debug information to AFLR3
# Please do not remove this - MCG
#AFLR3_FILES=$(wildcard $(AFLR)/src/aflr3/*.c)
#AFLR3_OBJS=$(patsubst %.c,%.o,$(AFLR3_FILES))

#VPATH      += $(AFLR)/src/ug
#UG_FILES=$(wildcard $(AFLR)/src/ug/*.c)
#UG_OBJS=$(patsubst %.c,%.o,$(UG_FILES))

#VPATH      += $(AFLR)/src/ug3
#UG3_FILES=$(wildcard $(AFLR)/src/ug3/*.c)
#UG3_OBJS=$(patsubst %.c,%.o,$(UG3_FILES))

#VPATH      += $(AFLR)/src/ug_io
#UG_IO_FILES = $(wildcard $(AFLR)/src/ug_io/*.c)
#UG_IO_OBJS  = $(subst $(AFLR)/src/ug_io,$(ODIR),$(UG_IO_FILES:.c=.o))

#VPATH      += $(AFLR)/src/ug_cpp
#UG_CPP_FILES = $(wildcard $(AFLR)/src/ug_cpp/*.cpp)
#UG_CPP_OBJS  = $(subst $(AFLR)/src/ug_cpp,$(ODIR),$(UG_CPP_FILES:.cpp=.o))

#VPATH    += $(AFLR)/src/otb
#OTB_FILES = $(wildcard $(AFLR)/src/otb/*.c)
#OTB_OBJS  = $(subst $(AFLR)/src/otb,$(ODIR),$(OTB_FILES:.c=.o))

#VPATH      += $(AFLR)/src/dgeom
#DGEOM_FILES = $(wildcard $(AFLR)/src/dgeom/*.c)
#DGEOM_OBJS  = $(subst $(AFLR)/src/dgeom,$(ODIR),$(DGEOM_FILES:.c=.o))

#VPATH      += $(AFLR)/src/aflr4
#AFLR4_FILES = $(wildcard $(AFLR)/src/aflr4/*.c)
#AFLR4_OBJS  = $(subst $(AFLR)/src/aflr4,$(ODIR),$(AFLR4_FILES:.c=.o))

#VPATH       += $(AFLR)/src/aflr43
#AFLR43_FILES = $(wildcard $(AFLR)/src/aflr43/*.c)
#AFLR43_OBJS  = $(subst $(AFLR)/src/aflr43,$(ODIR),$(AFLR43_FILES:.c=.o))

#VPATH      += $(AFLR)/src/anbl3
#ANBL3_FILES = $(wildcard $(AFLR)/src/anbl3/*.c)
#ANBL3_OBJS  = $(subst $(AFLR)/src/anbl3,$(ODIR),$(ANBL3_FILES:.c=.o))

OBJS+=$(AFLR3_OBJS) $(UG_OBJS) $(UG3_OBJS) $(UG_IO_OBJS) $(AFLR4_OBJS) $(AFLR43_OBJS) $(DGEOM_OBJS) $(OTB_OBJS) $(ANBL3_OBJS)
#OBJSP=$(UG_CPP_OBJS)
###############################################################################

$(LDIR)/aflr3AIM.so:	$(ODIR)/aflr3AIM.o $(ODIR)/aflr3_Interface.o \
			$(LDIR)/libutils.a  $(LDIR)/libaimUtil.a $(OBJS) $(OBJSP)
	touch $(LDIR)/aflr3AIM.so
	rm $(LDIR)/aflr3AIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/aflr3AIM.so $(ODIR)/aflr3AIM.o \
		$(ODIR)/aflr3_Interface.o $(OBJS) $(OBJSP) -L$(AFLR)/$(AFLR_ARCH)/lib \
		$(AFLR3_LIB) -L$(LDIR) -lutils \
		-laimUtil -locsm -legads -ludunits2 -ldl $(RPATH) $(CPPSLB) -lm

$(ODIR)/aflr3AIM.o:	$(ODIR)/aflr3_version aflr3AIM.c aflr3_Interface.h $(IDIR)/capsTypes.h \
			$(IDIR)/aimUtil.h $(IDIR)/aimMesh.h ../utils/meshUtils.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils -I$(AFLR)/src \
		aflr3AIM.c -o $(ODIR)/aflr3AIM.o

$(ODIR)/aflr3_Interface.o:	$(ODIR)/aflr3_version aflr3_Interface.c aflr3_Interface.h \
				../utils/meshUtils.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I../utils -I$(AFLR)/src \
		aflr3_Interface.c -o $(ODIR)/aflr3_Interface.o

# Compile and execute to check the version number
.DELETE_ON_ERROR:
$(ODIR)/aflr3_version:	aflr3_version.c $(AFLR)/src/aflr3/aflr3_version.h 
	$(CC) $(COPTS) $(DEFINE) -I$(IDIR) -I$(AFLR)/src aflr3_version.c \
		-L$(AFLR)/$(AFLR_ARCH)/lib -laflr3 \
		-o $(ODIR)/aflr3_version
	$(ODIR)/aflr3_version

# Build for each AFLR C source file
$(ODIR)/%.o: %.c $(EGADS_INC_FILES) | $(ODIR)/aflr3_version
	$(CC) -c $(COPTS) $(DEFINE) -D_ENABLE_EGADS_ -w -I$(IDIR) -I$(AFLR)/src $< -o $@

ifdef OBJSP
# Build for each AFLR C++ source file
$(ODIR)/%.o: %.cpp $(ODIR)/aflr3_version
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(AFLR)/src $< -o $@
endif

aflrlint:
	$(LINT) -I../../include -I$(IDIR) -I../utils -I$(AFLR)/src -I$(AFLR)/src/aflr3 $(wildcard *.c) -nullpass -exportlocal -bufferoverflowhigh

clean:
	-rm -f $(ODIR)/aflr3_Interface.o $(ODIR)/aflr3AIM.o $(ODIR)/aflr3_version $(OBJS)

cleanall:	clean
	-rm -f $(LDIR)/aflr3AIM.so

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen aflr3AIM_dox.cfg; cd latex; make; mv refman.pdf aflr3AIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f aflr3AIM.tag)
