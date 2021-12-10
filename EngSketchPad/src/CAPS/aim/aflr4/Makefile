#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
ifndef AFLR
$(error AFLR must be set -- Please fix the environment...)
endif
#
ifndef AFLR_ARCH
$(error AFLR_ARCH must be set -- Please fix the environment...)
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

AFLR4_LIB  = -laflr4 -laflr2c
AFLR4_LIB += -lsurf_auto -lsurf_auto_cpp -laflr4 -laflr3
AFLR4_LIB += -legen -lotb -lrec2 -lrec3 -ldftr3 -ldftr2 -ldgeom -lqtb
AFLR4_LIB += -lice2 -lice3 -lug3 -lug2 -lug -lug_mp

ifeq ("$(ESP_ARCH)","LINUX64")
AFLR4_LIB += -lrt
endif

# CAPS compiles egads_aflr4 library to use the latest egads.h file
VPATH             = $(AFLR)/src/egads_aflr4
EGADS_INC_FILES   = $(wildcard $(IDIR)/egads*.h)
# EGADS_AFLR4_FILES = $(wildcard $(AFLR)/src/egads_aflr4/*.c)
# EGADS_AFLR4_OBJS  = $(subst $(AFLR)/src/egads_aflr4,$(ODIR),$(EGADS_AFLR4_FILES:.c=.o))

EGADS_AFLR4_FILES = \
EG_getEffBodyTopos.c \
egads_aflr4_check_tess_coord.c \
egads_aflr4_check_tess_ordering.c \
egads_aflr4_get_bedge_coord.c \
egads_aflr4_get_edge_spacing.c \
egads_aflr4_get_face_spacing.c \
egads_aflr4_get_surf_grid.c \
egads_aflr4_get_tess.c \
egads_aflr4_version.c \
egads_auto_cad_geom_setup.c \
egads_cad_geom_add_attr.c \
egads_cad_geom_add_ff.c \
egads_cad_geom_add_ret_attr.c \
egads_cad_geom_bodyFlag.c \
egads_cad_geom_check_bc.c \
egads_cad_geom_data_cleanup.c \
egads_cad_geom_edgeData.c \
egads_cad_geom_edgeFace.c \
egads_cad_geom_edit.c \
egads_cad_geom_faceEdge.c \
egads_cad_geom_faceFace.c \
egads_cad_geom_faceFlag.c \
egads_cad_geom_faceLoop.c \
egads_cad_geom_file_read.c \
egads_cad_geom_file_write.c \
egads_cad_geom_idef2idef.c \
egads_cad_geom_loopEdge.c \
egads_cad_geom_nodeData.c \
egads_cad_geom_reset_attr.c \
egads_cad_geom_setup.c \
egads_cad_geom_write_topo_data.c \
egads_eval_arclen.c \
egads_eval_curv_at_uv.c \
egads_eval_edge_uv.c \
egads_eval_norm_at_uv.c \
egads_eval_uv_at_xyz.c \
egads_eval_uv_bounds.c \
egads_eval_xyz_at_u.c \
egads_eval_xyz_at_uv.c \
egads_set_ext_cad_data.c

EGADS_AFLR4_OBJS  = $(addprefix $(ODIR)/,$(EGADS_AFLR4_FILES:.c=.o))

###############################################################################
# Used to add debug information to AFLR4
# Please do not remove this - MCG
#VPATH      += $(AFLR)/src/aflr4
#AFLR4_FILES = $(wildcard $(AFLR)/src/aflr4/*.c)
#AFLR4_OBJS  = $(subst $(AFLR)/src/aflr4,$(ODIR),$(AFLR4_FILES:.c=.o))

#VPATH          += $(AFLR)/src/surf_auto
#SURF_AUTO_FILES = $(wildcard $(AFLR)/src/surf_auto/*.c)
#SURF_AUTO_OBJS  = $(subst $(AFLR)/src/surf_auto,$(ODIR),$(SURF_AUTO_FILES:.c=.o))

#VPATH    += $(AFLR)/src/ug3
#UG3_FILES = $(wildcard $(AFLR)/src/ug3/*.c)
#UG3_OBJS  = $(subst $(AFLR)/src/ug3,$(ODIR),$(UG3_FILES:.c=.o))

#VPATH      += $(AFLR)/src/dgeom
#DGEOM_FILES = $(wildcard $(AFLR)/src/dgeom/*.c)
#DGEOM_OBJS  = $(subst $(AFLR)/src/dgeom,$(ODIR),$(DGEOM_FILES:.c=.o))

#VPATH              += $(AFLR)/src/surf_auto_cpp
#SURF_AUTO_CPP_FILES = $(wildcard $(AFLR)/src/surf_auto_cpp/*.cpp)
#SURF_AUTO_CPP_OBJS  = $(subst $(AFLR)/src/surf_auto_cpp,$(ODIR),$(SURF_AUTO_CPP_FILES:.cpp=.o))

OBJS=$(EGADS_AFLR4_OBJS) $(AFLR4_OBJS) $(SURF_AUTO_OBJS) $(UG3_OBJS) $(DGEOM_OBJS)
#OBJSP=$(SURF_AUTO_CPP_OBJS)
###############################################################################

$(LDIR)/aflr4AIM.so:	$(ODIR)/aflr4_version $(ODIR)/aflr4AIM.o $(ODIR)/aflr4_Interface.o \
			$(LDIR)/libutils.a $(LDIR)/libaimUtil.a $(OBJS) $(OBJSP)
	touch $(LDIR)/aflr4AIM.so
	rm $(LDIR)/aflr4AIM.so
	$(CC) $(SOFLGS) -o $(LDIR)/aflr4AIM.so $(ODIR)/aflr4AIM.o \
		$(ODIR)/aflr4_Interface.o $(OBJS) $(OBJSP) -L$(AFLR)/$(AFLR_ARCH)/lib \
		$(AFLR4_LIB) -L$(LDIR) -lutils \
		-laimUtil -locsm -legads -ludunits2 -ldl $(RPATH) $(CPPSLB) -lm

$(ODIR)/aflr4AIM.o:	$(ODIR)/aflr4_version aflr4AIM.c $(IDIR)/capsTypes.h \
			../utils/meshUtils.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(AFLR)/src -I../utils aflr4AIM.c \
		-o $(ODIR)/aflr4AIM.o

$(ODIR)/aflr4_Interface.o:	$(ODIR)/aflr4_version aflr4_Interface.c ../utils/meshUtils.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) -I$(AFLR)/src \
		-I../utils aflr4_Interface.c -o $(ODIR)/aflr4_Interface.o

# Compile and execute to check the version number
.DELETE_ON_ERROR:
$(ODIR)/aflr4_version:	aflr4_version.c $(AFLR)/src/aflr4/aflr4_version.h 
	$(CC) $(COPTS) $(DEFINE) -I$(IDIR) -I$(AFLR)/src aflr4_version.c \
		-L$(AFLR)/$(AFLR_ARCH)/lib -laflr4 -lm \
		-o $(ODIR)/aflr4_version
	$(ODIR)/aflr4_version

# Build for each AFLR C source file
$(ODIR)/%.o: %.c $(ODIR)/aflr4_version $(EGADS_INC_FILES)
	$(CC) -c $(COPTS) $(DEFINE) -w -I$(IDIR) -I$(AFLR)/src $< -o $@

ifdef OBJSP
# Build for each AFLR C++ source file
$(ODIR)/%.o: %.cpp $(ODIR)/aflr4_version
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(AFLR)/src $< -o $@
endif

lint:
	$(LINT) -I../../include -I$(IDIR) $(wildcard *.c) -I../utils -I$(AFLR)/src -exportlocal

clean:
	-rm -f $(ODIR)/aflr4_Interface.o $(ODIR)/aflr4AIM.o $(ODIR)/aflr4_version $(OBJS) $(OBJSP)

cleanall:	clean
	-rm -f $(LDIR)/aflr4AIM.so

dox: doc/html/index.html
doc/html/index.html: $(wildcard *.c) $(wildcard ../utils/*.c) $(wildcard doc/*.cfg) $(wildcard doc/*.tex) $(wildcard doc/*.xml) $(wildcard doc/*.dox)
	(cd doc; doxygen aflr4AIM_dox.cfg; cd latex; make; mv refman.pdf aflr4AIM.pdf)

doxclean:
	(cd doc; rm -f INPUT; rm -rf html latex; rm -f aflr4AIM.tag)
