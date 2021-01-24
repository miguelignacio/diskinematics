# A prototype Makefile for H1
#
# last update: $Date: 2006/10/21 13:40:57 $
# by         : $Author: ozerov $
#
# The name of this module
MODULENAME = EventShapesAnalysis
#
# macro definitions
# DEFINES = 
#
# debug on or off (comment next line to disable debugger information)
DEBUG = -g -Wall
#
# The shared library to make (if any): e.g. libH1Tracks.so
# SLIB_TARGETS = 
SLIB = libEventShapes.so
#
# The static library to make (if any): e.g. libH1Tracks.a
# LIB_TARGETS = 
LIB = libEventShapes.a
#
# The binaries you would like to make, you have to provide the final linking
BIN = EventShapes
#clusters_ods dstar_mods jpsi_mods kaonfind_ods read_ods read_usertree \
rerun_finder steermanager empz_hat l1te_hat  \
lumicalc \
copyMyEvents write_eventlist create_eventlist read_eventlist dqhat_access \
test_binning cintsteering dbaccess

#ifeq ($(STANDALONE),yes)
#  # no rerun_rec
#else
#  BIN += rerun_rec 
#endif

#
# Objects you need for the libraries
OBJS = AnalysisBase.o
OBJS += AnalysisEventShapes.o
OBJS += H2020HistManager.o
OBJS += FidVolCut.o
OBJS += TDetectQedc.o
OBJS += elecCut.o
OBJS += JetTools.o
OBJS += EventshapeTools.o
OBJS += LumiFillList.o
OBJS += SpacLinearity.o
OBJS += Alignment.o
#OBJS += H1Boostedjets.o
OBJS += JetReweight.o
OBJS += Interpolate2d.o
OBJS += EventshapeTools.o

# OBJS += UserEvent.o UserFiller.o

# Dictionary objects, default is: DICTOBJS = $(OBJS:.o=Dict.o)
DICTOBJS = $(OBJS:.o=Dict.o) 

# Now your special libs and includes
# EXTRA_INCLUDEDIRS = ~/my_favourite_library/include
# EXTRA_LIBDIRS = ~/my_favourite_library/lib
# EXTRA_LINKFLAGS = -lmy_favourite_library
# Dummy objects
# DUMMYOBJS =

# Dummy dictionary objects
# DUMMYDICTOBJS = $(DUMMYOBJS:.o=Dict.o)

CXX_FLAGS= -g -Wall -fPIC -std=gnu++17
#extra files to be removed when calling 'make clean' or 'make clobber'
EXTRA_CLEAN_RULES=
EXTRA_DISTCLEAN_RULES=


#extra files for fjcontrib for H1Example 'centauro'
#FASTJET_VERS = $(FASTJET__HOME)/lib 
ifdef FJCONTRIB__HOME
  EXTRA_INCLUDEDIRS = $(FASTJET__HOME)/include $(FJCONTRIB__HOME)/include
endif


# defines LIBS
# defines Rules for .c .C .cpp .f Dict.C .so
ifeq ($(H1BASICMK),local)
  include  ../oo_tools/H1Basic.mk
else
  include  $(H1DIST)/releases/$(H1CURRENT)/oo_tools/H1Basic.mk
endif                                                             

# include the dependencies if created via 'make depend'
-include Makefile.dep

# the executables

EventShapes: outputdirs $(BUILDDIR)/main.o
	@echo ""
	@echo Linking programm $@: $@ ...
	$(LOADER) \
		$(BUILDDIR)/main.o \
		\
		$(OO_LIBS) -lEventShapes \
		\
		-lH1Geom -lH1Ods -lH1OOBanks -lH1Tracks -lH1Clusters -lH1Mods -lH1Hat -lH1Skeleton -lH1Steering \
		-lH1Pointers -lH1Arrays -lH1Tools -lH1Benchmarks -lH1Arrays -lH1Binning \
		-lH1PhysUtils  -lH1JetFinder -lH1Finder -lH1HadronicCalibration -lH1CaloTrigger \
		-lH1Calculator -lH1CalcPointers -lH1CalcWeights -lH1SoftLeptonId  -lH1Analysis \
		-lH1MagfieldOO -lH1Cuts	-lH1TrkFinder -lH1Selection -lH1ElecCalibration -lH1QCDFunc  -l${LHAPDFLIB} \
		\
		-L$(H1_LIBS) -L$(CERN_LIBS) \
		$(ROOTLIBS) -lEG -l${FORTRANLIB} \
		-L$(FJCONTRIB__HOME)/lib -lCentauro \
		-L$(FASTJET__HOME)/lib -lfastjet \
		\
		-o $(BINDIR)/$@
	@echo ""
	@echo "copying .pcm files to libdir"
	cp *.pcm $(LIBDIR)
	@echo ""
	@echo $@ created ...
	@echo Please run: $(BINDIR)/$@
	@echo ""

#
# to remove the generated *Ptr.[hC] files:
#clean::
#	rm  *Ptr.[Ch] 2>/dev/null
#	-cd $(BUILDDIR) &&\
#        rm Makefile.dep *.o *~ core *.bak 2>/dev/null

clean::
	rm ../bin/$(BUILDDIR)/EventShapes
	@echo "deleteing pcm-files from (local) libdir"
	rm $(LIBDIR)/*.pcm
