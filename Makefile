ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)PreSelDATA2017.o: $(INCLUDEDIR)src/PreSelDATA2017.C \
		$(INCLUDEDIR)src/HLTapply.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)PreSelDATA2017.o $<
$(OUTLIB)HLTapply.o: $(INCLUDEDIR)src/HLTapply.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HLTapply.o $<
#$(OUTLIB)TaPJpsiSelectionNaod.o: $(INCLUDEDIR)/src/TaPJpsiSelectionNaod.cc
#		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TaPJpsiSelectionNaod.o $<
#$(OUTLIB)TripletSelection.o: $(INCLUDEDIR)/src/TripletSelection.cc
#		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TripletSelection.o $<

# ==================== BParkApp =============================================
X3872App:$(INCLUDEDIR)/src/X3872Application.cc\
			$(OUTLIB)PreSelDATA2017.o \
			$(OUTLIB)HLTapply.o 
			$(CXX) $(CXXFLAGS) -ldl -o X3872App $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
X3872App.clean:
			 rm -f X3872App 

# ==================== reduced trees =============================================

clean:
		rm -f $(OUTLIB)*.o
		rm -f X3872App 

all:  X3872App


