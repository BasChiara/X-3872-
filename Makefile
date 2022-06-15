ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g -lGenVector # con errori tipo 'undefined reference to `ROOT::Math::GenVector::Throw(char const*)'
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew , $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)PreSelDATA2017.o: $(INCLUDEDIR)src/PreSelDATA2017.C 
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)PreSelDATA2017.o $<
$(OUTLIB)HLTapply.o: $(INCLUDEDIR)src/HLTapply.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HLTapply.o $<

# ==================== BParkApp =============================================
X3872App:$(INCLUDEDIR)/src/X3872Application.cc\
			$(OUTLIB)PreSelDATA2017.o \
			$(OUTLIB)HLTapply.o 
			$(CXX) $(CXXFLAGS) -ldl -o X3872App $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
X3872App.clean:
			 rm -f X3872App 
MVAnalysis:$(INCLUDEDIR)src/MVAnalysis.C
			  $(CXX) $(CXXFLAGS) -ldl -o MVAnalysis $(GLIBS) $(LDFLAGS) $ $<

# ==================== reduced trees =============================================

clean:
		rm -f $(OUTLIB)*.o
		rm -f X3872App 

mva: MVAnalysis

all:  X3872App


