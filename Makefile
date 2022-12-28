ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA -lRooFit -lRooFitCore -lMinuit

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
$(OUTLIB)OptimizerMVA.o: $(INCLUDEDIR)src/OptimizerMVA.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)OptimizerMVA.o $<
$(OUTLIB)CutterMVA.o: $(INCLUDEDIR)src/CutterMVA.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutterMVA.o $<
$(OUTLIB)BDTapply.o: $(INCLUDEDIR)src/BDTapply.C
		$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BDTapply.o $<

# ==================== X3872 App =============================================
X3872App:$(INCLUDEDIR)/src/X3872Application.cc\
			$(OUTLIB)PreSelDATA2017.o \
			$(OUTLIB)HLTapply.o 
			$(CXX) $(CXXFLAGS) -ldl -o X3872App $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
X3872App.clean:
			 rm -f X3872App 

# ==================== MVA optimization=============================================
MVAoptim:$(INCLUDEDIR)src/OptimizeMVA.cc\
			$(OUTLIB)OptimizerMVA.o
			$(CXX) $(CXXFLAGS) -ldl -o OptimizeMVA $(OUTLIB)/OptimizerMVA.o $(GLIBS) $(LDFLAGS) $ $<
MVAoptim.clean:
				rm -f OptimizeMVA 

# ==================== MVA-cut application =============================================
CUTapply:$(INCLUDEDIR)src/ApplyMVAcuts.cc\
			$(OUTLIB)OptimizerMVA.o
			$(CXX) $(CXXFLAGS) -ldl -o ApplyMVAcuts $(OUTLIB)/OptimizerMVA.o $(GLIBS) $(LDFLAGS) $ $<
CUTapply.clean:
				rm -f ApplyMVAcuts

# ==================== MVA-cut unblinded =============================================
CUTunbl: $(INCLUDEDIR)src/ApplyMVAcutsUNBL.cc\
			$(OUTLIB)CutterMVA.o\
			$(OUTLIB)BDTapply.o
			$(CXX) $(CXXFLAGS) -ldl -o ApplyMVAcutsUNBL $(OUTLIB)/CutterMVA.o $(OUTLIB)/BDTapply.o $(GLIBS) $(LDFLAGS) $ $<
CUTunbl.clean:
				rm -f ApplyMVAcutsUNBL

# ==================== MVA application =============================================
MVAnalysis:$(INCLUDEDIR)src/MVAnalysis.C
			  $(CXX) $(CXXFLAGS) -ldl -o MVAnalysis $(GLIBS) $(LDFLAGS) $ $<
MVAnalysis.clean:
				rm -f MVAnalysis


# ==================== reduced trees =============================================

clean:
		rm -f $(OUTLIB)*.o
		rm -f X3872App 
		rm -f MVAnalysis
		rm -f OptimizeMVA
		rm -f ApplyMVAcuts
		rm -f ApplyMVAcutsUNBL

mva: MVAnalysis

opt: MVAoptim

cut: CUTapply

unb: CUTunbl

hlt: X3872App
