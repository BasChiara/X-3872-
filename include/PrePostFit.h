#ifndef PrePostFit_h
#define PrePostFit_h

//#include "../src/B0toX3872K0s.C"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include "TTree.h"  
#include "TCanvas.h"  
#include "TStyle.h"  
#include <iostream>   
#include <vector>
#include <string>

#include <TMath.h>
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TLorentzVector.h"

class PrePostFit : public B0toX3872K0s{

  public:
    //CONSTRUCTOR - DECONSTRUCTOR
    PrePostFit(TTree *tree=0);
    virtual ~PrePostFit();
    // Loop() from B0toX3872K0s
    void Loop();

    //new functions
    void GenPartFillP4();
    void MCtruthMatching(const bool verbose = false);
    void ApplyTriggerSelection(const int Bidx, bool*isOkMuMu_L0, bool* isOkMuMu_L1, bool* isOkMuMu_L2, bool* isOkTrk_L0, bool*isOkTrk_L1);

    float DeltaPT(ROOT::Math::PtEtaPhiMVector genV, ROOT::Math::PtEtaPhiMVector recV);


  private:
    // outputs
    TFile* outFile_;
    TTree* outTree_;

    // PDG particle ID
    int isMum = 13, isPip = 211, isJPsi = 443, isRho = 113, isK0s = 310, isX3872 = 9920443, isB0 = 511;
    // PDG particle mass
    const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648; 
 
    // generator ptls 4-vectors
    ROOT::Math::PtEtaPhiMVector GenMumP4, GenMupP4; // muons
    ROOT::Math::PtEtaPhiMVector GenPimP4, GenPipP4; // pions
    ROOT::Math::PtEtaPhiMVector GenK0sP4;// K0s
    ROOT::Math::PtEtaPhiMVector GenRhoP4;// K0s

    // MC truth matching ptl 
    // ... idx
    int MCmatch_Mum_Idx, MCmatch_Mup_Idx;
    int MCmatch_Pim_Idx, MCmatch_Pip_Idx;
    int MCmatch_K0s_Idx;
    // ... 4 vectors

    ROOT::Math::PtEtaPhiMVector P4_Reco_Mu1, P4_Reco_Mu2; 
    ROOT::Math::PtEtaPhiMVector P4_Reco_Pi1, P4_Reco_Pi2;
    ROOT::Math::PtEtaPhiMVector P4_Reco_K0s;

    
    

}; //PrePostFit

#endif 
