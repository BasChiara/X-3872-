#ifndef SGNvsBKGvariables_h
#define SGNvsBKGvariables_h

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

class SGNvsBKGvariables : public B0toX3872K0s{

	public:
		//CONSTRUCTOR - DECONSTRUCTOR
		SGNvsBKGvariables(TTree *tree=0);
		virtual ~SGNvsBKGvariables();
		// Loop() from B0toX3872K0s
		void Loop();

		//new functions
		int ReadTree();
		int RecoPartFillP4(const int Bidx);
		float DeltaPhi_RestFrameB0(const int Bidx);


	private:
		// outputs
		TFile* outFileHisto_;
		TString inFilePath_;
		TFile* inFileTree_;
		TTree* inTree_;

		// PDG particle ID
		int isMum = 13, isPip = 211, isJPsi = 443, isRho = 113, isK0s = 310, isX3872 = 9920443, isB0 = 511;
		// PDG particle mass
		const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648; 


		// MC truth matching ptl-idx
		int MCmatch_Mum_Idx, MCmatch_Mup_Idx;
		int MCmatch_Pim_Idx, MCmatch_Pip_Idx;
		int MCmatch_K0s_Idx;
		int MCmatch_B0_Idx;
		
		// RECO tracks P4
		ROOT::Math::PtEtaPhiMVector P4_Reco_Mu1, P4_Reco_Mu2; 
		ROOT::Math::PtEtaPhiMVector P4_Reco_JPsi;
		ROOT::Math::PtEtaPhiMVector P4_Reco_Pi1, P4_Reco_Pi2;
		ROOT::Math::PtEtaPhiMVector P4_Reco_Rho;
		ROOT::Math::PtEtaPhiMVector P4_Reco_X3872;
		ROOT::Math::PtEtaPhiMVector P4_Reco_K0trk1, P4_Reco_K0trk2;
		ROOT::Math::PtEtaPhiMVector P4_Reco_K0s;
		ROOT::Math::PtEtaPhiMVector P4_Reco_B0;
		// B0 analysis result (TTree branches)
		Long64_t TriggerSel_event;
		UInt_t nBKG_B0 ;
		Int_t B0_SGN_idx, B0_BKG_idx[50];
		Bool_t B0_BKG_isTrueJPsi[50],B0_BKG_isTruePi1[50], B0_BKG_isTruePi2[50], B0_BKG_isTrueRho[50], B0_BKG_isTrueK0s[50];





}; //SGNvsBKGvariables

#endif 
