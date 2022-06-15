#ifndef TriggerSelection_h
#define TriggerSelection_h

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
#include <fstream>
#include <stddef.h>
#include <vector>
#include <string>

#include <TMath.h>
#include "Math/Vector4D.h"
#include <Math/GenVector/VectorUtil.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "TLorentzVector.h"

class TriggerSelection : public B0toX3872K0s{

	public:
		//CONSTRUCTOR - DECONSTRUCTOR
		TriggerSelection(TTree *tree=0);
		virtual ~TriggerSelection();
		// Loop() from B0toX3872K0s
		void Loop();

		//new functions
		void OutTree_setup();
		void GetFitParams();
		void GenPartFillP4();
		int  RecoPartFillP4(const int Bidx);

		void  MCtruthMatching(const bool verbose = false);
		float DeltaPT(ROOT::Math::PtEtaPhiMVector genV, ROOT::Math::PtEtaPhiMVector recV);
		bool  isMCmatchingFailed();
		void  WhichPtl_MCmissed(TH1* histo);
		void  B0cand_PTLmissed(const bool isMCmatchedJPsi, const bool isMCmatchedRho, const bool isMCmatchedK0s, TH1* histo);

		int  ApplyTriggerSelection_Muons(const int Bidx);
		int  ApplyTriggerSelection_Track(const int Bidx);
		
		bool inBKGregion(const TString ptl);
		void NK0s_per_Rho(TH1* histo);
		void NRho_per_K0s(TH1* histo);

		void PrintGenLevel(const int Bidx, const int ptlID, const bool isMCmatched_Pi1, const bool isMCmatched_Pi2);
		void DebugRho();

	private:
		// outputs
		std::string FileFitParB0, FileFitParK0s, FileFitParX;
		TFile* outFileHisto_;
		TFile* outFileTree_;
		TTree* outTree_;

		// PDG particle ID
		int isMum = 13, isPip = 211, isJPsi = 443, isRho = 113, isK0s = 310, isX3872 = 9920443, isB0 = 511;
		// PDG particle mass
		const float mMuon = 0.105658, mPion = 0.1395704, mK0s = 0.497648; 

		// generator ptls 4-vectors
		ROOT::Math::PtEtaPhiMVector GenMumP4, GenMupP4; // muons
		ROOT::Math::PtEtaPhiMVector GenPimP4, GenPipP4; // pions
		ROOT::Math::PtEtaPhiMVector GenK0sP4;// K0s
		ROOT::Math::PtEtaPhiMVector GenRhoP4;// Rho 

		// MC truth matching ptl-idx
		int MCmatch_Mum_Idx, MCmatch_Mup_Idx;
		int MCmatch_Pim_Idx, MCmatch_Pip_Idx;
		int MCmatch_K0s_Idx;
		int MCmatch_B0_Idx;
			// ... DRmin & DPt
		float MCmatch_Mum_DRmin, MCmatch_Mup_DRmin;
		float MCmatch_Mum_DpT, MCmatch_Mup_DpT;
		float MCmatch_Pim_DRmin, MCmatch_Pip_DRmin;
		float MCmatch_Pim_DpT, MCmatch_Pip_DpT;
		float MCmatch_K0s_DRmin;
		float MCmatch_K0s_DpT;
		float MCmatch_B0_DRmin;
		float MCmatch_B0_DpT;
		
		// RECO tracks P4
		ROOT::Math::PtEtaPhiMVector P4_Reco_Mu1, P4_Reco_Mu2; 
		ROOT::Math::PtEtaPhiMVector P4_Reco_Pi1, P4_Reco_Pi2;
		ROOT::Math::PtEtaPhiMVector P4_Reco_K0s;
		// X & K0s Signal Region
		double MX_nearLeft, MX_nearRight;
		double MK0s_nearLeft, MK0s_nearRight;

		// B0 analysis result (TTree branches)
		Long64_t TriggerSel_event;
		UInt_t nBKG_B0 ;
		Int_t B0_SGN_idx, B0_BKG_idx[50];
		Bool_t B0_BKG_isTrueJPsi[50], B0_BKG_isTruePi1[50], B0_BKG_isTruePi2[50], B0_BKG_isTrueRho[50], B0_BKG_isTrueK0s[50];





}; //TriggerSelection

#endif 
