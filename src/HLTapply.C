#include "../include/HLTapply.h" 

// CONSTRUCTOR-DESTRUCTOR
HLTapply::HLTapply(TTree *tree, const TString dataset) : PreSelDATA2017(tree){ 

	dataset_ = dataset;
	outFileHistoPath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/plots/HLTapply_" + dataset_ + ".root";
   outFileTreePath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BKG_" + dataset_ + ".root";
	
	FileFitParB0 = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";
	FileFitParK0s= "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/K0sparams.txt";
	FileFitParX  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/X3872params.txt";
   
	P4_Reco_Mu1.SetM(mMuon); P4_Reco_Mu2.SetM(mMuon);
   P4_Reco_Pi1.SetM(mPion); P4_Reco_Pi2.SetM(mPion);
   P4_Reco_K0s.SetM(mK0s);
}

HLTapply::~HLTapply(){
}





void HLTapply::Loop() {

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Long64_t Nbreak = nentries + 10, NPrint = (int)nentries/20.; 

   // ----- VARIABLES ----- //
   int Ncandidates_JPsi = 0, Ncandidates_Rho = 0, Ncandidates_X3872 = 0, Ncandidates_K0s = 0, Ncandidates_B0 = 0;
	int TOTcand_JPsi = 0, TOTcand_B0 = 0, TOTsidebands_B0 = 0 ;
   int NTriggeredEv = 0, NTrigSelEV = 0;
	int TriggeredTrkID;

   bool ToCountJPsi = true, ToCountRho = true, ToCountX = true, ToCountK0s = true;
	ROOT::Math::PtEtaPhiMVector Prev_P4_Reco_K0s(0,0,0,0);
   bool isTriggerON, isTrkQC, isSel_MuMu, isSel_Trk, isSel_Pim, isSel_Pip;

   int SelMu1_idx = -1, SelMu2_idx = -1;
	int prevMu1_idx = -1, prevMu2_idx = -1, prevPi1_idx = -1, prevPi2_idx = -1, prevK0s_idx = -1;

   float MuMu_M, PiPi_M, MuMuPi1_M, MuMuPi2_M, MuMuPiPi_M, MuMuPiPiK0s_M;
	const float MuMuPi_Mlow = 3.1, MuMuPi_Mhigh  = 3.9;
	float Nmin = 4, Nmax = 8;
		
	// ----- FIT PARAMETER FOR SIGNAL RANGE -----//
	GetFitParams();
	const float MB_nearLeft   = meanMB0  - Nmin * sigmaMB0,  MB_nearRight   = meanMB0  + Nmin * sigmaMB0;
	const float MB_farLeft    = meanMB0  - Nmax * sigmaMB0,  MB_farRight    = meanMB0  + Nmax * sigmaMB0;  
	std::cout << "B0 interval [" << MB_farLeft << "," << MB_nearLeft << "]" << " + [" << MB_nearRight << "," << MB_farRight << "]"  <<std::endl; 
	const float MX_nearLeft   = meanMX   - Nmin * sigmaMX,   MX_nearRight   = meanMX   + Nmin * sigmaMX;
	const float MK0s_nearLeft = meanMK0s - Nmin * sigmaMK0s, MK0s_nearRight = meanMK0s + Nmin * sigmaMK0s;
	bool B0_sidebans, X_SgnRegion, K0s_SgnRegion;
	// ----- OUTPUT TREE SETUP ----- //
	OutTree_setup();

	// ----- HISTOGRAMS ----- //

	// ---- TRIGGER INFO ----//
	int nbins = 50;
   TH1F h_Ncandidate_JPsi("Ncandidate_JPsi", "", 4, 0, 4);
   TH1F h_Ncandidate_Rho("Ncandidate_Rho", "", 10, 0, 10);
   TH1F h_Ncandidate_X3872("Ncandidate_X3872", "", 10, 0, 10);
   TH1F h_Ncandidate_K0s("Ncandidate_K0s", "", 5, 0, 5);
   TH1F h_Ncandidate_B0("Ncandidate_B0", "", 10, 0, 10);
	TH1F h_QC_PiTrks("QC_PiTrks", "", 2, 0, 2);
	TH1F h_QC_PiTrksTrig("QC_PiTrksTrig", "", 2, 0, 2);
   TH1F h_TriggerTrkMother("TriggerTrkMother", "", 2, 0, 2);
	TH1F h_N_K0sPerRho("N_K0sPerRho","", 3, 0., 3.);
	TH1F h_N_RhoPerK0s("N_RhoPerK0s","", 10, 0., 10);

	// ---- MASSES POST SELECTION ---- //
   nbins = 100;
   float Mlow, Mhigh;

	Mlow = 2.9; Mhigh = 3.3;
	TH1F h_VTXfit_JPsi_M("VTXfit_JPsi_M", "J\\Psi", nbins, Mlow, Mhigh);	

	Mlow = 0.4; Mhigh = 1.;
	TH1F h_PREfit_Rho_M("PREfit_Rho_M", "\\rho(770)", nbins, Mlow, Mhigh);	

	Mlow = 3.; Mhigh = 5.;
	TH1F h_PREfit_JPsiPi_M("PREfit_JPsiPi_M", "", nbins, Mlow, Mhigh);

	Mlow = 3.67; Mhigh = 4.07;
	TH1F h_PREfit_X3872_M("PREfit_X3872_M", "X(3872)", nbins, Mlow, Mhigh);	
	TH1F h_TOTfit_X3872_M("TOTfit_X3872_M", "X(3872)", nbins, Mlow, Mhigh);	

	Mlow = 0.45; Mhigh = 0.55;
	TH1F h_VTXfit_K0s_M("VTXfit_K0s_M", "\\ K_0^s", nbins, Mlow, Mhigh);	

	Mlow = 4.5; Mhigh = 6;
	TH1F h_TOTfit_B0_M("TOTfit_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	
	TH1F h_SideBands_B0_M("SideBand_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	

	// MVA vars sidebands
	nbins = 50;
	TH1F h_BKG_B0_pT("BKG_B0_pT", "", nbins, 0, 20);
	TH1F h_BKG_B0_LxySign("BKG_B0_LxySign", "", nbins, 0, 100);
	TH1F h_BKG_B0_SVchi2("BKG_B0_SVchi2", "", nbins, 0, 20);
	TH1F h_BKG_B0_SVp("BKG_B0_SVp", "", nbins, 0, 1.);
	TH1F h_BKG_B0_cosA("BKG_B0_cosA", "", nbins, 0.95, 1.);
	TH1F h_BKGb_Pi1_pT("BKGb_Pi1_pT", "", nbins, 0, .3);
	TH1F h_BKGb_Rho_pT("BKGb_Rho_pT", "", nbins, 0, .5);
	TH1F h_BKGb_DR_Pi1B0_Rho("BKGb_DR_Pi1B0_Rho", "", nbins,0, 1.);
	TH1F h_BKGb_Rho_D0("BKGb_Rho_D0", "", nbins, 0., 8.);
	

   // Loop on evens
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry +1  == Nbreak) break;    
      if ((jentry+1) % NPrint==0) std::cout << "--> " << Form("%3.1f",(float)(jentry+1)/nentries* 100.) << " \%"<< std::endl;

      // --> Check if TRIGGER is ON
      isTriggerON = (bool)HLT_DoubleMu4_JpsiTrk_Displaced;
      if ( !isTriggerON ) continue;
      NTriggeredEv++;

		// --> How many K0s per Rho? & How many Rho per K0s?
		//NK0s_per_Rho(&h_N_K0sPerRho);
		//NRho_per_K0s(&h_N_RhoPerK0s);

		// LOOP ON B0 CANDIDATES
		// ripristinate vars
		isSel_MuMu = false;
		prevMu1_idx = -1, prevMu2_idx = -1, prevPi1_idx = -1, prevPi1_idx = -1, prevK0s_idx = -1;
		Ncandidates_JPsi = 0, Ncandidates_Rho = 0, Ncandidates_X3872 = 0, Ncandidates_K0s = 0,  Ncandidates_B0 = 0 ; 
      ToCountJPsi = true, ToCountRho = true, ToCountX = true, ToCountK0s = true; 

		for (UInt_t b = 0; b < nB0; b++){

			isTrkQC = RecoPartFillP4(b);
			if(!isTrkQC)continue;

			PiPi_M = (P4_Reco_Pi1+P4_Reco_Pi2).M();
			MuMu_M = (P4_Reco_Mu1 + P4_Reco_Mu2).M();
			MuMuPi1_M = (P4_Reco_Mu1+P4_Reco_Mu2 + P4_Reco_Pi1).M(); MuMuPi2_M = (P4_Reco_Mu1+P4_Reco_Mu2 + P4_Reco_Pi2).M();
			MuMuPiPi_M = (P4_Reco_Mu1+P4_Reco_Mu2 + P4_Reco_Pi1+P4_Reco_Pi2).M();
			MuMuPiPiK0s_M = (P4_Reco_Mu1+P4_Reco_Mu2 + P4_Reco_Pi1+P4_Reco_Pi2 + P4_Reco_K0s).M();

	//==> TRIGGER SELECTION
      	// --> TRIGGER SELECTION on muons
			ToCountJPsi = (B0_mu1_idx[b] != prevMu1_idx) || (B0_mu2_idx[b] != prevMu2_idx);
			if (ToCountJPsi){
				isSel_MuMu = (bool)ApplyTriggerSelection_Muons(b);
				if (isSel_MuMu) Ncandidates_JPsi++;
				prevMu1_idx = B0_mu1_idx[b]; prevMu2_idx = B0_mu2_idx[b];
			}
			// --> TRIGGER SELECTION on tracks
			TriggeredTrkID = ApplyTriggerSelection_Track(b); // 1= Pi+-Rho 2= Pi--Rho 3=trkK0s
			isSel_Trk = (bool)TriggeredTrkID;
			if (TriggeredTrkID == 1)h_QC_PiTrksTrig.Fill(ProbeTracks_isMatchedToMuon[B0_pi1_idx[b]]);
			if (TriggeredTrkID == 2)h_QC_PiTrksTrig.Fill(ProbeTracks_isMatchedToMuon[B0_pi2_idx[b]]);
			 
			if( !(isSel_MuMu && isSel_Trk) ) continue;
			Ncandidates_B0++;

			//...JPsi
			if(ToCountJPsi)	h_VTXfit_JPsi_M.Fill(MuMu_M);

			//...Rho
			h_QC_PiTrks.Fill(ProbeTracks_isMatchedToMuon[B0_pi1_idx[b]]); h_QC_PiTrks.Fill(ProbeTracks_isMatchedToMuon[B0_pi2_idx[b]]);

			ToCountRho = (prevPi1_idx != B0_pi1_idx[b]) || (prevPi2_idx != B0_pi2_idx[b]) ;

			if (ToCountRho ){
				h_PREfit_Rho_M.Fill(PiPi_M);
				Ncandidates_Rho++;
			} 
			prevPi1_idx = B0_pi1_idx[b]; prevPi2_idx = B0_pi2_idx[b]; 

			//...X(3872)
			ToCountX = ToCountRho || ToCountJPsi;
			if (ToCountX){
				Ncandidates_X3872++;
				h_PREfit_X3872_M.Fill(MuMuPiPi_M);
				h_TOTfit_X3872_M.Fill(B0_finalFit_X_mass[b]);
			}  
			//... K0s
			ToCountK0s  = ROOT::Math::VectorUtil::DeltaR(Prev_P4_Reco_K0s, P4_Reco_K0s)  < 0.0001; 
			if (ToCountK0s){
				h_VTXfit_K0s_M.Fill(B0_K0s_nmcFitted_mass[b]);
				Ncandidates_K0s++;
			}
			Prev_P4_Reco_K0s = P4_Reco_K0s;

			//...B0

			M_B0 = B0_finalFit_mass[b];
			h_TOTfit_B0_M.Fill(M_B0);
			// --> SIDEBANDS
			B0_sidebans =  ((M_B0 > MB_farLeft) && (M_B0 < MB_nearLeft) )  || ( (M_B0 > MB_nearRight) && (M_B0 < MB_farRight));
			M_X3872 = B0_finalFit_X_mass[b];
			X_SgnRegion = (M_X3872 > MX_nearLeft) && (M_X3872 < MX_nearRight);
			M_K0s = B0_K0s_nmcFitted_mass[b];
			K0s_SgnRegion = (M_K0s > MK0s_nearLeft) && (M_K0s < MK0s_nearRight);
			if (B0_sidebans && X_SgnRegion && K0s_SgnRegion){ 
				TOTsidebands_B0++;
				h_SideBands_B0_M.Fill(M_B0);

				P4_Reco_B0.SetPt(B0_finalFit_pt[b]); P4_Reco_B0.SetEta(B0_finalFit_eta[b]); P4_Reco_B0.SetPhi(B0_finalFit_phi[b]), P4_Reco_B0.SetM(M_B0);
				P4_Reco_Pi1.SetPt(B0_finalFit_pi1_pt[b]); P4_Reco_Pi1.SetEta(B0_finalFit_pi1_eta[b]); P4_Reco_Pi1.SetPhi(B0_finalFit_pi1_phi[b]);	
				P4_Reco_Pi2.SetPt(B0_finalFit_pi2_pt[b]); P4_Reco_Pi2.SetEta(B0_finalFit_pi2_eta[b]); P4_Reco_Pi2.SetPhi(B0_finalFit_pi2_phi[b]);
				P4_Reco_Rho = P4_Reco_Pi1 + P4_Reco_Pi2;

				// TTree variables
				M_Rho = B0_finalFit_Rho_mass[b];
				M_mumu = B0_MuMu_fitted_mass[b];
				pTM_B0 = P4_Reco_B0.Pt()/M_B0;			
				LxySign_B0 = B0_lxySign_PV[b];
				SVprob = B0_svprob[b];
				CosAlpha_B0 = B0_cosAlpha_PV[b];
				pT_Rho = P4_Reco_Rho.Pt()/P4_Reco_B0.Pt();
				pT_Pi1 = P4_Reco_Pi1.Pt()/P4_Reco_B0.Pt();
				DR_Pi1B0 = ROOT::Math::VectorUtil::DeltaR(P4_Reco_Pi1, P4_Reco_B0);
				D0_Rho = B0_PiPi_pi1_d0sig[b];
				outTree_->Fill();
			}
			
		}// on B0



		TOTcand_B0 += Ncandidates_B0;
      TOTcand_JPsi += Ncandidates_JPsi;
		if (Ncandidates_B0 > 0 ){
			NTrigSelEV++;
			h_Ncandidate_JPsi.Fill(Ncandidates_JPsi);
			h_Ncandidate_Rho.Fill(Ncandidates_Rho);
			h_Ncandidate_X3872.Fill(Ncandidates_X3872);
			h_Ncandidate_K0s.Fill(Ncandidates_K0s);
			h_Ncandidate_B0.Fill(Ncandidates_B0);
		}	
	

   } // on events

   std::cout << "Triggered events: \t" << NTriggeredEv << std::endl;
   std::cout << "Events with at lesat 1 B0 candidate: \t" << NTrigSelEV << std::endl;
   std::cout << "Muon pairs passing Trigger Selection: \t" << TOTcand_JPsi << std::endl;
   std::cout << "Number of B0 candidates: \t" << TOTcand_B0 << std::endl;
   std::cout << "Number of B0 in the sidebands: \t" << TOTsidebands_B0 << std::endl;

	outFileTree_ = new TFile( outFileTreePath_, "RECREATE");	
	if (!outFileTree_->IsOpen()) std::cout << "	ERROR: cannot open Histo out-file " << outFileTreePath_ << std::endl;
	else std::cout << " OUTPUT --> " << outFileTreePath_ << std::endl;

	outFileTree_->cd();
	outTree_->Write();
	outFileTree_->Close();

	TFile* out_file = new TFile(outFileHistoPath_ , "RECREATE");
	if (!out_file->IsOpen()) std::cout << "	ERROR: cannot open Histo out-file " << outFileHistoPath_ << std::endl;
	else std::cout << " OUTPUT --> " << outFileHistoPath_ << std::endl;
	
	h_QC_PiTrks.Write(); 
	h_QC_PiTrksTrig.Write(); 
	h_N_K0sPerRho.Write();
	h_N_RhoPerK0s.Write();

	h_Ncandidate_JPsi.Write();
	h_Ncandidate_Rho.Write();
	h_Ncandidate_X3872.Write();
	h_Ncandidate_K0s.Write();
	h_Ncandidate_B0.Write();

	h_VTXfit_JPsi_M.Write(); 

	h_PREfit_Rho_M.Write();

	h_PREfit_X3872_M.Write();
	h_TOTfit_X3872_M.Write();

	h_VTXfit_K0s_M.Write(); 

	h_TOTfit_B0_M.Write(); 
	h_SideBands_B0_M.Write();

	out_file->Close();


}//Loop ()


void HLTapply::GetFitParams(){
	
	std::string line;
	int Nline = 0;

	char ParName[10];
	float M1, M2, S1, S2;

	// --- B0 FIT 
	std::ifstream inFileParB0(FileFitParB0);	
	if(!inFileParB0.is_open()) std::cout << "ERROR cannot open " << FileFitParB0 << std::endl;
	while(!inFileParB0.eof()){

		getline(inFileParB0, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 4) sscanf(line.c_str(), "%s %f", ParName, &M1);
		if(Nline == 5) sscanf(line.c_str(), "%s %f", ParName, &S1);
		if(Nline == 7) sscanf(line.c_str(), "%s %f", ParName, &M2);
		if(Nline == 8) sscanf(line.c_str(), "%s %f", ParName, &S2);
	}
	inFileParB0.close();
	meanMB0 = (M1 + M2)*0.5;
	sigmaMB0 = sqrt( S1*S1 + S2*S2);

	// --- K0s FIT 
	Nline = 0;
	std::ifstream inFileParK0s(FileFitParK0s);	
	if(!inFileParK0s.is_open()) std::cout << "ERROR cannot open " << FileFitParK0s << std::endl;
	while(!inFileParK0s.eof()){

		getline(inFileParK0s, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 4) sscanf(line.c_str(), "%s %f", ParName, &M1);
		if(Nline == 5) sscanf(line.c_str(), "%s %f", ParName, &S1);
		if(Nline == 7) sscanf(line.c_str(), "%s %f", ParName, &M2);
		if(Nline == 8) sscanf(line.c_str(), "%s %f", ParName, &S2);
	}
	inFileParK0s.close();
	meanMK0s = (M1 + M2)*0.5;
	sigmaMK0s = sqrt( S1*S1 + S2*S2);

	// --- X FIT 
	Nline = 0;
	std::ifstream inFileParX(FileFitParX);	
	if(!inFileParX.is_open()) std::cout << "ERROR cannot open " << FileFitParX << std::endl;
	while(!inFileParX.eof()){

		getline(inFileParX, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 4) sscanf(line.c_str(), "%s %f", ParName, &M1);
		if(Nline == 5) sscanf(line.c_str(), "%s %f", ParName, &S1);
		if(Nline == 7) sscanf(line.c_str(), "%s %f", ParName, &M2);
		if(Nline == 8) sscanf(line.c_str(), "%s %f", ParName, &S2);
	}
	inFileParX.close();
	meanMX = (M1 + M2)*0.5;
	sigmaMX = sqrt( S1*S1 + S2*S2);

	std::cout << " ---> MASS FIT PRESULTS " << std::endl;
	std::cout << "  <M_B0> = " << meanMB0 << "\tS(M_B0) = " << sigmaMB0 << std::endl;
	std::cout << "  <M_K0s> = " << meanMK0s << "\tS(M_K0s) = " << sigmaMK0s << std::endl;
	std::cout << "  <M_X> = " << meanMX << "\tS(M_X) = " << sigmaMX << std::endl;

	
}//GetFitParams()

int HLTapply::RecoPartFillP4(const int Bidx){
	int TrackQualityCheck = 1;

	//... muons P4
	if(!Muon_softId[B0_mu1_idx[Bidx]] || !Muon_softId[B0_mu2_idx[Bidx]]) TrackQualityCheck = 0; 
	P4_Reco_Mu1.SetPt(B0_MuMu_prefit_mu1_pt[Bidx]); P4_Reco_Mu1.SetEta(B0_MuMu_prefit_mu1_eta[Bidx]); P4_Reco_Mu1.SetPhi(B0_MuMu_prefit_mu1_phi[Bidx]);
	P4_Reco_Mu2.SetPt(B0_MuMu_prefit_mu2_pt[Bidx]); P4_Reco_Mu2.SetEta(B0_MuMu_prefit_mu2_eta[Bidx]); P4_Reco_Mu2.SetPhi(B0_MuMu_prefit_mu2_phi[Bidx]);
	//... tracks P4
	if(ProbeTracks_isMatchedToMuon[B0_pi1_idx[Bidx]] || ProbeTracks_isMatchedToMuon[B0_pi2_idx[Bidx]])TrackQualityCheck = 0;
	P4_Reco_Pi1.SetPt(B0_PiPi_prefit_pi1_pt[Bidx]); P4_Reco_Pi1.SetEta(B0_PiPi_prefit_pi1_eta[Bidx]); P4_Reco_Pi1.SetPhi(B0_PiPi_prefit_pi1_phi[Bidx]);
	P4_Reco_Pi2.SetPt(B0_PiPi_prefit_pi2_pt[Bidx]); P4_Reco_Pi2.SetEta(B0_PiPi_prefit_pi2_eta[Bidx]); P4_Reco_Pi2.SetPhi(B0_PiPi_prefit_pi2_phi[Bidx]);
	P4_Reco_K0s.SetPt(B0_K0s_mcFitted_pt[Bidx]); P4_Reco_K0s.SetEta(B0_K0s_mcFitted_eta[Bidx]); P4_Reco_K0s.SetPhi(B0_K0s_mcFitted_phi[Bidx]);

	return TrackQualityCheck;

}//RecoPartFillP4()




void HLTapply::OutTree_setup(){

	outTree_ = new TTree("B0sidebands", "B0sidebands");
	std::cout << " out tree setting up ... " << std::endl;
	
	outTree_->Branch("M_B0", &M_B0, "M_B0/F");
	outTree_->Branch("M_Rho", &M_Rho, "M_Rho/F");
	outTree_->Branch("M_X3872", &M_X3872, "M_X3872/F");
	outTree_->Branch("M_K0s", &M_K0s, "M_K0s/F");
	outTree_->Branch("M_mumu", &M_mumu, "M_mumu/F");

	outTree_->Branch("pTM_B0", &pTM_B0, "pTM_B0/F");
	outTree_->Branch("LxySign_B0", &LxySign_B0, "LxySign_B0/F");
	outTree_->Branch("SVprob", &SVprob, "SVprob/F");
	outTree_->Branch("CosAlpha_B0", &CosAlpha_B0, "CosAlpha_B0/F");
	outTree_->Branch("pT_Rho", &pT_Rho, "pT_Rho/F");
	outTree_->Branch("pT_Pi1", &pT_Pi1, "pT_Pi1/F");
	outTree_->Branch("DR_Pi1B0", &DR_Pi1B0, "DR_Pi1B0/F");
	outTree_->Branch("D0_Rho", &D0_Rho, "D0_Rho/F");


}//OutTree_setup()


int HLTapply::ApplyTriggerSelection_Muons(const int Bidx){
   // TRIGGER SETTINGS 
   const float Min_Mu_pT = 4.,Max_Mu_eta = 2.5, Max_Mu_dr = 2.;
   const float Min_MuMu_pT = 6.9, Low_MuMu_M = 2.9,  High_MuMu_M = 3.3, Max_MuMu_DCA = 0.5;
   const float Min_MuMu_LxyS = 3, Min_MuMu_cosAlpha = 0.9, Min_MuMu_SVp = 0.1;
	
	int mu1_idx, mu2_idx;
   bool isFiredMu1, isFiredMu2;
	bool isOK_mu1_step0 = false, isOK_mu2_step0 = false, MassCut = false, isOK_mumu_step1 = false, isOK_mumu_step2 = false; 

	int RETURN_VALUE = 0;

	mu1_idx = B0_mu1_idx[Bidx];
	mu2_idx = B0_mu2_idx[Bidx];

	// Fired Mu + muon tracks QUALITY CHECK
	isFiredMu1 = (bool)B0_MuMu_mu1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	isFiredMu2 = (bool)B0_MuMu_mu2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx]; 
	if ( (isFiredMu1 && isFiredMu2) && ( Muon_softId[mu1_idx] && Muon_softId[mu2_idx] )){ 
			// STEP 0
			isOK_mu1_step0 = true;
			if((P4_Reco_Mu1.Pt() < Min_Mu_pT) || (fabs(P4_Reco_Mu1.Eta()) > Max_Mu_eta)  || ( B0_MuMu_mu1_dr[Bidx]) > Max_Mu_dr ) isOK_mu1_step0 = false;
			isOK_mu2_step0 = true;
			if((P4_Reco_Mu2.Pt() < Min_Mu_pT) || (fabs(P4_Reco_Mu2.Eta()) > Max_Mu_eta)  || ( B0_MuMu_mu2_dr[Bidx]) > Max_Mu_dr ) isOK_mu2_step0 = false;

			if (isOK_mu1_step0 && isOK_mu2_step0){ 

				// STEP 1 
				isOK_mumu_step1 = true;
				MassCut = ( (P4_Reco_Mu1 + P4_Reco_Mu2).M() > Low_MuMu_M ) && ( (P4_Reco_Mu1 + P4_Reco_Mu2).M() < High_MuMu_M  );
				if ( !MassCut || ((P4_Reco_Mu1 + P4_Reco_Mu2).Pt() < Min_MuMu_pT ) || ( B0_MuMu_DCA[Bidx] > Max_MuMu_DCA )  )	isOK_mumu_step1 = false;
				// STEP 2	
				isOK_mumu_step2 = true;
				if((B0_MuMu_LxySign[Bidx] < Min_MuMu_LxyS) || (B0_MuMu_cosAlpha[Bidx] < Min_MuMu_cosAlpha ) || (B0_MuMu_sv_prob[Bidx] < Min_MuMu_SVp )) isOK_mumu_step2 = false;
			}
	}

	if (isOK_mu1_step0 && isOK_mu2_step0 && isOK_mumu_step1 && isOK_mumu_step2) RETURN_VALUE = 1; 

	return RETURN_VALUE;

}//ApplyTriggerSelection_Muons()


int HLTapply::ApplyTriggerSelection_Track(const int Bidx){
   //TRIGGER SETTINGS
   const float Min_Trk_pT = 1.2, Max_Trk_eta = 2.5, Min_Trk_D0S = 2.;
   bool isOK_trk_step0 = false, isOK_trk_step1 = false;

   int RETURN_VALUE = 0;

	bool isFired_RhoPi1 = (bool)B0_PiPi_p1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	bool isMatchedToMuon_Rho_Pi1	= ProbeTracks_isMatchedToMuon[B0_pi1_idx[Bidx]];
	bool isFired_RhoPi2 = (bool)B0_PiPi_p2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	bool isMatchedToMuon_Rho_Pi2 = ProbeTracks_isMatchedToMuon[B0_pi2_idx[Bidx]];

	bool isFired_K0sPi1 = (bool)B0_K0s_matchTrack1_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];
	bool isFired_K0sPi2 = (bool)B0_K0s_matchTrack2_fired_DoubleMu4_JpsiTrk_Displaced[Bidx];

	//LEVEL 0 
	isOK_trk_step0 =  (isFired_RhoPi1 || isFired_RhoPi2 || isFired_K0sPi1 || isFired_K0sPi2); // is fired at least 1

   // LEVEL 1
   isOK_trk_step1 = false;
   if (isFired_RhoPi1 && !isMatchedToMuon_Rho_Pi1){ 
		if( (P4_Reco_Pi1.Pt() > Min_Trk_pT) && (fabs(P4_Reco_Pi1.Eta()) < Max_Trk_eta) && (B0_PiPi_pi1_d0sig[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 1;
		}
   } 
   if (isFired_RhoPi2 && !isMatchedToMuon_Rho_Pi2){ 
      if( (P4_Reco_Pi2.Pt() > Min_Trk_pT) && (fabs(P4_Reco_Pi2.Eta()) < Max_Trk_eta) && (B0_PiPi_pi2_d0sig[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 2;
		}
   }
   if (isFired_K0sPi1){
      if( (B0_K0s_matchTrack1_pt[Bidx] > Min_Trk_pT) && (fabs(B0_K0s_matchTrack1_eta[Bidx]) < Max_Trk_eta) && (B0_K0s_matchTrack1_D0sign[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 3;
		}
   }
   if (isFired_K0sPi2){
      if( (B0_K0s_matchTrack2_pt[Bidx] > Min_Trk_pT) && (fabs(B0_K0s_matchTrack2_eta[Bidx]) < Max_Trk_eta) && (B0_K0s_matchTrack2_D0sign[Bidx] > Min_Trk_D0S) ){
			isOK_trk_step1 = true;
			RETURN_VALUE = 3;
		}
   }

	//if (isOK_trk_step0 && isOK_trk_step1) RETURN_VALUE = 1;	

	return RETURN_VALUE;

}//ApplyTriggerSelection_Tracks



void HLTapply::NK0s_per_Rho(TH1* histo){

	int NK0s = 1;
	float epsilon = 0.0001;
	int thisB0;
	int thisPi1, thisPi2;
	float DR_thisK_nextK;
	ROOT::Math::PtEtaPhiMVector P4_ThisK0s(0., 0., 0., mK0s);
	ROOT::Math::PtEtaPhiMVector P4_NextK0s(0., 0., 0., mK0s);

	for(unsigned int b = 0; b < nB0; b++){
		NK0s = 1;
		thisB0 = b;
		thisPi1 = B0_pi1_idx[b]; thisPi2 = B0_pi2_idx[b];
		P4_ThisK0s.SetPt(B0_K0s_mcFitted_pt[b]); P4_ThisK0s.SetPhi(B0_K0s_mcFitted_phi[b]); P4_ThisK0s.SetEta(B0_K0s_mcFitted_eta[b]);
		
		for(unsigned int next = thisB0; next < nB0; next++){
			
			if((thisPi1 == B0_pi1_idx[next]) && (thisPi2 == B0_pi2_idx[next])){	
				P4_NextK0s.SetPt(B0_K0s_mcFitted_pt[next]); P4_NextK0s.SetPhi(B0_K0s_mcFitted_phi[next]); P4_NextK0s.SetEta(B0_K0s_mcFitted_eta[next]);
				DR_thisK_nextK = ROOT::Math::VectorUtil::DeltaR(P4_NextK0s, P4_ThisK0s);
				if (DR_thisK_nextK > epsilon) NK0s++;
			}
		}// next B0 candidates
		histo->Fill(NK0s);
	}


}//NK0s_per_Rho()


void HLTapply::NRho_per_K0s(TH1* histo){

	int NRho = 1;
	float epsilon = 0.0001;
	int thisB0;
	int thisPi1, thisPi2;
	float DR_thisK_nextK;
	ROOT::Math::PtEtaPhiMVector P4_ThisK0s(0., 0., 0., mK0s);
	ROOT::Math::PtEtaPhiMVector P4_NextK0s(0., 0., 0., mK0s);

	for(unsigned int b = 0; b < nB0; b++){
		NRho = 1;
		thisB0 = b;
		thisPi1 = B0_pi1_idx[b]; thisPi2 = B0_pi2_idx[b];
		P4_ThisK0s.SetPt(B0_K0s_mcFitted_pt[b]); P4_ThisK0s.SetPhi(B0_K0s_mcFitted_phi[b]); P4_ThisK0s.SetEta(B0_K0s_mcFitted_eta[b]);
		
		for(unsigned int next = thisB0; next < nB0; next++){
			
			P4_NextK0s.SetPt(B0_K0s_mcFitted_pt[next]); P4_NextK0s.SetPhi(B0_K0s_mcFitted_phi[next]); P4_NextK0s.SetEta(B0_K0s_mcFitted_eta[next]);
			DR_thisK_nextK = ROOT::Math::VectorUtil::DeltaR(P4_NextK0s, P4_ThisK0s);
			
			if( DR_thisK_nextK < epsilon){
				if ((thisPi1 != B0_pi1_idx[next]) || (thisPi2 != B0_pi2_idx[next])) NRho++;
			}
		}// next B0 candidates
		histo->Fill(NRho);
	}
}//NRho_per_K0s()


