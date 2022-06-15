#include "../include/TriggerSelection.h" 

// CONSTRUCTOR-DESTRUCTOR
TriggerSelection::TriggerSelection(TTree *tree) : B0toX3872K0s(tree){ 

	FileFitParB0 = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";
	FileFitParK0s= "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/K0sparams.txt";
	FileFitParX  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/X3872params.txt";

   GenMumP4.SetM(mMuon); GenMupP4.SetM(mMuon);
   GenPimP4.SetM(mPion); GenPipP4.SetM(mPion);
   GenK0sP4.SetM(mK0s);
   
   P4_Reco_Mu1.SetM(mMuon); P4_Reco_Mu2.SetM(mMuon);
   P4_Reco_Pi1.SetM(mPion); P4_Reco_Pi2.SetM(mPion);
   P4_Reco_K0s.SetM(mK0s);
}

TriggerSelection::~TriggerSelection(){
}





void TriggerSelection::Loop() {

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Long64_t Nbreak = nentries + 10, NPrint = (int)nentries/20.; 

   // ----- VARIABLES ----- //
   int Nmatch_Mu = 0, Nmatch_JPsi = 0, Nmatch_Rho =0, Nmatch_X = 0, Nmatch_K0s = 0,  Nmatch_B0 = 0;
   int Ncandidates_JPsi = 0, Ncandidates_X3872 = 0, Ncandidates_B0 = 0, TOTcand_JPsi = 0, TOTcand_B0 = 0;
   int NTriggeredEv = 0, N_MCmatchedEV = 0, N_B0recoEV = 0;
	int TriggeredTrkID;

   bool ToCountJPsi = true, ToCountRho = true, ToCountX = true, ToCountK0s = true;
   bool isTriggerON, isTrkQC, isSel_MuMu, isSel_Trk, isSel_Pim, isSel_Pip;
	bool isMCmatched_JPsi, isMCmatched_Pi1, isMCmatched_Pi2, isMCmatched_Rho, isMCmatched_X3872, isMCmatched_K0s, isMCmatched_B0;	
	bool isSGNregion_X3872 = true, isSGNregion_K0s = true;
	bool inBKGregion_Pi1 = true, inBKGregion_Pi2 = true, inBKGregion_K0s = true;

   int SelMu1_idx = -1, SelMu2_idx = -1;
	int prevMu1_idx = -1, prevMu2_idx = -1, prevPi1_idx = -1, prevPi2_idx = -1, prevK0s_idx = -1;

   float MuMu_M, PiPi_M, MuMuPi1_M, MuMuPi2_M, MuMuPiPi_M, MuMuPiPiK0s_M;
	const float MuMuPi_Mlow = 3.1, MuMuPi_Mhigh  = 3.9;
	float DRpim, DRpip, DRk0s;
	
	

	// ----- OUTPUT TREE SETUP ----- //
	OutTree_setup();

	// ----- LOAD FIT PARAMS ----- //
	GetFitParams();

	// ----- HISTOGRAMS ----- //
	// ---- MC truth matching algorithm ----//
   int nbins = 100;
	TH2F h_DRminVSDpT_Mu("DRminVSDpT_Mu" , "\\mu", nbins, 0., 0.05, nbins, 0., 1.);		
	TH2F h_DRminVSDpT_Pi("DRminVSDpT_Pi" , "\\pi", nbins*2, 0., 0.5, nbins, 0., 2.);		
	TH2F h_DRminVSDpT_K0s("DRminVSDpT_K0s" , "\\ K_s^0", nbins*2, 0., 0.5, nbins, 0., 2.);		
	float M_aiut_B0, M_B0low = 5.25 , M_B0high = 5.35;
	TH2F h_DRminVSDpT_aiut_Pi("DRminVSDpT_aiut_Pi" , "Non-signal\\ \\pi", nbins*2, 0., 0.5, nbins, 0., 2.);		

	TH1F h_MCmatch_WhichPTLFail("MCmatch_WhichPTLFail","", 3 , 0, 3);
	TH1F h_B0cand_WhichPtlMCmiss("B0cand_WhichPtlMCmiss", "", 4, 0, 4);
	TH1F h_N_K0sPerRho("N_K0sPerRho","", 3, 0., 3.);
	TH1F h_N_RhoPerK0s("N_RhoPerK0s","", 10, 0., 10);

	// ---- TRIGGER INFO ----//
   nbins = 50;
   TH1F h_Ncandidate_JPsi("Ncandidate_JPsi", "", 4, 0, 4);
   TH1F h_Ncandidate_X3872("Ncandidate_X3872", "", 10, 0, 10);
   TH1F h_Ncandidate_B0("Ncandidate_B0", "", 10, 0, 10);
	TH1F h_QC_PiTrks("QC_PiTrks", "", 2, 0, 2);
	TH1F h_QC_PiTrksTrig("QC_PiTrksTrig", "", 2, 0, 2);
   TH1F h_TriggerTrkMother("TriggerTrkMother", "", 2, 0, 2);

	// ---- MASSES POST SELECTION ---- //
   nbins = 50;
   float Mlow, Mhigh;

	Mlow = 2.9; Mhigh = 3.3;
	TH1F h_PREfit_SGN_JPsi_M("PREfit_SGN_JPsi_M", "J\\Psi", nbins, Mlow, Mhigh);	
	TH1F h_PREfit_BKG_JPsi_M("PREfit_BKG_JPsi_M", "J\\Psi", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_SGN_JPsi_M("VTXfit_SGN_JPsi_M", "J\\Psi", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_BKG_JPsi_M("VTXfit_BKG_JPsi_M", "J\\Psi", nbins, Mlow, Mhigh);	

	Mlow = 0.4; Mhigh = 1.;
	TH1F h_PREfit_SGN_Rho_M("PREfit_SGN_Rho_M", "\\rho(770)", nbins, Mlow, Mhigh);	
	TH1F h_PREfit_BKG_Rho_M("PREfit_BKG_Rho_M", "\\rho(770)", nbins, Mlow, Mhigh);	

	Mlow = 3.; Mhigh = 5.;
	TH1F h_PREfit_SGN_JPsiPi_M("PREfit_SGN_JPsiPi_M", "", nbins, Mlow, Mhigh);
	TH1F h_PREfit_BKG_JPsiPi_M("PREfit_BKG_JPsiPi_M", "", nbins, Mlow, Mhigh);

   nbins = 100;
	Mlow = 3.7; Mhigh = 4.05;
	TH1F h_PREfit_SGN_X3872_M("PREfit_SGN_X3872_M", "X(3872)", nbins, Mlow, Mhigh);	
	TH1F h_PREfit_BKG_X3872_M("PREfit_BKG_X3872_M", "X(3872)", nbins, Mlow, Mhigh);	
	TH1F h_TOTfit_SGN_X3872_M("TOTfit_SGN_X3872_M", "X(3872)", nbins, Mlow, Mhigh);	
	TH1F h_TOTfit_BKG_X3872_M("TOTfit_BKG_X3872_M", "X(3872)", nbins, Mlow, Mhigh);	

	Mlow = 0.38; Mhigh = 0.62;
	TH1F h_PREfit_SGN_K0s_M("PREfit_SGN_K0s_M", "\\ K_0^s", nbins, Mlow, Mhigh);	
	TH1F h_PREfit_BKG_K0s_M("PREfit_BKG_K0s_M", "\\ K_0^s", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_SGN_K0s_M("VTXfit_SGN_K0s_M", "\\ K_0^s", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_BKG_K0s_M("VTXfit_BKG_K0s_M", "\\ K_0^s", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_BKG_B0picco_K0s_M("VTXfit_BKG_B0picco_K0s_M", "", nbins, Mlow, Mhigh);	

	nbins = 50;
	Mlow = 5.; Mhigh = 5.5;
	TH1F h_PREfit_SGN_B0_M("PREfit_SGN_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	
	TH1F h_PREfit_BKG_B0_M("PREfit_BKG_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_SGN_B0_M("VTXfit_SGN_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	
	TH1F h_VTXfit_BKG_B0_M("VTXfit_BKG_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	
	TH1F h_TOTfit_SGN_B0_M("TOTfit_SGN_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	
	TH1F h_TOTfit_BKG_B0_M("TOTfit_BKG_B0_M", "\\ B_0", nbins, Mlow, Mhigh);	



   // Loop on evens
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry +1  == Nbreak) break;    
      if ((jentry+1) % NPrint==0) cout << "--> " << Form("%3.0f", (float)(jentry+1)/nentries* 100.) << " \%"<< endl;

      // --> Check if TRIGGER is ON
      isTriggerON = (bool)HLT_DoubleMu4_JpsiTrk_Displaced;
      if ( !isTriggerON ) continue;
      NTriggeredEv++;

      // --> Find MC-truth
      GenPartFillP4();

      // --> Match MC-truth with reco tracks
      MCmatch_Mum_Idx = -1; MCmatch_Mup_Idx = -1;
      MCmatch_Pim_Idx = -1; MCmatch_Pip_Idx = -1;
      MCmatch_K0s_Idx = -1;
		MCmatch_B0_Idx = -1;
      MCtruthMatching();
		h_DRminVSDpT_Mu.Fill(MCmatch_Mum_DRmin, MCmatch_Mum_DpT); h_DRminVSDpT_Mu.Fill(MCmatch_Mup_DRmin, MCmatch_Mup_DpT);
		h_DRminVSDpT_Pi.Fill(MCmatch_Pim_DRmin, MCmatch_Pim_DpT); h_DRminVSDpT_Pi.Fill(MCmatch_Pip_DRmin, MCmatch_Pip_DpT);
		h_DRminVSDpT_K0s.Fill(MCmatch_K0s_DRmin, MCmatch_K0s_DpT);
		if ( !isMCmatchingFailed()) N_MCmatchedEV++;
		WhichPtl_MCmissed(&h_MCmatch_WhichPTLFail);

		// --> How many K0s per Rho? & How many Rho per K0s?
		NK0s_per_Rho(&h_N_K0sPerRho);
		NRho_per_K0s(&h_N_RhoPerK0s);

		// LOOP ON B0 CANDIDATES
		// ripristinate vars
		isSel_MuMu = false;
		prevMu1_idx = -1, prevMu2_idx = -1, prevPi1_idx = -1, prevPi1_idx = -1, prevK0s_idx = -1;
		Ncandidates_JPsi = 0, Ncandidates_X3872 = 0, Ncandidates_B0 = 0 ; 
      ToCountJPsi = true, ToCountRho = true, ToCountX = true, ToCountK0s = true; 

		B0_SGN_idx = -1;
		nBKG_B0 = 0;

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
			TriggerSel_event = jentry;

	//==> MC truth matching
			//...JPsi
			isMCmatched_JPsi = ((B0_mu1_idx[b] == MCmatch_Mum_Idx) && (B0_mu2_idx[b] == MCmatch_Mup_Idx)) 
									|| ((B0_mu1_idx[b] == MCmatch_Mup_Idx) && (B0_mu2_idx[b] == MCmatch_Mum_Idx));
			if (isMCmatched_JPsi && ToCountJPsi){
				h_PREfit_SGN_JPsi_M.Fill(MuMu_M);
				h_VTXfit_SGN_JPsi_M.Fill(B0_MuMu_fitted_mass[b]);
			}else if(ToCountJPsi){
				h_PREfit_BKG_JPsi_M.Fill(MuMu_M);
				h_VTXfit_BKG_JPsi_M.Fill(B0_MuMu_fitted_mass[b]);
			}

			//...Rho
			h_QC_PiTrks.Fill(ProbeTracks_isMatchedToMuon[B0_pi1_idx[b]]); h_QC_PiTrks.Fill(ProbeTracks_isMatchedToMuon[B0_pi2_idx[b]]);

			isMCmatched_Pi1 = (B0_pi1_idx[b] == MCmatch_Pim_Idx) || (B0_pi1_idx[b] == MCmatch_Pip_Idx); 
			if(!isMCmatched_Pi1 && !inBKGregion("Pi1")) continue;
			isMCmatched_Pi2 = (B0_pi2_idx[b] == MCmatch_Pim_Idx) || (B0_pi2_idx[b] == MCmatch_Pip_Idx);
			if(!isMCmatched_Pi2 && !inBKGregion("Pi2")) continue;
			isMCmatched_Rho = isMCmatched_Pi1 && isMCmatched_Pi2; 

			//PrintGenLevel(b, isPip, isMCmatched_Pi1, isMCmatched_Pi2 );	
			ToCountRho = (prevPi1_idx != B0_pi1_idx[b]) || (prevPi2_idx != B0_pi2_idx[b]) ;

			if(isMCmatched_Rho && ToCountRho){
				h_PREfit_SGN_Rho_M.Fill(PiPi_M);
				//if (PiPi_M < 0.7) DebugRho();
				Nmatch_Rho++;
			}else if (ToCountRho ){ 
				h_PREfit_BKG_Rho_M.Fill(PiPi_M);
				h_DRminVSDpT_aiut_Pi.Fill(ROOT::Math::VectorUtil::DeltaR(GenPimP4, P4_Reco_Pi1), DeltaPT(GenPimP4, P4_Reco_Pi1));
				h_DRminVSDpT_aiut_Pi.Fill(ROOT::Math::VectorUtil::DeltaR(GenPipP4, P4_Reco_Pi2), DeltaPT(GenPipP4, P4_Reco_Pi2));
			}  
			prevPi1_idx = B0_pi1_idx[b]; prevPi2_idx = B0_pi2_idx[b]; 

			//...X(3872)
			isMCmatched_X3872 = isMCmatched_Rho && isMCmatched_JPsi; 
			ToCountX = ToCountRho || ToCountJPsi;
			if( isMCmatched_X3872 && ToCountX){ 
				Nmatch_X++; Ncandidates_X3872++; 
				h_PREfit_SGN_X3872_M.Fill(MuMuPiPi_M);
				h_TOTfit_SGN_X3872_M.Fill(B0_finalFit_X_mass[b]);
			}else if (ToCountX){
				Ncandidates_X3872++;
				h_PREfit_BKG_X3872_M.Fill(MuMuPiPi_M);
				h_TOTfit_BKG_X3872_M.Fill(B0_finalFit_X_mass[b]);
			}  
			isSGNregion_X3872 = ( B0_finalFit_X_mass[b] > MX_nearLeft ) && ( B0_finalFit_X_mass[b] < MX_nearRight );

			//... K0s
			isMCmatched_K0s = (fabs(ROOT::Math::VectorUtil::DeltaR(GenK0sP4, P4_Reco_K0s) - MCmatch_K0s_DRmin) < 0.0001) && (MCmatch_K0s_Idx  >= 0); //(MCmatch_K0s_Idx == B0_k0short_idx[b]);
			ToCountK0s = B0_k0short_idx[b] != prevK0s_idx;
			if (isMCmatched_K0s && ToCountK0s){
				//std::cout << "DR @ idx "<< B0_k0short_idx[b] << "\t" << ROOT::Math::VectorUtil::DeltaR(GenK0sP4, P4_Reco_K0s) << "\t DRmin @ idx "<< MCmatch_K0s_Idx<< "\t" << MCmatch_K0s_DRmin << std::endl;
				Nmatch_K0s++;
				h_PREfit_SGN_K0s_M.Fill(P4_Reco_K0s.M());
				h_VTXfit_SGN_K0s_M.Fill(B0_K0s_nmcFitted_mass[b]);
			}else if (ToCountK0s){
				h_PREfit_BKG_K0s_M.Fill(P4_Reco_K0s.M());
				h_VTXfit_BKG_K0s_M.Fill(B0_K0s_nmcFitted_mass[b]);
			}
			prevK0s_idx = B0_k0short_idx[b];
			isSGNregion_K0s = ( B0_K0s_nmcFitted_mass[b] > MK0s_nearLeft) && ( B0_K0s_nmcFitted_mass[b] < MK0s_nearRight );
			

			//...B0
			if (isMCmatched_JPsi && isMCmatched_Rho && isMCmatched_K0s && isSGNregion_X3872 && isSGNregion_K0s){
				// TTree variables 
				B0_SGN_idx = b;
				Nmatch_B0++;
			
				h_PREfit_SGN_B0_M.Fill(MuMuPiPiK0s_M); 
				h_VTXfit_SGN_B0_M.Fill(B0_fitted_mass_womc[b]);
				h_TOTfit_SGN_B0_M.Fill(B0_finalFit_mass[b]);
			}else{
				// TTree variables 
				B0_BKG_idx[nBKG_B0] = b;
				B0_BKG_isTrueJPsi[nBKG_B0] = isMCmatched_JPsi;
				B0_BKG_isTruePi1[nBKG_B0] = isMCmatched_Pi1;
				B0_BKG_isTruePi2[nBKG_B0] = isMCmatched_Pi2;
				B0_BKG_isTrueRho[nBKG_B0] = isMCmatched_Rho;
				B0_BKG_isTrueK0s[nBKG_B0] = isMCmatched_K0s;
				nBKG_B0++;
				
				h_PREfit_BKG_B0_M.Fill(MuMuPiPiK0s_M);
				h_VTXfit_BKG_B0_M.Fill(B0_fitted_mass_womc[b]);
				h_TOTfit_BKG_B0_M.Fill(B0_finalFit_mass[b]);
				
				B0cand_PTLmissed(isMCmatched_JPsi, isMCmatched_Rho, isMCmatched_K0s, &h_B0cand_WhichPtlMCmiss);
			}
		}// on B0




      TOTcand_JPsi += Ncandidates_JPsi;
		if (Ncandidates_B0 > 0 ){
			outTree_->Fill();	
			N_B0recoEV++;
			h_Ncandidate_JPsi.Fill(Ncandidates_JPsi);
			h_Ncandidate_B0.Fill(Ncandidates_B0);
			h_Ncandidate_X3872.Fill(Ncandidates_X3872);
		}	
	

   } // on events

   std::cout << "Triggered events: " << NTriggeredEv << std::endl;
	std::cout << "MC matched events: " << N_MCmatchedEV << std::endl;
   std::cout << "Muon pairs passing Trigger Selection:" << TOTcand_JPsi << std::endl;
   std::cout << "Events with at lesat 1 B0 candidate:" << N_B0recoEV << std::endl;
   std::cout << "Rho matching MC thruth:" << Nmatch_Rho << std::endl;  
   std::cout << "X(3872) matching MC thruth:" << Nmatch_X << std::endl;
   std::cout << "K0s matching MC thruth:" << Nmatch_K0s << std::endl;
   std::cout << "B0 matching MC thruth:" << Nmatch_B0 << std::endl;


	outFileTree_ = new TFile("./results/B0_SGNvsBKG_anlaysis.root", "RECREATE");	
	outFileTree_->cd();
	outTree_->Write();
	outFileTree_->Close();


	TFile* out_file = new TFile("./plots/TriggerSelection.root", "RECREATE");

	h_DRminVSDpT_Mu.Write();
	h_DRminVSDpT_Pi.Write();
	h_DRminVSDpT_K0s.Write();
	h_DRminVSDpT_aiut_Pi.Write();
	
	h_MCmatch_WhichPTLFail.Write();
	h_B0cand_WhichPtlMCmiss.Write();
	h_N_K0sPerRho.Write();
	h_N_RhoPerK0s.Write();

	h_Ncandidate_JPsi.Write();
	h_QC_PiTrks.Write();
	h_QC_PiTrksTrig.Write();
	h_Ncandidate_X3872.Write();
	h_Ncandidate_B0.Write();

	h_PREfit_SGN_JPsi_M.Write();
	h_VTXfit_SGN_JPsi_M.Write();
	h_PREfit_BKG_JPsi_M.Write();
	h_VTXfit_BKG_JPsi_M.Write();

	h_PREfit_SGN_JPsiPi_M.Write();
	h_PREfit_BKG_JPsiPi_M.Write();

	h_PREfit_SGN_Rho_M.Write();
	h_PREfit_BKG_Rho_M.Write();

	h_VTXfit_SGN_K0s_M.Write();
	h_VTXfit_BKG_K0s_M.Write();
	h_VTXfit_BKG_B0picco_K0s_M.Write();

	h_PREfit_SGN_X3872_M.Write();
	h_TOTfit_SGN_X3872_M.Write();
	h_PREfit_BKG_X3872_M.Write();
	h_TOTfit_BKG_X3872_M.Write();


	h_PREfit_SGN_B0_M.Write(); 
	h_VTXfit_SGN_B0_M.Write(); 
	h_TOTfit_SGN_B0_M.Write(); 
	h_PREfit_BKG_B0_M.Write(); 
	h_VTXfit_BKG_B0_M.Write(); 
	h_TOTfit_BKG_B0_M.Write(); 

	out_file->Close();


}//Loop ()


void TriggerSelection::GetFitParams(){
	
	std::string line;
	int Nline = 0;

	char ParName[30];
	double err;
	/*// --- B0 FIT 
	std::ifstream inFileParB0(FileFitParB0);	
	if(!inFileParB0.is_open()) std::cout << "ERROR cannot open " << FileFitParB0 << std::endl;
	while(!inFileParB0.eof()){

		getline(inFileParB0, line); Nline++;
		//std::cout << Nline << "\t" << line << std::endl;
		if(line.c_str()[0] == '#') continue;
		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &MB_nearLeft);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &MB_nearRight);
		if(Nline == 5) sscanf(line.c_str(), "%s %lf", ParName, &MB_farRight);
		if(Nline == 6) sscanf(line.c_str(), "%s %lf", ParName, &MB_farLeft);
	}
	inFileParB0.close();*/

	// --- K0s FIT 
	Nline = 0;
	std::ifstream inFileParK0s(FileFitParK0s);	
	if(!inFileParK0s.is_open()) std::cout << "ERROR cannot open " << FileFitParK0s << std::endl;
	while(!inFileParK0s.eof()){

		getline(inFileParK0s, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &MK0s_nearLeft);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &MK0s_nearRight);
	}
	inFileParK0s.close();

	// --- X FIT 
	Nline = 0;
	std::ifstream inFileParX(FileFitParX);	
	if(!inFileParX.is_open()) std::cout << "ERROR cannot open " << FileFitParX << std::endl;
	while(!inFileParX.eof()){

		getline(inFileParX, line); Nline++;
		if(line.c_str()[0] == '#') continue;

		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &MX_nearLeft);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &MX_nearRight);
	}
	inFileParX.close();

	std::cout << " ---> MASS FIT RESULTS " << std::endl;
	//std::cout << "    B0 sidebands   [" << MB_farLeft << "," << MB_nearLeft << "]" << " + [" << MB_nearRight << "," << MB_farRight << "]"  <<std::endl; 
	std::cout << " X(3872) SGNregion [" << MX_nearLeft << "," << MX_nearRight << "]" << std::endl; 
	std::cout << "   K0s SGNregion   [" << MK0s_nearLeft << "," << MK0s_nearRight << "]" << std::endl; 
	
}//GetFitParams()

void TriggerSelection::GenPartFillP4(){
   UInt_t MumIdx = -1, MupIdx = -1, PimIdx = -1 , PipIdx = -1, K0sIdx = -1, RhoIdx = -1;
   for (UInt_t g = 0; g < nGenPart; g++){
      if( (GenPart_pdgId[g] ==  isMum) &&  (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi)) MumIdx = g;
      if( (GenPart_pdgId[g] == -isMum) &&  (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isJPsi)) MupIdx = g;
      if( (GenPart_pdgId[g] == -isPip) &&  
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isRho) &&
			 (GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isX3872)) PimIdx = g;
      if( (GenPart_pdgId[g] ==  isPip) &&  
			 (GenPart_pdgId[GenPart_genPartIdxMother[g]] == isRho) &&
			(GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]] == isX3872)	) PipIdx = g;
      if( (GenPart_pdgId[g] ==  isRho) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == isX3872)) RhoIdx = g;
      if( (GenPart_pdgId[g] ==  isK0s) &&  (abs(GenPart_pdgId[GenPart_genPartIdxMother[g]]) == abs(isB0))) K0sIdx = g;

   }//on generated ptl
   // mu- 				mu+
   if(MumIdx > 0){   
      GenMumP4.SetPt(GenPart_pt[MumIdx]);   
      GenMumP4.SetEta(GenPart_eta[MumIdx]);
      GenMumP4.SetPhi(GenPart_phi[MumIdx]); 
		//std::cout << "Mu - == pt " << GenMumP4.Pt() << " eta " << GenMumP4.Eta() << " mass " << GenPart_mass[MumIdx] << std::endl;
   } else std::cout << "No generated mu- is found" << std::endl;   

   if(MupIdx > 0){
      GenMupP4.SetPt(GenPart_pt[MupIdx]);
      GenMupP4.SetEta(GenPart_eta[MupIdx]);
      GenMupP4.SetPhi(GenPart_phi[MupIdx]);
		//std::cout << "Mu + == pt " << GenMupP4.Pt() << " eta " << GenMupP4.Eta() << " mass " << GenPart_mass[MupIdx] << std::endl;

   }else std::cout << "No generated mu+ is found" << std::endl;

   // pi- 				pi+   
   if(PimIdx > 0){
      GenPimP4.SetPt(GenPart_pt[PimIdx]);   
      GenPimP4.SetEta(GenPart_eta[PimIdx]); 
      GenPimP4.SetPhi(GenPart_phi[PimIdx]); 
		//std::cout << "Pi - == pt " << GenPimP4.Pt() << " eta " << GenPimP4.Eta() << " mass " << GenPart_mass[PimIdx] << std::endl;
   }else std::cout << "No generated pi- is found" << std::endl;
   if(PipIdx > 0){

      GenPipP4.SetPt(GenPart_pt[PipIdx]);
      GenPipP4.SetEta(GenPart_eta[PipIdx]);
      GenPipP4.SetPhi(GenPart_phi[PipIdx]);
		//std::cout << "Pi + == pt " << GenPipP4.Pt() << " eta " << GenPipP4.Eta() << " mass " << GenPart_mass[PipIdx] << std::endl;
   }else std::cout << "No generated pi+ is found" << std::endl;
	// Rho
   if(RhoIdx > 0 ){
      GenRhoP4.SetPt(GenPart_pt[RhoIdx]);
      GenRhoP4.SetEta(GenPart_eta[RhoIdx]);
      GenRhoP4.SetPhi(GenPart_phi[RhoIdx]);
      GenRhoP4.SetM(GenPart_mass[RhoIdx]);
   }else std::cout << "No generated Rho is found" << std::endl;

   // K0s
   if(K0sIdx > 0 ){
      GenK0sP4.SetPt(GenPart_pt[K0sIdx]);
      GenK0sP4.SetEta(GenPart_eta[K0sIdx]);
      GenK0sP4.SetPhi(GenPart_phi[K0sIdx]);
   }else std::cout << "No generated K0s is found" << std::endl;
}//GenPartFillP4()

int TriggerSelection::RecoPartFillP4(const int Bidx){
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


void TriggerSelection::MCtruthMatching(const bool verbose){

   const double DRmin_threshold = 0.03;
   const double DpT_threshold = 0.5;
   float DeltaPt;
 
   // ..... muons ..... //	
   ROOT::Math::PtEtaPhiMVector P4_RecoMu1, P4_RecoMu2;
   float DRminMum = 100., DR_gMum_rMu1, DR_gMum_rMu2, DRminMum_DpT = 0., DRminMup = 100., DR_gMup_rMu1, DR_gMup_rMu2, DRminMup_DpT = 0.;
   int DRminMum_Idx = -1, DRminMup_Idx = -1; 
   // ..... pions ..... //      
   ROOT::Math::PtEtaPhiMVector P4_RecoPi1, P4_RecoPi2;
   float DRminPim = 100., DR_gPim_rPi1, DR_gPim_rPi2, DRminPim_DpT = 0., DRminPip = 100., DR_gPip_rPi1, DR_gPip_rPi2, DRminPip_DpT = 0.;
   int DRminPim_Idx = -1, DRminPip_Idx = -1;


   // ..... K0s ..... //
   ROOT::Math::PtEtaPhiMVector P4_RecoK0s;
   float DRminK0s = 100., DR_gK0s_rK0s, DRminK0s_DpT =0;
   int DRminK0s_Idx = -1;

   for (UInt_t b = 0; b < nB0; b++){
      // ..... muons ..... //
		if (Muon_softId[B0_mu1_idx[b]] && Muon_softId[B0_mu2_idx[b]]){ // softId for both muons is required

			P4_RecoMu1.SetPt(B0_MuMu_prefit_mu1_pt[b]); P4_RecoMu1.SetEta(B0_MuMu_prefit_mu1_eta[b]); P4_RecoMu1.SetPhi(B0_MuMu_prefit_mu1_phi[b]);
			P4_RecoMu2.SetPt(B0_MuMu_prefit_mu2_pt[b]); P4_RecoMu2.SetEta(B0_MuMu_prefit_mu2_eta[b]); P4_RecoMu2.SetPhi(B0_MuMu_prefit_mu2_phi[b]);

			DR_gMum_rMu1 = ROOT::Math::VectorUtil::DeltaR(GenMumP4, P4_RecoMu1); // mu1
			DR_gMup_rMu1 = ROOT::Math::VectorUtil::DeltaR(GenMupP4, P4_RecoMu1);
			DR_gMum_rMu2 = ROOT::Math::VectorUtil::DeltaR(GenMumP4, P4_RecoMu2); // mu2
			DR_gMup_rMu2 = ROOT::Math::VectorUtil::DeltaR(GenMupP4, P4_RecoMu2);

			// DR(mu-, mu1) < DRmin(mu-) + m1 is nearer to mu- than mu2 + DRmin threshold
			// mu-
			if( (DR_gMum_rMu1 < DRminMum) && (DR_gMum_rMu1 < DR_gMum_rMu2) ){
				DRminMum = DR_gMum_rMu1;
				DRminMum_DpT = DeltaPT(GenMumP4, P4_RecoMu1);
				DRminMum_Idx = B0_mu1_idx[b];
			}
			if( (DR_gMum_rMu2 < DRminMum) && (DR_gMum_rMu2 < DR_gMum_rMu1) ){
				DRminMum = DR_gMum_rMu2;
				DRminMum_DpT = DeltaPT(GenMumP4, P4_RecoMu2);
				DRminMum_Idx = B0_mu2_idx[b];
			}
			// mu+
			if( (DR_gMup_rMu1 < DRminMup) && (DR_gMup_rMu1 < DR_gMup_rMu2) ){
				DRminMup = DR_gMup_rMu1;
				DRminMup_DpT = DeltaPT(GenMupP4, P4_RecoMu1);
				DRminMup_Idx = B0_mu1_idx[b];
			}
			if( (DR_gMup_rMu2 < DRminMup) && (DR_gMup_rMu2 < DR_gMup_rMu1) ){
				DRminMup = DR_gMup_rMu2;
				DRminMup_DpT = DeltaPT(GenMupP4, P4_RecoMu2);
				DRminMup_Idx = B0_mu2_idx[b];
			}

		}// ... muons //

		// ..... pions ..... //
		if (!ProbeTracks_isMatchedToMuon[B0_pi1_idx[b]] && !ProbeTracks_isMatchedToMuon[B0_pi2_idx[b]]){ // isMatchedToMuon must be false
			P4_RecoPi1.SetPt(B0_PiPi_prefit_pi1_pt[b]); P4_RecoPi1.SetEta(B0_PiPi_prefit_pi1_eta[b]); P4_RecoPi1.SetPhi(B0_PiPi_prefit_pi1_phi[b]);
			P4_RecoPi2.SetPt(B0_PiPi_prefit_pi2_pt[b]); P4_RecoPi2.SetEta(B0_PiPi_prefit_pi2_eta[b]); P4_RecoPi2.SetPhi(B0_PiPi_prefit_pi2_phi[b]);

			DR_gPim_rPi1 = ROOT::Math::VectorUtil::DeltaR(GenPimP4, P4_RecoPi1);//pi1
			DR_gPip_rPi1 = ROOT::Math::VectorUtil::DeltaR(GenPipP4, P4_RecoPi1);
			DR_gPim_rPi2 = ROOT::Math::VectorUtil::DeltaR(GenPimP4, P4_RecoPi2);//pi2
			DR_gPip_rPi2 = ROOT::Math::VectorUtil::DeltaR(GenPipP4, P4_RecoPi2);
			// pi-
			if( (DR_gPim_rPi1 < DRminPim) && (DR_gPim_rPi1 < DR_gPim_rPi2) ){
				DRminPim = DR_gPim_rPi1;
				DRminPim_DpT = DeltaPT(GenPimP4, P4_RecoPi1);
				DRminPim_Idx = B0_pi1_idx[b];
			}
			if( (DR_gPim_rPi2 < DRminPim) && (DR_gPim_rPi2 < DR_gPim_rPi1) ){
				DRminPim = DR_gPim_rPi2;
				DRminPim_DpT = DeltaPT(GenPimP4, P4_RecoPi2);
				DRminPim_Idx = B0_pi2_idx[b];
			}
			// pi+
			if( (DR_gPip_rPi1 < DRminPip) && (DR_gPip_rPi1 < DR_gPip_rPi2) ){
				DRminPip = DR_gPip_rPi1;
				DRminPip_DpT = DeltaPT(GenPipP4, P4_RecoPi1);
				DRminPip_Idx = B0_pi1_idx[b];
			}
			if( (DR_gPip_rPi2 < DRminPip) && (DR_gPip_rPi2 < DR_gPip_rPi1) ){
				DRminPip = DR_gPip_rPi2;
				DRminPip_DpT = DeltaPT(GenPipP4, P4_RecoPi2);
				DRminPip_Idx = B0_pi2_idx[b];
			}


		}// ... pions//

		// ..... K0s ..... //
		P4_RecoK0s.SetPt(B0_K0s_mcFitted_pt[b]); P4_RecoK0s.SetEta(B0_K0s_mcFitted_eta[b]); P4_RecoK0s.SetPhi(B0_K0s_mcFitted_phi[b]);
		DR_gK0s_rK0s = ROOT::Math::VectorUtil::DeltaR(GenK0sP4, P4_RecoK0s);
		if ( (DR_gK0s_rK0s < DRminK0s) ){ 
			DRminK0s = DR_gK0s_rK0s;
			DRminK0s_DpT = DeltaPT(GenK0sP4, P4_RecoK0s);
			DRminK0s_Idx = B0_k0short_idx[b];
		}


	}// on B0 candidates

	if(DRminMum < DRmin_threshold) MCmatch_Mum_Idx = DRminMum_Idx; 
	MCmatch_Mum_DRmin = DRminMum; MCmatch_Mum_DpT = DRminMum_DpT; 
	if(DRminMup < DRmin_threshold) MCmatch_Mup_Idx = DRminMup_Idx; 
	MCmatch_Mup_DRmin = DRminMup; MCmatch_Mup_DpT = DRminMup_DpT; 


	if( (DRminPim < DRmin_threshold) && (DRminPim_DpT < DpT_threshold)) MCmatch_Pim_Idx = DRminPim_Idx;
	MCmatch_Pim_DRmin = DRminPim; MCmatch_Pim_DpT = DRminPim_DpT;
	if( (DRminPip < DRmin_threshold) && (DRminPip_DpT < DpT_threshold)) MCmatch_Pip_Idx = DRminPip_Idx;
	MCmatch_Pip_DRmin = DRminPip; MCmatch_Pip_DpT = DRminPip_DpT; 

	if( (DRminK0s < DRmin_threshold) && (DRminK0s_DpT < DpT_threshold)) MCmatch_K0s_Idx = DRminK0s_Idx;
	MCmatch_K0s_DRmin = DRminK0s; MCmatch_K0s_DpT = DRminK0s_DpT;
 
	if (verbose){
		std::cout << "MC matching indices " << std::endl;
		std::cout << "mu- " << MCmatch_Mum_Idx << "\t mu+ " << MCmatch_Mup_Idx << std::endl;
		std::cout << "pi- " << MCmatch_Pim_Idx << "\t pi+ " << MCmatch_Pip_Idx << std::endl;
		std::cout << "K0s " << MCmatch_K0s_Idx << std::endl;
	}

}//MCtruthMatching()

float TriggerSelection::DeltaPT(ROOT::Math::PtEtaPhiMVector genV, ROOT::Math::PtEtaPhiMVector recV){
	return TMath::Abs(genV.Pt() - recV.Pt()) / genV.Pt();
}//DeltaPT()

bool TriggerSelection::isMCmatchingFailed(){

	bool isFailed = false;
	if (( MCmatch_Mum_Idx < 0 ) || ( MCmatch_Mup_Idx < 0)) isFailed = true;
	if (( MCmatch_Pim_Idx < 0 ) || ( MCmatch_Pip_Idx < 0)) isFailed = true;
	if (  MCmatch_K0s_Idx < 0 ) isFailed = true;

	return isFailed;
}//isMCmatchingFailed()


void TriggerSelection::OutTree_setup(){

	outTree_ = new TTree("B0_SGN_BKG", "B0_SGN_BKG");
	std::cout << " out tree setting up ... " << std::endl;
	
	outTree_->Branch("TriggerSel_event", &TriggerSel_event, "TriggerSel_event/L");
	outTree_->Branch("B0_SGN_idx", &B0_SGN_idx, "B0_SGN_idx/I");
	outTree_->Branch("nBKG_B0", &nBKG_B0, "nBKG_B0/i");
	outTree_->Branch("B0_BKG_idx", B0_BKG_idx, "B0_SGN_idx[nBKG_B0]/I");
	outTree_->Branch("B0_BKG_isTrueJPsi", B0_BKG_isTrueJPsi, "B0_BKG_isTrueJPsi[nBKG_B0]/O");
	outTree_->Branch("B0_BKG_isTruePi1", B0_BKG_isTruePi1, "B0_BKG_isTruePi1[nBKG_B0]/O");
	outTree_->Branch("B0_BKG_isTruePi2", B0_BKG_isTruePi2, "B0_BKG_isTruePi2[nBKG_B0]/O");
	outTree_->Branch("B0_BKG_isTrueRho", B0_BKG_isTrueRho, "B0_BKG_isTrueRho[nBKG_B0]/O");
	outTree_->Branch("B0_BKG_isTrueK0s", B0_BKG_isTrueK0s, "B0_BKG_isTrueK0s[nBKG_B0]/O");


}//OutTree_setup()


int TriggerSelection::ApplyTriggerSelection_Muons(const int Bidx){
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


int TriggerSelection::ApplyTriggerSelection_Track(const int Bidx){
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

void TriggerSelection::WhichPtl_MCmissed(TH1* histo){

	
	if (( MCmatch_Mum_Idx < 0 ) || ( MCmatch_Mup_Idx < 0)) histo->Fill(0.);
	if (( MCmatch_Pim_Idx < 0 ) || ( MCmatch_Pip_Idx < 0)) histo->Fill(1.); 
	if (  MCmatch_K0s_Idx < 0 ) histo->Fill(2.);

}//WhichPtl_MCmissed()


bool TriggerSelection::inBKGregion(const TString ptl){

	const double DRlowLim = 0.1;
	const double DpTlowLim = 0.8;
	ROOT::Math::PtEtaPhiMVector P4_Reco;
	bool isFar_Pim, isFar_Pip;
	int isFar = 0;

	if(ptl == "Pi1")P4_Reco = P4_Reco_Pi1; 
	if(ptl == "Pi2")P4_Reco = P4_Reco_Pi2; 

	isFar_Pim = (ROOT::Math::VectorUtil::DeltaR(GenPimP4 ,P4_Reco) > DRlowLim) && (DeltaPT(GenPimP4, P4_Reco));
	isFar_Pip = (ROOT::Math::VectorUtil::DeltaR(GenPipP4 ,P4_Reco) > DRlowLim) && (DeltaPT(GenPipP4, P4_Reco));

	if(isFar_Pim && isFar_Pip) isFar = 1;

	return isFar;
}//MCtruthProximity()


void TriggerSelection::B0cand_PTLmissed(const bool isMCmatchedJPsi, const bool isMCmatchedRho, const bool isMCmatchedK0s, TH1* histo){

	if(!isMCmatchedJPsi) histo->Fill(0.5);
	if(!isMCmatchedRho ) histo->Fill(1.5);
	if(!isMCmatchedK0s ) histo->Fill(2.5);
	if(!isMCmatchedRho && !isMCmatchedK0s) histo->Fill(3.5);


}//B0cand_PTLmissed()



void TriggerSelection::NK0s_per_Rho(TH1* histo){

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


void TriggerSelection::NRho_per_K0s(TH1* histo){

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




void TriggerSelection::PrintGenLevel(const int Bidx, const int ptlID, const bool isMCmatched_Pi1, const bool isMCmatched_Pi2 ){
	ROOT::Math::PtEtaPhiMVector P4_Gen(0., 0., 0., mPion);

   for (UInt_t g = 0; g < nGenPart; g++){
      if (abs(GenPart_pdgId[g]) ==  ptlID){ 

			P4_Gen.SetPt(GenPart_pt[g]); P4_Gen.SetEta(GenPart_eta[g]); P4_Gen.SetPhi(GenPart_phi[g]);
			if (!isMCmatched_Pi1){ 
				//DRhisto->Fill(ROOT::Math::VectorUtil::DeltaR(P4_Gen, P4_Reco_Pi1));
				if (ROOT::Math::VectorUtil::DeltaR(P4_Gen, P4_Reco_Pi1) < 0.05 ){
					std::cout << "RECO --> NON matching pi1\t" << "pT " << P4_Reco_Pi1.Pt() <<  "\teta " << P4_Reco_Pi1.Eta() <<  "\tphi " << P4_Reco_Pi1.Phi()<< "idx " << B0_pi1_idx[Bidx] << std::endl;
					std::cout << "GEN --> PTL mother " << GenPart_pdgId[GenPart_genPartIdxMother[g]] << "--> " << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[g]]]<< "\tpT " << GenPart_pt[g] << "\teta " << GenPart_eta[g] << "\tphi " << GenPart_phi[g] << std::endl;
					std::cout << "DR(pi1,.)" << ROOT::Math::VectorUtil::DeltaR(P4_Gen, P4_Reco_Pi1) << "\t vs DRmin  pi- " << MCmatch_Pim_DRmin << " idx " << MCmatch_Pim_Idx << "\t DRmin  pi+ "<< MCmatch_Pip_DRmin << " idx " << MCmatch_Pip_Idx<< std::endl;
					std::cout << "GEN mio pim " << GenPimP4.Pt() << " pip " << GenPipP4.Pt() << std::endl; 
				} 
			}
			if (!isMCmatched_Pi2){
				//DRhisto->Fill(ROOT::Math::VectorUtil::DeltaR(P4_Gen, P4_Reco_Pi2));
				if (ROOT::Math::VectorUtil::DeltaR(P4_Gen, P4_Reco_Pi2) < 0.05 ){
					std::cout << "RECO --> NON matching pi2\t" << "pT " << P4_Reco_Pi2.Pt() <<  "\teta " << P4_Reco_Pi2.Eta() <<  "\tphi " << P4_Reco_Pi2.Phi() << std::endl;
					std::cout << "GEN --> PTL mother " << GenPart_pdgId[GenPart_genPartIdxMother[g]]<< "--> " << GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother    [g]]] << "\tpT " << GenPart_pt[g] << "\teta " << GenPart_eta[g] << "\tphi " << GenPart_phi[g] << std::endl;
					std::cout << "DR(pi2,.)" << ROOT::Math::VectorUtil::DeltaR(P4_Gen, P4_Reco_Pi2) << "\t vs DRmin pi- " << MCmatch_Pim_DRmin << " idx " << MCmatch_Pim_Idx <<     "\t DRmin  pi+ "<< MCmatch_Pip_DRmin << " idx " << MCmatch_Pip_Idx<< std::endl; 
					std::cout << "GEN mio pim " << GenPimP4.Pt() << " pip " << GenPipP4.Pt() << std::endl; 
				}
			}
		}
	}
}//PrintGenLevel()


void TriggerSelection::DebugRho(){
		
		ROOT::Math::PtEtaPhiMVector P4_Rho(0., 0., 0., 0.);
		
		std::cout << "--- GENERATOR ---"<< std::endl;
		std::cout << "Pi - == pt " << GenPimP4.Pt() << " eta " << GenPimP4.Eta() << " mass " << GenPimP4.M() << std::endl;
		std::cout << "Pi + == pt " << GenPipP4.Pt() << " eta " << GenPipP4.Eta() << " mass " << GenPipP4.M() << std::endl;
		std::cout << "Rho  == pt " << GenRhoP4.Pt() << " eta " << GenRhoP4.Eta() << " mass " << GenRhoP4.M() << std::endl;
		std::cout << "--- RECO LEVEL ---"  << std::endl;
		std::cout << "Pi 1 == pt " << P4_Reco_Pi1.Pt() << " eta " << P4_Reco_Pi1.Eta() << " mass " << P4_Reco_Pi1.M() << std::endl;
		std::cout << "Pi 2 == pt " << P4_Reco_Pi2.Pt() << " eta " << P4_Reco_Pi2.Eta() << " mass " << P4_Reco_Pi2.M() << std::endl;
		std::cout << "DRmin con Pi- " << MCmatch_Pim_DRmin << "\t  DRmin con Pi+ " << MCmatch_Pip_DRmin << std::endl;
		P4_Rho = P4_Reco_Pi1 + P4_Reco_Pi2; 
		std::cout << "Rho == pt " << P4_Rho.Pt() << " eta " << P4_Rho.Eta() << " mass " << P4_Rho.M() << std::endl;




} // DebugRho()
