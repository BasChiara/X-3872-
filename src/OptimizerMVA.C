#include "../include/OptimizerMVA.h"

// RooFit
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooCurve.h"
using namespace RooFit;

OptimizerMVA::OptimizerMVA(){

	// SGN scale factor
	SGNfactor = 0.173451;

	// I/O files

	inFileSGN_path  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGN_MC.root";
	SGNtreeName     = "mcSIGNAL";
	outFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonSGN.root";

	outFileBKG_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonBKG.root";
	inFileBKG_path  = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/data/merged/BKG_MergData17.root";
	BKGtreeName     = "B0sidebands";

	inFileFitPar_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/SGNfit/B0params.txt";

	LoadFitRanges();

	SetUpInTreeSGN();
	SetUpInTreeBKG();

	// BDT model
	BDTweightPath_ = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/datasetBDT/weights/TMVAClassification_BDT_nT200_S25_D4_b03_nC35.weights.xml";


}//OptimizerMVA()

OptimizerMVA::~OptimizerMVA(){ 

	InFileSGN->Close();
	InFileBKG->Close();
}

//    ==== UTILITY ====     //
TFile* OptimizerMVA::GetFile(TString file){

   TFile *File(0);
	TString Fpath("");
	if (file == "InS") Fpath = inFileSGN_path; 
	if (file == "InB") Fpath = inFileBKG_path;
	if (file == "xS")  Fpath = outFileSGN_path;
	if (file == "xB")  Fpath = outFileBKG_path;

	// Open SIGNAL file
   if (!gSystem->AccessPathName(Fpath)) {
      File = TFile::Open(Fpath); // check if file in local directory exists
   }
   if (!File) {
      std::cout << "ERROR: could not open signal-data file" << std::endl;
      exit(1);
   }

	return File; 

}//GetFile()


void OptimizerMVA::SetUpInTreeSGN(){

   InFileSGN = GetFile("InS");
   // ===== SIGNAL TREE =====
	TreeSGN = (TTree*)InFileSGN->Get(SGNtreeName);

   TreeSGN->SetBranchAddress( "run", &run);
   TreeSGN->SetBranchAddress( "LumiBlock", &LumiBlock);
   TreeSGN->SetBranchAddress( "event", &event);

   TreeSGN->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN->SetBranchAddress( "pT_Rho", &pT_Rho);
	TreeSGN->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN->SetBranchAddress( "M_K0s", &M_K0s);

	 // friend-tree with BDT response
   TFile* outputBDT = GetFile("xS");
   TTree* Tree_xS= (TTree*)outputBDT->Get("TreeBDTx_S");
	Tree_xS->SetBranchAddress( "BDTx", &X);
	TreeSGN->AddFriend("TreeBDTx_S");


}//SetUpInTreeSGN()

void OptimizerMVA::SetUpInTreeBKG(){

	// ===== BACKGROUND TREE =====

   InFileBKG = GetFile("InB");
   TreeBKG = (TTree*)InFileBKG->Get(BKGtreeName);

   TreeBKG->SetBranchAddress( "run", &run);
   TreeBKG->SetBranchAddress( "LumiBlock", &LumiBlock);
   TreeBKG->SetBranchAddress( "event", &event);

   TreeBKG->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeBKG->SetBranchAddress( "SVprob", &SVprob);
   TreeBKG->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeBKG->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeBKG->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeBKG->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeBKG->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeBKG->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeBKG->SetBranchAddress( "M_Rho", &M_Rho);
   TreeBKG->SetBranchAddress( "M_B0", &M_B0);
   TreeBKG->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeBKG->SetBranchAddress( "M_X3872", &M_X3872);
   TreeBKG->SetBranchAddress( "M_K0s", &M_K0s);
	
   TFile* outputBDT = GetFile("xB");
   TTree* Tree_xB= (TTree*)outputBDT->Get("TreeBDTx_B");
	Tree_xB->SetBranchAddress( "BDTx", &X);
	TreeBKG->AddFriend("TreeBDTx_B");



}// SetUpInTreeBKG()



void OptimizerMVA::LoadFitRanges(){
	
	std::string line;
	int Nline = 0;

	char ParName[30];
	double err;
	// --- B0 FIT 
	std::ifstream inFileParB0(inFileFitPar_path);	
	if(!inFileParB0.is_open()) std::cout << "ERROR cannot open " << inFileFitPar_path << std::endl;
	while(!inFileParB0.eof()){

		getline(inFileParB0, line); Nline++;
		//std::cout << Nline << "\t" << line << std::endl;
		if(line.c_str()[0] == '#') continue;
		if(Nline == 3) sscanf(line.c_str(), "%s %lf", ParName, &BKGregNearL);
		if(Nline == 4) sscanf(line.c_str(), "%s %lf", ParName, &BKGregFarL);
		if(Nline == 5) sscanf(line.c_str(), "%s %lf", ParName, &BKGregNearR);
		if(Nline == 6) sscanf(line.c_str(), "%s %lf", ParName, &BKGregFarR);
		if(Nline == 9) sscanf(line.c_str(), "%s %lf", ParName, &SGNregL);
		if(Nline ==10) sscanf(line.c_str(), "%s %lf", ParName, &SGNregR);
	}
	inFileParB0.close();


}//LoadFitRanges()

int OptimizerMVA::GetSignalPreCut(){
	
	int S_precut = 0;
	int Nentries = TreeSGN->GetEntries();
	bool inSGNregion;

	for(int ien = 0; ien < Nentries; ien++){
		TreeSGN->GetEntry(ien);
		inSGNregion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if(inSGNregion)S_precut++;
	}
	
	return S_precut;
	
}//GetSignalPreCut()

void OptimizerMVA::GetCorrelation_BDT_MR(){

	double CorrF;
	const int NbinsX = 100; 
	const double Xlow = -0.8, Xhigh = 0.5;
	const int NbinsMR = 100; 
	const double MRlow = 0.4, MRhigh = .9;
	TH2F* h_CorrS = new TH2F("CorrS_XvsMrho", "", NbinsMR, MRlow, MRhigh, NbinsX, Xlow, Xhigh);	
	TH2F* h_CorrB = new TH2F("CorrB_XvsMrho", "", NbinsMR, MRlow, MRhigh, NbinsX, Xlow, Xhigh);	
	
	TreeSGN->Draw("TreeBDTx_S.BDTx:M_Rho>>CorrS_XvsMrho");
	h_CorrS->GetXaxis()->SetTitle("M(\\rho(770))\\ [GeV]");
	h_CorrS->GetYaxis()->SetTitle("BDT response [x]");
	h_CorrS->GetYaxis()->SetTitleOffset(1.05);
	h_CorrS->SetLineColor(kAzure +1);
	h_CorrS->SetFillColorAlpha(kAzure +1, 0.60);
	TreeBKG->Draw("TreeBDTx_B.BDTx:M_Rho>>CorrB_XvsMrho");
	h_CorrB->SetLineColor(kRed);
	h_CorrB->SetFillColorAlpha(kRed, 0.50);


	auto legend= new TLegend(0.525, 0.8,.89,.89);
	legend->SetTextSize(0.025);
	legend->SetBorderSize(0);
	
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);	
	gStyle->SetOptStat(0);
	h_CorrS->Draw("BOX");
	legend->AddEntry(h_CorrS, Form("SIGNAL (MC) Corr = %.3f", h_CorrS->GetCorrelationFactor()));
	h_CorrB->Draw("BOX SAME");
	legend->AddEntry(h_CorrB, Form("DATA 2017   Corr = %.3f", h_CorrB->GetCorrelationFactor()));
	legend->Draw();

	TString outPath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/";
	c1->SaveAs(outPath + "CorrelationBDTxMRho.png");

}

double OptimizerMVA::GetAvMultiplicity(const double Xcut, const double MRcut, TString dataset = "SGN"){
	
	TTree* inTree;
	if (dataset == "SGN") inTree = TreeSGN;
	if (dataset == "BKG") inTree = TreeBKG;

	double AvMultiplicity;
	int Nentries = inTree->GetEntries(), Ncand = 0, Nevents = 0;
	float prevRun = -1, prevLumiBlock = -1, prevEvent = -1;
	bool isNewEvent, passCut;

	for(int ientry = 0; ientry < Nentries; ientry++){
		inTree->GetEntry(ientry);
		passCut = (M_Rho > MRcut) && (X > Xcut );
		if (!passCut) continue;
		Ncand++;
		isNewEvent = !((run == prevRun) && (LumiBlock == prevLumiBlock ) && (event == prevEvent));
		prevRun = run; prevLumiBlock = LumiBlock; prevEvent = event;
		if ( isNewEvent ) Nevents++;

	}
	std::cout << "#cand vs #events\t" << Ncand << " " << Nevents << std::endl;

	return ((double)Ncand/Nevents);

}//GetAvMultiplicitySGN()



// ==== ANALYSIS TOOLS ==== //

int OptimizerMVA::ApplyBDT(){

// --> TMVA SET UP
   // This loads the library
   TMVA::Tools::Instance();
   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;
 
 
   // Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
 
   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t pTM_B0, SVprob, LxySign_B0, CosAlpha_B0;  
   Float_t DR_Pi1B0, pT_Pi1, pT_Rho, D0_Rho; 
   reader->AddVariable( "pTM_B0", &pTM_B0);
   reader->AddVariable( "SVprob", &SVprob);
   reader->AddVariable( "LxySign_B0", &LxySign_B0);
   reader->AddVariable( "CosAlpha_B0", &CosAlpha_B0);
   reader->AddVariable( "DR_Pi1B0", &DR_Pi1B0);
   reader->AddVariable( "pT_Pi1", &pT_Pi1);
   reader->AddVariable( "pT_Rho", &pT_Rho);
   reader->AddVariable( "D0_Rho", &D0_Rho);
 
   // Spectator variables declared in the training have to be added to the reader, too
   Float_t M_Rho; 
   reader->AddSpectator( "M_Rho", &M_Rho);
 
   // Book the BDT method
   TString methodName = TString("BDT method"); 
	reader->BookMVA( methodName, BDTweightPath_);

// --> SIGNAL & BACKGROUND files reading
	TFile* input_signal = GetFile("InS");
	std::cout << "     TMVAClassificationApp    : Using input file (SIGNAL): " << input_signal->GetName() << std::endl;
 
   TFile *input_background = GetFile("InB");
	std::cout << "     TMVAClassificationApp    : Using input file (BACKGROUND): " << input_background->GetName() << std::endl;

// --> EVENT LOOP 
 

	float M_B0, M_MuMu, M_X3872, M_K0s;
   // ===== SIGNAL TREE =====
   std::cout << "--- Select signal sample" << std::endl;
   TTree* TreeSGN = (TTree*)input_signal->Get(SGNtreeName);

   TreeSGN->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeSGN->SetBranchAddress( "SVprob", &SVprob);
   TreeSGN->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeSGN->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeSGN->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeSGN->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeSGN->SetBranchAddress( "pT_Rho", &pT_Rho);
	TreeSGN->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeSGN->SetBranchAddress( "M_Rho", &M_Rho);
   TreeSGN->SetBranchAddress( "M_B0", &M_B0);
   TreeSGN->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeSGN->SetBranchAddress( "M_X3872", &M_X3872);
   TreeSGN->SetBranchAddress( "M_K0s", &M_K0s);
	 // friend-tree with BDT response
	Float_t Xs;
	TString outFileSGN_path = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/BDTonSGN.root";
	TFile* output_signal = new TFile(outFileSGN_path, "RECREATE"); 
	TTree* TreeBDTx_S = new TTree("TreeBDTx_S", "friend with BDT response"); 
	TreeBDTx_S->Branch("BDTx", &Xs, "BDTx/F");


	// ===== BACKGROUND TREE =====

   std::cout << "--- Select background sample" << std::endl;
	TreeBKG = (TTree*)input_background->Get(BKGtreeName);

   TreeBKG->SetBranchAddress( "pTM_B0", &pTM_B0);
   TreeBKG->SetBranchAddress( "SVprob", &SVprob);
   TreeBKG->SetBranchAddress( "LxySign_B0", &LxySign_B0);
   TreeBKG->SetBranchAddress( "CosAlpha_B0", &CosAlpha_B0);
   TreeBKG->SetBranchAddress( "DR_Pi1B0", &DR_Pi1B0);
   TreeBKG->SetBranchAddress( "pT_Pi1", &pT_Pi1);
   TreeBKG->SetBranchAddress( "pT_Rho", &pT_Rho);
   TreeBKG->SetBranchAddress( "D0_Rho", &D0_Rho);
   TreeBKG->SetBranchAddress( "M_Rho", &M_Rho);
   TreeBKG->SetBranchAddress( "M_B0", &M_B0);
   TreeBKG->SetBranchAddress( "M_mumu", &M_MuMu);
   TreeBKG->SetBranchAddress( "M_X3872", &M_X3872);
   TreeBKG->SetBranchAddress( "M_K0s", &M_K0s);
	 // friend-tree with BDT response
	Float_t Xb;
	TFile* output_background = new TFile(outFileBKG_path, "RECREATE"); 
	TTree* TreeBDTx_B = new TTree("TreeBDTx_B", "friend with BDT response"); 
	TreeBDTx_B->Branch("BDTx", &Xb, "BDTx/F");
	


	// ---> SIGNAL evaluation 
   std::cout << "--- Processing: " << TreeSGN->GetEntries() << " signal events" << std::endl;
   TStopwatch sw;
   sw.Start();

	for (Long64_t ievt=0; ievt<TreeSGN->GetEntries();ievt++) {
 
		if (ievt%100 == 0) std::cout << "    ... Processing signal-event: " << ievt << std::endl;
		TreeSGN->GetEntry(ievt);
		Xs = reader->EvaluateMVA(methodName);

		TreeBDTx_S->Fill(); 
   }

	// ---> BACKGROUND evaluation 
   std::cout << "--- Processing: " << TreeBKG->GetEntries() << " background events" << std::endl;

	for (Long64_t ievt=0; ievt<TreeBKG->GetEntries();ievt++) {
 
		if (ievt%1000 == 0) std::cout << "    ... Processing background-event: " << ievt << std::endl;
		TreeBKG->GetEntry(ievt);
		Xb = reader->EvaluateMVA(methodName);

		TreeBDTx_B->Fill();
   }


   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();
	
	output_signal->cd();
	TreeBDTx_S->Write();

	output_background->cd();
	TreeBDTx_B->Write();

	// close files
	output_signal->Close();
	output_background->Close();
	input_signal->Close();
	input_background->Close();	



	return 0;

} //ApplyBDT()



double OptimizerMVA::BKG_NevExtraction(const double Xcut, const double MRcut){

	double Bs;

	// ==== ROOFIT SET UP ==== //
	
	RooRealVar B0m("B0m", "M(B_0)\\ [GeV]", 5.07947, 5.47934);

	B0m.setRange(   "SB1"      , BKGregFarL , BKGregNearL);       // SIDEBAND 1
	B0m.setRange(   "SB2"      , BKGregNearR, BKGregFarR);        // SIDEBAND 2
	B0m.setRange("BLINDregion" , BKGregNearL, BKGregNearR);        // SIDEBAND 2
	B0m.setRange("FULLregion"  , BKGregFarL , BKGregFarR);        // FULLregion RANGE 
	B0m.setRange("SGNregion"   , SGNregL    , SGNregR);           // SGNregion RANGE 

	RooDataSet B0m_ds("B0m_ds","B0m_ds", RooArgSet(B0m));

	// Expo + C	
	RooRealVar tau("tau", "", -15.0, -30.0 , -5.);
	RooExponential expo("expo", "", B0m, tau);
	RooProdPdf BKGdesc("BKGdesc", "BKGdesc", expo); 

	RooPolynomial Const("Const", "", B0m);
	RooProdPdf BKGbase("BKGbase", "BKGbase", Const); 

	
	// Fermi 
	RooRealVar Slope("Slope", "", 10.0, 1.0 , 100.);
	RooRealVar Flex("Flex", "", 4., 1. , 5.2);
	RooRealVar C("C", "", 0.5, 0. , 1.);
	
	RooGenericPdf Fermi("Fermi", "", "1./(1. + exp(( @0 - @1)*@2)) + @3", RooArgList(B0m, Flex, Slope, C));

	// Build PDF
	RooRealVar bkg_yield("nBKG", "", 1600, 100, 10000);
	RooRealVar f("f", "", 0., 1.);
	RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", Fermi, bkg_yield);
	//RooAddPdf BKGmodel("BKGmodel", "B0 sidebans shape", RooArgList(expo, Const), f);


	// ==== LOOP ON EVENTS ==== // 
	int N = TreeBKG->GetEntries();
	for (int i = 0; i < N; i++){

		TreeBKG->GetEntry(i);
		if( ( X < Xcut) || ( M_Rho < MRcut) ) continue;
		B0m = M_B0;
		B0m_ds.add(RooArgSet(B0m));

	}
	B0m_ds.Print("v");
	// ==== FIT PERFORM ==== //
	
	RooFitResult *ResBKGmodelSB2 = BKGmodel.fitTo(B0m_ds, Range("SB2"), Save());
	//ResBKGmodelSB2->Print();
	RooFitResult *ResBKGmodelSB1 = BKGmodel.fitTo(B0m_ds, Range("SB1"), Save());
	//ResBKGmodelSB1->Print();
	RooFitResult *ResBKGmodel= BKGmodel.fitTo(B0m_ds, Range("SB1,SB2"), Save());
	//ResBKGmodel->Print("v");

	// ==== SAVE RESULTS ==== //
	RooPlot* B0m_fr = B0m.frame(Bins(47));
	B0m_ds.plotOn(B0m_fr);
	B0m_fr->setInvisible("BLINDregion");
	BKGmodel.plotOn(B0m_fr, Range("FULLregion"), RooFit::NormRange("SB1,SB2"));
	BKGmodel.paramOn(B0m_fr, Layout(0.60));
	B0m_fr->SetTitle(Form("X-cut = %.2f Mrho-cut = %.3f", Xcut, MRcut));

	
	TText *TxtChi2= new TText(5.2, B0m_fr->GetMaximum()*0.85, Form("Chi^2 = %.3f", B0m_fr->chiSquare()));
   TxtChi2->SetTextSize(0.035);
   TxtChi2->SetTextColor(kRed);
   B0m_fr->addObject(TxtChi2);
	//std::cout << " ---> Chi^2 = " << B0m_fr->chiSquare() << std::endl;

	TH2 *h_Corr = ResBKGmodel->correlationHist();
	
	// ==== Bs ==== //
	RooAbsReal* IntSreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SGNregion")); //integrate sgn region
	double Is = IntSreg->getVal();
	RooAbsReal* IntBreg = BKGmodel.createIntegral(B0m, NormSet(B0m), Range("SB1,SB2"));   //integrate bkg region
	double Ib = IntBreg->getVal();

	Bs = B0m_ds.sumEntries() * Is/Ib; 

	// ==== SAVE ON FILE ==== //
	TString outFileName("/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/FitRes/FitBKG_X");
	if (Xcut < 0.) outFileName.Append(Form("m"));
	outFileName.Append(Form("%.0f_MR%.0f.root", fabs(Xcut)*100., MRcut*1000.));

	TFile* outFitFile = new TFile(outFileName, "RECREATE");
	B0m_fr->Write("FitPlot");
	ResBKGmodel->Write("FitResults");
	h_Corr->Write();

	outFitFile->Close();
	

	return Bs;

}//BKG_NevExtraction()




double OptimizerMVA::SGN_NevExtraction(const double Xcut, const double MRcut){

	int iSs = 0;
	double Ss;
	bool inSGNregion = false;
	int Nentries = TreeSGN->GetEntries();

	for(int ien = 0; ien < Nentries; ien++){
		TreeSGN->GetEntry(ien);
		inSGNregion = (M_B0 > SGNregL) && (M_B0 < SGNregR);
		if( (X < Xcut) || (M_Rho < MRcut)  || (!inSGNregion) ) continue;
		iSs ++;
	}
	Ss = iSs * SGNfactor;
	return Ss;

}//SGN_NevExtraction()
	



double OptimizerMVA::PunziSign(const double S, const double B, double* PSerr){

	double PunziSign, Sign_denom, Sign_denom_error;
	double SGNeff, SGNeff_error;
	const double b = 5.0;    // #sigmas corresp. to 5sigma significance level
	const double a = 2.0; // #sigmas corresp. to CL (90%--> 1.2816) (95% --> 1.6448) (CMStwiki --> 2.)

	SGNeff       = S / (GetSignalPreCut() * SGNfactor);
	SGNeff_error = SGNfactor / S; // error square
	
	Sign_denom        = b*b + 2.*a*sqrt(B) + b*sqrt(b*b + 4.*a*sqrt(B) + 4.*B );
	Sign_denom_error  = a/sqrt(B) + b * (a / sqrt(B) + 2. ) / sqrt(b*b + 4.*a*sqrt(B) + 4.*B );	
	Sign_denom_error /= Sign_denom;	

	PunziSign = SGNeff/Sign_denom;
	*PSerr = PunziSign*sqrt(SGNeff_error + Sign_denom_error*Sign_denom_error);


	return PunziSign;

}



int OptimizerMVA::makeSGNvsBKGplot(const double Xcut, const double MRcut){
	int Nbins = 70;
	double Mlow = 5.0, Mhigh = 5.6;
	TH1F* h_Data_B0 = new TH1F("Data_B0", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_B0  = new TH1F("SGN_B0" , "", Nbins, Mlow, Mhigh); 
	Nbins = 30;
	Mlow = 3.8, Mhigh = 3.95;
	TH1F* h_Data_X  = new TH1F("Data_X", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_X   = new TH1F("SGN_X" , "", Nbins, Mlow, Mhigh); 
	Nbins = 50;
	Mlow = 0.45, Mhigh = 0.55;
	TH1F* h_Data_K0s  = new TH1F("Data_K0s", "", Nbins, Mlow, Mhigh);
	TH1F* h_SGN_K0s   = new TH1F("SGN_K0s" , "", Nbins, Mlow, Mhigh); 

	// ==== LOOP ON EVENTS ==== // 
	// ...BKG
	int Nb = TreeBKG->GetEntries();
	for (int i = 0; i < Nb; i++){
		TreeBKG->GetEntry(i);
		if( ( X < Xcut) || ( M_Rho < MRcut) ) continue;
		h_Data_B0->Fill(M_B0);
		h_Data_X->Fill(M_X3872);
		h_Data_K0s->Fill(M_K0s);
	}
	int Ns = TreeSGN->GetEntries();
	for (int i = 0; i < Ns; i++){
		TreeSGN->GetEntry(i);
		if( ( X < Xcut) || ( M_Rho < MRcut) ) continue;
		h_SGN_B0->Fill(M_B0);
		h_SGN_X->Fill(M_X3872);
		h_SGN_K0s->Fill(M_K0s);
	}
	// Histo set-up
	// ---> B0
	h_SGN_B0->SetTitle(Form("X %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_B0->GetXaxis()->SetTitle("M(B_0)\\ [GeV]");
	h_SGN_B0->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_B0->GetXaxis()->GetBinWidth(1)));
	h_SGN_B0->SetLineWidth(3);
	h_SGN_B0->SetLineColor(kAzure + 1); h_SGN_B0->SetFillColorAlpha(kAzure + 1, 0.30);
	h_Data_B0->SetLineColor(kBlack);
	h_Data_B0->SetMarkerStyle(20);
	h_Data_B0->SetLineWidth(2);
	// ---> X(3872)
	h_SGN_X->SetTitle(Form("X %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_X->GetXaxis()->SetTitle("M(X)\\ [GeV]");
	h_SGN_X->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_X->GetXaxis()->GetBinWidth(1)));
	h_SGN_X->SetLineWidth(3);
	h_SGN_X->SetLineColor(kViolet + 1); h_SGN_X->SetFillColorAlpha(kViolet + 1, 0.30);
	h_Data_X->SetLineColor(kBlack);
	h_Data_X->SetMarkerStyle(20);
	h_Data_X->SetLineWidth(2);
	// ---> K0short 
	h_SGN_K0s->SetTitle(Form("K0s %.2f  Mrho %.3f", Xcut, MRcut));
	h_SGN_K0s->GetXaxis()->SetTitle("M(K_0^s)\\ [GeV]");
	h_SGN_K0s->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_SGN_K0s->GetXaxis()->GetBinWidth(1)));
	h_SGN_K0s->SetLineWidth(3);
	h_SGN_K0s->SetLineColor(kGreen + 1); h_SGN_K0s->SetFillColorAlpha(kGreen + 1, 0.30);
	h_Data_K0s->SetLineColor(kBlack);
	h_Data_K0s->SetMarkerStyle(20);
	h_Data_K0s->SetLineWidth(2);

	// get FIT
	TString inFileName = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/FitRes/FitBKG_X";
	if (Xcut < 0.) inFileName.Append(Form("m"));
	inFileName.Append(Form("%.0f_MR%.0f.root", fabs(Xcut)*100., MRcut*1000.));
	TFile* inFile = new TFile(inFileName);
	RooPlot* FitPlot = (RooPlot*)inFile->Get("FitPlot");
	RooCurve* FitCurve = (RooCurve*)FitPlot->getObject(1);
	FitCurve->SetLineColor(kRed);

	// Legend 
	auto legendB0 = new TLegend(0.60, 0.75,.89,.89);
	legendB0->SetTextSize(0.025);
	legendB0->SetBorderSize(0);
	auto legendX= new TLegend(0.60, 0.80,.89,.89);
	legendX->SetTextSize(0.025);
	legendX->SetBorderSize(0);
	auto legendK0s= new TLegend(0.60, 0.80,.89,.89);
	legendK0s->SetTextSize(0.025);
	legendK0s->SetBorderSize(0);

	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	h_SGN_B0->Scale(SGNfactor);
	h_SGN_B0->GetYaxis()->SetRangeUser(0., 1.2 * std::max(h_SGN_B0->GetMaximum(), h_Data_B0->GetMaximum()));
	h_SGN_B0->Draw("HIST");
	legendB0->AddEntry(h_SGN_B0, "SIGNAL (MC)");
	h_Data_B0->Draw("PE0 SAME");
	legendB0->AddEntry(h_Data_B0, "DATA 2017");
	FitCurve->Draw("SAME");
	legendB0->AddEntry(FitCurve, "SIDEBANDS-FIT");

	legendB0->Draw();
	TString outPath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/";
	c1->SaveAs(outPath + "B0massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	h_SGN_X->Scale(SGNfactor);
	h_SGN_X->GetYaxis()->SetRangeUser(0.,1.2 * std::max(h_SGN_X->GetMaximum(), h_Data_X->GetMaximum()));
	h_SGN_X->Draw("HIST");
	legendX->AddEntry(h_SGN_X, "SIGNAL (MC)");
	h_Data_X->Draw("PE0 SAME");
	legendX->AddEntry(h_Data_X, "DATA 2017");
	legendX->Draw();
	c1->SaveAs(outPath + "X3872massCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	h_SGN_K0s->Scale(SGNfactor);
	h_SGN_K0s->GetYaxis()->SetRangeUser(0.,1.2 * std::max(h_SGN_K0s->GetMaximum(), h_Data_K0s->GetMaximum()));
	h_SGN_K0s->Draw("HIST");
	legendK0s->AddEntry(h_SGN_K0s, "SIGNAL (MC)");
	h_Data_K0s->Draw("PE0 SAME");
	legendK0s->AddEntry(h_Data_K0s, "DATA 2017");
	legendK0s->Draw();
	c1->SaveAs(outPath + "K0smassCUT_" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	inFile->Close();
	delete h_Data_B0;
	delete h_SGN_B0;
	delete h_Data_X;
	delete h_SGN_X;
	delete h_Data_K0s;
	delete h_SGN_K0s;

	return 0;
}


int OptimizerMVA::makeMASSplot2D(const double Xcut, const double MRcut){

	int Nbins = 50;
	double Mlow_B0 = 5.20 , Mhigh_B0 = 5.40;
	double Mlow_X  = 3.80 , Mhigh_X  = 3.95;
	double Mlow_K0s= 0.45 , Mhigh_K0s= 0.55; 
	TH2F* h_XvsB0_m       = new TH2F( "XvsB0_m"     , "", Nbins, Mlow_B0 , Mhigh_B0,  Nbins, Mlow_X  , Mhigh_X  );
	TH2F* h_XvsB0_m_bkg   = new TH2F( "XvsB0_m_bkg" , "", Nbins, Mlow_B0 , Mhigh_B0,  Nbins, Mlow_X  , Mhigh_X  );
	TH2F* h_K0svsB0_m     = new TH2F("K0svsB0_m"    , "", Nbins, Mlow_B0 , Mhigh_B0,   Nbins, Mlow_K0s, Mhigh_K0s); 
	TH2F* h_K0svsB0_m_bkg = new TH2F("K0svsB0_m_bkg", "", Nbins, Mlow_B0 , Mhigh_B0,   Nbins, Mlow_K0s, Mhigh_K0s); 
	TH2F* h_XvsK0s_m      = new TH2F("XvsK0s_m"     , "", Nbins, Mlow_K0s, Mhigh_K0s , Nbins, Mlow_X, Mhigh_X); 
	TH2F* h_XvsK0s_m_bkg  = new TH2F("XvsK0s_m_bkg" , "", Nbins, Mlow_K0s, Mhigh_K0s , Nbins, Mlow_X, Mhigh_X); 

	TString CutOpt = Form("TreeBDTx_S.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	TString CutOptB = Form("TreeBDTx_B.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	TreeSGN->Draw("M_X3872:M_B0>>XvsB0_m",  CutOpt); 
	TreeBKG->Draw("M_X3872:M_B0>>XvsB0_m_bkg",  CutOptB); 
	h_XvsB0_m->SetFillColor(kAzure+1);  h_XvsB0_m_bkg->SetFillColor(kBlack);
	h_XvsB0_m->GetXaxis()->SetTitle("M(B_0)\\ [GeV]");
	h_XvsB0_m->GetYaxis()->SetTitle("M(X)\\ [GeV]");
	TreeSGN->Draw("M_K0s:M_B0>>K0svsB0_m", CutOpt); 
	TreeBKG->Draw("M_K0s:M_B0>>K0svsB0_m_bkg", CutOptB); 
	h_K0svsB0_m->SetFillColor(kGreen); h_K0svsB0_m_bkg->SetFillColor(kBlack);
	h_K0svsB0_m->GetXaxis()->SetTitle("M(B_0)\\ [GeV]");
	h_K0svsB0_m->GetYaxis()->SetTitle("M(K_0^s)\\ [GeV]");
	TreeSGN->Draw("M_X3872:M_K0s>>XvsK0s_m", CutOpt); 
	TreeBKG->Draw("M_X3872:M_K0s>>XvsK0s_m_bkg", CutOptB); 
	h_XvsK0s_m->SetFillColor(kViolet);  h_XvsK0s_m_bkg->SetFillColor(kBlack);
	h_XvsK0s_m->GetXaxis()->SetTitle("M(K_0^s)\\ [GeV]");
	h_XvsK0s_m->GetYaxis()->SetTitle("M(X)\\ [GeV]");

	
	TString outPath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/";
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	TPad* pad = new TPad("pad", "", 0.,0.,1., 1.);
	pad->SetLeftMargin(0.15); 
	pad->Draw();
	pad->cd();
	h_XvsB0_m->Draw("BOX"); h_XvsB0_m_bkg->Draw("BOX SAME");
	pad->Update();
	c1->SaveAs(outPath + "M_X3872vsM_B0" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	h_K0svsB0_m->Draw("BOX"); h_K0svsB0_m_bkg->Draw("BOX SAME");
	pad->Update();
	c1->SaveAs(outPath + "M_K0svsM_B0" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");
	h_XvsK0s_m->Draw("BOX"); h_XvsK0s_m_bkg->Draw("BOX SAME");
	pad->Update();
	c1->SaveAs(outPath + "M_X3872vsM_K0s" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	return 0;
}//makeMASSplot2D




void OptimizerMVA::GetControlPlots(const double Xcut, const double MRcut){

	int Nbins = 100;
	double Mlow_Rho = 0.4 , Mhigh_Rho= 0.9;
	double low_X = -0.80 , high_X  = 0.2;

	TH1F* h_MrhoCutS= new TH1F("MrhoCutS", "", Nbins, Mlow_Rho, Mhigh_Rho);
	TH1F* h_MrhoCutB= new TH1F("MrhoCutB", "", Nbins, Mlow_Rho, Mhigh_Rho);
	TH1F* h_BDTxCutS= new TH1F("BDTxCutS", "", Nbins, low_X, high_X);
	TH1F* h_BDTxCutB= new TH1F("BDTxCutB", "", Nbins, low_X, high_X);

	TString CutOptS = Form("TreeBDTx_S.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	TString CutOptB = Form("TreeBDTx_B.BDTx>%f && M_Rho>%f", Xcut, MRcut);
	// Mrho
	TreeSGN->Draw("M_Rho>>MrhoCutS",  CutOptS); 
	TreeBKG->Draw("M_Rho>>MrhoCutB",  CutOptB); 
	h_MrhoCutS->GetXaxis()->SetTitle("M(\\rho(770))\\ GeV");
	h_MrhoCutS->GetYaxis()->SetTitle(Form("Events/%.3f [GeV]", h_MrhoCutS->GetXaxis()->GetBinWidth(1)));
	h_MrhoCutS->SetLineWidth(3);	
	h_MrhoCutS->SetLineColor(kAzure +1); h_MrhoCutS->SetFillColorAlpha(kAzure +1, 0.30);
	h_MrhoCutB->SetLineColor(kRed); h_MrhoCutB->SetFillColorAlpha(kRed, 0.30);
	// X BDT
	TreeSGN->Draw("TreeBDTx_S.BDTx>>BDTxCutS",  CutOptS); 
	TreeBKG->Draw("TreeBDTx_B.BDTx>>BDTxCutB",  CutOptB); 
	h_BDTxCutS->GetXaxis()->SetTitle("BDT response [x]");
	h_BDTxCutS->GetYaxis()->SetTitle(Form("Events/%.3f ", h_BDTxCutS->GetXaxis()->GetBinWidth(1)));
	h_BDTxCutS->SetLineWidth(3);	
	h_BDTxCutS->SetLineColor(kAzure +1); h_BDTxCutS->SetFillColorAlpha(kAzure +1, 0.30);
	h_BDTxCutB->SetLineColor(kRed); h_BDTxCutB->SetFillColorAlpha(kRed, 0.30);


	TString outPath = "/afs/cern.ch/user/c/cbasile/CMSSW-10-6-20-Analysis/src/BParkNANO/B0toX3872K0s/results/CutOptimization/CutOptPlots/";
	TCanvas* c1 = new TCanvas("c1","canvas", 1024, 1024);
	gStyle->SetOptStat(0);
	h_MrhoCutS->Draw();
	h_MrhoCutB->Draw("same");

	c1->SaveAs(outPath + "M_RhoCUT" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");

	h_BDTxCutS->Draw();
	h_BDTxCutB->Draw("same");
	c1->SaveAs(outPath + "BDTxCUT" + Form("X%.0f_MR%.0f", fabs(Xcut)*100., MRcut*1000.) + ".png");



}//GetControlPlots()

