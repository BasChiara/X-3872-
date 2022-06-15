#include "./src/B0toX3872K0s.C"
#include "./src/PreSelDATA2017.C"
#include "./src/PrePostFit.C"
#include "./src/TriggerSelection.C"
#include "./src/HLTapply.C"
#include "./src/SGNvsBKGvariables.C"

#include <iostream>
#include <vector>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <dirent.h>
#include <sys/stat.h>
using namespace std;

std::string MCdata2017 = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/B0ToXKs_2022Apr29/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_BdToX3872Ks/220429_084035/";
std::string MCcentral2017 = "root://xrootd-cms.infn.it//store/user/crovelli/BdToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi_TuneCP5_13TeV-pythia8-evtgen/crab_central_signal_2017/BdToX3872Ks_X3872ToJPsiRho_JPsiToMuMu_RhoToPiPi_2017.root";
std::string Data2017_B = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/Xdata2017_2022May17/Charmonium/crab_data_Run2017B/";
std::string Data2017_C = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017C/Run2017C_";
std::string Data2017_D = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_";
std::string Data2017_E = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017E/Run2017E_";

int LxP_DirNum_MAP(const TString dataset){
	std::map <TString , int> LxP_DirNum{};
	
	LxP_DirNum["MC"] = 1;
	LxP_DirNum["data17B"] = 2;
	
	return LxP_DirNum[dataset];
}

int T2_DirNum_MAP(const TString dataset){
	std::map <TString , int> T2_DirNum{};
	
	T2_DirNum["signal17"] = 1;
	T2_DirNum["data17C"]  = 3;
	T2_DirNum["data17D"]  = 2;
	T2_DirNum["data17E"]  = 3;
	
	return T2_DirNum[dataset];
}


std::string Data2017_B_0001 = "/eos/cms/store/group/phys_bphys/crovelli/nanoaod_X/Xdata2017_2022May17/Charmonium/crab_data_Run2017B/220517_110819/0001";
std::string Data2017_D_0000 = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_0000.root";
std::string Data2017_D_0001 = "root://xrootd-cms.infn.it//store/user/crovelli/Charmonium/crab_data_Run2017D/Run2017D_0001.root";

TChain* LoadTree(TString dataset = "MC"){

	//bool LxP = false, T2 = false;	
	
	std::vector <string> LxPlusDirPath;
	std::vector <string> T2DirName; 
	string LxPlusDirName= "";
	std::vector<int> LxP_Ndir;
	std::vector<int> T2_Ndir;

	if (dataset ==    "MC"  ){
		LxPlusDirPath.push_back(MCdata2017); 
		LxP_Ndir.push_back(LxP_DirNum_MAP("MC"));
	}
	if (dataset == "data17B" or dataset == "data17"){ 
		LxPlusDirPath.push_back(Data2017_B); 
		LxP_Ndir.push_back(LxP_DirNum_MAP("data17B"));
	}
	if (dataset == "signal17"){
		T2DirName.push_back(MCcentral2017); 
		T2_Ndir.push_back(T2_DirNum_MAP("signal17")); 
	}
	if (dataset == "data17C" or dataset == "data17"){
		T2DirName.push_back(Data2017_C); 
		T2_Ndir.push_back(T2_DirNum_MAP("data17C")); 
	}
	if (dataset == "data17D" or dataset == "data17"){
		T2DirName.push_back(Data2017_D); 
		T2_Ndir.push_back(T2_DirNum_MAP("data17D")); 
	}
	if (dataset == "data17E" or dataset == "data17"){
		T2DirName.push_back(Data2017_E); 
		T2_Ndir.push_back(T2_DirNum_MAP("data17E")); 
	}
	

	int Nfiles = 0;
	TString tree_path, tree_name = "/Events";
	TChain* chain =new TChain("Events");

	
	for (unsigned int s = 0; s < LxPlusDirPath.size(); s++){

		struct dirent* file = NULL; // entry in the directory
		struct stat file_stats;
		const char* directory;
		DIR* dir_pointer = NULL;

		if (dataset == "MC") std::cout << " ===== READING MC data ===== " << std::endl;;
		if (dataset == "data17B" or dataset == "data17") std::cout << " ===== READING RUN2 2017(B) data ===== " << std::endl;
		
		for (int d = 0; d < LxP_Ndir[s]; d++ ){
			LxPlusDirName = LxPlusDirPath[s] + Form("%.4d", d);
			std::cout << LxPlusDirName << std::endl;
			directory = LxPlusDirName.c_str();
			dir_pointer = opendir(directory);//point to the directory

			while((file = readdir (dir_pointer))){
				if(file == NULL){
					std::cout << "ERROR null pointer to file" << std::endl;
					return NULL;
				}

				if (strcmp(file->d_name, "xNANO_") < 0) continue; // skip "." and ".." and "log"

				Nfiles ++;
				//std::cout << file->d_name << std::endl;
				tree_path = LxPlusDirName + "/" + file->d_name + tree_name; 
				chain->Add(tree_path);
			}
		}
	}
	if (dataset == "data17") cout<<" Number of events: " <<chain->GetEntries()<<std::endl;
	//================ LOADING FILES FROM T2 
	for (unsigned int s = 0; s < T2DirName.size(); s++){

		if (dataset == "data17D" or dataset == "data17")	std::cout << " ===== READING RUN2 2017(D) data ===== " << std::endl;
		if (dataset == "signal17"){
			tree_path =  T2DirName[s] + tree_name;
			chain->Add(tree_path);	
			std::cout << tree_path  << std::endl;
		}else{

			for (int d = 0; d < T2_Ndir[s]; d++ ){
				tree_path = T2DirName[s] + Form("%.4d", d) + ".root" + tree_name; 
				chain->Add(tree_path);	
				std::cout << tree_path  << std::endl;
			}
		}
	}
	std::cout << " ... LOADING " << Nfiles << " FILES ..." << std::endl;
	cout<<" Number of events: " <<chain->GetEntries()<<std::endl;

	return chain;

}//LoadTree()


int JustCheck(TString dataset = "signal17"){

	
	TChain* chain = LoadTree(dataset);
	return 0;

}


int RunPrePostFitAnalysis() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	PrePostFit tree( chain );

	tree.Loop();

	return 0;
}// RunPrePostFitAnalysis()


int RunTriggerSel() {

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	TriggerSelection tree( chain );

	tree.Loop();

	return 0;
}// RunTriggerSel()

int RunHLTapply(TString dataset ="MC") {

	TChain* chain = LoadTree(dataset);

	//================ Run analysis                                                                                                                         
	HLTapply tree( chain , dataset);

	tree.Loop();

	return 0;
}// RunHLTapply()

int RunSgnBkgAnalysis(){ 

	TChain* chain = LoadTree();

	//================ Run analysis                                                                                                                         
	SGNvsBKGvariables tree( chain );

	tree.Loop();

	return 0;
}//RunSgnBkgAnalysis() 
