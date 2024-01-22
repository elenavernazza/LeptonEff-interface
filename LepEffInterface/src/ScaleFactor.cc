#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"

// used for legacy /grid_mnt/data__data.polcms/cms/vernazza/FrameworkNanoAOD/hhbbtt-analysis/nanoaod_base_analysis/data/cmssw/CMSSW_12_3_0_pre6/src/HTT-utilities/LepEffInterface/data/Electron/Run2016BtoH/Electron_Ele24_eff.root
void ScaleFactor::init_ScaleFactor(TString inputRootFile) {
	etaIsAbsolute = true;
	TFile * fileIn = new TFile(inputRootFile, "read");
	// if root file not found
	if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from NTupleMaker/src/ScaleFactor.cc : File " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };
	
	std::string HistoBaseName = "ZMass";
	etaBinsH = (TH1D*)fileIn->Get("etaBinsH"); 
	std::string etaLabel, GraphName;
	int nEtaBins = etaBinsH->GetNbinsX();

 	for (int iBin=0; iBin<nEtaBins; iBin++){    
		etaLabel = etaBinsH->GetXaxis()->GetBinLabel(iBin+1);
		GraphName = HistoBaseName+etaLabel+"_Data";

		if (fileIn->GetListOfKeys()->Contains(TString(GraphName))){
			eff_data[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName)); 
			SetAxisBins(eff_data[etaLabel]);
		}
		else eff_data[etaLabel] = 0;

		GraphName = HistoBaseName+etaLabel+"_MC";
		if (fileIn->GetListOfKeys()->Contains(TString(GraphName))){
			eff_mc[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
			SetAxisBins(eff_mc[etaLabel]); 
		}
		else eff_mc[etaLabel] =0;

		if (eff_data[etaLabel] != 0 && eff_mc[etaLabel] != 0 ) {
			bool sameBinning = check_SameBinning(eff_data[etaLabel], eff_mc[etaLabel]);
			if (!sameBinning) {std::cout<< "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc . Can not proceed because ScaleFactor::check_SameBinning returned different pT binning for data and MC for eta label " << etaLabel << std::endl; exit(1); }; 
		}
	}
	
	return;
}


void ScaleFactor::init_EG_ScaleFactor(TString inputRootFile, bool isTrg){

  TFile * fileIn = new TFile(inputRootFile, "read");
  // if root file not found
  if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_EG_ScaleFactor(TString inputRootFile, bool isTrg) from LepEffInterface/src/ScaleFactor.cc : File " << inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); }
  TString hsfname = isTrg ? "SF2D":"EGamma_SF2D";
  TH2F *hSF = (TH2F*)fileIn->Get(hsfname);
  TH2F *heffdat = (TH2F*)fileIn->Get(isTrg ? "eff_data":"EGamma_SF2D");
  TH2F *heffmc = (TH2F*)fileIn->Get(isTrg ? "eff_mc":"EGamma_SF2D");

  // retrieve eta binning (ugly, but should work fine)
  const int nbin_eta = isTrg? hSF->GetNbinsY():hSF->GetNbinsX();
  TString eta_bins[nbin_eta] = {""};
  TString firstbinlabel = Form("EtaLt%.1f",isTrg?heffdat->GetYaxis()->GetBinUpEdge(1):hSF->GetXaxis()->GetBinUpEdge(1));
  TString lastbinlabel  = Form("EtaGt%.1f",isTrg?heffdat->GetYaxis()->GetBinLowEdge(nbin_eta):hSF->GetXaxis()->GetBinLowEdge(nbin_eta));
  firstbinlabel.ReplaceAll(".","p");
  lastbinlabel.ReplaceAll(".","p");

  // Check if there are negative eta values in the inputs
  etaIsAbsolute = (isTrg?heffdat->GetYaxis()->GetBinLowEdge(1):hSF->GetXaxis()->GetBinLowEdge(1)) >= 0.;

  //create etabinning Histo
  etaBinsH = new TH1D("etaBinsH","",nbin_eta, isTrg?heffdat->GetYaxis()->GetXbins()->GetArray():hSF->GetXaxis()->GetXbins()->GetArray());

  TString TetaLabel;
  TString GraphName;
  for (int iBin=1; iBin<=nbin_eta; iBin++) {
    if (iBin==1) {
      TetaLabel = firstbinlabel;
    }
    else if (iBin==nbin_eta) {
      TetaLabel = lastbinlabel;
    }
    else {
      TetaLabel = Form("Eta%.1fto%.1f",
		       isTrg?heffdat->GetYaxis()->GetBinLowEdge(iBin):hSF->GetXaxis()->GetBinLowEdge(iBin),
		       isTrg?heffdat->GetYaxis()->GetBinUpEdge(iBin):hSF->GetXaxis()->GetBinUpEdge(iBin));
      TetaLabel.ReplaceAll(".","p");
    }
    etaBinsH->GetXaxis()->SetBinLabel(iBin,TetaLabel);

    std::string etaLabel = (std::string)TetaLabel;
    //GraphName = TString(HistoBaseName)+"_"+etaLabel+"_Data";

    TH1F *hslice_data;
    TH1F *hslice_mc;
    if (!isTrg) {
      hslice_data = (TH1F*) hSF->ProjectionY("slicedata",iBin,iBin);
      hslice_mc = (TH1F*) hSF->ProjectionY("slicemc",iBin,iBin);
    } else {
      hslice_data = (TH1F*) heffdat->ProjectionX("slicedata",iBin,iBin);
      hslice_mc   = (TH1F*) heffmc->ProjectionX("slicemc",iBin,iBin);

    }
    const int nbin_pt = hslice_data->GetNbinsX();

    double data_pt_nom[nbin_pt] = {0};
    double data_eff_nom[nbin_pt] = {0};
    double data_pt_errlow[nbin_pt] = {0};
    double data_eff_errlow[nbin_pt] = {0};
    double data_pt_errhigh[nbin_pt] = {0};
    double data_eff_errhigh[nbin_pt] = {0};
    double mc_pt_nom[nbin_pt] = {0};
    double mc_eff_nom[nbin_pt] = {0};
    double mc_pt_errlow[nbin_pt] = {0};
    double mc_eff_errlow[nbin_pt] = {0};
    double mc_pt_errhigh[nbin_pt] = {0};
    double mc_eff_errhigh[nbin_pt] = {0};

    for (int iptbin=1; iptbin<=nbin_pt; iptbin++) {
	  //iptbin=0 is the underflow bin, which we do not want to keep. iptbin=nbin_pt+1 is the overflow bin
	  // TGraph uses difference between center value and bin edge as errlow and errhigh
      data_pt_nom[iptbin-1]      = hslice_data->GetXaxis()->GetBinCenter(iptbin);
      data_eff_nom[iptbin-1]     = hslice_data->GetBinContent(iptbin);
      data_pt_errlow[iptbin-1]   = hslice_data->GetXaxis()->GetBinCenter(iptbin) - hslice_data->GetXaxis()->GetBinLowEdge(iptbin);
      data_pt_errhigh[iptbin-1]  = hslice_data->GetXaxis()->GetBinUpEdge(iptbin) - hslice_data->GetXaxis()->GetBinCenter(iptbin);
      data_eff_errlow[iptbin-1]  = hslice_data->GetBinError(iptbin);
      data_eff_errhigh[iptbin-1] = hslice_data->GetBinError(iptbin);
      if (isTrg) {
		mc_pt_nom[iptbin-1]      = hslice_mc->GetXaxis()->GetBinCenter(iptbin);
		mc_eff_nom[iptbin-1]     = hslice_mc->GetBinContent(iptbin);
		mc_pt_errlow[iptbin-1]   = hslice_mc->GetXaxis()->GetBinCenter(iptbin)-hslice_mc->GetXaxis()->GetBinLowEdge(iptbin);
		mc_pt_errhigh[iptbin-1]  = hslice_mc->GetXaxis()->GetBinUpEdge(iptbin)-hslice_mc->GetXaxis()->GetBinCenter(iptbin);
		mc_eff_errlow[iptbin-1]  = hslice_mc->GetBinError(iptbin);
		mc_eff_errhigh[iptbin-1] = hslice_mc->GetBinError(iptbin);
      }
    }

    eff_data[etaLabel] = new TGraphAsymmErrors(nbin_pt, data_pt_nom, data_eff_nom, data_pt_errlow, data_pt_errhigh, data_eff_errlow, data_eff_errhigh);
	//eff_data[etaLabel]->Write(Form("%s",etaLabel));
	// std::cout << " Before: " << eff_data[etaLabel]->GetXaxis()->GetBinLowEdge(1) << ", " << eff_data[etaLabel]->GetXaxis()->GetBinUpEdge(nbin_pt) << std::endl;
    // ShiftAxisBins(eff_data[etaLabel]);
	// std::cout << " After: " << eff_data[etaLabel]->GetXaxis()->GetBinLowEdge(1) << ", " << eff_data[etaLabel]->GetXaxis()->GetBinUpEdge(nbin_pt) << std::endl;
    if (isTrg) {
      eff_mc[etaLabel] = new TGraphAsymmErrors(nbin_pt, mc_pt_nom, mc_eff_nom, mc_pt_errlow, mc_pt_errhigh, mc_eff_errlow, mc_eff_errhigh);
    //   ShiftAxisBins(eff_mc[etaLabel]);
    }
  }
}

void ScaleFactor::init_ScaleFactor(TString inputRootFile, std::string HistoBaseName) {

  TFile * fileIn = new TFile(inputRootFile, "read");
  // if root file not found
  if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc : File " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };

  if (HistoBaseName == "ZMass") { // efficiency file contains TGraphAsymmErrors + eta map (HTT group format for legacy)
	etaIsAbsolute = true;
    etaBinsH = (TH1D*)fileIn->Get("etaBinsH");
    std::string etaLabel, GraphName;
    int nEtaBins = etaBinsH->GetNbinsX();

    for (int iBin=0; iBin<nEtaBins; iBin++) {
      etaLabel = etaBinsH->GetXaxis()->GetBinLabel(iBin+1);
      GraphName = HistoBaseName+etaLabel+"_Data";

	  if (fileIn->GetListOfKeys()->Contains(TString(GraphName))) {
		// std::cout << " ### init_ScaleFactor" << std::endl;
		eff_data[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
		// std::cout << " Before: " << eff_data[etaLabel]->GetXaxis()->GetBinLowEdge(1) << ", " << eff_data[etaLabel]->GetXaxis()->GetBinUpEdge(eff_data[etaLabel]->GetN()) << std::endl;
		SetAxisBins(eff_data[etaLabel]);
		//eff_data[etaLabel]->Write(etaLabel.c_str());
		// std::cout << " After: " << eff_data[etaLabel]->GetXaxis()->GetBinLowEdge(1) << ", " << eff_data[etaLabel]->GetXaxis()->GetBinUpEdge(eff_data[etaLabel]->GetN()) << std::endl;

	  }
	  else eff_data[etaLabel] = 0;

	  GraphName = HistoBaseName+etaLabel+"_MC";
	  if (fileIn->GetListOfKeys()->Contains(TString(GraphName))) {
		eff_mc[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
		SetAxisBins(eff_mc[etaLabel]);
	  }
	  else eff_mc[etaLabel] =0;

	  if (eff_data[etaLabel] != 0 && eff_mc[etaLabel] != 0 ) {
		bool sameBinning = check_SameBinning(eff_data[etaLabel], eff_mc[etaLabel]);
		if (!sameBinning) {std::cout<< "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc . Can not proceed because ScaleFactor::check_SameBinning returned different pT binning for data and MC for eta label " << etaLabel << std::endl; exit(1); }; 
	  }
    }
  } 
  else { //efficiency maps in TH2F -> create eta map & TGraphAsymmErrors so that rest of the pipeline isn't disturbed

    // somewhat ugly as applying only to muon SFs
    TString effname_data = "_efficiencyData";
    TString effname_MC = "_efficiencyMC";

	etaIsAbsolute = true;
	TH2F *hSF = (TH2F*)fileIn->Get((TString)HistoBaseName);
    TH2F *heffdat = (TH2F*)fileIn->Get((TString)HistoBaseName+effname_data);
    TH2F *heffmc = (TH2F*)fileIn->Get((TString)HistoBaseName+effname_MC);

    // retrieve eta binning (ugly, but should work fine)
    const int nbin_eta = hSF->GetNbinsX();;
    TString eta_bins[nbin_eta] = {""};
    TString firstbinlabel = Form("EtaLt%.1f",hSF->GetXaxis()->GetBinUpEdge(1));
    TString lastbinlabel  = Form("EtaGt%.1f",hSF->GetXaxis()->GetBinLowEdge(nbin_eta));
    firstbinlabel.ReplaceAll(".","p");
    lastbinlabel.ReplaceAll(".","p");

    //create etabinning Histo
    etaBinsH = new TH1D("etaBinsH","",nbin_eta, hSF->GetXaxis()->GetXbins()->GetArray());

    TString TetaLabel;
    TString GraphName;
    for (int iBin=1; iBin<=nbin_eta; iBin++) {
      if (iBin==1) {
		TetaLabel = firstbinlabel;
      }
      else if (iBin==nbin_eta) {
		TetaLabel = lastbinlabel;
      }
      else {
		TetaLabel = Form("Eta%.1fto%.1f",
			hSF->GetXaxis()->GetBinLowEdge(iBin),
			hSF->GetXaxis()->GetBinUpEdge(iBin));
		TetaLabel.ReplaceAll(".","p");
      }
      etaBinsH->GetXaxis()->SetBinLabel(iBin,TetaLabel);

      std::string etaLabel = (std::string)TetaLabel;
      // GraphName = TString(HistoBaseName)+"_"+etaLabel+"_Data";
	  
      TH1F *hslice_data = (TH1F*)heffdat->ProjectionY("slicedata",iBin,iBin);
      TH1F *hslice_mc   = (TH1F*)heffmc->ProjectionY("slicemc",iBin,iBin);

      const int nbin_pt = hslice_data->GetNbinsX();

      double data_pt_nom[nbin_pt] = {0};
      double data_eff_nom[nbin_pt] = {0};
      double data_pt_errlow[nbin_pt] = {0};
      double data_eff_errlow[nbin_pt] = {0};
      double data_pt_errhigh[nbin_pt] = {0};
      double data_eff_errhigh[nbin_pt] = {0};

      double mc_pt_nom[nbin_pt] = {0};
      double mc_eff_nom[nbin_pt] = {0};
      double mc_pt_errlow[nbin_pt] = {0};
      double mc_eff_errlow[nbin_pt] = {0};
      double mc_pt_errhigh[nbin_pt] = {0};
      double mc_eff_errhigh[nbin_pt] = {0};

      for (int iptbin=1; iptbin<=nbin_pt; iptbin++) {
		data_pt_nom[iptbin-1]      = hslice_data->GetXaxis()->GetBinCenter(iptbin);
		data_eff_nom[iptbin-1]     = hslice_data->GetBinContent(iptbin);
		data_pt_errlow[iptbin-1]   = hslice_data->GetXaxis()->GetBinCenter(iptbin) - hslice_data->GetXaxis()->GetBinLowEdge(iptbin);
		data_pt_errhigh[iptbin-1]  = hslice_data->GetXaxis()->GetBinUpEdge(iptbin) - hslice_data->GetXaxis()->GetBinCenter(iptbin);
		data_eff_errlow[iptbin-1]  = hslice_data->GetBinError(iptbin);
		data_eff_errhigh[iptbin-1] = hslice_data->GetBinError(iptbin);

		mc_pt_nom[iptbin-1]      = hslice_mc->GetXaxis()->GetBinCenter(iptbin);
		mc_eff_nom[iptbin-1]     = hslice_mc->GetBinContent(iptbin);
		mc_pt_errlow[iptbin-1]   = hslice_mc->GetXaxis()->GetBinCenter(iptbin)-hslice_mc->GetXaxis()->GetBinLowEdge(iptbin);
		mc_pt_errhigh[iptbin-1]  = hslice_mc->GetXaxis()->GetBinUpEdge(iptbin)-hslice_mc->GetXaxis()->GetBinCenter(iptbin);
		mc_eff_errlow[iptbin-1]  = hslice_mc->GetBinError(iptbin);
		mc_eff_errhigh[iptbin-1] = hslice_mc->GetBinError(iptbin);
      }

      eff_data[etaLabel] = new TGraphAsymmErrors(nbin_pt, data_pt_nom, data_eff_nom, data_pt_errlow, data_pt_errhigh, data_eff_errlow, data_eff_errhigh);
      eff_mc[etaLabel]   = new TGraphAsymmErrors(nbin_pt, mc_pt_nom,   mc_eff_nom,   mc_pt_errlow,   mc_pt_errhigh,   mc_eff_errlow,   mc_eff_errhigh  );

	}
  }
  return;
}


void ScaleFactor::SetAxisBins(TGraphAsymmErrors* graph) {

	int NPOINTS = graph->GetN(); 
	double AXISBINS [NPOINTS+1] = {};
	for (int i=0; i<NPOINTS; i++) { AXISBINS[i] = (graph->GetX()[i] - graph->GetErrorXlow(i)); }
	AXISBINS[NPOINTS] = (graph->GetX()[NPOINTS-1] + graph->GetErrorXhigh(NPOINTS-1));
	graph->GetXaxis()->Set(NPOINTS, AXISBINS);
	return;
}
/* Not exactly sure what this is supposed to be doing. DO not use unless you want to debug it
void ScaleFactor::ShiftAxisBins(TGraphAsymmErrors* graph) {

	std::cout << "BEFORE : ";
	int NPOINTS_1 = graph->GetN(); 
	double AXISBINS_1 [NPOINTS_1+1] = {};
	for (int i=0; i<NPOINTS_1; i++) {
		std::cout << graph->GetX()[i] - graph->GetErrorXlow(i) << " - ";
		AXISBINS_1[i] = (graph->GetX()[i] - graph->GetErrorXlow(i)); 
	}
	AXISBINS_1[NPOINTS_1] = (graph->GetX()[NPOINTS_1-1] + graph->GetErrorXhigh(NPOINTS_1-1));
	std::cout << (graph->GetX()[NPOINTS_1-1] + graph->GetErrorXhigh(NPOINTS_1-1)) << std::endl;
	
	std::cout << "AFTER : 0";
	int NPOINTS = graph->GetN(); 
	double AXISBINS [NPOINTS+2] = {};
	AXISBINS[0] = 0.;
	for (int i = 0; i < NPOINTS; i++) {
		std::cout << " - " << graph->GetErrorXlow(i);
		AXISBINS[i+1] = graph->GetErrorXlow(i);
	}
	std::cout << " - " << graph->GetErrorXhigh(NPOINTS-1);;
	AXISBINS[NPOINTS+1] = graph->GetErrorXhigh(NPOINTS-1);
	std::cout << std::endl;
	graph->GetXaxis()->Set(NPOINTS, AXISBINS);
	return;
}
*/

bool ScaleFactor::check_SameBinning(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2){
	bool haveSameBins = false;
	int n1 = graph1->GetXaxis()->GetNbins();
	int n2 = graph2->GetXaxis()->GetNbins();
	if (n1 != n2 ) {return false;}
	else {
		haveSameBins = true;
		const int nbins = n1;
		double x1, x2;
		for (int i=0; i<nbins; i++){ 
			x1 = (graph1->GetXaxis()->GetXbins())->GetArray()[i];
			x2 = (graph2->GetXaxis()->GetXbins())->GetArray()[i]; 
			haveSameBins = haveSameBins and (x1== x2) ;
		}
	}

	return haveSameBins;
}


std::string ScaleFactor::FindEtaLabel(double Eta, std::string Which){
	if (etaIsAbsolute)
		Eta = fabs(Eta);
	int binNumber = etaBinsH->GetXaxis()->FindFixBin(Eta);
	std::string EtaLabel = etaBinsH->GetXaxis()->GetBinLabel(binNumber);
	std::map<std::string, TGraphAsymmErrors*>::iterator it;

	if (Which == "data"){
		it =  eff_data.find(EtaLabel);
		if ( it == eff_data.end()) { 
			std::string error = "ERROR in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : no object corresponding to eta label "
				+ EtaLabel + " for data ";
			std::cerr << error << std::endl;
			throw std::runtime_error(error); // sometimes throwing exception leads to segmentation fault, so we print the error first to stderr
		}
	}

	else if (Which == "mc"){
		it = eff_mc.find(EtaLabel);
		if (it == eff_mc.end()) { 
			std::string error = "ERROR in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : no object corresponding to eta label "
				+ EtaLabel + " for MC ";
			std::cerr << error << std::endl;
			throw std::runtime_error(error);
		}		
	}
	
     return EtaLabel;
}

/**
 * Computes the bin number for given pt. Output is in range [1; Npoints] (so shifted by 1 compared to points in TGraphAsymErrors)
*/
int ScaleFactor::FindPtBin( std::map<std::string, TGraphAsymmErrors *> eff_map, std::string EtaLabel, double Pt){
    int Npoints = eff_map[EtaLabel]->GetN();
	double ptMAX = (eff_map[EtaLabel]->GetX()[Npoints-1])+(eff_map[EtaLabel]->GetErrorXhigh(Npoints-1));
	double ptMIN = (eff_map[EtaLabel]->GetX()[0])-(eff_map[EtaLabel]->GetErrorXlow(0));
	// if pt is overflow, return last pt bin
 	if (Pt >= ptMAX ) return Npoints; 
	// if pt is underflow, return nonsense number and warning
	else if (Pt < ptMIN){ 
		std::cout<< "WARNING in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: pT too low (pt = " << Pt << "), min value is " << ptMIN << ". Returned efficiency =1. Weight will be 1. " << std::endl;
		return -99;
	}
	// if pt is in range
	else {
		// this does not work (it is a TGraph, not a histogram. )
		// return eff_map[EtaLabel]->GetXaxis()->FindFixBin(Pt);

		for (int graphBin=0; graphBin < Npoints; graphBin++) {
			if (Pt >= eff_map[EtaLabel]->GetPointX(graphBin) - eff_map[EtaLabel]->GetErrorXlow(graphBin)
			  	&& Pt < eff_map[EtaLabel]->GetPointX(graphBin) + eff_map[EtaLabel]->GetErrorXhigh(graphBin))
				return graphBin+1;
		}
		assert(false && "ScaleFactor::FindPtBin : TGraphAsymErrors has gaps : this should not happen");
	} 
}




double ScaleFactor::get_EfficiencyData(double pt, double eta){

    double eff;
	std::string label = FindEtaLabel(eta, "data");

	// std::cout << "  eta " << eta << ";  label  " << label << std::endl;

	int ptbin = FindPtBin(eff_data, label, pt); 

	// std::cout << "  pt " << pt << ";  ptbin  " << ptbin << "; old " << eff_data[label]->GetXaxis()->FindFixBin(pt) << std::endl;

	if (ptbin == -99) {
		eff = 1; // if pt is underflow 
	}
	else {
		eff = eff_data[label]->GetY()[ptbin-1];
	}

	if (eff > 1.) {std::cout<< "WARNING in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: Efficiency in data > 1. Set eff = 1." << std::endl; eff=1;} 
	if (eff < 0 ) {std::cout<< "WARNING in ScaleFactor::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: Negative efficiency in data. Set eff = 0." <<std::endl; eff=0;}

	return eff;
	
}


double ScaleFactor::get_EfficiencyMC(double pt, double eta) {

	double eff;		
	std::string label = FindEtaLabel(eta,"mc");

	int ptbin = FindPtBin(eff_mc, label, pt); 
	if (ptbin == -99){eff =1;} // if pt is underflow 
	else eff= eff_mc[label]->GetY()[ptbin-1];

	if (eff > 1.) {std::cout << "WARNING in ScaleFactor::get_EfficiencyMC(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : Efficiency in MC > 1. Set eff = 1." << std::endl; eff =1;} 		
	if (eff < 0 ) {std::cout << "WARNING in ScaleFactor::get_EfficiencyMC(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : Negative efficiency in MC. Set eff = 0." <<std::endl; eff =0;}

	return eff;

}



double ScaleFactor::get_ScaleFactor(double pt, double eta){
	
	double efficiency_data = get_EfficiencyData(pt, eta);
	double efficiency_mc = get_EfficiencyMC(pt, eta);
	double SF;

	if ( efficiency_mc != 0) {SF = efficiency_data/efficiency_mc;}
	else {
	SF=1.; std::cout << "WARNING in ScaleFactor::get_ScaleFactor(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : MC efficiency = 0. Scale Factor set to 1. ";
	}

	// std::cout << " ### DEBUG Pt = " << pt << " , Eta = " << eta << "		 SF = " << SF << std::endl;
	// std::cout << " eff data = " << efficiency_data << " ; eff mc = " << efficiency_mc << std::endl;

	return SF;	
	
}


double ScaleFactor::get_EfficiencyDataError(double pt, double eta) {

	double eff_error;
	std::string label = FindEtaLabel(eta,"data");
	int ptbin = FindPtBin(eff_data, label, pt); 
	if (ptbin == -99){eff_error =0.;} // if pt is underflow 
	else eff_error= eff_data[label]->GetErrorYhigh(ptbin-1); 
        // errors are supposed to be symmetric, can use GetErrorYhigh or GetErrorYlow

	double effData = get_EfficiencyData(pt,eta);
	if (eff_error > effData) eff_error = 0.5*effData;
	return eff_error;
}
	
	

double ScaleFactor::get_EfficiencyMCError(double pt, double eta){

	double eff_error;
	std::string label = FindEtaLabel(eta,"mc");
	int ptbin = FindPtBin(eff_mc, label, pt); 
	if (ptbin == -99){eff_error =0.;} // if pt is underflow 
	else eff_error= eff_mc[label]->GetErrorYhigh(ptbin-1); 
	// errors are supposed to be symmetric, can use GetErrorYhigh or GetErrorYlow

	double effMC = get_EfficiencyMC(pt,eta);
	if (eff_error > effMC ) eff_error = 0.5*effMC;
	return eff_error;
}

double ScaleFactor::get_ScaleFactorError(double pt, double eta){

	double SF_error = 0.;
	
	double effData = get_EfficiencyData(pt, eta);
	double effMC = get_EfficiencyMC(pt, eta);
	double errData = get_EfficiencyDataError(pt, eta);
	double errMC =  get_EfficiencyMCError(pt, eta);

	if (errData == 0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: uncertainty on data point = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (errMC ==0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: uncertainty on MC = 0, can not calculate uncerttainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (effData ==0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: efficiency in data = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (effMC ==0) {std::cout<<"WARNING in ScaleFactor::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: efficiency in MC = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	else {	
	SF_error = pow((errData/effData),2) + pow((errMC/effMC),2);
	SF_error = pow(SF_error, 0.5)*(effData/effMC);
	}
	return SF_error;
}
	

