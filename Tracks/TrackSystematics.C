#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
#include "iostream"
// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
// #include "./TrackMcQC_settings.h"
#include "./TrackSystematics_inputs.h"
#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramPlotting.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well


#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();


template <std::size_t N> 
std::vector<double> Get_SelTrackYield_Pt(TH1D* (&hYieldReturned)[N], bool isMC);
template <std::size_t N>
bool Get_systematics_trackSelectionVariation_pt_betterVersion(TH1D* (&hRelativeEfficiency)[N], TH1D* (&hRelativeEfficiency_PreBarlow)[N], bool isMC, std::string options);
void Draw_Systematics_trackSelectionVariation_pt_betterVersion(std::string options);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void TrackSystematics
() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();
  // TString* SaveAs_Title = new TString("");

  bool useSplit = false; //set to true if you want to see the influence of split
  float etaRange[2] = {-0.9, 0.9};
  float ptRange[2] = {0.15, 100};

  Draw_Systematics_trackSelectionVariation_pt_betterVersion("");
}
/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs() {
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  // gSystem->Load("libMinuit");
  // gSystem->Load("libSTEERBase");
  // gSystem->Load("libESD");
  // gSystem->Load("libAOD");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libCORRFW");
  // gSystem->Load("libPWGTools");
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLineScalePS(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}




// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextTrackDatasetComp(std::string options){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  //  if (options.find("track") == std::string::npos) {
  texcontextDatasetCompAndRadiusAndVarRange = *texDatasetsComparisonCommonDenominator;
  return texcontextDatasetCompAndRadiusAndVarRange;
}


// TString contextPtRange(float* PtRange){
//   std::stringstream ss;
//   ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
//   TString textContext((TString)ss.str());
//   // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
//   return textContext;
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <std::size_t N> 
std::vector<double> Get_SelTrackYield_Pt(TH1D* (&hYieldReturned)[N], bool isMC) {

  TFile* inputTFileArray[nDatasets];
  TString analysisWorkflowsArray[nDatasets];
  TString Datasets[nDatasets];

  TH2D* H1D_trackPtTrackSigma[nDatasets];
  TH2D* H1D_trackPtTrackSigma_ptHigh[nDatasets];

  TH1D* H1D_trackPt[nDatasets];
  TH1D* H1D_trackPt_ptHigh[nDatasets];

  TH1D* H1D_trackPt_rebinned[nDatasets];
  TH1D* H1D_trackPt_ptHigh_rebinned[nDatasets];

  TH1D* H1D_trackPt_yield_concatenated[nDatasets];

  int nBinsPt;
  std::vector<double> ptBinning;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    inputTFileArray[iDataset] = isMC ?  file_O2Analysis_list_MC[iDataset] : file_O2Analysis_list_Data[iDataset];
    analysisWorkflowsArray[iDataset] = isMC ? analysisWorkflowsMC[iDataset] : analysisWorkflowsData[iDataset];
    Datasets[iDataset] = isMC ? DatasetsMC[iDataset] : DatasetsData[iDataset];

    cout << "dataset " << iDataset << endl;
    H1D_trackPtTrackSigma[iDataset] = (TH2D*)((TH2D*)inputTFileArray[iDataset]->Get(analysisWorkflowsArray[iDataset]+"/h2_track_pt_track_sigmapt"))->Clone("Get_SelTrackYield_Pt"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H1D_trackPtTrackSigma_ptHigh[iDataset] = (TH2D*)((TH2D*)inputTFileArray[iDataset]->Get(analysisWorkflowsArray[iDataset]+"/h2_track_pt_high_track_sigmapt"))->Clone("Get_SelTrackYield_Pt"+Datasets[iDataset]+DatasetsNames[iDataset]);

    H1D_trackPt[iDataset] = (TH1D*)H1D_trackPtTrackSigma[iDataset]->ProjectionX("trackPt_tracks"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    H1D_trackPt_ptHigh[iDataset] = (TH1D*)H1D_trackPtTrackSigma_ptHigh[iDataset]->ProjectionX("trackPt_tracks_ptHigh"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");


    //tweaking the low-pt bins to make them larger close to 10GeV
    std::vector<double> xbinsVectorInitialLow = GetTH1Bins(H1D_trackPt[iDataset]);
    double* xbinsInitialLow = &xbinsVectorInitialLow[0];

    double ptBinsLowNew[500]; //500 to have a good margin
    double lastPt;
    int iBinPtNew = 0;
    int iBinPtInitialHisto = 0;
    int increment;
    float originalBinningEnd = 0.3; //above that we switch to custom binning
    while(iBinPtInitialHisto < H1D_trackPt[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialLow[iBinPtInitialHisto];
      if (lastPt < originalBinningEnd) { // top boundary for the unchanged binning
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt[iDataset]->GetXaxis()->GetBinWidth(1) /8;
      }
      ptBinsLowNew[iBinPtNew] = lastPt;
      iBinPtInitialHisto += increment;
      iBinPtNew += 1;
    }
    ptBinsLowNew[iBinPtNew-1] = 10; // replace last set bin edge by 10
    int nBinsLowNew = iBinPtNew - 1;
    if (nBinsLowNew >= 500){
      cout << "NEEDS TO RESERVE MORE MEMORY FOR ptBinsLowNew !!!!!!!!!!!!!!!!!!" << endl;
    }
    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsLowNew);


    //tweaking the high-pt bins to have them less prone to statistical fluctuations
    std::vector<double> xbinsVectorInitialHigh = GetTH1Bins(H1D_trackPt_ptHigh[iDataset]);
    double* xbinsInitialHigh = &xbinsVectorInitialHigh[0];

    double ptBinsHighNew[500]; //500 to have a good margin
    iBinPtNew = 0;
    iBinPtInitialHisto = 0;
    while(iBinPtInitialHisto < H1D_trackPt_ptHigh[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialHigh[iBinPtInitialHisto];
      if (lastPt < 10) { // was 25 before
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_ptHigh[iDataset]->GetXaxis()->GetBinWidth(1) ;
      }
      ptBinsHighNew[iBinPtNew] = lastPt;
      iBinPtInitialHisto += increment;
      iBinPtNew += 1;
    }
    ptBinsHighNew[iBinPtNew-1] = 100; // replace last set bin edge by 10
    int nBinsHighNew = iBinPtNew - 1;
    if (nBinsHighNew >= 500){
      cout << "NEEDS TO RESERVE MORE MEMORY FOR ptBinsHighNew !!!!!!!!!!!!!!!!!!" << endl;
    }
    H1D_trackPt_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_ptHigh_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset], ptBinsHighNew);


    // Merging high and low pt histograms:
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins(H1D_trackPt_rebinned[iDataset]);
    std::vector<double> xbinsVectorRight = GetTH1Bins(H1D_trackPt_ptHigh_rebinned[iDataset]);
    xbinsVectorRight.erase(xbinsVectorRight.begin());
    // cout << "xbinsVectorLeft.front() = " << xbinsVectorLeft.front() << ", xbinsVectorLeft.back() = " << xbinsVectorLeft.back() << ", xbinsVectorRight.front() = " << xbinsVectorRight.front() << ", xbinsVectorRight.back() = " << xbinsVectorRight.back() << endl;
    std::vector<double> xbinsVectorCombination = xbinsVectorLeft;
    xbinsVectorCombination.insert( xbinsVectorCombination.end(), xbinsVectorRight.begin(), xbinsVectorRight.end() );
    double* xbins_new = &xbinsVectorCombination[0];
    nBinsPt = xbinsVectorCombination.size()-1;

    if (iDataset == 0) {
      ptBinning = xbinsVectorCombination;
    }
    // cout << "xbinsVectorCombination.size() = " << xbinsVectorCombination.size() << endl;

    // making the new hist with concatenated bins
    TH1D H1D_trackPt_yield_concatenated_temp("H1D_trackPt_yield_concatenated_temp"+Datasets[iDataset]+DatasetsNames[iDataset], "H1D_trackPt_yield_concatenated_temp"+Datasets[iDataset]+DatasetsNames[iDataset], nBinsPt, xbins_new);
    H1D_trackPt_yield_concatenated[iDataset] = (TH1D*)H1D_trackPt_yield_concatenated_temp.Clone("H1D_trackPt_yield_concatenated"+Datasets[iDataset]+DatasetsNames[iDataset]);

    //filling the two new histograms
    for(int iBinX = 1; iBinX <= H1D_trackPt_rebinned[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_yield_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_yield_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_rebinned[iDataset]->GetBinError(iBinX));
    }
    for(int iBinX = 1; iBinX <= H1D_trackPt_ptHigh_rebinned[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_yield_concatenated[iDataset]->SetBinContent(H1D_trackPt_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_ptHigh_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_yield_concatenated[iDataset]->SetBinError(H1D_trackPt_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_ptHigh_rebinned[iDataset]->GetBinError(iBinX));
    }
    
    hYieldReturned[iDataset] = (TH1D*)((TH1D*)H1D_trackPt_yield_concatenated[iDataset]->Clone("hYieldReturned"+Datasets[iDataset]+DatasetsNames[iDataset]));
  }

  return ptBinning;
}



template <std::size_t N>
bool Get_systematics_trackSelectionVariation_pt_betterVersion(TH1D* (&hRelativeEfficiency)[N], TH1D* (&hRelativeEfficiency_PreBarlow)[N], bool isMC, std::string options) {
  // assumes nominal in 1st (0th) position
  // returns errors in MC or Data; not relative errors!
  int iNominal = 0;

  TString Datasets[nDatasets];
  
  // return histogram that has the systematics in its contents

  TH1D* H1D_trackPtYield[nDatasets];

  std::vector<double> ptBinning = Get_SelTrackYield_Pt(H1D_trackPtYield, isMC);

  // if (ptBinningMC.size() != ptBinningData.size()) {
  //   cout << "ptBinningMC and ptBinningData have different sizes, they both should be identical" << endl;
  // }
  int nBinsPt = ptBinning.size();
  double* ptBinningArray = &ptBinning[0];

  TH1D* H1D_trackPtYield_difference[nDatasets-1];
  for(int iDataset = 1; iDataset < nDatasets; iDataset++){
    Datasets[iDataset] = isMC ? DatasetsMC[iDataset] : DatasetsData[iDataset];

    H1D_trackPtYield_difference[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPtYield[iDataset]->Clone("H1D_trackPtYield_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]+(TString)(isMC ? "MC" : "Data")));
    H1D_trackPtYield_difference[iDataset-1]->Add(H1D_trackPtYield[iDataset], H1D_trackPtYield[iNominal], -1, 1); // 1*Nominal + -1*Tight
  }

  // cout << "Do I apply Barlow condition even though not param variation ? check paper again" << endl; YES, the subset of data thing is only shown for first demonstration, but barlow says it holds true even if that's not the cast

  TH1D* hTempYieldDifferenceToNominal[nDatasets]; // = new TH1D("hTempYieldDifferenceToNominal", "hTempYieldDifferenceToNominal", nBinsPt, ptBinningArray);
  TH1D* hTempYieldDifferenceToNominal_PreBarlow[nDatasets]; // = new TH1D("hTempYieldDifferenceToNominal_PreBarlow", "hTempYieldDifferenceToNominal_PreBarlow", nBinsPt, ptBinningArray);


  /////////////////
  // Barlow test //
  /////////////////

  TH1D* H1D_trackPtYield_nominal = H1D_trackPtYield[iNominal];
  double differenceToNominal;
  int id_SignalExtractionType_maxDeviation;
  double hSigmaBarlow[nBinsPt];
  TH1D* hYieldDifferenceToNominal[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hYieldDifferenceToNominal_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  for(int iDataset = 1; iDataset < nDatasets; iDataset++){ 
    hTempYieldDifferenceToNominal[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPtYield_difference[iDataset-1]->Clone("hTempYieldDifferenceToNominalPt"+Datasets[iDataset]+DatasetsNames[iDataset]+(TString)(isMC ? "MC" : "Data")));
    hTempYieldDifferenceToNominal_PreBarlow[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPtYield_difference[iDataset-1]->Clone("hTempYieldDifferenceToNominalPt_PreBarlow"+Datasets[iDataset]+DatasetsNames[iDataset]+(TString)(isMC ? "MC" : "Data")));
    hTempYieldDifferenceToNominal[iDataset-1]->Sumw2();
    hTempYieldDifferenceToNominal_PreBarlow[iDataset-1]->Sumw2();
    hTempYieldDifferenceToNominal[iDataset-1]->Reset("M");
    hTempYieldDifferenceToNominal_PreBarlow[iDataset-1]->Reset("M");

    for(int iBinPt = 1; iBinPt <= nBinsPt; iBinPt++){
      differenceToNominal = H1D_trackPtYield_difference[iDataset-1]->GetBinContent(iBinPt);

      // Barlow condition for systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
      Double_t StatUncertainty_REF = H1D_trackPtYield_nominal->GetBinError(iBinPt);
      Double_t StatUncertainty_MaxDeviationCase = H1D_trackPtYield[iDataset]->GetBinError(iBinPt);
  
      int PtArrayIterator = iBinPt - 1;
      hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
  
      hTempYieldDifferenceToNominal_PreBarlow[iDataset-1]->SetBinContent(iBinPt,differenceToNominal);
      hTempYieldDifferenceToNominal_PreBarlow[iDataset-1]->SetBinError(iBinPt,hSigmaBarlow[PtArrayIterator]);

      if (differenceToNominal > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
        hTempYieldDifferenceToNominal[iDataset-1]->SetBinContent(iBinPt,differenceToNominal);
        hTempYieldDifferenceToNominal[iDataset-1]->SetBinError(iBinPt,hSigmaBarlow[PtArrayIterator]);
      }
      else {
        hTempYieldDifferenceToNominal[iDataset-1]->SetBinContent(iBinPt,0.);
      }
    }

    // hYieldDifferenceToNominal = yield(nominal) - yield(tightCut)
    hYieldDifferenceToNominal[iDataset-1] = (TH1D*)((TH1D*)hTempYieldDifferenceToNominal[iDataset-1]->Clone("hYieldDifferenceToNominal_PtYield_"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
    hYieldDifferenceToNominal_PreBarlow[iDataset-1] = (TH1D*)((TH1D*)hTempYieldDifferenceToNominal_PreBarlow[iDataset-1]->Clone("hYieldDifferenceToNominal_PreBarlow_PtYield_"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
  }

  bool divideSuccess = true;
  bool divideSuccessPreBarlow = true;
  bool divideSuccessTemp;
  for(int iDataset = 1; iDataset < nDatasets; iDataset++){ 
    // hRelativeEfficiency = yield(tightCut)/yield(nominal) = (yield(nominal) - hYieldDifferenceToNominal)/yield(nominal)
    hRelativeEfficiency[iDataset-1] = (TH1D*)((TH1D*)hYieldDifferenceToNominal[iDataset-1]->Clone("hRelativeEfficiency"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
    hRelativeEfficiency[iDataset-1]->Scale(-1);
    hRelativeEfficiency[iDataset-1]->Add(H1D_trackPtYield_nominal);
    divideSuccessTemp = hRelativeEfficiency[iDataset-1]->Divide(H1D_trackPtYield_nominal);
    divideSuccess = divideSuccess && divideSuccessTemp;
    hRelativeEfficiency_PreBarlow[iDataset-1] = (TH1D*)((TH1D*)hYieldDifferenceToNominal_PreBarlow[iDataset-1]->Clone("hRelativeEfficiency_PreBarlow"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
    hRelativeEfficiency_PreBarlow[iDataset-1]->Scale(-1);
    hRelativeEfficiency_PreBarlow[iDataset-1]->Add(H1D_trackPtYield_nominal);
    divideSuccessTemp = hRelativeEfficiency_PreBarlow[iDataset-1]->Divide(H1D_trackPtYield_nominal);
    divideSuccessPreBarlow = divideSuccessPreBarlow && divideSuccessTemp;
  }
  return (divideSuccess && divideSuccessPreBarlow);
}

void Draw_Systematics_trackSelectionVariation_pt_betterVersion(std::string options) {

  TH1D* hRelativeEfficiency_Data[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hRelativeEfficiency_Data_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  TH1D* hRelativeEfficiency_MC[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hRelativeEfficiency_MC_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  bool isMC_Data = false;
  // Get_systematics_trackSelectionVariation_pt_betterVersion(hRelativeEfficiency_Data, hRelativeEfficiency_Data_PreBarlow, isMC_Data, options);
  bool divideSuccessData = Get_systematics_trackSelectionVariation_pt_betterVersion(hRelativeEfficiency_Data, hRelativeEfficiency_Data_PreBarlow, isMC_Data, options);
  if (!divideSuccessData) {
    cout << "divide failed for MC Get_systematics_trackSelectionVariation_pt_betterVersion() in Draw_Systematics_trackSelectionVariation_pt_betterVersion() Data" << endl;
  }

  bool isMC_MC = true;
  Get_systematics_trackSelectionVariation_pt_betterVersion(hRelativeEfficiency_MC, hRelativeEfficiency_MC_PreBarlow, isMC_MC, options);
  bool divideSuccessMC = Get_systematics_trackSelectionVariation_pt_betterVersion(hRelativeEfficiency_MC, hRelativeEfficiency_MC_PreBarlow, isMC_MC, options);
  if (!divideSuccessMC) {
    cout << "divide failed for Data Get_systematics_trackSelectionVariation_pt_betterVersion() in Draw_Systematics_trackSelectionVariation_pt_betterVersion() MC" << endl;
  }

  TString systematicsLegend[nDatasets-1];
  std::stringstream ss;

  TH1D* hSystematicUncertainty[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hSystematicUncertainty_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  TH1D* H1D_trackPt_efficiency[nDatasets];
  for(int iDataset = 1; iDataset < nDatasets; iDataset++){ 
    cout << "dataset " << iDataset << endl;
    hSystematicUncertainty[iDataset-1] = (TH1D*)hRelativeEfficiency_Data[iDataset-1]->Clone("hSystematicUncertainty"+DatasetsData[iDataset]+DatasetsNames[iDataset]);
    hSystematicUncertainty[iDataset-1]->Divide(hRelativeEfficiency_MC[iDataset-1]);
    hSystematicUncertainty_PreBarlow[iDataset-1] = (TH1D*)hRelativeEfficiency_Data_PreBarlow[iDataset-1]->Clone("hSystematicUncertainty_PreBarlow"+DatasetsData[iDataset]+DatasetsNames[iDataset]);
    hSystematicUncertainty_PreBarlow[iDataset-1]->Divide(hRelativeEfficiency_MC_PreBarlow[iDataset-1]);
   
    ss << DatasetsNames[iDataset];
    systematicsLegend[iDataset-1] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }


  TString* pdfName_RelativeEfficiency_Data = new TString("RelativeEfficiency_Data_pt");
  TString* pdfName_RelativeEfficiency_MC = new TString("RelativeEfficiency_MC_pt");

  TString* pdfName_systUncertainty = new TString("Systematics_Efficiency_pt");
  TString* pdfName_systUncertainty_PreBarlow = new TString("Systematics_Efficiency_pt_PreBarlow");

  TString* pdfName_systUncertainty_zoom = new TString("Systematics_Efficiency_pt_zoom");
  TString* pdfName_systUncertainty_PreBarlow_zoom = new TString("Systematics_Efficiency_pt_PreBarlow_zoom");

  TString textContext(contextCustomOneField(*texDatasetsComparisonType, ""));


  // std::array<std::array<float, 2>, 2> drawnWindowLog = {{{(float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(1)
  //                                                         +hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(2))/2, 
  //                                                         (float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(hSystematicUncertainty[0]->GetNbinsX())
  //                                                         +hSystematicUncertainty[0]->GetXaxis()->GetBinWidth(hSystematicUncertainty[0]->GetNbinsX()))}
  //                                                         , {0, 0.025}}};


  std::array<std::array<float, 2>, 2> drawnWindow_zoom = {{{(float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(hSystematicUncertainty[0]->GetNbinsX())
                                                          +hSystematicUncertainty[0]->GetXaxis()->GetBinWidth(hSystematicUncertainty[0]->GetNbinsX()))}
                                                          , {0.97, 1.06}}};

  Draw_TH1_Histograms(hRelativeEfficiency_Data, systematicsLegend, nDatasets-1, textContext, pdfName_RelativeEfficiency_Data, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hRelativeEfficiency_MC, systematicsLegend, nDatasets-1, textContext, pdfName_RelativeEfficiency_MC, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");


  Draw_TH1_Histograms(hSystematicUncertainty, systematicsLegend, nDatasets-1, textContext, pdfName_systUncertainty, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hSystematicUncertainty_PreBarlow, systematicsLegend, nDatasets-1, textContext, pdfName_systUncertainty_PreBarlow, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hSystematicUncertainty, systematicsLegend, nDatasets-1, textContext, pdfName_systUncertainty_zoom, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindow_zoom, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hSystematicUncertainty_PreBarlow, systematicsLegend, nDatasets-1, textContext, pdfName_systUncertainty_PreBarlow_zoom, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindow_zoom, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
}



// need option b on division probably