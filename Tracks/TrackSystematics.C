#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
#include "iostream"
#include "TFractionFitter.h"
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
#include <locale>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();


template <std::size_t N> 
std::vector<double> Get_SelTrackYield_Pt(TH1D* (&hYieldReturned)[N], bool isMC);
template <std::size_t N>
bool Get_systematics_trackSelectionVariation_pt(TH1D* (&hRelativeEfficiency)[N], TH1D* (&hRelativeEfficiency_PreBarlow)[N], bool isMC, std::string options);
void Draw_Systematics_trackSelectionVariation_pt(std::string options);


template <std::size_t N> 
void Get_PrimarySecondaryFractions_in_Data(double (&primaryFractionsInData)[N], double (&secondaryFractionsInData)[N], double* ptRange);
template <std::size_t N> 
void Get_DCAxy_yields(TH1D* (&hYieldReturned)[N], TFile* (&fileInput)[N], const TString (&analysisWorkflows)[N], std::string particleStatus, double* ptRange);
template <std::size_t N> 
void Get_ITSTPC_matching_efficiency_Data(TH1D* (&hMatchingEfficiencyReturned_inclusive)[N]);
template <std::size_t N> 
void Get_ITSTPC_matching_efficiency_MC(TH1D* (&hMatchingEfficiencyReturned_mcprimary)[N], TH1D* (&hMatchingEfficiencyReturned_mcsecondary)[N]);
template <std::size_t N> 
void Get_TrackYield_Pt_ITSorTPCorITSTPC_Data(TH1D* (&hYieldReturned_inclusive)[N], int detectorId);
template <std::size_t N> 
void Get_TrackYield_Pt_ITSorTPCorITSTPC_MC(TH1D* (&hYieldReturned_mcprimary)[N], TH1D* (&hYieldReturned_mcsecondary)[N], int detectorId);

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
  double ptRange[2] = {0.3, 100};

  // Draw_Systematics_trackSelectionVariation_pt("");
  double primaryFractionsInData[nDatasets];
  double secondaryFractionsInData[nDatasets];
  Get_PrimarySecondaryFractions_in_Data(primaryFractionsInData, secondaryFractionsInData, ptRange);
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Tracking Efficiency Systematics - quality cut variation ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
bool Get_systematics_trackSelectionVariation_pt(TH1D* (&hRelativeEfficiency)[N], TH1D* (&hRelativeEfficiency_PreBarlow)[N], bool isMC, std::string options) {
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

  // if (ptBinningMC.size() != ptBinningData.size()) {
  //   cout << "ptBinningMC and ptBinningData have different sizes, they both should be identical" << endl;
  // }

  TH1D* H1D_trackPtYield_difference[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  for(int iDataset = 1; iDataset < nDatasets; iDataset++){
    Datasets[iDataset] = isMC ? DatasetsMC[iDataset] : DatasetsData[iDataset];

    H1D_trackPtYield_difference[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPtYield[iDataset]->Clone("H1D_trackPtYield_ratios"+Datasets[iDataset]+DatasetsNames[iDataset]+(TString)(isMC ? "MC" : "Data")));
    H1D_trackPtYield_difference[iDataset-1]->Add(H1D_trackPtYield[iDataset], H1D_trackPtYield[iNominal], -1, 1); // 1*Nominal + -1*Tight
  }

  // cout << "Do I apply Barlow condition even though not param variation ? check paper again" << endl; YES, the subset of data thing is only shown for first demonstration, but barlow says it holds true even if that's not the cast

  TH1D* hTempYieldDifferenceToNominal[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hTempYieldDifferenceToNominal_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1


  /////////////////
  // Barlow test //
  /////////////////

  TH1D* H1D_trackPtYield_nominal = H1D_trackPtYield[iNominal];
  double differenceToNominal;
  int id_SignalExtractionType_maxDeviation;
  double hSigmaBarlow[nBinsPt];
  TH1D* hYieldDifferenceToNominal[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hYieldDifferenceToNominal_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  bool barlowTestSuccess[nDatasets][nBinsPt]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
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
        barlowTestSuccess[iDataset-1][iBinPt] = true;
      }
      else {
        hTempYieldDifferenceToNominal[iDataset-1]->SetBinContent(iBinPt,0.);
        barlowTestSuccess[iDataset-1][iBinPt]  = false;
      }
    }

    // hYieldDifferenceToNominal = yield(nominal) - yield(tightCut)
    hYieldDifferenceToNominal[iDataset-1] = (TH1D*)((TH1D*)hTempYieldDifferenceToNominal[iDataset-1]->Clone("hYieldDifferenceToNominal_PtYield_"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
    hYieldDifferenceToNominal_PreBarlow[iDataset-1] = (TH1D*)((TH1D*)hTempYieldDifferenceToNominal_PreBarlow[iDataset-1]->Clone("hYieldDifferenceToNominal_PreBarlow_PtYield_"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
  }

  // relative efficiency calculation
  bool divideSuccess = true;
  bool divideSuccessPreBarlow = true;
  bool divideSuccessTemp;
  for(int iDataset = 1; iDataset < nDatasets; iDataset++){ 
    // hRelativeEfficiency = yield(tightCut)/yield(nominal) = (yield(nominal) - hYieldDifferenceToNominal)/yield(nominal)
    hRelativeEfficiency[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPtYield[iDataset]->Clone("hRelativeEfficiency"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
    hRelativeEfficiency[iDataset-1]->Reset("M");
    divideSuccessTemp = hRelativeEfficiency[iDataset-1]->Divide(H1D_trackPtYield[iDataset], H1D_trackPtYield_nominal, 1, 1, "b"); // option b because numerator is subset of denominator
    divideSuccess = divideSuccess && divideSuccessTemp;
    for(int iBinPt = 1; iBinPt <= nBinsPt; iBinPt++){
      if (barlowTestSuccess[iDataset-1][iBinPt] == false) {
        // then difference between yield(tightCut) and yield(nominal) is not significant -> relative efficiency is set to 1 with error 0
        hRelativeEfficiency[iDataset-1]->SetBinContent(iBinPt, 1.);
        hRelativeEfficiency[iDataset-1]->SetBinError(iBinPt, 0.);
      }
    }

    hRelativeEfficiency_PreBarlow[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPtYield[iDataset]->Clone("hRelativeEfficiency_PreBarlow"+(TString)(isMC ? "MC" : "Data")+Datasets[iDataset]+DatasetsNames[iDataset]));
    hRelativeEfficiency_PreBarlow[iDataset-1]->Reset("M");
    divideSuccessTemp = hRelativeEfficiency_PreBarlow[iDataset-1]->Divide(H1D_trackPtYield[iDataset], H1D_trackPtYield_nominal, 1, 1, "b"); // option b because numerator is subset of denominator
    divideSuccessPreBarlow = divideSuccessPreBarlow && divideSuccessTemp;
  }
  return (divideSuccess && divideSuccessPreBarlow);
}

void Draw_Systematics_trackSelectionVariation_pt(std::string options) {

  TH1D* hRelativeEfficiency_Data[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hRelativeEfficiency_Data_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  TH1D* hRelativeEfficiency_MC[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
  TH1D* hRelativeEfficiency_MC_PreBarlow[nDatasets]; //only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1

  bool isMC_Data = false;
  // Get_systematics_trackSelectionVariation_pt(hRelativeEfficiency_Data, hRelativeEfficiency_Data_PreBarlow, isMC_Data, options);
  bool divideSuccessData = Get_systematics_trackSelectionVariation_pt(hRelativeEfficiency_Data, hRelativeEfficiency_Data_PreBarlow, isMC_Data, options);
  if (!divideSuccessData) {
    cout << "divide failed for MC Get_systematics_trackSelectionVariation_pt() in Draw_Systematics_trackSelectionVariation_pt() Data" << endl;
  }

  bool isMC_MC = true;
  Get_systematics_trackSelectionVariation_pt(hRelativeEfficiency_MC, hRelativeEfficiency_MC_PreBarlow, isMC_MC, options);
  bool divideSuccessMC = Get_systematics_trackSelectionVariation_pt(hRelativeEfficiency_MC, hRelativeEfficiency_MC_PreBarlow, isMC_MC, options);
  if (!divideSuccessMC) {
    cout << "divide failed for Data Get_systematics_trackSelectionVariation_pt() in Draw_Systematics_trackSelectionVariation_pt() MC" << endl;
  }

  TString systematicsLegend[nDatasets];//only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
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
  TString* pdfName_RelativeEfficiency_Data_zoom = new TString("RelativeEfficiency_Data_pt_zoom");
  TString* pdfName_RelativeEfficiency_MC_zoom = new TString("RelativeEfficiency_MC_pt_zoom");

  TString* pdfName_systUncertainty = new TString("Systematics_Efficiency_pt");
  TString* pdfName_systUncertainty_PreBarlow = new TString("Systematics_Efficiency_pt_PreBarlow");

  TString* pdfName_systUncertainty_zoom = new TString("Systematics_Efficiency_pt_zoom");
  TString* pdfName_systUncertainty_PreBarlow_zoom = new TString("Systematics_Efficiency_pt_PreBarlow_zoom");

  TString textContextData(contextCustomOneField(*texDatasetsComparisonType, ""));
  TString textContextMC(contextCustomOneField(*texDatasetsComparisonTypeMC, ""));


  // std::array<std::array<float, 2>, 2> drawnWindowLog = {{{(float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(1)
  //                                                         +hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(2))/2, 
  //                                                         (float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(hSystematicUncertainty[0]->GetNbinsX())
  //                                                         +hSystematicUncertainty[0]->GetXaxis()->GetBinWidth(hSystematicUncertainty[0]->GetNbinsX()))}
  //                                                         , {0, 0.025}}};


  std::array<std::array<float, 2>, 2> drawnWindow_uncertainty_zoom = {{{(float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(hSystematicUncertainty[0]->GetNbinsX())
                                                          +hSystematicUncertainty[0]->GetXaxis()->GetBinWidth(hSystematicUncertainty[0]->GetNbinsX()))}
                                                          , {0.97, 1.06}}};
  std::array<std::array<float, 2>, 2> drawnWindow_relativeefficiency_zoom = {{{(float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(hSystematicUncertainty[0]->GetXaxis()->GetBinLowEdge(hSystematicUncertainty[0]->GetNbinsX())
                                                          +hSystematicUncertainty[0]->GetXaxis()->GetBinWidth(hSystematicUncertainty[0]->GetNbinsX()))}
                                                          , {0.9, 1.05}}};

  Draw_TH1_Histograms(hRelativeEfficiency_Data, systematicsLegend, nDatasets-1, textContextData, pdfName_RelativeEfficiency_Data, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hRelativeEfficiency_MC, systematicsLegend, nDatasets-1, textContextMC, pdfName_RelativeEfficiency_MC, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hRelativeEfficiency_Data, systematicsLegend, nDatasets-1, textContextData, pdfName_RelativeEfficiency_Data_zoom, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindow_relativeefficiency_zoom, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hRelativeEfficiency_MC, systematicsLegend, nDatasets-1, textContextMC, pdfName_RelativeEfficiency_MC_zoom, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindow_relativeefficiency_zoom, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");


  Draw_TH1_Histograms(hSystematicUncertainty, systematicsLegend, nDatasets-1, textContextData, pdfName_systUncertainty, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hSystematicUncertainty_PreBarlow, systematicsLegend, nDatasets-1, textContextData, pdfName_systUncertainty_PreBarlow, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hSystematicUncertainty, systematicsLegend, nDatasets-1, textContextData, pdfName_systUncertainty_zoom, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindow_uncertainty_zoom, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
  Draw_TH1_Histograms(hSystematicUncertainty_PreBarlow, systematicsLegend, nDatasets-1, textContextData, pdfName_systUncertainty_PreBarlow_zoom, texPtReco, texRelativeEfficiencyRatio, texCollisionDataInfo, drawnWindow_uncertainty_zoom, legendPlacementAuto, contextPlacementAuto, "logx,ratioLine");
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// ITS-TPC matching efficiency systematic uncertainty ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <std::size_t N> 
void Get_TrackYield_Pt_ITSorTPCorITSTPC_MC(TH1D* (&hYieldReturned_mcprimary)[N], TH1D* (&hYieldReturned_mcsecondary)[N], int detectorId) {
  
  TH2D* H1D_trackPtTrackSigma_mcprimary[nDatasets];
  TH2D* H1D_trackPtTrackSigma_mcprimary_ptHigh[nDatasets];
  TH2D* H1D_trackPtTrackSigma_mcsecondary[nDatasets];
  TH2D* H1D_trackPtTrackSigma_mcsecondary_ptHigh[nDatasets];

  TH1D* H1D_trackPt_mcprimary[nDatasets];
  TH1D* H1D_trackPt_mcprimary_ptHigh[nDatasets];
  TH1D* H1D_trackPt_mcsecondary[nDatasets];
  TH1D* H1D_trackPt_mcsecondary_ptHigh[nDatasets];

  TH1D* H1D_trackPt_mcprimary_rebinned[nDatasets];
  TH1D* H1D_trackPt_mcprimary_ptHigh_rebinned[nDatasets];
  TH1D* H1D_trackPt_mcsecondary_rebinned[nDatasets];
  TH1D* H1D_trackPt_mcsecondary_ptHigh_rebinned[nDatasets];

  TH1D* H1D_trackPt_mcprimary_yield_concatenated[nDatasets];
  TH1D* H1D_trackPt_mcsecondary_yield_concatenated[nDatasets];

  int nBinsPt;
  std::vector<double> ptBinning;
  TString detector_suffix;
  if (detectorId == 0) {
    detector_suffix = "_ITS";
  } else if (detectorId == 1) {
    detector_suffix = "_TPC";
  } else if (detectorId == 2) {
    detector_suffix = "_ITSTPC";
  } else {
    cout << "detectorId wrong input, should be 0, 1 or 2" << endl;
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    cout << "dataset " << iDataset << endl;
    H1D_trackPtTrackSigma_mcprimary[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list_MC[iDataset]->Get(analysisWorkflowsMC[iDataset]+"/h2_track_pt_track_eta_mcprimary"+detector_suffix))->Clone("Get_TrackYield_Pt_ITSorTPCorITSTPC_MC_mcprimary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);
    H1D_trackPtTrackSigma_mcprimary_ptHigh[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list_MC[iDataset]->Get(analysisWorkflowsMC[iDataset]+"/h2_track_pt_high_track_sigmapt_mcprimary"+detector_suffix))->Clone("Get_TrackYield_Pt_ITSorTPCorITSTPC_MC_ptHigh_mcprimary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);
    H1D_trackPtTrackSigma_mcsecondary[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list_MC[iDataset]->Get(analysisWorkflowsMC[iDataset]+"/h2_track_pt_track_eta_mcsecondary"+detector_suffix))->Clone("Get_TrackYield_Pt_ITSorTPCorITSTPC_MC_mcsecondary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);
    H1D_trackPtTrackSigma_mcsecondary_ptHigh[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list_MC[iDataset]->Get(analysisWorkflowsMC[iDataset]+"/h2_track_pt_high_track_sigmapt_mcsecondary"+detector_suffix))->Clone("Get_TrackYield_Pt_ITSorTPCorITSTPC_MC_ptHigh_mcsecondary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);

    H1D_trackPt_mcprimary[iDataset] = (TH1D*)H1D_trackPtTrackSigma_mcprimary[iDataset]->ProjectionX("trackPt_tracks_mcprimary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, 0, -1, "e");
    H1D_trackPt_mcprimary_ptHigh[iDataset] = (TH1D*)H1D_trackPtTrackSigma_mcprimary_ptHigh[iDataset]->ProjectionX("trackPt_tracks_mcprimary_ptHigh"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, 0, -1, "e");
    H1D_trackPt_mcsecondary[iDataset] = (TH1D*)H1D_trackPtTrackSigma_mcsecondary[iDataset]->ProjectionX("trackPt_tracks_mcsecondary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, 0, -1, "e");
    H1D_trackPt_mcsecondary_ptHigh[iDataset] = (TH1D*)H1D_trackPtTrackSigma_mcsecondary_ptHigh[iDataset]->ProjectionX("trackPt_tracks_mcsecondary_ptHigh"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, 0, -1, "e");


    //tweaking the low-pt bins to make them larger close to 10GeV
    std::vector<double> xbinsVectorInitialLow = GetTH1Bins(H1D_trackPt_mcprimary[iDataset]);
    double* xbinsInitialLow = &xbinsVectorInitialLow[0];

    double ptBinsLowNew[500]; //500 to have a good margin
    double lastPt;
    int iBinPtNew = 0;
    int iBinPtInitialHisto = 0;
    int increment;
    float originalBinningEnd = 0.3; //above that we switch to custom binning
    while(iBinPtInitialHisto < H1D_trackPt_mcprimary[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialLow[iBinPtInitialHisto];
      if (lastPt < originalBinningEnd) { // top boundary for the unchanged binning
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_mcprimary[iDataset]->GetXaxis()->GetBinWidth(1) /8;
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
    H1D_trackPt_mcprimary_rebinned[iDataset] = (TH1D*)H1D_trackPt_mcprimary[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_mcprimary_rebinned"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, ptBinsLowNew);
    H1D_trackPt_mcsecondary_rebinned[iDataset] = (TH1D*)H1D_trackPt_mcsecondary[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_mcsecondary_rebinned"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, ptBinsLowNew);


    //tweaking the high-pt bins to have them less prone to statistical fluctuations
    std::vector<double> xbinsVectorInitialHigh = GetTH1Bins(H1D_trackPt_mcprimary_ptHigh[iDataset]);
    double* xbinsInitialHigh = &xbinsVectorInitialHigh[0];

    double ptBinsHighNew[500]; //500 to have a good margin
    iBinPtNew = 0;
    iBinPtInitialHisto = 0;
    while(iBinPtInitialHisto < H1D_trackPt_mcprimary_ptHigh[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialHigh[iBinPtInitialHisto];
      if (lastPt < 10) { // was 25 before
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_mcprimary_ptHigh[iDataset]->GetXaxis()->GetBinWidth(1) ;
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
    H1D_trackPt_mcprimary_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_mcprimary_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_mcprimary_ptHigh_rebinned"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, ptBinsHighNew);
    H1D_trackPt_mcsecondary_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_mcsecondary_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_mcsecondary_ptHigh_rebinned"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, ptBinsHighNew);


    // Merging high and low pt histograms:
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins(H1D_trackPt_mcprimary_rebinned[iDataset]);
    std::vector<double> xbinsVectorRight = GetTH1Bins(H1D_trackPt_mcprimary_ptHigh_rebinned[iDataset]);
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
    TH1D H1D_trackPt_mcprimary_yield_concatenated_temp("H1D_trackPt_mcprimary_yield_concatenated_temp"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, "H1D_trackPt_mcprimary_yield_concatenated_temp"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, nBinsPt, xbins_new);
    H1D_trackPt_mcprimary_yield_concatenated[iDataset] = (TH1D*)H1D_trackPt_mcprimary_yield_concatenated_temp.Clone("H1D_trackPt_mcprimary_yield_concatenated"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);
    TH1D H1D_trackPt_mcsecondary_yield_concatenated_temp("H1D_trackPt_mcsecondary_yield_concatenated_temp"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, "H1D_trackPt_mcsecondary_yield_concatenated_temp"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, nBinsPt, xbins_new);
    H1D_trackPt_mcsecondary_yield_concatenated[iDataset] = (TH1D*)H1D_trackPt_mcsecondary_yield_concatenated_temp.Clone("H1D_trackPt_mcsecondary_yield_concatenated"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);

    //filling the two new histograms
    for(int iBinX = 1; iBinX <= H1D_trackPt_mcprimary_rebinned[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_mcprimary_yield_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_mcprimary_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_mcprimary_yield_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_mcprimary_rebinned[iDataset]->GetBinError(iBinX));

      H1D_trackPt_mcsecondary_yield_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_mcsecondary_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_mcsecondary_yield_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_mcsecondary_rebinned[iDataset]->GetBinError(iBinX));
    }
    for(int iBinX = 1; iBinX <= H1D_trackPt_mcprimary_ptHigh_rebinned[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_mcprimary_yield_concatenated[iDataset]->SetBinContent(H1D_trackPt_mcprimary_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_mcprimary_ptHigh_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_mcprimary_yield_concatenated[iDataset]->SetBinError(H1D_trackPt_mcprimary_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_mcprimary_ptHigh_rebinned[iDataset]->GetBinError(iBinX));

      H1D_trackPt_mcsecondary_yield_concatenated[iDataset]->SetBinContent(H1D_trackPt_mcsecondary_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_mcsecondary_ptHigh_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_mcsecondary_yield_concatenated[iDataset]->SetBinError(H1D_trackPt_mcsecondary_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_mcsecondary_ptHigh_rebinned[iDataset]->GetBinError(iBinX));
    }
    
    hYieldReturned_mcprimary[iDataset] = (TH1D*)((TH1D*)H1D_trackPt_mcprimary_yield_concatenated[iDataset]->Clone("hYieldReturned_mcprimary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix));
    hYieldReturned_mcsecondary[iDataset] = (TH1D*)((TH1D*)H1D_trackPt_mcsecondary_yield_concatenated[iDataset]->Clone("hYieldReturned_mcsecondary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix));
  }
}

template <std::size_t N> 
void Get_TrackYield_Pt_ITSorTPCorITSTPC_Data(TH1D* (&hYieldReturned_inclusive)[N], int detectorId) {
  
  TH2D* H1D_trackPtTrackSigma_inclusive[nDatasets];
  TH2D* H1D_trackPtTrackSigma_inclusive_ptHigh[nDatasets];

  TH1D* H1D_trackPt_inclusive[nDatasets];
  TH1D* H1D_trackPt_inclusive_ptHigh[nDatasets];

  TH1D* H1D_trackPt_inclusive_rebinned[nDatasets];
  TH1D* H1D_trackPt_inclusive_ptHigh_rebinned[nDatasets];

  TH1D* H1D_trackPt_inclusive_yield_concatenated[nDatasets];

  int nBinsPt;
  std::vector<double> ptBinning;
  TString detector_suffix;
  if (detectorId == 0) {
    detector_suffix = "_ITS";
  } else if (detectorId == 1) {
    detector_suffix = "_TPC";
  } else if (detectorId == 2) {
    detector_suffix = "_ITSTPC";
  } else {
    cout << "detectorId wrong input, should be 0, 1 or 2" << endl;
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    cout << "dataset " << iDataset << endl;
    H1D_trackPtTrackSigma_inclusive[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list_MC[iDataset]->Get(analysisWorkflowsMC[iDataset]+"/h2_track_pt_track_eta_inclusive"+detector_suffix))->Clone("Get_TrackYield_Pt_ITSorTPCorITSTPC_Data"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);
    H1D_trackPtTrackSigma_inclusive_ptHigh[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list_MC[iDataset]->Get(analysisWorkflowsMC[iDataset]+"/h2_track_pt_high_track_sigmapt_inclusive"+detector_suffix))->Clone("Get_TrackYield_Pt_ITSorTPCorITSTPC_Data_ptHigh"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);

    H1D_trackPt_inclusive[iDataset] = (TH1D*)H1D_trackPtTrackSigma_inclusive[iDataset]->ProjectionX("trackPt_tracks_inclusive"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, 0, -1, "e");
    H1D_trackPt_inclusive_ptHigh[iDataset] = (TH1D*)H1D_trackPtTrackSigma_inclusive_ptHigh[iDataset]->ProjectionX("trackPt_tracks_inclusive_ptHigh"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, 0, -1, "e");

    //tweaking the low-pt bins to make them larger close to 10GeV
    std::vector<double> xbinsVectorInitialLow = GetTH1Bins(H1D_trackPt_inclusive[iDataset]);
    double* xbinsInitialLow = &xbinsVectorInitialLow[0];

    double ptBinsLowNew[500]; //500 to have a good margin
    double lastPt;
    int iBinPtNew = 0;
    int iBinPtInitialHisto = 0;
    int increment;
    float originalBinningEnd = 0.3; //above that we switch to custom binning
    while(iBinPtInitialHisto < H1D_trackPt_inclusive[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialLow[iBinPtInitialHisto];
      if (lastPt < originalBinningEnd) { // top boundary for the unchanged binning
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_inclusive[iDataset]->GetXaxis()->GetBinWidth(1) /8;
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
    H1D_trackPt_inclusive_rebinned[iDataset] = (TH1D*)H1D_trackPt_inclusive[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_inclusive_rebinned"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, ptBinsLowNew);


    //tweaking the high-pt bins to have them less prone to statistical fluctuations
    std::vector<double> xbinsVectorInitialHigh = GetTH1Bins(H1D_trackPt_inclusive_ptHigh[iDataset]);
    double* xbinsInitialHigh = &xbinsVectorInitialHigh[0];

    double ptBinsHighNew[500]; //500 to have a good margin
    iBinPtNew = 0;
    iBinPtInitialHisto = 0;
    while(iBinPtInitialHisto < H1D_trackPt_inclusive_ptHigh[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialHigh[iBinPtInitialHisto];
      if (lastPt < 10) { // was 25 before
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_inclusive_ptHigh[iDataset]->GetXaxis()->GetBinWidth(1) ;
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
    H1D_trackPt_inclusive_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_inclusive_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_inclusive_ptHigh_rebinned"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, ptBinsHighNew);


    // Merging high and low pt histograms:
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins(H1D_trackPt_inclusive_rebinned[iDataset]);
    std::vector<double> xbinsVectorRight = GetTH1Bins(H1D_trackPt_inclusive_ptHigh_rebinned[iDataset]);
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
    TH1D H1D_trackPt_inclusive_yield_concatenated_temp("H1D_trackPt_inclusive_yield_concatenated_temp"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, "H1D_trackPt_inclusive_yield_concatenated_temp"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix, nBinsPt, xbins_new);
    H1D_trackPt_inclusive_yield_concatenated[iDataset] = (TH1D*)H1D_trackPt_inclusive_yield_concatenated_temp.Clone("H1D_trackPt_inclusive_yield_concatenated"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix);

    //filling the two new histograms
    for(int iBinX = 1; iBinX <= H1D_trackPt_inclusive_rebinned[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_inclusive_yield_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_inclusive_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_inclusive_yield_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_inclusive_rebinned[iDataset]->GetBinError(iBinX));
    }
    for(int iBinX = 1; iBinX <= H1D_trackPt_inclusive_ptHigh_rebinned[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_inclusive_yield_concatenated[iDataset]->SetBinContent(H1D_trackPt_inclusive_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_inclusive_ptHigh_rebinned[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_inclusive_yield_concatenated[iDataset]->SetBinError(H1D_trackPt_inclusive_rebinned[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_inclusive_ptHigh_rebinned[iDataset]->GetBinError(iBinX));
    }
    
    hYieldReturned_inclusive[iDataset] = (TH1D*)((TH1D*)H1D_trackPt_inclusive_yield_concatenated[iDataset]->Clone("hYieldReturned_inclusive"+DatasetsMC[iDataset]+DatasetsNames[iDataset]+detector_suffix));
  }
}

template <std::size_t N> 
void Get_ITSTPC_matching_efficiency_MC(TH1D* (&hMatchingEfficiencyReturned_mcprimary)[N], TH1D* (&hMatchingEfficiencyReturned_mcsecondary)[N]) {
  
  TH1D* H1D_trackPt_mcprimary_yield_ITSTPC[nDatasets];
  TH1D* H1D_trackPt_mcprimary_yield_TPC[nDatasets];
  TH1D* H1D_trackPt_mcsecondary_yield_ITSTPC[nDatasets];
  TH1D* H1D_trackPt_mcsecondary_yield_TPC[nDatasets];

  int detectorIdITSTPC = 2;
  Get_TrackYield_Pt_ITSorTPCorITSTPC_MC(H1D_trackPt_mcprimary_yield_ITSTPC, H1D_trackPt_mcsecondary_yield_ITSTPC, detectorIdITSTPC);
  int detectorIdTPC = 1;
  Get_TrackYield_Pt_ITSorTPCorITSTPC_MC(H1D_trackPt_mcprimary_yield_TPC, H1D_trackPt_mcsecondary_yield_TPC, detectorIdTPC);
  
  bool divideSuccessMcPrimary, divideSuccessMcSecondary;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    cout << "dataset " << iDataset << endl;

    hMatchingEfficiencyReturned_mcprimary[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPt_mcprimary_yield_TPC[iDataset]->Clone("hMatchingEfficiencyReturned_mcprimary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]));
    hMatchingEfficiencyReturned_mcprimary[iDataset-1]->Reset("M");
    divideSuccessMcPrimary = hMatchingEfficiencyReturned_mcprimary[iDataset-1]->Divide(H1D_trackPt_mcprimary_yield_ITSTPC[iDataset], H1D_trackPt_mcprimary_yield_TPC[iDataset], 1, 1, "b"); // option b because numerator is subset of denominator
    if (!divideSuccessMcPrimary) {
      cout << "divideSuccessMcPrimary = false: divide failed in Get_ITSTPC_matching_efficiency_MC(), iDataset = " << iDataset << endl;
    }

    hMatchingEfficiencyReturned_mcsecondary[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPt_mcsecondary_yield_TPC[iDataset]->Clone("hMatchingEfficiencyReturned_mcsecondary"+DatasetsMC[iDataset]+DatasetsNames[iDataset]));
    hMatchingEfficiencyReturned_mcsecondary[iDataset-1]->Reset("M");
    divideSuccessMcSecondary = hMatchingEfficiencyReturned_mcsecondary[iDataset-1]->Divide(H1D_trackPt_mcsecondary_yield_ITSTPC[iDataset], H1D_trackPt_mcsecondary_yield_TPC[iDataset], 1, 1, "b"); // option b because numerator is subset of denominator
    if (!divideSuccessMcSecondary) {
      cout << "divideSuccessMcSecondary = false: divide failed in Get_ITSTPC_matching_efficiency_MC(), iDataset = " << iDataset << endl;
    }
  }
}

template <std::size_t N> 
void Get_ITSTPC_matching_efficiency_Data(TH1D* (&hMatchingEfficiencyReturned_inclusive)[N]) {
  
  TH1D* H1D_trackPt_inclusive_yield_ITSTPC[nDatasets];
  TH1D* H1D_trackPt_inclusive_yield_TPC[nDatasets];

  int detectorIdITSTPC = 2;
  Get_TrackYield_Pt_ITSorTPCorITSTPC_Data(H1D_trackPt_inclusive_yield_ITSTPC, detectorIdITSTPC);
  int detectorIdTPC = 1;
  Get_TrackYield_Pt_ITSorTPCorITSTPC_Data(H1D_trackPt_inclusive_yield_TPC, detectorIdTPC);
  
  bool divideSuccessData;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    cout << "dataset " << iDataset << endl;

    hMatchingEfficiencyReturned_inclusive[iDataset-1] = (TH1D*)((TH1D*)H1D_trackPt_inclusive_yield_TPC[iDataset]->Clone("hMatchingEfficiencyReturned_inclusive"+DatasetsData[iDataset]+DatasetsNames[iDataset]));
    hMatchingEfficiencyReturned_inclusive[iDataset-1]->Reset("M");
    divideSuccessData = hMatchingEfficiencyReturned_inclusive[iDataset-1]->Divide(H1D_trackPt_inclusive_yield_ITSTPC[iDataset], H1D_trackPt_inclusive_yield_TPC[iDataset], 1, 1, "b"); // option b because numerator is subset of denominator
    if (!divideSuccessData) {
      cout << "divideSuccessData = false: divide failed in Get_ITSTPC_matching_efficiency_Data(), iDataset = " << iDataset << endl;
    }
  }
}

template <std::size_t N> 
void Get_DCAxy_yields(TH1D* (&hYieldReturned)[N], TFile* (&fileInput)[N], const TString (&analysisWorkflows)[N], std::string particleStatus, double* ptRange) {
  TH2D* H2D_trackPt_trackDcaXY_ptLow[nDatasets];
  TH2D* H2D_trackPt_trackDcaXY_ptHigh[nDatasets];
  TH1D* hYieldReturned_temp[nDatasets];

  double ptLimit_lowHigh = 10;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    TString histSpecifier = (TString)DatasetsNames[iDataset]+(TString)particleStatus+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]";

    cout << "dataset " << iDataset << endl;
    H2D_trackPt_trackDcaXY_ptLow[iDataset] = (TH2D*)((TH2D*)fileInput[iDataset]->Get(analysisWorkflows[iDataset]+"/h_track_pt_track_dcaxy_"+(TString)particleStatus))->Clone("H2D_trackPt_trackDcaXY_ptLow"+(TString)histSpecifier);
    H2D_trackPt_trackDcaXY_ptLow[iDataset]->Sumw2();
    H2D_trackPt_trackDcaXY_ptHigh[iDataset] = (TH2D*)((TH2D*)fileInput[iDataset]->Get(analysisWorkflows[iDataset]+"/h_track_pt_high_track_dcaxy_"+(TString)particleStatus))->Clone("H2D_trackPt_trackDcaXY_ptHigh"+(TString)histSpecifier);
    H2D_trackPt_trackDcaXY_ptHigh[iDataset]->Sumw2();

    if (ptRange[1] < ptLimit_lowHigh) {
      int ibinPt_low = H2D_trackPt_trackDcaXY_ptLow[iDataset]->GetXaxis()->FindBin(ptRange[0]);
      int ibinPt_high = H2D_trackPt_trackDcaXY_ptLow[iDataset]->GetXaxis()->FindBin(ptRange[1]);
      hYieldReturned_temp[iDataset] = (TH1D*)H2D_trackPt_trackDcaXY_ptLow[iDataset]->ProjectionY("H1D_trackDcaXY"+histSpecifier, ibinPt_low, ibinPt_high, "e");

    } else if (ptLimit_lowHigh < ptRange[0]) {
      int ibinPt_low = H2D_trackPt_trackDcaXY_ptHigh[iDataset]->GetXaxis()->FindBin(ptRange[0]);
      int ibinPt_high = H2D_trackPt_trackDcaXY_ptHigh[iDataset]->GetXaxis()->FindBin(ptRange[1]);
      hYieldReturned_temp[iDataset] = (TH1D*)H2D_trackPt_trackDcaXY_ptHigh[iDataset]->ProjectionY("H1D_trackDcaXY_ptHigh"+histSpecifier, ibinPt_low, ibinPt_high, "e");
    } else {
      // ptLimit_lowHigh is inside the ptRange
      int ibinPt_low = H2D_trackPt_trackDcaXY_ptLow[iDataset]->GetXaxis()->FindBin(ptRange[0]);
      int ibinPt_high = H2D_trackPt_trackDcaXY_ptHigh[iDataset]->GetXaxis()->FindBin(ptRange[1]);
      hYieldReturned_temp[iDataset] = (TH1D*)H2D_trackPt_trackDcaXY_ptLow[iDataset]->ProjectionY("H1D_trackDcaXY"+histSpecifier, ibinPt_low, H2D_trackPt_trackDcaXY_ptLow[iDataset]->GetNbinsX(), "e");
      hYieldReturned_temp[iDataset]->Add((TH1D*)H2D_trackPt_trackDcaXY_ptHigh[iDataset]->ProjectionY("H1D_trackDcaXY_ptHigh"+histSpecifier, 1, ibinPt_high, "e"));
    }
    hYieldReturned[iDataset] = (TH1D*)hYieldReturned_temp[iDataset]->Rebin(1.,"hYieldReturned_temp"+histSpecifier);
  }
}

template <std::size_t N> 
void Get_PrimarySecondaryFractions_in_Data(double (&primaryFractionsInData)[N], double (&secondaryFractionsInData)[N], double* ptRange){
  // https://root.cern/doc/v628/classTFractionFitter.html
  // gotta be careful with errors using this package: https://arxiv.org/pdf/0803.2711

  TH1D *H1D_DCAxy_data[nDatasets];
  TH1D *H1D_DCAxy_mcprimary[nDatasets];
  TH1D *H1D_DCAxy_mcsecondary[nDatasets];

  std::string particleStatusInclusive("datainclusive");
  std::string particleStatusMcPrimaries("mcprimary");
  std::string particleStatusMcSecondaries("mcsecondary");
  Get_DCAxy_yields(H1D_DCAxy_data, file_O2Analysis_list_Data, analysisWorkflowsData, particleStatusInclusive, ptRange);
  Get_DCAxy_yields(H1D_DCAxy_mcprimary, file_O2Analysis_list_MC, analysisWorkflowsMC, particleStatusMcPrimaries, ptRange);
  Get_DCAxy_yields(H1D_DCAxy_mcsecondary, file_O2Analysis_list_MC, analysisWorkflowsMC, particleStatusMcSecondaries, ptRange);


  int iDataset = 0; // test 
  int minBin = H1D_DCAxy_data[iDataset]->FindFirstBinAbove(0);
  int maxBin = H1D_DCAxy_data[iDataset]->FindLastBinAbove(0);
  cout << "minBin = " << minBin << ", maxBin = " << maxBin << endl;

  // H1D_DCAxy_mcprimary[iDataset]->Scale(1./H1D_DCAxy_mcprimary[iDataset]->GetBinContent(minBin+5));
  // H1D_DCAxy_mcsecondary[iDataset]->Scale(1./H1D_DCAxy_mcsecondary[iDataset]->GetBinContent(minBin+5));
  // H1D_DCAxy_data[iDataset]->Scale(1./H1D_DCAxy_data[iDataset]->GetBinContent(minBin+5));

  H1D_DCAxy_data[iDataset]->Scale(1./H1D_DCAxy_data[iDataset]->GetEntries());
  H1D_DCAxy_mcprimary[iDataset]->Scale(1./H1D_DCAxy_mcprimary[iDataset]->GetEntries());
  H1D_DCAxy_mcsecondary[iDataset]->Scale(1./H1D_DCAxy_mcsecondary[iDataset]->GetEntries());

  // H1D_DCAxy_mcprimary[iDataset]->Scale(H1D_DCAxy_data[iDataset]->GetEntries()*1./H1D_DCAxy_mcprimary[iDataset]->GetEntries());
  // H1D_DCAxy_mcsecondary[iDataset]->Scale(H1D_DCAxy_data[iDataset]->GetEntries()*1./H1D_DCAxy_mcsecondary[iDataset]->GetEntries());


  TVirtualFitter::SetMaxIterations(10000);

  // retrieve histograms
  TObjArray* TObjArray_DCAxy_mc[nDatasets];
  TObjArray_DCAxy_mc[iDataset] = new TObjArray(2);        // MC histograms are put in this array
  TObjArray_DCAxy_mc[iDataset]->Add(H1D_DCAxy_mcprimary[iDataset]);
  TObjArray_DCAxy_mc[iDataset]->Add(H1D_DCAxy_mcsecondary[iDataset]);
  TFractionFitter* templateFit[nDatasets];
  templateFit[iDataset] = new TFractionFitter(H1D_DCAxy_data[iDataset], TObjArray_DCAxy_mc[iDataset]); // initialise
  // templateFit[iDataset]->Constrain(0, 0., 1.);               // constrain fraction 1 to be between 0 and 1
  templateFit[iDataset]->SetRangeX(minBin, maxBin);                   
  // templateFit[iDataset]->Constrain(1, 0., 1.);               // constrain fraction 1 to be between 0 and 1
  // ROOT::Fit::Fitter* fitter = fit->GetFitter();
  // fitter->Config().ParSettings(parameter #).Set(const std::string &name, double value, double step, double lower, double upper);

  // Int_t status = templateFit[iDataset]->Fit();               // perform the fit
  auto status = templateFit[iDataset]->Fit();               // perform the fit
   status = templateFit[iDataset]->Fit();               // perform the fit
   status = templateFit[iDataset]->Fit();               // perform the fit
   status = templateFit[iDataset]->Fit();               // perform the fit
   status = templateFit[iDataset]->Fit();               // perform the fit
   status = templateFit[iDataset]->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
  // if (status == 0) {                       // check on fit status
    auto matrix = status->GetCovarianceMatrix();
    matrix.Print();
    // the parameters of the fit are the fractions
    double fractionsArrayValue[2];
    double fractionsArrayError[2];
    templateFit[iDataset]->GetResult(0, fractionsArrayValue[0], fractionsArrayError[0]);
    templateFit[iDataset]->GetResult(1, fractionsArrayValue[1], fractionsArrayError[1]);
    cout << "fractionPrimaries = " << fractionsArrayValue[0] << ", fractionSecondaries = " << fractionsArrayValue[1] << endl;
    primaryFractionsInData[iDataset] = fractionsArrayValue[0];
    secondaryFractionsInData[iDataset] = fractionsArrayValue[1];

    TString systematicsLegend[4];//only 0 -> n-1 will be used, the n is only there (instead of n-1) so that this function is accepted when nDatasets=1
    systematicsLegend[0] = (TString)"data";
    systematicsLegend[1] = (TString)"mcprimaries";
    systematicsLegend[2] = (TString)"mcseconcaries";
    systematicsLegend[3] = (TString)"mcweightedsum";

    TString* pdfName = new TString("DcaTemplateFit_"+DatasetsNames[iDataset]);
    TString textContextData(contextCustomOneField(*texDatasetsComparisonType, ""));

    TH1D* DataAndMcTemplates[3]; // 0 data, 1 mc primaries, 2 mc secondaries
    DataAndMcTemplates[0] = (TH1D*)H1D_DCAxy_data[iDataset]->Clone("DataAndMcTemplatesData");
    // DataAndMcTemplates[0]->Scale(1./DataAndMcTemplates[0]->GetBinContent(minBin+5));
    DataAndMcTemplates[1] = (TH1D*)H1D_DCAxy_mcprimary[iDataset]->Clone("DataAndMcTemplatesMCPrim");
    // DataAndMcTemplates[1]->Scale(1./DataAndMcTemplates[1]->GetBinContent(minBin+5));
    // DataAndMcTemplates[1]->Scale(fractionsArrayValue[0]);
    DataAndMcTemplates[1]->Scale(0.6); // 0.85 and 0.15 works very well for 0.3 to 100gev
    DataAndMcTemplates[2] = (TH1D*)H1D_DCAxy_mcsecondary[iDataset]->Clone("DataAndMcTemplatesMCSec");
    // DataAndMcTemplates[2]->Scale(1./DataAndMcTemplates[2]->GetBinContent(minBin+5));
    // DataAndMcTemplates[2]->Scale(fractionsArrayValue[1]);
    DataAndMcTemplates[2]->Scale(0.4);

    DataAndMcTemplates[3] = (TH1D*)DataAndMcTemplates[1]->Clone("DataAndMcTemplatesMCSum");
    DataAndMcTemplates[3]->Add(DataAndMcTemplates[2]);

    TString* pdfName_2 = new TString("DcaTemplateFit_RatioDataToSum"+DatasetsNames[iDataset]);
    TH1D* DataToSumRatio = (TH1D*)DataAndMcTemplates[0]->Clone("DataToSumRatio");
    DataToSumRatio->Divide(DataAndMcTemplates[3]);
    // Data = MC1*[frac1*data.Integral()/mc1.Integral()] + MC2*[frac2*data.Integral()/mc2.Integral()] // if not initially normalised to integral of 1, should likely do this



    Draw_TH1_Histograms(DataAndMcTemplates, systematicsLegend, 4, textContextData, pdfName, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
    Draw_TH1_Histogram(DataToSumRatio, textContextData, pdfName_2, texPtReco, texRelativeEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy, ratioLine");
  // } else {
  //   cout << "template fit failed, status = " << status << endl;
  // }
  // this function scales those templates to get a fraction
}


//Minuit 2:5.1.2 fval(), edm(), nfcn()
// The method double FunctionMinimum::fval() returns the function value at the minimum, 
// the method double FunctionMinimum::edm() returns the Expected vertical Distance to the Minimum   
// and unsigned int FunctionMinimum::nfcn() returns the total number of function calls during the minimization.
