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
#include <RooUnfold.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"

//My Libraries
#include "./JetSpectrum_settings.h"
#include "./JetSpectrum_inputs.h"

#include "./JetSpectrum_ResponseMatrixFunctions.h"
#include "./JetSpectrum_ResponseMatrixFunctions.C"
#include "./JetSpectrum_SpectraGetters.h"
#include "./JetSpectrum_SpectraGetters.C"
#include "./JetSpectrum_Unfolding.h"
#include "./JetSpectrum_Unfolding.C"
#include "./JetSpectrum_EfficiencyPurityGetters.h"
#include "./JetSpectrum_EfficiencyPurityGetters.C"

#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" 
#include "../Utilities/HistogramUtilities.C"
#include "../Utilities/HistogramPlotting.C" 

#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();



void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);
void Draw_Systematics_TrackEff(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);


/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_systematics() {
  // Load necessary libraries
  LoadLibs_Systematics();
  // Set the default style
  SetStyle_Systematics();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // gathers the analysis options in a single char[]
  mcCollHistIsObsolete = inputMcCollHistIsObsolete;
  /// Do not run all functions together, select which one to run by commenting/uncommenting
  int iDataset = 0;
  int iRadius = 1;

  //######################################################### Unfolding method Systematics #####################################################
  // char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  // snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);
  // const int nUnfoldingMethods = 2;
  // char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes"}; // default is the first one in this list
  // int unfoldParameterInputList[2] = {8, 4};
  // Draw_Systematics_UnfoldMethod(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);

  //######################################################### Parameter variation Systematics #####################################################
  // char optionsAnalysis[100] = "";
  // snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", unfoldingPrior, unfoldingMethod);
  // int unfoldParameterInputMin = 7;
  // int unfoldParameterInputMax = 9;
  // int unfoldParameterInputStep = 1;
  // Draw_Systematics_parameterVariation(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);

  //######################################################### Track efficiency Systematics #####################################################
  char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);
  const int nUnfoldingMethods = 4;
  char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes", "Svd", "Bayes"}; // first two to be with nominal efficiency, last two with efficiency varied 
  int unfoldParameterInputList[4] = {7, 4, 4, 2}; // first two to be with nominal efficiency, last two with efficiency varied
  Draw_Systematics_TrackEff(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);




}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs_Systematics() {
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

void SetStyle_Systematics(Bool_t graypalette) {
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


void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {

  TH1D* hTempSystematicUncertainty = new TH1D("hTempSystematicUncertainty", "hTempSystematicUncertainty", nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]);
  TH1D* hTempSystematicUncertainty_PreBarlow = new TH1D("hTempSystematicUncertainty_PreBarlow", "hTempSystematicUncertainty_PreBarlow", nBinPtJetsGen[iRadius], ptBinsJetsGen[iRadius]);
  hTempSystematicUncertainty->Sumw2();
  hTempSystematicUncertainty_PreBarlow->Sumw2();
  hTempSystematicUncertainty->Reset("M");
  hTempSystematicUncertainty_PreBarlow->Reset("M");
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // return histogram that has the systematics in its contents
  TH1D* H1D_jetPt_unfolded[nUnfoldingMethods];
  TH1D* H1D_jetPt_unfolded_differences[nUnfoldingMethods-1];


  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    }
  }

  char optionsAnalysis_withUnfoldingMethod[100] = "";
  for(int iMethod = 0; iMethod < nUnfoldingMethods; iMethod++){
    snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDataset, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);

    if (iMethod != 0) {
      H1D_jetPt_unfolded_differences[iMethod-1] = (TH1D*)H1D_jetPt_unfolded[iMethod]->Clone("H1D_jetPt_unfolded_differences"+partialUniqueSpecifier);
      H1D_jetPt_unfolded_differences[iMethod-1]->Add(H1D_jetPt_unfolded[0],-1);
    }
    cout << "do I want the absolute value of the difference?" << endl;
  }

  // cout << "Do I apply Barlow condition even though not param variation ? check paper again" << endl; YES, the subset of data thing is only shown for first demonstration, but barlow says it holds true even if that's not the cast

  /////////////////
  // Barlow test //
  /////////////////

  TH1D* H1D_jetPt_unfolded_REF = H1D_jetPt_unfolded[0];
  double SystUncertainty;
  int id_SignalExtractionType_maxDeviation;
  double hSigmaBarlow[nBinPtJetsGen[iRadius]];
  for(int iBinPt = 1; iBinPt <= nBinPtJetsGen[iRadius]; iBinPt++){
    SystUncertainty = 0;
    for(int iMethod = 1; iMethod < nUnfoldingMethods; iMethod++){ // get maximum difference among the nUnfoldingMethods-1 ones, hold value with SystUncertainty, and the id of the method wîth id_SignalExtractionType_maxDeviation
      if (abs(H1D_jetPt_unfolded_differences[iMethod-1]->GetBinContent(iBinPt)) > SystUncertainty) {
        SystUncertainty = abs(H1D_jetPt_unfolded_differences[iMethod-1]->GetBinContent(iBinPt));
        id_SignalExtractionType_maxDeviation = iMethod;
      }
    }

    // Barlow condition for systematics (Systematic Errors: facts and fictions, by Roger Barlow, https://arxiv.org/abs/hep-ex/0207026)
    Double_t StatUncertainty_REF = H1D_jetPt_unfolded_REF->GetBinError(iBinPt);
    Double_t StatUncertainty_MaxDeviationCase = H1D_jetPt_unfolded[id_SignalExtractionType_maxDeviation]->GetBinError(iBinPt);
 
    int PtArrayIterator = iBinPt - 1;
    hSigmaBarlow[PtArrayIterator] = sqrt(abs(StatUncertainty_MaxDeviationCase*StatUncertainty_MaxDeviationCase - StatUncertainty_REF*StatUncertainty_REF)); //stat error of the difference in the case of subsample
 
    hTempSystematicUncertainty_PreBarlow->SetBinContent(iBinPt,SystUncertainty);
    hTempSystematicUncertainty_PreBarlow->SetBinError(iBinPt,hSigmaBarlow[PtArrayIterator]);

    if (SystUncertainty > N_SigmaBarlow*hSigmaBarlow[PtArrayIterator]) { //Could ask for 1Sigma, 4Sigma or whatever depending on how conservative we want to be; one suggested in PWGLF note is 2Sigma
      hTempSystematicUncertainty->SetBinContent(iBinPt,SystUncertainty);
      hTempSystematicUncertainty->SetBinError(iBinPt,hSigmaBarlow[PtArrayIterator]);
    }
    else {
      hTempSystematicUncertainty->SetBinContent(iBinPt,0.);
    }
  }

  hSystematicUncertainty = (TH1D*)hTempSystematicUncertainty->Clone("hSystematicUncertainty_UnfoldMethod"+partialUniqueSpecifier);
  hSystematicUncertainty_PreBarlow = (TH1D*)hTempSystematicUncertainty_PreBarlow->Clone("hSystematicUncertainty_PreBarlow_UnfoldMethod"+partialUniqueSpecifier);

}

void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {

  TH1D* hSystematicUncertainty;
  TH1D* hSystematicUncertainty_PreBarlow;
  Get_systematics_UnfoldMethod(hSystematicUncertainty, hSystematicUncertainty_PreBarlow, iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, options);

  TH1D* H1D_jetPt_unfolded;
  char optionsAnalysis_withUnfoldingMethod[100] = "";
  snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethod);

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    }
  }
  
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInputList[0], optionsAnalysis_withUnfoldingMethod); //SVD unfolding as reference
  hSystematicUncertainty->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  hSystematicUncertainty->Scale(100.0);
  hSystematicUncertainty_PreBarlow->Divide(H1D_jetPt_unfolded); //get it as a ratio of ref corrected yield
  hSystematicUncertainty_PreBarlow->Scale(100.0);


  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius])+"]_"+unfoldingMethodList[0]+"_kUnfold="+Form("%i", unfoldParameterInputList[0]);

  TString* pdfName = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier);
  TString* pdfName_PreBarlow = new TString("Systematics_UnfoldMethod_"+partialUniqueSpecifier+"_PreBarlow");

  // TString textContext("");
  TString textContext = Form(
    "#splitline{sys. unfolding method}"
    "{k_{svd} = %d, k_{bayes} = %d}",
    unfoldParameterInputList[0],
    unfoldParameterInputList[1]
  );

  TString* texRelativeErrPercent = new TString ("relative error (%)");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 25}}};
  Draw_TH1_Histogram(hSystematicUncertainty, textContext, pdfName, texPtJetRec, texRelativeErrPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSystematicUncertainty_PreBarlow, textContext, pdfName_PreBarlow, texPtJetRec, texRelativeErrPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
}

void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options) {
  cout << "########### Drawing systematics from parameter variation ###############" << endl;
  const int nUnfoldIteration = std::floor((unfoldIterationMax - unfoldIterationMin + 1)/step);

  TH1D* H1D_jetPt_unfolded[nUnfoldIteration];

  TString partialUniqueSpecifier;

  partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  int unfoldParameterInput = 0;

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    }
  }

  if (measuredInput == nullptr) {
    cout << "Error: measuredInput histogram is null!" << endl;
    return;
  }
  else {
    cout << "measuredInput histogram successfully retrieved." << endl;
  }

  for(int iUnfoldIteration = 0; iUnfoldIteration < nUnfoldIteration; iUnfoldIteration++){
    cout << "  entering the for loop "  << endl;
    unfoldParameterInput = unfoldIterationMax - iUnfoldIteration * step; 

    cout << "((((((((((((()))))))))))))" << endl;
    cout << "Iteration "<< iUnfoldIteration << endl;
    cout << "((((((((((((()))))))))))))" << endl;
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iUnfoldIteration], measuredInput, iDataset, iRadius, unfoldParameterInput, options);
  }

  int iNominal = nUnfoldIteration/2; // choose nominal (e.g. central iteration)
  TH1D* hNom = H1D_jetPt_unfolded[iNominal];

  TH1D* hSys_envelope = (TH1D*)hNom->Clone("hSys_envelope");
  hSys_envelope->Reset(); // will store absolute systematic (positive)

  int nBins = hNom->GetNbinsX();
  for (int ib = 1; ib <= nBins; ++ib) {
      double valNom = hNom->GetBinContent(ib);
      double maxAbs = 0.0;
      for (int i = 0; i < nUnfoldIteration; ++i) {
          double val = H1D_jetPt_unfolded[i]->GetBinContent(ib);
          double d = fabs(val - valNom);
          if (d > maxAbs) maxAbs = d;
      }
      hSys_envelope->SetBinContent(ib, maxAbs);
  }

  TH1D* hSys_rel = (TH1D*)hSys_envelope->Clone("hSys_relative");
  hSys_rel->Reset();

  for (int ib = 1; ib <= nBins; ++ib) {
      double absSys = hSys_envelope->GetBinContent(ib);
      double valNom = hNom->GetBinContent(ib);

      double rel = 0.0;
      if (valNom > 0) rel = absSys / valNom;

      hSys_rel->SetBinContent(ib, rel*100.0); // in percent
  }
  TString* pdfName_envolope = new TString("Systematics_deviation_ParameterVariation_"+partialUniqueSpecifier);
  TString* pdfName_relUnc = new TString("Systematics_RelativeUncertainty_ParameterVariation_"+partialUniqueSpecifier);
  TString textContext("SVD unfolding");
  TString* sigma = new TString("#sigma interation variation");
  TString* relativeErrors = new TString("realtive errors (%) ");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 1.6}}};
  // Draw_TH1_Histogram(hSys_envelope, textContext, pdfName_envolope, texPtJetRecX, sigma, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(hSys_rel, textContext, pdfName_relUnc, texPtJetRec, relativeErrors, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  
}

void Draw_Systematics_TrackEff(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options) {
  cout << "########### Drawing systematics from track efficiency variation ###############" << endl;
  TString partialUniqueSpecifier = Datasets[iDataset]+"_R="+Form("%.1f",arrayRadius[iRadius]);

  // return histogram that has the systematics in its contents
  TH1D* H1D_jetPt_unfolded[nUnfoldingMethods];

  TH1D* measuredInput;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options); 
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measuredInput, iDataset, iRadius, options);
    }
  } else{
    Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    if (useFineBinningTest) {
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measuredInput, iDataset, iRadius, options);
    }
  }

  char optionsAnalysis_withUnfoldingMethod[100] = "";
  for(int iMethod = 0; iMethod < nUnfoldingMethods; iMethod++){
    snprintf(optionsAnalysis_withUnfoldingMethod, sizeof(optionsAnalysis_withUnfoldingMethod), "%s,%s", options.c_str(), (const char*)unfoldingMethodList[iMethod]);
    int iDatasetMC = (iMethod == 0 || iMethod == 1) ? 0 : 1; // first two to be with nominal efficiency dataset 0, last two with efficiency varied dataset 1
    Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded[iMethod], measuredInput, iDatasetMC, iRadius, unfoldParameterInputList[iMethod], optionsAnalysis_withUnfoldingMethod);
  }

  // Get histograms
  TH1D* H1D_Nominal_SVD = (TH1D*) H1D_jetPt_unfolded[0]->Clone("H1D_Nominal_SVD");
  TH1D* H1D_Nominal_Bayes = (TH1D*) H1D_jetPt_unfolded[1]->Clone("H1D_Nominal_Bayes");
  TH1D* H1D_Reduced_SVD = (TH1D*) H1D_jetPt_unfolded[2]->Clone("H1D_Reduced_SVD");
  TH1D* H1D_Reduced_Bayes = (TH1D*) H1D_jetPt_unfolded[3]->Clone("H1D_Reduced_Bayes");

  // Create histograms for absolute differences and relative uncertainties
  TH1D *H1D_Delta_SVD = (TH1D*)H1D_Nominal_SVD->Clone("hDiff_NominalReduced_SVD");
  H1D_Delta_SVD->Reset();
  TH1D *H1D_RelativeUncertainty_SVD = (TH1D*)H1D_Nominal_SVD->Clone("H1D_RelativeUncertainty_SVD");
  H1D_RelativeUncertainty_SVD->Reset();

  TH1D *H1D_Delta_Bayes = (TH1D*)H1D_Nominal_Bayes->Clone("H1D_Delta_Bayes");
  H1D_Delta_Bayes->Reset();
  TH1D *H1D_RelativeUncertainty_Bayes = (TH1D*)H1D_Nominal_Bayes->Clone("H1D_RelativeUncertainty_Bayes");
  H1D_RelativeUncertainty_Bayes->Reset();

  // Compute absolute differences bin by bin and compute relative uncertainties
  for (int i = 1; i <= H1D_Nominal_SVD->GetNbinsX(); ++i) {
      H1D_Delta_SVD->SetBinContent(i, abs(H1D_Nominal_SVD->GetBinContent(i) - H1D_Reduced_SVD->GetBinContent(i)));
      H1D_Delta_SVD->SetBinError(i, 0.0);

      if (H1D_Nominal_SVD->GetBinContent(i) != 0) {
          double relUnc = (H1D_Delta_SVD->GetBinContent(i) / H1D_Nominal_SVD->GetBinContent(i)) * 100.0; // in percent
          H1D_RelativeUncertainty_SVD->SetBinContent(i, relUnc);
          H1D_RelativeUncertainty_SVD->SetBinError(i, 0.0);
      } else {
          H1D_RelativeUncertainty_SVD->SetBinContent(i, 0.0);
          H1D_RelativeUncertainty_SVD->SetBinError(i, 0.0);
      }
  }

  for (int i = 1; i <= H1D_Nominal_Bayes->GetNbinsX(); ++i) {
      H1D_Delta_Bayes->SetBinContent(i, abs(H1D_Nominal_Bayes->GetBinContent(i) - H1D_Reduced_Bayes->GetBinContent(i)));
      H1D_Delta_Bayes->SetBinError(i, 0.0);

      if (H1D_Nominal_Bayes->GetBinContent(i) != 0) {
          double relUnc = (H1D_Delta_Bayes->GetBinContent(i) / H1D_Nominal_Bayes->GetBinContent(i)) * 100.0; // in percent
          H1D_RelativeUncertainty_Bayes->SetBinContent(i, relUnc);
          H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
      } else {
          H1D_RelativeUncertainty_Bayes->SetBinContent(i, 0.0);
          H1D_RelativeUncertainty_Bayes->SetBinError(i, 0.0);
      }
  }
  TString textContextSVD("with 3% reduced track efficiency (SVD)");
  TString textContextBayes("with 3% reduced track efficiency (Bayes)");

  TString* pdfName_relUncSvd_logx = new TString("Relative_uncert_Svd_Data_train380686_Nominal_train533385_param7_ReducedTrackEff3perCent_train564527_param4_R=0.4_logx");
  TString* pdfName_relUncBayes_logx = new TString("Relative_uncert_Bayes_Data_train380686_Nominal_train533385_param4_ReducedTrackEff3perCent_train564527_param2_R=0.4_logx");
  TString* pdfName_relUncSvd = new TString("Relative_uncert_Svd_Data_train380686_Nominal_train533385_param7_ReducedTrackEff3perCent_train564527_param4_R=0.4");
  TString* pdfName_relUncBayes = new TString("Relative_uncert_Bayes_Data_train380686_Nominal_train533385_param4_ReducedTrackEff3perCent_train564527_param2_R=0.4");
  std::array<std::array<float, 2>, 2> drawnWindow = {{{5, 200}, {0, 80}}};

  Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd_logx, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes_logx, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "logx");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_SVD, textContextSVD, pdfName_relUncSvd, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  Draw_TH1_Histogram(H1D_RelativeUncertainty_Bayes, textContextBayes, pdfName_relUncBayes, texPtJetRec, texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, "");
  

}