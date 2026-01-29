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
#include <tuple>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();



void Get_systematics_UnfoldMethod(TH1D* &hSystematicUncertainty, TH1D* &hSystematicUncertainty_PreBarlow, int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_parameterVariation(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);
void Draw_Systematics_TrackEff(int iDataset, int iRadius, char** unfoldingMethodList, int* unfoldParameterInputList, int nUnfoldingMethods, std::string options);
void Draw_Systematics_SecondaryContamination(int iDataset, int iRadius, int unfoldParameterInput, std::string options);


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
  // char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  // snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);
  // const int nUnfoldingMethods = 4;
  // char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes", "Svd", "Bayes"}; // first two to be with nominal efficiency, last two with efficiency varied 
  // int unfoldParameterInputList[4] = {7, 4, 4, 2}; // first two to be with nominal efficiency, last two with efficiency varied
  // Draw_Systematics_TrackEff(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);

  //######################################################### Secondary tracks Systematics #####################################################
  char optionsAnalysis[100] = "";
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", unfoldingPrior, unfoldingMethod);
  int unfoldParameterInput = 7;
  Draw_Systematics_SecondaryContamination(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);



}

/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////


///////////// Fit functions //////////////////////////////////
// Convert a fitted function + fit result into a TGraphErrors that includes the 1σ confidence band of the fit.
TGraphErrors* getFunctionTGraphErrorsFromFitResult(double* xRangeFit, TF1* fitFunctionDrawn, TFitResultPtr fitResult, int nPointsGraph = 1000){
  std::vector<double> xAxisGraph= {};
  std::vector<double> yAxisGraph= {};
  std::vector<double> yAxisGraphErrors= {};
  // double* ;

  for(int iPoint = 0; iPoint < nPointsGraph; iPoint++){
    xAxisGraph.push_back(xRangeFit[0]+iPoint*1./nPointsGraph*(xRangeFit[1]-xRangeFit[0]));
    yAxisGraph.push_back(fitFunctionDrawn->Eval(xAxisGraph.back()));
    yAxisGraphErrors.push_back(0);
  }
  double oneSigmaInterval = 0.683;
  fitResult->GetConfidenceIntervals(nPointsGraph, 1, 1, &xAxisGraph[0], &yAxisGraphErrors[0], oneSigmaInterval, false);
  TGraphErrors* fitFunctionTGraphErrors = new TGraphErrors(nPointsGraph, &xAxisGraph[0], &yAxisGraph[0], nullptr, &yAxisGraphErrors[0]);
  return fitFunctionTGraphErrors;
}

// Fit a histogram with a double Tsallis-like function and return everything needed to propagate uncertainties
std::tuple<TF1*, TMatrixDSym, TFitResultPtr> FitDoubleTsallis(TH1D* &histogramInput, int nBinsX, double* binsX, double* xRangeFit) {
  TF1 *fitFunctionInit;
  TF1 *fitFunctionFinal;
  TF1 *fitFunctionDrawn; // drawn over the full range
  TFitResultPtr fFitResult;

  double parfitFunctionInit[4];
  double parfitFunctionFinal[4];
  // double parfitFunctionInit[8];
  // double parfitFunctionFinal[8];
  // const char* doubleTsallis = "([2]+[3]*x)*pow(1 + x/([0]*[1]), -[1]) + ([6]+[7]*x)*pow(1 + x/([4]*[5]), -[5])";
  const char* doubleTsallis = "([2]+[3]*x)*pow(1 + x/([0]*[1]), -[1])";


  ////////////////////////////////////////////////////////////////////
  //////////////////////////// Fit start /////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  fitFunctionInit = new TF1("fitFunctionInit_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  
  // Set parameter names
  fitFunctionInit->SetParName(0, "p0");
  fitFunctionInit->SetParName(1, "p1");
  fitFunctionInit->SetParName(2, "p2");
  fitFunctionInit->SetParName(3, "p3");
  // fitFunctionInit->SetParName(4, "p4");
  // fitFunctionInit->SetParName(5, "p5");
  // fitFunctionInit->SetParName(6, "p6");
  // fitFunctionInit->SetParName(7, "p7");

  // fitFunctionInit->SetParameters(0.5,  7,   50,  0,  1.2,  10,  300, 0);
  // //                             p0,   p1,   p2,  p3,  p4,  p5,  p6,  p7

  fitFunctionInit->SetParameters(0.5,  7,   0,  0);
  //                             p0,   p1,   p2,  p3

  fitFunctionInit->SetParLimits(0, 0.05, 1.0);
  fitFunctionInit->SetParLimits(1, 3.0, 10.0);
  fitFunctionInit->SetParLimits(2, -20.0, 2.0);
  fitFunctionInit->SetParLimits(3, -10.0, 70.0);
  // fitFunctionInit->SetParLimits(4, 0.05, 5.0);
  // fitFunctionInit->SetParLimits(5, 3.0, 30.0);
  // fitFunctionInit->SetParLimits(6, -50.0, 500.0);
  // fitFunctionInit->SetParLimits(7, -100.0, 100.0);

  histogramInput->Fit(fitFunctionInit, "R0Q"); // R = fit range, Q = quiet, L = likelihood
  fitFunctionInit->GetParameters(&parfitFunctionInit[0]); // Save initial parameters

  fitFunctionFinal = new TF1("fitFunctionFinal_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  
  for(int i=0; i<8; i++) fitFunctionFinal->SetParameter(i, parfitFunctionInit[i]);

  fFitResult = histogramInput->Fit(fitFunctionFinal, "RS");  
  fitFunctionFinal->GetParameters(&parfitFunctionFinal[0]);

  // Check covariance availability
  TMatrixDSym covMatrixFit; // default empty
  if (fFitResult && fFitResult->CovMatrixStatus() == 3) {
      covMatrixFit = fFitResult->GetCovarianceMatrix();
  } else {
      std::cout << "Warning: Covariance matrix not available!" << std::endl;
  }

  // TMatrixDSym covMatrixFit = fFitResult->GetCovarianceMatrix();

  // Double_t *pDataSmall = covMatrixFit.GetMatrixArray();
  // for (int i = 0; i < 2*2; i++) {
  //   cout << "i = " << i << ", covMatrixFit[i]" << pDataSmall[i] << endl;
  // }

  fitFunctionDrawn = new TF1("fitFunctionDrawn_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  for(int i=0; i<8; i++) fitFunctionDrawn->SetParameter(i, parfitFunctionFinal[i]);

  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> fitFunctionAndFitParams(fitFunctionDrawn, covMatrixFit, fFitResult);
  return fitFunctionAndFitParams;
}

// Use the fitted function to rebin a histogram and propagate fit uncertainties to the new bins
std::tuple<TH1D*, TGraphErrors*, TF1*> RebinWithDoubleTsallisFit(TH1D* &histogramInput, int nBinsX, double* binsX, double* xRangeFit) {
  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tsallisFitFunctionResult = FitDoubleTsallis(histogramInput, nBinsX, binsX, xRangeFit);
  TF1* fitFunctionDrawn = std::get<0>(tsallisFitFunctionResult);
  TFitResultPtr fitResult = std::get<2>(tsallisFitFunctionResult);
  TGraphErrors* fitFunctionTGraphErrors = getFunctionTGraphErrorsFromFitResult(xRangeFit, fitFunctionDrawn, fitResult);
  
  //////////////////////////// Rebin of input histogram /////////////////////////////

  TH1D* histogramRebinned = new TH1D("Unfolded: fit sampling", "Unfolded: fit sampling", nBinsX, binsX);
  for(int iBin = 0; iBin < nBinsX; iBin++){
    double xCenter = histogramRebinned->GetXaxis()->GetBinCenter(iBin);
    histogramRebinned->SetBinContent(iBin, fitFunctionDrawn->Eval(xCenter)); 
    double oneSigmaInterval = 0.683;
    double errorEval[1] = {0};
    double xEval[1] = {xCenter};
    fitResult->GetConfidenceIntervals(1, 1, 1, xEval, errorEval, oneSigmaInterval, false);
    histogramRebinned->SetBinError(iBin, errorEval[0]);
  }

  std::tuple<TH1D*, TGraphErrors*, TF1*> rebinResultAndFitFunction(histogramRebinned, fitFunctionTGraphErrors, fitFunctionDrawn);
  return rebinResultAndFitFunction;
}

//////////////////////////////////////////////////////////////

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

void Draw_Systematics_SecondaryContamination(int iDataset, int iRadius, int unfoldParameterInput, std::string options){
  cout << "########### Drawing systematics from secondary contamination variation ###############" << endl;
  TH1D* H1D_jetPt_unfolded;

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
  
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options); 

  // Define your bins and fit range
  int nBinsX = H1D_jetPt_unfolded->GetNbinsX();
  double* binsX = new double[nBinsX+1];
  for(int i=0; i<=nBinsX; i++) binsX[i] = H1D_jetPt_unfolded->GetBinLowEdge(i+1);

  double xRangeFit[2] = {5.0, 120.0}; // Fit range in GeV

  // Step 1: Rebin histogram using double Tsallis fit
  std::tuple<TH1D*, TGraphErrors*, TF1*> result = 
      RebinWithDoubleTsallisFit(H1D_jetPt_unfolded, nBinsX, binsX, xRangeFit);

  // Step 2: Extract outputs
  TH1D* hJetPtRebinned = std::get<0>(result);
  TGraphErrors* fitGraph = std::get<1>(result);
  TF1* fitFunctionDrawn = std::get<2>(result);

  // Step 3: Draw original histogram and rebinned fit
  TCanvas* c1 = new TCanvas("c1", "Double Tsallis Fit", 800, 600);
  H1D_jetPt_unfolded->SetMarkerStyle(20);
  H1D_jetPt_unfolded->SetMarkerColor(kBlack);
  H1D_jetPt_unfolded->Draw("E"); // original histogram with errors

  fitGraph->SetLineColor(kRed);
  fitGraph->SetLineWidth(2);
  fitGraph->Draw("L SAME"); // smooth fit with ±1σ band

  hJetPtRebinned->SetMarkerStyle(24);
  hJetPtRebinned->SetMarkerColor(kBlue);
  hJetPtRebinned->Draw("E SAME"); // rebinned histogram

  // c1->BuildLegend();
  // ================= Legend =================
  TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);   // transparent
  leg->SetTextSize(0.035);

  leg->AddEntry(H1D_jetPt_unfolded, "Unfolded data", "lep");
  leg->AddEntry(fitGraph, "Double Tsallis fit", "l");
  leg->AddEntry(hJetPtRebinned, "Fit sampling (rebinned)", "lep");

  leg->Draw();
  c1->Update();


  // TF1* ShiftTF1(TF1* f, double shift, const char* name="shifted") {
  //   return new TF1(name, [f, shift](double *x, double *){ return f->Eval(x[0] * shift); }, 
  //                  f->GetXmin(), f->GetXmax(), 0);
  // }

  // double shiftFactor = 0.005;
  // TF1* fUp   = ShiftTF1(fitFunctionDrawn, 1.0 + shiftFactor, "fUp");    // 1.005 * pT
  // TF1* fDown = ShiftTF1(fitFunctionDrawn, 1.0 - shiftFactor, "fDown");  // 0.995 * pT

  // TH1D* hRebinnedUp   = new TH1D("hRebinnedUp", "Rebinned Up (1.005x)", nBinsX, binsX);
  // TH1D* hRebinnedDown = new TH1D("hRebinnedDown", "Rebinned Down (0.995x)", nBinsX, binsX);

  // for(int iBin = 0; iBin < nBinsX; iBin++){
  //   double xCenter = hRebinnedUp->GetXaxis()->GetBinCenter(iBin);
  //   hRebinnedUp->SetBinContent(iBin, fUp->Eval(xCenter));     
  //   hRebinnedDown->SetBinContent(iBin, fDown->Eval(xCenter));
  // }

  // // --- Create histograms for absolute differences ---
  // TH1D* hDiffUp   = new TH1D("hDiffUp",   "Absolute difference Up",   nBinsX, binsX);
  // TH1D* hDiffDown = new TH1D("hDiffDown", "Absolute difference Down", nBinsX, binsX);

  // // --- Fill the difference histograms ---
  // for(int iBin = 0; iBin < nBinsX; iBin++){
  //     double nominal = hJetPtRebinned->GetBinContent(iBin);
  //     if(nominal == 0) nominal = 1e-12; // avoid division by zero

  //     double upVal   = hRebinnedUp->GetBinContent(iBin);
  //     double downVal = hRebinnedDown->GetBinContent(iBin);

  //     hRatioUp->SetBinContent(iBin,   fabs(upVal - nominal) / nominal * 100.0);
  //     hRatioDown->SetBinContent(iBin, fabs(downVal - nominal) / nominal * 100.0);
  // }

  // TH1D** deltaHistos = new TH1D*[2];
  // deltaHistos[0] = hRatioUp;    // Up variation
  // deltaHistos[1] = hRatioDown;  // Down variation

  // const TString Names[2] = {
  //   TString::Format("Up variation (+%.2f%%)", shiftFactor*100.0),
  //   TString::Format("Down variation (-%.2f%%)", shiftFactor*100.0)
  // };

  // TString pdfName = TString::Format("jet_spectrum_systematics_SecondaryTrackContamination_upDownVariation+%.3f.pdf", shiftFactor);
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{-25, 200},{0.85, 1.15}}};
  // std::array<std::array<float, 2>, 2> legendPlacement = {{{0.5, 0.7}, {0.75, 0.90}}}; // {{{x1, y1}, {x2, y2}}}
  // Draw_TH1_Histograms(deltaHistos, Names, 2, textContext, pdfNamePt, texPtJetRec , texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacementAuto, "");
}

/*
void Draw_Systematics_SecondaryContamination(int iDataset, int iRadius, int unfoldParameterInput, const double* xRangeFit, std::string options);


/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_systematics() {
  
  int iDataset = 0;
  int iRadius = 1;

  //######################################################### Secondary tracks Systematics #####################################################
  char optionsAnalysis[100] = "";
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s,%s", unfoldingPrior, unfoldingMethod);
  int unfoldParameterInput = 7;
  const double xRangeFit[2] = {5.0, 120.0}; // Fit range in GeV
  Draw_Systematics_SecondaryContamination(iDataset, iRadius, unfoldParameterInput, xRangeFit, optionsAnalysis);

}


///////////// Fit functions //////////////////////////////////
// Convert a fitted function + fit result into a TGraphErrors that includes the 1σ confidence band of the fit.
TGraphErrors* getFunctionTGraphErrorsFromFitResult(const double* xRangeFit, TF1* fitFunctionDrawn, TFitResultPtr fitResult, int nPointsGraph = 1000){
  std::vector<double> xAxisGraph= {};
  std::vector<double> yAxisGraph= {};
  std::vector<double> yAxisGraphErrors= {};
  // double* ;

  for(int iPoint = 0; iPoint < nPointsGraph; iPoint++){
    xAxisGraph.push_back(xRangeFit[0]+iPoint*1./nPointsGraph*(xRangeFit[1]-xRangeFit[0]));
    yAxisGraph.push_back(fitFunctionDrawn->Eval(xAxisGraph.back()));
    yAxisGraphErrors.push_back(0);
  }
  double oneSigmaInterval = 0.683;
  fitResult->GetConfidenceIntervals(nPointsGraph, 1, 1, &xAxisGraph[0], &yAxisGraphErrors[0], oneSigmaInterval, false);
  TGraphErrors* fitFunctionTGraphErrors = new TGraphErrors(nPointsGraph, &xAxisGraph[0], &yAxisGraph[0], nullptr, &yAxisGraphErrors[0]);
  return fitFunctionTGraphErrors;
}

// Fit a histogram with a double Tsallis-like function and return everything needed to propagate uncertainties
std::tuple<TF1*, TMatrixDSym, TFitResultPtr> FitDoubleTsallis(TH1D* &histogramInput, int nBinsX, double* binsX, const double* xRangeFit) {
  // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  // ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-4);
  const char* doubleTsallis = "([2]+[3]*x)*pow(1 + x/([0]*[1]), -[1]) + ([6]+[7]*x)*pow(1 + x/([4]*[5]), -[5])";

  TF1 *fitFunctionInit;
  TF1 *fitFunctionFinal;
  TF1 *fitFunctionDrawn; // drawn over the full range
  TFitResultPtr fFitResult;
  double parInit[8];
  double parFinal[8];

  // Initial Fit
  fitFunctionInit = new TF1("fitFunctionInit_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  // Set parameter names
  fitFunctionInit->SetParName(0, "p0");
  fitFunctionInit->SetParName(1, "p1");
  fitFunctionInit->SetParName(2, "p2");
  fitFunctionInit->SetParName(3, "p3");
  fitFunctionInit->SetParName(4, "p4");
  fitFunctionInit->SetParName(5, "p5");
  fitFunctionInit->SetParName(6, "p6");
  fitFunctionInit->SetParName(7, "p7");

  fitFunctionInit->SetParameters(0.5,  7,   -10,  50,  0.5,  5,  300, -70);
  //                             p0,   p1,   p2,  p3,  p4,  p5,  p6,  p7
  fitFunctionInit->SetParLimits(0, 0.05, 5.0);
  fitFunctionInit->SetParLimits(1, 3.0, 30.0);
  fitFunctionInit->SetParLimits(2, -20.0, 20.0);
  fitFunctionInit->SetParLimits(3, -10.0, 70.0);
  fitFunctionInit->SetParLimits(4, 0.05, 5.0);
  fitFunctionInit->SetParLimits(5, 3.0, 30.0);
  fitFunctionInit->SetParLimits(6, -50.0, 500.0);
  fitFunctionInit->SetParLimits(7, -100.0, 100.0);

  histogramInput->Fit(fitFunctionInit, "R0QL"); // R = fit range, Q = quiet, L = likelihood
  fitFunctionInit->GetParameters(&parInit[0]); // Save initial parameters

  // Final Fit
  fitFunctionFinal = new TF1("fitFunctionFinal_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  for(int i=0; i<8; i++) fitFunctionFinal->SetParameter(i, parInit[i]);
  
  fFitResult = histogramInput->Fit(fitFunctionFinal, "RSLH");  // S : return a TFitResultPtr (crucial for errors, covariance ...)
  fitFunctionFinal->GetParameters(&parFinal[0]);
  // Check covariance availability
  TMatrixDSym covMatrixFit; // default empty
  if (fFitResult && fFitResult->CovMatrixStatus() == 3) { // 0 : not calculated, 1 : approximated, 2 : forced pos. def., 3 : accurate
      covMatrixFit = fFitResult->GetCovarianceMatrix();
  } else {
      std::cout << "Warning: Covariance matrix not available!" << std::endl;
  }
  // TMatrixDSym covMatrixFit = fFitResult->GetCovarianceMatrix();
  // Double_t *pDataSmall = covMatrixFit.GetMatrixArray();
  // for (int i = 0; i < 2*2; i++) {
  //   cout << "i = " << i << ", covMatrixFit[i]" << pDataSmall[i] << endl;
  // }
  fitFunctionDrawn = new TF1("fitFunctionDrawn_", doubleTsallis, xRangeFit[0], xRangeFit[1]);
  for(int i=0; i<8; i++) fitFunctionDrawn->SetParameter(i, parFinal[i]);

  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> fitFunctionAndFitParams(fitFunctionDrawn, covMatrixFit, fFitResult);
  return fitFunctionAndFitParams;
}

// Use the fitted function to rebin a histogram and propagate fit uncertainties to the new bins
std::tuple<TH1D*, TGraphErrors*, TF1*> RebinWithDoubleTsallisFit(TH1D* &histogramInput, int nBinsX, double* binsX, const double* xRangeFit) {
  std::tuple<TF1*, TMatrixDSym, TFitResultPtr> tsallisFitFunctionResult = FitDoubleTsallis(histogramInput, nBinsX, binsX, xRangeFit);
  TF1* fitFunctionDrawn = std::get<0>(tsallisFitFunctionResult);
  TFitResultPtr fitResult = std::get<2>(tsallisFitFunctionResult);
  TGraphErrors* fitFunctionTGraphErrors = getFunctionTGraphErrorsFromFitResult(xRangeFit, fitFunctionDrawn, fitResult);
  
  //////////////////////////// Rebin of input histogram /////////////////////////////

  TH1D* histogramRebinned = new TH1D("Unfolded: fit sampling", "Unfolded: fit sampling", nBinsX, binsX);
  for(int iBin = 0; iBin < nBinsX; iBin++){
    double xCenter = histogramRebinned->GetXaxis()->GetBinCenter(iBin);
    histogramRebinned->SetBinContent(iBin, fitFunctionDrawn->Eval(xCenter)); 
    double oneSigmaInterval = 0.683;
    double errorEval[1] = {0};
    double xEval[1] = {xCenter};
    fitResult->GetConfidenceIntervals(1, 1, 1, xEval, errorEval, oneSigmaInterval, false); 
    histogramRebinned->SetBinError(iBin, errorEval[0]);
  }

  std::tuple<TH1D*, TGraphErrors*, TF1*> rebinResultAndFitFunction(histogramRebinned, fitFunctionTGraphErrors, fitFunctionDrawn);
  return rebinResultAndFitFunction;
}

TF1* ShiftTF1(const TF1* f, double shift, const char* name="shifted") {
    return new TF1(name, [f, shift](double *x, double *){ return f->Eval(x[0] * shift); }, 
                   f->GetXmin(), f->GetXmax(), 0);
  }
/////////////////////////////////////////////////////
void Draw_Systematics_SecondaryContamination(int iDataset, int iRadius, int unfoldParameterInput, const double* xRangeFit, std::string options){
  cout << "########### Drawing systematics from secondary contamination variation ###############" << endl;
  TH1D* H1D_jetPt_unfolded;

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
  
  Get_Pt_spectrum_unfolded(H1D_jetPt_unfolded, measuredInput, iDataset, iRadius, unfoldParameterInput, options); 

  // Define your bins and fit range
  int nBinsX = H1D_jetPt_unfolded->GetNbinsX();
  double* binsX = new double[nBinsX+1];
  for(int i=0; i<=nBinsX; i++) binsX[i] = H1D_jetPt_unfolded->GetBinLowEdge(i+1);

  // double xRangeFit[2] = {5.0, 120.0}; // Fit range in GeV
  cout << "########### Fit range: [" << xRangeFit[0] << " , " << xRangeFit[1] << "] GeV #############" << endl;
  // Step 1: Rebin histogram using double Tsallis fit
  std::tuple<TH1D*, TGraphErrors*, TF1*> result = RebinWithDoubleTsallisFit(H1D_jetPt_unfolded, nBinsX, binsX, xRangeFit);

  // Step 2: Extract outputs
  TH1D* hJetPtRebinned = std::get<0>(result);
  TGraphErrors* fitGraph = std::get<1>(result);
  TF1* fitFunctionDrawn = std::get<2>(result);

  // Step 3: Draw original histogram and rebinned fit
  TCanvas* c1 = new TCanvas("c1", "Double Tsallis Fit", 800, 600);
  H1D_jetPt_unfolded->SetMarkerStyle(20);
  H1D_jetPt_unfolded->SetMarkerColor(kBlack);
  H1D_jetPt_unfolded->Draw("E"); // original histogram with errors

  fitGraph->SetLineColor(kRed);
  fitGraph->SetLineWidth(2);
  fitGraph->Draw("L SAME"); // smooth fit with ±1σ band

  hJetPtRebinned->SetMarkerStyle(24);
  hJetPtRebinned->SetMarkerColor(kBlue);
  hJetPtRebinned->Draw("E SAME"); // rebinned histogram

  // c1->BuildLegend();
  // ================= Legend =================
  TLegend* leg = new TLegend(0.55, 0.65, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);   // transparent
  leg->SetTextSize(0.035);

  leg->AddEntry(H1D_jetPt_unfolded, "Unfolded data", "lep");
  leg->AddEntry(fitGraph, "Double Tsallis fit", "l");
  leg->AddEntry(hJetPtRebinned, "Fit sampling (rebinned)", "lep");

  leg->Draw();
  c1->Update();

  // double shiftFactor = 0.005;
  // TF1* fUp   = ShiftTF1(fitFunctionDrawn, 1.0 + shiftFactor, "fUp");    // 1.005 * pT
  // TF1* fDown = ShiftTF1(fitFunctionDrawn, 1.0 - shiftFactor, "fDown");  // 0.995 * pT

  // TString titleUp   = Form("Rebinned Up (+%.2f%%)",   shiftFactor * 100.0);
  // TString titleDown = Form("Rebinned Down (-%.2f%%)", shiftFactor * 100.0);

  // TH1D* hRebinnedUp   = new TH1D("hRebinnedUp",   titleUp,   nBinsX, binsX);
  // TH1D* hRebinnedDown = new TH1D("hRebinnedDown", titleDown, nBinsX, binsX);

  // for(int iBin = 1; iBin < nBinsX; iBin++){
  //   double xCenter = hRebinnedUp->GetXaxis()->GetBinCenter(iBin);
  //   hRebinnedUp->SetBinContent(iBin, fUp->Eval(xCenter));     
  //   hRebinnedDown->SetBinContent(iBin, fDown->Eval(xCenter));
  // }

  // // --- Create histograms for absolute differences ---
  // TH1D* hRatioUp   = new TH1D("hRatioUp",   "Ratio difference Up",   nBinsX, binsX);
  // TH1D* hRatioDown = new TH1D("hRatioDown", "Ratio difference Down", nBinsX, binsX);

  // // --- Fill the difference histograms ---
  // for(int iBin = 0; iBin < nBinsX; iBin++){
  //     double nominal = hJetPtRebinned->GetBinContent(iBin);
  //     if(nominal == 0) nominal = 1e-12; // avoid division by zero

  //     double upVal   = hRebinnedUp->GetBinContent(iBin);
  //     double downVal = hRebinnedDown->GetBinContent(iBin);

  //     hRatioUp->SetBinContent(iBin,   fabs(upVal - nominal) / nominal * 100.0);
  //     hRatioDown->SetBinContent(iBin, fabs(downVal - nominal) / nominal * 100.0);
  // }

  // TH1D** deltaHistos = new TH1D*[2];
  // deltaHistos[0] = hRatioUp;    // Up variation
  // deltaHistos[1] = hRatioDown;  // Down variation

  // const TString Names[2] = {
  //   TString::Format("Up variation (+%.2f%%)", shiftFactor*100.0),
  //   TString::Format("Down variation (-%.2f%%)", shiftFactor*100.0)
  // };

  // TString pdfName = TString::Format("jet_spectrum_systematics_SecondaryTrackContamination_upDownVariation+%.3f.pdf", shiftFactor);
  // std::array<std::array<float, 2>, 2> drawnWindow = {{{-25, 200},{0.85, 1.15}}};
  // std::array<std::array<float, 2>, 2> legendPlacement = {{{0.5, 0.7}, {0.75, 0.90}}}; // {{{x1, y1}, {x2, y2}}}
  // TString textContext = "";
  // TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV");
  // Draw_TH1_Histograms(deltaHistos, Names, 2, textContext, pdfName, texPtJetRec , texSystematicsPercent, texCollisionDataInfo, drawnWindow, legendPlacement, contextPlacementAuto, "");
}
  */