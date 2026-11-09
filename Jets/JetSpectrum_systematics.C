#if defined(__CLING__)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "TStyle.h"
#include "TFile.h"
#include "TDirectoryFile.h"
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
#include "./JetSpectrum_Systematics_input.h"

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
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle_Systematics(Bool_t graypalette=kFALSE);
void LoadLibs_Systematics();



void Get_systematics_UnfoldMethod(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, const char* const* methods, const int* parameters, int nMethods, std::string options, TH1D** reference = nullptr, bool drawSpectra = true);
void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, const char* const* methods, const int* parameters, int nMethods, std::string options, TH1D** relativeResult = nullptr);
void Get_systematics_SvdParameter(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, int nominalK, std::string options, TH1D** reference = nullptr, bool drawSpectra = true);
void Draw_Systematics_SvdParameter(int iDataset, int iRadius, int nominalK, std::string options, TH1D** relativeResult = nullptr);
void Get_systematics_BackgroundFluctuations(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, int unfoldParameter, std::string options, TH1D** reference = nullptr);
void Draw_Systematics_BackgroundFluctuations(int iDataset, int iRadius, int unfoldParameter, std::string options);
void Get_systematics_TrackingEfficiency(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, int unfoldParameter, std::string options, TH1D** reference = nullptr);
void Draw_Systematics_TrackingEfficiency(int iDataset, int iRadius, int unfoldParameter, std::string options);
void Get_systematics_PtRange(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius, int unfoldParameter, std::string options, bool varyMeasured, TH1D** reference = nullptr, bool drawSpectra = true);
void Draw_Systematics_PtRange(int iDataset, int iRadius, int unfoldParameter, std::string options, bool varyMeasured, TH1D** relativeResult = nullptr);
TH1D* DrawUnfoldingSys(int iDataset = 0, int iRadius = 0, int nominalSvdK = 11, int bayesIterations = 7);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_systematics() {
  cout << "Apply Barlow test: " << (applyBarlowTest ? "true" : "false") << endl;
  // Load necessary libraries
  LoadLibs_Systematics();
  // Set the default style
  SetStyle_Systematics();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");

  // gathers the analysis options in a single char[]

  int iDataset = 0;
  int iRadius = 0;

  char optionsAnalysis_withoutUnfoldingMethod[100] = "";
  snprintf(optionsAnalysis_withoutUnfoldingMethod, sizeof(optionsAnalysis_withoutUnfoldingMethod), "%s", unfoldingPrior);


  const int nUnfoldingMethods = 2;
  const char* unfoldingMethodList[nUnfoldingMethods] = {"Svd", "Bayes"}; // default is the first one in this list
  int unfoldParameterInputList[2] = {11, 7};

  // Draw_Systematics_UnfoldMethod(iDataset, iRadius, unfoldingMethodList, unfoldParameterInputList, nUnfoldingMethods, optionsAnalysis_withoutUnfoldingMethod);
  // Draw_Systematics_SvdParameter(iDataset, iRadius, unfoldParameterInputList[0], optionsAnalysis_withoutUnfoldingMethod);
  // Draw_Systematics_BackgroundFluctuations(iDataset, iRadius, unfoldParameterInputList[0], std::string(optionsAnalysis_withoutUnfoldingMethod) + "," + unfoldingMethodList[0]);
  // Draw_Systematics_TrackingEfficiency(iDataset, iRadius, unfoldParameterInputList[0], std::string(optionsAnalysis_withoutUnfoldingMethod) + "," + unfoldingMethodList[0]);
  // Draw_Systematics_PtRange(iDataset, iRadius, unfoldParameterInputList[0], std::string(optionsAnalysis_withoutUnfoldingMethod) + "," + unfoldingMethodList[0], true);
  // Draw_Systematics_PtRange(iDataset, iRadius, unfoldParameterInputList[0], std::string(optionsAnalysis_withoutUnfoldingMethod) + "," + unfoldingMethodList[0], false);
  // Calculate all four unfolding sources and draw the combined result.
  DrawUnfoldingSys(iDataset, iRadius, unfoldParameterInputList[0], unfoldParameterInputList[1]);
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

// First spectrum is nominal. Use the largest absolute variation in each bin.
// When Barlow is enabled, test each variation before taking the accepted envelope.
void Get_Systematics_Envelope(TH1D* &uncertainty, TH1D* &preBarlow, const std::vector<TH1D*>& spectra, TString name) {
  uncertainty = (TH1D*)spectra[0]->Clone("hSystematicUncertainty_" + name);
  uncertainty->Reset("ICES");
  preBarlow = nullptr;
  if (applyBarlowTest) {
    preBarlow = (TH1D*)uncertainty->Clone("hSystematicUncertainty_PreBarlow_" + name);
  }
   
  for (int bin = 1; bin <= spectra[0]->GetNbinsX(); ++bin) {
    double maximum = 0., acceptedMaximum = 0.;
    double maximumSigma = 0., acceptedSigma = 0.;
    for (size_t variation = 1; variation < spectra.size(); ++variation) {
      double difference = std::abs(spectra[variation]->GetBinContent(bin)
                                   - spectra[0]->GetBinContent(bin));
      if (!std::isfinite(difference)) throw std::runtime_error("Non-finite systematic variation");
      if (!applyBarlowTest) {
        maximum = std::max(maximum, difference);
        continue; // No Barlow error calculation or threshold test.
      }
      double error = spectra[variation]->GetBinError(bin);
      double referenceError = spectra[0]->GetBinError(bin);
      double sigma = std::sqrt(std::abs(error*error - referenceError*referenceError));
      if (!std::isfinite(sigma)) throw std::runtime_error("Non-finite Barlow uncertainty");
      if (difference > maximum) {
        maximum = difference;
        maximumSigma = sigma;
      }
      if (difference > N_SigmaBarlow*sigma && difference > acceptedMaximum) {
        acceptedMaximum = difference;
        acceptedSigma = sigma;
      }
    }
    uncertainty->SetBinContent(bin, applyBarlowTest ? acceptedMaximum : maximum);
    uncertainty->SetBinError(bin, applyBarlowTest ? acceptedSigma : 0.);
    if (applyBarlowTest) {
      preBarlow->SetBinContent(bin, maximum);
      preBarlow->SetBinError(bin, maximumSigma);
    }
  }
}

TH1D* Get_Systematics_MeasuredInput(int iDataset, int iRadius, std::string options) {
  TH1D* measured = nullptr;
  if (!normGenAndMeasByNEvtsForUnfoldingInput) {
    if (useFineBinningTest)
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(measured, iDataset, iRadius, options);
    else
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(measured, iDataset, iRadius, options);
  } else {
    if (useFineBinningTest)
      Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(measured, iDataset, iRadius, options);
    else
      Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(measured, iDataset, iRadius, options);
  }
  return measured;
}

void Draw_Systematics_Spectra(std::vector<TH1D*>& spectra, const std::vector<TString>& legends,
                             TString name, int iRadius) {
  TString* pdfName = new TString(name);
  TString* yTitle = normaliseUnfoldingResultsAtEnd ? texJet_d2Ndptdeta_EventNorm : new TString("d^{2}#it{N}_{jet}/d#it{p}_{T}d#it{#eta} (GeV/#it{c})^{-1}");
  std::array<std::array<float, 2>, 2> legendPlacement = {{{0.53, 0.59}, {0.94, 0.79}}};
  if (spectra.size() > 3) legendPlacement = {{{0.53, 0.47}, {0.94, 0.79}}};
  Draw_TH1_Histograms(spectra.data(), legends.data(), spectra.size(), contextJetRadius(arrayRadius[iRadius]), pdfName, texPtJetGen, yTitle, texCollisionDataInfo, drawnWindowUnfoldedMeasurement, legendPlacement, contextPlacementAuto, "logy");
}

// Convert the absolute envelope to a fraction of the nominal spectrum.
// No independent ratio-error propagation: nominal and variations share the data.
void Draw_Systematics_RelativeUncertainty(TH1D* uncertainty, TH1D* preBarlow,
                                         TH1D* reference, TString name, int iRadius) {
  auto relativeWindow = drawnWindowUnfoldedMeasurement;
  relativeWindow[1] = {relativeSystematicsRange[0], relativeSystematicsRange[1]};
  TH1D* histograms[2] = {uncertainty, preBarlow};
  int nHistograms = applyBarlowTest ? 2 : 1;
  for (int h = 0; h < nHistograms; ++h) {
    for (int bin = 1; bin <= reference->GetNbinsX(); ++bin) {
      double denominator = std::abs(reference->GetBinContent(bin));
      if (!std::isfinite(denominator)) throw std::runtime_error("Non-finite nominal spectrum");
      if (denominator == 0.) {
        if (histograms[h]->GetBinContent(bin) != 0.)
          throw std::runtime_error("Cannot form relative uncertainty with zero nominal yield");
        histograms[h]->SetBinContent(bin, 0.);
        histograms[h]->SetBinError(bin, 0.);
      } else {
        histograms[h]->SetBinContent(bin, histograms[h]->GetBinContent(bin)/denominator);
        histograms[h]->SetBinError(bin, histograms[h]->GetBinError(bin)/denominator);
      }
    }
    TString suffix = !applyBarlowTest ? "_NoBarlow" : (h == 0 ? "_Barlow" : "_PreBarlow");
    TString* pdfName = new TString(name + suffix);
    TString* yTitle = new TString("Relative systematic uncertainty");
    Draw_TH1_Histogram(histograms[h], contextJetRadius(arrayRadius[iRadius]), pdfName,
                       texPtJetGen, yTitle, texCollisionDataInfo, relativeWindow,
                       legendPlacementAuto, contextPlacementAuto, "");
    cout << name << suffix << ": relative maximum = " << histograms[h]->GetMaximum() << endl;
  }
}

void Get_systematics_UnfoldMethod(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset, int iRadius,
                                  const char* const* methods, const int* parameters, int nMethods,
                                  std::string options, TH1D** reference, bool drawSpectra) {
  if (nMethods < 2) throw std::runtime_error("At least two unfolding methods are required");
  TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
  std::vector<TH1D*> spectra(nMethods);
  std::vector<TString> legends(nMethods);
  TString name = "UnfoldMethod_" + Datasets[iDataset] + Form("_R=%.1f", arrayRadius[iRadius]);
  for (int method = 0; method < nMethods; ++method) {
    Get_Pt_spectrum_unfolded(spectra[method], measured, iDataset, iRadius, parameters[method],
                             options + "," + methods[method]);
    legends[method] = TString(methods[method]) == "Svd"
      ? Form("SVD, k = %d", parameters[method])
      : Form("%s, iter. = %d", methods[method], parameters[method]);
    name += Form("_%s%d", methods[method], parameters[method]);
  }
  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  if (drawSpectra) Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_UnfoldMethod(int iDataset, int iRadius, const char* const* methods,
                                   const int* parameters, int nMethods, std::string options, TH1D** relativeResult) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_UnfoldMethod(uncertainty, preBarlow, iDataset, iRadius, methods, parameters,
                               nMethods, options, &reference);
  TString name = "Systematics_UnfoldMethod_" + Datasets[iDataset]
                 + Form("_R=%.1f", arrayRadius[iRadius]);
  for (int method = 0; method < nMethods; ++method)
    name += Form("_%s%d", methods[method], parameters[method]);
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
  if (relativeResult) *relativeResult = uncertainty;
}

void Get_systematics_SvdParameter(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset,
                                 int iRadius, int nominalK, std::string options, TH1D** reference, bool drawSpectra) {
  if (automaticBestSvdParameter)
    throw std::runtime_error("Disable automaticBestSvdParameter for an explicit k variation");
  int nBins = useFineBinningTest ? nBinPtJetsFine[iRadius]
                               : std::min(nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius]);
  if (nominalK < 2 || nominalK + 1 > nBins)
    throw std::runtime_error("SVD k-1 and k+1 must lie within the response bin count");
  TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
  const int parameters[3] = {nominalK, nominalK - 1, nominalK + 1};
  std::vector<TH1D*> spectra(3);
  std::vector<TString> legends(3);
  TString name = "SvdParameter_" + Datasets[iDataset]
                 + Form("_R=%.1f_k%d", arrayRadius[iRadius], nominalK);
  for (int variation = 0; variation < 3; ++variation) {
    Get_Pt_spectrum_unfolded(spectra[variation], measured, iDataset, iRadius, parameters[variation],
                             options + ",Svd");
    legends[variation] = Form("SVD, k = %d%s", parameters[variation], variation == 0 ? " (nominal)" : "");
  }
  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  if (drawSpectra) Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_SvdParameter(int iDataset, int iRadius, int nominalK, std::string options, TH1D** relativeResult) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_SvdParameter(uncertainty, preBarlow, iDataset, iRadius, nominalK, options, &reference);
  TString name = "Systematics_SvdParameter_" + Datasets[iDataset]
                 + Form("_R=%.1f_k%d", arrayRadius[iRadius], nominalK);
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
  if (relativeResult) *relativeResult = uncertainty;
}

void Get_systematics_BackgroundFluctuations(TH1D* &uncertainty, TH1D* &preBarlow,
                                          int iDataset, int iRadius, int unfoldParameter,
                                          std::string options, TH1D** reference) {
  if (!useFactorisedMatrix || options.find("useIdentityForFluctResp") != std::string::npos)
    throw std::runtime_error("Background systematics require the fluctuation response to be enabled");
  if (options.find("Svd") != std::string::npos && automaticBestSvdParameter)
    throw std::runtime_error("Use a fixed SVD parameter for background variations");
  if (nBackgroundVariations < 1)
    throw std::runtime_error("At least one background variation is required");

  TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
  std::vector<TH1D*> spectra(nBackgroundVariations + 1);
  std::vector<TString> legends(nBackgroundVariations + 1);
  TString method = options.find("Svd") != std::string::npos ? "Svd" : "Bayes";
  TString name = "BackgroundFluctuations_" + Datasets[iDataset]
                 + Form("_R=%.1f_%s%d", arrayRadius[iRadius], method.Data(), unfoldParameter);

  // Nominal background remains the file/workflow/randomConeType in nominal inputs.
  Get_Pt_spectrum_unfolded(spectra[0], measured, iDataset, iRadius, unfoldParameter, options);
  legends[0] = "Nominal RC";
  for (int variation = 0; variation < nBackgroundVariations; ++variation) {
    TH2D* fluctuations = nullptr;
    Get_PtResponseMatrix_Fluctuations(fluctuations,
      backgroundVariationRandomConeTypes[variation], iDataset, iRadius, "",
      file_O2Analysis_Bkg_Systematics_list[iDataset], backgroundVariationWorkflows[variation]);
    fluctuations->SetName(Form("BackgroundVariation%d", variation + 1));
    fluctuations->SetDirectory(nullptr);
    cout << "Background variation " << variation + 1 << ": "
         << file_O2Analysis_Bkg_Systematics_list[iDataset]->GetName() << ", "
         << backgroundVariationWorkflows[variation] << ", randomConeType="
         << backgroundVariationRandomConeTypes[variation] << endl;
    // Explicit matrix input: no changes to nominal file pointers or settings.
    Get_Pt_spectrum_unfolded(spectra[variation + 1], measured, iDataset, iRadius,
                             unfoldParameter, options, false, false, fluctuations);
    delete fluctuations;
    legends[variation + 1] = backgroundVariationLabels[variation];
  }
  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_BackgroundFluctuations(int iDataset, int iRadius, int unfoldParameter,
                                           std::string options) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_BackgroundFluctuations(uncertainty, preBarlow, iDataset, iRadius,
                                        unfoldParameter, options, &reference);
  TString method = options.find("Svd") != std::string::npos ? "Svd" : "Bayes";
  TString name = "Systematics_BackgroundFluctuations_" + Datasets[iDataset]
                 + Form("_R=%.1f_%s%d", arrayRadius[iRadius], method.Data(), unfoldParameter);
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
}

// All MC-dependent getters read this file pointer: response, prior, fakes and efficiency.
// Restore nominal input on both normal return and a C++ exception.
void Get_Systematics_UnfoldedWithMC(TH1D* &spectrum, TH1D* &measured, int iDataset,
                                   int iRadius, int unfoldParameter, std::string options,
                                   TFile* mcFile, TString mcWorkflow) {
  if (!mcFile || mcFile->IsZombie())
    throw std::runtime_error("Cannot open tracking-efficiency MC file");
  if (mcWorkflow != analysisWorkflowMC)
    throw std::runtime_error("MC variation must use the nominal workflow name");

  TString matched = "h2_jet_pt_mcd_jet_pt_mcp_matchedgeo";
  matched += doBkgSubtractionInMC ? "_rhoareasubtracted"
    : (etaCutOnMatchedJetsIsObsoleteVersion ? "_mcpetaconstraint" : "_mcpmcdetaconstraint");
  const TString required[] = {
    doBkgSubtractionInMC ? "h_jet_pt_rhoareasubtracted" : "h_jet_pt",
    doBkgSubtractionInMC ? "h_jet_pt_part_rhoareasubtracted" : "h_jet_pt_part",
    matched,
    mcIsWeighted ? "h_collisions_weighted" : "h_collisions",
    mcCollHistIsObsolete
      ? (mcIsWeighted ? "h_mcColl_counts_weight" : "h_mcColl_counts")
      : (mcIsWeighted ? "h_mccollisions_weighted" : "h_mccollisions")
  };
  for (const TString& histogram : required) {
    if (!dynamic_cast<TH1*>(mcFile->Get(mcWorkflow + "/" + histogram)))
      throw std::runtime_error(std::string("Missing MC histogram: ")
                               + mcFile->GetName() + ":" + mcWorkflow.Data() + "/" + histogram.Data());
  }

  TFile* nominalMC = file_O2Analysis_MCfile_GeneralResponse[iDataset];
  file_O2Analysis_MCfile_GeneralResponse[iDataset] = mcFile;
  try {
    Get_Pt_spectrum_unfolded(spectrum, measured, iDataset, iRadius, unfoldParameter,
                             options, false, false);
  } catch (...) {
    file_O2Analysis_MCfile_GeneralResponse[iDataset] = nominalMC;
    throw;
  }
  file_O2Analysis_MCfile_GeneralResponse[iDataset] = nominalMC;
}

void Get_systematics_TrackingEfficiency(TH1D* &uncertainty, TH1D* &preBarlow,
                                       int iDataset, int iRadius, int unfoldParameter,
                                       std::string options, TH1D** reference) {
  if (!file_O2Analysis_MCfile_TrackingEfficiency_Systematics[iDataset])
    throw std::runtime_error("Tracking-efficiency MC input is not configured");
  if (options.find("Svd") != std::string::npos && automaticBestSvdParameter)
    throw std::runtime_error("Use a fixed SVD parameter for tracking-efficiency variations");
  TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
  std::vector<TH1D*> spectra(2);
  std::vector<TString> legends = {"Nominal MC", "Tracking eff. varied"};
  TString method = options.find("Svd") != std::string::npos ? "Svd" : "Bayes";
  TString name = "TrackingEfficiency_" + Datasets[iDataset]
                 + Form("_R=%.1f_%s%d", arrayRadius[iRadius], method.Data(), unfoldParameter);

  Get_Pt_spectrum_unfolded(spectra[0], measured, iDataset, iRadius, unfoldParameter, options);
  spectra[0]->SetName("Nominal_" + name);
  spectra[0]->SetDirectory(nullptr);
  cout << "Tracking efficiency variation: "
       << file_O2Analysis_MCfile_TrackingEfficiency_Systematics[iDataset]->GetName()
       << ", " << analysisWorkflowMC_TrackingEfficiency_Systematics << endl;
  Get_Systematics_UnfoldedWithMC(spectra[1], measured, iDataset, iRadius, unfoldParameter, options,
    file_O2Analysis_MCfile_TrackingEfficiency_Systematics[iDataset],
    analysisWorkflowMC_TrackingEfficiency_Systematics);
  spectra[1]->SetName("Variation_" + name);
  spectra[1]->SetDirectory(nullptr);

  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_TrackingEfficiency(int iDataset, int iRadius, int unfoldParameter,
                                       std::string options) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_TrackingEfficiency(uncertainty, preBarlow, iDataset, iRadius,
                                    unfoldParameter, options, &reference);
  TString method = options.find("Svd") != std::string::npos ? "Svd" : "Bayes";
  TString name = "Systematics_TrackingEfficiency_" + Datasets[iDataset]
                 + Form("_R=%.1f_%s%d", arrayRadius[iRadius], method.Data(), unfoldParameter);
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
}

// Keep all nominal interior edges. Insert only the requested new endpoints.
std::vector<double> Get_Systematics_PtRangeEdges(const double* nominal, int nBins,
                                               int capacity, double minimum, double maximum,
                                               int iRadius) {
  if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum >= maximum)
    throw std::runtime_error("Invalid pT-range variation");
  for (double boundary : {minimum, maximum}) {
    bool onFineEdge = false;
    for (int edge = 0; edge <= nBinPtJetsFine[iRadius]; ++edge)
      if (std::abs(boundary - ptBinsJetsFine[iRadius][edge]) < 1.e-8) onFineEdge = true;
    if (!onFineEdge)
      throw std::runtime_error("pT-range endpoints must be edges of the fine response binning");
  }
  std::vector<double> edges = {minimum};
  for (int edge = 0; edge <= nBins; ++edge) {
    if (edge && nominal[edge] <= nominal[edge-1])
      throw std::runtime_error("Nominal pT bin edges must be strictly increasing");
    if (nominal[edge] > minimum && nominal[edge] < maximum) edges.push_back(nominal[edge]);
  }
  edges.push_back(maximum);
  if (int(edges.size()) > capacity)
    throw std::runtime_error("pT-range variation exceeds the bin-array capacity");
  return edges;
}

// Change only one row and its count for the selected radius, then restore both rows.
// Re-read measured data and rebuild response/prior/efficiencies for every variation.
void Get_Systematics_UnfoldedWithPtRange(TH1D* &spectrum, int iDataset, int iRadius,
                                        int unfoldParameter, std::string options,
                                        bool varyMeasured, double minimum, double maximum) {
  if (useFineBinningTest)
    throw std::runtime_error("Disable useFineBinningTest for rec/gen pT-range variations");
  if (options.find("Svd") != std::string::npos && automaticBestSvdParameter)
    throw std::runtime_error("Use a fixed SVD parameter for pT-range variations");
  const int recCapacity = sizeof(ptBinsJetsRec[iRadius])/sizeof(double);
  const int genCapacity = sizeof(ptBinsJetsGen[iRadius])/sizeof(double);
  const int nominalRecCount = nBinPtJetsRec[iRadius], nominalGenCount = nBinPtJetsGen[iRadius];
  std::vector<double> nominalRec(ptBinsJetsRec[iRadius], ptBinsJetsRec[iRadius] + recCapacity);
  std::vector<double> nominalGen(ptBinsJetsGen[iRadius], ptBinsJetsGen[iRadius] + genCapacity);
  std::vector<double> edges = Get_Systematics_PtRangeEdges(
    varyMeasured ? nominalRec.data() : nominalGen.data(),
    varyMeasured ? nominalRecCount : nominalGenCount,
    varyMeasured ? recCapacity : genCapacity, minimum, maximum, iRadius);
  auto restore = [&]() {
    std::copy(nominalRec.begin(), nominalRec.end(), ptBinsJetsRec[iRadius]);
    std::copy(nominalGen.begin(), nominalGen.end(), ptBinsJetsGen[iRadius]);
    nBinPtJetsRec[iRadius] = nominalRecCount;
    nBinPtJetsGen[iRadius] = nominalGenCount;
  };
  double* selectedEdges = varyMeasured ? ptBinsJetsRec[iRadius] : ptBinsJetsGen[iRadius];
  std::copy(edges.begin(), edges.end(), selectedEdges);
  (varyMeasured ? nBinPtJetsRec[iRadius] : nBinPtJetsGen[iRadius]) = edges.size() - 1;
  try {
    if (options.find("Svd") != std::string::npos &&
        (unfoldParameter < 1 || unfoldParameter > std::min(nBinPtJetsRec[iRadius], nBinPtJetsGen[iRadius])))
      throw std::runtime_error("SVD parameter is outside the varied response bin count");
    TH1D* measured = Get_Systematics_MeasuredInput(iDataset, iRadius, options);
    cout << (varyMeasured ? "Measured" : "Unfolded") << " pT range: "
         << minimum << "-" << maximum << ", rec bins=" << nBinPtJetsRec[iRadius]
         << ", gen bins=" << nBinPtJetsGen[iRadius] << endl;
    Get_Pt_spectrum_unfolded(spectrum, measured, iDataset, iRadius, unfoldParameter, options);
    spectrum->SetName(Form("%sPtRange_%g_%g_%s_R%g", varyMeasured ? "Measured" : "Unfolded",
                           minimum, maximum, Datasets[iDataset].Data(), arrayRadius[iRadius]));
    spectrum->SetDirectory(nullptr);
  } catch (...) {
    restore();
    throw;
  }
  restore();
}

// Compare identical physical bins in the final display range only: do not divide
// by bin index or interpolate a truncated density bin (140-180 versus 140-200).
void Align_Systematics_CommonPtBins(std::vector<TH1D*>& spectra, TString name) {
  std::vector<int> nominalBins;
  std::vector<double> edges;
  for (int bin = 1; bin <= spectra[0]->GetNbinsX(); ++bin) {
    double low = spectra[0]->GetXaxis()->GetBinLowEdge(bin);
    double high = spectra[0]->GetXaxis()->GetBinUpEdge(bin);
    if (high <= ptWindowDisplay[0] + 1.e-8 || low >= ptWindowDisplay[1] - 1.e-8) continue;
    bool common = true;
    for (TH1D* spectrum : spectra) {
      int matched = spectrum->GetXaxis()->FindFixBin((low + high)/2.);
      if (matched < 1 || matched > spectrum->GetNbinsX() ||
          std::abs(spectrum->GetXaxis()->GetBinLowEdge(matched) - low) > 1.e-8 ||
          std::abs(spectrum->GetXaxis()->GetBinUpEdge(matched) - high) > 1.e-8) common = false;
    }
    if (!common) continue;
    if (!nominalBins.empty() && bin != nominalBins.back() + 1)
      throw std::runtime_error("Common pT comparison bins are not contiguous");
    if (edges.empty()) edges.push_back(low);
    nominalBins.push_back(bin);
    edges.push_back(high);
  }
  if (nominalBins.empty()) throw std::runtime_error("No common physical pT bins");
  if (edges.front() > ptWindowDisplay[0] + 1.e-8 || edges.back() < ptWindowDisplay[1] - 1.e-8)
    throw std::runtime_error("Common pT bins do not cover the final display range");
  for (size_t variation = 0; variation < spectra.size(); ++variation) {
    TH1D* aligned = new TH1D(Form("CommonBins_%s_%zu", name.Data(), variation), "",
                             nominalBins.size(), edges.data());
    aligned->SetDirectory(nullptr);
    for (size_t bin = 0; bin < nominalBins.size(); ++bin) {
      int matched = spectra[variation]->GetXaxis()->FindFixBin((edges[bin] + edges[bin+1])/2.);
      aligned->SetBinContent(bin + 1, spectra[variation]->GetBinContent(matched));
      aligned->SetBinError(bin + 1, spectra[variation]->GetBinError(matched));
    }
    spectra[variation] = aligned;
  }
  cout << name << ": common comparison bins " << edges.front() << "-" << edges.back()
       << ", " << nominalBins.size() << " bins" << endl;
}

void Get_systematics_PtRange(TH1D* &uncertainty, TH1D* &preBarlow, int iDataset,
                            int iRadius, int unfoldParameter, std::string options,
                            bool varyMeasured, TH1D** reference, bool drawSpectra) {
  int nVariations = varyMeasured ? nMeasuredPtRangeVariations : nUnfoldedPtRangeVariations;
  const double (*ranges)[2] = varyMeasured ? measuredPtRangeVariations : unfoldedPtRangeVariations;
  if (nVariations < 1) throw std::runtime_error("At least one pT-range variation is required");
  std::vector<TH1D*> spectra(nVariations + 1);
  std::vector<TString> legends(nVariations + 1);
  TString source = varyMeasured ? "MeasuredPtRange" : "UnfoldedPtRange";
  TString method = options.find("Svd") != std::string::npos ? "Svd" : "Bayes";
  TString name = source + "_" + Datasets[iDataset]
                 + Form("_R=%.1f_%s%d", arrayRadius[iRadius], method.Data(), unfoldParameter);
  double nominalMinimum = varyMeasured ? ptBinsJetsRec[iRadius][0] : ptBinsJetsGen[iRadius][0];
  double nominalMaximum = varyMeasured ? ptBinsJetsRec[iRadius][nBinPtJetsRec[iRadius]]
                                       : ptBinsJetsGen[iRadius][nBinPtJetsGen[iRadius]];
  Get_Systematics_UnfoldedWithPtRange(spectra[0], iDataset, iRadius, unfoldParameter,
                                     options, varyMeasured, nominalMinimum, nominalMaximum);
  legends[0] = Form("Nominal %g-%g", nominalMinimum, nominalMaximum);
  for (int variation = 0; variation < nVariations; ++variation) {
    Get_Systematics_UnfoldedWithPtRange(spectra[variation+1], iDataset, iRadius,
      unfoldParameter, options, varyMeasured, ranges[variation][0], ranges[variation][1]);
    legends[variation+1] = Form("%s %g-%g", varyMeasured ? "Measured" : "Unfolded",
                               ranges[variation][0], ranges[variation][1]);
  }
  Align_Systematics_CommonPtBins(spectra, name);
  Get_Systematics_Envelope(uncertainty, preBarlow, spectra, name);
  if (drawSpectra) Draw_Systematics_Spectra(spectra, legends, "Systematics_Spectra_" + name, iRadius);
  if (reference) *reference = spectra[0];
}

void Draw_Systematics_PtRange(int iDataset, int iRadius, int unfoldParameter,
                             std::string options, bool varyMeasured, TH1D** relativeResult) {
  TH1D *uncertainty, *preBarlow, *reference;
  Get_systematics_PtRange(uncertainty, preBarlow, iDataset, iRadius,
                          unfoldParameter, options, varyMeasured, &reference);
  TString source = varyMeasured ? "MeasuredPtRange" : "UnfoldedPtRange";
  TString method = options.find("Svd") != std::string::npos ? "Svd" : "Bayes";
  TString name = "Systematics_" + source + "_" + Datasets[iDataset]
                 + Form("_R=%.1f_%s%d", arrayRadius[iRadius], method.Data(), unfoldParameter);
  // These small relative uncertainties need longer decimal tick labels.
  float originalLabelSize = gStyle->GetLabelSize("Y");
  gStyle->SetLabelSize(0.035, "Y");
  Draw_Systematics_RelativeUncertainty(uncertainty, preBarlow, reference, name, iRadius);
  gStyle->SetLabelSize(originalLabelSize, "Y");
  if (relativeResult) *relativeResult = uncertainty;
}

// Calculate all four sources internally and draw their bin-by-bin maximum.
// A private ROOT directory isolates intermediate histogram names from existing plots.
TH1D* DrawUnfoldingSys(int iDataset, int iRadius, int nominalSvdK, int bayesIterations) {
  if (iDataset < 0 || iDataset >= nDatasets || iRadius < 0 || iRadius >= nRadius)
    throw std::runtime_error("Invalid dataset or radius for unfolding systematics");
  const char* methods[2] = {"Svd", "Bayes"};
  const int parameters[2] = {nominalSvdK, bayesIterations};
  if (bayesIterations < 1) throw std::runtime_error("Bayes iterations must be positive");
  std::string options = std::string(unfoldingPrior) + ",noUnfoldingDiagnostics";
  TH1D* relative[4] = {nullptr, nullptr, nullptr, nullptr};
  cout << "DrawUnfoldingSys: SVD k=" << nominalSvdK << ", Bayes iterations="
       << bayesIterations << ", Barlow=" << applyBarlowTest << endl;
  {
    TDirectoryFile calculationDirectory("UnfoldingSysCalculation", "Temporary calculation");
    TDirectory::TContext directoryContext(&calculationDirectory);
    TH1D* absolute[4] = {nullptr, nullptr, nullptr, nullptr};
    TH1D* reference[4] = {nullptr, nullptr, nullptr, nullptr};
    TH1D* preBarlow = nullptr;
    Get_systematics_UnfoldMethod(absolute[0], preBarlow, iDataset, iRadius,
      methods, parameters, 2, options, &reference[0], false);
    Get_systematics_SvdParameter(absolute[1], preBarlow, iDataset, iRadius,
      nominalSvdK, options, &reference[1], false);
    Get_systematics_PtRange(absolute[2], preBarlow, iDataset, iRadius,
      nominalSvdK, options + ",Svd", true, &reference[2], false);
    Get_systematics_PtRange(absolute[3], preBarlow, iDataset, iRadius,
      nominalSvdK, options + ",Svd", false, &reference[3], false);
    for (int source = 0; source < 4; ++source) {
      relative[source] = (TH1D*)absolute[source]->Clone(Form("CombinedRelativeInput%d", source));
      relative[source]->SetDirectory(nullptr);
      for (int bin = 1; bin <= relative[source]->GetNbinsX(); ++bin) {
        double denominator = std::abs(reference[source]->GetBinContent(bin));
        if (!std::isfinite(denominator)) throw std::runtime_error("Non-finite nominal spectrum");
        if (denominator == 0. && relative[source]->GetBinContent(bin) != 0.)
          throw std::runtime_error("Cannot form relative uncertainty with zero nominal yield");
        relative[source]->SetBinContent(bin, denominator == 0. ? 0.
          : relative[source]->GetBinContent(bin)/denominator);
        relative[source]->SetBinError(bin, denominator == 0. ? 0.
          : relative[source]->GetBinError(bin)/denominator);
      }
    }
  }
  std::vector<TH1D*> sources = {relative[0], relative[1], relative[2], relative[3]};
  for (TH1D* source : sources)
    if (!source) throw std::runtime_error("Missing relative unfolding-systematics input");
  if (!std::isfinite(relativeSystematicsRange[0]) || !std::isfinite(relativeSystematicsRange[1])
      || relativeSystematicsRange[0] >= relativeSystematicsRange[1])
    throw std::runtime_error("Invalid relative-systematics display range");
  TString name = "UnfoldingCombined_" + Datasets[iDataset]
                 + Form("_R=%.1f", arrayRadius[iRadius])
                 + (applyBarlowTest ? "_Barlow" : "_NoBarlow");
  Align_Systematics_CommonPtBins(sources, name);
  for (TH1D* source : relative) delete source;
  TH1D* combined = (TH1D*)sources[0]->Clone("hSystematicUncertainty_" + name);
  combined->SetDirectory(nullptr);
  combined->Reset("ICES");
  for (int bin = 1; bin <= combined->GetNbinsX(); ++bin) {
    double maximum = 0.;
    for (TH1D* source : sources) {
      double value = source->GetBinContent(bin);
      if (!std::isfinite(value) || value < 0.)
        throw std::runtime_error("Invalid relative unfolding systematic uncertainty");
      maximum = std::max(maximum, value);
    }
    combined->SetBinContent(bin, maximum);
    combined->SetBinError(bin, 0.);
  }

  TCanvas* canvas = new TCanvas("canvas_" + name, "Unfolding systematics", 800, 800);
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.10);
  canvas->SetBottomMargin(0.15);
  canvas->SetTopMargin(0.10);
  canvas->SetGridx();
  canvas->SetGridy();
  TH1* frame = canvas->DrawFrame(ptWindowDisplay[0], relativeSystematicsRange[0],
                                 ptWindowDisplay[1], relativeSystematicsRange[1]);
  frame->SetXTitle(texPtJetGen->Data());
  frame->SetYTitle("Relative systematic uncertainty");
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.3);
  const int colors[4] = {kRed+1, kOrange+7, kGreen+2, kMagenta+1};
  const char* labels[4] = {"Unfolding algorithm", "SVD regularisation",
                            "Measured p_{T} range", "Unfolded p_{T} range"};
  TLegend* legend = new TLegend(0.49, 0.59, 0.93, 0.83);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  for (int source = 0; source < 4; ++source) {
    sources[source]->SetLineColor(colors[source]);
    sources[source]->SetStats(false);
    sources[source]->SetLineWidth(2);
    sources[source]->Draw("HIST SAME");
    legend->AddEntry(sources[source], labels[source], "l");
  }
  combined->SetLineColor(kBlack);
  combined->SetStats(false);
  combined->SetLineWidth(3);
  combined->SetLineStyle(2);
  combined->Draw("HIST SAME");
  legend->AddEntry(combined, "Unfolding total (max)", "l");
  legend->Draw();
  TLatex* context = new TLatex();
  context->SetNDC();
  context->SetTextFont(42);
  context->SetTextSize(0.035);
  context->DrawLatex(0.18, 0.85, texCollisionDataInfo->Data());
  context->DrawLatex(0.18, 0.80, contextJetRadius(arrayRadius[iRadius]).Data());
  for (const char* format : {"pdf", "png", "eps"}) {
    TString folder = TString(format) + "Folder";
    gSystem->mkdir(folder, true);
    canvas->SaveAs(folder + "/Systematics_" + name + "." + format);
  }
  cout << name << ": combined relative maximum = " << combined->GetMaximum() << endl;
  return combined;
}
