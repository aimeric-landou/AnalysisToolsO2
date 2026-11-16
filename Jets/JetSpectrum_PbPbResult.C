////////////////////////////////////////////////////////////////////////////////////////////////////
// STANDALONE FINAL-RESULT PLOTTING MACRO -- NOT PART OF THE RunMacro EXECUTION FLOW
//
// Purpose:
//   Read completed analysis results and produce the final preliminary jet-spectrum plots,
//   Rcp plots, Rcp systematic-uncertainty summaries, and the ALICE--ATLAS Rcp comparison.
//   The plotting style follows the original JetUnfoldResults.C preliminary figures.
//   This macro does not run unfolding or generate the required centrality input files.
//
// REQUIRED BEFORE USE:
//   First run JetSpectrum_RunMacro.C separately for each of the four centrality classes:
//   0-10%, 10-30%, 30-50%, and 50-80%, with the matching inputs and unfolding settings.
//   Enable saveAnalysisResults, choose a different output filename for each centrality,
//   and call DrawAndSave_JetSpectrumResultsAndSystematics(...) to save both the nominal spectra and the complete systematic uncertainties and variation spectra.
//
//   All four completed ROOT files must be available in pbpbResultInputFolder, for example:
//     JetSpectrum_results_lead5_25kHz_cent0010.root
//     JetSpectrum_results_lead5_25kHz_cent1030.root
//     JetSpectrum_results_lead5_25kHz_cent3050.root
//     JetSpectrum_results_lead5_25kHz_cent5080.root
//   Each file must contain HJetSpectra_measured, HJetSpectra_Svdkdefault,
//   HJetSpectra_refolded, Systematics_uncertainty/, and Systematics_variation/.
//   The local ATLAS reference ROOT file configured in pbpbRcpAtlasFile is also required
//   when running Draw_JetSpectrumPbPbRcp(), which automatically draws the ATLAS comparison.
//
// WHEN TO RUN:
//   Use this macro only after unfolding and systematic checks are complete and all four centrality result files have been saved successfully; it is the final plotting stage.
//   Check the input paths, Ncoll values, and plotting settings below, then run independently
//   from the parent directory of AnalysisToolsO2:
//     root -l AnalysisToolsO2/Jets/JetSpectrum_PbPbResult.C
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef JETSPECTRUM_PBPBRESULT_C
#define JETSPECTRUM_PBPBRESULT_C

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TLine.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

// Function declarations
void JetSpectrum_PbPbResult();
TCanvas* Draw_JetSpectrumPbPbSpectra();
int Find_JetSpectrumPbPbBin(const TH1* histogram, double low, double high);
double Get_JetSpectrumPbPbMeasuredRelativeError(const TH1* measured, double low, double high);
TH1* Read_JetSpectrumPbPbHistogram(TFile* file, TString key);
TH1* Find_JetSpectrumPbPbVariation(TFile* file, TString prefix);
TH1D* Get_JetSpectrumPbPbRcp(TFile* centralFile, TFile* peripheralFile, TH1* centralSpectrum, TH1* peripheralSpectrum, int centrality, bool useMeasuredStatistics = false);
TH1D* Calculate_JetSpectrumPbPbRcpSystematics(int centrality, TH1D** sourceResults = nullptr);
TH1D* Draw_JetSpectrumPbPbRcpSystematics(int centrality);
void Draw_JetSpectrumPbPbRcpUncertainty(int centrality, TH1D** sources, TH1D* total);
TCanvas* Draw_JetSpectrumPbPbRcp();
void Save_JetSpectrumPbPbCanvas(TCanvas* canvas, TString name);

//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Global settings ////////////////////////////////////////
//////////// Independent of the single-centrality RunMacro settings and input files.//////////
//////////////////////////////////////////////////////////////////////////////////////////////

TString pbpbResultInputFolder = "rootfile_folder";
const int nPbPbResultCentralities = 4;
const TString pbpbResultInputFiles[nPbPbResultCentralities] = {
  "JetSpectrum_results_lead5_25kHz_cent0010.root",
  "JetSpectrum_results_lead5_25kHz_cent1030.root",
  "JetSpectrum_results_lead5_25kHz_cent3050.root",
  "JetSpectrum_results_lead5_25kHz_cent5080.root"
};
// Exactly the active Ncoll values in JetUnfoldResults.C (not its commented alternatives).
double pbpbResultNcoll[nPbPbResultCentralities] = {1597.73, 785.13, 258.14, 47.38};
const double pbpbResultMinimumPt[nPbPbResultCentralities] = {40., 30., 30., 30.};
const double pbpbResultMaximumPt = 140.;
bool pbpbResultUseMeasuredRelativeStatistics = true; // original preliminary convention
bool savePbPbResultPlots = true;
TString pbpbResultPlotName = "Corrected_Jetpt_lead5_centrality_Ncoll";

const int pbpbResultColors[nPbPbResultCentralities] = {kRed, kOrange+1, kBlue, kGreen+2};
const int pbpbResultMarkers[nPbPbResultCentralities] = {kFullCircle, kFullSquare, kFullCross, kFullCrossX};
const double pbpbResultMarkerSizes[nPbPbResultCentralities] = {1.2, 1.0, 1.2, 1.2};
const char* pbpbResultCentralityLabels[nPbPbResultCentralities] = {
  "0#font[122]{-}10%", "10#font[122]{-}30%",
  "30#font[122]{-}50%", "50#font[122]{-}80%"
};

bool drawPbPbRcpSystematics = true;
const int pbpbRcpPeripheralIndex = 3; // 50-80%
const TString pbpbResultCentralityCodes[4] = {"0010", "1030", "3050", "5080"};
double pbpbRcpNcollUncertaintyPercent[3] = {5., 5., 5.}; // separate from the returned total and Rcp bands
double pbpbRcpSystematicsRange[2] = {0., 30.};
TString pbpbRcpPlotName = "CorrRcp_lead5_centrality";
TString pbpbRcpAtlasFile = "/Users/fengwenhui/public_plots/JetSizeDependence_JetSuppression_PbPb276_ATLAS/HEPData-ins1126965-v1-root.root";
TString pbpbRcpAtlasPlotName = "CorrRcp_lead5_cent0010_5080_ATLASCom";

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Main macro //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void JetSpectrum_PbPbResult() {
  Draw_JetSpectrumPbPbSpectra();
  Draw_JetSpectrumPbPbRcp();
}


// Function implementations
int Find_JetSpectrumPbPbBin(const TH1* histogram, double low, double high) {
  int bin = histogram->GetXaxis()->FindFixBin((low + high)/2.);
  if (bin < 1 || bin > histogram->GetNbinsX() || std::abs(histogram->GetXaxis()->GetBinLowEdge(bin)-low) > 1e-6 || std::abs(histogram->GetXaxis()->GetBinUpEdge(bin)-high) > 1e-6)
    throw std::runtime_error("Incompatible physical pT bin edges in saved results");
  return bin;
}

double Get_JetSpectrumPbPbMeasuredRelativeError(const TH1* measured, double low, double high) {
  double yield = 0., variance = 0., next = low;
  for (int bin = 1; bin <= measured->GetNbinsX(); ++bin) {
    double left = measured->GetXaxis()->GetBinLowEdge(bin);
    double right = measured->GetXaxis()->GetBinUpEdge(bin);
    if (right <= low+1e-6 || left >= high-1e-6) continue;
    if (left < low) left = low;
    if (right > high) right = high;
    if (std::abs(left-next) > 1e-6)
      throw std::runtime_error("Measured bins do not cover unfolded pT bin");
    double width = right-left;
    yield += measured->GetBinContent(bin)*width;
    variance += std::pow(measured->GetBinError(bin)*width, 2);
    next = right;
  }
  if (std::abs(next-high) > 1e-6 || !std::isfinite(yield) || yield <= 0. || !std::isfinite(variance)) throw std::runtime_error("Invalid measured statistical input");
  return std::sqrt(variance)/yield;
}

TCanvas* Draw_JetSpectrumPbPbSpectra() {
  TH1D* spectra[4];
  TGraphAsymmErrors* bands[4];
  double ncollRelative = -1.;
  bool commonNcollUncertainty = true;

  for (int cent = 0; cent < 4; ++cent) {
    TString inputPath = pbpbResultInputFolder + "/" + pbpbResultInputFiles[cent];
    TFile* inputFile = TFile::Open(inputPath, "READ");
    if (!inputFile || inputFile->IsZombie())
      throw std::runtime_error("Cannot open PbPb result file");

    TH1* measured = (TH1*)inputFile->Get("HJetSpectra_measured");
    TH1D* unfolded = (TH1D*)inputFile->Get("HJetSpectra_Svdkdefault");
    TH1* relativeSystematics = (TH1*)inputFile->Get("Systematics_uncertainty/h1_syst_total");
    TH1* relativeNcoll = (TH1*)inputFile->Get("Systematics_uncertainty/h1_syst_ncoll");
    if (!measured || !unfolded || !relativeSystematics || !relativeNcoll)
      throw std::runtime_error("Missing histogram in PbPb result file");
    if (!std::isfinite(pbpbResultNcoll[cent]) || pbpbResultNcoll[cent] <= 0.)
      throw std::runtime_error("Ncoll must be positive");

    spectra[cent] = (TH1D*)unfolded->Clone(Form("jetSpectra_cent%d_Ncoll", cent));
    spectra[cent]->SetDirectory(nullptr);
    spectra[cent]->Reset("ICES");
    spectra[cent]->SetStats(false);
    spectra[cent]->SetMarkerColor(pbpbResultColors[cent]);
    spectra[cent]->SetLineColor(pbpbResultColors[cent]);
    spectra[cent]->SetMarkerStyle(pbpbResultMarkers[cent]);
    spectra[cent]->SetMarkerSize(pbpbResultMarkerSizes[cent]);
    bands[cent] = new TGraphAsymmErrors();
    bands[cent]->SetName(Form("Graph_sys_cent%d", cent));
    for (int bin = 1; bin <= spectra[cent]->GetNbinsX(); ++bin) {
      double low = spectra[cent]->GetXaxis()->GetBinLowEdge(bin);
      double high = spectra[cent]->GetXaxis()->GetBinUpEdge(bin);
      if (low < pbpbResultMinimumPt[cent]-1e-6 || high > pbpbResultMaximumPt+1e-6) continue;
      int sysBin = Find_JetSpectrumPbPbBin(relativeSystematics, low, high);
      int ncollBin = Find_JetSpectrumPbPbBin(relativeNcoll, low, high);
      double value = unfolded->GetBinContent(bin);
      double stat = unfolded->GetBinError(bin);
      double relative = relativeSystematics->GetBinContent(sysBin);
      double normalisation = relativeNcoll->GetBinContent(ncollBin);
      if (!std::isfinite(value) || value <= 0. || !std::isfinite(stat) || stat < 0.
          || !std::isfinite(relative) || relative < 0.
          || !std::isfinite(normalisation) || normalisation < 0.)
        throw std::runtime_error("Invalid spectrum or fractional uncertainty");
      if (pbpbResultUseMeasuredRelativeStatistics)
        stat = value * Get_JetSpectrumPbPbMeasuredRelativeError(measured, low, high);
      value = value / pbpbResultNcoll[cent];
      stat = stat / pbpbResultNcoll[cent];
      spectra[cent]->SetBinContent(bin, value);
      spectra[cent]->SetBinError(bin, stat);
      int point = bands[cent]->GetN();
      bands[cent]->SetPoint(point, (low+high)/2., value);
      bands[cent]->SetPointError(point, (high-low)/2., (high-low)/2., relative*value, relative*value);
      if (ncollRelative < 0.) ncollRelative = normalisation;
      else if (std::abs(normalisation-ncollRelative) > 1e-8) commonNcollUncertainty = false;
    }
    bands[cent]->SetFillColorAlpha(pbpbResultColors[cent], 0.3);
    bands[cent]->SetLineWidth(2);
    bands[cent]->SetLineColor(pbpbResultColors[cent]);
    std::cout << "PbPb spectrum " << pbpbResultCentralityLabels[cent] << ": " << bands[cent]->GetN()
              << " points, Ncoll=" << pbpbResultNcoll[cent] << std::endl;
    inputFile->Close();
    delete inputFile;
  }

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas* canvas = new TCanvas("Cspectra_cent", "Cspectra_cent", 1000, 1000);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.02);
  canvas->SetLeftMargin(0.185);
  canvas->SetRightMargin(0.05);
  canvas->SetLogy();
  TH1F* frame = canvas->DrawFrame(30, 3e-10, 140, 3e-5);
  frame->SetTitle(" ");
  frame->GetXaxis()->SetTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
  frame->GetYaxis()->SetTitle("#frac{1}{#LT #it{N}_{coll} #GT} #frac{1}{#it{N}_{evt}} #frac{d^{2}#it{N}}{d#it{p}_{T, ch jet} d#it{#eta}_{jet} } (GeV/#it{c})^{#font[122]{-}1}");
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetXaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.95);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->Draw("");
  for (int cent = 0; cent < 4; ++cent) {
    bands[cent]->Draw("2 same");
  }
  for (int cent = 0; cent < 4; ++cent) {
    spectra[cent]->Draw("lep same");
  }
  TGraphAsymmErrors* legendBand = new TGraphAsymmErrors();
  legendBand->SetName("Graph_sys_legend");
  legendBand->SetFillColorAlpha(kGray+2, 0.5);
  legendBand->SetLineWidth(0);
  legendBand->SetLineColor(kGray+2);
  TLegend* legend = new TLegend(0.55, 0.60, 0.88, 0.78);
  legend->AddEntry((TObject*)nullptr, "Centrality", "");
  for (int cent = 0; cent < 4; ++cent) {
    legend->AddEntry(spectra[cent], pbpbResultCentralityLabels[cent], "lep");
  }
  legend->AddEntry(legendBand, "syst. uncertainty", "f");
  legend->SetBorderSize(0);
  legend->SetTextSize(0.033);
  legend->Draw("same");
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.035);
  text.SetTextFont(42);
  text.DrawLatex(0.21, 0.92, "ALICE Preliminary");
  text.DrawLatex(0.21, 0.87, "Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  text.DrawLatex(0.21, 0.83, "Ch-particle jets, anti-#it{k}_{T}, #it{R} = 0.2, |#it{#eta}_{jet}| < 0.7");
  text.DrawLatex(0.21, 0.78, "#it{p}_{T, track}^{lead.} > 5 GeV/#it{c}");
  text.SetTextSize(0.03);
  if (commonNcollUncertainty) {
    text.DrawLatex(0.57, 0.57, Form("#LT #it{N}_{coll} #GT uncertainty #pm %g%%", 100.*ncollRelative));
  } else {
    text.DrawLatex(0.57, 0.57, "#LT #it{N}_{coll} #GT uncertainty included");
  }

  if (savePbPbResultPlots) {
    gSystem->mkdir("pngFolder/PbPbResult", true);
    gSystem->mkdir("pdfFolder/PbPbResult", true);
    gSystem->mkdir("epsFolder/PbPbResult", true);
    canvas->SaveAs("pngFolder/PbPbResult/" + pbpbResultPlotName + ".png");
    canvas->SaveAs("pdfFolder/PbPbResult/" + pbpbResultPlotName + ".pdf");
    canvas->SaveAs("epsFolder/PbPbResult/" + pbpbResultPlotName + ".eps");
  }
  return canvas;
}

TH1* Read_JetSpectrumPbPbHistogram(TFile* file, TString key) {
  TH1* histogram = dynamic_cast<TH1*>(file->Get(key));
  if (!histogram || histogram->GetDimension() != 1)
    throw std::runtime_error(std::string("Missing histogram: ") + file->GetName() + ":" + key.Data());
  return histogram;
}

TH1* Find_JetSpectrumPbPbVariation(TFile* file, TString prefix) {
  TDirectory* directory = file->GetDirectory("Systematics_variation");
  if (!directory) throw std::runtime_error("Missing Systematics_variation directory");
  TString fullPrefix = "h1_syst_variation_" + prefix;
  TString matchedKey;
  int matches = 0;
  TIter next(directory->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    if (TString(key->GetName()).BeginsWith(fullPrefix)) {
      matchedKey = key->GetName();
      ++matches;
    }
  }
  if (matches != 1)
    throw std::runtime_error(std::string("Missing or ambiguous variation: ") + file->GetName() + ":" + fullPrefix.Data());
  return Read_JetSpectrumPbPbHistogram(file, "Systematics_variation/" + matchedKey);
}

TH1D* Get_JetSpectrumPbPbRcp(TFile* centralFile, TFile* peripheralFile, TH1* centralSpectrum, TH1* peripheralSpectrum, int centrality, bool useMeasuredStatistics) {
  if (centrality < 0 || centrality >= 3) throw std::runtime_error("Rcp numerator must be 0-10, 10-30 or 30-50%");
  double edges[256];
  int nBins = 0;
  if (centralSpectrum->GetNbinsX() > 255) throw std::runtime_error("Too many Rcp bins");
  for (int bin = 1; bin <= centralSpectrum->GetNbinsX(); ++bin) {
    double low = centralSpectrum->GetXaxis()->GetBinLowEdge(bin);
    double high = centralSpectrum->GetXaxis()->GetBinUpEdge(bin);
    if (low < pbpbResultMinimumPt[centrality]-1e-6 || high > pbpbResultMaximumPt+1e-6) continue;
    if (nBins > 0 && std::abs(edges[nBins]-low) > 1e-6) throw std::runtime_error("Non-contiguous Rcp bins");
    edges[nBins] = low;
    edges[++nBins] = high;
  }
  if (nBins == 0) throw std::runtime_error("No Rcp bins in final pT range");
  if (pbpbResultNcoll[centrality] <= 0. || pbpbResultNcoll[pbpbRcpPeripheralIndex] <= 0.)
    throw std::runtime_error("Ncoll must be positive");
  TH1D* ratio = new TH1D("h1_Rcp_cent" + pbpbResultCentralityCodes[centrality], "", nBins, edges);
  ratio->SetDirectory(nullptr);
  TH1* measuredCentral = nullptr;
  TH1* measuredPeripheral = nullptr;
  if (useMeasuredStatistics) {
    measuredCentral = Read_JetSpectrumPbPbHistogram(centralFile, "HJetSpectra_measured");
    measuredPeripheral = Read_JetSpectrumPbPbHistogram(peripheralFile, "HJetSpectra_measured");
  }
  double normalisation = pbpbResultNcoll[pbpbRcpPeripheralIndex] / pbpbResultNcoll[centrality];
  for (int bin = 1; bin <= nBins; ++bin) {
    double low = edges[bin-1];
    double high = edges[bin];
    int centralBin = Find_JetSpectrumPbPbBin(centralSpectrum, low, high);
    int peripheralBin = Find_JetSpectrumPbPbBin(peripheralSpectrum, low, high);
    double central = centralSpectrum->GetBinContent(centralBin);
    double peripheral = peripheralSpectrum->GetBinContent(peripheralBin);
    if (!std::isfinite(central) || !std::isfinite(peripheral) || central <= 0. || peripheral <= 0.)
      throw std::runtime_error("Invalid Rcp numerator or denominator");
    double centralRelativeError = centralSpectrum->GetBinError(centralBin) / central;
    double peripheralRelativeError = peripheralSpectrum->GetBinError(peripheralBin) / peripheral;
    if (useMeasuredStatistics) {
      centralRelativeError = Get_JetSpectrumPbPbMeasuredRelativeError(measuredCentral, low, high);
      peripheralRelativeError = Get_JetSpectrumPbPbMeasuredRelativeError(measuredPeripheral, low, high);
    }
    double value = central / peripheral * normalisation;
    double error = value * std::sqrt(centralRelativeError*centralRelativeError + peripheralRelativeError*peripheralRelativeError);
    if (!std::isfinite(value) || !std::isfinite(error)) throw std::runtime_error("Non-finite Rcp value or error");
    ratio->SetBinContent(bin, value);
    ratio->SetBinError(bin, error);
  }
  return ratio;
}

const int nPbPbRcpSources = 9;
const TString pbpbRcpSourceNames[9] = {
  "UnfoldingMethod", "SvdParameter", "MeasuredPtRange", "UnfoldedPtRange",
  "Prior", "BackgroundFluctuations", "TrackingEfficiency", "Vertex", "Centrality"
};
const char* pbpbRcpSourceLabels[11] = {
  "Unfolding algorithm", "SVD regularisation", "Measured p_{T} range", "Unfolded p_{T} range",
  "Prior", "Background", "Tracking efficiency", "Vertex selection", "Centrality selection",
  "Unfolding (max)", "#LTN_{coll}#GT"
};
const int pbpbRcpSourceColors[11] = {
  kRed+1, kOrange+7, kGreen+2, kMagenta+1,
  kGreen+2, kBlue, kCyan+2, kViolet+2, kMagenta, kRed, kOrange-1
};
const int pbpbRcpVariationCounts[9] = {1, 2, 3, 2, 2, 3, 1, 1, 1};
const TString pbpbRcpVariationPrefixes[9][3] = {
  {"UnfoldingMethod_Bayes_ite", "", ""},
  {"SvdParameter_kMinus1_k", "SvdParameter_kPlus1_k", ""},
  {"MeasuredPtRange_Var1_", "MeasuredPtRange_Var2_", "MeasuredPtRange_Var3_"},
  {"UnfoldedPtRange_Var1_", "UnfoldedPtRange_Var2_", ""},
  {"Prior_Var1_", "Prior_Var2_", ""},
  {"BackgroundFluctuations_jet_background_analysis_task_DeltaR15_RCtype1",
   "BackgroundFluctuations_jet_background_analysis_task_RCtype2",
   "BackgroundFluctuations_jet_background_analysis_task_RCtype3"},
  {"TrackingEfficiency_Varied", "", ""},
  {"Vertex_ZVertexMax8cm", "", ""},
  {"Centrality_FT0CVariant1", "", ""}
};

TH1D* Calculate_JetSpectrumPbPbRcpSystematics(int centrality, TH1D** sourceResults) {
  if (centrality < 0 || centrality >= 3) throw std::runtime_error("Invalid Rcp centrality");
  if (!std::isfinite(pbpbRcpNcollUncertaintyPercent[centrality]) || pbpbRcpNcollUncertaintyPercent[centrality] < 0.)
    throw std::runtime_error("Invalid Rcp Ncoll uncertainty");
  TFile* centralFile = TFile::Open(pbpbResultInputFolder + "/" + pbpbResultInputFiles[centrality], "READ");
  TFile* peripheralFile = TFile::Open(pbpbResultInputFolder + "/" + pbpbResultInputFiles[pbpbRcpPeripheralIndex], "READ");
  if (!centralFile || centralFile->IsZombie() || !peripheralFile || peripheralFile->IsZombie())
    throw std::runtime_error("Cannot open Rcp inputs");
  TH1* central = Read_JetSpectrumPbPbHistogram(centralFile, "HJetSpectra_Svdkdefault");
  TH1* peripheral = Read_JetSpectrumPbPbHistogram(peripheralFile, "HJetSpectra_Svdkdefault");
  TH1D* nominal = Get_JetSpectrumPbPbRcp(centralFile, peripheralFile, central, peripheral, centrality);
  TH1D* sources[11];
  for (int source = 0; source < 11; ++source) {
    TString sourceName;
    if (source < 9) sourceName = pbpbRcpSourceNames[source];
    else if (source == 9) sourceName = "Unfolding";
    else sourceName = "Ncoll";
    sources[source] = (TH1D*)nominal->Clone("h1_syst_Rcp_" + sourceName + "_cent" + pbpbResultCentralityCodes[centrality]);
    sources[source]->SetDirectory(nullptr);
    sources[source]->Reset("ICES");
  }
  for (int source = 0; source < nPbPbRcpSources; ++source) {
    TH1D* variations[3] = {nullptr, nullptr, nullptr};
    for (int variation = 0; variation < pbpbRcpVariationCounts[source]; ++variation) {
      TString prefix = pbpbRcpVariationPrefixes[source][variation];
      TH1* variedCentral = Find_JetSpectrumPbPbVariation(centralFile, prefix);
      TH1* variedPeripheral = Find_JetSpectrumPbPbVariation(peripheralFile, prefix);
      std::cout << "Rcp cent" << pbpbResultCentralityCodes[centrality] << " "
                << variedCentral->GetName() << " / " << variedPeripheral->GetName() << std::endl;
      variations[variation] = Get_JetSpectrumPbPbRcp(centralFile, peripheralFile,
        variedCentral, variedPeripheral, centrality);
      variations[variation]->SetName(Form("h1_Rcp_%s_var%d_cent", pbpbRcpSourceNames[source].Data(), variation+1)
        + pbpbResultCentralityCodes[centrality]);
      for (int bin = 1; bin <= nominal->GetNbinsX(); ++bin) {
        double low = nominal->GetXaxis()->GetBinLowEdge(bin);
        double high = nominal->GetXaxis()->GetBinUpEdge(bin);
        int variedBin = Find_JetSpectrumPbPbBin(variations[variation], low, high);
        double difference = std::abs(variations[variation]->GetBinContent(variedBin) - nominal->GetBinContent(bin));
        double relative = difference / nominal->GetBinContent(bin);
        if (relative > sources[source]->GetBinContent(bin)) sources[source]->SetBinContent(bin, relative);
      }
    }
    for (int variation = 0; variation < pbpbRcpVariationCounts[source]; ++variation) delete variations[variation];
  }
  TH1D* total = (TH1D*)nominal->Clone("h1_syst_Rcp_total_cent" + pbpbResultCentralityCodes[centrality]);
  total->SetDirectory(nullptr);
  total->Reset("ICES");
  for (int bin = 1; bin <= nominal->GetNbinsX(); ++bin) {
    double unfolding = 0.;
    for (int source = 0; source < 4; ++source)
      if (sources[source]->GetBinContent(bin) > unfolding) unfolding = sources[source]->GetBinContent(bin);
    sources[9]->SetBinContent(bin, unfolding);
    sources[10]->SetBinContent(bin, pbpbRcpNcollUncertaintyPercent[centrality]/100.);
    double variance = unfolding * unfolding;
    for (int source = 4; source < 9; ++source) variance += std::pow(sources[source]->GetBinContent(bin), 2);
    // Ncoll is displayed separately as in the original preliminary Rcp figure.
    total->SetBinContent(bin, std::sqrt(variance));
  }
  if (sourceResults) {
    for (int source = 0; source < 11; ++source) sourceResults[source] = sources[source];
  }
  centralFile->Close();
  peripheralFile->Close();
  delete centralFile;
  delete peripheralFile;
  delete nominal;
  return total;
}

void Save_JetSpectrumPbPbCanvas(TCanvas* canvas, TString name) {
  if (!savePbPbResultPlots) return;
  gSystem->mkdir("pngFolder/PbPbResult", true);
  gSystem->mkdir("pdfFolder/PbPbResult", true);
  gSystem->mkdir("epsFolder/PbPbResult", true);
  canvas->SaveAs("pngFolder/PbPbResult/" + name + ".png");
  canvas->SaveAs("pdfFolder/PbPbResult/" + name + ".pdf");
  canvas->SaveAs("epsFolder/PbPbResult/" + name + ".eps");
}

void Draw_JetSpectrumPbPbRcpUncertainty(int centrality, TH1D** sources, TH1D* total) {
  TString name = "Rcp_TotalSystematics_cent" + pbpbResultCentralityCodes[centrality];
  TCanvas* canvas = new TCanvas("canvas_" + name, name, 1000, 1000);
  canvas->SetBottomMargin(.12);
  canvas->SetTopMargin(.03);
  canvas->SetLeftMargin(.14);
  canvas->SetRightMargin(.05);
  canvas->SetGridx();
  canvas->SetGridy();
  TH1F* frame = canvas->DrawFrame(30, pbpbRcpSystematicsRange[0], 140, pbpbRcpSystematicsRange[1]);
  frame->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
  frame->SetYTitle("R_{CP} relative systematic uncertainty (%)");
  frame->GetYaxis()->SetTitleSize(.035);
  frame->GetYaxis()->SetTitleOffset(1.6);
  TLegend* legend = new TLegend(.56, .70, .90, .95);
  legend->SetBorderSize(0);
  legend->SetTextSize(.024);
  for (int source = 4; source < 11; ++source) {
    TH1D* percent = (TH1D*)sources[source]->Clone(TString(sources[source]->GetName()) + "_Percent");
    percent->SetDirectory(nullptr);
    percent->Scale(100.);
    percent->SetStats(false);
    percent->SetLineWidth(2);
    percent->SetLineColor(pbpbRcpSourceColors[source]);
    percent->Draw("hist same");
    legend->AddEntry(percent, pbpbRcpSourceLabels[source], "l");
  }
  TH1D* percent = (TH1D*)total->Clone(TString(total->GetName()) + "_Percent");
  percent->SetDirectory(nullptr);
  percent->Scale(100.);
  percent->SetStats(false);
  percent->SetLineColor(kBlack);
  percent->SetLineWidth(2);
  percent->Draw("hist same");
  legend->AddEntry(percent, "Total (excl. N_{coll})", "l");
  legend->Draw();
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(.03);
  text.DrawLatex(.17, .92, "R_{CP} systematic uncertainty");
  text.DrawLatex(.17, .87, TString("(") + pbpbResultCentralityLabels[centrality] + ")/(50#font[122]{-}80%)");
  Save_JetSpectrumPbPbCanvas(canvas, name);
}

TH1D* Draw_JetSpectrumPbPbRcpSystematics(int centrality) {
  TH1D* sources[11];
  TH1D* total = Calculate_JetSpectrumPbPbRcpSystematics(centrality, sources);
  Draw_JetSpectrumPbPbRcpUncertainty(centrality, sources, total);
  std::cout << "Rcp cent" << pbpbResultCentralityCodes[centrality]
            << " maximum systematic = " << 100.*total->GetMaximum() << "% (excluding Ncoll)" << std::endl;
  return total;
}

TCanvas* Draw_JetSpectrumPbPbRcp() {
  TH1D* ratios[3];
  TGraphAsymmErrors* bands[3];
  TGraphAsymmErrors* normBoxes[3];
  for (int cent = 0; cent < 3; ++cent) {
    TFile* centralFile = TFile::Open(pbpbResultInputFolder + "/" + pbpbResultInputFiles[cent], "READ");
    TFile* peripheralFile = TFile::Open(pbpbResultInputFolder + "/" + pbpbResultInputFiles[pbpbRcpPeripheralIndex], "READ");
    if (!centralFile || centralFile->IsZombie() || !peripheralFile || peripheralFile->IsZombie())
      throw std::runtime_error("Cannot open nominal Rcp spectra");
    TH1* central = Read_JetSpectrumPbPbHistogram(centralFile, "HJetSpectra_Svdkdefault");
    TH1* peripheral = Read_JetSpectrumPbPbHistogram(peripheralFile, "HJetSpectra_Svdkdefault");
    ratios[cent] = Get_JetSpectrumPbPbRcp(centralFile, peripheralFile, central, peripheral, cent, pbpbResultUseMeasuredRelativeStatistics);
    ratios[cent]->SetLineColor(pbpbResultColors[cent]);
    ratios[cent]->SetMarkerColor(pbpbResultColors[cent]);
    ratios[cent]->SetMarkerStyle(pbpbResultMarkers[cent]);
    ratios[cent]->SetMarkerSize(pbpbResultMarkerSizes[cent]);
    TH1D* uncertainty;
    if (drawPbPbRcpSystematics) uncertainty = Draw_JetSpectrumPbPbRcpSystematics(cent);
    else uncertainty = Calculate_JetSpectrumPbPbRcpSystematics(cent);
    bands[cent] = new TGraphAsymmErrors();
    bands[cent]->SetName("Graph_Rcp_sys_cent" + pbpbResultCentralityCodes[cent]);
    for (int bin = 1; bin <= ratios[cent]->GetNbinsX(); ++bin) {
      double x = ratios[cent]->GetBinCenter(bin);
      double y = ratios[cent]->GetBinContent(bin);
      double width = ratios[cent]->GetBinWidth(bin)/2.;
      bands[cent]->SetPoint(bin-1, x, y);
      bands[cent]->SetPointError(bin-1, width, width,
        y*uncertainty->GetBinContent(bin), y*uncertainty->GetBinContent(bin));
    }
    bands[cent]->SetFillColorAlpha(pbpbResultColors[cent], .3);
    bands[cent]->SetLineColor(pbpbResultColors[cent]);
    bands[cent]->SetLineWidth(cent == 0 ? 0 : 2);
    normBoxes[cent] = new TGraphAsymmErrors();
    normBoxes[cent]->SetName(Form("Graph_Rcp_Ncoll_cent%d", cent));
    normBoxes[cent]->SetPoint(0, 35, 1.);
    normBoxes[cent]->SetPointError(0, 3., 3., pbpbRcpNcollUncertaintyPercent[cent]/100., pbpbRcpNcollUncertaintyPercent[cent]/100.);
    normBoxes[cent]->SetFillColor(kGray+2);
    normBoxes[cent]->SetFillStyle(3345);
    normBoxes[cent]->SetLineWidth(0);
    centralFile->Close();
    peripheralFile->Close();
    delete centralFile;
    delete peripheralFile;
  }
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas* canvas = new TCanvas("Crcp_alice", "Crcp_alice", 1000, 1000);
  canvas->SetBottomMargin(.12);
  canvas->SetTopMargin(.02);
  canvas->SetLeftMargin(.14);
  canvas->SetRightMargin(.05);
  TH1F* frame = canvas->DrawFrame(25, .3, 145, 1.4);
  frame->SetTitle(" ");
  frame->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
  frame->SetYTitle("#it{R}_{CP}");
  frame->GetYaxis()->SetTitleSize(.045);
  frame->GetXaxis()->SetTitleSize(.04);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetLabelSize(.035);
  frame->GetXaxis()->SetLabelSize(.035);
  frame->GetXaxis()->SetNdivisions(508);
  TLine* unity = new TLine(35, 1., 145, 1.);
  unity->SetLineColor(kBlack);
  unity->SetLineStyle(9);
  unity->SetLineWidth(2);
  unity->Draw("same");
  for (int cent = 0; cent < 3; ++cent) {
    bands[cent]->Draw("2 same");
    ratios[cent]->Draw("lep same");
  }
  normBoxes[0]->Draw("2 same");
  bool sameNormalisation = true;
  for (int cent = 1; cent < 3; ++cent)
    if (pbpbRcpNcollUncertaintyPercent[cent] != pbpbRcpNcollUncertaintyPercent[0]) sameNormalisation = false;
  if (!sameNormalisation) {
    for (int cent = 1; cent < 3; ++cent) {
      normBoxes[cent]->SetPoint(0, 35.+7.*cent, 1.);
      normBoxes[cent]->SetFillColor(pbpbResultColors[cent]);
      normBoxes[cent]->Draw("2 same");
    }
  }
  TGraphAsymmErrors* legendBand = new TGraphAsymmErrors();
  legendBand->SetFillColorAlpha(kGray+2, .3);
  legendBand->SetLineWidth(0);
  TLegend* legend = new TLegend(.54, .72, .90, .92);
  for (int cent = 0; cent < 3; ++cent)
    legend->AddEntry(ratios[cent], TString("(") + pbpbResultCentralityLabels[cent] + ")/(50#font[122]{-}80%)", "lep");
  legend->AddEntry(legendBand, "syst. uncertainty", "f");
  legend->AddEntry(normBoxes[0], "#LT #it{N}_{coll} #GT uncertainty", "f");
  legend->SetBorderSize(0);
  legend->SetTextSize(.035);
  legend->Draw();
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(.035);
  text.DrawLatex(.17, .92, "ALICE Preliminary");
  text.DrawLatex(.17, .88, "Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  text.DrawLatex(.17, .835, "Ch-particle jets");
  text.DrawLatex(.17, .79, "anti-#it{k}_{T}, #it{R} = 0.2, |#it{#eta}_{jet}| < 0.7");
  text.DrawLatex(.17, .75, "#it{p}_{T, track}^{lead.} > 5 GeV/#it{c}");
  Save_JetSpectrumPbPbCanvas(canvas, pbpbRcpPlotName);

  TFile* atlasFile = TFile::Open(pbpbRcpAtlasFile, "READ");
  if (!atlasFile || atlasFile->IsZombie()) throw std::runtime_error("Cannot open ATLAS Rcp file");
  TH1* atlas = (TH1*)Read_JetSpectrumPbPbHistogram(atlasFile, "Table 2/Hist1D_y1")->Clone("ATLAS_RcpR02_cent0010");
  atlas->SetDirectory(nullptr);
  TH1* atlasStat = Read_JetSpectrumPbPbHistogram(atlasFile, "Table 2/Hist1D_y1_e1");
  TH1* atlasE2Up = Read_JetSpectrumPbPbHistogram(atlasFile, "Table 2/Hist1D_y1_e2plus");
  TH1* atlasE2Down = Read_JetSpectrumPbPbHistogram(atlasFile, "Table 2/Hist1D_y1_e2minus");
  TH1* atlasE3Up = Read_JetSpectrumPbPbHistogram(atlasFile, "Table 2/Hist1D_y1_e3plus");
  TH1* atlasE3Down = Read_JetSpectrumPbPbHistogram(atlasFile, "Table 2/Hist1D_y1_e3minus");
  TGraphAsymmErrors* atlasE2 = new TGraphAsymmErrors();
  TGraphAsymmErrors* atlasE3 = new TGraphAsymmErrors();
  atlasE2->SetName("Graph_ATLAS_Rcp_FullyCorrelated");
  atlasE3->SetName("Graph_ATLAS_Rcp_PartiallyCorrelated");
  for (int bin = 1; bin <= atlas->GetNbinsX(); ++bin) {
    double x = atlas->GetBinCenter(bin);
    double y = atlas->GetBinContent(bin);
    double ex = atlas->GetBinWidth(bin)/2.;
    // Preserve the error conversion and asymmetric ordering of JetUnfoldResults.C.
    atlas->SetBinError(bin, atlasStat->GetBinContent(bin)*y);
    atlasE2->SetPoint(bin-1, x, y);
    atlasE2->SetPointError(bin-1, ex, ex,
      std::abs(atlasE2Up->GetBinContent(bin))*y, std::abs(atlasE2Down->GetBinContent(bin))*y);
    atlasE3->SetPoint(bin-1, x, y);
    atlasE3->SetPointError(bin-1, ex, ex,
      std::abs(atlasE3Up->GetBinContent(bin))*y, std::abs(atlasE3Down->GetBinContent(bin))*y);
  }
  atlasFile->Close();
  delete atlasFile;

  TCanvas* comparison = new TCanvas("crcp_cent0", "crcp_cent0", 1000, 1000);
  comparison->SetBottomMargin(.12);
  comparison->SetTopMargin(.02);
  comparison->SetLeftMargin(.14);
  comparison->SetRightMargin(.05);
  TH1F* comparisonFrame = comparison->DrawFrame(35, 0, 145, 1.2);
  comparisonFrame->SetTitle(" ");
  comparisonFrame->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
  comparisonFrame->SetYTitle("#it{R}_{CP}");
  comparisonFrame->GetYaxis()->SetTitleSize(.045);
  comparisonFrame->GetXaxis()->SetTitleSize(.04);
  comparisonFrame->GetYaxis()->SetTitleOffset(1.2);
  comparisonFrame->GetXaxis()->SetTitleOffset(1.2);
  comparisonFrame->GetYaxis()->SetLabelSize(.035);
  comparisonFrame->GetXaxis()->SetLabelSize(.035);
  comparisonFrame->GetXaxis()->SetNdivisions(508);
  TLine* comparisonUnity = new TLine(35, 1., 145, 1.);
  comparisonUnity->SetLineColor(kBlack);
  comparisonUnity->SetLineStyle(9);
  comparisonUnity->SetLineWidth(2);
  comparisonUnity->Draw("same");
  gStyle->SetHatchesLineWidth(3);
  TGraphAsymmErrors* aliceNorm = (TGraphAsymmErrors*)normBoxes[0]->Clone("Graph_ALICE_ATLASComparison_Ncoll");
  aliceNorm->SetPoint(0, 137, 1.);
  aliceNorm->SetPointError(0, 2., 2., pbpbRcpNcollUncertaintyPercent[0]/100., pbpbRcpNcollUncertaintyPercent[0]/100.);
  aliceNorm->Draw("2 same");

  atlas->SetStats(false);
  atlas->SetMarkerColor(kGray+2);
  atlas->SetLineColor(kGray+2);
  atlas->SetMarkerStyle(53);
  atlas->SetMarkerSize(.8);
  atlasE2->SetLineColor(kGray+2);
  atlasE2->SetLineWidth(2);
  atlasE2->SetFillColorAlpha(0, 0);
  atlasE2->SetMarkerSize(0);
  atlasE3->SetFillColorAlpha(kBlue, .5);
  atlasE3->SetLineWidth(2);
  atlasE3->SetLineColor(kBlue);
  atlasE2->Draw("5 same");
  atlasE3->Draw("2 same");
  atlas->Draw("lep same");
  bands[0]->Draw("2 same");
  TH1D* alice = (TH1D*)ratios[0]->Clone("ALICE_RcpR02_cent0010_ATLASComparison");
  alice->SetDirectory(nullptr);
  alice->SetMarkerSize(1.);
  alice->Draw("lep same");

  TLegend* aliceLegend = new TLegend(.4, .58, .88, .76);
  aliceLegend->AddEntry(alice, "ALICE Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.36 TeV", "lep");
  aliceLegend->AddEntry((TObject*)nullptr, "(0#font[122]{-}10%)/(50#font[122]{-}80%), |#it{#eta}_{jet}| < 0.7", "");
  aliceLegend->AddEntry(bands[0], "ALICE syst. uncertainty", "f");
  aliceLegend->AddEntry(aliceNorm, "ALICE #LT #it{N}_{coll} #GT uncertainty", "f");
  aliceLegend->SetBorderSize(0);
  aliceLegend->SetTextSize(.028);
  aliceLegend->Draw();
  TLegend* atlasLegend = new TLegend(.4, .14, .88, .37);
  atlasLegend->AddEntry(atlas, "ATLAS Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "lep");
  atlasLegend->AddEntry((TObject*)nullptr, "(0#font[122]{-}10%)/(60#font[122]{-}80%)", "");
  atlasLegend->AddEntry((TObject*)nullptr, "Charged + Neutral jets, |#it{#eta}_{jet}| < 2.1", "");
  atlasLegend->AddEntry(atlasE2, "ATLAS fully correlated syst.", "f");
  atlasLegend->AddEntry(atlasE3, "ATLAS partially correlated syst.", "f");
  atlasLegend->AddEntry((TObject*)nullptr, "#it{Phys. Lett. B} 719, 220#font[122]{-}241(2013)", "");
  atlasLegend->SetBorderSize(0);
  atlasLegend->SetTextSize(.028);
  atlasLegend->Draw();
  TLatex comparisonText;
  comparisonText.SetNDC();
  comparisonText.SetTextFont(42);
  comparisonText.SetTextSize(.035);
  comparisonText.DrawLatex(.18, .91, "ALICE Preliminary");
  comparisonText.DrawLatex(.18, .86, "Ch-particle jets, anti-#it{k}_{T}, #it{R} = 0.2");
  Save_JetSpectrumPbPbCanvas(comparison, pbpbRcpAtlasPlotName);
  return canvas;
}

#endif
