#ifndef JETSPECTRUM_OUTPUT_H
#define JETSPECTRUM_OUTPUT_H

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "JetSpectrum_settings.h"
#include <initializer_list>
#include <utility>
#include <memory>
#include <cmath>
#include <stdexcept>

// Start a fresh result file once per macro run, only when saving is enabled.
inline void Create_JetSpectrumResultsFile() {
  if (!saveAnalysisResults) return;
  if (analysisResultsFileName.IsNull()) throw std::runtime_error("Set analysisResultsFileName before saving");
  TDirectory::TContext context;
  TFile output(analysisResultsFileName, "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create analysis result file");
}

// Never transfer plotted histograms to the output file. Write detached copies
// and restore gDirectory, including on exceptions. Fixed keys work for every centrality.
inline void Save_JetSpectrumHistograms(
    std::initializer_list<std::pair<const char*, TH1*>> histograms,
    const char* directoryName = "") {
  if (!saveAnalysisResults) return;
  if (analysisResultsFileName.IsNull()) throw std::runtime_error("Set analysisResultsFileName before saving");
  TDirectory::TContext context;
  TFile output(analysisResultsFileName, "UPDATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot open analysis result file for writing");
  TDirectory* directory = &output;
  if (TString(directoryName).Length()) {
    directory = output.GetDirectory(directoryName);
    if (!directory) directory = output.mkdir(directoryName);
    if (!directory) throw std::runtime_error("Cannot create result directory");
  }
  directory->cd();
  for (const auto& entry : histograms) {
    if (!entry.second) throw std::runtime_error("Missing histogram to save");
    for (int bin = 0; bin <= entry.second->GetNbinsX()+1; ++bin)
      if (!std::isfinite(entry.second->GetBinContent(bin)) || !std::isfinite(entry.second->GetBinError(bin)))
        throw std::runtime_error("Non-finite histogram to save");
    std::unique_ptr<TH1> copy((TH1*)entry.second->Clone(entry.first));
    copy->SetDirectory(nullptr);
    copy->GetXaxis()->SetRange(0, 0);
    copy->SetMinimum(); copy->SetMaximum();
    if (TString(directoryName) == "Systematics_uncertainty")
      copy->GetYaxis()->SetTitle("Relative systematic uncertainty (fraction)");
    if (copy->Write(entry.first, TObject::kOverwrite) <= 0)
      throw std::runtime_error("Failed to write result histogram");
  }
}

// Encode numerical parameters in ROOT keys without rounding distinct settings.
inline TString JetSpectrumVariationNumber(double value) {
  TString token = Form("%.12g", value);
  token.ReplaceAll(".", "p");
  token.ReplaceAll("-", "m");
  token.ReplaceAll("+", "");
  return token;
}

// Store spectra, not fractional uncertainties or percent plotting copies.
// Called before Barlow/envelope calculations and before display-range cropping.
inline void Save_JetSpectrumVariation(TH1* spectrum, const TString& variationName) {
  TString key = "h1_syst_variation_" + variationName;
  Save_JetSpectrumHistograms({{key.Data(), spectrum}}, "Systematics_variation");
}

#endif
