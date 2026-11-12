// Run inside the ALICE RooUnfold environment:
// root -l Jets/JetSpectrum_RunMacro.C
#if defined(__CLING__)
R__LOAD_LIBRARY(libRooUnfold)
#endif

#include "./JetSpectrum_DrawingFunctions.h"
#include "./JetSpectrum_DrawingFunctions.C"
#include "./JetSpectrum_systematics.C"

void SetStyle(Bool_t graypalette=kFALSE);

void JetSpectrum_RunMacro() {
  SetStyle();
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s", unfoldingPrior, unfoldingMethod);
  Create_JetSpectrumResultsFile();

  int iDataset = 0;
  int iRadius = 0;
  int unfoldParameterInputList[2] = {12, 6}; // nominal SVD k, alternative Bayes iterations
  int unfoldParameterInput = unfoldParameterInputList[0];

  // Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  // Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  // Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius, optionsAnalysis);
  Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// For calculate systematics uncertainty ///////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  const char* methods[2] = {"Svd", "Bayes"};
  std::string priorOptions = unfoldingPrior;
  std::string svdOptions = priorOptions + ",Svd";

  // Unfolding algorithm：SVD k=11 vs Bayes 7
  // Draw_Systematics_UnfoldMethod(iDataset, iRadius, methods, unfoldParameterInputList, 2, priorOptions);

  // SVD k±1
  // Draw_Systematics_SvdParameter(iDataset, iRadius, unfoldParameterInput, priorOptions);

  // Measured pT range
  // Draw_Systematics_PtRange(iDataset, iRadius, unfoldParameterInput, svdOptions, true);

  // Unfolded pT range
  // Draw_Systematics_PtRange(iDataset, iRadius, unfoldParameterInput, svdOptions, false);

  // Background fluctuations
  // Draw_Systematics_BackgroundFluctuations(iDataset, iRadius, unfoldParameterInput, svdOptions);

  // Tracking efficiency
  // Draw_Systematics_TrackingEfficiency(iDataset, iRadius, unfoldParameterInput, svdOptions);

  // Model prior
  // Draw_Systematics_Prior(iDataset, iRadius, unfoldParameterInput);

  // Vertex
  // Draw_Systematics_Vertex(iDataset, iRadius, unfoldParameterInput, svdOptions);

  // Centrality variations
  // Draw_Systematics_Centrality(iDataset, iRadius, unfoldParameterInput, svdOptions);

  DrawTotalSys(iDataset, iRadius, unfoldParameterInputList[0], unfoldParameterInputList[1]);
  // TH1D* totalUncertainty = DrawTotalSys(iDataset, iRadius, unfoldParameterInputList[0], unfoldParameterInputList[1]);
}
