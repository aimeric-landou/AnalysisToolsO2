#include "./JetSpectrum_DrawingFunctions.h"
#include "./JetSpectrum_DrawingFunctions.C"

using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void JetSpectrum_RunMacro() {
  // Set the default style
  SetStyle();

  // gathers the analysis options in a single char[]
  snprintf(optionsAnalysis, sizeof(optionsAnalysis), "%s,%s", unfoldingPrior, unfoldingMethod);
  cout << "Analysis options are: " << optionsAnalysis << endl;

  mcCollHistIsObsolete = inputMcCollHistIsObsolete;

  int iDataset = 0;
  int iRadius = 0;

  Draw_ResponseMatrices_Fluctuations(iDataset, iRadius);
  Draw_ResponseMatrices_detectorResponse(iDataset, iRadius);
  Draw_ResponseMatrices_DetectorAndFluctuationsCombined(iDataset, iRadius, optionsAnalysis);

  // Draw_ResponseMatrices_Fluctuations(1, iRadius);
  // Draw_ResponseMatrices_detectorResponse(1, iRadius);
  // Draw_ResponseMatrices_DetectorAndFluctuationsCombined(1, iRadius, optionsAnalysis);

  // // Draw_Pt_spectrum_unfolded_FluctResponseOnly(iDataset, iRadius, optionsAnalysis); // NOT FIXED YET - result meaningless
  // Draw_Pt_spectrum_raw(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_raw(iDataset, iRadius, optionsAnalysis+(std::string)"noEventNormNorBinWidthScaling");
  // Draw_Pt_spectrum_mcp(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcp(iDataset, iRadius, optionsAnalysis+(std::string)"noEventNormNorBinWidthScaling");
  // Draw_Pt_spectrum_mcdMatched(iDataset, iRadius, optionsAnalysis);
  // Draw_Pt_spectrum_mcdMatched(iDataset, iRadius, optionsAnalysis+(std::string)"noEventNormNorBinWidthScaling");

  // Draw_Pt_efficiency_jets(iRadius, optionsAnalysis);
  // Draw_kinematicEfficiency(iRadius, optionsAnalysis);
  // Draw_FakeRatio(iRadius, optionsAnalysis);

  // int unfoldParameterInput = 5;
  // Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput, optionsAnalysis);
  // int unfoldParameterInput1 = 6;
  // Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput1, optionsAnalysis);
  // int unfoldParameterInput2 = 8;
  // Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput2, optionsAnalysis);
  // Draw_Pt_spectrum_unfolded_datasetComparison(iRadius, unfoldParameterInput2, optionsAnalysis);
  int unfoldParameterInput3 = 10;
  Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput3, optionsAnalysis);
  // int unfoldParameterInput4 = 12;
  // Draw_Pt_spectrum_unfolded_singleDataset(iDataset, iRadius, unfoldParameterInput4, optionsAnalysis);

  // int unfoldParameterInputMin = 7;
  // int unfoldParameterInputMax = 12;
  // int unfoldParameterInputStep = 2;
  // Draw_Pt_spectrum_unfolded_parameterVariation_singleDataset(iDataset, iRadius, unfoldParameterInputMin, unfoldParameterInputMax, unfoldParameterInputStep, optionsAnalysis);

}

