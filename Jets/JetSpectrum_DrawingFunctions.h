#ifndef JETSPECTRUM_DRAWINGFUNCTIONS_H
#define JETSPECTRUM_DRAWINGFUNCTIONS_H


void IterationLegend(TString* iterationLegend, int unfoldIterationMin, int unfoldIterationMax, int step);

void Draw_ResponseMatrices_Fluctuations(int iDataset, int iRadius);
void Draw_ResponseMatrices_detectorResponse(int iDataset, int iRadius);
void Draw_ResponseMatrices_DetectorAndFluctuationsCombined(int iDataset, int iRadius, std::string options);

void Draw_Pt_spectrum_raw(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_mcp(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_mcdMatched(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_FluctResponseOnly(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_singleDataset(int iDataset, int iRadius, int unfoldParameterInput, std::string options);
void Draw_Pt_spectrum_unfolded_datasetComparison(int iRadius, int unfoldParameterInput, std::string options);
void Draw_Pt_TestSpectrum_unfolded(int iDataset, int iRadius, std::string options);
void Draw_Pt_spectrum_unfolded_parameterVariation_singleDataset(int iDataset, int iRadius, int unfoldIterationMin, int unfoldIterationMax, int step, std::string options);

void Draw_Pt_efficiency_jets(int iRadius, std::string options);
void Draw_kinematicEfficiency(int iRadius, std::string options);
void Draw_FakeRatio(int iRadius, std::string options);

#endif