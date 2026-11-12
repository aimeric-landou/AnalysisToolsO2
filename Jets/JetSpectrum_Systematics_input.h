#ifndef JETSPECTRUM_SYSTEMATICS_INPUT_H
#define JETSPECTRUM_SYSTEMATICS_INPUT_H

#include "TFile.h"
#include "TString.h"
#include "JetSpectrum_inputs.h"

// 0-10% Pb-Pb systematics inputs on 2023 Pb-Pb Analysis HP2026 preliminary.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// 00-10%: cent0010 ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////// Background fluctuations //////////////////////////////////////////////////////
TFile* file_O2Analysis_Bkg_Systematics_list[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23zzh_bkg_DelatRs_HY683074/AnalysisResults.root")
};
const TString analysisWorkflowBkg_Systematics = "jet-background-analysis-task";
const TString analysisWorkflowBkg_DeltaR08 = "jet-background-analysis-task_DeltaR08";
const TString analysisWorkflowBkg_DeltaR15 = "jet-background-analysis-task_DeltaR15";

// Same three background variations used in the May analysis.
// randomConeType: 1 = no leading jet, 2 = randomised tracks,
// 3 = randomised tracks without the leading jet (see JetSpectrum_settings.h).
const int nBackgroundVariations = 3;
const TString backgroundVariationWorkflows[nBackgroundVariations] = {
  analysisWorkflowBkg_DeltaR15,
  analysisWorkflowBkg_Systematics,
  analysisWorkflowBkg_Systematics
};
const int backgroundVariationRandomConeTypes[nBackgroundVariations] = {1, 2, 3};
const TString backgroundVariationLabels[nBackgroundVariations] = {
  "RC, #DeltaR > 1.5", "Randomised tracks", "Randomised, no lead"
};

///////////////////////////////////////////////////////// Tracking efficiency //////////////////////////////////////////////////////
TFile* file_O2Analysis_MCfile_TrackingEfficiency_Systematics[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_cent0010_new_TrackingEff_sys_HY688109/AnalysisResults.root")
};
const TString analysisWorkflowMC_TrackingEfficiency_Systematics = "jet-spectra-charged_lead5";

////////////////////////////////////////////////// Event selection: |zVertex| < 8 cm ///////////////////////////////////////////////
TFile* file_O2Analysis_ZVertex8_Systematics_list[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPbGolden_zVertex8_cent0010_1030_3050_HY669452/AnalysisResults_25kHz.root")
};
const TString analysisWorkflowData_ZVertex8_Systematics = "jet-spectra-charged_cent0010";

///////////////////////////////////////////////// Event selection: centrality FT0C Variant1 ////////////////////////////////////////
TFile* file_O2Analysis_Centrality_Systematics_list[nDatasets] = {
  new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPb_gloden_centFT0C_systematics_HY687336/AnalysisResults_centFT0CVariant1_25kHz_train687336.root")
};
const TString analysisWorkflowData_Centrality_Systematics = "jet-spectra-charged_cent0010";

 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// 10-30%: cent1030 ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TFile* file_O2Analysis_Bkg_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23zzh_bkg_DelatRs_HY683074/AnalysisResults.root")
// };
// const TString analysisWorkflowBkg_Systematics = "jet-background-analysis-task";
// const TString analysisWorkflowBkg_DeltaR08 = "jet-background-analysis-task_DeltaR08";
// const TString analysisWorkflowBkg_DeltaR15 = "jet-background-analysis-task_DeltaR15";

// // Same three background variations used in the May analysis.
// // randomConeType: 1 = no leading jet, 2 = randomised tracks,
// // 3 = randomised tracks without the leading jet (see JetSpectrum_settings.h).
// const int nBackgroundVariations = 3;
// const TString backgroundVariationWorkflows[nBackgroundVariations] = {
//   analysisWorkflowBkg_DeltaR15,
//   analysisWorkflowBkg_Systematics,
//   analysisWorkflowBkg_Systematics
// };
// const int backgroundVariationRandomConeTypes[nBackgroundVariations] = {1, 2, 3};
// const TString backgroundVariationLabels[nBackgroundVariations] = {
//   "RC, #DeltaR > 1.5", "Randomised tracks", "Randomised, no lead"
// };
// ///////////////////////////////////////////////////////// Tracking efficiency //////////////////////////////////////////////////////
// TFile* file_O2Analysis_MCfile_TrackingEfficiency_Systematics[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_cent1030_TrackingEff_sys_HY687298/AnalysisResults.root")
// };
// const TString analysisWorkflowMC_TrackingEfficiency_Systematics = "jet-spectra-charged_lead5";
// ////////////////////////////////////////////////// Event selection: |zVertex| < 8 cm ///////////////////////////////////////////////
// TFile* file_O2Analysis_ZVertex8_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPbGolden_zVertex8_cent0010_1030_3050_HY669452/AnalysisResults_25kHz.root")
// };
// const TString analysisWorkflowData_ZVertex8_Systematics = "jet-spectra-charged_cent1030";
// ///////////////////////////////////////////////// Event selection: centrality FT0C Variant1 ////////////////////////////////////////
// TFile* file_O2Analysis_Centrality_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPb_gloden_centFT0C_systematics_HY687336/AnalysisResults_centFT0CVariant1_25kHz_train687336.root")
// };
// const TString analysisWorkflowData_Centrality_Systematics = "jet-spectra-charged_cent1030";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// 30-50%: cent3050 ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TFile* file_O2Analysis_Bkg_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23zzh_bkg_DelatRs_HY683074/AnalysisResults.root")
// };
// const TString analysisWorkflowBkg_Systematics = "jet-background-analysis-task";
// const TString analysisWorkflowBkg_DeltaR08 = "jet-background-analysis-task_DeltaR08";
// const TString analysisWorkflowBkg_DeltaR15 = "jet-background-analysis-task_DeltaR15";

// // Same three background variations used in the May analysis.
// // randomConeType: 1 = no leading jet, 2 = randomised tracks,
// // 3 = randomised tracks without the leading jet (see JetSpectrum_settings.h).
// const int nBackgroundVariations = 3;
// const TString backgroundVariationWorkflows[nBackgroundVariations] = {
//   analysisWorkflowBkg_DeltaR15,
//   analysisWorkflowBkg_Systematics,
//   analysisWorkflowBkg_Systematics
// };
// const int backgroundVariationRandomConeTypes[nBackgroundVariations] = {1, 2, 3};
// const TString backgroundVariationLabels[nBackgroundVariations] = {
//   "RC, #DeltaR > 1.5", "Randomised tracks", "Randomised, no lead"
// };

// ///////////////////////////////////////////////////////// Tracking efficiency //////////////////////////////////////////////////////
// TFile* file_O2Analysis_MCfile_TrackingEfficiency_Systematics[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_cent3050_new_trackingEff_sys_HY688111/AnalysisResults.root")
// };
// const TString analysisWorkflowMC_TrackingEfficiency_Systematics = "jet-spectra-charged_lead5";

// ////////////////////////////////////////////////// Event selection: |zVertex| < 8 cm ///////////////////////////////////////////////
// TFile* file_O2Analysis_ZVertex8_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPbGolden_zVertex8_cent0010_1030_3050_HY669452/AnalysisResults_25kHz.root")
// };
// const TString analysisWorkflowData_ZVertex8_Systematics = "jet-spectra-charged_cent3050";

// ///////////////////////////////////////////////// Event selection: centrality FT0C Variant1 ////////////////////////////////////////
// TFile* file_O2Analysis_Centrality_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPb_gloden_centFT0C_systematics_HY687336/AnalysisResults_centFT0CVariant1_25kHz_train687336.root")
// };
// const TString analysisWorkflowData_Centrality_Systematics = "jet-spectra-charged_cent3050";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////// 50-80%: cent5080 //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TFile* file_O2Analysis_Bkg_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23zzh_bkg_DelatRs_HY683074/AnalysisResults.root")
// };
// const TString analysisWorkflowBkg_Systematics = "jet-background-analysis-task";
// const TString analysisWorkflowBkg_DeltaR08 = "jet-background-analysis-task_DeltaR08";
// const TString analysisWorkflowBkg_DeltaR15 = "jet-background-analysis-task_DeltaR15";

// // Same three background variations used in the May analysis.
// // randomConeType: 1 = no leading jet, 2 = randomised tracks,
// // 3 = randomised tracks without the leading jet (see JetSpectrum_settings.h).
// const int nBackgroundVariations = 3;
// const TString backgroundVariationWorkflows[nBackgroundVariations] = {
//   analysisWorkflowBkg_DeltaR15,
//   analysisWorkflowBkg_Systematics,
//   analysisWorkflowBkg_Systematics
// };
// const int backgroundVariationRandomConeTypes[nBackgroundVariations] = {1, 2, 3};
// const TString backgroundVariationLabels[nBackgroundVariations] = {
//   "RC, #DeltaR > 1.5", "Randomised tracks", "Randomised, no lead"
// };

// ///////////////////////////////////////////////////////// Tracking efficiency //////////////////////////////////////////////////////
// TFile* file_O2Analysis_MCfile_TrackingEfficiency_Systematics[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_jetjetMC/input/JJMC_cent5080_trackingeffSys_2timesRun2_HY684498/AnalysisResults.root")
// };
// const TString analysisWorkflowMC_TrackingEfficiency_Systematics = "jet-spectra-charged_lead5";

// ////////////////////////////////////////////////// Event selection: |zVertex| < 8 cm ///////////////////////////////////////////////
// TFile* file_O2Analysis_ZVertex8_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPbGolden_Zvertex8_cent5080_HY680479/AnalysisResults_25kHz.root")
// };
// const TString analysisWorkflowData_ZVertex8_Systematics = "jet-spectra-charged_cent5080";

// ///////////////////////////////////////////////// Event selection: centrality FT0C Variant1 ////////////////////////////////////////
// TFile* file_O2Analysis_Centrality_Systematics_list[nDatasets] = {
//   new TFile("/Users/fengwenhui/run_hyperloop/pbpb_apass4_spectra/input/23PbPb_gloden_centFT0C_systematics_HY687336/AnalysisResults_centFT0CVariant1_25kHz_train687336.root")
// };
// const TString analysisWorkflowData_Centrality_Systematics = "jet-spectra-charged_cent5080";


#endif
