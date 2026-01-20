
// This is a template. To use the JetSpectrum_DrawingMacro.C, rename this file to JetSpectrum_inputs.h and edit it how you want.


#ifndef JETSPECTRUM_INPUTS_H
#define JETSPECTRUM_INPUTS_H

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////       file access choice       ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
TFile* file_O2Analysis_run2ComparisonFileHannaBossiLaura = new TFile("Datasets/Run2_Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange/Unfolding_AreaBased_HannahMethod_R020_Nominal_ExtendedPtRange.root");
TFile* file_O2Analysis_run2ComparisonFileMLPaper = new TFile("Datasets/Run2_Unfolding_MachineLearningMethod_R020/Ch-jetSuppression_PbPb502TeV.root");


//////// -------- LHC23zzh pass 4 with - pp sim anchored to PbPb 10% - lead05 ///////
TString* texEnergyPbPb = new TString("#sqrt{#it{s}_{NN}} = 5.36 TeV"); 
TString* texEnergy = new TString("pp, #sqrt{#it{s}} = 5.36 TeV"); 
TString* texCollisionDataType = new TString("0#font[122]{-}10% Pb#font[122]{-}Pb"); 
TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+", "+(TString)*texEnergyPbPb); 
TString* texCollisionMCType = new TString("PYTHIA + GEANT4"); 
TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+", "+(TString)*texEnergy); 
const TString* texDatasetsComparisonType = new TString("0#font[122]{-}10% cent.");
const TString* texDatasetsComparisonCommonDenominator = new TString("ALICE performance");

const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_train420439"
                                    };
const TString DatasetsNames[nDatasets] = {"ppAnchorPbPb"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
                                      };
TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_train420439/AnalysisResults.root")
                                                    };
TFile* file_O2Analysis_MCfile_UnfoldingControl_input[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_0010centEffCorrection_train571118_MCUnfoldingInput/AnalysisResults.root")
                                                                        }; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

TFile* file_O2Analysis_MCfile_UnfoldingControl_response[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_0010centEffCorrection_train571118_MCResponse/AnalysisResults.root")
                                                                        }; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)


const TString trainIdData = "";
const TString analysisWorkflowData = "jet-spectra-charged_lead_05_100"+trainIdData;

const TString trainIdBkg = "";
const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
const TString trainIdUnfoldingControl = "";
const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged_lead_05_100"+trainIdUnfoldingControl;

const TString trainIdMC = "";
const TString analysisWorkflowMC = "jet-spectra-charged_lead_05_100"+trainIdMC;
const bool etaCutOnMatchedJetsIsObsoleteVersion = false;
bool inputMcCollHistIsObsolete = true;


#endif
