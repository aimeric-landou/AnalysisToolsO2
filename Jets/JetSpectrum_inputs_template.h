
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
// TString* texCollisionDataType = new TString("0#minus10% Pb#minusPb"); 
TString* texCollisionDataType = new TString("0#font[122]{-}10% Pb#font[122]{-}Pb"); 
TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+", "+(TString)*texEnergyPbPb); 
TString* texCollisionMCType = new TString("PYTHIA + GEANT4"); 
TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+", "+(TString)*texEnergy); 
const TString* texDatasetsComparisonType = new TString("0#font[122]{-}10% cent.");
// const TString* texDatasetsComparisonType = new TString("50-70% centrality");
const TString* texDatasetsComparisonCommonDenominator = new TString("ALICE performance");

const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_train420439"
                                    };
const TString DatasetsNames[nDatasets] = {"ppAnchorPbPb"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
                                      };
TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_train420439/AnalysisResults.root")
                                                    };
TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_train420439/AnalysisResults.root")
                                                                        }; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)


// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23_PbPb_pass4_goldenRuns_occupancy01000_train372068"
//                                     };
// const TString DatasetsNames[nDatasets] = {"noTrackEffCorrection"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_FullStats_R02_noTrackEffCorrection_train395261/AnalysisResults.root")
//                                                     };
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_FullStats_R02_noTrackEffCorrection_train395261/AnalysisResults.root")
//                                                                         }; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// const TString trainIdData = "_id26156"; // lead05
// const TString analysisWorkflowData = "jet-spectra-charged_central"+trainIdData;
const TString trainIdData = ""; // lead05
const TString analysisWorkflowData = "jet-spectra-charged_lead_05_100"+trainIdData;

const TString trainIdBkg = "";
const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
const TString trainIdUnfoldingControl = "";
const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged_lead_05_100"+trainIdUnfoldingControl;

const TString trainIdMC = "";
const TString analysisWorkflowMC = "jet-spectra-charged_lead_05_100"+trainIdMC;
const bool etaCutOnMatchedJetsIsObsoleteVersion = false;



// //////// -------- LHC23zzh pass 4 with - pp sim anchored to PbPb 10% - lead03 ///////
// TString* texEnergyPbPb = new TString("#sqrt{#it{s}_{NN}} = 5.36 TeV"); 
// TString* texEnergy = new TString("pp, #sqrt{#it{s}} = 5.36 TeV"); 
// // TString* texCollisionDataType = new TString("0#minus10% Pb#minusPb"); 
// TString* texCollisionDataType = new TString("0#font[122]{-}10% Pb#font[122]{-}Pb"); 
// TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+", "+(TString)*texEnergyPbPb); 
// TString* texCollisionMCType = new TString("PYTHIA + GEANT4"); 
// TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+", "+(TString)*texEnergy); 
// const TString* texDatasetsComparisonType = new TString("0#font[122]{-}10% centrality");
// // const TString* texDatasetsComparisonType = new TString("50-70% centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("ALICE performance");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23_PbPb_pass4_goldenRuns_occupancy01000_train372068"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_FullStats_R02_train388528_leadTrack03/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_FullStats_R02_train388528_leadTrack03/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// // const TString trainIdData = "_id26156"; // lead05
// const TString trainIdData = "_id26155"; // lead03
// const TString analysisWorkflowData = "jet-spectra-charged_central"+trainIdData;
// // const TString analysisWorkflowData = "jet-spectra-charged_peripheral"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = false;




// //////// -------- pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o pass7");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC22o_pass7_train261733"};
// const TString DatasetsNames[nDatasets] = {"LHC22o_pass7"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/LHC24f3b_train261768/half1/AnalysisResults.root");

// // TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC24f3_sel8MC_train240962/AnalysisResults.root")};
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC24f3b_train261768/half2/AnalysisResults.root")};

// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/LHC24b1b_sel8Full_train239409/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4_weighted/AnalysisResults.root");

// const TString trainId = "";
// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";
// // const TString analysisWorkflowMC = "jet-finder-charged-qa_global_CollMatch";
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;





// //////// -------- pp spectrum with Joonsuk files tests -------- ////////
// TString* texEnergyPbPb = new TString(""); 
// TString* texEnergy = new TString("#sqrt{#it{s}} = 13.6 TeV"); 
// TString* texCollisionDataType = new TString("Data pp"); 
// TString* texCollisionDataInfo = new TString((TString)*texCollisionDataType+" "+(TString)*texEnergy); 
// TString* texCollisionMCType = new TString("MC pp PYTHIA"); 
// TString* texCollisionMCInfo = new TString((TString)*texCollisionMCType+" "+(TString)*texEnergy); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"Data_LHC22o_minBias_small_Joonsuk"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = {new TFile("Datasets/ppSim_LHC24f3b_withTrackTuner_Joonsuk/AnalysisResults.root")};
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/ppSim_LHC24f3b_withTrackTuner_Joonsuk/AnalysisResults.root")};

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;

// const TString trainId = "";
// const TString analysisWorkflowData = "jet-finder-charged-qa"+trainId;
// const TString analysisWorkflow_unfoldingControl = "jet-finder-charged-qa"+trainId;

// const TString analysisWorkflowMC = "jet-finder-charged-qa";
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;







// //////// -------- Pb-Pb spectrum with Wenhui Angantyr files tests -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"Data_halfMCAngantyr_train361349"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/MC_halfMCAngantyr_train361349/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/MC_halfMCAngantyr_train361349/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdData = "";
// const TString analysisWorkflowData = "jet-spectra-charged"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;








// //////// -------- LHC25b6 - pp sim anchored to PbPb 10%  - just looking at eff and purity ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV - Pb-Pb anchor"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_10percent_train366181_jetSpectrumWorkflow"};
// const TString DatasetsNames[nDatasets] = {"LHC25b6"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_10percent_train366181_jetSpectrumWorkflow/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC25b6_pp_sim_PbPbAnchor_10percent_train366181_jetSpectrumWorkflow/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdData = "";
// const TString analysisWorkflowData = "jet-spectra-charged"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = true;




// //////// -------- test LHC25b6 efficiency finderQA vs spectraCharged ///////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC25b6_localWorkflow"};
// const TString DatasetsNames[nDatasets] = {"spectraCharged"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/ppSim_LHC23d4_weighted_withLeadingTrackCut/AnalysisResults.root");
// // TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/MC_halfMCAngantyr/AnalysisResults.root");
// TFile* file_O2Analysis_MCfileForMatrix[nDatasets] = new TFile("Datasets/LHC25b6_localWorkflow/AnalysisResults.root");
// TFile* file_O2Analysis_ppSimDetectorEffect_unfoldingControl[nDatasets] = {new TFile("Datasets/LHC25b6_localWorkflow/AnalysisResults.root")}; // use this MC file as input to unfolding (with h_jet_pt_rhoareasubtracted distrib on file) and as comparison to gen (with h_jet_pt_part distrib on file)

// // const TString trainId = "_id12832";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_5090_lead5"+trainId;
// // const TString trainId = "_id12436";
// // const TString analysisWorkflowData = "jet-finder-charged-qa_central_0010_lead5"+trainId;
// const TString trainIdData = "_id26156";
// const TString analysisWorkflowData = "jet-spectra-charged_central"+trainIdData;
// const TString trainIdBkg = "";
// const TString analysisWorkflowBkg = "jet-background-analysis"+trainIdBkg;
// const TString trainIdUnfoldingControl = "";
// const TString analysisWorkflow_unfoldingControl = "jet-spectra-charged_noOccupancyCut"+trainIdUnfoldingControl;

// const TString trainIdMC = "";
// const TString analysisWorkflowMC = "jet-spectra-charged"+trainIdMC;
// const bool etaCutOnMatchedJetsIsObsoleteVersion = false;


#endif
