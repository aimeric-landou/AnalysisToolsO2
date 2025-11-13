// This is a template. To use the TrackMcQC.C, rename this file to TrackMcQC_inputs.h and edit it how you want.

// //////// -------- Track Sel check rough -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_GlobalWithOptionalTpc", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalWithOptionalTpc", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";




// //////// -------- Track Sel check Detailed-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_intermTrackSel2", "LHC23zzh_intermTrackSel1bis", "LHC23zzh_intermTrackSel1", "LHC23zzh_GlobalTpcOptional", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "intermTrackSel2", "intermTrackSel1bis", "intermTrackSel1", "GlobalTpcOptional", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Track Sel check small -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMTracks", "LHC23zzh_GlobalTracks"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";


// //////// -------- Track Sel check three -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMTracks", "LHC23zzh_GlobalTpcOptional", "LHC23zzh_GlobalTracks"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTpcOptional", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";



// //////// -------- Apass vs Cpasses comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Cpass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass1", "LHC23zzh_cpass11", "LHC23zzh_cpass10", "LHC23zzh_cpass8", "LHC23zzh_cpass1", "LHC23zzh_cpass0"};
// const TString DatasetsNames[nDatasets] = {"apass1", "cpass11", "cpass10", "cpass8", "cpass1", "cpass0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Apass - single -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Apass2");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_run544122"};
// const TString DatasetsNames[nDatasets] = {"apass2"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Track efficeincy calculation comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23k6d_AllChecksAndSel8Global", "LHC23k6d_AllChecksAndSel8"};
// const TString DatasetsNames[nDatasets] = {"globalTracks", "MMsel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";

// //////// -------- Track efficeincy calculation comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_MMsel_Cent0010", "LHC23zzh_MMsel_Cent5090", "LHC23zzh_Global_Cent0010", "LHC23zzh_Global_Cent5090"};
// const TString DatasetsNames[nDatasets] = {"MMsel 0-10%", "MMsel 50-90%", "globalTracks 0-10%", "globalTracks 50-90%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "track-efficiency";



// //////// -------- Track efficeincy calculation comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24d2", "LHC24d2"};
// const TString DatasetsNames[nDatasets] = {"LHC24d2 - globalTracks", "LHC24d2 - uniformTracks"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency_uniformTracks"
//                                           };


// //////// -------- Track efficiency calculation centrality comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23k6d");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23k6d_Global_cent0010", "LHC23k6d_Global_cent5090", "LHC23k6d_Uniform_cent0010", "LHC23k6d_Uniform_cent5090"};
// const TString DatasetsNames[nDatasets] = {"Global 0-10%", "Global 50-90%", "Uniform 0-10%", "Uniform 50-90%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency"
//                                           };


// //////// -------- LHC24b1 vs ppJetAnchored -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh_apass3");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"AnalysisResults_LHC24b1_noTrackSel", "AnalysisResults_ppJets_likeLHC24b1_noTrackSel", "AnalysisResults_unanchored_ppJets"};
// const TString DatasetsNames[nDatasets] = {"LHC24b1", "ppJetsAnchored", "ppJetsUnanchored"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+".root"),
//                                         new TFile(Datasets[1]+".root"),
//                                         new TFile(Datasets[2]+".root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"qa-event-track/Tracks",
//                                           "qa-event-track/Tracks",
//                                           "qa-event-track/Tracks"
//                                           };

// //////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut 2cm-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 15;
// const TString Datasets[nDatasets] = {"LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_2", "LHC24d2_1"};
// const TString DatasetsNames[nDatasets] = {"uniform DCAz 0.005", "uniform DCAz 0.01", "uniform DCAz 0.03", "uniform DCAz 0.05", "uniform DCAz 0.1", "uniform DCAz 0.2", "uniform DCAz 0.5", "uniform DCAz 1.0", "uniform DCAz 2.0", "global DCAz 0.005", "global DCAz 0.01", "global DCAz 0.03", "global DCAz 0.05", "global DCAz 1.0", "global DCAz 2.0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[10]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[11]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[12]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[13]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[14]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAz0005",
//                                             "track-efficiency_uniformTracks_DCAxyz001",
//                                             "track-efficiency_uniformTracks_DCAz003",
//                                             "track-efficiency_uniformTracks_DCAxy005_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy01_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy02_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy05_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy1_id13159",
//                                             "track-efficiency_id13159",
//                                             "track-efficiency_globalTracks_DCAz0005",
//                                             "track-efficiency_globalTracks_DCAz001",
//                                             "track-efficiency_globalTracks_DCAz003",
//                                             "track-efficiency_globalTracks_DCAxy005_id13171",
//                                             "track-efficiency_globalTracks_DCAxy1_id13171",
//                                             "track-efficiency_id13171",
//                                           };


// //////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut sans 2cm-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 13;
// const TString Datasets[nDatasets] = {"LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_1", "LHC24d2_2", "LHC24d2_3", "LHC24d2_3", "LHC24d2_3", "LHC24d2_2", "LHC24d2_2"};
// const TString DatasetsNames[nDatasets] = {"uniform DCAz 0.005", "uniform DCAz 0.01", "uniform DCAz 0.03", "uniform DCAz 0.05", "uniform DCAz 0.1", "uniform DCAz 0.2", "uniform DCAz 0.5", "uniform DCAz 1.0", "global DCAz 0.005", "global DCAz 0.01", "global DCAz 0.03", "global DCAz 0.05", "global DCAz 1.0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[10]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[11]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[12]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAz0005",
//                                             "track-efficiency_uniformTracks_DCAxyz001",
//                                             "track-efficiency_uniformTracks_DCAz003",
//                                             "track-efficiency_uniformTracks_DCAxy005_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy01_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy02_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy05_id13159",
//                                             "track-efficiency_uniformTracks_DCAxy1_id13159",
//                                             "track-efficiency_globalTracks_DCAz0005",
//                                             "track-efficiency_globalTracks_DCAz001",
//                                             "track-efficiency_globalTracks_DCAz003",
//                                             "track-efficiency_globalTracks_DCAxy005_id13171",
//                                             "track-efficiency_globalTracks_DCAxy1_id13171"
//                                           };

                                        

// //////// -------- Track efficiency calculation dcaz cut comparison - avec DCAz cut sans 2cm-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("EfficiencyCalc");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24d2_2", "LHC24d2_2"};
// const TString DatasetsNames[nDatasets] = {"uniform", "global"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_uniformTracks_DCAxy1_id13159",
//                                             "track-efficiency_globalTracks_DCAxy1_id13171"
//                                           };


// //////// -------- Run 3 MC local -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("global tracks");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2 - global tracks");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC24d2_global"};
// const TString DatasetsNames[nDatasets] = {"global"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_globalTracks_DCAz20"
//                                           };



// //////// -------- MC pp anchored PbPB vs pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("anchor");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet sim");
// const int nDatasets = 3;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_local"};
// const TString DatasetsNames[nDatasets] = {"PbPbAnchored", "ppAnchored", "LHC24f3 pp"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                             "track-efficiency",
//                                             "track-efficiency"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true, false};

// // jaime testJetAnchoredGapGen-20240801-155530

// //////// -------- MC pp anchored PbPB vs pp -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("includePartonEvent");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet sim");
// const int nDatasets = 2;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetPbPb_includePartonEvent"};
// const TString DatasetsNames[nDatasets] = {"includePartonEvent off", "includePartonEvent on"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true};

// // jaime testJetAnchoredGapGen-20240801-155530


// //////// -------- PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2 - anchored to apass3");
// const int nDatasets = 4;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"LHC24d2_PbPbMCapass4Anchored_train258203", "LHC24d2_PbPbMCapass4Anchored_train258203", "LHC24d2_PbPbMCapass4Anchored_train258203", "LHC24d2_PbPbMCapass4Anchored_train258203"};
// const TString DatasetsNames[nDatasets] = {"0-10%", "10-30%", "30-50%", "50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010",
//                                             "track-efficiency_central1030",
//                                             "track-efficiency_central3050",
//                                             "track-efficiency_central5070"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};

// // jaime testJetAnchoredGapGen-20240801-155530



// //////// -------- accept split vs not accepting split collisions - PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("splitCollisions");
// const TString* texDatasetsComparisonCommonDenominator = new TString("noSplit");
// const int nDatasets = 8;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3", "LHC24g3"};
// const TString DatasetsNames[nDatasets] = {"noSplit 0-10%", "noSplit 10-30%", "noSplit 30-50%", "noSplit 50-70%", "okSplit 0-10%", "okSplit 10-30%", "okSplit 30-50%", "okSplit 50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010_id16386",
//                                             "track-efficiency_central1030_id16386",
//                                             "track-efficiency_central3050_id16840",
//                                             "track-efficiency_central5070_id16840",
//                                             "track-efficiency_central0010_id16837",
//                                             "track-efficiency_central1030_id16837",
//                                             "track-efficiency_central3050_id16841",
//                                             "track-efficiency_central5070_id16841",
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false, false, false, false};

// // jaime testJetAnchoredGapGen-20240801-155530


//////// -------- accept split vs not accepting split collisions - PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("splitCollisions");
// const TString* texDatasetsComparisonCommonDenominator = new TString("noSplit");
// const int nDatasets = 1;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"LHC24g3"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"noSplit 0-10%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults.root")};
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010_id16386"};

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";

// //////// -------- local pp test sim anchored to PbPb comparison to pp anchor -------- ////////
// TString* texCollisionDataInfo = new TString("PbPb MC #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("ppAnchorPbPb");
// const TString* texDatasetsComparisonCommonDenominator = new TString("sim");
// const int nDatasets = 1;
// // const TString Datasets[nDatasets] = {"LHC23zzh_apass3_global", "LHC23zzh_apass3_uniform", "LHC23zzh_apass3_itsOnly", "LHC23zzh_apass3_globalWithoutTpcCrossedRows"};
// // const TString Datasets[nDatasets] = {"anchoredJetJetPbPb", "anchoredJetJetpp", "LHC24f3_train240962", "LHC24f3_local", "unanchoredJetJet_train230486"};
// const TString Datasets[nDatasets] = {"pp_sim_anchored_to_PbPb_5360GeV"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"PbPb anchor (5.36TeV)"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                           };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency"
//                                             };

// //////// -------- accept split vs not accepting split collisions - PbPb MC anchored to LHC23zzh apass4 -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24g3 - sel8FullPbPb");
// const int nDatasets = 8;
// // const TString Datasets[nDatasets] = {"LHC24g3","LHC24d2"};
// // const TString Datasets[nDatasets] = {"LHC24g3"};
//  const TString Datasets[nDatasets] = {"LHC24g3","LHC24d2","LHC24g3","LHC24d2","LHC24g3","LHC24d2","LHC24g3","LHC24d2"};
// const TString DatasetsNames[nDatasets] = {"0-10% pass4","0-10% pass3","10-30% pass4","10-30% pass3","30-50% pass4","30-50% pass3","50-70% pass4","50-70% pass3"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[1]+"/AnalysisResults_LHC24d2_centrality.root"),
//                                         new TFile(Datasets[2]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[3]+"/AnalysisResults_LHC24d2_centrality.root"),
//                                         new TFile(Datasets[4]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[5]+"/AnalysisResults_LHC24d2_centrality.root"),
//                                         new TFile(Datasets[4]+"/AnalysisResults_LHC24g3_centrality.root"),
//                                         new TFile(Datasets[5]+"/AnalysisResults_LHC24d2_centrality.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"track-efficiency","track-efficiency"};
// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central0010","track-efficiency_central0010","track-efficiency_central1030","track-efficiency_central1030","track-efficiency_central3050","track-efficiency_central3050","track-efficiency_central5070","track-efficiency_central5070"};

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false,false,false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";



// //////// -------- pp sim test Anchor Pb-Pb gap comparison-------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("gap value");
// const TString* texDatasetsComparisonCommonDenominator = new TString("MC anchored to Pb-Pb");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"pp_sim_anchored_to_PbPb_5360GeV_gap2", "pp_sim_anchored_to_PbPb_5360GeV_gap3", "pp_sim_anchored_to_PbPb_5360GeV_gap4", "pp_sim_anchored_to_PbPb_5360GeV_gap5"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"gap2", "gap3", "gap4", "gap5"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency",
//                                           "track-efficiency"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- LHC25b6 - pp sim anchored to PbPb 10% ////////
// TString* texCollisionDataInfo = new TString("0-10% Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC25b6_pp_sim_PbPbAnchor_10percent_train370115", "LHC25b4a_pp_ref_Tracks_train371362"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet MC Pb-Pb anchor", "ppRef gen.purp. MC"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;



// //////// -------- Efficiency differences between pp and PbPb 0-10% ////////
// TString* texCollisionDataInfo = new TString("0-10% #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24g3_PbPb_sim_FullStats_R02_train367457_IR25kHz", "LHC25b6_pp_sim_PbPbAnchor_FullStats_R02_train397699_IR25kHz_bugfixed"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"Pb-Pb gen.purp.", "jet-jet Pb-Pb anchor"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_central5070",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false, true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;






// //////// -------- Efficiency differences between pp sims ////////
// TString* texCollisionDataInfo = new TString("0-10% #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("e");
// const TString* texDatasetsComparisonCommonDenominator = new TString("eee");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV", "jetjet_PbPbAnchorMC_5TeV_testWeights", "jetjet_ppAnchorMC_13TeV", "ppRefGenPurposeMC_5TeV", "genPurp_PbPb_sim_train367457_IRall"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet Pb-Pb anchor", "jet-jet Pb-Pb anchor weighted", "jet-jet pp anchor", "pp ref gen.purp.", "Pb-Pb gen.purp."};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency_central5070"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true, true, false, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;



// //////// -------- Efficiency differences between pp and PbPb 0-10% ////////
// TString* texCollisionDataInfo = new TString("0-10% #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("e");
// const TString* texDatasetsComparisonCommonDenominator = new TString("eee");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet Pb-Pb anchor"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;




// //////// -------- Efficiency differences between different ptHat bins ////////
// TString* texCollisionDataInfo = new TString("#sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("ptHatBin");
// const TString* texDatasetsComparisonCommonDenominator = new TString("nonWeighted");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_005050_ptHat", "jetjet_PbPbAnchorMC_5TeV_050100_ptHat", "jetjet_PbPbAnchorMC_5TeV_100200_ptHat", "jetjet_PbPbAnchorMC_5TeV_200300_ptHat"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jj Pb-Pb anchor 5-50 pTHat", "jj Pb-Pb anchor 50-100 pTHat", "jj Pb-Pb anchor 100-200 pTHat", "jj Pb-Pb anchor 200-300 pTHat"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;




// //////// -------- Efficiency differences between different ptHat bins ////////
// TString* texCollisionDataInfo = new TString("#sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("ptHatBin");
// const TString* texDatasetsComparisonCommonDenominator = new TString("weighted");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_005050_ptHat_weighted", "jetjet_PbPbAnchorMC_5TeV_050100_ptHat_weighted", "jetjet_PbPbAnchorMC_5TeV_100200_ptHat_weighted", "jetjet_PbPbAnchorMC_5TeV_200300_ptHat_weighted"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jj Pb-Pb anchor 5-50 pTHat", "jj Pb-Pb anchor 50-100 pTHat", "jj Pb-Pb anchor 100-200 pTHat", "jj Pb-Pb anchor 200-300 pTHat"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true, true, true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;




// //////// -------- Efficiency differences between pp Pb-Pb sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_train428451", "genPurp_PbPb_sim_train367457_IRall", "genPurp_PbPb_sim_train367457_IRall"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"pp jet-jet Pb-Pb anchor", "Pb-Pb gen.purp. 00-10%", "Pb-Pb gen.purp. 50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_id29229",
//                                              "track-efficiency_central0010",
//                                              "track-efficiency_central5070"
//                                           };



// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Efficiency differences between pp sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"jetjet_useTrueTrackWeight", "jetjet_useFalseTrackWeight"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jetjet_useTrueTrackWeight", "jetjet_useFalseTrackWeight" };
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency",
//                                              "track-efficiency"
//                                           };

// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Efficiency differences between pp sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jetjet pp PbPbAnchor 5.36TeV" };
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_id29229"
//                                           };

// const TString wagonId[nDatasets] = {""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Efficiency differences between pp only sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("pp simulation");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_train428451", "ppRefGenPurposeMC_5TeV_train429398", "jetjet_ppAnchorMC_13TeV_train429402", "ppGenPurposeMC_13TeV_train428885"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet 5.36TeV fullstats", "ref gen.purp. 5.36TeV", "jet-jet 13.6TeV", "gen.purp. 13.6TeV"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_id29229",
//                                              "track-efficiency_id31072",
//                                              "track-efficiency_id30280",
//                                              "track-efficiency_id29252"
//                                           };
// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false, true, false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs"; //
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Efficiency differences between pp only sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("pp simulation");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"ppRefGenPurposeMC_5TeV_train429398", "ppGenPurposeMC_13TeV_train428885"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"5.36TeV pp ref gen.purp.", "13.6TeV gen.purp."};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {
//                                              "track-efficiency_id31072",
//                                              "track-efficiency_id29252"
//                                           };
// const TString wagonId[nDatasets] = {"",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false, false};
// const std::string histDatasetComparisonStructure = ""; //
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Efficiency differences between different settings pp PbPb anchor ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_train427391", "jetjet_PbPbAnchorMC_5TeV_train427391", "jetjet_PbPbAnchorMC_5TeV_train427391"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"trueWeight/ptHatCut", "noTrueWeight/ptHatCut", "trueWeight/noPtHatCut"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults_train424259.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_useTrueTrackWeight_id29229",
//                                              "track-efficiency_id29229",
//                                              "track-efficiency_id29229"
//                                           };
// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true, true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Efficiency differences between pp Pb-Pb sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("ptHatFractionMax");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_train428451", "jetjet_PbPbAnchorMC_5TeV_train428451", "jetjet_PbPbAnchorMC_5TeV_train428451", "jetjet_PbPbAnchorMC_5TeV_train428451"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"fraction 4", "fraction 3.5", "fraction 3", "fraction 2.5"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults_pthat4.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_id29229",
//                                              "track-efficiency_ptHatFraction35",
//                                              "track-efficiency_ptHatFraction30",
//                                              "track-efficiency_ptHatFraction25"
//                                           };



// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true, true, true};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;



// //////// -------- Efficiency differences between pp only sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("ptHatFractionMax");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet pp simulation");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"jetjet_ppAnchorMC_13TeV_train429402", "jetjet_ppAnchorMC_13TeV_train429402", "jetjet_ppAnchorMC_13TeV_train429402", "jetjet_ppAnchorMC_13TeV_train429402"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet 2.5 ptHatFractionMax", "jet-jet 3.0 ptHatFractionMax", "jet-jet 3.5 ptHatFractionMax", "jet-jet 4.0 ptHatFractionMax"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_ptHatFraction25_id30280",
//                                              "track-efficiency_ptHatFraction30_id30280",
//                                              "track-efficiency_ptHatFraction35_id30280",
//                                              "track-efficiency_id30280"
//                                           };
// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true, true, true};
// const std::string histDatasetComparisonStructure = ""; //
// const bool datasetsAreSubsetsofId0 = false;


// //////// -------- Variation of track cuts ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("trackCutVariation");
// const TString* texDatasetsComparisonCommonDenominator = new TString("jet-jet pp simulation");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"jetGlobalTrack", "globalRecreationPartial_noChange", "globalRecreationPartial_ITSHits", "globalRecreationPartial_ITSHitsDCAz", "globalRecreationPartial_ITSHitsDCAzChi2ITS", "globalRecreationPartial_ITSHitsDCAzChi2ITSTPC", "globalRecreationPartial_ITSHitsDCAzChi2ITSTPCCrossedRowsOverFindable", "globalRecreationPartial_ITSHitsDCAzChi2ITSTPCCrossedRowsOverFindableCrossedRows"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet globalTrack", "noChange", "ITSHits", "ITSHitsDCAz", "ITSHitsDCAzChi2ITS", "ITSHitsDCAzChi2ITSTPC", "ITSHitsDCAzChi2ITSTPC+CrossedRowsOverFindable", "ITSHitsDCAzChi2ITSTPC+CrossedRowsOverFindable+CrossedRows"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"), 
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency","track-efficiency"
//                                           };
// const TString wagonId[nDatasets] = {"","","","","","","",""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, true, true, true, true, true, true, true};
// const std::string histDatasetComparisonStructure = ""; //
// const bool datasetsAreSubsetsofId0 = false;



// //////// -------- Efficiency differences between pp only sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("pp simulation");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_train428451", "ppRefGenPurposeMC_5TeV_train429398", "jetjet_ppAnchorMC_13TeV_train429402", "ppGenPurposeMC_13TeV_train428885"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"jet-jet 5.36TeV", "ref gen.purp. 5.36TeV", "jet-jet 13.6TeV", "gen.purp. 13.6TeV"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_id29229",
//                                              "track-efficiency_id31072",
//                                              "track-efficiency_id30280",
//                                              "track-efficiency_id29252"
//                                           };
// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false, true, false};
// const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs"; //
// const bool datasetsAreSubsetsofId0 = false;



// //////// -------- Efficiency differences between pp Pb-Pb sims ////////
// TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("simType");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"jetjet_PbPbAnchorMC_5TeV_train428451", "genPurp_PbPb_sim_train367457_IRall", "genPurp_PbPb_sim_train367457_IRall"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"pp jet-jet Pb-Pb anchor", "Pb-Pb gen.purp. 00-10%", "Pb-Pb gen.purp. 50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };

// const TString analysisWorkflow[nDatasets] = {"track-efficiency_id29229",
//                                              "track-efficiency_central0010",
//                                              "track-efficiency_central5070"
//                                           };
// const TString wagonId[nDatasets] = {"",
//                                     "",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {true, false, false};
// const std::string histDatasetComparisonStructure = "";
// const bool datasetsAreSubsetsofId0 = false;


//////// -------- test cleanUpCollHists ////////
TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
const TString* texDatasetsComparisonType = new TString("simType");
const TString* texDatasetsComparisonCommonDenominator = new TString("");
const int nDatasets = 2;
const TString Datasets[nDatasets] = {"jetjet_ppAnchoredPbPb_5TeV_localNewCleanUpCollHists", "ppRefGenPurposeMC_5TeV_localNewCleanUpCollHists"};
// const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
const TString DatasetsNames[nDatasets] = {"pp jet-jet Pb-Pb anchor", "pp ref MB MC"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
                                        };

const TString analysisWorkflow[nDatasets] = {"track-efficiency",
                                             "track-efficiency"
                                          };
const TString wagonId[nDatasets] = {"",
                                    ""
                                    };
const bool isDatasetWeighted[nDatasets] = {true, false};
const std::string histDatasetComparisonStructure = "";
const bool datasetsAreSubsetsofId0 = false;
