// This is a template. To use the TrackMcQC.C, rename this file to TrackMcQC_inputs.h and edit it how you want.

//////// -------- tracking systematics - track sel variation  ////////
TString* texCollisionDataInfo = new TString("PYTHIA MC #sqrt{#it{s}_{(NN)}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Pb-Pb Angantyr 0-10%");
const TString* texDatasetsComparisonType = new TString("Pb-Pb Angantyr 50-80%");
const TString* texDatasetsComparisonCommonDenominator = new TString("");

const int nDatasets = 6;
const TString DatasetsMC[nDatasets] = {"LHC24g3fix_medium_MC_wrongconfigData_train676663",
                                      "LHC24g3fix_medium_MC_wrongconfigData_train676663",
                                      "LHC24g3fix_medium_MC_wrongconfigData_train676663",
                                      "LHC24g3fix_medium_MC_wrongconfigData_train676663",
                                      "LHC24g3fix_medium_MC_wrongconfigData_train676663",
                                      "LHC24g3fix_medium_MC_wrongconfigData_train676663"};
                                      
const TString DatasetsData[nDatasets] = {"LHC23zzh_pass4_Data_train675029",
                                      "LHC23zzh_pass4_Data_train675029",
                                      "LHC23zzh_pass4_Data_train675029",
                                      "LHC23zzh_pass4_Data_train675029",
                                      "LHC23zzh_pass4_Data_train675029",
                                      "LHC23zzh_pass4_Data_train675029"};

const TString DatasetsNames[nDatasets] = {"nominal", "tightXRows120-5/pt", "tightXRows80", "tightXRows/Findable", "chi2TPC", "chi2ITS"};

TFile* file_O2Analysis_list_MC[nDatasets] = {new TFile("Datasets/"+DatasetsMC[0]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsMC[1]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsMC[2]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsMC[3]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsMC[4]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsMC[5]+"/AnalysisResults.root")
                                          // new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
                                        };
TFile* file_O2Analysis_list_Data[nDatasets] = {new TFile("Datasets/"+DatasetsData[0]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsData[1]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsData[2]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsData[3]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsData[4]+"/AnalysisResults.root"),
                                          new TFile("Datasets/"+DatasetsData[5]+"/AnalysisResults.root")
                                          // new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
                                        };

// IDs: 
// cent0010 data -> id53962
// cent5080 data -> id53963
// cent0010 MC -> id53962 but because wrong train config (ran data instead of MC settings)
// cent5080 MC -> id53963 but because wrong train config (ran data instead of MC settings)
// const TString analysisWorkflowsMC[nDatasets] = {"track-efficiency_id53962",
//                                               "track-efficiency_crossedRows120times5overPt_id53962",
//                                               "track-efficiency_crossedRows080_id53962",
//                                               "track-efficiency_crossedRowsOverFindable_id53962",
//                                               "track-efficiency_chi2TPC_id53962",
//                                               "track-efficiency_chi2ITS_id53962"
//                                           };
// const TString analysisWorkflowsData[nDatasets] = {"track-efficiency_id53962",
//                                               "track-efficiency_crossedRows120times5overPt_id53962",
//                                               "track-efficiency_crossedRows080_id53962",
//                                               "track-efficiency_crossedRowsOverFindable_id53962",
//                                               "track-efficiency_chi2TPC_id53962",
//                                               "track-efficiency_chi2ITS_id53962"
//                                           };

const TString analysisWorkflowsMC[nDatasets] = {"track-efficiency_id53963",
                                              "track-efficiency_crossedRows120times5overPt_id53963",
                                              "track-efficiency_crossedRows080_id53963",
                                              "track-efficiency_crossedRowsOverFindable_id53963",
                                              "track-efficiency_chi2TPC_id53963",
                                              "track-efficiency_chi2ITS_id53963"
                                          };
const TString analysisWorkflowsData[nDatasets] = {"track-efficiency_id53963",
                                              "track-efficiency_crossedRows120times5overPt_id53963",
                                              "track-efficiency_crossedRows080_id53963",
                                              "track-efficiency_crossedRowsOverFindable_id53963",
                                              "track-efficiency_chi2TPC_id53963",
                                              "track-efficiency_chi2ITS_id53963"
                                          };
                                          
// const bool isDatasetWeightedMC[nDatasets] = {false, false, false, false, false, false};
// const bool isDatasetWeightedData[nDatasets] = {false, false, false, false, false, false};

const std::string histDatasetComparisonStructure = "";
const bool datasetsAreSubsetsofId0 = true;

