// This is a template. To use the OccupancyQC.C, rename this file to OccupancyQC_settings.h and edit it how you want.

// // // Bin edge control
// // float GLOBAL_epsilon = 0.00001;

// // Analysis settings
// const int nJetFinderQaType = 3;
// const TString jetFinderQaHistType[nJetFinderQaType] = {"", "_rhoareasubtracted", "_eventwiseconstituentsubtracted"};
// const int nJetType = 3;
// const TString jetType[nJetType] = {"charged", "neutral", "full"};
// const int nJetLevel = 3;
// const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
// // const int nRadius = 3;
// // const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
// // float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// // const float areaDisplayMax[nRadius] = {0.5, 1, 1.5};
// // const int nRadius = 1;
// // const TString RadiusLegend[nRadius] = {"R = 0.4"};
// // float arrayRadius[nRadius] = {0.4};
// // const float areaDisplayMax[nRadius] = {0.5};
// const int nRadius = 9;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
// double arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
// const float areaDisplayMax[nRadius] = {0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.5};

// // Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
// const int iJetType = 0;
// const int iJetLevel = 0;

// // Choice of jet QA type (uncorrected jets, background corrected jet (rho area version), background corrected jet (rho area version) with area cut)
// const int iJetFinderQaType = 0;

// // Choice of Random Cone method:
// const TString methodRandomConeHistName = "withoutleadingjet"; 
// // hist names:                 {"",                   "withoutleadingjet", "randomtrackdirection"};
// // correspond to               {"Random Cones (RC)",  "RC w/o leadJet",    "RC rand(#eta,#phi)"};
// // Default window for random cone:
// std::array<std::array<float, 2>, 2> drawnWindowRCdefault = {{{-30, 60}, {5E-7, 20}}}; // {{xmin, xmax}, {ymin, ymax}}



// // const int nCentralityBins = 6;
// // const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 50, 90};
// // const int nCentralityBins = 3;
// // const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 50, 80};
// const int nCentralityBins = 4;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 30, 50, 70};
// // const int nCentralityBins = 7;
// // const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 60, 70, 80};

 
// const int nTracksBins = 11;
// const float arrayNTracksBinning[nTracksBins+1] = {0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 3000};

// // TFile* file_AliAnalysis = new TFile("../AnalysisResults_Run2_merged_Jaime.root");
// TFile* file_AliAnalysis;




// //////// -------- local check PbPb sel8FullPbPb - all centralities-------- ////////
// TString* texCollisionDataInfo = new TString("PbPb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_sel8FullPbPb", "LHC23zzh_sel8", "LHC24g3_MC_PbPb_WenhuiLocal_selMCFullPbPb", "LHC24g3_MC_PbPb_WenhuiLocal_selMC"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"Data sel8FullPbPb", "Data sel8", "MC selMCFullPbPb", "MC selMC"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"), 
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"), 
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"), 
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-occupancy-qa", 
//                                             "jet-occupancy-qa", 
//                                             "occupancy-qa", 
//                                             "occupancy-qa"
//                                             };

// const TString wagonId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const std::string histDatasetComparisonStructure = "pairs";

// //////// -------- local check PbPb sel8FullPbPb - 0070% centralities-------- ////////
// TString* texCollisionDataInfo = new TString("PbPb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh 00-70%");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_sel8_0070Cent", "LHC23zzh_sel8FullPbPb_0070Cent", "LHC24g3_MC_PbPb_WenhuiLocal_selMCFullPbPb_0070Cent", "LHC24g3_MC_PbPb_WenhuiLocal_selMC_0070Cent"}; //DatasetFiles
// const TString DatasetsNames[nDatasets] = {"Data sel8", "Data sel8FullPbPb", "MC selMCFullPbPb", "MC selMC"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"), 
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"), 
//                                           new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"), 
//                                           new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                           };
// const TString analysisWorkflow[nDatasets] = {"occupancy-qa", 
//                                             "occupancy-qa"
//                                             }; 

// const TString wagonId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};
// const std::string histDatasetComparisonStructure = "pairs";


// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_newSel8FullPbPb_train346468_lowCent","LHC23zzh_newSel8FullPbPb_train346468_highCent"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"00-10%","50-70%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_id24141",
//                                           "jet-finder-charged-qa_id24142"
//                                           };

// const TString wagonId[nDatasets] = {"_id24141",
//                                     "_id24142"
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false,false};
// const std::string histDatasetComparisonStructure = "";

//////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
const TString* texDatasetsComparisonType = new TString("centrality");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
const int nDatasets = 6;
const TString Datasets[nDatasets] = {"LHC23zzh_newSel8FullPbPb_train346468_lowCent", "LHC24g3_newSelMCFullPbPb_train346481_lowCent","LHC23zzh_newSel8FullPbPb_train346468_highCent", "LHC24g3_newSelMCFullPbPb_train346481_highCent", "LHC23zzh_newSel8FullPbPb_train346468_inclCent", "LHC24g3_newSelMCFullPbPb_train346481_inclCent"};
// const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
const TString DatasetsNames[nDatasets] = {"00-10% data","00-10% MC","50-70% data","50-70% MC", "allCent data", "allCent MC"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
                                        };
const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_id24141",
                                          "jet-finder-charged-qa_id24234",
                                          "jet-finder-charged-qa_id24142",
                                          "jet-finder-charged-qa_cent5070_id24234",
                                          "jet-finder-charged-qa_id9300",
                                          "jet-finder-charged-qa_id24422"
                                          };

const TString wagonId[nDatasets] = {"_id24141",
                                    "_id24234",
                                    "_id24142",
                                    "_id24234",
                                    "_id9300",
                                    "_id24422"
                                    };
const bool isDatasetWeighted[nDatasets] = {false,false,false,false,false,false};
const std::string histDatasetComparisonStructure = "twoByTwoDatasetPairs";


// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("allCent");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_newSel8FullPbPb_train346468_lowCent"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"allCent"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_id9300"
//                                           };

// const TString wagonId[nDatasets] = {"_id9300"
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDatasetComparisonStructure = "";

// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centrality");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_newSel8FullPbPb_train346468_lowCent","LHC23zzh_newSel8FullPbPb_train346468_highCent"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"allCent","allCent"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_id9300",
//                                           "jet-finder-charged-qa_id24142"
//                                           };

// const TString wagonId[nDatasets] = {"_id9300",
//                                     "_id24142"
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false,false};
// const std::string histDatasetComparisonStructure = "";

// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_newSel8FullPbPb_train346468_lowCent","LHC24g3_newSelMCFullPbPb_train346481_inclCent"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"Data","MC"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_id9300",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString wagonId[nDatasets] = {"_id9300",
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false,false};
// const std::string histDatasetComparisonStructure = "";


// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts - 1000 occupancy cut with new event sel -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_newSel8FullPbPb_train348332_IgorVersion"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"Data Igor"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {
//                                           "jet-finder-charged-qa"
//                                           };

// const TString wagonId[nDatasets] = {
//                                     ""
//                                     };
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDatasetComparisonStructure = "";