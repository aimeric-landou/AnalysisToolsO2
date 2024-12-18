
// // Bin edge control
// float GLOBAL_epsilon = 0.00001;

// Analysis settings
const int nJetFinderQaType = 3;
const TString jetFinderQaHistType[nJetFinderQaType] = {"", "_rhoareasubtracted", "_eventwiseconstituentsubtracted"};
const int nJetType = 3;
const TString jetType[nJetType] = {"charged", "neutral", "full"};
const int nJetLevel = 3;
const TString jetLevel[nJetLevel] = {"data", "mcd", "mcp"};
// const int nRadius = 3;
// const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.4", "R = 0.6"};
// float arrayRadius[nRadius] = {0.2, 0.4, 0.6};
// const float areaDisplayMax[nRadius] = {0.5, 1, 1.5};
// const int nRadius = 1;
// const TString RadiusLegend[nRadius] = {"R = 0.4"};
// float arrayRadius[nRadius] = {0.4};
// const float areaDisplayMax[nRadius] = {0.5};
const int nRadius = 9;
const TString RadiusLegend[nRadius] = {"R = 0.2", "R = 0.25", "R = 0.3", "R = 0.35", "R = 0.4", "R = 0.45", "R = 0.5", "R = 0.55", "R = 0.6"};
double arrayRadius[nRadius] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
const float areaDisplayMax[nRadius] = {0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.5};

// Choice of jet type (charged, neutral, full) and level (data, detector level, particle level)
const int iJetType = 0;
const int iJetLevel = 0;

// Choice of jet QA type (uncorrected jets, background corrected jet (rho area version), background corrected jet (rho area version) with area cut)
const int iJetFinderQaType = 0;

// Choice of Random Cone method:
const TString methodHistName = "withoutleadingjet"; 
// hist names:                 {"",                   "withoutleadingjet", "randomtrackdirection"};
// correspond to               {"Random Cones (RC)",  "RC w/o leadJet",    "RC rand(#eta,#phi)"};
// Default window for random cone:
std::array<std::array<float, 2>, 2> drawnWindowRCdefault = {{{-30, 60}, {5E-7, 20}}}; // {{xmin, xmax}, {ymin, ymax}}



// const int nCentralityBins = 6;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 50, 90};
// const int nCentralityBins = 3;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 50, 80};
const int nCentralityBins = 4;
const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 30, 50, 70};
// const int nCentralityBins = 7;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 40, 60, 70, 80};

 
const int nTracksBins = 11;
const float arrayNTracksBinning[nTracksBins+1] = {0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 3000};

// TFile* file_AliAnalysis = new TFile("../AnalysisResults_Run2_merged_Jaime.root");
TFile* file_AliAnalysis;


// Options to be set:
// //////// -------- Full cpass0 Analysis -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("cpass0 periods");
// const int nDatasets = 10;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass0","LHC23zzi_cpass0", "LHC23zzf_cpass0", "LHC23zzg_cpass0", "LHC23zx_cpass0", "LHC23zy_cpass0", "LHC23zz_cpass0", "LHC23zza_cpass0", "LHC23zzb_cpass0", "LHC23zze_cpass0"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh","LHC23zzi", "LHC23zzf", "LHC23zzg", "LHC23zx", "LHC23zy", "LHC23zz", "LHC23zza", "LHC23zzb", "LHC23zze"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[8]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[9]+"/AnalysisResults.root")
//                                       };
// TString analysisWorkflow = "jet-finder-"+jetType[iJetType]+"-qa";

// //////// -------- Flat Phi Periods -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("good Phi cpass0 periods");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass0", "LHC23zzf_cpass0", "LHC23zzg_cpass0", "LHC23zzi_cpass0"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh", "LHC23zzf", "LHC23zzg", "LHC23zzi"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Hole Phi Periods -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Period"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("bad Phi cpass0 periods");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zx_cpass0", "LHC23zy_cpass0", "LHC23zz_cpass0", "LHC23zza_cpass0", "LHC23zzb_cpass0", "LHC23zze_cpass0"};
// const TString DatasetsNames[nDatasets] = {"LHC23zx", "LHC23zy", "LHC23zz", "LHC23zza", "LHC23zzb", "LHC23zze"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Cpass comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Cpass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass11", "LHC23zzh_cpass10", "LHC23zzh_cpass8", "LHC23zzh_cpass1", "LHC23zzh_cpass0"};
// const TString DatasetsNames[nDatasets] = {"cpass11", "cpass10", "cpass8", "cpass1", "cpass0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Run Comparison - LHC23zzh -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Run"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh_cpass0");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"LHC23zzh_run544122","LHC23zzh_run544126","LHC23zzh_run544124","LHC23zzh_run544121","LHC23zzh_run544116","LHC23zzh_run544098","LHC23zzh_run544095","LHC23zzh_run544091"};
// const TString DatasetsNames[nDatasets] = {"run 544122","run 544126","run 544124","run 544121","run 544116","run 544098","run 544095","run 544091"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
//                                       };
// TString analysisWorkflow = "jet-finder-"+jetType[iJetType]+"-qa";

// //////// -------- Track Eta cut comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("#eta_{track} cut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh cpass 11");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"LHC23zzh_eta0.9", "LHC23zzh_eta0.8", "LHC23zzh_eta0.7", "LHC23zzh_eta0.6", "LHC23zzh_eta0.5"};
// const TString DatasetsNames[nDatasets] = {"-0.9 < #eta_{Tracks} < 0.9", "-0.8 < #eta_{Tracks} < 0.8", "-0.7 < #eta_{Tracks} < 0.7", "-0.6 < #eta_{Tracks} < 0.6", "-0.5 < #eta_{Tracks} < 0.5"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Jet Radius comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh cpass 11");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_cpass11"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- TRD tracks check -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Track Sel check small - Anton First-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Track Sel check small - Global First-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_GlobalDefault", "LHC23zzh_AntonTrackSel"};
// const TString DatasetsNames[nDatasets] = {"GlobalTrackSel", "MMTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Track Sel check single -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Track Sel check Detailed-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 7;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel", "LHC23zzh_AntonInAcceptanceTrackSel", "LHC23zzh_intermTrackSel3", "LHC23zzh_intermTrackSel2", "LHC23zzh_intermTrackSel1bis", "LHC23zzh_intermTrackSel1", "LHC23zzh_GlobalDefault"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel", "AntonInAcceptance", "intermTrackSel3", "intermTrackSel2", "intermTrackSel1bis", "intermTrackSel1", "GlobalTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Area float vs double check -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("Area double vs float");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_doubleArea", "LHC23zzh_floatArea"};
// const TString DatasetsNames[nDatasets] = {"areaDouble", "areaFloat"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow = "jet-finder-charged-qa";


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





// //////// -------- Check trackSel single -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TrackSel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_AntonTrackSel"};
// const TString DatasetsNames[nDatasets] = {"MMTrackSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";



// //////// -------- Leading Track Cut -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_LeadCut0", "LHC23zzh_LeadCut1", "LHC23zzh_LeadCut2", "LHC23zzh_LeadCut3", "LHC23zzh_LeadCut4", "LHC23zzh_LeadCut5"};
// const TString DatasetsNames[nDatasets] = {"LeadCut 0 GeV", "LeadCut 1 GeV", "LeadCut 2 GeV", "LeadCut 3 GeV", "LeadCut 4 GeV", "LeadCut 5 GeV"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";



// //////// -------- Leading Track Cut - single -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_LeadCut5"};
// const TString DatasetsNames[nDatasets] = {"LeadCut 5 GeV"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";



// //////// -------- Apasses comparison - apass 2 first -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Apass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_run544122", "LHC23zzh_apass1_run544095"};
// const TString DatasetsNames[nDatasets] = {"apass2 run 544122", "apass1 run 544095"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Apasses comparison - apass 1 first -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Apass");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass1_run544095", "LHC23zzh_apass2_run544122"};
// const TString DatasetsNames[nDatasets] = {"apass1 run 544095", "apass2 run 544122"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";



// //////// -------- Leading Track Max Cut -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackUpperCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 1 run 544095");
// const int nDatasets = 7;
// const TString Datasets[nDatasets] = {"LHC23zzh_LeadCutNone", "LHC23zzh_LeadCut10", "LHC23zzh_LeadCut5", "LHC23zzh_LeadCut4", "LHC23zzh_LeadCut3", "LHC23zzh_LeadCut2", "LHC23zzh_LeadCut1"};
// const TString DatasetsNames[nDatasets] = {"LeadTrack < 999 GeV", "LeadTrack < 10 GeV", "LeadTrack < 5 GeV", "LeadTrack < 4 GeV", "LeadTrack < 3 GeV", "LeadTrack < 2 GeV", "LeadTrack < 1 GeV" };
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Leading Track Max Cut - single -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544095");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_LeadCut5"};
// const TString DatasetsNames[nDatasets] = {"LeadCut 5 GeV"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Test jet leading pt vs Rho -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackUpperCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 1 run 544095");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_LeadCutNone"};
// const TString DatasetsNames[nDatasets] = {""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Rho Fix (leading jets removed) vs pre Rho fix -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackUpperCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_postRhoFix", "LHC23zzh_preRhoFix"};
// const TString DatasetsNames[nDatasets] = {"postRhoFix", "preRhoFix"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- MC vs Data -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("LeadTrackUpperCut");
// const TString* texDatasetsComparisonCommonDenominator = new TString("anchor apass2");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_MC_LHC23k6d", "LHC23zzh_Data_LHC23zzhRun544122"};
// const TString DatasetsNames[nDatasets] = {"MC_LHC23k6d", "Data_LHC23zzhRun544122"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";



// //////// -------- JetValidationFramework -------- ////////
// TString* texCollisionDataInfo = new TString("Run 2, pp 13 TeV, LHC18b, real"); 
// const TString* texDatasetsComparisonType = new TString("");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"AliPhysics", "O2Physics"};
// const TString DatasetsNames[nDatasets] = {"AliPhysics", "O2Physics"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                           new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- RandomCone Fix (leading jets removed) vs pre RandomCone fix -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("RandomConeLeadFix");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_preRandomConeFix_withAreaCut", "LHC23zzh_apass2_postRandomConeFix_noAreaCut", "LHC23zzh_apass2_postRandomConeFix_withAreaCut"};
// const TString DatasetsNames[nDatasets] = {"preRCfix_withAreaCut", "postRCfix_noAreaCut", "postRCfix_withAreaCut"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- Timeframe Border Cut comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("TimeframeBorder");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_postTimeFrameBorderCut", "LHC23zzh_apass2_preTimeFrameBorderCut"};
// const TString DatasetsNames[nDatasets] = {"postTimeBorderCut", "preTimeBorderCut"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";


// //////// -------- Run 2 vs Run 3 comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Run2Run3");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_Area04_leadConstCut5", "LHC23zzh_apass2_Area04_leadConstCut0", "LHC23zzh_apass2_Area06_leadConstCut0"};
// const TString DatasetsNames[nDatasets] = {"Area04_leadConstCut5", "Area04_leadConstCut0", "Area06_leadConstCut0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-finder-charged-qa";

// //////// -------- stdRho centrality comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("CentralityWindows");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122 - stdRho");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_stdRho_withCentralityCutWindows", "LHC23zzh_apass2_stdRho_withCentralityCutWindows", "LHC23zzh_apass2_stdRho_withCentralityCutWindows", "LHC23zzh_apass2_stdRho_withCentralityCutWindows", "LHC23zzh_apass2_stdRho_withCentralityCutWindows", "LHC23zzh_apass2_stdRho_withCentralityCutWindows"};
// const TString DatasetsNames[nDatasets] = {"central0010", "central1020", "central2030", "central3040", "central4050", "central5090"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_1020",
//                                           "jet-finder-charged-qa_central_2030",
//                                           "jet-finder-charged-qa_central_3040",
//                                           "jet-finder-charged-qa_central_4050",
//                                           "jet-finder-charged-qa_central_5090",
//                                           };


// //////// -------- sparseRho centrality comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("CentralityWindows");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122 - sparseRho");
// const int nDatasets = 6;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sparseRho_withCentralityCutWindows", "LHC23zzh_apass2_sparseRho_withCentralityCutWindows", "LHC23zzh_apass2_sparseRho_withCentralityCutWindows", "LHC23zzh_apass2_sparseRho_withCentralityCutWindows", "LHC23zzh_apass2_sparseRho_withCentralityCutWindows", "LHC23zzh_apass2_sparseRho_withCentralityCutWindows"};
// const TString DatasetsNames[nDatasets] = {"central0010", "central1020", "central2030", "central3040", "central4050", "central5090"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_1020",
//                                           "jet-finder-charged-qa_central_2030",
//                                           "jet-finder-charged-qa_central_3040",
//                                           "jet-finder-charged-qa_central_4050",
//                                           "jet-finder-charged-qa_central_5090",
//                                           };



// //////// -------- Sparse vs Standard Rho comparison - central 0-10% events -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("sparseVsStd");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sparseRho_withCentralityCutWindows", "LHC23zzh_apass2_stdRho_withCentralityCutWindows"};
// const TString DatasetsNames[nDatasets] = {"sparseRho", "stdRho"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010"
//                                           };


// //////// -------- Standard Rho - central 0-10% events -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("UpperVsLowerArm");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_inclusive", "LHC23zzh_apass2_lowerArm"};
// const TString DatasetsNames[nDatasets] = {"inclusive", "lowerArm"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa"
//                                           };




// //////// -------- Run Comparison - LHC23zzh apass 2 - central 0-10% -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Run"); 
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh_apass2 - 0-10% - areaCut 0.4");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_run544124","LHC23zzh_apass2_run544123","LHC23zzh_apass2_run544122","LHC23zzh_apass2_run544121","LHC23zzh_apass2_run544116","LHC23zzh_apass2_run544098","LHC23zzh_apass2_run544095","LHC23zzh_apass2_run544091"};
// const TString DatasetsNames[nDatasets] = {"run 544124","run 544123","run 544122","run 544121","run 544116","run 544098","run 544095","run 544091"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                       new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root")
//                                       };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_0010"
//                                           };


// //////// -------- ITSROF and BC and GoodZvtxFT0PV cuts vs before -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("newEvtCuts");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8ITSROFuBCuZvtxCut", "LHC23zzh_apass2_sel8ITSROFuBC", "LHC23zzh_apass2_sel8ITSROF", "LHC23zzh_apass2_sel8", "LHC23zzh_apass2_noSel"};
// const TString DatasetsNames[nDatasets] = {"sel8ITSROFuBCuZvtxCut", "sel8ITSROFuBC", "sel8ITSROF", "sel8", "noSel"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// //////// -------- sel8 vs sel8Full - ITSROF and BC and GoodZvtxFT0PV cuts vs before -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("newEvtCuts");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8", "LHC23zzh_apass2_sel8Full"};
// const TString DatasetsNames[nDatasets] = {"sel8", "sel8Full"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// //////// -------- sel8Full - apass2 - centrality comp -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centralityWindow");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8Full", "LHC23zzh_apass2_sel8Full"};
// const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
//                                           "jet-finder-charged-qa_central_5090"
//                                           };

// //////// -------- apass3 vs apass 2 with sel8Full -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("apass3v2");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh run 544122");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8Full", "LHC23zzh_apass3_sel8Full"};
// const TString DatasetsNames[nDatasets] = {"apass2", "apass3"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };



// //////// -------- Run 2 vs Run 3 comparison -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV or 5.02 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Run2Run3");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass3_sel8Full", "LHC18q"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh", "LHC18q"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// //////// -------- sel8Full vs GoodZvtxFT0PV vs sel8 -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("newEvtCuts");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass 2 run 544122");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_sel8", "LHC23zzh_apass2_sel8uZvtxCut", "LHC23zzh_apass2_sel8uBC", "LHC23zzh_apass2_sel8Full"};
// const TString DatasetsNames[nDatasets] = {"sel8", "sel8uZvtxCut", "sel8uBC", "sel8Full"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };


// //////// -------- sel8Full - apass 3 - centrality comp -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("centralityWindow");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass3_sel8Full", "LHC23zzh_apass3_sel8Full"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"", ""};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };


// //////// -------- sel8Full - apass 3 - centrality comp -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("efficiencyDiscardingUpdate");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass2_preEffUpdate", "LHC23zzh_apass2_postEffUpdate_eps=1.0", "LHC23zzh_apass2_postEffUpdate_eps=0.5"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"preEffUpdate", "postEffUpdate eps=1.0", "postEffUpdate eps=0.5"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };


// //////// -------- sel8Full - apass 3 - centrality comp -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_MC_LHC23k6d"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23k6d"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
                                        


// //////// -------- sel8Full - apass 3 - centrality comp -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24d2");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24d2_sel8_and_sel8Full", "LHC24d2_sel8_and_sel8Full"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"sel8Full", "sel8"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa_sel8"
//                                           };

// const TString trainId = "";


// //////// -------- LHC23zzh - apass 4 -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_areaLeadingCuts_centralityWindowsAvailable/run544116", "LHC23zzh_apass4_areaLeadingCuts_centralityWindowsAvailable/run544123"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"run544116", "run544123"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false};




// //////// -------- pp MC sel8MC vs sel8Full apass6 anchorr? -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV");
// const TString* texDatasetsComparisonType = new TString("MC");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC24b1b");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC24b1b_sel8MC_train239181", "LHC24b1b_sel8Full_train239409"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"sel8MC", "sel8Full"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";



// //////// -------- LHC23zzh - apass 4 with area and leadingtrackpt cuts-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train254040/merged"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh_apass4_areaLeadingCuts"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDrawColorsOption = "";


// //////// -------- LHC23zzh - apass 4 withOUT area and leadingtrackpt cuts-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train254040/merged"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh_apass4_NoAreaNorLeadingCuts"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_noAreaLeadingJetCut"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDrawColorsOption = "";



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
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
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
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {true, true};

// // jaime testJetAnchoredGapGen-20240801-155530



// //////// -------- LHC23zzh - apass 4 withOUT area and leadingtrackpt cuts VS LHC22o pp-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_areaLeadingCuts_AND_NoAreaNorLeadingCuts", "LHC22o_pass6_train238827"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh apass4", "LHC22o pass6"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_noAreaLeadingJetCut",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false};


// //////// -------- LHC23zzh - apass 4 occupancy analysis sel8FullPbPb-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8FullPbPb");
// const int nDatasets = 5;
// const TString Datasets[nDatasets] = {"sel8FullPbPb/LHC23zzh_apass4_occupancy01000", "sel8FullPbPb/LHC23zzh_apass4_occupancy02000", "sel8FullPbPb/LHC23zzh_apass4_occupancy06000", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-2k occupancy", "2k-6k occupancy", "6k-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false};
// const std::string histDrawColorsOption = "";


// //////// -------- LHC23zzh - apass 4 occupancy analysis sel8Full-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8Full");
// const int nDatasets = 4;
// const TString Datasets[nDatasets] = {"sel8Full/LHC23zzh_apass4_occupancy01000", "sel8Full/LHC23zzh_apass4_occupancy05000", "sel8Full/LHC23zzh_apass4_occupancy10000", "sel8Full/LHC23zzh_apass4_occupancy20000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k occupancy", "1k-5k occupancy", "5k-10k occupancy", "10-20k occupancy"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false};




// //////// -------- LHC23zzh - apass 4 occupancy analysis - sel8Full vs sel8FullPbPb -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("sel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"sel8Full/LHC23zzh_apass4_occupancy01000", "sel8FullPbPb/LHC23zzh_apass4_occupancy01000", "sel8Full/LHC23zzh_apass4_occupancy05000", "sel8FullPbPb/LHC23zzh_apass4_occupancy05000", "sel8Full/LHC23zzh_apass4_occupancy10000", "sel8FullPbPb/LHC23zzh_apass4_occupancy10000", "sel8Full/LHC23zzh_apass4_occupancy20000", "sel8FullPbPb/LHC23zzh_apass4_occupancy20000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k sel8Full", "0k-1k sel8FullPbPb", "1k-5k sel8Full", "1k-5k sel8FullPbPb", "5k-10k sel8Full", "5k-10k sel8FullPbPb", "10-20k sel8Full", "10-20k sel8FullPbPb"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false, false, false, false};
// // const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true, true, true, true, true, true};

// const std::string histDrawColorsOption = "colorPairs";


// //////// -------- LHC23zzh - apass 4 train train253451 vs train247272 - bkgFluct unexpected change-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train253451", "LHC23zzh_apass4_train245750"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"apass4", "apass3"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false};



// //////// -------- Run 3 data apass4 - run by run -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 3;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train254040/run544116", "LHC23zzh_apass4_train254040/run544121", "LHC23zzh_apass4_train254040/run544123"};
// const TString DatasetsNames[nDatasets] = {"run 544116", "run 544121", "run 544123"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                             "jet-finder-charged-qa",
//                                             "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false};
// const std::string histDrawColorsOption = "";

// //////// -------- Run 3 data apass4 - merged runs -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV"); 
// const TString* texDatasetsComparisonType = new TString("run");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4_train254040/merged"};
// const TString DatasetsNames[nDatasets] = {"mergedRuns"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDrawColorsOption = "";



// //////// -------- LHC23zzh - apass 4 occupancy vs run analysis -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("sel");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4");
// const int nDatasets = 8;
// const TString Datasets[nDatasets] = {"run544123/LHC23zzh_apass4_occupancy01000", "run544121/LHC23zzh_apass4_occupancy01000", "run544116/LHC23zzh_apass4_occupancy01000", "run544095/LHC23zzh_apass4_occupancy01000", "run544123/LHC23zzh_apass4_occupancy10000", "run544121/LHC23zzh_apass4_occupancy10000", "run544116/LHC23zzh_apass4_occupancy10000", "run544095/LHC23zzh_apass4_occupancy10000"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"0k-1k run544123", "0k-1k run544121", "0k-1k run544116", "0k-1k run544095", "5k-10k run544123", "5k-10k run544121", "5k-10k run544116", "5k-10k run544095"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[2]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[3]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[4]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[5]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[6]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[7]+"/AnalysisResults.root"),
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa",
//                                           "jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false, false, false, false};
// // const bool trackHistsObsoleteVersion[nDatasets] = {true, true, true, true, true, true, true, true};

// const std::string histDrawColorsOption = "";





// //////// -------- LHC23zzh - apass 4 - random cone methods -------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh apass4");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh_apass4"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDrawColorsOption = "";


// //////// -------- LHC23zzh - apass 4 withOUT area and leadingtrackpt cuts-------- ////////
// TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
// const TString* texDatasetsComparisonType = new TString("Data");
// const TString* texDatasetsComparisonCommonDenominator = new TString("");
// const int nDatasets = 1;
// const TString Datasets[nDatasets] = {"LHC23zzh_apass4"};
// // const TString DatasetsNames[nDatasets] = {"0-10%", "50-90%"};
// const TString DatasetsNames[nDatasets] = {"LHC23zzh_apass4_NoAreaNorLeadingCuts"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root")
//                                         };
// // const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_central_0010",
// //                                           "jet-finder-charged-qa_central_5090"
// //                                           };
// const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa_noAreaLeadingJetCut"
//                                           };

// const TString trainId = "";
// const bool isDatasetWeighted[nDatasets] = {false};
// const std::string histDrawColorsOption = "";

//////// -------- LHC23zzh - apass 4 occupancy analysis sel8FullPbPb-------- ////////
TString* texCollisionDataInfo = new TString("Pb-Pb #sqrt{#it{s}} = 5.36 TeV");
const TString* texDatasetsComparisonType = new TString("Occupancy analysis");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC23zzh pass4 - sel8FullPbPb");
const int nDatasets = 1;
const TString Datasets[nDatasets] = {"LHC24g3"};
const TString DatasetsNames[nDatasets] = {"test_bonjour"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile(Datasets[0]+"/AnalysisResults.root")
                                        };
const TString analysisWorkflow[nDatasets] = {"jet-finder-charged-qa"
                                          };

const TString trainId = "";
const bool isDatasetWeighted[nDatasets] = {false, false, false, false, false};
const std::string histDrawColorsOption = "";