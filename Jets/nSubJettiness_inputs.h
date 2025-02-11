// To avoid your local _inputs.h files to be replaced by the template in the git repository, you should ask git to ignore them with git update-index --assume-unchanged (for more details see https://stackoverflow.com/questions/3319479/can-i-git-commit-a-file-and-ignore-its-content-changes)

// TFile* file_AliAnalysis = new TFile("../AnalysisResults_Run2_merged_Jaime.root");


// //////// -------- Inclusive vs D0 -------- ////////
// TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
// const TString* texDatasetsComparisonType = new TString("Inclusive/HF");
// const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o apass6");
// const int nDatasets = 2;
// const TString Datasets[nDatasets] = {"Inclusive/JE_LHC22o_pass6_small_C_R6_18_F_R6_8", "HF_D0/JE_HF_LHC22_pass6_medium_2P3PDstar_D0C_R4_4"};
// const TString DatasetsNames[nDatasets] = {"inclusive", "D0"};
// TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
//                                         new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
//                                         };
// const TString analysisWorkflow = "jet-nsubjettiness";


//////// -------- Inclusive - before/after cleaning of workflow -------- ////////
TString* texCollisionDataInfo = new TString("pp #sqrt{#it{s}} = 13.6 TeV"); 
const TString* texDatasetsComparisonType = new TString("Inclusive/HF");
const TString* texDatasetsComparisonCommonDenominator = new TString("LHC22o apass6");
const int nDatasets = 2;
const TString Datasets[nDatasets] = {"Inclusive/JE_LHC22o_pass6_small_C_R6_18_F_R6_8", "Inclusive/cleaningOfWorkflow"};
const TString DatasetsNames[nDatasets] = {"beforeCleaning", "afterCleaning"};
TFile* file_O2Analysis_list[nDatasets] = {new TFile("Datasets/"+Datasets[0]+"/AnalysisResults.root"),
                                        new TFile("Datasets/"+Datasets[1]+"/AnalysisResults.root")
                                        };
const TString analysisWorkflow = "jet-nsubjettiness";
