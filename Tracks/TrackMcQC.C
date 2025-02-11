#include "TStyle.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRatioPlot.h"
#include "TLegend.h"
#include "TH1.h"
#include "iostream"
// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./TrackMcQC_settings.h"
#include "./TrackMcQC_inputs.h"
#include "../Settings/AxisTitles.h"
#include "../Settings/GlobalSettings.h"
#include "../Utilities/AnalysisUtilities.h"
#include "../Utilities/HistogramUtilities.h"
#include "../Utilities/HistogramPlotting.h"
#include "../Utilities/AnalysisUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramUtilities.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well
#include "../Utilities/HistogramPlotting.C" // bizarre but only including the .h fils doesn't work for the standard 'root macro.C+' method, I need to include the .C as well


#include<array>
#include <iomanip>
#include <sstream>
#include <string.h>
using namespace std;

// Misc utilities
void SetStyle(Bool_t graypalette=kFALSE);
void LoadLibs();

// // Plot Utilities
// TString contextDatasetComp(std::string options);

//////////// QC plot functions
// Dataset comparison

void Draw_Efficiency_Pt_DatasetComparison(float* etaRange, bool useSplit);
void Draw_Efficiency_Eta_DatasetComparison(float* PtRange, bool useSplit);
void Draw_Efficiency_Phi_DatasetComparison(float* ptRange, float* etaRange, bool useSplit);
void Draw_Efficiency_Eta_PtRangeComparison(float* ptRange, int nPtRanges, int iDataset, std::string options);
void Draw_Efficiency_Phi_PtRangeComparison(float* ptRange, int nPtRanges, float* etaRange, int iDataset, std::string options);
void Draw_Efficiency_Phi_DatasetComparison_finerPhi(float* ptRange, float* etaRange);

void Draw_Efficiency_Pt_DatasetComparison_legacy(float* etaRange);
void Draw_Efficiency_Eta_DatasetComparison_legacy(float* PtRange);
void Draw_Efficiency_Phi_DatasetComparison_legacy(float* ptRange, float* etaRange);
void Draw_Efficiency_Eta_PtRangeComparison_legacy(float* ptRange, int nPtRanges, int iDataset, std::string options);
void Draw_Efficiency_Phi_PtRangeComparison_legacy(float* ptRange, int nPtRanges, float* etaRange, int iDataset, std::string options);
void Draw_Efficiency_Phi_DatasetComparison_finerPhi_legacy(float* ptRange, float* etaRange);
void Draw_Efficiency_Pt_ratio_etaNeg_etaPos_DatasetComparison(float* etaRange, bool useSplit);


void Draw_Purity_Pt_DatasetComparison(float* etaRange, bool useSplit);
void Draw_Purity_Eta_DatasetComparison(float* etaRange, bool useSplit);
void Draw_Purity_Phi_DatasetComparison(float* etaRange, bool useSplit);
void Draw_Purity_Pt_ratio_etaNeg_etaPos_DatasetComparison(float* etaRange, bool useSplit);


void Draw_mcPt_DatasetComparison(int particleStatusOption);
void Draw_Vx_DatasetComparison(int particleStatusOption);
void Draw_Vy_DatasetComparison(int particleStatusOption);
void Draw_Vz_DatasetComparison(int particleStatusOption);
void Draw_y_iu_DatasetComparison();



void Draw_Pt_gen_DatasetComparison(float* etaRange, std::string options);
void Draw_Eta_gen_DatasetComparison(float* ptRange, std::string options);
void Draw_Phi_gen_DatasetComparison(float* ptRange, float* etaRange, std::string options);



void Draw_Pt_gen_DatasetComparison_H2CentVersion(std::string options);
void Draw_Eta_gen_DatasetComparison_H2CentVersion(std::string options);
void Draw_Phi_gen_DatasetComparison_H2CentVersion(std::string options);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void TrackMcQC() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();
  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");
  bool useSplit = true; //set to true if you want to see the influence of split
  float etaRange[2] = {-0.9, 0.9};
  //  Draw_Efficiency_Pt_DatasetComparison(etaRange,useSplit);
  float ptRange1[2] = {0.15, 100};
  // Draw_Efficiency_Eta_DatasetComparison(ptRange1,useSplit);
  // Draw_Efficiency_Phi_DatasetComparison(ptRange1, etaRange, useSplit);

  // Draw_Efficiency_Pt_ratio_etaNeg_etaPos_DatasetComparison(etaRange, useSplit);
  // // Draw_Efficiency_Phi_DatasetComparison_finerPhi(ptRange1, etaRange); // only works with very specific datasets created locally

  // Draw_Purity_Pt_DatasetComparison(etaRange, useSplit);
  // Draw_Purity_Eta_DatasetComparison(etaRange, useSplit);
  // Draw_Purity_Phi_DatasetComparison(etaRange, useSplit);
  Draw_Purity_Pt_ratio_etaNeg_etaPos_DatasetComparison(etaRange, useSplit);

  // int nPtRanges = 4;
  // float ptRange[5] = {0, 0.15, 0.5, 1, 100};
  // for(int iDataset = 0; iDataset < nDatasets; iDataset++){
  //   Draw_Efficiency_Eta_PtRangeComparison(ptRange, nPtRanges, iDataset, "");
  //   Draw_Efficiency_Eta_PtRangeComparison(ptRange, nPtRanges, iDataset, "scaled");
  //   Draw_Efficiency_Phi_PtRangeComparison(ptRange, nPtRanges, etaRange, iDataset, "");
  //   Draw_Efficiency_Phi_PtRangeComparison(ptRange, nPtRanges, etaRange, iDataset, "scaled");
  // }


  // FOR SIMULATION CHECKS, don't use in standard mc files
  // int particleStatusOption = 2;
  // Draw_mcPt_DatasetComparison(particleStatusOption);
  // Draw_Vx_DatasetComparison(particleStatusOption);
  // Draw_Vy_DatasetComparison(particleStatusOption);
  // Draw_Vz_DatasetComparison(particleStatusOption);
  // Draw_y_iu_DatasetComparison();

  // float ptRange[2] = {0.150, 100};
  // float etaRange[2] = {-0.9, 0.9};
  // Draw_Pt_gen_DatasetComparison(etaRange, "primaries, ratio");
  // Draw_Eta_gen_DatasetComparison(ptRange, "primaries, ratio");
  // Draw_Phi_gen_DatasetComparison(ptRange, etaRange, "primaries, ratio");
  // Draw_Pt_gen_DatasetComparison(etaRange, "secondaries, ratio");
  // Draw_Eta_gen_DatasetComparison(ptRange, "secondaries, ratio");
  // Draw_Phi_gen_DatasetComparison(ptRange, etaRange, "secondaries, ratio");
  // Draw_Pt_gen_DatasetComparison(etaRange, "primaries,secondaries, ratio");
  // Draw_Eta_gen_DatasetComparison(ptRange, "primaries,secondaries, ratio");
  // Draw_Phi_gen_DatasetComparison(ptRange, etaRange, "primaries,secondaries, ratio");


  // Draw_Pt_gen_DatasetComparison_H2CentVersion("primaries,secondaries, ratio");
  // Draw_Eta_gen_DatasetComparison_H2CentVersion("primaries,secondaries, ratio");
  // Draw_Phi_gen_DatasetComparison_H2CentVersion("primaries,secondaries, ratio");
}
/////////////////////////////////////////////////////
/////////////////// Misc utilities //////////////////
/////////////////////////////////////////////////////

void LoadLibs() {
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libVMC");
  // gSystem->Load("libTree");
  // gSystem->Load("libMinuit");
  // gSystem->Load("libSTEERBase");
  // gSystem->Load("libESD");
  // gSystem->Load("libAOD");
  // gSystem->Load("libANALYSIS");
  // gSystem->Load("libANALYSISalice");
  // gSystem->Load("libCORRFW");
  // gSystem->Load("libPWGTools");
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLineScalePS(1);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}




// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////// Context Utilities /////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TString contextTrackDatasetComp(std::string options){
  TString texcontextDatasetCompAndRadiusAndVarRange;
  //  if (options.find("track") == std::string::npos) {
  texcontextDatasetCompAndRadiusAndVarRange = *texDatasetsComparisonCommonDenominator;
  return texcontextDatasetCompAndRadiusAndVarRange;
}


// TString contextPtRange(float* PtRange){
//   std::stringstream ss;
//   ss << PtRange[0] << " < #it{p}_{T} < " << PtRange[1];
//   TString textContext((TString)ss.str());
//   // TString texDataset(Form("%.0f", PtRange[0])+" < #it{p}_{T} < "+Form("%.0f", PtRange[1]));
//   return textContext;
// }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Efficiency_Pt_DatasetComparison(float* etaRange, bool useSplit) {
  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_mcparticles_ptHigh[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_ptHigh[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split_ptHigh[nDatasets];
  // TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPt_mcparticles[nDatasets];
  TH1D* H1D_trackPt_mcparticles_rebinned[nDatasets];
  TH1D* H1D_trackPt_mcparticles_ptHigh[nDatasets];
  TH1D* H1D_trackPt_mcparticles_ptHigh_rebinned[nDatasets];
  TH1D* H1D_trackPt_assoctracks[nDatasets];
  TH1D* H1D_trackPt_assoctracks_rebinned[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split_rebinned[nDatasets];
  TH1D* H1D_trackPt_assoctracks_ptHigh[nDatasets];
  TH1D* H1D_trackPt_assoctracks_ptHigh_rebinned[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split_ptHigh[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split_ptHigh_rebinned[nDatasets];
  // TH1D* H1D_trackPt_assoctracks_nonSplitSecondary[nDatasets];
  TH1D* H1D_trackPt_efficiency[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPt_efficiency_ptHigh[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected_ptHigh[nDatasets];
  TH1D* H1D_trackPt_efficiency_split[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitNotCorrected_ptHigh[nDatasets];

  // TH1D* H1D_trackPt_efficiency_splitAndSecondaryCorrected[nDatasets];
  TH1D* H1D_trackPt_efficiency_concatenated[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected_concatenated[nDatasets];
  TH1D* H1D_trackPt_efficiency_split_concatenated[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_ptHigh = false;
  bool divideSuccess_split = false;
  bool divideSuccess_split_ptHigh = false;
  bool divideSuccess_split_notCorrected = false;
  bool divideSuccess_split__notCorrectedPtHigh = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Pt_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_mcparticles_ptHigh[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_high_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Pt_DatasetComparison_mcparticles_ptHigh"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks_ptHigh[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelColl_ptHigh"+Datasets[iDataset]);

    if (useSplit == true) {
      H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
      H3D_trackPtEtaPhi_assoctracks_split_ptHigh[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_high_particle_eta_particle_phi_associatedtrack_split_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollSplit_ptHigh"+Datasets[iDataset]);
    }
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon + deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon -deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon );
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon);

    H1D_trackPt_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionX("trackPt_mcparticles"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    H1D_trackPt_mcparticles_ptHigh[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles_ptHigh[iDataset]->ProjectionX("trackPt_mcparticles_ptHigh"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");

    H1D_trackPt_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionX("trackPt_tracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    H1D_trackPt_assoctracks_ptHigh[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_ptHigh[iDataset]->ProjectionX("trackPt_tracks_ptHigh"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");

    if (useSplit == true) {
        H1D_trackPt_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionX("trackPt_assoctracks_split"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
        H1D_trackPt_assoctracks_split_ptHigh[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split_ptHigh[iDataset]->ProjectionX("trackPt_assoctracks_split_ptHigh"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    }
    
    //tweaking the low-pt bins to make them larger close to 10GeV
    std::vector<double> xbinsVectorInitialLow = GetTH1Bins(H1D_trackPt_mcparticles[iDataset]);
    double* xbinsInitialLow = &xbinsVectorInitialLow[0];

    double ptBinsLowNew[500]; //500 to have a good margin
    double lastPt;
    int iBinPtNew = 0;
    int iBinPtInitialHisto = 0;
    int increment;
    while(iBinPtInitialHisto < H1D_trackPt_mcparticles[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialLow[iBinPtInitialHisto];
      if (lastPt < 0.300) {
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_mcparticles[iDataset]->GetXaxis()->GetBinWidth(1) /8;
      }
      ptBinsLowNew[iBinPtNew] = lastPt;
      iBinPtInitialHisto += increment;
      iBinPtNew += 1;
    }
    ptBinsLowNew[iBinPtNew-1] = 10; // replace last set bin edge by 10
    int nBinsLowNew = iBinPtNew - 1;
    if (nBinsLowNew >= 500){
      cout << "NEEDS TO RESERVE MORE MEMORY FOR ptBinsLowNew !!!!!!!!!!!!!!!!!!" << endl;
    }
    H1D_trackPt_mcparticles_rebinned[iDataset] = (TH1D*)H1D_trackPt_mcparticles[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_mcparticles_rebinned"+Datasets[iDataset], ptBinsLowNew);
    H1D_trackPt_assoctracks_rebinned[iDataset] = (TH1D*)H1D_trackPt_assoctracks[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_assoctracks_rebinned"+Datasets[iDataset], ptBinsLowNew);
    if (useSplit == true) {
      H1D_trackPt_assoctracks_split_rebinned[iDataset] = (TH1D*)H1D_trackPt_assoctracks_split[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_assoctracks_split_rebinned"+Datasets[iDataset], ptBinsLowNew);
    }

    //tweaking the high-pt bins to have them less prone to statistical fluctuations
    std::vector<double> xbinsVectorInitialHigh = GetTH1Bins(H1D_trackPt_mcparticles_ptHigh[iDataset]);
    double* xbinsInitialHigh = &xbinsVectorInitialHigh[0];

    double ptBinsHighNew[500]; //500 to have a good margin
    iBinPtNew = 0;
    iBinPtInitialHisto = 0;
    while(iBinPtInitialHisto < H1D_trackPt_mcparticles_ptHigh[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialHigh[iBinPtInitialHisto];
      if (lastPt < 25) {
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_mcparticles_ptHigh[iDataset]->GetXaxis()->GetBinWidth(1) ;
      }
      ptBinsHighNew[iBinPtNew] = lastPt;
      iBinPtInitialHisto += increment;
      iBinPtNew += 1;
    }
    ptBinsHighNew[iBinPtNew-1] = 100; // replace last set bin edge by 10
    int nBinsHighNew = iBinPtNew - 1;
    if (nBinsHighNew >= 500){
      cout << "NEEDS TO RESERVE MORE MEMORY FOR ptBinsHighNew !!!!!!!!!!!!!!!!!!" << endl;
    }
    H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_mcparticles_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_mcparticles_ptHigh_rebinned"+Datasets[iDataset], ptBinsHighNew);
    H1D_trackPt_assoctracks_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_assoctracks_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_assoctracks_ptHigh_rebinned"+Datasets[iDataset], ptBinsHighNew);
    if (useSplit == true) {
      H1D_trackPt_assoctracks_split_ptHigh_rebinned[iDataset] = (TH1D*)H1D_trackPt_assoctracks_split_ptHigh[iDataset]->Rebin(nBinsHighNew, "H1D_trackPt_assoctracks_split_ptHigh_rebinned"+Datasets[iDataset], ptBinsHighNew);
    }
    // getting the efficiency
    H1D_trackPt_efficiency[iDataset] = (TH1D*)H1D_trackPt_mcparticles_rebinned[iDataset]->Clone("trackPt_efficiency"+Datasets[iDataset]);
    H1D_trackPt_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_efficiency[iDataset]->Divide(H1D_trackPt_assoctracks_rebinned[iDataset], H1D_trackPt_mcparticles_rebinned[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    H1D_trackPt_efficiency_ptHigh[iDataset] = (TH1D*)H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset]->Clone("trackPt_efficiency_ptHigh"+Datasets[iDataset]);
    H1D_trackPt_efficiency_ptHigh[iDataset]->Reset("M");
    divideSuccess_ptHigh = H1D_trackPt_efficiency_ptHigh[iDataset]->Divide(H1D_trackPt_assoctracks_ptHigh_rebinned[iDataset], H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    // correcting for split tracks
    if (useSplit == true) {
      //getting the efficiency : for split corrected 
      H1D_trackPt_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPt_mcparticles_rebinned[iDataset]->Clone("trackPt_efficiency_splitCorrected"+Datasets[iDataset]);
      H1D_trackPt_efficiency_splitCorrected[iDataset]->Reset("M");
      divideSuccess_split = H1D_trackPt_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPt_assoctracks_split_rebinned[iDataset], H1D_trackPt_mcparticles_rebinned[iDataset]);
      H1D_trackPt_efficiency_splitCorrected[iDataset]->Scale(-1.);
      H1D_trackPt_efficiency_splitCorrected[iDataset]->Add(H1D_trackPt_efficiency[iDataset],1.);
      H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset] = (TH1D*)H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset]->Clone("trackPt_efficiency_splitCorrected_ptHigh"+Datasets[iDataset]);
      H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset]->Reset("M");
      divideSuccess_split_ptHigh = H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset]->Divide(H1D_trackPt_assoctracks_split_ptHigh_rebinned[iDataset], H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
      H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset]->Scale(-1.);
      H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset]->Add(H1D_trackPt_efficiency_ptHigh[iDataset],1.);
      //getting the efficiency : for split only
      H1D_trackPt_efficiency_split[iDataset] = (TH1D*)H1D_trackPt_mcparticles_rebinned[iDataset]->Clone("trackPt_efficiency_split"+Datasets[iDataset]);
      H1D_trackPt_efficiency_split[iDataset]->Reset("M");
      divideSuccess_split_notCorrected = H1D_trackPt_efficiency_split[iDataset]->Divide(H1D_trackPt_assoctracks_split_rebinned[iDataset], H1D_trackPt_mcparticles_rebinned[iDataset]);
      H1D_trackPt_efficiency_splitNotCorrected_ptHigh[iDataset] = (TH1D*)H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset]->Clone("trackPt_efficiency_splitNotCorrected_ptHigh"+Datasets[iDataset]);
      H1D_trackPt_efficiency_splitNotCorrected_ptHigh[iDataset]->Reset("M");
      divideSuccess_split__notCorrectedPtHigh = H1D_trackPt_efficiency_splitNotCorrected_ptHigh[iDataset]->Divide(H1D_trackPt_assoctracks_split_ptHigh_rebinned[iDataset], H1D_trackPt_mcparticles_ptHigh_rebinned[iDataset], 1., 1., "b");
    }

    // Merging high and low pt histograms:
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins(H1D_trackPt_efficiency[iDataset]);
    std::vector<double> xbinsVectorRight = GetTH1Bins(H1D_trackPt_efficiency_ptHigh[iDataset]);
    xbinsVectorRight.erase(xbinsVectorRight.begin());
    // cout << "xbinsVectorLeft.front() = " << xbinsVectorLeft.front() << ", xbinsVectorLeft.back() = " << xbinsVectorLeft.back() << ", xbinsVectorRight.front() = " << xbinsVectorRight.front() << ", xbinsVectorRight.back() = " << xbinsVectorRight.back() << endl;
    std::vector<double> xbinsVectorCombination = xbinsVectorLeft;
    xbinsVectorCombination.insert( xbinsVectorCombination.end(), xbinsVectorRight.begin(), xbinsVectorRight.end() );
    double* xbins_new = &xbinsVectorCombination[0];
    // cout << "xbinsVectorCombination.size() = " << xbinsVectorCombination.size() << endl;

    // making the new hist with concatenated bins
    TH1D H1D_trackPt_efficiency_concatenated_temp("H1D_trackPt_efficiency_concatenated_temp"+Datasets[iDataset], "H1D_trackPt_efficiency_concatenated_temp"+Datasets[iDataset], xbinsVectorCombination.size()-1, xbins_new);
    H1D_trackPt_efficiency_concatenated[iDataset] = (TH1D*)H1D_trackPt_efficiency_concatenated_temp.Clone("H1D_trackPt_efficiency_concatenated"+Datasets[iDataset]);
    if (useSplit == true) {
      H1D_trackPt_efficiency_splitCorrected_concatenated[iDataset] = (TH1D*)H1D_trackPt_efficiency_concatenated_temp.Clone("H1D_trackPt_efficiency_splitCorrected_concatenated"+Datasets[iDataset]);
      H1D_trackPt_efficiency_split_concatenated[iDataset] = (TH1D*)H1D_trackPt_efficiency_concatenated_temp.Clone("H1D_trackPt_efficiency_split_concatenated"+Datasets[iDataset]);
    }
    //filling the two new histograms
    for(int iBinX = 1; iBinX <= H1D_trackPt_efficiency[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_efficiency_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_efficiency[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_efficiency_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_efficiency[iDataset]->GetBinError(iBinX));
      if (useSplit == true) {
        H1D_trackPt_efficiency_splitCorrected_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_efficiency_splitCorrected[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_efficiency_splitCorrected_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_efficiency_splitCorrected[iDataset]->GetBinError(iBinX));
        H1D_trackPt_efficiency_split_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_efficiency_split[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_efficiency_split_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_efficiency_split[iDataset]->GetBinError(iBinX));
      }
    }
    for(int iBinX = 1; iBinX <= H1D_trackPt_efficiency_ptHigh[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_efficiency_concatenated[iDataset]->SetBinContent(H1D_trackPt_efficiency[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_efficiency_ptHigh[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_efficiency_concatenated[iDataset]->SetBinError(H1D_trackPt_efficiency[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_efficiency_ptHigh[iDataset]->GetBinError(iBinX));
      if (useSplit == true) {
        H1D_trackPt_efficiency_splitCorrected_concatenated[iDataset]->SetBinContent(H1D_trackPt_efficiency_splitCorrected[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_efficiency_splitCorrected_concatenated[iDataset]->SetBinError(H1D_trackPt_efficiency_splitCorrected[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_efficiency_splitCorrected_ptHigh[iDataset]->GetBinError(iBinX));
        H1D_trackPt_efficiency_split_concatenated[iDataset]->SetBinContent(H1D_trackPt_efficiency_split[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_efficiency_splitNotCorrected_ptHigh[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_efficiency_split_concatenated[iDataset]->SetBinError(H1D_trackPt_efficiency_split[iDataset]->GetNbinsX()+iBinX, H1D_trackPt_efficiency_splitNotCorrected_ptHigh[iDataset]->GetBinError(iBinX));
      }
    }
    
  }
  TH1D* H1D_trackPt_efficiency_ratios[nDatasets];
  TString DatasetsNamesPairRatio[nDatasets];
  int nHistPairRatio = (int)nDatasets / 2;
  bool divideSuccessRatio;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
      
    if (histDrawColorsOption.find("colorPairs") != std::string::npos) { //colorPairs assumes a Datasets array like so: {pair1_element1, pair1_element2, pair2_element1, pair2_element2..., pairN_element1, pairN_element2}
    cout<<"test1,1"<<endl;
      if (iDataset < nHistPairRatio) {
        DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
        H1D_trackPt_efficiency_ratios[iDataset] = (TH1D*)H1D_trackPt_efficiency_concatenated[2*iDataset]->Clone("H1D_trackPt_efficiency_concatenated_ratio"+Datasets[iDataset]);
        H1D_trackPt_efficiency_ratios[iDataset]->Reset("M");
        cout<<"test1"<<endl;
        divideSuccessRatio = H1D_trackPt_efficiency_ratios[iDataset]->Divide(H1D_trackPt_efficiency_concatenated[2*iDataset], H1D_trackPt_efficiency_concatenated[2*iDataset+1]);
        cout<<"test2"<<endl;
      }
    } else {
      cout<<"test3"<<endl;
      H1D_trackPt_efficiency_ratios[iDataset] = (TH1D*)H1D_trackPt_efficiency_concatenated[iDataset]->Clone("H1D_trackPt_efficiency_concatenated_ratio"+Datasets[iDataset]);
      H1D_trackPt_efficiency_ratios[iDataset]->Reset("M");
      divideSuccessRatio = H1D_trackPt_efficiency_ratios[iDataset]->Divide(H1D_trackPt_efficiency_concatenated[iDataset], H1D_trackPt_efficiency_concatenated[0]);
    }
  }

  TString* pdfNameEntriesNorm = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  TString* pdfNameEntriesNorm_split = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_split");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));

  std::array<std::array<float, 2>, 2> drawnWindowLog = {{{(float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX())
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinWidth(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX()))}
                                                          , {0.8, 1.2}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnWindowLogZoom = {{{(float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX())
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinWidth(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX()))}
                                                          , {0.98, 1.02}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> legendPlacementratio = {{{0.55, 0.55}, {0.9, 0.9}}}; // {{xmin, xmax}, {ymin, ymax}}                                                   
                                                    
  cout << "drawnWindowLog x: {" << drawnWindowLog[0][0] << ", " << drawnWindowLog[0][1] << "}"<< endl;
  cout << "lowEdge(0) = " << H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(0) << endl;
  cout << "lowEdge(1) = " << H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(1) << endl;
  // std::array<std::array<float, 2>, 2> drawnWindowLog = {{{-999, -999}
  //                                                         , {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  std::array<float, 2> contextPlacementEfficiency = {{0.33, 0.85}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> legendPlacementEfficiency = {{{0.55, 0.2}, {0.85, 0.6}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnWindowLogEff = {{{(float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX())
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinWidth(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX()))}
                                                          , {0, 1}}};
  std::array<std::array<float, 2>, 2> drawnWindowLogEff_SplitCorrected = {{{-999,-999},{0, 1}}};
  std::array<std::array<float, 2>, 2> drawnWindowLogEffRatio = {{{(float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(1)
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(2))/2, 
                                                          (float)(H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinLowEdge(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX())
                                                          +H1D_trackPt_efficiency_concatenated[0]->GetXaxis()->GetBinWidth(H1D_trackPt_efficiency_concatenated[0]->GetNbinsX()))}
                                                          , {0.7, 2.2}}};
  std::array<std::array<float, 2>, 2> legendPlacementEfficiencyRatio = {{{0.35, 0.5}, {0.7, 0.7}}}; 
  std::array<float, 2> contextPlacementEfficiencyRatio = {{0.35, 0.85}};     

  if (divideSuccess == true && divideSuccess_ptHigh == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_concatenated, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowLogEff, legendPlacementEfficiency, contextPlacementEfficiency, "logx,efficiency,150MevLine"+histDrawColorsOption);
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Pt_DatasetComparison" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_split == true && divideSuccess_split_ptHigh == true ) {
      Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_splitCorrected_concatenated, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowLogEff_SplitCorrected, legendPlacementEfficiency, contextPlacementEfficiency, "logx,efficiency,150MevLine"+histDrawColorsOption);
    }
    if (divideSuccess_split_notCorrected == true && divideSuccess_split__notCorrectedPtHigh ==true) {
      Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_split_concatenated, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_split, texPtMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowLogEff_SplitCorrected, legendPlacementEfficiency, contextPlacementEfficiency, "logx,efficiency,150MevLine"+histDrawColorsOption);
    }
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Pt_DatasetComparison" << endl;
  }
  TString* pdfName_ratio = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_ratio");
  if (divideSuccessRatio == true) {
    if (histDrawColorsOption.find("colorPairs") != std::string::npos) {
      Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texPtX, texRatio, texCollisionDataInfo, drawnWindowLogEffRatio, legendPlacementEfficiencyRatio, contextPlacementEfficiencyRatio, "logx,150MevLine,standardratio,ratioLine");
      TString* pdfName_ratio_zoom = new TString("track_Pt_efficiency"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_ratio_zoom");
      //Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio_zoom, texPtX, texRatio, texCollisionDataInfo, drawnWindowLogZoom, legendPlacementratio, contextPlacementAuto, "logx,150MevLine,ratioLine");
    } else {
    Draw_TH1_Histograms_in_one(H1D_trackPt_efficiency_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowLog, legendPlacementratio, contextPlacementAuto, "logx,150MevLine,standardratio,noMarkerFirst,ratioLine"+histDrawColorsOption);
    }
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}

void Draw_Efficiency_Eta_DatasetComparison(float* ptRange, bool useSplit) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackEta_mcparticles[nDatasets];

  TH1D* H1D_trackEta_assoctracks[nDatasets];
  TH1D* H1D_trackEta_assoctracks_split[nDatasets];
  TH1D* H1D_trackEta_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackEta_efficiency[nDatasets];
  // TH1D* H1D_trackEta_efficiency_split[nDatasets];
  TH1D* H1D_trackEta_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackEta_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Eta_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    if (useSplit == true) {
      H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    }
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Eta_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
    // int ibincutlow = H1D_trackEta[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackEta[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackEta[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }


    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[0]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[1]);
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);

    H1D_trackEta_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionY("trackEta_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionY("trackEta_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackEta_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionY("trackEta_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    }
    // H1D_trackEta_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionY("trackEta_assoctracks_nonSplitSecondary"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_efficiency[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency"+Datasets[iDataset]);
    H1D_trackEta_efficiency[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_efficiency[iDataset]->Divide(H1D_trackEta_assoctracks[iDataset], H1D_trackEta_mcparticles[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    // H1D_trackEta_efficiency_split[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_split"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_split[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackEta_efficiency_split[iDataset]->Divide(H1D_trackEta_assoctracks_split[iDataset], H1D_trackEta_mcparticles[iDataset]);

    // H1D_trackEta_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackEta_efficiency[iDataset]->Clone("trackEta_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_splitCorrected[iDataset]->Add(H1D_trackEta_efficiency_split[iDataset],-1.);

    if (useSplit == true) {
      H1D_trackEta_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_split"+Datasets[iDataset]);
      H1D_trackEta_efficiency_splitCorrected[iDataset]->Reset("M");
      divideSuccess_split = H1D_trackEta_efficiency_splitCorrected[iDataset]->Divide(H1D_trackEta_assoctracks_split[iDataset], H1D_trackEta_mcparticles[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
      H1D_trackEta_efficiency_splitCorrected[iDataset]->Scale(-1.);
      H1D_trackEta_efficiency_splitCorrected[iDataset]->Add(H1D_trackEta_efficiency[iDataset],1.);
    }
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackEta_mcparticles[iDataset]->Clone("trackEta_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackEta_assoctracks_nonSplitSecondary[iDataset], H1D_trackEta_mcparticles[iDataset]);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackEta_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackEta_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitTracks");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitCorrected");
  std::array<std::array<float, 2>, 2> legendPlacementEta = {{{0.55, 0.55}, {0.85, 0.85}}};
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Eta_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextPtRange(ptRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texEtaMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementEta, contextPlacementAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_split == true) {
      Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texEtaMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementEta, contextPlacementAuto, "efficiency");
    }
  }
  // Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texEtaMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
  
}
void Draw_Efficiency_Phi_DatasetComparison(float* ptRange, float* etaRange, bool useSplit) {

  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[nDatasets];

  TH1D* H1D_trackPhi_mcparticles[nDatasets];

  TH1D* H1D_trackPhi_assoctracks[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_split[nDatasets];
  TH1D* H1D_trackPhi_assoctracks_nonSplitSecondary[nDatasets];
  
  TH1D* H1D_trackPhi_efficiency[nDatasets];
  // TH1D* H1D_trackPhi_efficiency_split[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitCorrected[nDatasets];
  TH1D* H1D_trackPhi_efficiency_splitAndSecondaryCorrected[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Phi_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    if (useSplit == true) {
      H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
      }
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrackSelCollNonSplitNonPrimary"))->Clone("Draw_Efficiency_Phi_DatasetComparison_associatedtrackSelCollNonSplitSecondary"+Datasets[iDataset]);
    // int ibincutlow = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackPhi[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }


    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[0]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->FindBin(ptRange[1]);
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon + deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon - deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon);
    // H3D_trackPtEtaPhi_mcparticles[iDataset]->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    // H3D_trackPtEtaPhi_assoctracks[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_split[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);
    // H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionZ("trackPhi_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionZ("trackPhi_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    if (useSplit == true) {
      H1D_trackPhi_assoctracks_split[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionZ("trackPhi_assoctracks_split"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
      }
    // H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_nonSplitSecondary[iDataset]->ProjectionZ("trackPhi_assoctracks_nonSplitSecondary"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency"+Datasets[iDataset]);
    H1D_trackPhi_efficiency[iDataset]->Reset("M");
    cout << "  - " << DatasetsNames[iDataset] << ":" << endl;
    cout << "     phi eff #### tracks count = " << H1D_trackPhi_assoctracks[iDataset]->GetEntries() << ", part count = " << H1D_trackPhi_mcparticles[iDataset]->GetEntries() << endl;
    divideSuccess = H1D_trackPhi_efficiency[iDataset]->Divide(H1D_trackPhi_assoctracks[iDataset], H1D_trackPhi_mcparticles[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    // H1D_trackPhi_efficiency_split[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_split[iDataset]->Reset("M");
    // divideSuccess = H1D_trackPhi_efficiency_split[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset]);

    // H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_efficiency[iDataset]->Clone("trackPhi_efficiency_splitCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency_split[iDataset],-1.);
    if (useSplit == true) {
      H1D_trackPhi_efficiency_splitCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_split"+Datasets[iDataset]);
      H1D_trackPhi_efficiency_splitCorrected[iDataset]->Reset("M");
      divideSuccess_split = H1D_trackPhi_efficiency_splitCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_split[iDataset], H1D_trackPhi_mcparticles[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
      H1D_trackPhi_efficiency_splitCorrected[iDataset]->Scale(-1.);
      H1D_trackPhi_efficiency_splitCorrected[iDataset]->Add(H1D_trackPhi_efficiency[iDataset],1.);
    }
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset] = (TH1D*)H1D_trackPhi_mcparticles[iDataset]->Clone("trackPhi_efficiency_splitAndSecondaryCorrected"+Datasets[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Reset("M");
    // divideSuccess_split = H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Divide(H1D_trackPhi_assoctracks_nonSplitSecondary[iDataset], H1D_trackPhi_mcparticles[iDataset]);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Scale(-1.);
    // H1D_trackPhi_efficiency_splitAndSecondaryCorrected[iDataset]->Add(H1D_trackPhi_efficiency_splitCorrected[iDataset],1.);
  }

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_split");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  // TString* pdfNameEntriesNorm_splitAndSecondaryCorrected = new TString("track_Phi_efficiency"+dummyName[0]+"_@pt["+Form("%.2f", ptRange[0])+","+Form("%.1f", ptRange[1])+"]"+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitAndSecondaryCorrected");

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, "#splitline{"+contextPtRange(ptRange)+"}{"+contextEtaRange(etaRange)+"}", ""));
  std::array<std::array<float, 2>, 2> legendPlacementPhi = {{{0.55, 0.55}, {0.85, 0.85}}};


  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPhiMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementPhi, contextPlacementAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_split == true) {
      Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPhiMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementPhi, contextPlacementAuto, "efficiency");
    }
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
  // Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency_splitAndSecondaryCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitAndSecondaryCorrected, texPhiMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
}


void Draw_Efficiency_Eta_PtRangeComparison(float* ptRange, int nPtRanges, int iDataset, std::string options) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Eta_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Eta_PtCutComparison"+Datasets[iDataset]);

  TH1D* H1D_trackEta_mcparticles[nPtRanges];
  TH1D* H1D_trackEta_assoctracks[nPtRanges];
  
  TH1D* H1D_trackEta_efficiency[nPtRanges];

  bool divideSuccess = false;


  TString PtCutsLegend[nPtRanges];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
    // int ibincutlow = H1D_trackEta[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackEta[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackEta[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackEta[iDataset]->SetBinContent(ibin, 0);
    // }


    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt+1]);
    H3D_trackPtEtaPhi_mcparticles->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_trackPtEtaPhi_assoctracks->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);

    H1D_trackEta_mcparticles[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_mcparticles->ProjectionY("trackEta_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    H1D_trackEta_assoctracks[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_assoctracks->ProjectionY("trackEta_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, 0, -1, "e");
    // H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(1.,"trackEta_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackEta_rebinned[iDataset]);
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(1, H1D_trackEta_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackEta_rebinned[iDataset]->Scale(1./H1D_trackEta_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackEta_efficiency[iBinPt] = (TH1D*)H1D_trackEta_mcparticles[iBinPt]->Clone("trackEta_efficiency"+Datasets[iBinPt]);
    H1D_trackEta_efficiency[iBinPt]->Reset("M");
    divideSuccess = H1D_trackEta_efficiency[iBinPt]->Divide(H1D_trackEta_assoctracks[iBinPt], H1D_trackEta_mcparticles[iBinPt], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
    
    if (options.find("scaled") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackEta_efficiency[iBinPt]);
    }

    ss << "pt[" << ptRange[iBinPt] << ", " << ptRange[iBinPt+1] << "]";
    PtCutsLegend[iBinPt] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
  // for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
  //   H1D_trackEta_efficiency[iBinPt]->Scale(H1D_trackEta_efficiency[nPtRanges-1]->GetEntries()/H1D_trackEta_efficiency[iBinPt]->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
  // }

  TString optDelimiter;
  options.find("scaled") == std::string::npos ? optDelimiter = "_" : optDelimiter = "";

  TString* pdfNameEntriesNorm = new TString("track_Eta_efficiency_ptRangeComp"+Datasets[iDataset]+optDelimiter+(TString)options);

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, DatasetsNames[iDataset],  ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texEtaMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Eta_DatasetComparison" << endl;
  }
}


void Draw_Efficiency_Phi_PtRangeComparison(float* ptRange, int nPtRanges, float* etaRange, int iDataset, std::string options) {

  TH3D* H3D_trackPtEtaPhi_mcparticles;
  TH3D* H3D_trackPtEtaPhi_assoctracks;

  H3D_trackPtEtaPhi_mcparticles = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_mcparticles"+Datasets[iDataset]);
  H3D_trackPtEtaPhi_assoctracks = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Phi_PtRangeComparison_associatedtrackSelColl"+Datasets[iDataset]);
  // H3D_jetRjetPtjetEta = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_jet_r_jet_pt_jet_eta"+jetFinderQaHistType[iJetFinderQaType]))->Clone("Draw_Phi_PtCutComparison"+Datasets[iDataset]);

  TH1D* H1D_trackPhi_mcparticles[nPtRanges];
  TH1D* H1D_trackPhi_assoctracks[nPtRanges];
  
  TH1D* H1D_trackPhi_efficiency[nPtRanges];

  bool divideSuccess = false;


  TString PtCutsLegend[nPtRanges];
  std::stringstream ss;

  TString* yAxisLabel;

  for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
    // int ibincutlow = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(-0.8);
    // int ibincuthigh = H1D_trackPhi[iDataset]->GetXaxis()->FindBin(0.8);
    // for (int ibin = 1; ibin < ibincutlow; ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }
    // for (int ibin = ibincuthigh; ibin < H1D_trackPhi[iDataset]->GetNbinsX(); ibin++){
    //   H1D_trackPhi[iDataset]->SetBinContent(ibin, 0);
    // }

    int ibinPt_low = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt]);
    int ibinPt_high = H3D_trackPtEtaPhi_mcparticles->GetXaxis()->FindBin(ptRange[iBinPt+1]);
    H3D_trackPtEtaPhi_mcparticles->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    H3D_trackPtEtaPhi_assoctracks->GetXaxis()->SetRange(ibinPt_low,ibinPt_high);
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon + deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon - deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon);
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon);
    H3D_trackPtEtaPhi_mcparticles->GetXaxis()->SetRange(ibinEta_low_mcpart,ibinEta_high_mcpart);
    H3D_trackPtEtaPhi_assoctracks->GetXaxis()->SetRange(ibinEta_low_track,ibinEta_high_track);

    H1D_trackPhi_mcparticles[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_mcparticles->ProjectionZ("trackPhi_mcparticles"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_mcpart, ibinEta_high_mcpart, "e");
    H1D_trackPhi_assoctracks[iBinPt] = (TH1D*)H3D_trackPtEtaPhi_assoctracks->ProjectionZ("trackPhi_assoctracks"+Datasets[iDataset]+"_@pt["+Form("%.2f", ptRange[iBinPt])+","+Form("%.1f", ptRange[iBinPt+1])+"]", ibinPt_low, ibinPt_high, ibinEta_low_track, ibinEta_high_track, "e");
    // H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(1.,"trackPhi_rebinned"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(1, H1D_trackPhi_rebinned[iDataset]->GetNbinsX()),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
    // H1D_trackPhi_rebinned[iDataset]->Scale(1./H1D_trackPhi_rebinned[iDataset]->Integral(ibincutlow, ibincuthigh),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.

    H1D_trackPhi_efficiency[iBinPt] = (TH1D*)H1D_trackPhi_mcparticles[iBinPt]->Clone("trackPhi_efficiency"+Datasets[iBinPt]);
    H1D_trackPhi_efficiency[iBinPt]->Reset("M");
    divideSuccess = H1D_trackPhi_efficiency[iBinPt]->Divide(H1D_trackPhi_assoctracks[iBinPt], H1D_trackPhi_mcparticles[iBinPt], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
    
    if (options.find("scaled") != std::string::npos) {
        NormaliseYieldToIntegral(H1D_trackPhi_efficiency[iBinPt]);
    }

    ss << "pt[" << ptRange[iBinPt] << ", " << ptRange[iBinPt+1] << "]";
    PtCutsLegend[iBinPt] = (TString)ss.str();
    ss.str("");
    ss.clear();
  }
//   for(int iBinPt = 0; iBinPt < nPtRanges; iBinPt++){
//     H1D_trackPhi_efficiency[iBinPt]->Scale(H1D_trackPhi_efficiency[nPtRanges-1]->GetEntries()/H1D_trackPhi_efficiency[iBinPt]->GetEntries(),"width"); // If option contains "width" the bin contents and errors are divided by the bin width.
//   }

  TString optDelimiter;
  options.find("scaled") == std::string::npos ? optDelimiter = "_" : optDelimiter = "";

  TString* pdfNameEntriesNorm = new TString("track_Phi_efficiency_ptRangeComp"+Datasets[iDataset]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]"+optDelimiter+(TString)options);

  // TString textContext(contextDatasetComp(""));
  TString textContext(contextCustomThreeFields(*texDatasetsComparisonCommonDenominator, DatasetsNames[iDataset], contextEtaRange(etaRange), ""));

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_efficiency, PtCutsLegend, nPtRanges, textContext, pdfNameEntriesNorm, texPhiMC, texTrackEfficiency, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Phi_DatasetComparison" << endl;
  }
}



void Draw_mcPt_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vx[nDatasets];
  TH1D* H1D_vx_rebinned[nDatasets];

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vx[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"ptmc"))->Clone("Draw_mcPt_DatasetComparison"+Datasets[iDataset]);
    H1D_vx_rebinned[iDataset] = (TH1D*)H1D_vx[iDataset]->Rebin(1.,"ptmc_rebinned"+Datasets[iDataset]);

    NormaliseYieldToNEntries(H1D_vx_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Pt_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vx_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtMC, texTrackPtYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
}

void Draw_Vx_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vx[nDatasets];
  TH1D* H1D_vx_rebinned[nDatasets];

  TH1D* dummy = new TH1D("h1_2", "h1 title_2", 5000, -1, 1);

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vx[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"vxmc"))->Clone("Draw_Vx_DatasetComparison"+Datasets[iDataset]);
    if (particleStatusOption == 1) {
      std::vector<double> O2H1DYbinsVector = GetTH1Bins(dummy);
      double* O2yBins = &O2H1DYbinsVector[0];
      H1D_vx_rebinned[iDataset] = (TH1D*)H1D_vx[iDataset]->Rebin(5000,"vxmc_rebinned"+Datasets[iDataset], O2yBins);
    }
    else {
      H1D_vx_rebinned[iDataset] = (TH1D*)H1D_vx[iDataset]->Rebin(50.,"vxmc_rebinned"+Datasets[iDataset]);
    }

    NormaliseYieldToNEntries(H1D_vx_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Vx_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vx_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texVx, texPartVxYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
}

void Draw_Vy_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vy[nDatasets];
  TH1D* H1D_vy_rebinned[nDatasets];

  TH1D* dummy = new TH1D("h1_2", "h1 title_2", 5000, -1, 1);

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vy[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"vymc"))->Clone("Draw_Vy_DatasetComparison"+Datasets[iDataset]);
    if (particleStatusOption == 1) {
      std::vector<double> O2H1DYbinsVector = GetTH1Bins(dummy);
      double* O2yBins = &O2H1DYbinsVector[0];
      H1D_vy_rebinned[iDataset] = (TH1D*)H1D_vy[iDataset]->Rebin(5000,"vymc_rebinned"+Datasets[iDataset], O2yBins);
    }
    else {
      H1D_vy_rebinned[iDataset] = (TH1D*)H1D_vy[iDataset]->Rebin(50.,"vymc_rebinned"+Datasets[iDataset]);
    }

    NormaliseYieldToNEntries(H1D_vy_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Vy_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vy_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texVy, texPartVyYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
}

void Draw_Vz_DatasetComparison(int particleStatusOption) {
  TH1D* H1D_vz[nDatasets];
  TH1D* H1D_vz_rebinned[nDatasets];

  TString Folder;
  TString Name;
  if (particleStatusOption == 0) {
    Folder = "";
    Name = "";
  }
  if (particleStatusOption == 1) {
    Folder = "Primaries/";
    Name = "Primaries_";
  }
  if (particleStatusOption == 2) {
    Folder = "TransportSecondaries/";
    Name = "TransportSecondaries_";
  }

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H1D_vz[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/Kine/"+Folder+"vzmc"))->Clone("Draw_Vz_DatasetComparison"+Datasets[iDataset]);
    H1D_vz_rebinned[iDataset] = (TH1D*)H1D_vz[iDataset]->Rebin(50.,"vzmc_rebinned"+Datasets[iDataset]);

    NormaliseYieldToNEntries(H1D_vz_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_"+Name+"Vz_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_vz_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texVz, texPartVzYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logy");
}


void Draw_y_iu_DatasetComparison() {
  TH1D* H1D_IUy[nDatasets];
  TH1D* H1D_IUy_rebinned[nDatasets];

  TH1D* dummy = new TH1D("h1", "h1 title", 500, -20, 20);

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    H1D_IUy[iDataset] = (TH1D*)((TH1D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/IU/y"))->Clone("Draw_Y_IU_DatasetComparison"+Datasets[iDataset]);
    // H1D_IUy[iDataset]->GetXaxis()->SetRangeUser(-20,20);
    std::vector<double> O2H1DYbinsVector = GetTH1Bins(dummy);
    double* O2yBins = &O2H1DYbinsVector[0];
    H1D_IUy_rebinned[iDataset] = (TH1D*)H1D_IUy[iDataset]->Rebin(500, "track_y_iu_rebinned"+Datasets[iDataset], O2yBins);
    // H1D_IUy_rebinned[iDataset]->GetXaxis()->SetRangeUser(-20,20);

    // int ibinYlow = H1D_IUy_rebinned[iDataset]->GetXaxis()->FindBin(-20);
    // int ibinYhigh = H1D_IUy_rebinned[iDataset]->GetXaxis()->FindBin(20);

    NormaliseYieldToNEntries(H1D_IUy_rebinned[iDataset]);
  }

  TString* pdfName = new TString("track_IU_y_DataComp");

  TString textContext("");

  Draw_TH1_Histograms_in_one(H1D_IUy_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texIUy, texPartYYield_EntriesNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
}



void Draw_Purity_Pt_DatasetComparison(float* etaRange, bool useSplit) {
  //low pT
  TH3D* H3D_trackPtEtaPhi_primaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_secondaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_splitPrimaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_splitSecondaryTracks[nDatasets];
  //high pT 
  TH3D* H3D_trackPtEtaPhi_high_primaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_high_secondaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_high_splitPrimaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_high_splitSecondaryTracks[nDatasets];

  //low pT
  TH1D* H1D_trackPt_primary[nDatasets];
  TH1D* H1D_trackPt_secondary[nDatasets];
  TH1D* H1D_trackPt_splitPrimaryTracks[nDatasets];
  TH1D* H1D_trackPt_splitSecondaryTracks[nDatasets];

  TH1D* H1D_trackPt_primary_rebinned[nDatasets] ;
  TH1D* H1D_trackPt_secondary_rebinned[nDatasets] ;
  TH1D* H1D_trackPt_splitPrimaryTracks_rebinned[nDatasets] ;
  TH1D* H1D_trackPt_splitSecondaryTracks_rebinned[nDatasets];
  
  //high pT
  TH1D* H1D_trackPtHigh_primary[nDatasets];
  TH1D* H1D_trackPtHigh_secondary[nDatasets];
  TH1D* H1D_trackPtHigh_splitPrimaryTracks[nDatasets];
  TH1D* H1D_trackPtHigh_splitSecondaryTracks[nDatasets];

  TH1D* H1D_trackPtHigh_primary_rebinned[nDatasets] ;
  TH1D* H1D_trackPtHigh_secondary_rebinned[nDatasets] ;
  TH1D* H1D_trackPtHigh_splitPrimaryTracks_rebinned[nDatasets];
  TH1D* H1D_trackPtHigh_splitSecondaryTracks_rebinned[nDatasets];

  //low pT histos for declaring ratio
  TH1D* H1D_trackPt_primaryPurity[nDatasets];
  TH1D* H1D_primaryPurity_denominator[nDatasets];
  TH1D* H1D_trackPt_split[nDatasets];
  TH1D* H1D_trackPt_nonsplit[nDatasets];

  TH1D* H1D_trackPt_purity_concatenated[nDatasets];
  TH1D* H1D_trackPt_purity_concatenatedSplitCorrected[nDatasets];
  TH1D* H1D_trackPt_purity_concatenatedSplit[nDatasets];
  
  TH1D* H1D_trackPtHigh_primaryPurity[nDatasets];
  TH1D* H1D_primaryPurity_High_denominator[nDatasets];

  TH1D* H1D_numerator_corrected[nDatasets];
  TH1D* H1D_numerator_corrected_high[nDatasets];

  TH1D* H1D_trackPt_primaryPurity_SplitCorrected[nDatasets];
  TH1D* H1D_trackPtHigh_primaryPurity_SplitCorrected[nDatasets];
  TH1D* H1D_trackPt_primaryPurity_Split[nDatasets];
  TH1D* H1D_trackPtHigh_primaryPurity_Split[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_High = false;
  bool divideSuccess_SplitCorrected = false ;
  bool divideSuccess_SplitCorrected_High = false ;
  bool divideSuccess_Split = false ;
  bool divideSuccess_Split_High = false ;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    ////////////////////////// low pT 

    //primary (split + nonsplit)
    H3D_trackPtEtaPhi_primaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_primary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_primary"+Datasets[iDataset]);
    //non primary (split + non split)
    H3D_trackPtEtaPhi_secondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_nonprimary"+Datasets[iDataset]);
    if (useSplit == true) {
      //split primary 
      H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_splitPrimary"+Datasets[iDataset]);
      //split non primary 
      H3D_trackPtEtaPhi_splitSecondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_split_nonprimary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_splitnonPrimary"+Datasets[iDataset]); 
    }
    int ibinEta_low_track = H3D_trackPtEtaPhi_secondaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon);
    int ibinEta_high_track = H3D_trackPtEtaPhi_secondaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon);
    
    H1D_trackPt_primary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_primaryTracks[iDataset]->ProjectionX("trackPt_primary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    H1D_trackPt_secondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_secondaryTracks[iDataset]->ProjectionX("trackPt_secondary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPt_splitPrimaryTracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset]->ProjectionX("trackPt_splitPrimaryTracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
      H1D_trackPt_splitSecondaryTracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_splitSecondaryTracks[iDataset]->ProjectionX("trackPt_splitNonPrimaryTracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    }
    
     //////////////////////////// high pT

    H3D_trackPtEtaPhi_high_primaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_high_track_eta_track_phi_associatedtrack_primary"))->Clone("Draw_Purity_PtHigh_DatasetComparison_associatedtrack_primary"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_high_secondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_high_track_eta_track_phi_associatedtrack_nonprimary"))->Clone("Draw_Purity_PtHigh_DatasetComparison_associatedtrack_nonprimary"+Datasets[iDataset]);

    if (useSplit == true) {
      //split primary 
      H3D_trackPtEtaPhi_high_splitPrimaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_high_track_eta_track_phi_associatedtrack_split_primary"))->Clone("Draw_Purity_PtHigh_DatasetComparison_associatedtrack_splitPrimary"+Datasets[iDataset]);
      //non split primary 
      H3D_trackPtEtaPhi_high_splitSecondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_high_track_eta_track_phi_associatedtrack_split_nonprimary"))->Clone("Draw_Purity_PtHigh_DatasetComparison_associatedtrack_splitnonPrimary"+Datasets[iDataset]); //new
    }

    H1D_trackPtHigh_primary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_high_primaryTracks[iDataset]->ProjectionX("trackPtHigh_primary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    H1D_trackPtHigh_secondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_high_secondaryTracks[iDataset]->ProjectionX("trackPtHigh_secondary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPtHigh_splitPrimaryTracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_high_splitPrimaryTracks[iDataset]->ProjectionX("trackPtHigh_splitPrimaryTracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
      H1D_trackPtHigh_splitSecondaryTracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_high_splitSecondaryTracks[iDataset]->ProjectionX("trackPtHihg_splitNonPrimaryTracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    }

    //tweaking the low-pt bins to make them larger close to 10GeV
    std::vector<double> xbinsVectorInitialLow = GetTH1Bins(H1D_trackPt_primary[iDataset]);
    double* xbinsInitialLow = &xbinsVectorInitialLow[0];

    double ptBinsLowNew[500]; //500 to have a good margin
    double lastPt;
    int iBinPtNew = 0;
    int iBinPtInitialHisto = 0;
    int increment;
    while(iBinPtInitialHisto < H1D_trackPt_primary[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialLow[iBinPtInitialHisto];
      if (lastPt < 0.300) {
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPt_primary[iDataset]->GetXaxis()->GetBinWidth(1) /8;
      }
      ptBinsLowNew[iBinPtNew] = lastPt;
      iBinPtInitialHisto += increment;
      iBinPtNew += 1;
    }
    ptBinsLowNew[iBinPtNew-1] = 10; // replace last set bin edge by 10
    int nBinsLowNew = iBinPtNew - 1;
    if (nBinsLowNew >= 500){
      cout << "NEEDS TO RESERVE MORE MEMORY FOR ptBinsLowNew !!!!!!!!!!!!!!!!!!" << endl;
    }
    H1D_trackPt_primary_rebinned[iDataset] = (TH1D*)H1D_trackPt_primary[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_primary_rebinned"+Datasets[iDataset], ptBinsLowNew);
    H1D_trackPt_secondary_rebinned[iDataset] = (TH1D*)H1D_trackPt_secondary[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_secondary_rebinned"+Datasets[iDataset], ptBinsLowNew);
    if (useSplit == true) {
      H1D_trackPt_splitPrimaryTracks_rebinned[iDataset] = (TH1D*)H1D_trackPt_splitPrimaryTracks[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_splitPrimaryTracks_rebinned"+Datasets[iDataset], ptBinsLowNew);
      H1D_trackPt_splitSecondaryTracks_rebinned[iDataset] = (TH1D*)H1D_trackPt_splitSecondaryTracks[iDataset]->Rebin(nBinsLowNew, "H1D_trackPt_splitNonPrimaryTracks_rebinned"+Datasets[iDataset], ptBinsLowNew);
    }
    //tweaking the high-pt bins to have them less prone to statistical fluctuations
    std::vector<double> xbinsVectorInitialHigh = GetTH1Bins(H1D_trackPtHigh_primary[iDataset]);
    double* xbinsInitialHigh = &xbinsVectorInitialHigh[0];

    double ptBinsHighNew[500]; //500 to have a good margin
    iBinPtNew = 0;
    iBinPtInitialHisto = 0;
    while(iBinPtInitialHisto < H1D_trackPtHigh_primary[iDataset]->GetNbinsX()){
      lastPt = xbinsInitialHigh[iBinPtInitialHisto];
      if (lastPt < 25) {
        increment = 1;
      } else {
        increment = 1+ lastPt/H1D_trackPtHigh_primary[iDataset]->GetXaxis()->GetBinWidth(1) ;
      }
      ptBinsHighNew[iBinPtNew] = lastPt;
      iBinPtInitialHisto += increment;
      iBinPtNew += 1;
    }
    ptBinsHighNew[iBinPtNew-1] = 100; // replace last set bin edge by 10
    int nBinsHighNew = iBinPtNew - 1;
    if (nBinsHighNew >= 500){
      cout << "NEEDS TO RESERVE MORE MEMORY FOR ptBinsHighNew !!!!!!!!!!!!!!!!!!" << endl;
    }
    H1D_trackPtHigh_primary_rebinned[iDataset] = (TH1D*)H1D_trackPtHigh_primary[iDataset]->Rebin(nBinsHighNew, "H1D_trackPtHigh_primary_rebinned"+Datasets[iDataset], ptBinsHighNew);
    H1D_trackPtHigh_secondary_rebinned[iDataset] = (TH1D*)H1D_trackPtHigh_secondary[iDataset]->Rebin(nBinsHighNew, "H1D_trackPtHigh_secondary_rebinned"+Datasets[iDataset], ptBinsHighNew);
    if (useSplit == true) {
      H1D_trackPtHigh_splitPrimaryTracks_rebinned[iDataset] = (TH1D*)H1D_trackPtHigh_splitPrimaryTracks[iDataset]->Rebin(nBinsHighNew, "H1D_trackPtHigh_assoctracks_split_rebinned"+Datasets[iDataset], ptBinsHighNew);
      H1D_trackPtHigh_splitSecondaryTracks_rebinned[iDataset] = (TH1D*)H1D_trackPtHigh_splitSecondaryTracks[iDataset]->Rebin(nBinsHighNew, "H1D_trackPtHigh_assoctracks_split_rebinned"+Datasets[iDataset], ptBinsHighNew);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////// declaring the ratios ///////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // % of primary
    H1D_trackPt_primaryPurity[iDataset] = (TH1D*)H1D_trackPt_primary_rebinned[iDataset]->Clone("trackPt_primaryPurity_rebinned"+Datasets[iDataset]);
    H1D_trackPt_primaryPurity[iDataset]->Reset("M");
    
    H1D_primaryPurity_denominator[iDataset] = (TH1D*)H1D_trackPt_primary_rebinned[iDataset]->Clone("trackPt_primaryPurity_denominator_rebinned1"+Datasets[iDataset]);
    H1D_primaryPurity_denominator[iDataset]->Add(H1D_trackPt_secondary_rebinned[iDataset],1.);
    divideSuccess = H1D_trackPt_primaryPurity[iDataset]->Divide(H1D_trackPt_primary_rebinned[iDataset], H1D_primaryPurity_denominator[iDataset]);

    H1D_trackPtHigh_primaryPurity[iDataset] = (TH1D*)H1D_trackPtHigh_primary_rebinned[iDataset]->Clone("trackPtHigh_primaryPurity_rebinned"+Datasets[iDataset]);
    H1D_trackPtHigh_primaryPurity[iDataset]->Reset("M");
    
    H1D_primaryPurity_High_denominator[iDataset] = (TH1D*)H1D_trackPtHigh_primary_rebinned[iDataset]->Clone("trackPtHigh_primaryPurity_denominator_rebinned1"+Datasets[iDataset]);
    H1D_primaryPurity_High_denominator[iDataset]->Add(H1D_trackPtHigh_secondary_rebinned[iDataset],1.);
    divideSuccess_High = H1D_trackPtHigh_primaryPurity[iDataset]->Divide(H1D_trackPtHigh_primary_rebinned[iDataset], H1D_primaryPurity_High_denominator[iDataset]);

    
    if (useSplit == true) {

    // % of primary with split corrected
      H1D_trackPt_primaryPurity_SplitCorrected[iDataset] = (TH1D*)H1D_trackPt_primary_rebinned[iDataset]->Clone("trackPt_primaryPurityCorrected_rebinned"+Datasets[iDataset]);
      H1D_trackPt_primaryPurity_SplitCorrected[iDataset]->Reset("M");

      H1D_numerator_corrected[iDataset] = (TH1D*)H1D_trackPt_primary_rebinned[iDataset]->Clone("trackPt_primaryPurity_numerator_corrected"+Datasets[iDataset]);
      H1D_numerator_corrected[iDataset]->Add(H1D_trackPt_splitPrimaryTracks_rebinned[iDataset], -1);
    
      divideSuccess_SplitCorrected = H1D_trackPt_primaryPurity_SplitCorrected[iDataset]->Divide(H1D_numerator_corrected[iDataset], H1D_primaryPurity_denominator[iDataset]);

      H1D_trackPtHigh_primaryPurity_SplitCorrected[iDataset] = (TH1D*)H1D_trackPtHigh_primary_rebinned[iDataset]->Clone("trackPtHigh_primaryPurityCorrected_rebinned"+Datasets[iDataset]);
      H1D_trackPtHigh_primaryPurity_SplitCorrected[iDataset]->Reset("M");

      H1D_numerator_corrected_high[iDataset] = (TH1D*)H1D_trackPtHigh_primary_rebinned[iDataset]->Clone("trackPtHigh_primaryPurity_numerator_corrected"+Datasets[iDataset]);
      H1D_numerator_corrected_high[iDataset]->Add(H1D_trackPtHigh_splitPrimaryTracks_rebinned[iDataset], -1);
    
      divideSuccess_SplitCorrected_High = H1D_trackPtHigh_primaryPurity_SplitCorrected[iDataset]->Divide(H1D_numerator_corrected_high[iDataset], H1D_primaryPurity_High_denominator[iDataset]);

      // % of split primary 
      H1D_trackPt_primaryPurity_Split[iDataset] = (TH1D*)H1D_trackPt_splitPrimaryTracks_rebinned[iDataset]->Clone("trackPt_primaryPuritySplit_rebinned"+Datasets[iDataset]);
      H1D_trackPt_primaryPurity_Split[iDataset]->Reset("M");
    
      divideSuccess_Split = H1D_trackPt_primaryPurity_Split[iDataset]->Divide(H1D_trackPt_splitPrimaryTracks_rebinned[iDataset], H1D_primaryPurity_denominator[iDataset]);

      H1D_trackPtHigh_primaryPurity_Split[iDataset] = (TH1D*)H1D_trackPtHigh_splitPrimaryTracks_rebinned[iDataset]->Clone("trackPtHigh_primaryPuritySplit_rebinned"+Datasets[iDataset]);
      H1D_trackPtHigh_primaryPurity_Split[iDataset]->Reset("M");
    
      divideSuccess_Split_High = H1D_trackPtHigh_primaryPurity_Split[iDataset]->Divide(H1D_trackPtHigh_splitPrimaryTracks_rebinned[iDataset], H1D_primaryPurity_High_denominator[iDataset]);

    }

    // Merging high and low pt histograms:
    // x-axis
    std::vector<double> xbinsVectorLeft = GetTH1Bins(H1D_trackPt_primaryPurity[iDataset]);
    std::vector<double> xbinsVectorRight = GetTH1Bins(H1D_trackPtHigh_primaryPurity[iDataset]);
    xbinsVectorRight.erase(xbinsVectorRight.begin());
    // cout << "xbinsVectorLeft.front() = " << xbinsVectorLeft.front() << ", xbinsVectorLeft.back() = " << xbinsVectorLeft.back() << ", xbinsVectorRight.front() = " << xbinsVectorRight.front() << ", xbinsVectorRight.back() = " << xbinsVectorRight.back() << endl;
    std::vector<double> xbinsVectorCombination = xbinsVectorLeft;
    xbinsVectorCombination.insert( xbinsVectorCombination.end(), xbinsVectorRight.begin(), xbinsVectorRight.end() );
    double* xbins_new = &xbinsVectorCombination[0];
    // cout << "xbinsVectorCombination.size() = " << xbinsVectorCombination.size() << endl;

    // making the new hist with concatenated bins
    TH1D H1D_trackPt_purity_concatenated_temp("H1D_trackPt_purity_concatenated_temp"+Datasets[iDataset], "H1D_trackPt_purity_concatenated_temp"+Datasets[iDataset], xbinsVectorCombination.size()-1, xbins_new);
    H1D_trackPt_purity_concatenated[iDataset] = (TH1D*)H1D_trackPt_purity_concatenated_temp.Clone("H1D_trackPt_efficiency_concatenated"+Datasets[iDataset]);
    if (useSplit == true) {
      H1D_trackPt_purity_concatenatedSplitCorrected[iDataset] = (TH1D*)H1D_trackPt_purity_concatenated_temp.Clone("H1D_trackPt_purity_concatenatedSplitCorrected"+Datasets[iDataset]);
      H1D_trackPt_purity_concatenatedSplit[iDataset] = (TH1D*)H1D_trackPt_purity_concatenated_temp.Clone("H1D_trackPt_purity_concatenatedSplit"+Datasets[iDataset]);
    }
   //filling the two new histograms
    for(int iBinX = 1; iBinX <= H1D_trackPt_primaryPurity[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_purity_concatenated[iDataset]->SetBinContent(iBinX, H1D_trackPt_primaryPurity[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_purity_concatenated[iDataset]->SetBinError(iBinX, H1D_trackPt_primaryPurity[iDataset]->GetBinError(iBinX));
      if (useSplit == true) {
        H1D_trackPt_purity_concatenatedSplitCorrected[iDataset]->SetBinContent(iBinX, H1D_trackPt_primaryPurity_SplitCorrected[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_purity_concatenatedSplitCorrected[iDataset]->SetBinError(iBinX, H1D_trackPt_primaryPurity_SplitCorrected[iDataset]->GetBinError(iBinX));
        H1D_trackPt_purity_concatenatedSplit[iDataset]->SetBinContent(iBinX, H1D_trackPt_primaryPurity_Split[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_purity_concatenatedSplit[iDataset]->SetBinError(iBinX, H1D_trackPt_primaryPurity_Split[iDataset]->GetBinError(iBinX));
      }
    }
    for(int iBinX = 1; iBinX <= H1D_trackPtHigh_primaryPurity[iDataset]->GetNbinsX(); iBinX++){
      H1D_trackPt_purity_concatenated[iDataset]->SetBinContent(H1D_trackPt_primaryPurity[iDataset]->GetNbinsX()+iBinX, H1D_trackPtHigh_primaryPurity[iDataset]->GetBinContent(iBinX));
      H1D_trackPt_purity_concatenated[iDataset]->SetBinError(H1D_trackPt_primaryPurity[iDataset]->GetNbinsX()+iBinX, H1D_trackPtHigh_primaryPurity[iDataset]->GetBinError(iBinX));
      if (useSplit == true) {
        H1D_trackPt_purity_concatenatedSplitCorrected[iDataset]->SetBinContent(H1D_trackPt_primaryPurity_SplitCorrected[iDataset]->GetNbinsX()+iBinX, H1D_trackPtHigh_primaryPurity_SplitCorrected[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_purity_concatenatedSplitCorrected[iDataset]->SetBinError(H1D_trackPt_primaryPurity_SplitCorrected[iDataset]->GetNbinsX()+iBinX, H1D_trackPtHigh_primaryPurity_SplitCorrected[iDataset]->GetBinError(iBinX));
        H1D_trackPt_purity_concatenatedSplit[iDataset]->SetBinContent(H1D_trackPt_primaryPurity_Split[iDataset]->GetNbinsX()+iBinX, H1D_trackPtHigh_primaryPurity_Split[iDataset]->GetBinContent(iBinX));
        H1D_trackPt_purity_concatenatedSplit[iDataset]->SetBinError(H1D_trackPt_primaryPurity_Split[iDataset]->GetNbinsX()+iBinX, H1D_trackPtHigh_primaryPurity_Split[iDataset]->GetBinError(iBinX));
      }
    }
  }
    TH1D* H1D_trackPt_purity_ratios[nDatasets];
    TString DatasetsNamesPairRatio[nDatasets];
    int nHistPairRatio = (int)nDatasets / 2;
    bool divideSuccessRatio;
    for(int iDataset = 0; iDataset < nDatasets; iDataset++){
      if (histDrawColorsOption.find("colorPairs") != std::string::npos) { //colorPairs assumes a Datasets array like so: {pair1_element1, pair1_element2, pair2_element1, pair2_element2..., pairN_element1, pairN_element2}
        if (iDataset < nHistPairRatio) {
          DatasetsNamesPairRatio[iDataset] = DatasetsNames[2*iDataset]+(TString)"/"+DatasetsNames[2*iDataset+1];
          H1D_trackPt_purity_ratios[iDataset] = (TH1D*)H1D_trackPt_purity_concatenated[2*iDataset]->Clone("H1D_trackPt_purity_concatenated_ratio"+Datasets[iDataset]);
          H1D_trackPt_purity_ratios[iDataset]->Reset("M");
          divideSuccessRatio = H1D_trackPt_purity_ratios[iDataset]->Divide(H1D_trackPt_purity_concatenated[2*iDataset], H1D_trackPt_purity_concatenated[2*iDataset+1]);
        }
      } else {
        H1D_trackPt_purity_ratios[iDataset] = (TH1D*)H1D_trackPt_purity_concatenated[iDataset]->Clone("H1D_trackPt_purity_concatenated_ratio"+Datasets[iDataset]);
        H1D_trackPt_purity_ratios[iDataset]->Reset("M");
        divideSuccessRatio = H1D_trackPt_purity_ratios[iDataset]->Divide(H1D_trackPt_purity_concatenated[iDataset], H1D_trackPt_purity_concatenated[0]);
      }
    }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////// drawing histograms //////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  TString* pdfNameEntriesNorm = new TString("track_Pt_primaryPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_splitCorrected = new TString("track_Pt_purity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");
  TString* pdfNameEntriesNorm_split = new TString("track_Pt_purity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_split");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));

  std::array<std::array<float, 2>, 2> drawnWindowPurity = {{{-999, -999} , {0.6, 1.2}}}; // {{xmin, xmax}, {ymin, ymax}} //0,6 et 1.2
  std::array<std::array<float, 2>, 2> drawnWindowPuritySplit = {{{-999, -999} , {0, 1.2}}}; 
  std::array<float, 2> contextPlacementPurity = {{0.20, 0.85}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<float, 2> contextPlacementPurityRatio = {{0.2, 0.85}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> legendPlacementPurity = {{{0.45, 0.2}, {0.8, 0.55}}}; // {{xmin, xmax}, {ymin, ymax}}

  if (divideSuccess == true && divideSuccess_High ==true ) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_purity_concatenated, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPtRec, texTrackPurity, texCollisionDataInfo, drawnWindowPurity, legendPlacementPurity, contextPlacementPurity, "logx,efficiency,150MevLine"+histDrawColorsOption);
  }
  else {
    cout << "Divide failed in Draw_Purity_Pt_DatasetComparison_legacy" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_SplitCorrected == true && divideSuccess_SplitCorrected_High == true) {
      Draw_TH1_Histograms_in_one(H1D_trackPt_purity_concatenatedSplitCorrected, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_splitCorrected, texPtRec, texTrackPurity, texCollisionDataInfo, drawnWindowPurity, legendPlacementPurity, contextPlacementAuto, "logx,efficiency,150MevLine"+histDrawColorsOption);
    }
    else {
    cout << "Divide failed in Draw_Purity_Pt_DatasetComparison_legacy" << endl;
    }
    if (divideSuccess_Split == true && divideSuccess_Split_High == true) {
      Draw_TH1_Histograms_in_one(H1D_trackPt_purity_concatenatedSplit, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_split, texPtRec, texTrackPurity, texCollisionDataInfo, drawnWindowPuritySplit, legendPlacementPurity, contextPlacementAuto, "logx,efficiency,150MevLine"+histDrawColorsOption);
    }
    else {
    cout << "Divide failed in Draw_Purity_Pt_DatasetComparison_legacy" << endl;
    }
  }


  std::array<std::array<float, 2>, 2> drawnWindowLog = {{{-999, -999} , {0.8, 1.2}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnWindowLogZoom = {{{-999, -999} , {0.98, 1.18}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> legendPlacementPurityRatio = {{{0.2, 0.2}, {0.8, 0.4}}}; // {{xmin, xmax}, {ymin, ymax}}
  TString* pdfName_ratio = new TString("track_Pt_purity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_ratio");
  if (divideSuccessRatio == true) {
    if (histDrawColorsOption.find("colorPairs") != std::string::npos) {
      Draw_TH1_Histograms_in_one(H1D_trackPt_purity_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio, texPtX, texRatio, texCollisionDataInfo, drawnWindowLog, legendPlacementPurityRatio, contextPlacementPurityRatio, "logx,150MevLine,standardratio,ratioLine");
      TString* pdfName_ratio_zoom = new TString("track_Pt_purity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_ratio_zoom");
      //Draw_TH1_Histograms_in_one(H1D_trackPt_purity_ratios, DatasetsNamesPairRatio, nHistPairRatio, textContext, pdfName_ratio_zoom, texPtX, texRatio, texCollisionDataInfo, drawnWindowLogZoom, legendPlacementAuto, contextPlacementAuto, "logx,150MevLine,ratioLine");
    } else {
    Draw_TH1_Histograms_in_one(H1D_trackPt_purity_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtX, texRatioDatasets, texCollisionDataInfo, drawnWindowLog, legendPlacementAuto, contextPlacementAuto, "logx,150MevLine,standardratio,noMarkerFirst,ratioLine"+histDrawColorsOption);
    }
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}



void Draw_Purity_Eta_DatasetComparison(float* etaRange, bool useSplit) {
  TH3D* H3D_trackPtEtaPhi_primaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_secondaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_splitPrimaryTracks[nDatasets];

  TH1D* H1D_trackEta_primary[nDatasets];
  TH1D* H1D_trackEta_secondary[nDatasets];
  TH1D* H1D_trackEta_splitPrimaryTracks[nDatasets];
  
  TH1D* H1D_trackEta_primaryPurity[nDatasets];
  TH1D* H1D_primaryPurity_denominator[nDatasets];
  TH1D* H1D_trackEta_splitPurity[nDatasets];
  TH1D* H1D_splitPurity_denominator[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_primaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_primary"))->Clone("Draw_Purity_Eta_DatasetComparison_associatedtrack_primary"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_secondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"))->Clone("Draw_Purity_Eta_DatasetComparison_associatedtrack_nonprimary"+Datasets[iDataset]);
    if (useSplit == true) {
      H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"))->Clone("Draw_Purity_Eta_DatasetComparison_associatedtrack_splitPrimary"+Datasets[iDataset]);
    }
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon + deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon - deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_secondaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon);
    int ibinEta_high_track = H3D_trackPtEtaPhi_secondaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon);

    H1D_trackEta_primary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_primaryTracks[iDataset]->ProjectionY("trackEta_primary"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    H1D_trackEta_secondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_secondaryTracks[iDataset]->ProjectionY("trackEta_secondary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackEta_splitPrimaryTracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset]->ProjectionY("trackEta_splitPrimaryTracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    }   
    // naming the histogram and resetting it to have a chosen name
    H1D_trackEta_primaryPurity[iDataset] = (TH1D*)H1D_trackEta_primary[iDataset]->Clone("trackEta_primaryPurity"+Datasets[iDataset]);
    H1D_trackEta_primaryPurity[iDataset]->Reset("M");

    H1D_primaryPurity_denominator[iDataset] = (TH1D*)H1D_trackEta_primary[iDataset]->Clone("trackEta_primaryPurity_denominator"+Datasets[iDataset]);
    H1D_primaryPurity_denominator[iDataset]->Add(H1D_trackEta_secondary[iDataset],1.);
    divideSuccess = H1D_trackEta_primaryPurity[iDataset]->Divide(H1D_trackEta_primary[iDataset], H1D_primaryPurity_denominator[iDataset]);

    // naming the histogram and resetting it to have a chosen name
    if (useSplit == true) {
      H1D_trackEta_splitPurity[iDataset] = (TH1D*)H1D_trackEta_primary[iDataset]->Clone("trackEta_splitPurity"+Datasets[iDataset]);
      H1D_trackEta_splitPurity[iDataset]->Reset("M");

      H1D_splitPurity_denominator[iDataset] = (TH1D*)H1D_trackEta_primary[iDataset]->Clone("trackEta_splitPurity_denominator"+Datasets[iDataset]);
      H1D_splitPurity_denominator[iDataset]->Add(H1D_trackEta_splitPrimaryTracks[iDataset],1.);
      divideSuccess_split = H1D_trackEta_splitPurity[iDataset]->Divide(H1D_trackEta_primary[iDataset], H1D_splitPurity_denominator[iDataset]);
    }
  }

  TString* pdfNameEntriesNorm = new TString("track_Eta_primaryPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_zoom = new TString("track_Eta_primaryPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_zoom");
  TString* pdfNameEntriesNorm_split = new TString("track_Eta_splitPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split_zoom = new TString("track_Eta_splitPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_zoom");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));


  std::array<std::array<float, 2>, 2> drawnZoom_prim = {{{(float)H1D_trackEta_primaryPurity[0]->GetXaxis()->GetXmin(), (float)H1D_trackEta_primaryPurity[0]->GetXaxis()->GetXmax()}, 
                                                    {0.9, 1}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnZoom_split = {{{(float)H1D_trackEta_primaryPurity[0]->GetXaxis()->GetXmin(), (float)H1D_trackEta_primaryPurity[0]->GetXaxis()->GetXmax()}, 
                                                    {0.996, 1.002}}}; // {{xmin, xmax}, {ymin, ymax}}

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_primaryPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texEtaMC, texTrackPurity, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
    Draw_TH1_Histograms_in_one(H1D_trackEta_primaryPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_zoom, texEtaMC, texTrackPurity, texCollisionDataInfo, drawnZoom_prim, legendPlacementAuto, contextPlacementAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Purity_Eta_DatasetComparison" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_split == true) {
      cout<<"test5"<<endl;
      Draw_TH1_Histograms_in_one(H1D_trackEta_splitPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_split, texEtaMC, texTrackPurity, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
      Draw_TH1_Histograms_in_one(H1D_trackEta_splitPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_split_zoom, texEtaMC, texTrackPurity, texCollisionDataInfo, drawnZoom_split, legendPlacementAuto, contextPlacementAuto, "efficiency");
    }
    else {
      cout << "Divide failed in Draw_Purity_Eta_DatasetComparison" << endl;
    }
  }
}


void Draw_Purity_Phi_DatasetComparison(float* etaRange, bool useSplit) {
  TH3D* H3D_trackPtEtaPhi_primaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_secondaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_splitPrimaryTracks[nDatasets];

  TH1D* H1D_trackPhi_primary[nDatasets];
  TH1D* H1D_trackPhi_secondary[nDatasets];
  TH1D* H1D_trackPhi_splitPrimaryTracks[nDatasets];
  
  TH1D* H1D_trackPhi_primaryPurity[nDatasets];
  TH1D* H1D_primaryPurity_denominator[nDatasets];
  TH1D* H1D_trackPhi_splitPurity[nDatasets];
  TH1D* H1D_splitPurity_denominator[nDatasets];

  bool divideSuccess = false;
  bool divideSuccess_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_primaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_primary"))->Clone("Draw_Purity_Eta_DatasetComparison_associatedtrack_primary"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_secondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"))->Clone("Draw_Purity_Eta_DatasetComparison_associatedtrack_nonprimary"+Datasets[iDataset]);
    if (useSplit == true) {
      H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"))->Clone("Draw_Purity_Eta_DatasetComparison_associatedtrack_splitPrimary"+Datasets[iDataset]);
    }
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon + deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon - deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_secondaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[0] + GLOBAL_epsilon);
    int ibinEta_high_track = H3D_trackPtEtaPhi_secondaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[1] - GLOBAL_epsilon);

    H1D_trackPhi_primary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_primaryTracks[iDataset]->ProjectionZ("trackPhi_primary"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    H1D_trackPhi_secondary[iDataset] = (TH1D*)H3D_trackPtEtaPhi_secondaryTracks[iDataset]->ProjectionZ("trackPhi_secondary"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPhi_splitPrimaryTracks[iDataset] = (TH1D*)H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset]->ProjectionZ("trackPhi_splitPrimaryTracks"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 0, -1, "e");
    }
    // naming the histogram and resetting it to have a chosen name
    H1D_trackPhi_primaryPurity[iDataset] = (TH1D*)H1D_trackPhi_primary[iDataset]->Clone("trackPhi_primaryPurity"+Datasets[iDataset]);
    H1D_trackPhi_primaryPurity[iDataset]->Reset("M");

    H1D_primaryPurity_denominator[iDataset] = (TH1D*)H1D_trackPhi_primary[iDataset]->Clone("trackPhi_primaryPurity_denominator"+Datasets[iDataset]);
    H1D_primaryPurity_denominator[iDataset]->Add(H1D_trackPhi_secondary[iDataset],1.);
    divideSuccess = H1D_trackPhi_primaryPurity[iDataset]->Divide(H1D_trackPhi_primary[iDataset], H1D_primaryPurity_denominator[iDataset]);

    // naming the histogram and resetting it to have a chosen name
    if (useSplit == true) {
      H1D_trackPhi_splitPurity[iDataset] = (TH1D*)H1D_trackPhi_primary[iDataset]->Clone("trackPhi_splitPurity"+Datasets[iDataset]);
      H1D_trackPhi_splitPurity[iDataset]->Reset("M");

      H1D_splitPurity_denominator[iDataset] = (TH1D*)H1D_trackPhi_primary[iDataset]->Clone("trackPhi_splitPurity_denominator"+Datasets[iDataset]);
      H1D_splitPurity_denominator[iDataset]->Add(H1D_trackPhi_splitPrimaryTracks[iDataset],1.);
      divideSuccess_split = H1D_trackPhi_splitPurity[iDataset]->Divide(H1D_trackPhi_primary[iDataset], H1D_splitPurity_denominator[iDataset]);
    }
  }

  TString* pdfNameEntriesNorm = new TString("track_Phi_primaryPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_zoom = new TString("track_Phi_primaryPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_zoom");
  TString* pdfNameEntriesNorm_split = new TString("track_Phi_splitPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfNameEntriesNorm_split_zoom = new TString("track_Phi_splitPurity"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_zoom");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));


  std::array<std::array<float, 2>, 2> drawnZoom_prim = {{{(float)H1D_trackPhi_primaryPurity[0]->GetXaxis()->GetXmin(), (float)H1D_trackPhi_primaryPurity[0]->GetXaxis()->GetXmax()}, 
                                                    {0.9, 1}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 2>, 2> drawnZoom_split = {{{(float)H1D_trackPhi_primaryPurity[0]->GetXaxis()->GetXmin(), (float)H1D_trackPhi_primaryPurity[0]->GetXaxis()->GetXmax()}, 
                                                    {0.995, 1.005}}}; // {{xmin, xmax}, {ymin, ymax}}

  if (divideSuccess == true) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_primaryPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm, texPhiRec, texTrackPurity, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
    Draw_TH1_Histograms_in_one(H1D_trackPhi_primaryPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_zoom, texPhiRec, texTrackPurity, texCollisionDataInfo, drawnZoom_prim, legendPlacementAuto, contextPlacementAuto, "efficiency");
  }
  else {
    cout << "Divide failed in Draw_Purity_Phi_DatasetComparison" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_split == true) {
      Draw_TH1_Histograms_in_one(H1D_trackPhi_splitPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_split, texPhiRec, texTrackPurity, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "efficiency");
      Draw_TH1_Histograms_in_one(H1D_trackPhi_splitPurity, DatasetsNames, nDatasets, textContext, pdfNameEntriesNorm_split_zoom, texPhiRec, texTrackPurity, texCollisionDataInfo, drawnZoom_split, legendPlacementAuto, contextPlacementAuto, "efficiency");
    }
    else {
      cout << "Divide failed in Draw_Purity_Phi_DatasetComparison" << endl;
    }
  }
}



void Draw_Efficiency_Pt_ratio_etaNeg_etaPos_DatasetComparison(float* etaRange, bool useSplit) {
  TH3D* H3D_trackPtEtaPhi_mcparticles[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_assoctracks_split[nDatasets];

  TH1D* H1D_trackPt_mcparticles_pos[nDatasets];

  TH1D* H1D_trackPt_assoctracks_pos[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split_pos[nDatasets];
  
  TH1D* H1D_trackPt_efficiency_pos[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected_pos[nDatasets];


  TH1D* H1D_trackPt_mcparticles_neg[nDatasets];

  TH1D* H1D_trackPt_assoctracks_neg[nDatasets];
  TH1D* H1D_trackPt_assoctracks_split_neg[nDatasets];
  
  TH1D* H1D_trackPt_efficiency_neg[nDatasets];
  TH1D* H1D_trackPt_efficiency_splitCorrected_neg[nDatasets];

  TH1D* H1D_efficiency_etaRightLeftRatio[nDatasets];
  TH1D* H1D_efficiencySplitCorrected_etaRightLeftRatio[nDatasets];


  bool divideSuccess = false;
  bool divideSuccess_split = false;
  bool divideSuccess_etaNegPos = false;
  bool divideSuccess_etaNegPos_splitCorrected = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_mcparticles[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Efficiency_Pt_DatasetComparison_mcparticles"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_assoctracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelColl"+Datasets[iDataset]);
    if (useSplit == true) {
      H3D_trackPtEtaPhi_assoctracks_split[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_associatedtrack_split_primary"))->Clone("Draw_Efficiency_Pt_DatasetComparison_associatedtrackSelCollSplit"+Datasets[iDataset]);
    }
    int ibinEta_low_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0] + deltaEtaMcVsTrackEfficiency);
    int ibinEta_zero_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(0. + GLOBAL_epsilon + deltaEtaMcVsTrackEfficiency);
    int ibinEta_high_mcpart = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1] -deltaEtaMcVsTrackEfficiency);
    int ibinEta_low_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_zero_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(0. + GLOBAL_epsilon );
    int ibinEta_high_track = H3D_trackPtEtaPhi_mcparticles[iDataset]->GetYaxis()->FindBin(etaRange[1]);

    //// positive etas
    H1D_trackPt_mcparticles_pos[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionX("trackPt_mcparticles_pos"+Datasets[iDataset], ibinEta_zero_mcpart, ibinEta_high_mcpart, 0, -1, "e");
    H1D_trackPt_assoctracks_pos[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionX("trackPt_tracks_split_pos"+Datasets[iDataset], ibinEta_zero_track, ibinEta_high_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPt_assoctracks_split_pos[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionX("trackPt_assoctracks_split_pos"+Datasets[iDataset], ibinEta_zero_track, ibinEta_high_track, 0, -1, "e");
    }

    // naming the histogram and resetting it to have a chosen name
    H1D_trackPt_efficiency_pos[iDataset] = (TH1D*)H1D_trackPt_mcparticles_pos[iDataset]->Clone("trackPt_efficiency_pos"+Datasets[iDataset]);
    H1D_trackPt_efficiency_pos[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_efficiency_pos[iDataset]->Divide(H1D_trackPt_assoctracks_pos[iDataset], H1D_trackPt_mcparticles_pos[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    // first naming the histogram and resetting it to have a chosen unique name
    if (useSplit == true) {
      H1D_trackPt_efficiency_splitCorrected_pos[iDataset] = (TH1D*)H1D_trackPt_mcparticles_pos[iDataset]->Clone("trackPt_efficiency_split_pos"+Datasets[iDataset]);
      H1D_trackPt_efficiency_splitCorrected_pos[iDataset]->Reset("M");
      divideSuccess_split = H1D_trackPt_efficiency_splitCorrected_pos[iDataset]->Divide(H1D_trackPt_assoctracks_split_pos[iDataset], H1D_trackPt_mcparticles_pos[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
      H1D_trackPt_efficiency_splitCorrected_pos[iDataset]->Scale(-1.);
      H1D_trackPt_efficiency_splitCorrected_pos[iDataset]->Add(H1D_trackPt_efficiency_pos[iDataset],1.);
    }

    //// negative etas
    H1D_trackPt_mcparticles_neg[iDataset] = (TH1D*)H3D_trackPtEtaPhi_mcparticles[iDataset]->ProjectionX("trackPt_mcparticles_neg"+Datasets[iDataset], ibinEta_low_mcpart, ibinEta_zero_mcpart, 0, -1, "e");
    H1D_trackPt_assoctracks_neg[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks[iDataset]->ProjectionX("trackPt_tracks_split_neg"+Datasets[iDataset], ibinEta_low_track, ibinEta_zero_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPt_assoctracks_split_neg[iDataset] = (TH1D*)H3D_trackPtEtaPhi_assoctracks_split[iDataset]->ProjectionX("trackPt_assoctracks_split_neg"+Datasets[iDataset], ibinEta_low_track, ibinEta_zero_track, 0, -1, "e");
    }

    // naming the histogram and resetting it to have a chosen name
    H1D_trackPt_efficiency_neg[iDataset] = (TH1D*)H1D_trackPt_mcparticles_neg[iDataset]->Clone("trackPt_efficiency_neg"+Datasets[iDataset]);
    H1D_trackPt_efficiency_neg[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_efficiency_neg[iDataset]->Divide(H1D_trackPt_assoctracks_neg[iDataset], H1D_trackPt_mcparticles_neg[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    // first naming the histogram and resetting it to have a chosen unique name
    if (useSplit == true) {
      H1D_trackPt_efficiency_splitCorrected_neg[iDataset] = (TH1D*)H1D_trackPt_mcparticles_neg[iDataset]->Clone("trackPt_efficiency_split_neg"+Datasets[iDataset]);
      H1D_trackPt_efficiency_splitCorrected_neg[iDataset]->Reset("M");
      divideSuccess_split = H1D_trackPt_efficiency_splitCorrected_neg[iDataset]->Divide(H1D_trackPt_assoctracks_split_neg[iDataset], H1D_trackPt_mcparticles_neg[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
      H1D_trackPt_efficiency_splitCorrected_neg[iDataset]->Scale(-1.);
      H1D_trackPt_efficiency_splitCorrected_neg[iDataset]->Add(H1D_trackPt_efficiency_neg[iDataset],1.);
    }

    //// ratio eta right / eta left
    H1D_efficiency_etaRightLeftRatio[iDataset] = (TH1D*)H1D_trackPt_efficiency_neg[iDataset]->Clone("H1D_efficiency_etaRightLeftRatio"+Datasets[iDataset]);
    H1D_efficiency_etaRightLeftRatio[iDataset]->Reset("M");
    divideSuccess_etaNegPos = H1D_efficiency_etaRightLeftRatio[iDataset]->Divide(H1D_trackPt_efficiency_pos[iDataset], H1D_trackPt_efficiency_neg[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency

    if (useSplit == true) {
      H1D_efficiencySplitCorrected_etaRightLeftRatio[iDataset] = (TH1D*)H1D_trackPt_efficiency_splitCorrected_neg[iDataset]->Clone("H1D_efficiencySplitCorrected_etaRightLeftRatio"+Datasets[iDataset]);
      H1D_efficiencySplitCorrected_etaRightLeftRatio[iDataset]->Reset("M");
      divideSuccess_etaNegPos_splitCorrected = H1D_efficiencySplitCorrected_etaRightLeftRatio[iDataset]->Divide(H1D_trackPt_efficiency_splitCorrected_pos[iDataset], H1D_trackPt_efficiency_splitCorrected_neg[iDataset], 1., 1., "b"); // option b for binomial because efficiency: https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSTRANGENESSEfficiency
    }
  }
  TString* pdfName = new TString("track_Pt_efficiency_etaRightLeftRatio"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfName_split = new TString("track_Pt_efficiency_etaRightLeftRatio"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]_splitCorrected");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));

  if (divideSuccess_etaNegPos == true) {
    Draw_TH1_Histograms_in_one(H1D_efficiency_etaRightLeftRatio, DatasetsNames, nDatasets, textContext, pdfName, texPtMC, texTrackEfficiencyRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,standardratio,150MevLine");
  }
  else {
    cout << "Divide failed in Draw_Efficiency_Pt_ratio_etaNeg_etaPos_DatasetComparison: divideSuccess_etaNegPos" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_etaNegPos_splitCorrected == true) {
      Draw_TH1_Histograms_in_one(H1D_efficiencySplitCorrected_etaRightLeftRatio, DatasetsNames, nDatasets, textContext, pdfName_split, texPtMC, texTrackEfficiencyRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,standardratio,150MevLine");
    }
    else {
      cout << "Divide failed in Draw_Efficiency_Pt_ratio_etaNeg_etaPos_DatasetComparison: divideSuccess_etaNegPos_splitCorrected" << endl;
    }
  }  
}



void Draw_Purity_Pt_ratio_etaNeg_etaPos_DatasetComparison(float* etaRange, bool useSplit) {
  TH3D* H3D_trackPtEtaPhi_primaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_secondaryTracks[nDatasets];
  TH3D* H3D_trackPtEtaPhi_splitPrimaryTracks[nDatasets];

  TH1D* H1D_trackPt_primary_pos[nDatasets];
  TH1D* H1D_trackPt_secondary_pos[nDatasets];
  TH1D* H1D_trackPt_splitPrimaryTracks_pos[nDatasets];
  
  TH1D* H1D_trackPt_primaryPurity_pos[nDatasets];
  TH1D* H1D_primaryPurity_denominator_pos[nDatasets];
  TH1D* H1D_trackPt_splitPurity_pos[nDatasets];
  TH1D* H1D_splitPurity_denominator_pos[nDatasets];

  TH1D* H1D_trackPt_primary_neg[nDatasets];
  TH1D* H1D_trackPt_secondary_neg[nDatasets];
  TH1D* H1D_trackPt_splitPrimaryTracks_neg[nDatasets];
  
  TH1D* H1D_trackPt_primaryPurity_neg[nDatasets];
  TH1D* H1D_primaryPurity_denominator_neg[nDatasets];
  TH1D* H1D_trackPt_splitPurity_neg[nDatasets];
  TH1D* H1D_splitPurity_denominator_neg[nDatasets];

  TH1D* H1D_primaryPurityPt_etaRightLeftRatio[nDatasets];
  TH1D* H1D_splitPurityPt_etaRightLeftRatio[nDatasets];


  bool divideSuccess = false;
  bool divideSuccess_split = false;
  bool divideSuccess_etaNegPos_prim = false;
  bool divideSuccess_etaNegPos_split = false;

  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_trackPtEtaPhi_primaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_primary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_primary"+Datasets[iDataset]);
    H3D_trackPtEtaPhi_secondaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_nonprimary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_nonprimary"+Datasets[iDataset]);
    if (useSplit == true) {
      H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_track_pt_track_eta_track_phi_associatedtrack_split_primary"))->Clone("Draw_Purity_Pt_DatasetComparison_associatedtrack_splitPrimary"+Datasets[iDataset]);
    }

    int ibinEta_low_track = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_zero_track = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(0.);
    int ibinEta_high_track = H3D_trackPtEtaPhi_primaryTracks[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    cout << "ibinEta_low_track, ibinEta_zero_track, ibinEta_high_track" << ibinEta_low_track << ", " << ibinEta_zero_track << ", " << ibinEta_high_track << endl;

    // positive etas
    H1D_trackPt_primary_pos[iDataset] = (TH1D*)H3D_trackPtEtaPhi_primaryTracks[iDataset]->ProjectionX("trackPt_primary_pos"+Datasets[iDataset], ibinEta_low_track, ibinEta_zero_track, 0, -1, "e");
    H1D_trackPt_secondary_pos[iDataset] = (TH1D*)H3D_trackPtEtaPhi_secondaryTracks[iDataset]->ProjectionX("trackPt_secondary_pos"+Datasets[iDataset], ibinEta_low_track, ibinEta_zero_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPt_splitPrimaryTracks_pos[iDataset] = (TH1D*)H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset]->ProjectionX("trackPt_splitPrimaryTracks_pos"+Datasets[iDataset], ibinEta_low_track, ibinEta_zero_track, 0, -1, "e");
    }

    // naming the histogram and resetting it to have a chosen name
    H1D_trackPt_primaryPurity_pos[iDataset] = (TH1D*)H1D_trackPt_primary_pos[iDataset]->Clone("trackPt_primaryPurity_pos"+Datasets[iDataset]);
    H1D_trackPt_primaryPurity_pos[iDataset]->Reset("M");

    H1D_primaryPurity_denominator_pos[iDataset] = (TH1D*)H1D_trackPt_primary_pos[iDataset]->Clone("trackPt_primaryPurity_denominator_pos"+Datasets[iDataset]);
    H1D_primaryPurity_denominator_pos[iDataset]->Add(H1D_trackPt_secondary_pos[iDataset],1.);
    divideSuccess = H1D_trackPt_primaryPurity_pos[iDataset]->Divide(H1D_trackPt_primary_pos[iDataset], H1D_primaryPurity_denominator_pos[iDataset]);

    // naming the histogram and resetting it to have a chosen name
    if (useSplit == true) {
      H1D_trackPt_splitPurity_pos[iDataset] = (TH1D*)H1D_trackPt_primary_pos[iDataset]->Clone("trackPt_splitPurity_pos"+Datasets[iDataset]);
      H1D_trackPt_splitPurity_pos[iDataset]->Reset("M");

      H1D_splitPurity_denominator_pos[iDataset] = (TH1D*)H1D_trackPt_primary_pos[iDataset]->Clone("trackPt_splitPurity_denominator_pos"+Datasets[iDataset]);
      H1D_splitPurity_denominator_pos[iDataset]->Add(H1D_trackPt_splitPrimaryTracks_pos[iDataset],1.);
      divideSuccess_split = H1D_trackPt_splitPurity_pos[iDataset]->Divide(H1D_trackPt_primary_pos[iDataset], H1D_splitPurity_denominator_pos[iDataset]);
    }
    // negative etas
    H1D_trackPt_primary_neg[iDataset] = (TH1D*)H3D_trackPtEtaPhi_primaryTracks[iDataset]->ProjectionX("trackPt_primary_neg"+Datasets[iDataset], ibinEta_zero_track, ibinEta_high_track, 0, -1, "e");
    H1D_trackPt_secondary_neg[iDataset] = (TH1D*)H3D_trackPtEtaPhi_secondaryTracks[iDataset]->ProjectionX("trackPt_secondary_neg"+Datasets[iDataset], ibinEta_zero_track, ibinEta_high_track, 0, -1, "e");
    if (useSplit == true) {
      H1D_trackPt_splitPrimaryTracks_neg[iDataset] = (TH1D*)H3D_trackPtEtaPhi_splitPrimaryTracks[iDataset]->ProjectionX("trackPt_splitPrimaryTracks_neg"+Datasets[iDataset], ibinEta_zero_track, ibinEta_high_track, 0, -1, "e");
    }
    // naming the histogram and resetting it to have a chosen name
    H1D_trackPt_primaryPurity_neg[iDataset] = (TH1D*)H1D_trackPt_primary_neg[iDataset]->Clone("trackPt_primaryPurity_neg"+Datasets[iDataset]);
    H1D_trackPt_primaryPurity_neg[iDataset]->Reset("M");

    H1D_primaryPurity_denominator_neg[iDataset] = (TH1D*)H1D_trackPt_primary_neg[iDataset]->Clone("trackPt_primaryPurity_denominator_neg"+Datasets[iDataset]);
    H1D_primaryPurity_denominator_neg[iDataset]->Add(H1D_trackPt_secondary_neg[iDataset],1.);
    divideSuccess = H1D_trackPt_primaryPurity_neg[iDataset]->Divide(H1D_trackPt_primary_neg[iDataset], H1D_primaryPurity_denominator_neg[iDataset]);

    // naming the histogram and resetting it to have a chosen name
    if (useSplit == true) {
      H1D_trackPt_splitPurity_neg[iDataset] = (TH1D*)H1D_trackPt_primary_neg[iDataset]->Clone("trackPt_splitPurity_neg"+Datasets[iDataset]);
      H1D_trackPt_splitPurity_neg[iDataset]->Reset("M");

      H1D_splitPurity_denominator_neg[iDataset] = (TH1D*)H1D_trackPt_primary_neg[iDataset]->Clone("trackPt_splitPurity_denominator_neg"+Datasets[iDataset]);
      H1D_splitPurity_denominator_neg[iDataset]->Add(H1D_trackPt_splitPrimaryTracks_neg[iDataset],1.);
      divideSuccess_split = H1D_trackPt_splitPurity_neg[iDataset]->Divide(H1D_trackPt_primary_neg[iDataset], H1D_splitPurity_denominator_neg[iDataset]);
    }
    //// ratio eta right / eta left
    H1D_primaryPurityPt_etaRightLeftRatio[iDataset] = (TH1D*)H1D_primaryPurity_denominator_neg[iDataset]->Clone("H1D_primaryPurityPt_etaRightLeftRatio"+Datasets[iDataset]);
    H1D_primaryPurityPt_etaRightLeftRatio[iDataset]->Reset("M");
    divideSuccess_etaNegPos_prim = H1D_primaryPurityPt_etaRightLeftRatio[iDataset]->Divide(H1D_trackPt_primaryPurity_pos[iDataset], H1D_trackPt_primaryPurity_neg[iDataset]);
    
    if (useSplit == true) {
      H1D_splitPurityPt_etaRightLeftRatio[iDataset] = (TH1D*)H1D_splitPurity_denominator_neg[iDataset]->Clone("H1D_splitPurityPt_etaRightLeftRatio"+Datasets[iDataset]);
      H1D_splitPurityPt_etaRightLeftRatio[iDataset]->Reset("M");
      divideSuccess_etaNegPos_split = H1D_splitPurityPt_etaRightLeftRatio[iDataset]->Divide(H1D_trackPt_splitPurity_pos[iDataset], H1D_trackPt_splitPurity_neg[iDataset]);
    }
  }
  TString* pdfName = new TString("track_Pt_primaryPurity_etaRightLeftRatio"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");
  TString* pdfName_split = new TString("track_Pt_splitPurity_etaRightLeftRatio"+dummyName[0]+"_@eta["+Form("%.1f", etaRange[0])+","+Form("%.1f", etaRange[1])+"]");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, contextEtaRange(etaRange), ""));

  if (divideSuccess_etaNegPos_prim == true) {
    Draw_TH1_Histograms_in_one(H1D_primaryPurityPt_etaRightLeftRatio, DatasetsNames, nDatasets, textContext, pdfName, texPtRec, texTrackPurityRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,standardratio,150MevLine");
  }
  else {
    cout << "Divide failed in Draw_Purity_Pt_ratio_etaNeg_etaPos_DatasetComparison: divideSuccess_etaNegPos_prim" << endl;
  }
  if (useSplit == true) {
    if (divideSuccess_etaNegPos_split == true) {
      Draw_TH1_Histograms_in_one(H1D_splitPurityPt_etaRightLeftRatio, DatasetsNames, nDatasets, textContext, pdfName_split, texPtRec, texTrackPurityRatioEtaComparison, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,standardratio,150MevLine");
    }
    else {
     cout << "Divide failed in Draw_Purity_Pt_ratio_etaNeg_etaPos_DatasetComparison: divideSuccess_etaNegPos_split" << endl;
    }
  }
}




void Draw_Pt_gen_DatasetComparison(float* etaRange, std::string options) {
  TH3D* H3D_ptetaphi_track[nDatasets];

  TH1D* H1D_trackPt[nDatasets];
  TH1D* H1D_trackPt_rebinned[nDatasets];
  
  TH1D* H1D_trackPt_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_ptetaphi_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]);
    H3D_ptetaphi_track[iDataset]->Reset("M");
    if (options.find("primaries") != std::string::npos) {
      H3D_ptetaphi_track[iDataset]->Add((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"));
    }
    if (options.find("secondaries") != std::string::npos) {
      H3D_ptetaphi_track[iDataset]->Add((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary"));
    }


    int ibinEta_low_track = H3D_ptetaphi_track[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_ptetaphi_track[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    
    H1D_trackPt[iDataset] = (TH1D*)H3D_ptetaphi_track[iDataset]->ProjectionX("trackPt_"+Datasets[iDataset], ibinEta_low_track, ibinEta_high_track, 1, H3D_ptetaphi_track[iDataset]->GetNbinsZ(), "e");

    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(2.,"trackPt_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);

    if (isDatasetWeighted[iDataset]) {
      Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
    } else {
      Nevents = GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]);
    }
    NormaliseYieldToNEvents(H1D_trackPt_rebinned[iDataset], Nevents);
    cout << "Dataset " << iDataset << ": Nevents = " << Nevents << endl;

    H1D_trackPt_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_trackPt_rebinned[0]);
  }

  TString trackComposition = "";
  if (options.find("primaries") != std::string::npos) {
    trackComposition = (TString)"primaries";
  }
  if (options.find("secondaries") != std::string::npos) {
    trackComposition = (TString)"secondaries";
  }
  if (options.find("primaries") != std::string::npos && options.find("secondaries") != std::string::npos) {
    trackComposition = (TString)"primaries+secondaries";
  }

  TString* pdfName = new TString("track_Pt_gen_DataComp"+trackComposition);
  TString* pdfName_ratio = new TString("track_Pt_gen_DataComp_ratio"+trackComposition);

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtMC, texTrackPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
  if (divideSuccess == true && options.find("ratio") != std::string::npos) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtMC, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoratio,noMarkerFirst");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}


void Draw_Eta_gen_DatasetComparison(float* ptRange, std::string options) {
  TH3D* H3D_ptetaphi_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){


    H3D_ptetaphi_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]);
    H3D_ptetaphi_track[iDataset]->Reset("M");
    if (options.find("primaries") != std::string::npos) {
      H3D_ptetaphi_track[iDataset]->Add((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"));
    }
    if (options.find("secondaries") != std::string::npos) {
      H3D_ptetaphi_track[iDataset]->Add((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary"));
    }

    int ibinPt_low_track = H3D_ptetaphi_track[iDataset]->GetXaxis()->FindBin(ptRange[0]);
    int ibinPt_high_track = H3D_ptetaphi_track[iDataset]->GetXaxis()->FindBin(ptRange[1]);

    H1D_trackEta[iDataset] = (TH1D*)H3D_ptetaphi_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset], ibinPt_low_track, ibinPt_high_track, 1, H3D_ptetaphi_track[iDataset]->GetNbinsZ(), "e");

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(5.,"trackEta_rebinned"+Datasets[iDataset]);


    if (isDatasetWeighted[iDataset]) {
      Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
    } else {
      Nevents = GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]);
    }
    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iDataset], Nevents);

    H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0]);
  }

  TString trackComposition = "";
  if (options.find("primaries") != std::string::npos) {
    trackComposition = (TString)"primaries";
  }
  if (options.find("secondaries") != std::string::npos) {
    trackComposition = (TString)"secondaries";
  }
  if (options.find("primaries") != std::string::npos && options.find("secondaries") != std::string::npos) {
    trackComposition = (TString)"primaries+secondaries";
  }

  TString* pdfNameEventNorm = new TString("track_Eta_gen_DataComp_EventNorm"+trackComposition);
  TString* pdfNameEventNorm_ratio = new TString("track_Eta_gen_DataComp_EventNorm_ratio"+trackComposition);


  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfNameEventNorm, texEtaMC, texTrackEtaYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  if (divideSuccess == true && options.find("ratio") != std::string::npos) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfNameEventNorm_ratio, texEtaMC, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoratio,noMarkerFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}


void Draw_Phi_gen_DatasetComparison(float* ptRange, float* etaRange, std::string options) { 
  TH3D* H3D_ptetaphi_track[nDatasets];

  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_ptetaphi_track[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]);
    H3D_ptetaphi_track[iDataset]->Reset("M");
    if (options.find("primaries") != std::string::npos) {
      H3D_ptetaphi_track[iDataset]->Add((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest"));
    }
    if (options.find("secondaries") != std::string::npos) {
      H3D_ptetaphi_track[iDataset]->Add((TH3D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h3_particle_pt_particle_eta_particle_phi_mcpart_nonprimary"));
    }

    int ibinPt_low_track = H3D_ptetaphi_track[iDataset]->GetXaxis()->FindBin(ptRange[0]);
    int ibinPt_high_track = H3D_ptetaphi_track[iDataset]->GetXaxis()->FindBin(ptRange[1]);
    int ibinEta_low_track = H3D_ptetaphi_track[iDataset]->GetYaxis()->FindBin(etaRange[0]);
    int ibinEta_high_track = H3D_ptetaphi_track[iDataset]->GetYaxis()->FindBin(etaRange[1]);
    
    H1D_trackPhi[iDataset] = (TH1D*)H3D_ptetaphi_track[iDataset]->ProjectionZ("trackPhi_"+Datasets[iDataset], ibinPt_low_track, ibinPt_high_track, ibinEta_low_track, ibinEta_high_track, "e");
    H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(5.,"trackPhi_rebinned_"+Datasets[iDataset]);
    cout << "H1D_trackPhi_rebinned[iDataset]->GetEntries() preNorm = " << H1D_trackPhi_rebinned[iDataset]->Integral() << endl;

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    if (isDatasetWeighted[iDataset]) {
      Nevents = GetNEventsSelected_JetFramework_weighted(file_O2Analysis_list[iDataset]);
    } else {
      Nevents = GetNEventsSelected_TrackEffWorkflow(file_O2Analysis_list[iDataset]);
    }
    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iDataset], Nevents);
    cout << "H1D_trackPhi_rebinned[iDataset]->GetEntries() postNorm = " << H1D_trackPhi_rebinned[iDataset]->Integral() << ", Nevents = " << Nevents << endl;

    H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[iDataset], H1D_trackPhi_rebinned[0]);
  }

  TString trackComposition = "";
  if (options.find("primaries") != std::string::npos) {
    trackComposition = (TString)"_primaries";
  }
  if (options.find("secondaries") != std::string::npos) {
    trackComposition = (TString)"_secondaries";
  }
  if (options.find("primaries") != std::string::npos && options.find("secondaries") != std::string::npos) {
    trackComposition = (TString)"_primaries+secondaries";
  }

  TString* pdfName = new TString("track_Phi_gen_DataComp"+trackComposition);
  TString* pdfName_ratio = new TString("track_Phi_gen_DataComp_ratio"+trackComposition);

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiMC, texTrackPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  if (divideSuccess == true && options.find("ratio") != std::string::npos) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiMC, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoratio,noMarkerFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}









void Draw_Pt_gen_DatasetComparison_H2CentVersion(std::string options) {
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackPt[nDatasets];
  TH1D* H1D_trackPt_rebinned[nDatasets];
  
  TH1D* H1D_trackPt_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_particle_pt"))->Clone("Draw_Pt_DatasetComparison"+Datasets[iDataset]);
    
    H1D_trackPt[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPt_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");


    H1D_trackPt_rebinned[iDataset] = (TH1D*)H1D_trackPt[iDataset]->Rebin(2.,"trackPt_rebinned_"+Datasets[iDataset]);

    // NormaliseYieldToNEntries(H1D_trackPt_rebinned[iDataset]);

    if (isDatasetWeighted[iDataset]) {
      Nevents = GetNEventsSelected_TrackEffWorkflow_gen_weighted(file_O2Analysis_list[iDataset]);
    } else {
      Nevents = GetNEventsSelected_TrackEffWorkflow_gen(file_O2Analysis_list[iDataset]);
    }
    NormaliseYieldToNEvents(H1D_trackPt_rebinned[iDataset], Nevents);


    H1D_trackPt_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPt_rebinned[iDataset]->Clone("trackPt_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPt_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPt_rebinned_ratios[iDataset]->Divide(H1D_trackPt_rebinned[iDataset], H1D_trackPt_rebinned[0]);
  }

  TString* pdfName = new TString("track_Pt_gen_DataComp");
  TString* pdfName_ratio = new TString("track_Pt_gen_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPtMC, texTrackPtYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "logx,logy");
  if (divideSuccess == true && options.find("ratio") != std::string::npos) {
    Draw_TH1_Histograms_in_one(H1D_trackPt_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPtMC, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoratio,noMarkerFirst,logx");
  }
  else {
    cout << "Divide failed in Draw_Pt_DatasetComparison" << endl;
  }
}


void Draw_Eta_gen_DatasetComparison_H2CentVersion(std::string options) {
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackEta[nDatasets];
  TH1D* H1D_trackEta_rebinned[nDatasets];
  
  TH1D* H1D_trackEta_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents; 
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){


    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_particle_eta"))->Clone("Draw_Eta_DatasetComparison"+Datasets[iDataset]);

    H1D_trackEta[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackEta_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");

    H1D_trackEta_rebinned[iDataset] = (TH1D*)H1D_trackEta[iDataset]->Rebin(5.,"trackEta_rebinned"+Datasets[iDataset]);


    if (isDatasetWeighted[iDataset]) {
      Nevents = GetNEventsSelected_TrackEffWorkflow_gen_weighted(file_O2Analysis_list[iDataset]);
    } else {
      Nevents = GetNEventsSelected_TrackEffWorkflow_gen(file_O2Analysis_list[iDataset]);
    }
    NormaliseYieldToNEvents(H1D_trackEta_rebinned[iDataset], Nevents);

    H1D_trackEta_rebinned_ratios[iDataset] = (TH1D*)H1D_trackEta_rebinned[iDataset]->Clone("trackEta_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackEta_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackEta_rebinned_ratios[iDataset]->Divide(H1D_trackEta_rebinned[iDataset], H1D_trackEta_rebinned[0]);
  }

  TString* pdfNameEventNorm = new TString("track_Eta_gen_DataComp_EventNorm");
  TString* pdfNameEventNorm_ratio = new TString("track_Eta_gen_DataComp_EventNorm_ratio");


  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned, DatasetsNames, nDatasets, textContext, pdfNameEventNorm, texEtaMC, texTrackEtaYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  if (divideSuccess == true && options.find("ratio") != std::string::npos) {
    Draw_TH1_Histograms_in_one(H1D_trackEta_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfNameEventNorm_ratio, texEtaMC, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoratio,noMarkerFirst");
  }
  else {
    cout << "Divide failed in Draw_Eta_DatasetComparison" << endl;
  }
}


void Draw_Phi_gen_DatasetComparison_H2CentVersion(std::string options) { 
  TH2D* H2D_centrality_track[nDatasets];

  TH1D* H1D_trackPhi[nDatasets];
  TH1D* H1D_trackPhi_rebinned[nDatasets];
  
  TH1D* H1D_trackPhi_rebinned_ratios[nDatasets];

  bool divideSuccess = false;

  double Nevents;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H2D_centrality_track[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_centrality_particle_phi"))->Clone("Draw_Phi_DatasetComparison"+Datasets[iDataset]);

    H1D_trackPhi[iDataset] = (TH1D*)H2D_centrality_track[iDataset]->ProjectionY("trackPhi_"+Datasets[iDataset], 1, H2D_centrality_track[iDataset]->GetNbinsX(), "e");
    H1D_trackPhi_rebinned[iDataset] = (TH1D*)H1D_trackPhi[iDataset]->Rebin(5.,"trackPhi_rebinned_"+Datasets[iDataset]);
    cout << "H1D_trackPhi_rebinned[iDataset]->GetEntries() preNorm = " << H1D_trackPhi_rebinned[iDataset]->Integral() << endl;

    // NormaliseYieldToNEntries(H1D_trackPhi_rebinned[iDataset]);
    if (isDatasetWeighted[iDataset]) {
      Nevents = GetNEventsSelected_TrackEffWorkflow_gen_weighted(file_O2Analysis_list[iDataset]);
    } else {
      Nevents = GetNEventsSelected_TrackEffWorkflow_gen(file_O2Analysis_list[iDataset]);
    }
    NormaliseYieldToNEvents(H1D_trackPhi_rebinned[iDataset], Nevents);
    cout << "H1D_trackPhi_rebinned[iDataset]->GetEntries() postNorm = " << H1D_trackPhi_rebinned[iDataset]->Integral() << ", Nevents = " << Nevents << endl;

    H1D_trackPhi_rebinned_ratios[iDataset] = (TH1D*)H1D_trackPhi_rebinned[iDataset]->Clone("trackPhi_rebinned_ratios"+Datasets[iDataset]);
    H1D_trackPhi_rebinned_ratios[iDataset]->Reset("M");
    divideSuccess = H1D_trackPhi_rebinned_ratios[iDataset]->Divide(H1D_trackPhi_rebinned[iDataset], H1D_trackPhi_rebinned[0]);
  }

  TString* pdfName = new TString("track_Phi_gen_DataComp");
  TString* pdfName_ratio = new TString("track_Phi_gen_DataComp_ratio");

  TString textContext(contextTrackDatasetComp(""));

  Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned, DatasetsNames, nDatasets, textContext, pdfName, texPhiMC, texTrackPhiYield_EventNorm, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "");
  if (divideSuccess == true && options.find("ratio") != std::string::npos) {
    Draw_TH1_Histograms_in_one(H1D_trackPhi_rebinned_ratios, DatasetsNames, nDatasets, textContext, pdfName_ratio, texPhiMC, texRatioDatasets, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, "autoratio,noMarkerFirst");
  }
  else {
    cout << "Divide failed in Draw_Phi_DatasetComparison" << endl;
  }
}

