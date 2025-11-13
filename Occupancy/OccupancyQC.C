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
#include "TPolyLine.h"

// #include <RooUnfold.h>
// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"

//My Libraries
#include "./OccupancyQC_settings.h"
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
// TString contextDataset1D(int iDataset, float* variableRange, std::string options);
// TString contextDatasetCompAndRadiusAndVarRange(float jetRadius, float* variableRange, std::string options);
// TString contextPtRange(float* PtRange);
// TString contextEtaRange(float* PtRange);
// TString contextJetRadius(float jetRadius);

//////////// QC plot functions


// Dataset comparison
void Draw_Ntracks_vs_Occupancy_DatasetComparison_JetVersion(std::string options);
void Draw_Ntracks_vs_Occupancy_DatasetComparison_Igor(std::string options);

/////////////////////////////////////////////////////
///////////////////// Main Macro ////////////////////
/////////////////////////////////////////////////////

void OccupancyQC() {
  // Load necessary libraries
  LoadLibs();
  // Set the default style
  SetStyle();

  // TString* SaveAs_Title = new TString("");
  TString* texXtitle = new TString("");
  TString* texYtitle = new TString("");
  // TString* Extra = new TString("");


  Draw_Ntracks_vs_Occupancy_DatasetComparison_JetVersion("");
  // Draw_Ntracks_vs_Occupancy_DatasetComparison_Igor("");
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



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////// QC  plot functions /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Draw_Ntracks_vs_Occupancy_DatasetComparison_JetVersion(std::string options) {

  TH2D* H2D_occupancy_ntracksselptetacuts_postcollsel[nDatasets];
  TH2D* H2D_occupancy_ntracksselptetacuts_precollsel[nDatasets];
  TH2D* H2D_occupancy_ntracksselptetacuts_postcollsel_rebinned[nDatasets];
  TH2D* H2D_occupancy_ntracksselptetacuts_precollsel_rebinned[nDatasets];
  TH1D* H1D_occupancy_averageNTracksselptetacuts_postcollsel_withProfile[nDatasets];
  TH1D* H1D_occupancy_averageNTracksselptetacuts_precollsel_withProfile[nDatasets];

  // int nBinsOccupancy, nBinsNTracks, numeratorAverageNTracks, denominatorAverageNTracks, nEventsAtOccupancyNTracks;

  // TString* yAxisLabel;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    cout << "test1" << endl;
    H2D_occupancy_ntracksselptetacuts_postcollsel[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_occupancy_ntracksselptetacuts_postsel"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_ntracksselptetacuts_postcollsel"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_occupancy_ntracksselptetacuts_postcollsel[iDataset]->Sumw2();
    H2D_occupancy_ntracksselptetacuts_precollsel[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_occupancy_ntracksselptetacuts_presel"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_ntracksselptetacuts_precollsel"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_occupancy_ntracksselptetacuts_precollsel[iDataset]->Sumw2();

    // H2D_occupancy_ntrackssel_postcollsel[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_occupancy_ntrackssel_postcollsel"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_ntrackssel_postcollsel"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_occupancy_ntrackssel_postcollsel[iDataset]->Sumw2();
    // H2D_occupancy_ntrackssel_precollsel[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_occupancy_ntrackssel_precollsel"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_ntrackssel_precollsel"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_occupancy_ntrackssel_precollsel[iDataset]->Sumw2();

    // H2D_occupancy_ntracksall_postcollsel[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_occupancy_ntracksall_postcollsel"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_ntracksall_postcollsel"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_occupancy_ntracksall_postcollsel[iDataset]->Sumw2();
    // H2D_occupancy_ntracksall_precollsel[iDataset] = (TH2D*)((TH2D*)file_O2Analysis_list[iDataset]->Get(analysisWorkflow[iDataset]+"/h2_occupancy_ntracksall_precollsel"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_ntracksall_precollsel"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_occupancy_ntracksall_precollsel[iDataset]->Sumw2();


    H2D_occupancy_ntracksselptetacuts_postcollsel_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntracksselptetacuts_postcollsel[iDataset]->RebinX(1.,"H2D_occupancy_ntracksselptetacuts_postcollsel_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H2D_occupancy_ntracksselptetacuts_precollsel_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntracksselptetacuts_precollsel[iDataset]->RebinX(1.,"H2D_occupancy_ntracksselptetacuts_precollsel_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);

    // H2D_occupancy_ntrackssel_postcollsel_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntrackssel_postcollsel[iDataset]->RebinX(5.,"H2D_occupancy_ntrackssel_postcollsel_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_occupancy_ntrackssel_precollsel_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntrackssel_precollsel[iDataset]->RebinX(5.,"H2D_occupancy_ntrackssel_precollsel_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);

    // H2D_occupancy_ntracksall_postcollsel_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntracksall_postcollsel[iDataset]->RebinX(5.,"H2D_occupancy_ntracksall_postcollsel_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);
    // H2D_occupancy_ntracksall_precollsel_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntracksall_precollsel[iDataset]->RebinX(5.,"H2D_occupancy_ntracksall_precollsel_rebinned"+Datasets[iDataset]+DatasetsNames[iDataset]);


    H1D_occupancy_averageNTracksselptetacuts_postcollsel_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntracksselptetacuts_postcollsel_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTracksselptetacuts_postcollsel_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    H1D_occupancy_averageNTracksselptetacuts_precollsel_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntracksselptetacuts_precollsel_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTracksselptetacuts_precollsel_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");

    // H1D_occupancy_averageNTrackssel_postcollsel_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntrackssel_postcollsel_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTrackssel_postcollsel_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    // H1D_occupancy_averageNTrackssel_precollsel_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntrackssel_precollsel_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTrackssel_precollsel_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");

    // H1D_occupancy_averageNTracksall_postcollsel_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntracksall_postcollsel_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTracksall_postcollsel_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    // H1D_occupancy_averageNTracksall_precollsel_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntracksall_precollsel_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTracksall_precollsel_withProfile"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    cout << "test2" << endl;
  }

  TString* pdfName_trackselptetacuts_precollsel = new TString("occupancy_averageNTracks_DataComp_preCollSel");
  TString* pdfName_trackselptetacuts_postcollsel = new TString("occupancy_averageNTracks_DataComp_postCollSel");

  // TString* pdfName_tracksel_precollsel = new TString("occupancy_averageNTracks_DataComp_preCollSel");
  // TString* pdfName_tracksel_postcollsel = new TString("occupancy_averageNTracks_DataComp_postCollSel");

  // TString* pdfName_trackall_precollsel = new TString("occupancy_averageNTracks_DataComp_preCollSel");
  // TString* pdfName_trackall_postcollsel = new TString("occupancy_averageNTracks_DataComp_postCollSel");



  TString* pdfName3 = new TString("occupancy_NTracks_DataComp_JetVersion_postCollSel");
  TString* pdfName3zoom = new TString("occupancy_NTracks_DataComp_JetVersion_postCollSel_zoom");
  TString* pdfName4 = new TString("occupancy_NTracks_DataComp_JetVersion_preCollSel");
  TString* pdfName4zoom = new TString("occupancy_NTracks_DataComp_JetVersion_preCollSel_zoom");
  // TString* pdfName2 = new TString("occupancy_averageNTracks_DataComp_testWithoutProfile");

  // TString textContext(contextDatasetCompAndRadiusAndVarRange(jetRadius, etaRange, "eta"));
  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, (TString)"globalTracks", ""));
  std::array<std::array<float, 2>, 2> drawnWindow = {{{-999, -999}, {0, 2500}}}; // {{xmin, xmax}, {ymin, ymax}}
  std::array<std::array<float, 3>, 3> drawnWindow2D = {{{000, 2000}, {0, 6000}, {-999, -999}}}; // {{xmin, xmax}, {ymin, ymax}}

  const std::array<float, 2> contextPlacementCustom = {{-999, 0.4}}; // {{x_topleft, y_topleft}}

  Draw_TH1_Histograms(H1D_occupancy_averageNTracksselptetacuts_precollsel_withProfile, DatasetsNames, nDatasets, textContext, pdfName_trackselptetacuts_precollsel, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementCustom, ""+histDatasetComparisonStructure);
  Draw_TH1_Histograms(H1D_occupancy_averageNTracksselptetacuts_postcollsel_withProfile, DatasetsNames, nDatasets, textContext, pdfName_trackselptetacuts_postcollsel, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementCustom, ""+histDatasetComparisonStructure);
  Draw_TH2_Histograms(H2D_occupancy_ntracksselptetacuts_postcollsel, DatasetsNames, nDatasets, textContext, pdfName3zoom, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH2_Histograms(H2D_occupancy_ntracksselptetacuts_postcollsel, DatasetsNames, nDatasets, textContext, pdfName3, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH2_Histograms(H2D_occupancy_ntracksselptetacuts_precollsel, DatasetsNames, nDatasets, textContext, pdfName4zoom, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow2D, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH2_Histograms(H2D_occupancy_ntracksselptetacuts_precollsel, DatasetsNames, nDatasets, textContext, pdfName4, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow2D, th2ContoursNone, contourNumberNone, "logz");
}

void Draw_Ntracks_vs_Occupancy_DatasetComparison_Igor(std::string options) {

  TH3D* H3D_occupancy_ntrackssel_postcollsel8[nDatasets];
  TH2D* H2D_occupancy_ntrackssel_postcollsel8[nDatasets];
  TH2D* H2D_occupancy_ntrackssel_postcollsel8_rebinned[nDatasets];
  TH1D* H1D_occupancy_averageNTrackssel_postcollsel8_withProfile[nDatasets];

  int nBinsOccupancy, nBinsNTracks, numeratorAverageNTracks, denominatorAverageNTracks, nEventsAtOccupancyNTracks;

  // TString* yAxisLabel;
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){

    H3D_occupancy_ntrackssel_postcollsel8[iDataset] = (TH3D*)((TH3D*)file_O2Analysis_list[iDataset]->Get("detector-occupancy-qa-task/nTracksGlobal_vs_nPV_vs_occup_pure"))->Clone("Draw_Ntracks_vs_Occupancy_DatasetComparison_Igor"+Datasets[iDataset]+DatasetsNames[iDataset]);
    H3D_occupancy_ntrackssel_postcollsel8[iDataset]->Sumw2();

    H2D_occupancy_ntrackssel_postcollsel8[iDataset] = (TH2D*)H3D_occupancy_ntrackssel_postcollsel8[iDataset]->Project3D((TString)Form("%d", iDataset)+"_Igor_ntracksoccupancy_e_yz"); //can't use letter D in this or it seems to replace the histogram in current pad (see documentation of ProjectionX function. Isn't mentioned in project3D sadly)

    H2D_occupancy_ntrackssel_postcollsel8_rebinned[iDataset] = (TH2D*)H2D_occupancy_ntrackssel_postcollsel8[iDataset]->RebinX(1.,"H2D_occupancy_ntrackssel_postcollsel8_rebinned_Igor"+Datasets[iDataset]+DatasetsNames[iDataset]);

    H1D_occupancy_averageNTrackssel_postcollsel8_withProfile[iDataset] = (TH1D*)H2D_occupancy_ntrackssel_postcollsel8_rebinned[iDataset]->ProfileX("H1D_occupancy_averageNTrackssel_postcollsel8_withProfile_Igor"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");

    // check that TProfile works well
    // H1D_occupancy_averageNTrackssel_postcollsel8[iDataset] = (TH1D*)H2D_occupancy_ntrackssel_postcollsel8[iDataset]->ProjectionX("H1D_occupancy_averageNTrackssel_postcollsel8_ntrackssel_postcollsel8"+Datasets[iDataset]+DatasetsNames[iDataset], 0, -1, "e");
    // H1D_occupancy_averageNTrackssel_postcollsel8[iDataset]->Reset("M");

    // nBinsOccupancy = H2D_occupancy_ntrackssel_postcollsel8[iDataset]->GetNbinsX();
    // nBinsNTracks = H2D_occupancy_ntrackssel_postcollsel8[iDataset]->GetNbinsY();
    // for(int iOccupancy = 0; iOccupancy < nBinsOccupancy; iOccupancy++){
    //   numeratorAverageNTracks = 0;
    //   denominatorAverageNTracks = 0;
    //   for(int iNTracks = 0; iNTracks < nBinsNTracks; iNTracks++){
    //     nEventsAtOccupancyNTracks = H2D_occupancy_ntrackssel_postcollsel8[iDataset]->GetBinContent(iOccupancy, iNTracks);
    //     denominatorAverageNTracks += nEventsAtOccupancyNTracks;
    //     numeratorAverageNTracks += iNTracks * nEventsAtOccupancyNTracks;
    //   }
    //   H1D_occupancy_averageNTrackssel_postcollsel8[iDataset]->SetBinContent(iOccupancy, numeratorAverageNTracks / denominatorAverageNTracks);
    // }

  }

  TString* pdfName = new TString("2D_occupancy_NTracks_DataComp_postcollsel_IgorVERSION");
  TString* pdfNamepostcollsel = new TString("occupancy_averageNTracks_DataComp_postcollsel_IgorVERSION");

  TString textContext(contextCustomTwoFields(*texDatasetsComparisonCommonDenominator, (TString)"globalTracks", ""));
  std::array<std::array<float, 2>, 2> drawnWindow = {{{-999, -999}, {0, 1200}}}; // {{xmin, xmax}, {ymin, ymax}}

  // Draw_TH1_Histograms(H1D_occupancy_averageNTrackssel_postcollsel8, DatasetsNames, nDatasets, textContext, pdfName, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindowAuto, legendPlacementAuto, contextPlacementAuto, ""+histDatasetComparisonStructure);
  Draw_TH2_Histograms(H2D_occupancy_ntrackssel_postcollsel8, DatasetsNames, nDatasets, textContext, pdfName, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow2DAuto, th2ContoursNone, contourNumberNone, "logz");
  Draw_TH1_Histograms(H1D_occupancy_averageNTrackssel_postcollsel8_withProfile, DatasetsNames, nDatasets, textContext, pdfNamepostcollsel, texWeightOccupancy, texMeanNtracks, texCollisionDataInfo, drawnWindow, legendPlacementAuto, contextPlacementAuto, ""+histDatasetComparisonStructure);

}
