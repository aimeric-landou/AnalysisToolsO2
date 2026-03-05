// This is a template. To use the TrackQC.C, rename this file to TrackQC_settings.h and edit it how you want.

// Bin edge control
// float GLOBAL_epsilon = 0.00001;

// Analysis settings
const int nCollSystems = 1;
const TString collSystems[nCollSystems] = {"pp"};

const int iCollSystem = 1; // if 1 starting from 0 instead of 1 

// const int nCentralityBins = 7;
// const float arrayCentralityBinning[nCentralityBins+1] = {0, 10, 20, 30, 50, 70, 100};
const int nCentralityBins = 1;
const float arrayCentralityBinning[nCentralityBins+1] = {-10, 0};


