#ifndef JETSPECTRUM_EFFICIENCYPURITYGETTERS_H
#define JETSPECTRUM_EFFICIENCYPURITYGETTERS_H

bool  Get_Pt_JetEfficiency(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, std::string options, bool controlMC);
bool  Get_Pt_JetEfficiency_fineBinning(TH1D* &H1D_jetEfficiency, int iDataset, int iRadius, std::string options, bool controlMC);
bool Get_ResponseMatrix_Pt_KinematicEffiency(TH1D* &H1D_kinematicEfficiency, TH2D* H2D_jetPtResponseMatrix, TString name_H1D_kinematicEfficiency, int iRadius);
bool Get_Pt_JetFakes(TH1D* &H1D_jetFakes, int iDataset, int iRadius, std::string options, bool controlMC);
bool Get_Pt_JetFakes_fineBinning(TH1D* &H1D_jetFakes, int iDataset, int iRadius, std::string options, bool controlMC);


#endif