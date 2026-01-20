#ifndef JETSPECTRUM_SPECTRAGETTERS_H
#define JETSPECTRUM_SPECTRAGETTERS_H

void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_Pt_spectrum_bkgCorrected_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_Pt_spectrum_bkgCorrected_recBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options, bool controlMC);

void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_Pt_spectrum_bkgCorrected_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_Pt_spectrum_bkgCorrected_genBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options, bool controlMC);

void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_Pt_spectrum_bkgCorrected_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt, int iDataset, int iRadius, std::string options, bool controlMC);
void Get_Pt_spectrum_bkgCorrected_fineBinning(TH1D* &H1D_jetPt_raw, int iDataset, int iRadius, std::string options, bool controlMC);

void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcp_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcp_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcp_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcp_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcp_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcp_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcdMatched_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcdMatched_genBinning(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcpMatched_genBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcpMatched_genBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcd_fineBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcd_fineBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEndAndEvtNorm(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcd_recBinning_preWidthScalingAtEnd(TH1D* &H1D_jetPt_mcd, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcd_recBinning(TH1D* &H1D_jetPt_mcp, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);

void Get_Pt_spectrum_mcpMatched_fineBinning(TH1D* &H1D_jetPt_mcpMatched, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcdMatched_recBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);
void Get_Pt_spectrum_mcdMatched_fineBinning(TH1D* &H1D_jetPt_mcdMatched, int iDataset, int iRadius, std::string options, bool controlMC, bool controlMC_useResponseSplit);


#endif