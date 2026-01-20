#if !defined( __CINT__) || defined(__MAKECINT__)

#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include <iostream>
#include "TPaveText.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH3F.h"
#include "THn.h"
#include  "TDatime.h" 

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h" 
#include "RooUnfoldSvd.h" 
#include "TSVDUnfold.h"
#include "RooUnfoldBinByBin.h" 

#include "UnfoldHelper.h"

#endif

using namespace std;

TString strOutDir = "./outSys_varyResp";

//TString baseDir = "/hera/alice/obusch";
TString baseDir = "~/";


TString fnameResp(Form("%s/work/FFprep/miniJets/unfold/files/LHC10f6a/lego733/AnalysisResults.root",baseDir.Data()));
TString strInFileData(Form("%s/work/FFprep/miniJets/unfold/files/data/lego141/out10de.root",baseDir.Data())); 


TString strIDResp 
= "clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODMC_tpc1_tof3_cl0"; 

TString strIDFFRec
= "clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODb_tpc1_tof3_cl0";

TString strID_data 
= "clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_noGenJets_trackTypeUndef_jetTypeUndef_tpc1_tof3_cl0";


// TString strIDFFGen
// = "clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODMCb_tpc1_tof3_cl0"; 

// TString strIDFFRec
// = "clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00_clustersAODMC2_ANTIKT04_B0_Filter00272_Cut00150_Skip00_AODMC2b_AODb_tpc1_tof3_cl0";

UnfoldHelper *UFH; 

// ---------------------------------------------------

TList* getList(TString name = ""){

  TList* list = 0x0; 

    
  TString dirName  = "PWGJE_FragmentationFunction_" + name;
  TString listName = "fracfunc_" + name;    

  gDirectory->cd(dirName);

  list = (TList*) gDirectory->Get(listName);

  return list;
}

// -------------------------------------------------

void unfoldSpec(RooUnfoldResponse* respJetPt, 
		TH1D* h1JetPtRec, TH1D** hSpecUnfolded, TH2D** h2Cov, Int_t nIterSpec, Bool_t doBayes = kTRUE, Int_t errTreatment = 0){
  
  // ----------------
  // unfold spectrum 
  
  cout<<" unfold spectrum "<<endl;
    
  RooUnfold* unfoldSpec; 
  if(doBayes){
    unfoldSpec = new RooUnfoldBayes(respJetPt, h1JetPtRec, nIterSpec); 
    cout<<" Bayes unfolding, nIter "<<nIterSpec<<endl;
  }
  else if(nIterSpec){
    unfoldSpec = new RooUnfoldSvd(respJetPt, h1JetPtRec, nIterSpec); 
    unfoldSpec->IncludeSystematics (1); // TEST !!! error due to response matrix variations via pseudodata variations 
    cout<<" SVD unfolding, nIter "<<nIterSpec<<endl;
  }
  else{
    unfoldSpec = new RooUnfoldBinByBin(respJetPt, h1JetPtRec); 
    cout<<" use Bin-by-Bin correction "<<endl;
  }
    

  if(errTreatment == -1)     *hSpecUnfolded = (TH1D*) unfoldSpec->Hreco(RooUnfold::kCovariance);
  else if(errTreatment == 1) *hSpecUnfolded = (TH1D*) unfoldSpec->Hreco(RooUnfold::kErrors);
  else if(errTreatment == 2) *hSpecUnfolded = (TH1D*) unfoldSpec->Hreco(RooUnfold::kNoError);
  else if(errTreatment == 3) *hSpecUnfolded = (TH1D*) unfoldSpec->Hreco(RooUnfold::kCovToy);
  else                       *hSpecUnfolded = (TH1D*) unfoldSpec->Hreco();

  cout<<" hSpecUnfolded bin 5 error "<<(*hSpecUnfolded)->GetBinError(5)<<endl;


  // covariance matrix
  *h2Cov = UFH->getCovariance(unfoldSpec);

  delete unfoldSpec;
}

// ------------------------------------------------------------------------
void smearResponse(THnSparse* hnResponseOrg, THnSparse* hnResponseSmeared){

  const Int_t axisM0 = 0; 
  const Int_t axisT0 = 1;
 
  Int_t nBinsM0 = hnResponseOrg->GetAxis(axisM0)->GetNbins();
  Int_t nBinsT0 = hnResponseOrg->GetAxis(axisT0)->GetNbins();

  cout<<" smearResponse, nBinsM0 "<<nBinsM0<<endl;
  cout<<" smearResponse, nBinsT0 "<<nBinsT0<<endl;

  for(Int_t iT0=1; iT0<=nBinsT0; iT0++){

    if(!(iT0%100)) cout<<" iT0 "<<iT0<<endl;

    for(Int_t iM0=1; iM0<=nBinsM0; iM0++){
      
      Double_t coordM0 = hnResponseOrg->GetAxis(axisM0)->GetBinCenter(iM0);
      Double_t coordT0 = hnResponseOrg->GetAxis(axisT0)->GetBinCenter(iT0);
      
      Double_t binCoord[] = {coordM0,coordT0};
      
      Long64_t binIndex = hnResponseOrg->GetBin(binCoord);
      
      Double_t resp = hnResponseOrg->GetBinContent(binIndex); 
      Double_t err  = hnResponseOrg->GetBinError(binIndex); 
      Double_t relErr = 0; 
      if(resp) relErr = err/resp; 
      
      Double_t contSmeared = gRandom->Gaus(resp,err);
      Double_t errSmeared  = relErr*contSmeared; // was: elErr * resp ???

      if(contSmeared < 0){ // if smeared entry < 0, keep original
	contSmeared = resp;
	errSmeared  = err;
      }
      
      Long64_t binIndexSmeared = hnResponseSmeared->GetBin(binCoord);
	
      // if(resp) cout<<" smearResponse: coordM0 "<<coordM0<<" coordT0 "<<coordT0
      // 		   <<" resp "<<resp<<" err "<<err<<" relErr "<<err/resp<<" contSmeared "
      // 		   <<contSmeared<<" errSmeared "<<errSmeared<<endl;
      
      if(binIndexSmeared != binIndex){
	cout<<" smearResponse: mismatch binIndex "<<binIndex<<" binIndexSmeared "<<binIndexSmeared<<endl;
	exit(0);
      }
      
      hnResponseSmeared->SetBinContent(binIndex,contSmeared); 
      hnResponseSmeared->SetBinError(binIndex,errSmeared); 
    }
  } 
}

// ------------------------------------------------------------------------
void resetErrors(THnSparse* hnResponse){
  
  const Int_t axisM0 = 0; 
  const Int_t axisT0 = 1;
  
  Int_t nBinsM0 = hnResponse->GetAxis(axisM0)->GetNbins();
  Int_t nBinsT0 = hnResponse->GetAxis(axisT0)->GetNbins();

  cout<<" resetErrors, nBinsM0 "<<nBinsM0<<endl;
  cout<<" resetErrors, nBinsT0 "<<nBinsT0<<endl;

  for(Int_t iT0=1; iT0<=nBinsT0; iT0++){
    for(Int_t iM0=1; iM0<=nBinsM0; iM0++){
	
      Double_t coordM0 = hnResponse->GetAxis(axisM0)->GetBinCenter(iM0);
      Double_t coordT0 = hnResponse->GetAxis(axisT0)->GetBinCenter(iT0);
      
      Double_t binCoord[] = {coordM0,coordT0};
      
      Long64_t binIndex = hnResponse->GetBin(binCoord);
      
      hnResponse->SetBinError(binIndex,0); 
    }
  }
}

// -------------------------------------------------------------

void reweightHistG(TH1D* histG, TH1D* histOrg, TH1D* histSmear){

  Int_t nBins0 = histG->GetXaxis()->GetNbins();

  for(Int_t bin0=0; bin0<=nBins0; bin0++){ 

    Double_t contG     = histG->GetBinContent(bin0);
    Double_t errG      = histG->GetBinError(bin0);
      
    Double_t contOrg   = histOrg->GetBinContent(bin0);
      
    Double_t contSmear = histSmear->GetBinContent(bin0);
      
    if(contOrg){

      Double_t weight = contSmear / contOrg;

      Double_t contNew = contG * weight;
      Double_t errNew  = errG * weight;

      cout<<" bin0 "<<bin0<<" weight "<<weight<<" contG "<<contG<<" contNew "<<contNew<<" errG "<<errG<<" errNew "<<errNew<<endl;

      histG->SetBinContent(bin0,contNew);
      histG->SetBinError(bin0,errNew);
    }
  }
}

// ---------------------------------------------------

void GetResponse(RooUnfoldResponse** respJetPt, Int_t nBinsX, Double_t* binsX, Bool_t smearResp, Bool_t resetErr = kFALSE){
  
  TFile f(fnameResp); 

  TList* list = getList(strIDResp);

  THnSparse* hnRespJetPt_org      = (THnSparse*) list->FindObject("hnResponseJetPt"); 
  THnSparse* hnRespJetPtHistG_org = (THnSparse*) list->FindObject("hnRespJetPtHistG"); 
  THnSparse* hnRespJetPtHistM_org = (THnSparse*) list->FindObject("hnRespJetPtHistM"); // rec jets rec pt
  
  f.Close();
  
  // ---------------------  
  // rebin before smearing

  THnSparse* hnRespJetPtHistG = UFH->rebinTHn(hnRespJetPtHistG_org,nBinsX,binsX);
  THnSparse* hnRespJetPtHistM = UFH->rebinTHn(hnRespJetPtHistM_org,nBinsX,binsX);
  THnSparse* hnRespJetPt      = UFH->rebinTHn(hnRespJetPt_org,nBinsX,binsX,nBinsX,binsX);

  // -------------------------------
  // projections for RooUResponse
  
  TH1D* h1RespJetPtHistG = (TH1D*) hnRespJetPtHistG->Projection(0);  
  TH1D* h1RespJetPtHistM = (TH1D*) hnRespJetPtHistM->Projection(0); 

  TH1D* h1RespJetPt_unsmeared_projGen = hnRespJetPt->Projection(1);

  TH1D* h1RespJetPt_smeared_projRec  = 0;
  TH1D* h1RespJetPt_smeared_projGen  = 0;
  TH1D* h1RespJetPtHistG_rew = 0;
  
  if(smearResp){

    // smear response
    THnSparse* hnRespJetPt_tmp = (THnSparse*) hnRespJetPt->Clone("hnRespJetPt_tmp");
    smearResponse(hnRespJetPt_tmp, hnRespJetPt); // args: org, smeared
    delete hnRespJetPt_tmp;

    // project smeared response on rec axis to replace HistM (avoid spurious 'fakes' correction)
    h1RespJetPt_smeared_projRec = hnRespJetPt->Projection(0);
    h1RespJetPt_smeared_projRec->SetName("h1RespJetPt_smeared_projRec");

    cout<<" projection response entries "<<h1RespJetPt_smeared_projRec->GetEntries()<<endl;
    cout<<"  h1RespJetPtHistM entries "<<h1RespJetPtHistM->GetEntries()<<endl;

    // project smeared response on gen axis to replace HistG (avoid spurious 'efficiency' correction)
    h1RespJetPt_smeared_projGen = hnRespJetPt->Projection(1);
    h1RespJetPtHistG_rew = (TH1D*) h1RespJetPtHistG->Clone("h1RespJetPtHistG_rew");    
    reweightHistG(h1RespJetPtHistG_rew,h1RespJetPt_unsmeared_projGen,h1RespJetPt_smeared_projGen);
  }

  if(resetErr) resetErrors(hnRespJetPt);
  
  RooUnfoldResponse* respJetPtp;
  if(smearResp) respJetPtp = UFH->RooUResponseFromTHn(hnRespJetPt, h1RespJetPt_smeared_projRec, h1RespJetPtHistG_rew, 1, 1, 1, 1);
  else          respJetPtp = UFH->RooUResponseFromTHn(hnRespJetPt, h1RespJetPtHistM, h1RespJetPtHistG, 1, 1, 1, 1);

  *respJetPt = respJetPtp;

  delete hnRespJetPt_org;      
  delete hnRespJetPtHistG_org; 
  delete hnRespJetPtHistM_org;
  
  delete hnRespJetPt;      
  delete hnRespJetPtHistG; 
  delete hnRespJetPtHistM; 

  delete h1RespJetPtHistG;
  delete h1RespJetPtHistM;

  delete h1RespJetPt_smeared_projRec;
  delete h1RespJetPt_smeared_projGen;
  delete h1RespJetPtHistG_rew;
  
}

// ------------------------------------------------------------

void smearTH1(TH1D* hist, TH1D* hSmeared){

  // smear TH1 entries  according to stat error

  hSmeared->Reset();

  for(Int_t binx=1; binx<=hist->GetNbinsX(); binx++){
    
    Double_t cont = hist->GetBinContent(binx); 
    Double_t err  = hist->GetBinError(binx); 
      
    Double_t contSmeared = 0; 
    Double_t errCorr     = 0;

    if(cont){
      contSmeared = gRandom->Gaus(cont,err);
      errCorr = err*contSmeared/cont;
    }

    if(contSmeared < 0){ // if smeared entry < 0, keep original
      contSmeared = cont;
      errCorr     = err;
    }

    cout<<" hist name "<<hist->GetName()<<" binx "<<binx
	<<" cent x "<<hist->GetXaxis()->GetBinCenter(binx)
	<<" cont "<<cont<<" err "<<err<<" relErr "<<err/cont<<" contSmeared "<<contSmeared<<" errCorr "<<errCorr<<endl;

    hSmeared->SetBinContent(binx,contSmeared);
    hSmeared->SetBinError(binx,errCorr);
  }
}


// -------------------------------------------

void FillVariations(TH1D* hist, TProfile* histP){

  if(hist->GetNbinsX() != histP->GetNbinsX()) { cout<<" discrepancy hist/histP "<<endl; exit(0); }
  if(hist->GetXaxis()->GetXmin() != histP->GetXaxis()->GetXmin()) { cout<<" discrepancy hist/histP "<<endl; exit(0); }
  if(hist->GetXaxis()->GetXmax() != histP->GetXaxis()->GetXmax()) { cout<<" discrepancy hist/histP "<<endl; exit(0); }
  
  for(Int_t binx=1; binx<hist->GetNbinsX()+1; binx++){
  
    Double_t cent = hist->GetXaxis()->GetBinCenter(binx);
    Double_t cont = hist->GetBinContent(binx); 
    histP->Fill(cent,cont);
  }
} 


// ---------------------------------------------------

void unfoldSpec_statVariations(Int_t nRepeats = 100, Bool_t doBayes = kTRUE, Bool_t varySpec = kTRUE, Bool_t varyResp = kFALSE, Int_t errTreatment = 0){
  
  // 1D unfolding of spec 
  // test RooUnfold errors 
  // varyFFResp: smear response matrix
  // varyFFSpec: smear input spectrum histo
  // for pure spectrum variations 100 repeats take ~ 2-3 seconds, but response matrix variations very slow (leaky !) 
  
  UFH = new UnfoldHelper(fnameResp);
  
  Double_t binsX[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,38,44,50,58,66,76,86,100,120,150,200};
  Int_t nBinsX     = 29;

  Int_t nIterSpecBayes = 4;
  Int_t nIterSpecSVD   = 11;

  Int_t nIterSpec = doBayes ? nIterSpecBayes : nIterSpecSVD;
  
  // ------------------------------

  TFile f2(strInFileData);

  TList* list = getList(strID_data);

  TH1D* h1JetPtRec_org = (TH1D*) list->FindObject("hJetSpecIncRec");
  
  h1JetPtRec_org->SetDirectory(0);

  f2.Close();
 
  // ---------------------------
  // bias corr

  Bool_t biasCorr = kTRUE;

  if(biasCorr){

    TH1D* fh1JetPtRecMatchBC;
    TH1D* fh1JetPtRecIncBC;

    TH2D* fh2FFZRecMatchBC;
    TH2D* fh2FFZRecIncBC;

    TH1D* hRatioJetPtRecMatch  = UFH->ratioFFJetPtMatchInc(strIDResp, strIDFFRec, &fh1JetPtRecMatchBC, &fh1JetPtRecIncBC, kTRUE);

    h1JetPtRec_org->Multiply(hRatioJetPtRecMatch);
  }

  // ---------------------------
  // rebin
  
  TH1D* h1JetPtRec = (TH1D*) h1JetPtRec_org->Rebin(nBinsX,"h1JetPtRec_reb",binsX);

  // ---------------------------
  // response matrix 

  RooUnfoldResponse* respJetPt = 0;

  Bool_t smearResp = kFALSE;
  Bool_t resetErr  = kFALSE; // don't remember why errors were reset ? (errTreatment == 0) ? kFALSE : kTRUE;

  GetResponse(&respJetPt, nBinsX, binsX, smearResp, resetErr);

  TH1D* hSpecUnfolded  = 0;
  TH2D* h2Cov = 0;
  
  unfoldSpec(respJetPt, h1JetPtRec, &hSpecUnfolded, &h2Cov, nIterSpec, doBayes, errTreatment);

  //return; // TEST!!!
 
  // ---------------------------
  // smear + unfold

  TProfile* tpJetPtVar = new TProfile("tpJetPtVar","",nBinsX,binsX,"S");  
    
  for(Int_t rep=0; rep<nRepeats; rep++){

    cout<<"unfold FF rep "<<rep<<endl;

    TH1D* h1JetPtRec_rep =  new TH1D(*h1JetPtRec);
    h1JetPtRec_rep->SetName(Form("h1JetPtRec_rep%d",rep));
    
    RooUnfoldResponse* respJetPt_rep = 0;
    
    TH1D* hUnfolded_rep = 0; 
    TH2D* h2Cov         = 0; 
    
    if(varyResp){
      Bool_t smearResp = kTRUE;
      GetResponse(&respJetPt_rep, nBinsX, binsX, smearResp);
    }

    // delete respJetPt_rep; // TEST !!!
    // continue; // TEST !!!
    
    if(varySpec){
      h1JetPtRec_rep->Reset();
      smearTH1(h1JetPtRec,h1JetPtRec_rep);
    }

    if(varyResp) unfoldSpec(respJetPt_rep, h1JetPtRec_rep,  &hUnfolded_rep, &h2Cov, nIterSpec, doBayes, errTreatment);
    else         unfoldSpec(respJetPt,     h1JetPtRec_rep,  &hUnfolded_rep, &h2Cov, nIterSpec, doBayes, errTreatment);    
    
    FillVariations(hUnfolded_rep,tpJetPtVar); 

    cout<<" rep "<<rep<<" hUnfolded_rep bin Cont 5 "<<hUnfolded_rep->GetBinContent(5)<<endl;
    
    delete hUnfolded_rep;
    delete respJetPt_rep;
    delete h1JetPtRec_rep;		    
  }

  // --------------------------------------

  // unfolded spec errors as assigned by RooUnfold
  TH1F* hSpecUnfoldedErr = (TH1F*) hSpecUnfolded->Clone("hSpecUnfoldedErr");
  hSpecUnfoldedErr->Reset();
  
  for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
    if(hSpecUnfolded->GetBinContent(bin)){
      hSpecUnfoldedErr->SetBinContent(bin,hSpecUnfolded->GetBinError(bin)/hSpecUnfolded->GetBinContent(bin));
      hSpecUnfoldedErr->SetBinError(bin,0);
    }
  }

  // poissonian errors of unfolded spec bin content 
  TH1F* hSpecUnfoldedPoissErr = (TH1F*) hSpecUnfolded->Clone("hSpecUnfoldedPoissErr");
  hSpecUnfoldedPoissErr->Reset();
  
  for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
    if(hSpecUnfolded->GetBinContent(bin)){
      Double_t cont = hSpecUnfolded->GetBinContent(bin);
      Double_t poissRelErr = 0;
      if(cont > 0){
	poissRelErr = 1/TMath::Sqrt(cont);
      }
      hSpecUnfoldedPoissErr->SetBinContent(bin,poissRelErr);
      hSpecUnfoldedPoissErr->SetBinError(bin,0);
    }
  }

  // poissonian errors of raw spec bin content 
  TH1F* hSpecRawErr = (TH1F*) h1JetPtRec->Clone("hSpecRawErr");
  hSpecRawErr->Reset();
  
  for(Int_t bin=1; bin<=h1JetPtRec->GetNbinsX(); bin++){ 
    if(h1JetPtRec->GetBinContent(bin)){
      hSpecRawErr->SetBinContent(bin,h1JetPtRec->GetBinError(bin)/h1JetPtRec->GetBinContent(bin));
      hSpecRawErr->SetBinError(bin,0);
    }
  }

  // error from covariance
  TH1F* hCovDiagErr = (TH1F*) hSpecUnfolded->Clone("hCovDiagErr");
  hCovDiagErr->Reset();

  for(Int_t bin=1; bin<=hSpecUnfolded->GetNbinsX(); bin++){ 
    Double_t cont = hSpecUnfolded->GetBinContent(bin);
    if(cont){
      Double_t err = TMath::Sqrt(h2Cov->GetBinContent(bin,bin));
      hCovDiagErr->SetBinContent(bin,err/cont);
      hCovDiagErr->SetBinError(bin,0);
    }
  }

  
  // error from variations 
  TH1F* hVarErr = (TH1F*) hSpecUnfolded->Clone("hVarErr");
  hVarErr->Reset();
  
  for(Int_t bin=1; bin<=hVarErr->GetNbinsX(); bin++){ 

    //     //cout<<" jbin "<<jbin<<" bin "<<bin<<" hVar err "<<hVar->GetBinError(bin)<<endl;
    
    if(tpJetPtVar->GetBinContent(bin)){
      hVarErr->SetBinContent(bin,tpJetPtVar->GetBinError(bin)/tpJetPtVar->GetBinContent(bin));
      hVarErr->SetBinError(bin,0);
    }
  }

  // -----------------
  // scale by binwidth
  hSpecUnfolded->Scale(1,"width");
  TH1D* hJetPtVar = tpJetPtVar->ProjectionX("hJetPtVar"); // option 'width' doesn't work for TProfile, so first project to TH1
  hJetPtVar->Scale(1,"width");
  
  // ---------------------------------------
  // plot 

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","",610,470);
  c1->Divide(2,1);
   
  c1->cd(1);
  gPad->SetLogy();
  hSpecUnfolded->SetTitle("unfolded jet spectrum");
  hSpecUnfolded->SetLineColor(2);
  hSpecUnfolded->GetXaxis()->SetRangeUser(0,100);
  hSpecUnfolded->SetXTitle("p_{T}^{jet}"); 
  hSpecUnfolded->SetYTitle("1/nJets dN/dp_{T}"); 
  hSpecUnfolded->DrawCopy();
  
  // --

  c1->cd(2);
  gPad->SetLogy();
  hJetPtVar->SetTitle(Form("unfolded pseudodata, %d variations, error = spread",nRepeats));
  hJetPtVar->SetLineColor(2);
  hJetPtVar->GetXaxis()->SetRangeUser(0,100);
  hJetPtVar->SetXTitle("p_{T}^{jet}"); 
  hJetPtVar->SetYTitle("1/nJets dN/dp_{T}"); 
  hJetPtVar->DrawCopy();
  
  // --
  
  TCanvas *c2 = new TCanvas("c2","",460,560);
  c2->Divide(1,1);
   
  c2->cd(1);
  TString strTit;
  if(doBayes) strTit = "relative error, Bayes unfolding";
  else        strTit = "relative error, SVD unfolding";
  hVarErr->SetTitle(strTit);
  hVarErr->SetLineColor(2);
  hVarErr->GetXaxis()->SetRangeUser(0,100);
  hVarErr->GetYaxis()->SetRangeUser(0,0.2);
  hVarErr->SetMarkerStyle(20);
  hVarErr->SetMarkerColor(4);
  hVarErr->SetXTitle("p_{T}^{jet} (GeV/c)"); 
  hVarErr->SetYTitle("relative error");
  hVarErr->GetYaxis()->SetTitleOffset(1.2);
  hVarErr->Draw("PM");

  hSpecUnfoldedPoissErr->SetMarkerStyle(24);
  hSpecUnfoldedPoissErr->SetMarkerColor(4);
  hSpecUnfoldedPoissErr->Draw("PM same");

  
  hSpecUnfoldedErr->SetMarkerStyle(20);
  hSpecUnfoldedErr->SetMarkerColor(2);
  hSpecUnfoldedErr->Draw("PM same");
   
  TLegend* leg1 = new TLegend(0.15,0.70,0.61,0.87);
  leg1->SetTextSize(0.02);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(1);
  if(errTreatment == -1)     leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kCovariance, from cov. matrix","P");
  else if(errTreatment == 1) leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kErrors, from diag. elements of cov. matrix","P");
  else if(errTreatment == 2) leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kNoError","P");
  else if(errTreatment == 3) leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment kCovToy","P");
  else 	 leg1->AddEntry(hSpecUnfoldedErr,"RooUnfold: errTreatment default ","P");
  leg1->AddEntry(hVarErr,"unfolded pseudodata: spread","P");
  leg1->AddEntry(hSpecUnfoldedPoissErr,"Poisson error unfolded spectrum","P");
  leg1->Draw("");


  TCanvas *c3 = new TCanvas("c3","",460,560);
  c3->Divide(1,1);
   
  c3->cd(1);
  TString strTit2;
  if(doBayes) strTit2 = "relative error cov matrix, Bayes unfolding";
  else        strTit2 = "relative error cov matrix, SVD unfolding";
  hCovDiagErr->SetTitle(strTit2);
  hCovDiagErr->SetLineColor(2);
  hCovDiagErr->GetXaxis()->SetRangeUser(0,100);
  hCovDiagErr->GetYaxis()->SetRangeUser(0,0.2);
  hCovDiagErr->SetMarkerStyle(20);
  hCovDiagErr->SetMarkerColor(4);
  hCovDiagErr->SetXTitle("p_{T}^{jet} (GeV/c)"); 
  hCovDiagErr->SetYTitle("relative error");
  hCovDiagErr->GetYaxis()->SetTitleOffset(1.2);
  hCovDiagErr->Draw("PM");

   
  // --

  TString strApp;
  if(doBayes) strApp.Form("Bayes_statVariations_nRep%d",nRepeats);
  else        strApp.Form("SVD_statVariations_nRep%d",nRepeats);
  if(errTreatment == -1) strApp += "_errkCov";
  if(errTreatment ==  1) strApp += "_errkErr";
  if(errTreatment ==  2) strApp += "_errkNoError";
  if(errTreatment ==  3) strApp += "_errkCovToy";
  if(varySpec) strApp += "_varySpec";
  if(varyResp) strApp += "_varyResp";
  
  // c1->SaveAs(Form("spectra_%s.pdf",strApp.Data()));
  c2->SaveAs(Form("errors_%s.pdf",strApp.Data()));

  // ---- 

  delete hSpecUnfolded;
}

