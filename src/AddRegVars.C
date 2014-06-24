#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include "TLorentzVector.h"
using namespace TMVA;

void AddRegVars(char*, char*, char*, int);

int main(int argc, char *argv[]){
  if(argc==1){
    printf("AddRegVars [treeName] [splitFile]\n");
    return 1;
  }

  char treeName[80];
  int splitFile=-1;//-1 for entire file, 0 for even events only, 1 for odd events only
  if(argc>1) sprintf(treeName,"%s",argv[1]);
  if(argc>2) splitFile=atoi(argv[2]);

  char inputFile[160];
  sprintf(inputFile,"files/tree_Signal.root");

  char outputFile[160];
  sprintf(outputFile,"files/tree_Signal_regression.root");

  AddRegVars(inputFile, outputFile, treeName, splitFile);

  return 0;
}

void AddRegVars(char* inputFileName, char* outputFileName, char* treeName, int splitFile){

  TFile *inFile = new TFile(inputFileName,"READ");
  TTree *inTree = (TTree*)inFile->Get(treeName);
  TFile *outFile = new TFile(outputFileName,"UPDATE");//for the first iteration, change UPDATE to RECREATE
  outFile->Delete(Form("%s;*",treeName));
  TTree *outTree = new TTree(treeName,"Regression data with results");

  const int maxP = 15;
  int event;
  float rho25, MET, METphi;
  float jet_pt[maxP], jet_eta[maxP], jet_phi[maxP], jet_e[maxP], jet_secVtxPt[maxP], jet_secVtxM[maxP], jet_secVtx3dL[maxP], jet_secVtx3deL[maxP], jet_softLeptPt[maxP], jet_softLeptPtRel[maxP], jet_softLeptDR[maxP], jet_emfrac[maxP], jet_hadfrac[maxP], jet_chadfrac[maxP], jet_nhadfrac[maxP], jet_phofrac[maxP], jet_mufrac[maxP], jet_elefrac[maxP], jet_JECUnc[maxP], jet_leadTrackPt[maxP], jet_genPt[maxP], jet_csvBtag[maxP];
  int jet_bgenMatched[maxP], jet_flavour[maxP], jet_softLeptIdLooseMu[maxP], jet_softLeptIdEle95[maxP], jet_nNeutrals[maxP], jet_nCharged[maxP];

  inTree->SetBranchAddress("event", &event);
  inTree->SetBranchAddress("rho", &rho25);
  inTree->SetBranchAddress("met_corr_pfmet", &MET);
  inTree->SetBranchAddress("met_corr_phi_pfmet", &METphi);

  for(int i=0; i<maxP; ++i){
    inTree->SetBranchAddress(TString::Format("j%i_pt",i+1),&(jet_pt[i]));
    inTree->SetBranchAddress(TString::Format("j%i_phi",i+1),&jet_phi[i]);
    inTree->SetBranchAddress(TString::Format("j%i_eta",i+1),&jet_eta[i]);
    inTree->SetBranchAddress(TString::Format("j%i_e",i+1),&jet_e[i]);
    inTree->SetBranchAddress(TString::Format("j%i_bgenMatched",i+1),&jet_bgenMatched[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtxPt",i+1),&jet_secVtxPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtxM",i+1),&jet_secVtxM[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtx3dL",i+1),&jet_secVtx3dL[i]);
    inTree->SetBranchAddress(TString::Format("j%i_secVtx3deL",i+1),&jet_secVtx3deL[i]);
    inTree->SetBranchAddress(TString::Format("j%i_emfrac",i+1),&jet_emfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_hadfrac",i+1),&jet_hadfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_nNeutrals",i+1),&jet_nNeutrals[i]);
    inTree->SetBranchAddress(TString::Format("j%i_nCharged",i+1),&jet_nCharged[i]);
    inTree->SetBranchAddress(TString::Format("j%i_chadfrac",i+1),&jet_chadfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_nhadfrac",i+1),&jet_nhadfrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_phofrac",i+1),&jet_phofrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_mufrac",i+1),&jet_mufrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_elefrac",i+1),&jet_elefrac[i]);
    inTree->SetBranchAddress(TString::Format("j%i_JECUnc",i+1),&jet_JECUnc[i]);
    inTree->SetBranchAddress(TString::Format("j%i_leadTrackPt",i+1),&jet_leadTrackPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptPt",i+1),&jet_softLeptPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptPtRel",i+1),&jet_softLeptPtRel[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptDR",i+1),&jet_softLeptDR[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptIdLooseMu",i+1),&jet_softLeptIdLooseMu[i]);
    inTree->SetBranchAddress(TString::Format("j%i_softLeptIdEle95",i+1),&jet_softLeptIdEle95[i]);
    inTree->SetBranchAddress(TString::Format("j%i_genPt",i+1),&jet_genPt[i]);
    inTree->SetBranchAddress(TString::Format("j%i_flavour",i+1),&jet_flavour[i]);
    inTree->SetBranchAddress(TString::Format("j%i_csvBtag",i+1),&jet_csvBtag[i]);
  }


  float bjet_pt, bjet_phi, bjet_eta, bjet_e, bjet_secVtxPt, bjet_secVtxM, bjet_secVtx3dL, bjet_secVtx3deL, bjet_emfrac, bjet_hadfrac, bjet_chadfrac, bjet_nhadfrac, bjet_phofrac, bjet_mufrac, bjet_elefrac, bjet_JECUnc, bjet_leadTrackPt, bjet_softLeptPt, bjet_softLeptPtRel, bjet_softLeptDR, bjet_genPt;
  int bjet_bgenMatched, bjet_flavour, bjet_softLeptIdLooseMu, bjet_softLeptIdEle95, bjet_nNeutrals, bjet_nCharged;

  int j1_hasSoftLept, j2_hasSoftLept;
  float j1_pt, j2_pt, j1_genPt, j2_genPt, j1_ptcorr, j2_ptcorr, j1_DRMETJet, j2_DRMETJet, Mjj, Mjjcorr;

  outTree->Branch("event", &event, "event/I");
  outTree->Branch("j1_pt", &j1_pt, "j1_pt/F");
  outTree->Branch("j2_pt", &j2_pt, "j2_pt/F");
  outTree->Branch("j1_genPt", &j1_genPt, "j1_genPt/F");
  outTree->Branch("j2_genPt", &j2_genPt, "j2_genPt/F");
  outTree->Branch("j1_ptcorr", &j1_ptcorr, "j1_ptcorr/F");
  outTree->Branch("j2_ptcorr", &j2_ptcorr, "j2_ptcorr/F");
  outTree->Branch("j1_hasSoftLept", &j1_hasSoftLept, "j1_hasSoftLept/I");
  outTree->Branch("j2_hasSoftLept", &j2_hasSoftLept, "j2_hasSoftLept/I");
  outTree->Branch("j1_DRMETJet", &j1_DRMETJet, "j1_DRMETJet/F");
  outTree->Branch("j2_DRMETJet", &j2_DRMETJet, "j2_DRMETJet/F");
  outTree->Branch("Mjj", &Mjj, "Mjj/F");
  outTree->Branch("Mjjcorr", &Mjjcorr, "Mjjcorr/F");

  float var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,var18,var19,var20,var21;
  Reader *readerRegres = new Reader( "!Color:!Silent" );

  readerRegres->AddVariable( "bjet_pt", &var1);
  readerRegres->AddVariable( "bjet_eta", &var2);
  readerRegres->AddVariable( "bjet_m", &var3);
  readerRegres->AddVariable( "bjet_phofrac", &var4);
  readerRegres->AddVariable( "bjet_nhadfrac", &var5);
  readerRegres->AddVariable( "(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptPt) : (-99)", &var6);
  readerRegres->AddVariable( "(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptPtRel) : (-99)", &var7);
  readerRegres->AddVariable( "bjet_secVtxM", &var8);
  readerRegres->AddVariable( "bjet_secVtx3deL", &var9);
  readerRegres->AddVariable( "MET", &var10);
  readerRegres->AddVariable( "(abs(bjet_phi-METphi)>3.14159265 ) ? (2*3.14159265-abs(bjet_phi-METphi)) : (abs(bjet_phi-METphi))", &var11);
  readerRegres->AddVariable( "bjet_leadTrackPt", &var12);
  readerRegres->AddVariable( "bjet_nCharged+bjet_nNeutrals", &var13);
  readerRegres->AddVariable( "rho", &var14);
  readerRegres->AddSpectator( "bjet_mufrac", &var15);
  readerRegres->AddSpectator( "bjet_elefrac", &var16);
  readerRegres->AddSpectator( "bjet_chadfrac", &var17);
  readerRegres->AddSpectator( "(bjet_softLeptIdLooseMu==1 || bjet_softLeptIdEle95==1) ? (bjet_softLeptDR) : (-99)", &var18);
  readerRegres->AddSpectator( "bjet_secVtxPt", &var19);
  readerRegres->AddSpectator( "bjet_secVtx3dL", &var20);
  readerRegres->AddSpectator( "bjet_JECUnc", &var21);

  if(strcmp(treeName,"ggHH_8TeV")==0){
    printf("SM training applied\n");
    readerRegres->BookMVA("BDTG",TString::Format("weights/TMVARegression_SM_BDTG.weights.xml"));
  }
  else{
    printf("resonant training applied\n");
    readerRegres->BookMVA("BDTG",TString::Format("weights/TMVARegression_resonant_BDTG.weights.xml"));
  }

  unsigned nEntries = inTree->GetEntries();
  printf("\nProcessing %i entries of %i (splitFile=%i)\n\n", inTree->GetEntries(Form("event%2==%i || %i==-1",splitFile,splitFile)), inTree->GetEntries(), splitFile);

  for(unsigned i=0; i<nEntries; ++i){
    inTree->GetEntry(i);
    if(splitFile!=-1 && splitFile!=event%2) continue;

    std::vector<int> bjets;
    for(int j=0; j<maxP && jet_pt[j]>0; ++j){
      if(jet_csvBtag[j]>0.679) bjets.push_back(j);
    }
    if(bjets.size()>2){
      int tmp1=-1, tmp2=-1; float ptH=-1;
      for(int j=0; j<bjets.size(); ++j){
	for(int k=j+1; k<bjets.size(); ++k){
	  TLorentzVector tmpj1, tmpj2;
	  tmpj1.SetPtEtaPhiE(jet_pt[bjets[j]],jet_eta[bjets[j]],jet_phi[bjets[j]],jet_e[bjets[j]]);
	  tmpj2.SetPtEtaPhiE(jet_pt[bjets[k]],jet_eta[bjets[k]],jet_phi[bjets[k]],jet_e[bjets[k]]);
	  if(ptH < (tmpj1+tmpj2).Pt() ){
	    ptH=(tmpj1+tmpj2).Pt();
	    tmp1=j; tmp2=k;
	  }
	}
      }
      bjets.clear();
      bjets.push_back(tmp1);
      bjets.push_back(tmp2);
    }
    else if(bjets.size()<2) continue;
    //else if(bjets.size()==1){
    //  int tmp2=-1; float ptH=-1;
    //  TLorentzVector tmpj1;
    //  tmpj1.SetPtEtaPhiE(jet_pt[bjets[0]],jet_eta[bjets[0]],jet_phi[bjets[0]],jet_e[bjets[0]]);
    //  for(int j=0; j<maxP && jet_pt[j]>0; ++j){
    //	if(j==bjets[0]) continue;
    //	TLorentzVector tmpj2;
    //	tmpj2.SetPtEtaPhiE(jet_pt[j],jet_eta[j],jet_phi[j],jet_e[j]);
    //	if(ptH < (tmpj1+tmpj2).Pt() ){
    //	  ptH=(tmpj1+tmpj2).Pt();
    //	  tmp2=j;
    //	}
    //  }
    //  bjets.push_back(tmp2);
    //}

    var1=jet_pt[bjets[0]];
    var2=jet_eta[bjets[0]];
    var3=sqrt(pow(jet_e[bjets[0]],2)-pow(jet_pt[bjets[0]]*(1+sinh(jet_eta[bjets[0]])),2));
    var4=jet_phofrac[bjets[0]];
    var5=jet_nhadfrac[bjets[0]];
    var6=(jet_softLeptIdLooseMu[bjets[0]]==1 || jet_softLeptIdEle95[bjets[0]]==1) ? (jet_softLeptPt[bjets[0]]) : (-99);
    var7=(jet_softLeptIdLooseMu[bjets[0]]==1 || jet_softLeptIdEle95[bjets[0]]==1) ? (jet_softLeptPtRel[bjets[0]]) : (-99);
    var8=jet_secVtxM[bjets[0]];
    var9=jet_secVtx3deL[bjets[0]];
    var10=MET;
    var11=(fabs(jet_phi[bjets[0]]-METphi)>3.14159265 ) ? (2*3.14159265-fabs(jet_phi[bjets[0]]-METphi)) : (fabs(jet_phi[bjets[0]]-METphi));
    var12=jet_leadTrackPt[bjets[0]];
    var13=jet_nCharged[bjets[0]]+jet_nNeutrals[bjets[0]];
    var14=rho25;
    var15=jet_mufrac[bjets[0]];
    var16=jet_elefrac[bjets[0]];
    var17=jet_chadfrac[bjets[0]];
    var18=(jet_softLeptIdLooseMu[bjets[0]]==1 || jet_softLeptIdEle95[bjets[0]]==1) ? (jet_softLeptDR[bjets[0]]) : (-99);
    var19=jet_secVtxPt[bjets[0]];
    var20=jet_secVtx3dL[bjets[0]];
    var21=jet_JECUnc[bjets[0]];
    j1_pt=jet_pt[bjets[0]];
    j1_genPt=jet_genPt[bjets[0]];
    j1_ptcorr=readerRegres->EvaluateRegression("BDTG")[0];
    j1_DRMETJet=var11;
    j1_hasSoftLept=(jet_softLeptIdLooseMu[bjets[0]]==1 || jet_softLeptIdEle95[bjets[0]]==1);

    var1=jet_pt[bjets[1]];
    var2=jet_eta[bjets[1]];
    var3=sqrt(pow(jet_e[bjets[1]],2)-pow(jet_pt[bjets[1]]*(1+sinh(jet_eta[bjets[1]])),2));
    var4=jet_phofrac[bjets[1]];
    var5=jet_nhadfrac[bjets[1]];
    var6=(jet_softLeptIdLooseMu[bjets[1]]==1 || jet_softLeptIdEle95[bjets[1]]==1) ? (jet_softLeptPt[bjets[1]]) : (-99);
    var7=(jet_softLeptIdLooseMu[bjets[1]]==1 || jet_softLeptIdEle95[bjets[1]]==1) ? (jet_softLeptPtRel[bjets[1]]) : (-99);
    var8=jet_secVtxM[bjets[1]];
    var9=jet_secVtx3deL[bjets[1]];
    var10=MET;
    var11=(fabs(jet_phi[bjets[1]]-METphi)>3.14159265 ) ? (2*3.14159265-fabs(jet_phi[bjets[1]]-METphi)) : (fabs(jet_phi[bjets[1]]-METphi));
    var12=jet_leadTrackPt[bjets[1]];
    var13=jet_nCharged[bjets[1]]+jet_nNeutrals[bjets[1]];
    var14=rho25;
    var15=jet_mufrac[bjets[1]];
    var16=jet_elefrac[bjets[1]];
    var17=jet_chadfrac[bjets[1]];
    var18=(jet_softLeptIdLooseMu[bjets[1]]==1 || jet_softLeptIdEle95[bjets[1]]==1) ? (jet_softLeptDR[bjets[1]]) : (-99);
    var19=jet_secVtxPt[bjets[1]];
    var20=jet_secVtx3dL[bjets[1]];
    var21=jet_JECUnc[bjets[1]];
    j2_pt=jet_pt[bjets[1]];
    j2_genPt=jet_genPt[bjets[1]];
    j2_ptcorr=readerRegres->EvaluateRegression("BDTG")[0];
    j2_DRMETJet=var11;
    j2_hasSoftLept=(jet_softLeptIdLooseMu[bjets[1]]==1 || jet_softLeptIdEle95[bjets[1]]==1);

    TLorentzVector bjet41,bjet42;
    bjet41.SetPtEtaPhiE(j1_pt,jet_eta[bjets[0]],jet_phi[bjets[0]],jet_e[bjets[0]]);
    bjet42.SetPtEtaPhiE(j2_pt,jet_eta[bjets[1]],jet_phi[bjets[1]],jet_e[bjets[1]]);
    Mjj = (bjet41+bjet42).M();
    bjet41.SetPtEtaPhiE(j1_ptcorr,jet_eta[bjets[0]],jet_phi[bjets[0]],jet_e[bjets[0]]*j1_ptcorr/j1_pt);
    bjet42.SetPtEtaPhiE(j2_ptcorr,jet_eta[bjets[1]],jet_phi[bjets[1]],jet_e[bjets[1]]*j2_ptcorr/j2_pt);
    Mjjcorr = (bjet41+bjet42).M();

    outTree->Fill();
  }
  outTree->Write();
  outFile->Close();
  inFile->Close();
}
