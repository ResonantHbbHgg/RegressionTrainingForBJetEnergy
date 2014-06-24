#include "TROOT.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "THStack.h"
#include "TTree.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TStyle.h"
#include "TCut.h"
#include "TError.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <utility>
#include <string>
#include <sstream>
#include "/afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/src/tdrstyle.C"

// RooFit headers
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"

using namespace std;
using namespace RooFit;

int doAltSample=0;
bool doHistApprox=false;
void RegVal(int, int, char*, char*);
void plotParameters(TCanvas *c, int Rmass, RooArgList *param, RooArgList *paramCorr);
void plotParameters(RooArgList *param, RooArgList *paramCorr, TCanvas *c, int Rmass, bool isVoigtian, bool isRes, double KSprob);
void plotParameters(TCanvas *c, int Rmass, double *arrayUncorr, int nEventsUncorr, double *arrayCorr, int nEventsCorr, RooRealVar *sig1, RooRealVar *sig2 );
void sortArray(double*, const int);
double sigma_eff(double*, int, double&);
double sigma_eff(TH1*);

int main(int argc, char *argv[])
{
  if(argc==1){
    printf("RegVal2 [optimIndex] [Rmass] [fileNameOdd] [fileNameEven=\"\"]\n");
    return 1;
  }

  int optimIndex=445;
  int optimIndex2=53;
  int mass=300;
  bool splitFile=1;
  char *newFileBase="";
  if(argc>1) optimIndex=atoi(argv[1]);
  if(argc>2) mass=atoi(argv[2]);
  if(argc>3) splitFile=atoi(argv[3]);
  if(argc>4) newFileBase=argv[4];


  char inputFile[80], inputFileEven[80], inputFileOdd[80];
  sprintf(inputFile,"RegOutputFiles/tree_M%i_%i.root",mass,optimIndex);
  sprintf(inputFileEven,"RegOutputFiles/tree_M%i_%i_even.root",mass,optimIndex);
  sprintf(inputFileOdd,"RegOutputFiles/tree_M%i_%i_odd.root",mass,optimIndex);

  char newInputFile[80];
  sprintf(newInputFile,"RegOutputFiles/%s_regOutput_%i.root",newFileBase,optimIndex);//4MB
  if(strstr(newFileBase,"ZH_ZToLL_HToBB_M-125")!=NULL) doAltSample=1;////181MB
  else if(strstr(newFileBase,"WH_ZH_HToGG_M-125")!=NULL) doAltSample=2;//1.1MB
  else if(strstr(newFileBase,"ZH_ZToBB_HToInv_M-125")!=NULL) doAltSample=3;//90MB
  else if(strstr(newFileBase,"ZZ_ZToLL_ZToBB")!=NULL) doAltSample=4;//481MB

  if(newFileBase[0]!='\0')
    RegVal(optimIndex, mass, newInputFile, "");
  else if(splitFile)
    RegVal(optimIndex, mass, inputFileOdd, inputFileEven);
  else
    RegVal(optimIndex, mass, inputFile, "");

  return 0;

}

void RegVal(int optimIndex, int Rmass, char *fileNameOdd, char *fileNameEven)
{
  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".L /afs/cern.ch/work/h/hebda/HggHbb/CMSSW_5_2_5/src/ggAnalysis/regression/src/tdrstyle.C");
  setTDRStyle();
  int mass=125;
  bool doEvenFile=1;
  if(fileNameEven[0]=='\0') doEvenFile=0;

  //TCut thecut = "hJet_pt[0]>20. && hJet_pt[1]>20. && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_csv[0]>0. && hJet_csv[1]>0. && hJet_ptLeadTrack[0]<1500. && hJet_ptLeadTrack[1]<1500. && hJet_genJetPt[0]>0. && hJet_genJetPt[1]>0. && hJet_puJetIdL[0]>0.0 && hJet_puJetIdL[1]>0.0";
  TCut trainingCut = "hJet_pt[0]>20. && hJet_pt[1]>20. && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_genJetPt[0]>0. && hJet_genJetPt[1]>0. && hJet_dRJetGenJet[0]<0.4 && hJet_dRJetGenJet[1]<0.4";
  TCut testingCut = "hJet_pt[0]>20. && hJet_pt[1]>20. && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_genJetPt[0]>0. && hJet_genJetPt[1]>0. && hJet_dRJetGenJet[0]<0.4 && hJet_dRJetGenJet[1]<0.4";

 TFile *fOdd, *fEven;
  TTree *tOdd, *tEven, *tOddCut, *tEvenCut;
  fOdd = new TFile(fileNameOdd);
  tOdd = (TTree*)fOdd->Get("Events");
  if(doEvenFile){
    fEven = new TFile(fileNameEven);
    tEven = (TTree*)fEven->Get("Events");
  } 

  TFile *temp = new TFile(Form("temp_%i.root",optimIndex),"RECREATE");//avoids errors
  tOddCut = tOdd->CopyTree(testingCut);
  if(doEvenFile) tEvenCut = tEven->CopyTree(trainingCut);
  else tEvenCut = new TTree();

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TGaxis::SetMaxDigits(3);

  // Observables
  double lowVal=40, highVal=200;
  //double lowVal=60., highVal=180.;
  if(doAltSample>1){ lowVal=20.; highVal=170; }
  double lowValVar=0., highValVar=9999.;
  int nbinsVar=(highValVar-lowValVar)/3, nbins=(highVal-lowVal)/3;
  if(doAltSample==1 || doAltSample==3 || doAltSample==4) {
    nbinsVar=1.0*(highValVar-lowValVar);
    nbins=1.0*(highVal-lowVal);
    doHistApprox=true;
  }

  //restrict the range
  lowValVar=lowVal; highValVar=highVal;
  nbinsVar=nbins;
  double fitRangeLow=85, fitRangeHigh=165;
  if(Rmass==300) {fitRangeLow=85, fitRangeHigh=165; }
  else if(Rmass==500) {fitRangeLow=100, fitRangeHigh=145; }
  else if(Rmass==700 ) {fitRangeLow=100, fitRangeHigh=150; }
  else if(Rmass==1000) {fitRangeLow=110, fitRangeHigh=160; }
  fitRangeLow=60; fitRangeHigh=180;
  TCut rangeCut = TString::Format("Hmass>%.0f && Hmass<%.0f",fitRangeLow,fitRangeHigh).Data();
  TCut rangeCutCorr = TString::Format("HmassCorr>%.0f && HmassCorr<%.0f",fitRangeLow,fitRangeHigh).Data();

  RooRealVar Hmass("Hmass", "M_{jj}", lowValVar, highValVar, "GeV");
  Hmass.setBins(nbinsVar);
  RooRealVar HmassCorr("HmassCorr", "M_{jj}", lowValVar, highValVar, "GeV");
  HmassCorr.setBins(nbinsVar);
  RooDataSet datasetOdd("radion", "radion", tOddCut, RooArgList(Hmass,HmassCorr));
  RooDataSet datasetEven("radion", "radion", tEvenCut, RooArgList(Hmass,HmassCorr));

  double mean_jj = datasetOdd.mean(Hmass);
  double mean_jjCorr = datasetOdd.mean(HmassCorr);
  double rms_jj = datasetOdd.rmsVar(Hmass)->getVal();
  double rms_jjCorr = datasetOdd.rmsVar(HmassCorr)->getVal();
  double mean_jjErr = datasetOdd.meanVar(Hmass)->getError();//rms_jj/sqrt(datasetOdd.numEntries());
  double mean_jjCorrErr = datasetOdd.meanVar(HmassCorr)->getError();//rms_jjCorr/sqrt(datasetOdd.numEntries());
  double rms_jjErr = datasetOdd.rmsVar(Hmass)->getError();
  double rms_jjCorrErr = datasetOdd.rmsVar(HmassCorr)->getError();
  double newMean, newMeanErr, resChangeV, resChangeErrV, resChangeG, resChangeErrG, resChangeS, resChangeErrS, KSprob;

  bool GAUSS = true;
  bool VOIGT = true;
  bool doGen = false;

  TH1D* HmassUncorrHist = new TH1D("h1","h1",nbinsVar,lowValVar,highValVar);
  TH1D* HmassCorrHist = new TH1D("h2","h2",nbinsVar,lowValVar,highValVar);
  if(doHistApprox){
    tOddCut->Draw("Hmass>>h1",testingCut);
    tOddCut->Draw("HmassCorr>>h2",testingCut);
  }

  if(doEvenFile){
    int na = tOdd->Draw("HmassCorr",trainingCut+rangeCutCorr);
    double* a = tOdd->GetV1();
    sortArray(a,na);

    int nb = tEven->Draw("HmassCorr",trainingCut+rangeCutCorr);
    double* b = tEven->GetV1();
    sortArray(b,nb);

    KSprob = TMath::KolmogorovTest(na,a,nb,b,"");
  }
  else
    KSprob=-99;

  tOddCut->Draw("Hmass",testingCut+rangeCut);
  double* HmassUncorrArray = tOddCut->GetV1();
  int sizeOfHmassUncorrArray = tOddCut->GetEntries(testingCut+rangeCut);
  tOddCut->Draw("HmassCorr",testingCut+rangeCutCorr);
  double* HmassCorrArray = tOddCut->GetV1();
  int sizeOfHmassCorrArray = tOddCut->GetEntries(testingCut+rangeCutCorr);
  if(!doHistApprox){
    sortArray(HmassUncorrArray,sizeOfHmassUncorrArray);
    sortArray(HmassCorrArray,sizeOfHmassCorrArray);
  }

  if(VOIGT)
    {
      RooRealVar mu_voigt("mu_voigt", "mean (uncorr.)", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
      RooRealVar width_voigt("width_voigt", "width (uncorr.)", rms_jj, 0.1*rms_jj, 5.*rms_jj, "GeV");
      RooRealVar sigma_voigt("sigma_voigt", "sigma (uncorr.)", 5., 0.0001, 20., "GeV");
      RooVoigtian voigt("voigt", "voigt", Hmass, mu_voigt, width_voigt, sigma_voigt);

      RooRealVar mu_voigtCorr("mu_voigtCorr", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
      RooRealVar width_voigtCorr("width_voigtCorr", "width (corr.)", rms_jjCorr, 0.1*rms_jjCorr, 5.*rms_jjCorr, "GeV");
      RooRealVar sigma_voigtCorr("sigma_voigtCorr", "sigma (corr.)", 5., 0.0001, 20., "GeV");
      RooVoigtian voigtCorr("voigtCorr", "voigtCorr", HmassCorr, mu_voigtCorr, width_voigtCorr, sigma_voigtCorr);

      RooPlot * jj_frame = Hmass.frame(lowVal,highVal);
      RooPlot * jjCorr_frame = HmassCorr.frame(lowVal,highVal);
      datasetOdd.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      datasetOdd.plotOn(jjCorr_frame, LineColor(kBlue+2), MarkerColor(kBlue), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      RooFitResult *f = voigt.fitTo(datasetOdd, Save(), Range(fitRangeLow,fitRangeHigh));
      voigt.plotOn(jj_frame, LineColor(kRed+2), LineWidth(3));
      RooFitResult * fCorr = voigtCorr.fitTo(datasetOdd, Save(), Range(fitRangeLow,fitRangeHigh));
      voigtCorr.plotOn(jjCorr_frame, LineColor(kBlue+2), LineWidth(3));
      RooArgList p = f->floatParsFinal();
      RooArgList pCorr = fCorr->floatParsFinal();
      jjCorr_frame->GetYaxis()->SetTitleOffset(1.25);
      jjCorr_frame->SetTitle("");
      jjCorr_frame->Draw();
      jj_frame->Draw("same");
      double fg = 2. * sigma_voigt.getVal() * sqrt(2. * log(2));
      double fl = 2. * width_voigt.getVal();
      double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
      double fgCorr = 2. * sigma_voigtCorr.getVal() * sqrt(2. * log(2));
      double flCorr = 2. * width_voigtCorr.getVal();
      double fvCorr = 0.5346 * flCorr + sqrt(0.2166 * pow(flCorr, 2.) + pow(fgCorr, 2.));
      double res= fv / mu_voigt.getVal() * 100. / (2. * sqrt(2. * log(2.)));
      double resErr= res * sqrt(pow(mu_voigt.getError()/mu_voigt.getVal(),2) + pow(sigma_voigt.getError()/sigma_voigt.getVal(),2) + pow(width_voigt.getError()/width_voigt.getVal(),2));
      double resCorr= fvCorr / mu_voigtCorr.getVal() * 100. / (2. * sqrt(2. * log(2.)));
      double resCorrErr= res * sqrt(pow(mu_voigtCorr.getError()/mu_voigtCorr.getVal(),2) + pow(sigma_voigtCorr.getError()/sigma_voigtCorr.getVal(),2) + pow(width_voigtCorr.getError()/width_voigtCorr.getVal(),2));
      newMean = mu_voigtCorr.getVal();
      newMeanErr = mu_voigtCorr.getError();
      resChangeV = (res-resCorr)/res * 100.;
      resChangeErrV = fabs(resChangeV) * sqrt(pow(resCorrErr/resCorr,2) + pow(resErr/res,2));
      double resS = doHistApprox ? sigma_eff(HmassUncorrHist)/mean_jj * 100. : sigma_eff(HmassUncorrArray,sizeOfHmassUncorrArray,mean_jj)/mean_jj * 100.;
      double resSerr = fabs(resS)*sqrt(pow(mean_jjErr/mean_jj,2)+pow(rms_jjErr/(resS*mean_jj/100.),2));
      double resSCorr = doHistApprox ? sigma_eff(HmassCorrHist)/mean_jjCorr * 100. : sigma_eff(HmassCorrArray,sizeOfHmassCorrArray,mean_jjCorr)/mean_jjCorr * 100.;

      double resSCorrerr = fabs(resSCorr)*sqrt(pow(mean_jjCorrErr/mean_jjCorr,2)+pow(rms_jjCorrErr/(resSCorr*mean_jjCorr/100.),2));
      resChangeS = (resS-resSCorr)/resS * 100.;
      resChangeErrS = fabs(resChangeS) * sqrt(pow(resSerr/resS,2) + pow(resSCorrerr/resSCorr,2));
      //plotParameters( c1, Rmass, HmassUncorrArray, sizeOfHmassUncorrArray, HmassCorrArray, sizeOfHmassCorrArray, datasetOdd.rmsVar(Hmass), datasetOdd.rmsVar(HmassCorr));
      plotParameters( c1, Rmass, &p, &pCorr);

      //c1->Print(Form("plots/resImp_voigt_M-%i_%i.png",Rmass,optimIndex));
      c1->Print(Form("plots/resImp_M-%i_%i.png",Rmass,optimIndex));
      c1->Clear();

      if(doEvenFile){
	RooRealVar mu_voigtTrain("mu_voigtTrain", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
	RooRealVar width_voigtTrain("width_voigtTrain", "width (corr.)", rms_jjCorr, 0.1*rms_jjCorr, 5.*rms_jjCorr, "GeV");
	RooRealVar sigma_voigtTrain("sigma_voigtTrain", "sigma (corr.)", 5., 0.0001, 20., "GeV");
	RooVoigtian voigtTrain("voigtTrain", "voigtTrain", HmassCorr, mu_voigtTrain, width_voigtTrain, sigma_voigtTrain);

	RooPlot * jjTrain_frame = HmassCorr.frame(lowVal,highVal);
	datasetEven.plotOn(jjTrain_frame, LineColor(kGreen+2), MarkerColor(kGreen), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
	RooFitResult *fTrain = voigtTrain.fitTo(datasetEven, Save(),Range(fitRangeLow,fitRangeHigh));
	voigtTrain.plotOn(jjTrain_frame, LineColor(kGreen+2), LineWidth(3));
	RooArgList pTrain = fTrain->floatParsFinal();
	jjCorr_frame->Draw();
	jjTrain_frame->Draw("same");
	jj_frame->Draw("same");

	plotParameters( &pTrain, &pCorr, c1, Rmass, 1, 0, KSprob);
	//c1->Print(Form("plots/ksTest_voigt_M-%i_%i.png",Rmass,optimIndex));
	c1->Print(Form("plots/ksTest_M-%i_%i.png",Rmass,optimIndex));
	c1->Clear();
      }
    }
  else
    {
      RooRealVar mu_gauss("mu_gauss", "mean (uncorr.)", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
      RooRealVar sigma1_gauss("sigma1_gauss", "sigma1 (uncorr.)", rms_jj, 0.01*rms_jj, 1000.*rms_jj, "GeV");
      RooRealVar sigma2_gauss("sigma2_gauss", "sigma2 (uncorr.)", rms_jj, 0.01*rms_jj, 1000.*rms_jj, "GeV");
      RooRealVar sigma3_gauss("sigma3_gauss", "sigma3 (uncorr.)", rms_jj, 0.01*rms_jj, 1000.*rms_jj, "GeV");
      RooGaussian gauss1("gauss1", "gauss1", Hmass, mu_gauss, sigma1_gauss);
      RooGaussian gauss2("gauss2", "gauss2", Hmass, mu_gauss, sigma2_gauss);
      RooGaussian gauss3("gauss3", "gauss3", Hmass, mu_gauss, sigma3_gauss);
      RooRealVar b1("b1","b1",1.0,0.0,sizeOfHmassCorrArray*100);
      RooRealVar b2("b2","b2",1.0,0.0,sizeOfHmassCorrArray*100);
      RooRealVar b3("b3","b3",1.0,0.0,sizeOfHmassCorrArray*100);
      RooAddPdf doubleGauss("dg","dg",RooArgList(gauss1,gauss2,gauss3),RooArgList(b1,b2,b3));

      RooRealVar mu_gaussCorr("mu_gaussCorr", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
      RooRealVar sigma1_gaussCorr("sigma1_gaussCorr", "sigma1 (corr.)", rms_jjCorr, 0.01*rms_jjCorr, 1000.*rms_jjCorr, "GeV");
      RooRealVar sigma2_gaussCorr("sigma2_gaussCorr", "sigma2 (corr.)", rms_jjCorr, 0.01*rms_jjCorr, 1000.*rms_jjCorr, "GeV");
      RooRealVar sigma3_gaussCorr("sigma3_gaussCorr", "sigma3 (corr.)", rms_jjCorr, 0.01*rms_jjCorr, 1000.*rms_jjCorr, "GeV");
      RooGaussian gauss1Corr("gauss1Corr", "gauss1Corr", HmassCorr, mu_gaussCorr, sigma1_gaussCorr);
      RooGaussian gauss2Corr("gauss2Corr", "gauss2Corr", HmassCorr, mu_gaussCorr, sigma2_gaussCorr);
      RooGaussian gauss3Corr("gauss3Corr", "gauss3Corr", HmassCorr, mu_gaussCorr, sigma3_gaussCorr);
      RooRealVar d1("d1","d1",1.0,0.0,sizeOfHmassCorrArray*100);
      RooRealVar d2("d2","d2",1.0,0.0,sizeOfHmassCorrArray*100);
      RooRealVar d3("d3","d3",1.0,0.0,sizeOfHmassCorrArray*100);
      RooAddPdf doubleGaussCorr("dgCorr","dgCorr",RooArgList(gauss1Corr,gauss2Corr,gauss3Corr),RooArgList(d1,d2,d3));

      RooPlot * jj_frame = Hmass.frame(lowVal,highVal);
      RooPlot * jjCorr_frame = HmassCorr.frame(lowVal,highVal);
      datasetOdd.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      datasetOdd.plotOn(jjCorr_frame, LineColor(kBlue+2), MarkerColor(kBlue), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      RooFitResult *f = doubleGauss.fitTo(datasetOdd, Save(), Range(lowVal,highVal));
      doubleGauss.plotOn(jj_frame, LineColor(kRed+2), LineWidth(3));
      RooFitResult * fCorr = doubleGaussCorr.fitTo(datasetOdd, Save(), Range(lowVal,highVal));
      doubleGaussCorr.plotOn(jjCorr_frame, LineColor(kBlue+2), LineWidth(3));
      RooArgList p = f->floatParsFinal();
      RooArgList pCorr = fCorr->floatParsFinal();
      jjCorr_frame->GetYaxis()->SetTitleOffset(1.25);
      jjCorr_frame->SetTitle("");
      jjCorr_frame->Draw();
      jj_frame->Draw("same");
      double f1 = ((RooRealVar*)f->floatParsFinal().at(f->floatParsFinal().index("b1")))->getVal();
      double f2 = ((RooRealVar*)f->floatParsFinal().at(f->floatParsFinal().index("b2")))->getVal();
      double f3 = ((RooRealVar*)f->floatParsFinal().at(f->floatParsFinal().index("b3")))->getVal();
      double f1Err = ((RooRealVar*)f->floatParsFinal().at(f->floatParsFinal().index("b1")))->getError();
      double f2Err = ((RooRealVar*)f->floatParsFinal().at(f->floatParsFinal().index("b2")))->getError();
      double f3Err = ((RooRealVar*)f->floatParsFinal().at(f->floatParsFinal().index("b3")))->getError();
      double g1 = ((RooRealVar*)fCorr->floatParsFinal().at(fCorr->floatParsFinal().index("d1")))->getVal();
      double g2 = ((RooRealVar*)fCorr->floatParsFinal().at(fCorr->floatParsFinal().index("d2")))->getVal();
      double g3 = ((RooRealVar*)fCorr->floatParsFinal().at(fCorr->floatParsFinal().index("d3")))->getVal();
      double g1Err = ((RooRealVar*)fCorr->floatParsFinal().at(fCorr->floatParsFinal().index("d1")))->getError();
      double g2Err = ((RooRealVar*)fCorr->floatParsFinal().at(fCorr->floatParsFinal().index("d2")))->getError();
      double g3Err = ((RooRealVar*)fCorr->floatParsFinal().at(fCorr->floatParsFinal().index("d3")))->getError();
      double res= (f1*sigma1_gauss.getVal()+f2*sigma2_gauss.getVal()+f3*sigma3_gauss.getVal())/mu_gauss.getVal() * 100.0;
      double resErr= res * sqrt(pow(mu_gauss.getError()/mu_gauss.getVal(),2) + pow(sigma1_gauss.getError()/sigma1_gauss.getVal(),2) + pow(sigma2_gauss.getError()/sigma2_gauss.getVal(),2) + pow(sigma3_gauss.getError()/sigma3_gauss.getVal(),2) + pow(f1Err/f1,2) + pow(f2Err/f2,2) + pow(f3Err/f3,2));
      double resCorr= (g1*sigma1_gaussCorr.getVal()+g2*sigma2_gaussCorr.getVal()+g3*sigma3_gaussCorr.getVal())/mu_gaussCorr.getVal() * 100.0;
      double resCorrErr= resCorr * sqrt(pow(mu_gaussCorr.getError()/mu_gaussCorr.getVal(),2) + pow(sigma1_gaussCorr.getError()/sigma1_gaussCorr.getVal(),2) + pow(sigma2_gaussCorr.getError()/sigma2_gaussCorr.getVal(),2) + pow(sigma3_gaussCorr.getError()/sigma3_gaussCorr.getVal(),2) + pow(g1Err/g1,2) + pow(g2Err/g2,2) + pow(g3Err/g3,2));
      newMean = mu_gaussCorr.getVal();
      newMeanErr = mu_gaussCorr.getError();
      resChangeV = (res-resCorr)/res * 100.;
      resChangeErrV = fabs(resChangeV) * sqrt(pow(resCorrErr/resCorr,2) + pow(resErr/res,2));
      double resS = doHistApprox ? sigma_eff(HmassUncorrHist)/mean_jj * 100. : sigma_eff(HmassUncorrArray,sizeOfHmassUncorrArray,mean_jj)/mean_jj * 100.;
      double resSerr = fabs(resS)*sqrt(pow(mean_jjErr/mean_jj,2)+pow(rms_jjErr/(resS*mean_jj/100.),2));
      double resSCorr = doHistApprox ? sigma_eff(HmassCorrHist)/mean_jjCorr * 100. : sigma_eff(HmassCorrArray,sizeOfHmassCorrArray,mean_jjCorr)/mean_jjCorr * 100.;
      double resSCorrerr = fabs(resSCorr)*sqrt(pow(mean_jjCorrErr/mean_jjCorr,2)+pow(rms_jjCorrErr/(resSCorr*mean_jjCorr/100.),2));
      resChangeS = (resS-resSCorr)/resS * 100.;
      resChangeErrS = fabs(resChangeS) * sqrt(pow(resSerr/resS,2) + pow(resSCorrerr/resSCorr,2));
      plotParameters( c1, Rmass, HmassUncorrArray, sizeOfHmassUncorrArray, HmassCorrArray, sizeOfHmassCorrArray, datasetOdd.rmsVar(Hmass), datasetOdd.rmsVar(HmassCorr));
      //c1->Print(Form("plots/resImp_voigt_M-%i_%i.png",Rmass,optimIndex));
      c1->Print(Form("plots/resImp_M-%i_%i.png",Rmass,optimIndex));
      c1->Clear();

      if(doEvenFile){
	RooRealVar mu_voigtTrain("mu_voigtTrain", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
	RooRealVar width_voigtTrain("width_voigtTrain", "width (corr.)", rms_jjCorr, 0.1*rms_jjCorr, 5.*rms_jjCorr, "GeV");
	RooRealVar sigma_voigtTrain("sigma_voigtTrain", "sigma (corr.)", 5., 0.0001, 20., "GeV");
	RooVoigtian voigtTrain("voigtTrain", "voigtTrain", HmassCorr, mu_voigtTrain, width_voigtTrain, sigma_voigtTrain);

	RooPlot * jjTrain_frame = HmassCorr.frame(lowVal,highVal);
	datasetEven.plotOn(jjTrain_frame, LineColor(kGreen+2), MarkerColor(kGreen), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
	RooFitResult *fTrain = voigtTrain.fitTo(datasetEven, Save());
	voigtTrain.plotOn(jjTrain_frame, LineColor(kGreen+2), LineWidth(3));
	RooArgList pTrain = fTrain->floatParsFinal();
	jjCorr_frame->Draw();
	jjTrain_frame->Draw("same");
	jj_frame->Draw("same");

	plotParameters( &pTrain, &pCorr, c1, Rmass, 1, 0, KSprob);
	//c1->Print(Form("plots/ksTest_voigt_M-%i_%i.png",Rmass,optimIndex));
	c1->Print(Form("plots/ksTest_M-%i_%i.png",Rmass,optimIndex));
	c1->Clear();
      }
    }
 

  //switch(Rmass){
  //case 300 : mean_jj=mean_jjCorr=120.4; rms_jj=rms_jjCorr=15.0; break;
  //case 500 : mean_jj=mean_jjCorr=125.4; rms_jj=rms_jjCorr=15.0; break;
  //case 700 : mean_jj=mean_jjCorr=129.5; rms_jj=rms_jjCorr=15.0; break;
  //case 1000: mean_jj=mean_jjCorr=134.2; rms_jj=rms_jjCorr=15.0; break;
  //}
  
  if(GAUSS)
    {
      RooRealVar mu_gauss("mu_gauss", "mean (uncorr.)", mean_jj, mean_jj-rms_jj, mean_jj+rms_jj, "GeV");
      RooRealVar sigma_gauss("sigma_gauss", "sigma (uncorr.)", rms_jj, .01*rms_jj, 5.*rms_jj, "GeV");
      RooGaussian gauss("gauss", "gauss", Hmass, mu_gauss, sigma_gauss);

      RooRealVar mu_gaussCorr("mu_gaussCorr", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
      RooRealVar sigma_gaussCorr("sigma_gaussCorr", "sigma (corr.)", rms_jjCorr, .01*rms_jjCorr, 5.*rms_jjCorr, "GeV");
      RooGaussian gaussCorr("gaussCorr", "gaussCorr", HmassCorr, mu_gaussCorr, sigma_gaussCorr);

      RooPlot * jj_frame = Hmass.frame(lowVal,highVal);
      RooPlot * jjCorr_frame = HmassCorr.frame(lowVal,highVal);
      datasetOdd.plotOn(jj_frame, LineColor(kRed+2), MarkerColor(kRed), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      datasetOdd.plotOn(jjCorr_frame, LineColor(kBlue+2), MarkerColor(kBlue), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
      RooFitResult *f = gauss.fitTo(datasetOdd, Save(), Range(mean_jj-25, mean_jj+20));
      gauss.plotOn(jj_frame, LineColor(kRed+2), LineWidth(3));
      RooFitResult * fCorr = gaussCorr.fitTo(datasetOdd, Save(), Range(mean_jjCorr-25, mean_jjCorr+20));
      gaussCorr.plotOn(jjCorr_frame, LineColor(kBlue+2), LineWidth(3));
      RooArgList p = f->floatParsFinal();
      RooArgList pCorr = fCorr->floatParsFinal();
      jjCorr_frame->GetYaxis()->SetTitleOffset(1.25);
      jjCorr_frame->SetTitle("");
      jjCorr_frame->Draw();
      jj_frame->Draw("same");
      double res = sigma_gauss.getVal() / mu_gauss.getVal() * 100.;
      double resErr= res * sqrt(pow(mu_gauss.getError()/mu_gauss.getVal(),2) + pow(sigma_gauss.getError()/sigma_gauss.getVal(),2));
      double resCorr = sigma_gaussCorr.getVal() / mu_gaussCorr.getVal() * 100.;
      double resCorrErr= res * sqrt(pow(mu_gaussCorr.getError()/mu_gaussCorr.getVal(),2) + pow(sigma_gaussCorr.getError()/sigma_gaussCorr.getVal(),2));

      resChangeG = (res-resCorr)/res * 100.;
      resChangeErrG = fabs(resChangeG) * sqrt(pow(resCorrErr/resCorr,2) + pow(resErr/res,2));
      plotParameters( &p, &pCorr, c1, Rmass, 0, 1, -1);
      //c1->Print(Form("plots/resImp_gauss_M-%i_%i.png",Rmass,optimIndex));
      c1->Clear();

      if(doEvenFile){
	RooRealVar mu_gaussTrain("mu_gaussTrain", "mean (corr.)", mean_jjCorr, mean_jjCorr-rms_jjCorr, mean_jjCorr+rms_jjCorr, "GeV");
	RooRealVar sigma_gaussTrain("sigma_gaussTrain", "sigma (corr.)", rms_jjCorr, .01*rms_jjCorr, 5.*rms_jjCorr, "GeV");
	RooGaussian gaussTrain("gaussTrain", "gaussTrain", HmassCorr, mu_gaussTrain, sigma_gaussTrain);

	RooPlot * jjTrain_frame = HmassCorr.frame(lowVal,highVal);
	datasetEven.plotOn(jjTrain_frame, LineColor(kGreen+2), MarkerColor(kGreen), MarkerStyle(8), MarkerSize(1.0), Binning(nbins,lowVal,highVal));
	RooFitResult * fTrain = gaussTrain.fitTo(datasetEven, Save(), Range(mean_jjCorr-25, mean_jjCorr+20));
	gaussTrain.plotOn(jjTrain_frame, LineColor(kGreen+2), LineWidth(3));
	RooArgList pTrain = fTrain->floatParsFinal();
	jjCorr_frame->Draw();
	jjTrain_frame->Draw("same");
	jj_frame->Draw("same");
	plotParameters( &pTrain, &pCorr, c1, Rmass, 0, 0, KSprob);
	//c1->Print(Form("plots/ksTest_gauss_M-%i_%i.png",Rmass,optimIndex));
	c1->Clear();
      }
    }

  if(doGen){
    TH1F *shape_gen = new TH1F("shape_gen","shape_gen",nbins,lowVal,highVal);
    tOdd->Project("shape_gen","HmassGen",testingCut);
    shape_gen->Sumw2();
    tOddCut->Draw("HmassGen",testingCut);
    double* HmassGenArray = tOddCut->GetV1();
    int sizeOfHmassGenArray = tOddCut->GetEntries(testingCut+"HmassGen>60 && HmassGen<180");
    double meanHmassGen=-1;
    sortArray(HmassGenArray,sizeOfHmassGenArray);
    double sigmaEffHmassGen = sigma_eff(HmassGenArray,sizeOfHmassGenArray,meanHmassGen);

    char *xLabel = (char*) "M(b#bar{b}) [GeV]";
    char *yLabel = (char*) "Events / (3 GeV)";

    gStyle->SetOptStat(0000);
    gStyle->SetOptTitle(0);

    shape_gen->Draw("HIST");
    shape_gen->GetXaxis()->SetTitle(xLabel);
    shape_gen->GetXaxis()->SetLabelSize(0.045);
    shape_gen->GetXaxis()->SetTitleSize(0.045);
    shape_gen->GetXaxis()->SetTitleOffset(1.1);
    shape_gen->GetYaxis()->SetTitle(yLabel);
    shape_gen->GetYaxis()->SetLabelSize(0.04);
    shape_gen->GetYaxis()->SetTitleSize(0.042);
    shape_gen->GetYaxis()->SetTitleOffset(1.55);

    TLatex *l11 = new TLatex;
    l11->SetNDC();
    l11->SetTextAlign(12);
    l11->SetTextSize(0.045);
    l11->DrawLatex(0.17,0.89,"CMS Simulation");
    l11->SetTextSize(0.035);
    l11->DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
    switch(doAltSample){
    case 0: 
      l11->DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
      l11->DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
      break;
    case 1: l11->DrawLatex(0.17,0.80,"ZH#rightarrow #ell#bar{#ell} b#bar{b}"); break;
    case 2: l11->DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
    case 3: l11->DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
    case 4: l11->DrawLatex(0.17,0.80,"ZZ#rightarrow #ell#bar{#ell} b#bar{b}"); break;
    }

    TLatex *l12 = new TLatex;
    l12->SetNDC();
    l12->SetTextAlign(12);
    l12->SetTextSize(0.025);
    l12->SetTextColor(kBlue+2);
    l12->DrawLatex(0.17,0.70,"Generated Jet M(b#bar{b}):");
    //l12->DrawLatex(0.17,0.67,Form("  Mean = %.1f #pm %.1f",shape_gen->GetMean(),shape_gen->GetMeanError()));
    //l12->DrawLatex(0.17,0.64,Form("  Sigma = %.1f #pm %.1f",shape_gen->GetRMS(),shape_gen->GetRMSError()));
    l12->DrawLatex(0.17,0.67,Form("  Mean = %.1f #pm %.1f",meanHmassGen,shape_gen->GetMeanError()));
    l12->DrawLatex(0.17,0.64,Form("  Sigma = %.1f #pm %.1f",sigmaEffHmassGen,shape_gen->GetRMSError()));
    double resGen = sigmaEffHmassGen/meanHmassGen * 100.;
    double resGenErr = resGen * sqrt(pow(shape_gen->GetRMSError()/shape_gen->GetRMS(),2) + pow(shape_gen->GetMeanError()/shape_gen->GetMean(),2));
    l12->DrawLatex(0.17,0.61,Form("  Res. = %.1f%% #pm %.1f%%",resGen,resGenErr));

    l11->Draw();
    l12->Draw();
    c1->SaveAs(Form("plots/genMbb_M-%i.png",Rmass));
  }

  cout << "\n\n";
  cout << "Resolution change from SigmaEff: " << resChangeS << " +/- " << resChangeErrS << endl; 
  cout << "Resolution change from Voigtian: " << resChangeV << " +/- " << resChangeErrV << endl; 
  cout << "Resolution change from Gaussian: " << resChangeG << " +/- " << resChangeErrG << endl; 
  cout << "New mean: " << newMean << " +/- " << newMeanErr << endl;
  if(doEvenFile) cout << "KS probability: " << KSprob << endl; 

  FILE* output = fopen("results.txt","a");
  //fprintf(output,"%3i %4i %5.2f +/- %4.2f :: %5.2f +/- %4.2f :: %5.2f +/- %4.2f :: %5.1f +/- %3.1f :: %9.7f \n",optimIndex,Rmass,resChangeS,resChangeErrS,resChangeV,resChangeErrV,resChangeG,resChangeErrG,newMean,newMeanErr,KSprob);
  fprintf(output,"%3i %4i %5.2f +/- %4.2f :: %5.1f +/- %3.1f :: %9.7f :: %i \n",optimIndex,Rmass,resChangeV,resChangeErrV,newMean,newMeanErr,KSprob,tOddCut->GetEntries("HmassCorr>500"));
  fclose(output);

  fOdd->Close();
  if(doEvenFile) fEven->Close();
  temp->Close();

  return;
}

void plotParameters(RooArgList *param, RooArgList *paramCorr, TCanvas *c, int Rmass, bool isVoigtian, bool isRes, double KSprob )
{
  c->cd(1);

  double mean1, mean2, sigma1, sigma2, width1, width2, res1, res2;
  double mean1err, mean2err, sigma1err, sigma2err, width1err, width2err, res1err, res2err;

  RooRealVar* obj = new RooRealVar();
  TIterator *it = (TIterator*) param->createIterator();
  obj = (RooRealVar*)it->Next();
  mean1 = obj->getVal();
  mean1err = obj->getError();
  obj = (RooRealVar*)it->Next();
  sigma1 = obj->getVal();
  sigma1err = obj->getError();
  if(isVoigtian){
    obj = (RooRealVar*)it->Next();
    width1 = obj->getVal();
    width1err = obj->getError();
  }
  it = (TIterator*) paramCorr->createIterator();
  obj = (RooRealVar*)it->Next();
  mean2 = obj->getVal();
  mean2err = obj->getError();
  obj = (RooRealVar*)it->Next();
  sigma2 = obj->getVal();
  sigma2err = obj->getError();
  if(isVoigtian){
    obj = (RooRealVar*)it->Next();
    width2 = obj->getVal();
    width2err = obj->getError();
  }

  if(isVoigtian){
    double fg = 2. * sigma1 * sqrt(2. * log(2));
    double fl = 2. * width1;
    double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
    res1 = fv / mean1 * 100. / (2. * sqrt(2. * log(2.)));
    res1err = res1 * sqrt(pow(mean1err/mean1,2) + pow(sigma1err/sigma1,2) + pow(width1err/width1,2));

    fg = 2. * sigma2 * sqrt(2. * log(2));
    fl = 2. * width2;
    fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
    res2 = fv / mean2 * 100. / (2. * sqrt(2. * log(2.)));
    res2err = res2 * sqrt(pow(mean2err/mean2,2) + pow(sigma2err/sigma2,2) + pow(width2err/width2,2));
  }
  else{
    res1 = sigma1/mean1 * 100.;
    res1err = res1 * sqrt(pow(mean1err/mean1,2) + pow(sigma1err/sigma1,2));
    res2 = sigma2/mean2 * 100.;
    res2err = res2 * sqrt(pow(mean2err/mean2,2) + pow(sigma2err/sigma2,2));
  }

  double resChange = fabs(res1-res2)/res1 * 100.;
  double resChangeErr = resChange * sqrt(pow(res1err/res1,2) + pow(res2err/res2,2));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.17,0.89,"CMS Simulation");
  latexLabel.SetTextSize(0.035);
  latexLabel.DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
  switch(doAltSample){
  case 0:
    latexLabel.DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
    latexLabel.DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
    break;
  case 1: latexLabel.DrawLatex(0.17,0.80,"ZH#rightarrow l#bar{l} b#bar{b}"); break;
  case 2: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
  case 3: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
  case 4: latexLabel.DrawLatex(0.17,0.80,"ZZ#rightarrow l#bar{l} b#bar{b}"); break;
  }

  int pos=0;
  float topPos=0.89, spacing=0.03;
  latexLabel.SetTextSize(0.025);
  if(isRes){
    latexLabel.SetTextColor(kRed+2);
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Before correction:");
  }
  else{
    latexLabel.SetTextColor(kGreen+2);
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction, training:");    
  }
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean1,mean1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma1,sigma1err));
  if(isVoigtian) latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width1,width1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res1,res1err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  if(isRes){
    latexLabel.SetTextColor(kBlue+2);
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction:");
  }
  else{
    latexLabel.SetTextColor(kBlue+2);
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction, testing:");    
  }
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean2,mean2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma2,sigma2err));
  if(isVoigtian) latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width2,width2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res2,res2err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlack);
  if(isRes){
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Resolution Change:");
    if(res1>res2)
      latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  + %.1f%% #pm %.1f%%",resChange,resChangeErr));
    else
      latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  - %.1f%% #pm %.1f%%",resChange,resChangeErr));
  }
  else{
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"KS Probability:");
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  %.7f%%",KSprob));
  }

  latexLabel.Draw();

  return;
}

void plotParameters(TCanvas *c, int Rmass, RooArgList *param, RooArgList *paramCorr)
{
  c->cd(1);

  double mean1, mean2, sigma1, sigma2, width1, width2, res1, res2;
  double mean1err, mean2err, sigma1err, sigma2err, width1err, width2err, res1err, res2err;

  RooRealVar* obj = new RooRealVar();
  TIterator *it = (TIterator*) param->createIterator();
  obj = (RooRealVar*)it->Next();
  mean1 = obj->getVal();
  mean1err = obj->getError();
  obj = (RooRealVar*)it->Next();
  sigma1 = obj->getVal();
  sigma1err = obj->getError();
  obj = (RooRealVar*)it->Next();
  width1 = obj->getVal();
  width1err = obj->getError();
  
  it = (TIterator*) paramCorr->createIterator();
  obj = (RooRealVar*)it->Next();
  mean2 = obj->getVal();
  mean2err = obj->getError();
  obj = (RooRealVar*)it->Next();
  sigma2 = obj->getVal();
  sigma2err = obj->getError();
  obj = (RooRealVar*)it->Next();
  width2 = obj->getVal();
  width2err = obj->getError();

  double fg = 2. * sigma1 * sqrt(2. * log(2));
  double fl = 2. * width1;
  double fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
  res1 = fv / mean1 * 100. / (2. * sqrt(2. * log(2.)));
  res1err = res1 * sqrt(pow(mean1err/mean1,2) + pow(sigma1err/sigma1,2) + pow(width1err/width1,2));

  fg = 2. * sigma2 * sqrt(2. * log(2));
  fl = 2. * width2;
  fv = 0.5346 * fl + sqrt(0.2166 * pow(fl, 2.) + pow(fg, 2.));
  res2 = fv / mean2 * 100. / (2. * sqrt(2. * log(2.)));
  res2err = res2 * sqrt(pow(mean2err/mean2,2) + pow(sigma2err/sigma2,2) + pow(width2err/width2,2));

  double resChange = fabs(res1-res2)/res1 * 100.;
  double resChangeErr = resChange * sqrt(pow(res1err/res1,2) + pow(res2err/res2,2));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.17,0.89,"CMS Simulation");
  latexLabel.SetTextSize(0.035);
  latexLabel.DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
  switch(doAltSample){
  case 0:
    latexLabel.DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
    latexLabel.DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
    break;
  case 1: latexLabel.DrawLatex(0.17,0.80,"ZH#rightarrow l#bar{l} b#bar{b}"); break;
  case 2: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
  case 3: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
  case 4: latexLabel.DrawLatex(0.17,0.80,"ZZ#rightarrow l#bar{l} b#bar{b}"); break;
  }

  int pos=0;
  float topPos=0.89, spacing=0.03;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kRed+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Before correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean1,mean1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma1,sigma1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width1,width1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res1,res1err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlue+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean2,mean2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sigma = (%.1f #pm %.1f) GeV",sigma2,sigma2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Width = (%.1f #pm %.1f) GeV",width2,width2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res2,res2err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Resolution Change:");
  if(res1>res2)
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  + %.1f%% #pm %.1f%%",resChange,resChangeErr));
  else
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  - %.1f%% #pm %.1f%%",resChange,resChangeErr));

  latexLabel.Draw();

  return;
}

void plotParameters(TCanvas *c, int Rmass, double *arrayUncorr, int nEventsUncorr, double *arrayCorr, int nEventsCorr, RooRealVar *sig1, RooRealVar *sig2 )
{
  c->cd(1);

  double mean1=-1, mean2=-1, sigma1, sigma2, res1, res2;
  double mean1err, mean2err, sigma1err, sigma2err, res1err, res2err;

  if(doHistApprox){
    sigma1 = sig1->getVal();
    sigma2 = sig2->getVal();
    for(int i=0; i<nEventsUncorr; ++i)
      mean1+=arrayUncorr[i];
    for(int i=0; i<nEventsCorr; ++i)
      mean2+=arrayCorr[i];
    mean1/=nEventsUncorr;
    mean2/=nEventsCorr;
  }
  else{
    sigma1 = sigma_eff(arrayUncorr, nEventsUncorr, mean1);
    sigma2 = sigma_eff(arrayCorr, nEventsCorr, mean2);
  }

  res1 = sigma1/mean1 * 100.;
  res2 = sigma2/mean2 * 100.;

  double tmp1=0.0, tmp2=0.0;
  for(int i=0; i<nEventsUncorr; ++i)
    tmp1+=pow(arrayUncorr[i]-mean1,2);
  for(int i=0; i<nEventsCorr; ++i)
    tmp2+=pow(arrayCorr[i]-mean2,2);
  mean1err = sqrt(tmp1)/(nEventsUncorr-1);
  mean2err = sqrt(tmp2)/(nEventsCorr-1);
  sigma1err = sig1->getError();
  sigma2err = sig2->getError();

  res1err = fabs(res1) * sqrt(pow(sigma1err/sigma1,2) + pow(mean1err/mean1,2));
  res2err = fabs(res2) * sqrt(pow(sigma2err/sigma2,2) + pow(mean2err/mean2,2));

  double resChange = fabs(res1-res2)/res1 * 100.;
  double resChangeErr = resChange * sqrt(pow(res1err/res1,2) + pow(res2err/res2,2));

  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.045);
  latexLabel.SetTextAlign(12);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.17,0.89,"CMS Simulation");
  latexLabel.SetTextSize(0.035);
  latexLabel.DrawLatex(0.17,0.84,"#sqrt{s} = 8 TeV");
  switch(doAltSample){
  case 0:
    latexLabel.DrawLatex(0.17,0.80,"R#rightarrow H(#gamma#gamma)H(b#bar{b})");
    latexLabel.DrawLatex(0.17,0.75,TString::Format("M_{R} = %i GeV",Rmass));
    break;
  case 1: latexLabel.DrawLatex(0.17,0.80,"ZH#rightarrow l#bar{l} b#bar{b}"); break;
  case 2: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #gamma#gamma b#bar{b}"); break;
  case 3: latexLabel.DrawLatex(0.17,0.80,"HZ#rightarrow #nu#bar{#nu} b#bar{b}"); break;
  case 4: latexLabel.DrawLatex(0.17,0.80,"ZZ#rightarrow l#bar{l} b#bar{b}"); break;
  }

  int pos=0;
  float topPos=0.89, spacing=0.03;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kRed+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Before correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean1,mean1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sig. Eff. = (%.1f #pm %.1f) GeV",sigma1,sigma1err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res1,res1err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlue+2);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"After correction:");
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Mean = (%.1f #pm %.1f) GeV",mean2,mean2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Sig. Eff. = (%.1f #pm %.1f) GeV",sigma2,sigma2err));
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  Res. = %.1f%% #pm %.1f%%",res2,res2err));
  
  pos++;
  latexLabel.SetTextSize(0.025);
  latexLabel.SetTextColor(kBlack);
  latexLabel.DrawLatex(0.67,topPos-spacing*pos++,"Resolution Change:");
  if(res1>res2)
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  + %.1f%% #pm %.1f%%",resChange,resChangeErr));
  else
    latexLabel.DrawLatex(0.67,topPos-spacing*pos++,Form("  - %.1f%% #pm %.1f%%",resChange,resChangeErr));

  latexLabel.Draw();

  return;
}

void sortArray(double* a, const int len){
  double tmpVal;
  int tmpIndex;

  for(int i=0; i<len; ++i){
    tmpVal=a[i];
    tmpIndex=i;
    for(int j=i+1; j<len; ++j){
      if(a[j]<tmpVal){
	tmpVal=a[j];
	tmpIndex=j;
      }
    }
    if(tmpIndex!=i){
      a[tmpIndex]=a[i];
      a[i]=tmpVal;
    }
  }

  cout<<"!!! "<<len<<' '<<a[0]<<' '<<a[len-1]<<' '<<endl;
}

double sigma_eff(double *data, int nEvents, double &mean){
  cout<<"@@@ "<<nEvents<<' '<<data[0]<<' '<<data[nEvents-1]<<' '<<endl;
  sortArray(data,nEvents);
  cout<<"@@@ "<<nEvents<<' '<<data[0]<<' '<<data[nEvents-1]<<' '<<endl;

  if(mean<0){
    mean=0.0;
    for(int i=0; i<nEvents; ++i) mean+=data[i];
    mean/=nEvents;
  }

  double width=1.0, step=0.01, frac=0.0;
  do{
    width+=step;
    int nEventsPass=0;
    for(int i=0; i<nEvents && data[i]<mean+width/2.0; ++i){
      if(data[i]>mean-width/2.0 && data[i]<mean+width/2.0) ++nEventsPass;
    }
    frac=(double)nEventsPass/nEvents;
  }
  while(frac<0.6827);

  return width/2.0;
}

double sigma_eff(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

