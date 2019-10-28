#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

using namespace std;

//R. K. Verma
//Sat Jul 14 14:47:08 IST 2018

///////////////////////////////////////////
//CHANNEL
bool isMuChannel = false;
bool isEleChannel = true;

//SAVE HISTOS ON DISK
bool isSaveHisto = true;

TFile *fDYdd = new TFile("all_DY_dd.root","RECREATE");
///////////////////////////////////////////


//INPUT FILES
TFile* fMuData  = TFile::Open("inputFiles/all_muData.root");
TFile* fEleData = TFile::Open("inputFiles/all_eleData.root");
TFile* fBkg	= TFile::Open("inputFiles/all_Bkg.root");
TFile* fVV	    = TFile::Open("inputFiles/all_VV.root");
TFile* fDY	    = TFile::Open("inputFiles/all_DY.root");
TFile* fTT	    = TFile::Open("inputFiles/all_TT.root");
TFile* fST	    = TFile::Open("inputFiles/all_ST.root");
TFile* fWJ	    = TFile::Open("inputFiles/all_WJets.root");
TFile *fSigMuMuZ250= TFile::Open("inputFiles/all_ExLepMuMuZ_M250.root");
TFile *fSigEEZ250  = TFile::Open("inputFiles/all_ExLepEEZ_M250.root");
TFile *fSigMuMuZ1000= TFile::Open("inputFiles/all_ExLepMuMuZ_M1000.root");
TFile *fSigEEZ1000  = TFile::Open("inputFiles/all_ExLepEEZ_M1000.root");
TFile *fSigMuMuZ2000= TFile::Open("inputFiles/all_ExLepMuMuZ_M2000.root");
TFile *fSigEEZ2000  = TFile::Open("inputFiles/all_ExLepEEZ_M2000.root");
TFile *fSigMuMuZ5000= TFile::Open("inputFiles/all_ExLepMuMuZ_M5000.root");
TFile *fSigEEZ5000  = TFile::Open("inputFiles/all_ExLepEEZ_M5000.root");
TFile *fNonDYBkg= TFile::Open("inputFiles/all_NonDYBkg.root");

TString chDir = "baseLowMET";
class MyExLepStackHisto{
  public :
	//get histogram from root file. Return empty hist, if the hist does not exist.
	TH1F* getHisto(TFile *inRootFile, TString chDir, TString baseDir, TString histDir, TString histName, bool isMerge=false, int mergeFrom=10, int nbins=10);
	//rebin histo
	TH1F* binnnedHisto(TH1F* h, int mergeFrom, int nbins, vector<int>binRange);
	vector<int> pickBins(TH1F* h, int mergeFrom, int nbins);
        //decorate histogram
        TH1F* decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit);
	//function to stack histos
    	TH1F* stackHisto(TFile *inRootFile, TString chDir, TString baseDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, bool isMerge, int mergeFrom, int nbins);
	//function to add histograms
	TH1F*  addHistoForUnc(TString chDir, TString baseDir, TString histDir, TString histName);
	//Up/down error in the unc band
	double errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus);
	double errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus);
	//unc graph
	TGraphAsymmErrors *UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, bool isFullGraph = false, bool isRatioGraph = false);
        TPaveText *paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font );
        //qcd from data
    	vector<double> getTransFactDY(TFile* fData, TFile* fTT, TFile* fWJ, TFile* fVV, TString chDir, TString baseDir, TString histDir, TString histName, TString xTitle, double xmin, double xmax);
    	vector<double> getTransFactDY2(TFile* fData, TFile* fNonDYBkg, TString chDir, TString baseDir, TString histDir, TString histName, TString xTitle, double xmin, double xmax);
    	TH1F* getDataDrivenDY(TString chDir, TString baseDir, TString histDir, TString histName, double transFact=1.0, double errorTransFact = 0.0);
    	double getStatUnc(TH1F* hCentral, double sError = 0.0);
    	void makeHistoPositive(TH1F* hist, bool setErrorZero = false);
    	//function to overlap histograms
    	void overlapHisto(TH1F *h1, TH1F *h2, bool isRatio, TString histDir, TString histName, TString xTitle, double xmin, double xmax);
  private :
	int dont_use ;
};

//--------------------------------------------//
//define various functions
//--------------------------------------------//
TH1F* MyExLepStackHisto::binnnedHisto(TH1F* h, int mergeFrom, int nbins, vector<int>binRange){
  int finalBins = mergeFrom + binRange.size();
  Double_t xbins[nbins];
  TAxis *axis = h->GetXaxis();
  int i;
  for(i=0;i<=mergeFrom;i++)  xbins[i] = axis->GetBinLowEdge(1+i);
  for(int k = 0; k<int(binRange.size()); k++){
    xbins[k+mergeFrom+1] = axis->GetBinLowEdge(binRange[k]+1);
  }
  //for (i=0;i<=finalBins;i++) printf("xbins[%d]=%g\n",i,xbins[i]);
  TH1F *hnew = (TH1F*)h->Rebin(finalBins, "hnew",xbins);
  for (i=0;i<int(binRange.size());i++) {
    int binWidth = 1.0;
    if(binRange[i] == mergeFrom){
      cout<<"ERORR: change the xMax of "<<h->GetName()<<endl;
      exit(1);
    } 
    binWidth = binRange[i] - mergeFrom;
    int newBin = mergeFrom+i+1;
    double newBinCont = hnew->GetBinContent(newBin)/binWidth;
    double newBinErr = hnew->GetBinError(newBin)/binWidth;
    hnew->SetBinContent(newBin, newBinCont);
    hnew->SetBinError(newBin,   newBinErr);
  }
  return hnew;
}

vector<int> MyExLepStackHisto:: pickBins(TH1F* h, int mergeFrom, int nbins){
  //int nbins = h->GetNbinsX();
  int remainingBins = nbins - mergeFrom;
  //10%, 15%, 25%, 50%
  vector<int>percentBin;
  percentBin.push_back(int(remainingBins*0.10));
  percentBin.push_back(int(remainingBins*0.15));
  percentBin.push_back(int(remainingBins*0.25));
  int lastBin = remainingBins - (
	  int(remainingBins*0.10)+
	  int(remainingBins*0.15)+
	  int(remainingBins*0.25));
  percentBin.push_back(lastBin);
  vector<int>binRange;
  int count = mergeFrom;
  for(int i=0; i<int(percentBin.size()); i++){
    count = count + percentBin[i];
    binRange.push_back(count);
  }
  return binRange;
}
TH1F*  MyExLepStackHisto:: getHisto(TFile *inRootFile, TString chDir, TString baseDir, TString histDir, TString histName, bool isMerge=false, int mergeFrom=10, int nbins=10){
  TH1F* hist;
  TString fullPath = chDir+baseDir+histDir+histName;
  //TString fullPath = string(inRootFile->GetName())+"/"+chDir+baseDir+histDir+histName;
  string exception_msg ("The histogram path, "+fullPath+", does not exist");
  try{
    if(!(inRootFile->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    //cout<<"WARNING:"<<e<<endl;
  }
  try{
    if(!(fTT->Get(fullPath)))
       throw  exception_msg.c_str();
  }catch (const char *e){
    cout<<"\033[01;31mERROR: \033[0m"<<e<< endl;
    exit(0);
  }
  if(!(inRootFile->Get(fullPath))){
    hist = (TH1F*)(fTT->Get(fullPath))->Clone(histName);
    hist->Reset();
  }else hist = (TH1F*)(inRootFile->Get(fullPath))->Clone(histName);
  //merge bins
  TH1F * hNew;
  if(isMerge){
  vector<int>binRange = pickBins(hist, mergeFrom, nbins);
  hNew = (TH1F*) binnnedHisto(hist, mergeFrom, nbins, binRange);
  }
  if(isMerge) return hNew;
  else return hist;
}


TH1F* MyExLepStackHisto:: decorateHisto(TH1F* hist, TString myTit, TString xTit, TString yTit){
  hist->SetTitle(myTit);
  hist->GetXaxis()->SetTitle(xTit);
  hist->GetYaxis()->SetTitle(yTit);
  hist->GetYaxis()->SetTitleOffset(1.00);
  hist->GetXaxis()->SetTitleOffset(1.00);
  hist->GetYaxis()->SetTitleSize(0.10);
  hist->GetXaxis()->SetTitleSize(0.10);
  hist->GetXaxis()->SetLabelSize(0.10);
  hist->GetYaxis()->SetLabelSize(0.10);
  hist->GetXaxis()->SetTickLength(0.05);
  hist->GetXaxis()->SetNdivisions(10);
  hist->GetYaxis()->SetNdivisions(5);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTickLength(0.04);
  return hist;
}

TH1F *  MyExLepStackHisto:: stackHisto(TFile *inRootFile, TString chDir, TString baseDir, TString histDir, TString histName, int color, double scale, THStack* MuptStack, TH1F* hMC, bool isMerge, int mergeFrom, int nbins){
  TH1F* hist = getHisto(inRootFile, chDir, baseDir, histDir, histName, isMerge, mergeFrom,  nbins);
  hist->Scale(scale);
  hist->SetFillColor(color);
  //leg->AddEntry(hist,lable,"F");
  MuptStack->Add(hist);
  hMC->Add(hist);
  return hist;
}

double MyExLepStackHisto:: errBandUp(int iBin, TH1F *hCentral, TH1F *hJESPlus, TH1F *hJERPlus, TH1F *bTagPlus){
  double errUp = sqrt(pow(fabs(hJESPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) +
		  pow(fabs(hJERPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) +
		  pow(fabs(bTagPlus->GetBinContent(iBin+1) - hCentral->GetBinContent(iBin+1)),2) +
		  pow(hCentral->GetBinError(iBin+1),2));
  return errUp;
}

double MyExLepStackHisto:: errBandDown(int iBin, TH1F *hCentral, TH1F *hJESMinus, TH1F *hJERMinus, TH1F *bTagMinus){
  double errDown =sqrt(pow(fabs(hCentral->GetBinContent(iBin+1) - hJESMinus->GetBinContent(iBin+1)),2) +
		  pow(fabs(hCentral->GetBinContent(iBin+1) - hJERMinus->GetBinContent(iBin+1)),2) +
		  pow(fabs(hCentral->GetBinContent(iBin+1) - bTagMinus->GetBinContent(iBin+1)),2) +
		  pow(hCentral->GetBinError(iBin+1),2));
  return errDown;
}

TGraphAsymmErrors * MyExLepStackHisto:: UNCGRAPH(TH1F *hCentral, TH1F *hJESPlus, TH1F *hJESMinus, TH1F *hJERPlus, TH1F *hJERMinus, TH1F *bTagPlus, TH1F *bTagMinus, bool isFullGraph = false, bool isRatioGraph = false){
  TGraphAsymmErrors *gr;
  int n1 = hCentral->GetNbinsX();
  double *Yval, *errorU, *errorD, *XerrorU, *XerrorD, *Xval ;
  Yval = new double[n1]; errorU = new double[n1]; errorD = new double[n1];
  XerrorU=new double[n1]; XerrorD=new double[n1]; Xval=new double[n1];
  //cout << "No. of bins= " << n1 << endl;
  for(int i=0; i<n1; i++){
    if(isFullGraph){
    Yval[i]   = hCentral->GetBinContent(i+1);
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus);
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus);
    }
    if(isRatioGraph){
    Yval[i]   = 1;
    errorU[i] = errBandUp(i, hCentral, hJESPlus, hJERPlus, bTagPlus);
    errorD[i] = errBandDown(i, hCentral, hJESMinus, hJERMinus, bTagMinus);
    //cout<<"bin = "<<i<<endl;
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<hCentral->GetBinContent(i+1)<<endl;
    errorU[i] = errorU[i]/hCentral->GetBinContent(i+1);
    errorD[i] = errorD[i]/hCentral->GetBinContent(i+1);
    //cout<<Yval[i]<<"\t"<<errorU[i]<<"\t"<<hCentral->GetBinContent(i+1)<<endl;
    }
    Xval[i]   = hCentral->GetBinCenter(i+1);
    XerrorU[i]= hCentral->GetBinWidth(i+1)/2;
    XerrorD[i]= hCentral->GetBinWidth(i+1)/2;
  }
  gr = new TGraphAsymmErrors(n1, Xval, Yval, XerrorD, XerrorU, errorD, errorU);
  return gr;
  delete [] Yval; delete [] errorU; delete [] errorD; delete [] XerrorU; delete [] XerrorD; delete [] Xval;
}

TH1F*  MyExLepStackHisto:: addHistoForUnc(TString chDir, TString baseDir, TString histDir, TString histName){
  TH1F* hVV =   	getHisto(fVV,   chDir, baseDir, histDir, histName);
  TH1F* hDY =   	getHisto(fDY,    chDir, baseDir, histDir, histName);
  TH1F* hWJ =   	getHisto(fWJ,   chDir, baseDir, histDir, histName);
  TH1F* hTT =   	getHisto(fTT,   chDir, baseDir, histDir, histName);
  TH1F* hAll = (TH1F*)hVV->Clone("hAllMC");
  hAll->Add(hDY);
  hAll->Add(hWJ);
  hAll->Add(hTT);
  return hAll;
}

TPaveText * MyExLepStackHisto:: paveText(double minX, double minY, double maxX, double maxY, int lineColor, int fillColor, int size, int style, int font ){
  TPaveText *pt = new TPaveText(minX, minY, maxX, maxY, "brNDC"); // good_v1
  pt->SetBorderSize(size);
  pt->SetFillColor(fillColor);
  pt->SetFillStyle(style);
  pt->SetLineColor(lineColor);
  pt->SetTextFont(font);
  return pt;
}

void MyExLepStackHisto::makeHistoPositive(TH1F* hist, bool setErrorZero = false){
  for(int ibin=1; ibin<hist->GetNbinsX(); ibin++){
    double binCont = hist->GetBinContent(ibin);
    if(binCont<0){
      hist->SetBinContent(ibin, 0);
      if(setErrorZero) hist->SetBinError(ibin, 0);
    }
  }
}

double MyExLepStackHisto::getStatUnc(TH1F* hCentral, double sError = 0.0){
  double  norm = hCentral->IntegralAndError(1, hCentral->GetNbinsX(), sError);
  //double statUnc = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;
  double statUnc = sError;
  return statUnc;
}
void MyExLepStackHisto::overlapHisto(TH1F *h1, TH1F *h2, bool isRatio, TString histDir, TString histName, TString xTitle, double xmin, double xmax){
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);
  const float xpad[2] = {0.0, 1.0};
  const float ypad[3] = {0.0, 0.30,0.98};
  TCanvas *canv = new TCanvas();
  canv->Divide(1, 2);
  canv->cd(1);
  gPad->SetPad(xpad[0],ypad[1],xpad[1],ypad[2]);
  //gPad->SetTopMargin(1.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.0);
  //legend box
  TLegend* leg = new TLegend(0.7518792,0.6261504,0.9512081,0.9198861,NULL,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetFillColor(10); leg->SetTextSize(0.07);
  //pave text CMS box
  TPaveText *pt = new TPaveText(0.11,0.9354,0.90,0.9362, "brNDC"); // good_v1
  pt->SetBorderSize(1); pt->SetFillColor(19);
  pt->SetFillStyle(0); pt->SetTextSize(0.08);
  pt->SetLineColor(0); pt->SetTextFont(132);
  TText *text = pt->AddText(histDir+": #sqrt{s}=13 TeV, 35.9 fb^{-1}; ");
  text->SetTextAlign(11);
  //pave text channel box
  TPaveText *ch = new TPaveText(1.00,0.9154898,0.7510067,0.9762187,"brNDC");
  ch->SetFillColor(19); ch->SetFillStyle(0);
  ch->SetLineColor(0); ch->SetTextSize(0.08);
  ch->SetBorderSize(1);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //data-MC from isolated region
  leg->AddEntry(h1,"Region-C","LP");
  h1->SetMarkerColor(kRed);
  h1->SetTitle("");
  h1->SetLineColor(kRed);
  h1->Scale(1/h1->Integral());
  cout<<h1->GetMaximum()<<endl;
  h1->GetYaxis()->SetRangeUser(-0.01,  1.2*h1->GetMaximum());
  h1->GetXaxis()->SetRangeUser(xmin, xmax);
  h1->GetYaxis()->SetTitleOffset(0.90);
  h1->GetYaxis()->SetTitle("#splitline{Data - non DY background}{(normalized to 1)}");
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetTickLength(0.04);
  h1->GetYaxis()->SetNdivisions(5);
  h1->GetYaxis()->SetRangeUser(-0.5, 1.0);
  h1->Draw("e1");
  //data-MC from non-isolated region
  leg->AddEntry(h2,"Region-D","LP");
  h2->SetMarkerColor(kGreen);
  h2->SetLineColor(kGreen);
  h2->Scale(1/h2->Integral());
  h2->Draw("SAME");
  pt->Draw();
  leg->Draw();
  ch->Draw();
  canv->cd(2);
  gPad->SetPad(xpad[0],ypad[0],xpad[1], ypad[1]);
  gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0.5); gPad->SetGridy();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  TH1F *h12 = (TH1F*)h1->Clone("h12");
  h12->Reset();
  h12->Add(h1);
  h12->SetTitle("");
  h12->Divide(h2); h12->SetMarkerStyle(20); h12->SetMarkerSize(0.8);
  h12->SetMarkerColor(kBlack); h12->SetLineColor(kBlack); h12->GetYaxis()->SetRangeUser(-10, 10);
  //h12->GetXaxis()->SetRangeUser(xmin, xmax);
  h12->GetXaxis()->SetTickLength(0.13);
  h12->GetYaxis()->SetTickLength(0.04);
  h12->GetXaxis()->SetTitle(xTitle);
  h12->GetYaxis()->SetTitleOffset(0.50);
  h12->GetXaxis()->SetTitleOffset(1.10);
  h12->GetYaxis()->SetTitle("#frac{Region-C}{Region-D}");
  h12->GetYaxis()->CenterTitle();
  h12->GetYaxis()->SetTitleSize(0.12);
  h12->GetXaxis()->SetTitleSize(0.15);
  h12->GetXaxis()->SetLabelSize(0.10);
  h12->GetYaxis()->SetLabelSize(0.15);
  h12->GetYaxis()->SetNdivisions(5);
  h12->Draw("E"); // use "P" or "AP"
  //base line
  TF1 *baseLine = new TF1("baseLine","1", -100, 2000);
  baseLine->SetLineColor(kCyan+1);
  baseLine->Draw("SAME");
  if(isSaveHisto){
    TString outDir = "$PWD/QCD/"+histDir;
    TString histDir_ = "";
    if(histDir.Contains("BTag")) histDir_ = "BTag";
    if(histDir.Contains("KinFit")) histDir_ = "KinFit";
    mkdir("QCD/", S_IRWXU);
    mkdir("QCD/"+histDir, S_IRWXU);
    TString outFile(outDir);
    if(isMuChannel)outFile += "mu_"+histDir_+"_"+histName+".png";
    if(isEleChannel)outFile += "ele_"+histDir_+"_"+histName+".png";
    canv->SaveAs(outFile);
    //canv->Close();
  }
}

vector<double> MyExLepStackHisto::getTransFactDY(TFile* fData, TFile* fTT, TFile* fWJ, TFile* fVV, TString chDir, TString baseDir, TString histDir, TString histName, TString xTitle, double xmin, double xmax){
  //RegionC = LowMET, Iso
  //if(histName.Contains("CTagEx")) histName="mjj_kfit";
  TH1F* hVV_RegC = getHisto(fVV,   chDir, baseDir, histDir, "C_"+histName);//Reg = Region
  TH1F* hWJ_RegC = getHisto(fWJ,   chDir, baseDir, histDir, "C_"+histName);
  TH1F* hTT_RegC = getHisto(fTT,   chDir, baseDir, histDir, "C_"+histName);
  TH1F* hMC_RegC = (TH1F*)hVV_RegC->Clone("hAllMC");
  hMC_RegC->Add(hWJ_RegC);
  hMC_RegC->Add(hTT_RegC);
  TH1F* hData_RegC= (TH1F*) getHisto(fData, chDir, baseDir, histDir, "C_"+histName);
  //RegionD = LowMET, NonIso
  TH1F* hVV_RegD = getHisto(fVV,   chDir, baseDir, histDir, "D_"+histName);
  TH1F* hWJ_RegD = getHisto(fWJ,   chDir, baseDir, histDir, "D_"+histName);
  TH1F* hTT_RegD = getHisto(fTT,   chDir, baseDir, histDir, "D_"+histName);
  TH1F* hMC_RegD = (TH1F*)hVV_RegD->Clone("hAllMC");
  hMC_RegD->Add(hWJ_RegD);
  hMC_RegD->Add(hTT_RegD);
  TH1F* hData_RegD=  getHisto(fData, chDir, baseDir, histDir, "D_"+histName);
  TH1F* hDiffC = (TH1F*)hData_RegC->Clone("hDiffC");
  hDiffC->Add(hMC_RegC, -1);
  TH1F* hDiffD = (TH1F*)hData_RegD->Clone("hDiffD");
  hDiffD->Add(hMC_RegD, -1);
  //If binContent < 0, set it to 0
  makeHistoPositive(hDiffC, false);
  makeHistoPositive(hDiffD, false);
  double intDiffC   = hDiffC->Integral();
  double errDiffC   = getStatUnc(hDiffC, 0.0);
  double intDiffD   = hDiffD->Integral();
  double errDiffD   = getStatUnc(hDiffD, 0.0);
  //Ratio of (Data-MC) from RegionD and RegionC
  double ratioDiffCD = intDiffC/intDiffD;
  double tmpC = errDiffC/intDiffC;
  double tmpD = errDiffD/intDiffD;
  double errDiffCD = ratioDiffCD*sqrt(tmpD*tmpD + tmpC*tmpC);
  cout<<"\n--------------------------"<<endl;
  printf("intC =  %f, intD = %f", intDiffC, intDiffD);
  vector<double>sfAndErr;
  sfAndErr.push_back(ratioDiffCD);
  sfAndErr.push_back(errDiffCD);
  //overlap two histo
  overlapHisto(hDiffC, hDiffD, true, histDir, histName, xTitle, xmin, xmax);
  return sfAndErr;
}

vector<double> MyExLepStackHisto::getTransFactDY2(TFile* fData, TFile* fNonDYBkg, TString chDir, TString baseDir, TString histDir, TString histName, TString xTitle, double xmin, double xmax){
  TH1F* hMC_RegC = getHisto(fNonDYBkg,   chDir, baseDir, histDir, "C_"+histName);
  TH1F* hData_RegC= (TH1F*) getHisto(fData, chDir, baseDir, histDir, "C_"+histName);
  //RegionD = LowMET, NonIso
  TH1F* hMC_RegD = getHisto(fNonDYBkg,   chDir, baseDir, histDir, "D_"+histName);
  TH1F* hData_RegD=  getHisto(fData, chDir, baseDir, histDir, "D_"+histName);
  TH1F* hDiffC = (TH1F*)hData_RegC->Clone("hDiffC");
  hDiffC->Add(hMC_RegC, -1);
  TH1F* hDiffD = (TH1F*)hData_RegD->Clone("hDiffD");
  hDiffD->Add(hMC_RegD, -1);
  //If binContent < 0, set it to 0
  makeHistoPositive(hDiffC, false);
  makeHistoPositive(hDiffD, false);
  double intDiffC   = hDiffC->Integral();
  double errDiffC   = getStatUnc(hDiffC, 0.0);
  double intDiffD   = hDiffD->Integral();
  double errDiffD   = getStatUnc(hDiffD, 0.0);
  //Ratio of (Data-MC) from RegionD and RegionC
  double ratioDiffCD = intDiffC/intDiffD;
  double tmpC = errDiffC/intDiffC;
  double tmpD = errDiffD/intDiffD;
  double errDiffCD = ratioDiffCD*sqrt(tmpD*tmpD + tmpC*tmpC);
  cout<<"\n--------------------------"<<endl;
  printf("intC =  %f, intD = %f", intDiffC, intDiffD);
  vector<double>sfAndErr;
  sfAndErr.push_back(ratioDiffCD);
  sfAndErr.push_back(errDiffCD);
  //overlap two histo
  overlapHisto(hDiffC, hDiffD, true, histDir, histName, xTitle, xmin, xmax);
  return sfAndErr;
}

TH1F* MyExLepStackHisto:: getDataDrivenDY(TString chDir, TString baseDir, TString histDir, TString histName, double transFact=1.0, double errorTransFact = 0.0){
  TH1F* hVV =   getHisto(fVV,   chDir, baseDir, histDir, "B_"+histName);
  TH1F* hWJ =   getHisto(fWJ,   chDir, baseDir, histDir, "B_"+histName);
  TH1F* hTT =   getHisto(fTT,   chDir, baseDir, histDir, "B_"+histName);
  TH1F* hData;
  if(chDir.Contains("Muon"))
      hData = getHisto(fMuData, chDir, baseDir, histDir, "B_"+histName);
  if(chDir.Contains("Electron"))
      hData = getHisto(fEleData, chDir, baseDir, histDir, "B_"+histName);
  TH1F* hOtherMC = (TH1F*)hVV->Clone("hOtherMC");
  hOtherMC->Add(hWJ);
  hOtherMC->Add(hTT);
  TH1F* hDataDrivenDY = (TH1F*)hData->Clone("A_"+histName);
  hDataDrivenDY->Add(hOtherMC, -1);
  //If binContent < 0, set it to 0
  makeHistoPositive(hDataDrivenDY, false);
  cout<<histDir<<"B_"+histName<<endl;
  double sError = 0.0;
  double  norm = hDataDrivenDY->IntegralAndError(1, hDataDrivenDY->GetNbinsX(), sError);
  cout<<"intB = "<<norm<<", intB_err = "<<sError<<endl;
  cout<<"qcdSF = "<<transFact<<", qcdSF_err = "<<errorTransFact<<endl;
  double tot_bin_cont = 0.0;
  double tot_bin_err = 0.0;
  for(int ibin=1; ibin<hDataDrivenDY->GetNbinsX(); ibin++){
      double bin_cont = hDataDrivenDY->GetBinContent(ibin);
      double bin_err = hDataDrivenDY->GetBinError(ibin);
      double new_bin_cont = transFact*bin_cont;
      double new_bin_err = sqrt(pow(bin_cont*errorTransFact, 2) + pow(bin_err* transFact, 2));
      tot_bin_cont = tot_bin_cont + new_bin_cont;
      tot_bin_err = tot_bin_err + new_bin_err*new_bin_err;
      hDataDrivenDY->SetBinContent(ibin, new_bin_cont);
      hDataDrivenDY->SetBinError(ibin, new_bin_err);
    }
  cout<<"tot_bin_cont= "<<tot_bin_cont<<", tot_bin_err = "<<sqrt(tot_bin_err)<<endl;
  return hDataDrivenDY;
}

