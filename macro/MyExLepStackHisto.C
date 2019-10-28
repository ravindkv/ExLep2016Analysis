#include "MyExLepStackHisto.h"

void plotStackedHisto(TString chDir, TString baseDir, TString histDir, TString histName, TString xTitle,bool isData=false, bool isSig=false, double xmin=0, double xmax=10, double unc=false, bool isDYdd=false, bool isMerge=false){
  MyExLepStackHisto MyELSHisto;
  string hist_name (histName);
  //Pad
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(3);
  TCanvas *c1 = new TCanvas(histName, histName, 700, 700);
  //TCanvas *c1 = new TCanvas("c1", "Data_MC", 400, 600);
  const float xpad[2] = {0.,1};
  const float ypad[4] = {0.,0.30,0.30,0.98};
  if(isData){
    //c1->Divide(1, 2); c1->cd(1);
    c1->Divide(1, 2, 0, 0); c1->cd(1);
    gPad->SetPad(xpad[0],ypad[2],xpad[1],ypad[3]);
    if(isData) gPad->SetLogy(true);
    if(hist_name.find("mjj") != string::npos) gPad->SetLogy(false);
  }
  //-------------------------------
  //merge tailing bins
  double mergeThres = 0.01; // 1%
  //-------------------------------
  TH1F* hBkg = MyELSHisto.getHisto(fBkg, chDir, baseDir, histDir, histName);
  int mergeFrom = 0.0;
  double integralData = hBkg->Integral();
  int totBins = int(hBkg->GetNbinsX());
  for(int bin=1; bin<totBins+1; bin++){ 
    double sError = 0.0;
    double partIntegral = hBkg->IntegralAndError(bin, totBins, sError);
    if(double(partIntegral/integralData) <= mergeThres){
      mergeFrom = bin;
      break;
    }
  }
  double binWidth = hBkg->GetBinWidth(5);
  double nbins = (xmax-xmin)/binWidth;
  cout<<"=================="<<endl;
  cout<<"mergeFrom = "<<mergeFrom<<endl;
  cout<<"binWidth = "<<binWidth<<endl;
  cout<<"nbins = "<<nbins<<endl;

  bool isMuChannel = false;
  bool isEleChannel = false;
  if(chDir.Contains("Muon")) isMuChannel = true;
  if(chDir.Contains("Electron")) isEleChannel = true;
  //-------------------------------
  // stack MC Bkg histo
  //-------------------------------
  THStack* hStack = new THStack("hStack","");
  TH1F* hST = MyELSHisto.getHisto(fST, chDir, baseDir, histDir, histName, isMerge, mergeFrom, nbins);
  TH1F* hMC = (TH1F*)hST->Clone("hMC");
  int col_depth =0;
  hST->SetFillColor(kBlue + col_depth);
  hStack->Add(hST);
  TH1F* hWJ = MyELSHisto.stackHisto(fWJ, chDir, baseDir, histDir, histName, kOrange+col_depth, 1,   hStack, hMC, isMerge, mergeFrom, nbins);
  TH1F* hVV = MyELSHisto.stackHisto(fVV, chDir, baseDir, histDir, histName, kViolet +col_depth , 1,   hStack, hMC, isMerge, mergeFrom, nbins);
  TH1F* hTT = MyELSHisto.stackHisto(fTT, chDir, baseDir, histDir, histName, kCyan+col_depth, 1,   hStack, hMC, isMerge, mergeFrom, nbins);

  // trim the histDir string
  std::string histDir_str;
  std::string histDir_orig(histDir);
  std::remove_copy(histDir_orig.begin(), histDir_orig.end(), std::back_inserter(histDir_str), '/');
  TString histDir_(histDir_str);
  //
  //-------------------------------
  // DY from Data
  //-------------------------------
  //qcd scale factors for data-driven DY
  double dySF = 1.0;
  double dyErr = 0.0;
  if(isDYdd){
    vector<double> sfAndErr;
    if(isMuChannel) sfAndErr = MyELSHisto.getTransFactDY2(fMuData, fNonDYBkg, chDir, baseDir, histDir, histName, xTitle, xmin, xmax);
    if(isEleChannel) sfAndErr = MyELSHisto.getTransFactDY2(fEleData, fNonDYBkg, chDir, baseDir, histDir, histName, xTitle, xmin, xmax);
    //vector<double> sfAndErr = MyELSHisto.getTransFactDY(fData, fTT, fWJ, fVV, chDir, baseDir, histDir, histName, xTitle, xmin, xmax);
    dySF = sfAndErr[0];
    dyErr = sfAndErr[1];
  }
  TH1F * hDYdd = MyELSHisto.getHisto(fDY, chDir, baseDir, histDir, histName);
  hDYdd->Reset(); // initialize empty hist
  if(isDYdd){
    hDYdd = MyELSHisto.getDataDrivenDY(chDir, baseDir, histDir, histName,  dySF,  dyErr);
    hDYdd->SetFillColor(kGreen+col_depth);
    hDYdd->GetXaxis()->SetRangeUser(xmin,xmax);
    //create same dir to the data driven qcd file
    std::string histPath = std::string(chDir+baseDir+histDir_);
    TDirectory *d = fDYdd->GetDirectory(histPath.c_str());
    if(!d) fDYdd->mkdir(histPath.c_str());
    fDYdd->cd(histPath.c_str());
    //hDY->Draw();
    hDYdd->Write();
    hStack->Add(hDYdd);
    hMC->Add(hDYdd);
  }
  else hDYdd = MyELSHisto.stackHisto(fDY, chDir, baseDir, histDir, histName, kGreen +col_depth, 1,   hStack, hMC, isMerge, mergeFrom, nbins); 

  if(isData) c1->cd(1);
  else c1->cd();
  gPad->SetTopMargin(0.10);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  hStack->Draw("HIST");
  hStack->SetMinimum(1.0);
  hStack->SetMinimum(0.01);
  //if(histDir.Contains("ZTag")) hStack->SetMinimum(0.1);
  hStack->GetXaxis()->SetRangeUser(xmin, xmax);
  //cout<<hStack->GetMaximum()<<endl;
  if(isData){
    hStack->GetYaxis()->SetTitleOffset(0.70);
    hStack->GetYaxis()->SetTitleSize(0.10);   
    hStack->GetYaxis()->SetLabelSize(0.07);   
    hStack->GetYaxis()->SetTickLength(0.04); 
    hStack->GetYaxis()->SetTitle("Events");
    hStack->GetXaxis()->SetTitleOffset(1.20);
  }
  else{
  hStack->GetYaxis()->SetTitle("Events");
  hStack->GetXaxis()->SetTitle(xTitle);
  hStack->GetXaxis()->SetTitleSize(0.07);
  hStack->GetXaxis()->SetLabelSize(0.05);   
  hStack->GetXaxis()->SetTickLength(0.05); 
  //hStack->GetYaxis()->SetNdivisions(5);
  hStack->GetYaxis()->SetTitleSize(0.08);   
  hStack->GetYaxis()->SetLabelSize(0.05);   
  hStack->GetYaxis()->SetTickLength(0.04); 
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.15);
  hStack->GetYaxis()->SetTitleOffset(0.80);
  hStack->GetXaxis()->SetTitleOffset(0.90);
  }

  //-------------------------------------///
  //unc band
  //-------------------------------------///
  TGraphAsymmErrors *UncBand;
  if(unc){
  UncBand = MyELSHisto.UNCGRAPH(
            MyELSHisto.addHistoForUnc("base/", 	 baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JESPlus/",      baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JESMinus/",     baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JERPlus/",      baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JERMinus/",     baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("bTagPlus/",     baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("bTagMinus/",    baseDir, histDir, histName),
	    true, false);
  UncBand->SetFillColor(17);
  UncBand->SetFillStyle(3008);
  UncBand->Draw(" E2 same");
  }

  TH1F * hData;
  if(isMuChannel) hData = MyELSHisto.getHisto(fMuData, chDir, baseDir, histDir, histName, isMerge, mergeFrom, nbins);
  if(isEleChannel) hData = MyELSHisto.getHisto(fEleData, chDir, baseDir, histDir, histName, isMerge, mergeFrom, nbins);
  ///MyELSHisto.decorateHisto(hData, "", xTitle, "Events");
  hData->SetFillColor(kBlack);
  hData->SetMarkerStyle(20); 
  hData->SetMarkerSize(1.2);
  if(isData)hData->Draw("Esame"); 

  //-------------------------------
  //Signal 
  //-------------------------------
  TH1F* hSig250;
  if(isMuChannel) hSig250 = MyELSHisto.getHisto(fSigMuMuZ250, chDir, baseDir, histDir, histName);
  if(isEleChannel) hSig250 = MyELSHisto.getHisto(fSigEEZ250, chDir, baseDir, histDir, histName);
  ///MyELSHisto.decorateHisto(hSig250, "", xTitle, "Events");
  hSig250->SetLineColor(kRed); hSig250->SetLineStyle(2);
  hSig250->SetLineWidth(3); hSig250->SetFillColor(0);
  if(isSig)hSig250->Draw("HISTSAME"); 

  TH1F* hSig1000;
  if(isMuChannel) hSig1000 = MyELSHisto.getHisto(fSigMuMuZ1000, chDir, baseDir, histDir, histName);
  if(isEleChannel) hSig1000 = MyELSHisto.getHisto(fSigEEZ1000, chDir, baseDir, histDir, histName);
  ///MyELSHisto.decorateHisto(hSig1000, "", xTitle, "Events");
  hSig1000->SetLineColor(kGreen); hSig1000->SetLineStyle(2);
  hSig1000->SetLineWidth(3); hSig1000->SetFillColor(0);
  if(isSig)hSig1000->Draw("HISTSAME"); 

  TH1F* hSig2000;
  if(isMuChannel) hSig2000 = MyELSHisto.getHisto(fSigMuMuZ2000, chDir, baseDir, histDir, histName);
  if(isEleChannel) hSig2000 = MyELSHisto.getHisto(fSigEEZ2000, chDir, baseDir, histDir, histName);
  ///MyELSHisto.decorateHisto(hSig2000, "", xTitle, "Events");
  hSig2000->SetLineColor(kOrange); hSig2000->SetLineStyle(2);
  hSig2000->SetLineWidth(3); hSig2000->SetFillColor(0);
  if(isSig)hSig2000->Draw("HISTSAME"); 

  TH1F* hSig5000;
  if(isMuChannel) hSig5000 = MyELSHisto.getHisto(fSigMuMuZ5000, chDir, baseDir, histDir, histName);
  if(isEleChannel) hSig5000 = MyELSHisto.getHisto(fSigEEZ5000, chDir, baseDir, histDir, histName);
  ///MyELSHisto.decorateHisto(hSig5000, "", xTitle, "Events");
  hSig5000->SetLineColor(kOrange+2); hSig5000->SetLineStyle(2);
  hSig5000->SetLineWidth(3); hSig5000->SetFillColor(0);
  if(isSig)hSig5000->Draw("HISTSAME"); 
  //-------------------------------
  //Legends
  //-------------------------------
  TLegend* leg = new TLegend(0.7218792,0.4061504,0.9212081,0.8798861,NULL,"brNDC");
  //TLegend* leg = new TLegend(0.7618792,0.3061504,0.9712081,0.8798861,NULL,"brNDC");
  if(hist_name.find("pt") != string::npos || hist_name.find("mt") != string::npos || hist_name.find("Fit") != string::npos ||hist_name.find("RelIso") != string::npos){
    leg = new TLegend(0.5018792,0.6061504,0.9512081,0.8898861,NULL,"brNDC");
    leg->SetNColumns(2);
  }
  if(hist_name.find("mlZ") != string::npos){
    leg = new TLegend(0.6518792,0.3061504,0.8512081,0.8798861,NULL,"brNDC");
  }
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(kBlack);
  leg->SetTextFont(42);
  leg->SetTextAngle(0);
  leg->SetTextSize(0.04);
  leg->SetTextAlign(12);
  //leg->AddEntry(hSig250, "#splitline{Signal}{M_{H^{+}} = 120 GeV}","L");
  if(isData)leg->AddEntry(hData,"Data","PE");
  leg->AddEntry(hDYdd,"DY","F");
  leg->AddEntry(hTT,"t#bar{t} + jets","F");
  leg->AddEntry(hVV,"VV","F");
  leg->AddEntry(hWJ,"W + jets","F");
  leg->AddEntry(hST,"Single t","F");
  //if(isSig)leg->AddEntry(hSig250, "Signal","L");
  if(unc)leg->AddEntry(UncBand, "Uncertainty","F");
  if(isSig)leg->AddEntry(hSig250,  "M_{l^{*}} = 250 GeV","L");
  if(isSig)leg->AddEntry(hSig1000, "M_{l^{*}} = 1 TeV","L");
  if(isSig)leg->AddEntry(hSig2000, "M_{l^{*}} = 2 TeV","L");
  if(isSig)leg->AddEntry(hSig5000, "M_{l^{*}} = 5 TeV","L");
  //if(isSig)leg->AddEntry(hSig250, "Signal","L");
  //if(isSig)leg->AddEntry((TObject*)0, "(M_{l^{*}} = 250 GeV)","");
  leg->Draw();

  double yMax = 0;
  if(hData->GetMaximum() > hSig250->GetMaximum()) yMax = hData->GetMaximum();
  else yMax = hSig250->GetMaximum();
  if(yMax < hMC->GetMaximum()) yMax = hMC->GetMaximum();

  if(isData) hStack->SetMaximum(4.0*hStack->GetMaximum());
  else hStack->SetMaximum(1.1*yMax);
  TPaveText *cct = MyELSHisto.paveText(0.30,0.8454,0.40,0.8462, 0, 19, 1, 0, 132);
  cct->SetTextSize(0.055);
  cct->AddText("CMS, Preliminary");

  //-------------------------------------///
  //  Draw Pave Text 
  //-------------------------------------///
  //hist name
  TPaveText *hLable = MyELSHisto.paveText(0.6513423,0.7754898,0.6010067,0.8962187, 0, 19, 1, 0, 132);
  hLable->SetTextSize(0.07);
  hLable->AddText(xTitle);
  
  //channel
  TPaveText *ch = MyELSHisto.paveText(0.823,0.9154898,0.9210067,0.9762187, 0, 19, 1, 0, 132);
  ch->SetTextSize(0.10);
  if(isMuChannel) ch->AddText("#mu + jets");
  if(isEleChannel) ch->AddText("e + jets");
  //CMS prili
  TPaveText *pt = MyELSHisto.paveText(0.01,0.9554,0.82,0.9562, 0, 19, 1, 0, 132);
  if(isData) pt->SetTextSize(0.080);
  else pt->SetTextSize(0.05);
  pt->AddText(histDir_+": 35.9 fb^{-1} (13 TeV)");
  //TText *text = pt->AddText(histDir+": CMS Preliminary, #sqrt{s} = 13 TeV, 35.9 fb^{-1}");
  pt->Draw();
  if(isSig) cct->Draw();
  ch->Draw();
  //hLable->Draw();
  gPad->RedrawAxis();
  c1->Update();
  
  //-------------------------------------///
  // Ratio = DATA/Bkg
  //-------------------------------------///
  if(isData){
    c1->cd(2);
    gPad->SetTopMargin(0); gPad->SetBottomMargin(0.5); //gPad->SetGridy();
    if(histDir=="") gPad->SetBottomMargin(0.55);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.05);
    gPad->SetPad(xpad[0],ypad[0],xpad[1],ypad[2]);
    TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
    hRatio->Reset();
    hRatio->Add(hData);
    hRatio->Divide(hMC); 
    MyELSHisto.decorateHisto(hRatio, "", xTitle, "#frac{Data}{Bkg}");
    hRatio->SetFillColor(kBlack);
    if(histDir.Contains("ZTag")) hRatio->GetYaxis()->SetRangeUser(0, 2);
    else hRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
    //else hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetXaxis()->SetRangeUser(xmin, xmax);
    hRatio->GetYaxis()->SetTitleOffset(0.40);
    hRatio->GetXaxis()->SetTitleOffset(0.90);
    hRatio->SetMarkerStyle(20); hRatio->SetMarkerSize(1.2);
    hRatio->GetYaxis()->SetTitleSize(0.15); 
    hRatio->GetXaxis()->SetTitleSize(0.15);
    hRatio->GetXaxis()->SetLabelSize(0.10); 
    hRatio->GetYaxis()->SetLabelSize(0.10); 
    if(hist_name.find("mjj") != string::npos){
      hRatio->GetXaxis()->SetTitleSize(0.05); 
      hRatio->GetXaxis()->SetTitleOffset(1.40);
    }
    //lable x-axis, for cutflow
    if(histName=="cutflow"){
      vector<string >cut_label;
      if(isEleChannel){
        cut_label.push_back("Ele trigger");
      }
      if(isMuChannel){
        cut_label.push_back("Mu trigger");
      }
      cut_label.push_back("Control Sel");
      cut_label.push_back("b-jet veto");
      cut_label.push_back("PreSel");
      cut_label.push_back("ZTag Sel");
      for(int istep=0; istep<cut_label.size(); istep++ ){
       hRatio->GetXaxis()->SetBinLabel(istep+1, cut_label[istep].c_str());
      }
      hRatio->GetXaxis()->LabelsOption("v");
      hRatio->GetXaxis()->SetTickLength(0.08); 
      hRatio->GetXaxis()->SetLabelOffset(0.03);
      hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
    }
    //unc band
    hRatio->Draw("E"); // use "P" or "AP"
    if(unc){
    TGraphAsymmErrors *UncBand_Ratio;
    UncBand_Ratio = MyELSHisto.UNCGRAPH(
	    MyELSHisto.addHistoForUnc("base/", 	baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JESPlus/",     baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JESMinus/",    baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JERPlus/",     baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("JERMinus/",    baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("bTagPlus/",    baseDir, histDir, histName),
      	    MyELSHisto.addHistoForUnc("bTagMinus/",   baseDir, histDir, histName),
	    false, true);
    UncBand_Ratio->SetFillColor(17);
    //UncBand_Ratio->SetFillStyle(19);
    UncBand_Ratio->Draw("E2 same");
    }
    hRatio->Draw("E same"); // use "P" or "AP"
    //base line at 1
    TF1 *baseLine = new TF1("baseLine","1", -100, 2000); 
    baseLine->SetLineColor(kBlack);
    baseLine->Draw("SAME");
    c1->Update();
  }
  if(isSaveHisto){
    mkdir("outputFiles", S_IRWXU);
    mkdir("outputFiles/"+chDir, S_IRWXU);
    mkdir("outputFiles/"+chDir+histDir, S_IRWXU);
    //mkdir(histDir_, S_IRWXU);
    TString outFile("$PWD/");
    outFile += "outputFiles/"+chDir+histDir+histName;
    if(isMuChannel) outFile += "_mu"+histDir_+".pdf";
    if(isEleChannel) outFile += "_ele"+histDir_+".pdf";
    c1->SaveAs(outFile);
    //c1->Close();
  }
}
void stackAllHisto(TString chDir, TString histDir, bool isDYdd){
  TString baseDir = "base/";
  bool isDataMjj= false;
  bool isData = true;
  //flags
  bool isSig = true;
  bool isUnc = false;
  bool isDDtmp = false;
  //---------------------------------------
  bool isMuChannel = false;
  bool isEleChannel = false;
  if(chDir.Contains("Muon")) isMuChannel = true;
  if(chDir.Contains("Electron")) isEleChannel = true;
  //---------------------------------------
  //plotStackedHisto(chDir, baseDir, "", "cutflow","cutflow", isData,  isSig,  0.5, 7.5, false);
  if(!isDYdd){
    plotStackedHisto(chDir, baseDir, histDir, "pt_jet", "Pt^{jets} [GeV]", isData, isSig,  0, 1400.0, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "MET", "MET [GeV]", isData, isSig,  0, 1000.0, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "eta_jet", "#eta^{jets}", isData, isSig,  -3.0, 5.0, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "final_multi_jet", "N^{jets}", isData, isSig,  0, 7 , isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "nvtx", "N^{vertex}", isData, isSig,  0, 70, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "rhoAll", "#rho", isData, isSig,  0, 70, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "ak8Pmass", "Pruned mass of jets [GeV]", isData, isSig,  0.0, 600, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "ak8Tau21", "Jet-subjetiness (#tau_{21})", isData, isSig,  0, 1.5, isUnc);
    TString ch =""; 
    if(isMuChannel)  ch = "#mu";
    if(isEleChannel) ch = "e";
    plotStackedHisto(chDir, baseDir, histDir, "dR1", "#Delta R (jets, "+ch+"_{1})", isData, isSig,  0, 8, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "dR2", "#Delta R (jets, "+ch+"_{2})", isData, isSig,  0, 8, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "pt1stLep", "Pt^{"+ch+"_{1}} [GeV]", isData, isSig,  0.0, 1000.0, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "eta1stLep", "#eta^{"+ch+"_{1}}", isData, isSig,  -3.0, 5.0, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "pt2ndLep", "Pt^{"+ch+"_{2}} [GeV]", isData, isSig,  0.0, 1000.0, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "eta2ndLep", "#eta^{"+ch+"_{2}}", isData, isSig,  -3.0, 5.0, isUnc);
    //plotStackedHisto(chDir, baseDir, histDir, "multi_mu", "N^{"+ch+"}", isData, isSig,  1, 5, isUnc);
    plotStackedHisto(chDir, baseDir, histDir, "mll", "M^{"+ch+"_{1}"+ch+"_{2}} [GeV]", isData, isSig,  0, 1400, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "pt_Z", "Pt^{"+ch+"_{1}"+ch+"_{2}}[GeV]", isData, isSig,  0, 1000, isUnc, isDDtmp, true);
    plotStackedHisto(chDir, baseDir, histDir, "eta_Z", "#eta^{"+ch+"_{1}"+ch+"_{2}}", isData, isSig,  -3.0, 4.5, isUnc);
    if(histDir=="ZTag/"){
      plotStackedHisto(chDir, baseDir, histDir, "mlZ_max","M^{"+ch+" Z}_{max}", isDataMjj, isSig,  200, 3000, isUnc, isDDtmp, true);
      plotStackedHisto(chDir, baseDir, histDir, "mlZ_min","M^{"+ch+" Z}_{min}", isDataMjj, isSig,  50, 1500, isUnc, isDDtmp, true);
    }
    plotStackedHisto(chDir, baseDir, "ControlP/", "pfCISV","b-discriminator", isData, isSig,  0.0, 1.5, isUnc, isDDtmp);
    plotStackedHisto(chDir, baseDir, "ControlP/", "multi_bjet","N^{b-jet}", isData, isSig,  -0.5, 3.5, isUnc, isDDtmp);
  }
  else{
   vector<string>dirVec;
   dirVec.push_back("ZTag1");
   dirVec.push_back("ZTag2");
   dirVec.push_back("ZTag3");
   dirVec.push_back("ZTag4");
   dirVec.push_back("ZTag5");
   dirVec.push_back("ZTag6");
   dirVec.push_back("ZTag7");
   dirVec.push_back("ZTag8");
   dirVec.push_back("ZTag9");
   dirVec.push_back("ZTag10");
   vector<string> sigMass;
   sigMass.push_back("250");   
   sigMass.push_back("500");   
   sigMass.push_back("750");   
   sigMass.push_back("1000");  
   sigMass.push_back("1250");  
   sigMass.push_back("1500");  
   sigMass.push_back("1750");  
   sigMass.push_back("2000");  
   sigMass.push_back("2500");  
   sigMass.push_back("3000");  
   sigMass.push_back("3500");  
   sigMass.push_back("4000");  
   sigMass.push_back("4500");  
   sigMass.push_back("5000");  
   for(unsigned int l =0; l<sigMass.size(); l++){
     TString mass = sigMass[l];
     plotStackedHisto(chDir, baseDir, "ZTag3/", "mlZ_max_sig"+mass, "A_mlZ_max_sig_ZTag3_"+mass, isData, isSig, 20, 2000, isUnc, isDYdd);
     /*
     for(unsigned int d = 0; d<dirVec.size(); d++){
       TString dir = dirVec[d];
       plotStackedHisto(chDir, baseDir, dir+"/", "mlZ_max_sig"+mass, "A_mlZ_max_sig_"+dir+"_"+mass, isData, isSig, 20, 2000, isUnc, isDYdd);
     }
     */
   }
  }
}

void MyExLepStackHisto(){
  bool isDYdd = false;
  //stackAllHisto("Muon/", "ControlP/",       isDYdd);
  //stackAllHisto("Muon/", "PreSel/",         isDYdd);
  //stackAllHisto("Muon/", "ZTag/",           isDYdd);
  stackAllHisto("Electron/", "ControlP/",   isDYdd);
  /*
  stackAllHisto("Muon/", "ZTag/",           isDYdd);

  stackAllHisto("Electron/", "ControlP/",   isDYdd);
  stackAllHisto("Electron/", "PreSel/",     isDYdd);
  stackAllHisto("Electron/", "ZTag/",       isDYdd);
  */
}

