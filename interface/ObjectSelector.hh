#ifndef _objectselector_h_
#define _objectselector_h_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include <exception>

#ifdef _STANDALONE
#include "Reader.h"
#else
#include "interface/Reader.h"
#endif
#include "interface/UncertaintyComputer.hh"
#endif

class ObjectSelector : public UncertaintyComputer
{
public : 
  ObjectSelector() : UncertaintyComputer()
  {
    setDefaultSelection();
  }
  virtual ~ObjectSelector(){}
  ///~ObjectSelector(){}

  void setDefaultSelection(){
    defaultSelection_=true;
  }

  // preselection of objects
  void preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex);
  void preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isData=false);
  void preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ, int jes=0, int jer=0);

  //Loose Lepton veto
  bool looseElectronVeto(unsigned long selectedElectron, const vector<MyElectron> & vE, MyVertex & vertex);
  bool looseMuonVeto( int selectedMuon, const vector<MyMuon> & vM);
  
  //HighPt Muon ID
  bool isHighPtMuon(const MyMuon * m);
  bool heepElectronID_HEEPV70(const MyElectron *e, MyVertex & vertex);
  // object cleaning
  void ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, double DR );

  void JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE, vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, double DR);
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV); 
    
  
private :
  bool defaultSelection_;
  ClassDef(ObjectSelector, 1)
};
#endif
