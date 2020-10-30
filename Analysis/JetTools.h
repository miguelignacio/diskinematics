#ifndef __JETTOOLS_H
#define __JETTOOLS_H

#include "H1Mods/H1PartJet.h"
#include <TLorentzVector.h>
#include "H1Arrays/H1ArrayI.h"
#include "H1Arrays/H1ArrayB.h"
#include "H1Arrays/H1ArrayD.h"

class JetTools
{

 public:

  JetTools();
  ~JetTools();

  H1ArrayI MatchModsJets(Float_t R0);
  H1ArrayD UnmatchedRecJetsPt( H1ArrayI MatchedGenJets );
  H1ArrayB MatchedRecJetsArray( H1ArrayI MatchedGenJets );
  //H1ArrayI MatchJetsUnfold(Float_t R0 = 0.6 , Int_t DontMatchGenI = -1 , Int_t DontMatchRecI = -1 );
  H1ArrayI MatchJetsUnfold( Bool_t* GenJetsPSArr, Bool_t* RecJetsPSArr , Int_t DontMatchGenI = -1 , Int_t DontMatchRecI = -1);
  H1ArrayI MatchJetsUnfold2( Bool_t* GenJetsPSArr, Bool_t* RecJetsPSArr );
  static Double_t Dist(TLorentzVector* jet1, TLorentzVector* jet2);  

  // Temporarily   
  static void SetJetMatchingRadius(Float_t r){fJetMatchingRadius = r;}  
  static Float_t GetJetMatchingRadius(){return fJetMatchingRadius;}  
   
  TLorentzVector CalcElecForBoost(Float_t q2, Float_t y, Float_t phi);  
  Float_t GenElecPhotDist();
  Float_t GetInvMassFromJetParts(H1PartJet* jet);
   
private:  
  static Float_t fJetMatchingRadius;

};

#endif
