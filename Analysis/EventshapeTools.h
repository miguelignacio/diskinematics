#ifndef __EVENTSHAPETOOLS_H
#define __EVENTSHAPETOOLS_H

#include <TLorentzVector.h>
#include "H1PhysUtils/H1BoostedJets.h"

class EventshapeTools : public TObject {

public:

  EventshapeTools();
  ~EventshapeTools();

  TLorentzVector CalcScatElec(double q2, double y, double x, double phi, double Ep);
  H1Boost BoostToBreitFrame(double q2, double y, double x, double phi);
  H1Boost BoostToLabFrame(H1Boost boost);

/*
  H1ArrayI MatchModsJets(Float_t R0);
  H1ArrayD UnmatchedRecJetsPt( H1ArrayI MatchedGenJets );
  H1ArrayB MatchedRecJetsArray( H1ArrayI MatchedGenJets );

  // Temporarily   
  static void SetJetMatchingRadius(Float_t r){fJetMatchingRadius = r;}  
  static Float_t GetJetMatchingRadius(){return fJetMatchingRadius;}  
   
  TLorentzVector CalcElecForBoost(Float_t q2, Float_t y, Float_t phi);  
  Float_t GenElecPhotDist();
  Float_t GetInvMassFromJetParts(H1PartJet* jet);
*/
protected:
  ClassDef(EventshapeTools,0)
//private:  

};

#endif
