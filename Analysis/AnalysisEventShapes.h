#ifndef __ANALYSISEVENTSHAPES_H
#define __ANALYSISEVENTSHAPES_H

#include "H1PhysUtils/H1BoostedJets.h"
#include "AnalysisBase.h"
#include "fastjet/PseudoJet.hh"
class H1Part;

using namespace std;
//using namespace fastjet;

// ________________________________________________________________ //
//!
//!  AnalysisEventShapes : AnalysisBase
//!
//!  Analysis class for the event shape analysis
//!
//!  AnalysisEventShapes fills
//!     + all control histograms
//!     + applies final kinematic cuts
//!       basic cuts are implemented in Base
//!     + fill histograms for final cross sections
//!
//!  Function starting with 'Do' are called by main
//!  virtual 'Do'-functions must be implemented in inherited class
//!
class AnalysisEventShapes : public AnalysisBase { //: public TObject {
   
public:

   //! Constructors
   AnalysisEventShapes(TString chain);
   //AnalysisEventShapes(const TString& chain);
   ~AnalysisEventShapes();

   virtual void DoReset() override;
   virtual void DoInitialSettings() override;
   virtual bool DoAnalysisCutsGen() override;
   virtual bool DoAnalysisCutsRec() override;
   virtual void DoCrossSectionObservablesGen() override;
   virtual void DoCrossSectionObservablesRec() override;
   virtual void DoControlPlotsGen() override;
   virtual void DoControlPlotsRec() override;
   virtual void DoControlPlotsGenRec() override;
   virtual void DoCrossSectionsGenRec() override;

protected:   

   // void SetSystematics();

    void ClassicalEventShapes (const string& hm, const vector<TLorentzVector> BoostedHFS);
   
 

protected:   

   // --- event classification
   bool fAnalysisCutsGen = false;
   bool fAnalysisCutsRec = false;

   struct CrossSectionQuantities {
     //similar to a class, but all variables are public
      double wgt                = 0 ;     //!<  event weight 
      bool IsGood               = false;  //!<  all cuts are fulfilles
      double Q2                 = 0 ;     //!<  Q2 for cross sections
      double Y                  = 0 ;     //!<  y for cross sections 
      double X                  = 0 ;     //!<  x for cross sections
      double tau_zQ             = 0 ;     //!<  Definition of tau_zQ from ...
      double tau1b              = 0 ;     //!<  Definition of tau_1^b from https://arxiv.org/pdf/1303.6952.pdf
      double tau_zP             = 0 ;     //!<
      double sumpz              = 0 ;     //!<
      vector<TLorentzVector> breit_current; //!< 4-vectors of all particles in the current hem. in the Breit frame
      std::vector<fastjet::PseudoJet>     genjets;      
      H1Boost BoostToBreit;                 //!< Lorentz boost to breit frame ( q + 2*x*P = 0 )
      H1Boost BoostToLab;                   //!< Lorentz boost back to the lab frame
   };

   CrossSectionQuantities   fRec;
   CrossSectionQuantities   fGen;

};

#endif
