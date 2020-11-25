#ifndef __ANALYSISEVENTSHAPES_H
#define __ANALYSISEVENTSHAPES_H

#include "H1PhysUtils/H1BoostedJets.h"
#include "AnalysisBase.h"
class H1Part;

using namespace std;


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

   // histogram fillers
   void FillBasicNCDISHists(const string& hm, double Q2, double Y, double X, double weight);
   void FillBasicH1CalculatorGen(const string& hm, double weight);
   void FillBasicH1CalculatorRec(const string& hm, double weight);
   void FillTrackPlots(const string& hm, const std::vector<H1Part*>& parts , double weight);
   void FillTrackPlots(const string& hm, TObjArray* parts , double weight);


   //primitive variable einfach so ( bsp weight)
   //objekte als pointer  ( darf geaendert werden) / reference ( const ref: darf nicht editiert werden)
   //   vector<TLorentzVector> BoostParticleArray(const string& hm, const vector<H1PartCand*> partarray, H1Boost& myboost, double Q2, double X, double Y);

   vector<TLorentzVector> BoostParticleArray (const string& hm, const vector<H1PartCand*> partarray, H1Boost& myboost, double Q2, double X, double Y);
   vector<TLorentzVector> CalculateEventShape_tauzQ (const string& hm, const vector<TLorentzVector> partarray, H1Boost& myboost, double Q2);
   void PlotKinematicVaribles( const string& hm,double Q2, double X, double Y);
   void ClassicalEventShapes (const string& hm, const vector<TLorentzVector> BoostedHFS);
   TLorentzVector BoostToBreitFrame ( double Q2, double y, double E0, double Phi );
   H1Boost CalcBoost(double q2, double y, double x, double phi, double Ep);   

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
      double tau_zQ             = 99;     //!<  Definition of tau_zQ from ...
      double tau1b              = 99;     //!<  Definition of tau_1^b from https://arxiv.org/pdf/1303.6952.pdf
      double tau_zP             = 99;     //!<
      double sumpz              = 0;      //!<
      vector<TLorentzVector> breit_current; //!<
   };

   CrossSectionQuantities   fRec;
   CrossSectionQuantities   fGen;

};

#endif
