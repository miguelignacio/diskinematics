#ifndef __ANALYSISBASE_H
#define __ANALYSISBASE_H

// stl includes

// H1 includes
class  H1DetectorStatus;
class  H1RunList;

// Analysis includes
#include "SpacLinearity.h"
#include "Alignment.h"
#include "H2020HistManager.h"
class  FidVolCut;
class TTree;
using namespace std;


// ________________________________________________________________ //
//!
//!  AnalysisBase
//!
//!  Base class for the event shape analysis
//!
//!  AnalysisBase provides 'interface' for derived analysis class
//!  and also applies *cuts*
//!
//!  Cuts are presently implemented from JetsAtHighQ2 and reproduce
//!  this analysis exactly.
//!
//!  Function starting with 'Do' are called by main
//!  virtual 'Do'-functions must be implemented in inherited class
//!  
class AnalysisBase { //: public TObject {
   
public:

   //! Constructors
   AnalysisBase(TString chain);
   virtual ~AnalysisBase();

   void DoBaseReset();
   virtual void DoReset() = 0;
   void DoBaseInitialSettings();
   virtual void DoInitialSettings() = 0;
   bool DoBasicCutsRec();
   bool DoBasicCutsGen();
   virtual bool DoAnalysisCutsGen() = 0;
   virtual bool DoAnalysisCutsRec() = 0;
   virtual void DoCrossSectionObservablesGen() = 0;
   virtual void DoCrossSectionObservablesRec() = 0;
   virtual void DoControlPlotsGen() = 0;
   virtual void DoControlPlotsRec() = 0;
   virtual void DoControlPlotsGenRec() = 0;
   virtual void DoCrossSectionsGenRec() = 0;
   
   void InitMiniTree(); //!< init mini tree for azimuthal correlation analysis
   void FillMiniTree();  //!< fil mini tree for azimuthal correlation analysis
   void WriteMiniTree(); //!< wrie mini tree for azimuthal correlation analysis
   //void PrintEventInfo();

   void DoWriteHistograms();

protected:   

   double CalcDist( const TLorentzVector& jet1, const TLorentzVector& jet2 );

   // void SetSystematics();

   //! convert a TObjArray into a std::vector, using proper type_cast
   template<class TObj, class TArr>
   std::vector<TObj> to_vector(TArr* array) { 
      std::vector<TObj> ret(array->GetEntries());
      for ( int i = 0 ; i<array->GetEntries() ; i++ ) 
         ret[i] = static_cast<TObj>(array->At(i));
      return ret;
   }

protected:   

   // --- Run parameters
   bool fGenOnly         =  false;  // set externally
   bool fNoRadMC         =  false;  // set externally
   bool IsMC             =  true;   // init
   bool IsBkgMC          =  false;  // set externally
   bool IsHQ             =  true;   // high-Q2 or low-Q2 MC

   // --- extra objects
   H1RunList* fRunList;               //!< run list
   H1DetectorStatus* fDectStatus;     //!< Detector status
   FidVolCut* fFidVolCut;             //!< fiducial volume cut

   // --- chain identification
   TString fChain;
   TString fChainName;
   double  fLumiMC;
   double  fLumiData;

   // --- event/run specific parameters
   bool fBasicCutsRec = false;
   bool fBasicCutsGen = false;
   Double_t f_VtxZMin, f_VtxZMax;

   // --- mini tree
   TTree* fMiniTree = NULL;

};

#endif
