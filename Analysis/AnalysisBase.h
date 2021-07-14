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
#include<vector>
class  FidVolCut;
class TTree;
using namespace std;

// ________________________________________________________________ //
//!
//!  TreeVariables
//!
//!  Helper class to keep the variables to be stored into the mini-tree
//!
//!  
class TreeVariables {
public:
   float event_weight =1 ; //weight

   float event_y_e = 0;        //  y            
   float event_Q2_e = 0;       //  Q2           
   float event_y_sigma = 0;        //  y                                                                                                                                                                
      
   float event_Q2_sigma = 0;       //  Q2

   float event_y_da = 0;        // 
   float event_Q2_da = 0;        // 
   float event_y_esigma = 0;   //
   float event_Q2_esigma = 0;   //

   float event_Q2_h = 0 ; //
   float event_y_h  = 0; //
         
   float gen_event_x =0;     // gen x
   float gen_event_y =0;     // gen y
   float gen_event_Q2 =0;    // gen Q2

   float vertex_z = 0;       //  vertex_z     
   float ptmiss = 0;         //  ptmiss   
   float pth   = 0;           // pth    
   float ptratio_da = 0;        //  ptratio with respect to double angle
   float ptratio_ele = 0;      // ptratio with respect to electron      
   float acoplanarity = 0;   //  acoplanarity 
   float Empz = 0;           //  Empz         
   float e_px = 0;           //  e_pt         
   float e_py = 0;          //  e_phi        
   float e_pz = 0;          //  e_rap        
   float HFS_E =0;
   float HFS_px = 0;         //
   float HFS_py = 0;         //
   float HFS_pz = 0;         //
   float gene_px = 0;          //  e_eta        
   float gene_py = 0;            //  e_p          
   float gene_pz = 0;        //  e_theta      
  
};


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
   void SetSysShift(int sys = -9999);
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
   void FillMiniTree();  //!< fill mini tree for azimuthal correlation analysis
   void WriteMiniTree(); //!< wrie mini tree for azimuthal correlation analysis
   //void PrintEventInfo();

   void DoWriteHistograms();

   const TString& GetChainName() const { return fChainName;} //!< get chain name
   
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
   int fSys = -9999; 
   // --- mini tree
   TTree* fMiniTree = NULL;
   TreeVariables fTreeVar;

};


#endif
