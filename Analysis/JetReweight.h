#ifndef H_JETREWEIGHT
#define H_JETREWEIGHT

#include <map>
#include <string>
#include "TH1D.h"
#include "TH2D.h"
#include "TH1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "Interpolate2d.h"

using namespace std;
 
class OneReweight {
public:
   OneReweight(string InterpolFile, double* genvar1, double* genvar2) ;
   OneReweight(TF2* fun2, double* genvar1, double* genvar2) ;
   OneReweight(TGraph2D* g2, double* genvar1, double* genvar2) ;
   ~OneReweight() {};
   double GetRW();
   double GetRW(double x, double y=0);
   void SetLogx(bool logx=true) { fLogx=logx;}
   void SetLogy(bool logy=true) { fLogy=logy;}
   void SetExtrapolation(double ex) {fExtrapol=ex;};

protected:
   Interpolate2d fInterpolation;
   TF2* fFun2;
   TGraph2D* fg2;
   double* fVar1;
   double* fVar2;
   double fExtrapol;
   bool fLogx;
   bool fLogy;
};

class JetReweight  {
 public:
   JetReweight() {};
   ~JetReweight() {};
   OneReweight* AddReweight(TString InterpolFile, Double_t* genvar1, Double_t* genvar2=NULL);
   OneReweight* AddReweight(TF2* fun2, Double_t* genvar1, Double_t* genvar2=NULL);
   OneReweight* AddReweight(TGraph2D* g2, Double_t* genvar1, Double_t* genvar2=NULL);
   double GetWeight();
   vector<OneReweight> fReweights;
};

#endif
