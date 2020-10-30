#include "JetReweight.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"

OneReweight::OneReweight(string InterpolFile, double* genvar1, double* genvar2) : 
   fFun2(NULL),
   fg2(NULL),
   fVar1(genvar1), 
   fVar2(genvar2),
   fExtrapol(1),
   fLogx(false),
   fLogy(false)
{
   fInterpolation.ReadData(InterpolFile.c_str());
}
OneReweight::OneReweight(TF2* fun2, double* genvar1, double* genvar2) : 
   fFun2(fun2),
   fg2(NULL),
   fVar1(genvar1), 
   fVar2(genvar2),
   fExtrapol(1),
   fLogx(false),
   fLogy(false)
{
   
}

OneReweight::OneReweight(TGraph2D* g2, double* genvar1, double* genvar2) : 
   fFun2(NULL),
   fg2(g2),
   fVar1(genvar1), 
   fVar2(genvar2),
   fExtrapol(1),
   fLogx(false),
   fLogy(false)
{
   
}



double OneReweight::GetRW(){
   // get reweight from pointers
   return GetRW(*fVar1,*fVar2);
}

double OneReweight::GetRW(double x, double y){
   // get reweight 
   if ( x==0 || y==0 )  return 1;
   else if ( fFun2 != NULL ) return fFun2->Eval(x,y);
   else if ( fg2 != NULL )   {
      if ( x<fg2->GetXmin() ) x=fg2->GetXmin();
      if ( x>fg2->GetXmax() ) x=fg2->GetXmax();
      if ( y<fg2->GetYmin() ) y=fg2->GetYmin();
      if ( y>fg2->GetYmax() ) y=fg2->GetYmax();
      return fg2->Interpolate(x,y);
   }
   else return fInterpolation.Interpolate(x,y,fLogx,fLogy,fExtrapol);
}


// _________________________________________________________________________________________ //
//                                   JetReweight                                             //
// _________________________________________________________________________________________ //

OneReweight* JetReweight::AddReweight(TString InterpolFile, Double_t* genvar1, Double_t* genvar2) {
   
   // open files
   string file = InterpolFile.Data();
   fReweights.push_back(OneReweight(file,genvar1,genvar2));
   return &(fReweights.back());
   
}

OneReweight* JetReweight::AddReweight(TF2* fun2, Double_t* genvar1, Double_t* genvar2) {
   // open files
   fReweights.push_back(OneReweight(fun2,genvar1,genvar2));
   return &(fReweights.back());
   
}

OneReweight* JetReweight::AddReweight(TGraph2D* g2, Double_t* genvar1, Double_t* genvar2) {
   // open files
   fReweights.push_back(OneReweight(g2,genvar1,genvar2));
   return &(fReweights.back());
   
}

double JetReweight::GetWeight(){
   // get weight!
   double wgt = 1;
   for ( unsigned int iw = 0 ; iw<fReweights.size() ; iw++ ){
      double w = fReweights[iw].GetRW();
      wgt*=w;
      //cout<<"w="<<w<<"\twgt="<<wgt<<endl;
   }
   return wgt;
}

