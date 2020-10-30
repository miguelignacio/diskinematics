#ifndef H_SPACLINEARITY
#define H_SPACLINEARITY

#include <TVector3.h>
#include "H1Steering/H1Steer.h"

class SpacLinearity : public H1Steer {
 public:
   SpacLinearity(void);
   void SetAlpha(const char *in);
   void SetSuperX(const char *in);
   void SetSuperY(const char *in);
   void SetSuperXY(int i,double dx, double dy);
   void SetAlphaByIndex(int i, Float_t a);
   inline Float_t GetAlpha(int k) const { return fAlpha[k]; }
   inline Float_t GetSuperX(int i) const {return fSuperX[i]; }
   inline Float_t GetSuperY(int i) const {return fSuperY[i]; }
   TVector3 GetPosCorr(TVector3 const &pos,
                       Float_t const *grid,Float_t const *offset) const;
   void PrintSteer(std::ostream &out,bool forcePrint=false) const;
 protected:
   static void PrintN(std::ostream &out,char const *name,Float_t const *f,int n,
                      int npl=0);
   mutable bool alreadyPrinted;
   Float_t fAlpha[5];
   Float_t fSuperX[81];
   Float_t fSuperY[81];
   ClassDef(SpacLinearity,0)
};

#endif
