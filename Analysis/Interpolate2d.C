#include "Interpolate2d.h"
#include <TMatrixD.h>
#include <TMath.h>
#include <fstream>

using namespace std;

void Interpolate2d::ImportHistogram(const TH2 *histogram,bool useLogZ) {
   fLogZ=useLogZ;
   int nx=histogram->GetNbinsX();
   int ny=histogram->GetNbinsY();
   fX.ResizeTo(nx+1);
   fY.ResizeTo(ny+1);
   for(int ix=1;ix<=nx;ix++) {
      fX(ix-1)=histogram->GetXaxis()->GetBinLowEdge(ix);
      fX(ix)=histogram->GetXaxis()->GetBinUpEdge(ix);
   }
   for(int iy=1;iy<=ny;iy++) {
      fY(iy-1)=histogram->GetYaxis()->GetBinLowEdge(iy);
      fY(iy)=histogram->GetYaxis()->GetBinUpEdge(iy);
   }
   fZ.ResizeTo(nx,ny);
   for(int iy=1;iy<=ny;iy++) {
      for(int ix=1;ix<=nx;ix++) {
         if(fLogZ) {
            fZ(ix-1,iy-1)=TMath::Log(histogram->GetBinContent(ix,iy));
         } else {
            fZ(ix-1,iy-1)=histogram->GetBinContent(ix,iy);
         }
      }
   }
}

double Interpolate2d::Interpolate(double x,double y,
                                  bool useLogX,bool useLogY,
                                  double extrapol) {
   int nx=fX.GetNrows()-1;
   int ny=fY.GetNrows()-1;
   int ix,iy;
   for(ix=0;ix<=nx;ix++) {
      if(x<fX(ix)) break;
   }
   for(iy=0;iy<=ny;iy++) {
      if(y<fY(iy)) break;
   }
   if(ix<1) {
      ix=1;
      x=fX(0);
   } else if(ix>nx) {
      ix=nx;
      x=fX(nx);
   }
   if(iy<1) {
      iy=1;
      y=fY(0);
   } else if(iy>ny) {
      iy=ny;
      y=fY(ny);
   }
   double dx,dy;
   if(useLogX) {
      dx=TMath::Log(x*x/(fX(ix-1)*fX(ix)))/
         TMath::Log(fX(ix)/fX(ix-1));
   } else {
      dx=(2.*x-fX(ix-1)-fX(ix))/(fX(ix)-fX(ix-1));
   }
   if(useLogY) {
      dy=TMath::Log(y*y/(fY(iy-1)*fY(iy)))/
         TMath::Log(fY(iy)/fY(iy-1));
   } else {
      dy=(2.*y-fY(iy-1)-fY(iy))/(fY(iy)-fY(iy-1));
   }
   double z00=fZ(ix-1,iy-1);
   double z10,z01,z11;
   int idx=(dx<0)?-2:0;
   int idy=(dy<0)?-2:0;
   double fx=1.0;
   double fy=1.0;
   if((ix+idx<0)||(ix+idx>=nx)) {
      idx=-idx-2;
      fx= -extrapol;
   }
   if((iy+idy<0)||(iy+idy>=ny)) {
      idy=-idy-2;
      fy= -extrapol;
   }
   z10=fx*fZ(ix+idx,iy-1)+(1.-fx)*z00;
   z01=fy*fZ(ix-1,iy+idy)+(1.-fy)*z00;
   z11=fy*(fx*fZ(ix+idx,iy+idy)+(1.-fx)*fZ(ix-1,iy+idy))
      +(1.-fy)*(fx*fZ(ix+idx,iy-1)+(1.-fx)*z00);
   double zP=0.25*(z00+z01+z10+z11);
   double r=z00;
   dx=fabs(dx);
   dy=fabs(dy);
   if(dx>dy) {
      r+=0.5*dx*(z10-z00)+dy*(zP-0.5*(z10+z00));
   } else {
      r+=0.5*dy*(z01-z00)+dx*(zP-0.5*(z01+z00));
   }
   if(fLogZ) r=TMath::Exp(r);
   return r;
}

bool Interpolate2d::ReadData(const char *coeffFile) {
   std::ifstream in(coeffFile);
   if ( !in ) {
      cout<<"Error opening file: "<<coeffFile<<endl;
      exit(1);
   }

   return ReadData(in);
}

bool Interpolate2d::WriteData(const char *coeffFile) {
   std::ofstream out(coeffFile);
   return WriteData(out);
}

bool Interpolate2d::ReadData(std::istream &in) {
   int nx=0,ny=0;
   in>>nx;
   if((!in.fail())&&(nx>0)) {
      fX.ResizeTo(nx+1);
      for(int i=0;i<=nx;i++) {
         in>>fX(i);
      }
      if(in.fail()) nx=0;
   }
   in>>ny;
   if((!in.fail())&&(ny>0)) {
      fY.ResizeTo(ny+1);
      for(int i=0;i<=ny;i++) {
         in>>fY(i);
      }
      if(in.fail()) ny=0;
   }
   in>>fLogZ;
   if((!in.fail())&&(nx>0)&&(ny>0)) {
      fZ.ResizeTo(nx,ny);
      for(int j=0;j<ny;j++) {
         for(int i=0;i<nx;i++) {
            in>>fZ(i,j);
         }
      }
   }
   if((!in.fail())&&(nx>0)&&(ny>0)) {
      return true;
   }
   cout<<"Interplate2D. Error reading stream."<<endl;
   exit(1);
   return false;
}

bool Interpolate2d::WriteData(std::ostream &out) {
   int nx=fX.GetNrows()-1;
   int ny=fY.GetNrows()-1;
   out<<nx;
   for(int i=0;i<=nx;i++) {
      out<<" "<<fX(i);
   }
   out<<"\n"<<ny;
   for(int i=0;i<=ny;i++) {
      out<<" "<<fY(i);
   }
   out<<"\n"<<fLogZ<<"\n";
   for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) {
         out<<" "<<fZ(i,j);
      }
      out<<"\n";
   }
   return out.fail();
}
