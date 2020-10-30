#include <iostream>
#include <cmath>
#include "LumiFillList.h"
#include "H1Tools/H1RunList.h"

using namespace std;

void LumiRange::Add(int fill,int run,float lumi) {
  if((fFirstRun==0)&&(fLastRun==0)) {
    fFirstFill=fLastFill=fill;
    fFirstRun=fLastRun=run;
    fLumi=lumi;
    //fGood=good;
    runlist[run]=lumi;
  } else {
    fFirstFill=min(fFirstFill,fill);
    fLastFill=max(fLastFill,fill);
    fFirstRun=min(fFirstRun,run);
    fLastRun=max(fLastRun,run);
    fLumi += lumi;
    runlist[run]=lumi;
    /* if(fGood !=good) {
      std::cout<<"LumiRange::Add run "
               <<run<<" to range "<<fFirstRun<<"-"<<fLastRun
               <<" with good="<<good<<"!="<<fGood<<"\n"
               <<" this is considered as a fatal error\n";
      exit(0);
      } */
  }
}

void LumiRange::Add(LumiRange const &lfr) {
  if((fFirstRun==0)&&(fLastRun==0)) {
    *this=lfr;
  } else {
    fFirstFill=min(fFirstFill,lfr.FirstFill());
    fLastFill=max(fLastFill,lfr.LastFill());
    fFirstRun=min(fFirstRun,lfr.FirstRun());
    fLastRun=max(fLastRun,lfr.LastRun());
    fLumi += lfr.GetLumi();
    for(std::map<int,double>::const_iterator i=lfr.runlist.begin();
        i!=lfr.runlist.end();i++) {
      runlist[(*i).first]=(*i).second;
    }
  }
}

LumiFillList::LumiFillList(H1RunList *rl0) {
  int n0=rl0->GetNrOfRuns();
  //n1=rl1->GetNrOfRuns();
  int i0=0 /* ,i1=0 */;
  // set up a map of lumi fills
  int last_run=0;
  int last_fill=0;
  // bool last_selected=false;
  double sum_l0=0. /* ,sum_l1=0.0 */;
  while((i0<n0)/* ||(i1<n1) */) {
    H1RunLumi const *l0=((i0<n0)? rl0->At(i0) : 0);
    // H1RunLumi const *l1=((i1<n1)? rl1->At(i1) : 0);
    int run0=l0 ? l0->GetRunNumber() : 999999;
    // int run1=l1 ? l1->GetRunNumber() : 999999;
    int fill=(l0 ? l0->GetLumiFillNr() : 0);
    // cout<<" === "<<i0<<" "<<run0<<" === "<<i1<<" "<<run1<<" fill="<<fill<<"\n";
    int run=run0;
    float lumi=(l0 ? l0->GetCorrectedLumi() : 0.0);
    //bool selected=true;
    //if(run0==run1) {
      // good run
      i0++;
      //i1++;
      sum_l0 += lumi;
      //sum_l1 += lumi;
      /* } else {
      // bad run
      selected=false;
      if(run0<run1) {
        i0++;
      } else {
         cout<<"LumiFillList: error, run missing in list0\n";
        i1++;
        run=run1;
        fill=l1->GetLumiFillNr();
        lumi=l1->GetCorrectedLumi();
      }
      sum_l0 += lumi;
    }
    if(last_selected != selected) {
      cout<<"LumiFillList: run="<<run0<<" selected="<<selected<<"\n";
      last_selected = selected;
      } */
    if(fill<last_fill) {
       // cout<<"LumiFillList: fill not in order "<<fill<<" -> "<<last_fill<<"\n";
      fill=last_fill;
    }
    if((run<=last_run)||(fill<last_fill)) {
      cout<<"LumiFillList: run not in increasing order "
          <<run<<" "<<last_run<<" "<<fill<<" "<<last_fill<<"\n";
    }
    last_run=run;
    last_fill=fill;
    (*this)[fill].Add(fill,run,lumi /* ,selected */);
  }
  /* 
  cout<<"LumiFillList while filling: "<<sum_l0<<" "<<sum_l1<<"\n";
  sum_l0=0.;
  sum_l1=0.0;
  for(const_iterator i=begin();i!=end();i++) {
     LumiRange const &lr=(*i).second;
     if(lr.IsGood()) {
        sum_l1 += lr.GetLumi();
     }
     sum_l0 += lr.GetLumi();
  }
  cout<<"LumiFillList from fills: "<<sum_l0<<" "<<sum_l1<<"\n";
  sum_l0=0.;
  sum_l1=0.0;
  for(const_iterator i=begin();i!=end();i++) {
     LumiRange const &lr=(*i).second;
     Double_t l=0.0;
     for(std::map<int,double>::const_iterator j=lr.BeginRun();j !=lr.EndRun();
         j++) {
        l += (*j).second;
     }
     if(lr.IsGood()) {
        sum_l1 += l;
     }
     sum_l0 += l;
  }
  cout<<"LumiFillList from runs: "<<sum_l0<<" "<<sum_l1<<"\n";
  */
}
