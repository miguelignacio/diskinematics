#ifndef __LumiFillList__
#define __LumiFillList__

#include <map>

class H1RunList;

class LumiRange {
 public:
  LumiRange(void) { Clear(); }
  void Clear(void) {
    fFirstRun=0;fLastRun=0;fFirstFill=0;fLastFill=0;fLumi=0.;
    // fGood=false;
    fBin=0;
  }
  void Add(int fill,int run, float lumi /* ,bool good */);
  void Add(LumiRange const&);
  inline int FirstRun(void) const { return fFirstRun; }
  inline int LastRun(void) const { return fLastRun; }
  inline int FirstFill(void) const { return fFirstFill; }
  inline int LastFill(void) const { return fLastFill; }
  inline double GetLumi(void) const { return fLumi; }
  //inline bool IsGood(void) const { return fGood; }
  inline void SetBin(int i) { fBin=i; }
  inline int GetBin(void) const { return fBin; }
  inline std::map<int,double>::const_iterator BeginRun(void) const {
    return runlist.begin(); }
  inline std::map<int,double>::const_iterator EndRun(void) const {
    return runlist.end(); }
  inline unsigned int Size(void) const { return runlist.size(); }
 protected:
  int fFirstFill,fLastFill;
  int fFirstRun,fLastRun;
  int fBin;
  double fLumi;
  std::map<int,double> runlist;
  //bool fGood;
};

class LumiFillList : public std::map<int,LumiRange> {
 public:
  inline LumiFillList(void) {}
  LumiFillList(H1RunList *rl);
};

#endif
