//----------APS Analysis
//----------Author's Name:Jingbo WANG
//----------Copyright:Those valid for ANL
//----------Modified:22/07/2014
#ifndef UVEVENT
#define UVEVENT

#include <TClonesArray.h>
#include <TVirtualFFT.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <vector>

using namespace std;

class UVevent : public TObject {

private:
  float WavSampleRate;
  int WavDimSize;
  float WavDeltat;

public:
  UVevent();
  virtual ~UVevent();
  void Clear(Option_t *option = "");

  int evtno;
  TClonesArray *fwav; //->array with all waveforms (not used yet)
  int fnwav;          // Number of waveforms
  int ClusterSize;    // Number of fired strips

  static TClonesArray *fgWaveforms;

  int cutflag; //
  vector<float> t;
  vector<float> fz;
  float ttrans; //
  float tdiff;
  float x; // Parallel postion
  float y; // Transverse position

  void Initialize();
  int Setevtno(int fevtno) {
    evtno = fevtno;
    return evtno;
  }
  float SetWavSampleRate(float frate) {
    WavSampleRate = frate;
    return WavSampleRate;
  }
  int SetWavDimSize(int fsize) {
    WavDimSize = fsize;
    return WavDimSize;
  }
  int SetCutflag(int flag) {
    cutflag = flag;
    return cutflag;
  }
  float SetTransitTime(float ftrans) {
    ttrans = ftrans;
    return ttrans;
  }
  float SetDifferentialTime(float ftdiff) {
    tdiff = ftdiff;
    return tdiff;
  }
  int SetWavNumber(int fn) {
    fnwav = fn;
    return fnwav;
  }
  int Getevtno() const { return evtno; }
  TClonesArray *GetWaveforms() const { return fwav; }
  int GetWavNumber() const { return fnwav; }
  void UVeventAnalysis();

  ClassDef(UVevent, 1)
};


#endif
