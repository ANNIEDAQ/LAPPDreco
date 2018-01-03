#ifndef TEMPLATEFIT
#define TEMPLATEFIT

#include <TObject.h>
#include <TFile.h>
#include <TF1.h>
#include <TROOT.h>
#include <TH1F.h>
#include "TClonesArray.h"
#include "TMinuit.h"
#include "TSplineFit.h"

class TemplateFit: public TObject { 

private: 

  TH1D* _ptemplate;
  TH1D* _phist;
  TH1D* _fithist;
  bool first_template;
  int tempcount;
  int instancecount;
  double _noise;
  TMinuit* fMinuitScaleAndOffset;

public: 

  //  static void fit_template_fcn(Int_t&, Double_t*, Double_t &f, Double_t* par, Int_t);

  TemplateFit(TString tname);
  TemplateFit(); 
  virtual ~TemplateFit(); 
  void Clear(Option_t *option ="");
  void AddTemplate(TH1D* ntemp);
  void AddTemplateNoNorm(TH1D* ntemp);
  void AddSmoothedTemplate(TH1D* stemp);	
  void SetTemplate(TH1D* ntemp);
  void Combine(TH1D* ohist,double mshift,double mscale);
  void SetPulseHist(TH1D*);
  void SetFitHists(TH1D* sfhist, TH1D* sthist);
  void SetFitRange(Double_t lr, Double_t hr);
  void SetTempRange(Double_t lr, Double_t hr);
  void SetNoise(double tnoise);
  void SetStartingValues(double tstart, double scalestart);
  void CalcFixedDegFree(TH1D* ttemp);
  bool IsFirstTemp();
  Double_t CalculateFCN(Double_t &f, Double_t *par);
  TH1D* ShiftAndScale(TH1D* ohist, double mshift, double mscale);
  TH1D* ReturnTemplate(bool NormIt);
  TH1D* ReturnGTemplate();
  TH1D* ReturnPulseHist();
  void FitWithTemplate(TH1D* fhist, TH1D* temphist, double &moffset, double &mscale, double &mchi2);
  void FitWithTemplate(TH1D* fhist, double &moffset, double &mscale, double &mchi2);


  ClassDef(TemplateFit,1) 

}; 

extern TSplineFit* tshist;
extern TH1D* gtemplate;
extern TH1D* ghist;
extern bool gusefitrange;
extern bool gusetemprange;
extern double glrange;
extern double ghrange;
extern int glbin;
extern int ghbin;
extern double gtlrange;
extern double gthrange;
extern double gchi2,gdegfree;
extern int gtlbin;
extern int gthbin;
extern int gdegfree_fixed;

#endif  
