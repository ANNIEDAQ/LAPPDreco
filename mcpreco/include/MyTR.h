//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 11 11:16:41 2016 by ROOT version 5.34/05
// from TTree lappd/lappd
// found on file: blabla.root
//////////////////////////////////////////////////////////

#ifndef MyTR_h
#define MyTR_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxfwav = 4;

class MyTR {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //UVevent         *event;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Float_t         WavSampleRate;
   Int_t           WavDimSize;
   Float_t         WavDeltat;
   Int_t           evtno;
   Int_t           fwav_;
   UInt_t          fwav_fUniqueID[kMaxfwav];   //[fwav_]
   UInt_t          fwav_fBits[kMaxfwav];   //[fwav_]
   Int_t           fwav_DoFFT[kMaxfwav];   //[fwav_]
   Int_t           fwav_Samplerate[kMaxfwav];   //[fwav_]
   Float_t         fwav_Deltat[kMaxfwav];   //[fwav_]
   Int_t           fwav_DimSize[kMaxfwav];   //[fwav_]
   Float_t         fwav_baseline[kMaxfwav];   //[fwav_]
   Float_t         fwav_bnoise[kMaxfwav];   //[fwav_]
   Float_t         fwav_qfast[kMaxfwav];   //[fwav_]
   Float_t         fwav_qtot[kMaxfwav];   //[fwav_]
   Float_t         fwav_gain[kMaxfwav];   //[fwav_]
   Float_t         fwav_amp[kMaxfwav];   //[fwav_]
   Float_t         fwav_FWHM[kMaxfwav];   //[fwav_]
   Float_t         fwav_FW20pM[kMaxfwav];   //[fwav_]
   Double_t        fwav_time[kMaxfwav];   //[fwav_]
   Double_t        fwav_risingtime[kMaxfwav];   //[fwav_]
   Int_t           fwav_MinimumTot[kMaxfwav];   //[fwav_]
   Float_t         fwav_TotThreshold[kMaxfwav];   //[fwav_]
   Int_t           fwav_npeaks[kMaxfwav];   //[fwav_]
   Int_t           fwav_evtno[kMaxfwav];   //[fwav_]
   Int_t           fwav_chno[kMaxfwav];   //[fwav_]
   Int_t           fwav_WavID[kMaxfwav];   //[fwav_]
   Double_t        fwav_gmean[kMaxfwav];   //[fwav_]
   Double_t        fwav_gsigma[kMaxfwav];   //[fwav_]
   Double_t        fwav_gpeak[kMaxfwav];   //[fwav_]
   Double_t        fwav_gtime[kMaxfwav];   //[fwav_]
   Double_t        fwav_grisetime[kMaxfwav];   //[fwav_]
   Double_t        fwav_gchi2[kMaxfwav];   //[fwav_]
   Double_t        fwav_gdegfree[kMaxfwav];   //[fwav_]
   Double_t        fwav_goffset[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggmean[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggsigma[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggpeak[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggmean2[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggsigma2[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggpeak2[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggtime[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggtime2[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggrisetime[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggchi2[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggdegfree[kMaxfwav];   //[fwav_]
   Double_t        fwav_ggoffset[kMaxfwav];   //[fwav_]
   Double_t        fwav_ttime[kMaxfwav];   //[fwav_]
   Double_t        fwav_tscale[kMaxfwav];   //[fwav_]
   Double_t        fwav_tamp[kMaxfwav];   //[fwav_]
   Double_t        fwav_tcharge[kMaxfwav];   //[fwav_]
   Double_t        fwav_tchi2[kMaxfwav];   //[fwav_]
   Double_t        fwav_tdegfree[kMaxfwav];   //[fwav_]
   Double_t        fwav_tNS[kMaxfwav];   //[fwav_]
   Char_t          fwav_WavName[kMaxfwav][100];   //[fwav_]
   Char_t          fwav_SplineName[kMaxfwav][100];   //[fwav_]
   Char_t          fwav_SplineTitle[kMaxfwav][100];   //[fwav_]
   Int_t           fwav_PointsPerSpline[kMaxfwav];   //[fwav_]
   Char_t          fwav_CFDName[kMaxfwav][100];   //[fwav_]
   Char_t          fwav_CFDTitle[kMaxfwav][100];   //[fwav_]
   Char_t          fwav_CFDSplineName[kMaxfwav][100];   //[fwav_]
   Char_t          fwav_CFDSplineTitle[kMaxfwav][100];   //[fwav_]
   Float_t         fwav_TimeStamp[kMaxfwav];   //[fwav_]
   vector<float>   fwav_vol[kMaxfwav];
   vector<float>   fwav_vol_fft[kMaxfwav];
   vector<float>   fwav_re_fft[kMaxfwav];
   vector<float>   fwav_im_fft[kMaxfwav];
   Int_t           fwav_IfDynamicWindow[kMaxfwav];   //[fwav_]
   Double_t        fwav_FitWindow_min[kMaxfwav];   //[fwav_]
   Double_t        fwav_FitWindow_max[kMaxfwav];   //[fwav_]
   Double_t        fwav_DisWindow_min[kMaxfwav];   //[fwav_]
   Double_t        fwav_DisWindow_max[kMaxfwav];   //[fwav_]
   Double_t        fwav_GaussRange_min[kMaxfwav];   //[fwav_]
   Double_t        fwav_GaussRange_max[kMaxfwav];   //[fwav_]
   Float_t         fwav_CutoffFrequency[kMaxfwav];   //[fwav_]
   Double_t        fwav_AmpThreshold[kMaxfwav];   //[fwav_]
   Double_t        fwav_Fraction_CFD[kMaxfwav];   //[fwav_]
   Int_t           fwav_Delay_CFD[kMaxfwav];   //[fwav_]
   char*           fwav_AnalysisOption[kMaxfwav];   //[fwav_]
   Int_t           fnwav;
   Int_t           ClusterSize;
   Int_t           cutflag;
   vector<float>   t;
   vector<float>   fz;
   Float_t         ttrans;
   Float_t         tdiff;
   Float_t         x;
   Float_t         y;

   // List of branches
   TBranch        *b_event_fUniqueID;   //!
   TBranch        *b_event_fBits;   //!
   TBranch        *b_event_WavSampleRate;   //!
   TBranch        *b_event_WavDimSize;   //!
   TBranch        *b_event_WavDeltat;   //!
   TBranch        *b_event_evtno;   //!
   TBranch        *b_event_fwav_;   //!
   TBranch        *b_fwav_fUniqueID;   //!
   TBranch        *b_fwav_fBits;   //!
   TBranch        *b_fwav_DoFFT;   //!
   TBranch        *b_fwav_Samplerate;   //!
   TBranch        *b_fwav_Deltat;   //!
   TBranch        *b_fwav_DimSize;   //!
   TBranch        *b_fwav_baseline;   //!
   TBranch        *b_fwav_bnoise;   //!
   TBranch        *b_fwav_qfast;   //!
   TBranch        *b_fwav_qtot;   //!
   TBranch        *b_fwav_gain;   //!
   TBranch        *b_fwav_amp;   //!
   TBranch        *b_fwav_FWHM;   //!
   TBranch        *b_fwav_FW20pM;   //!
   TBranch        *b_fwav_time;   //!
   TBranch        *b_fwav_risingtime;   //!
   TBranch        *b_fwav_MinimumTot;   //!
   TBranch        *b_fwav_TotThreshold;   //!
   TBranch        *b_fwav_npeaks;   //!
   TBranch        *b_fwav_evtno;   //!
   TBranch        *b_fwav_chno;   //!
   TBranch        *b_fwav_WavID;   //!
   TBranch        *b_fwav_gmean;   //!
   TBranch        *b_fwav_gsigma;   //!
   TBranch        *b_fwav_gpeak;   //!
   TBranch        *b_fwav_gtime;   //!
   TBranch        *b_fwav_grisetime;   //!
   TBranch        *b_fwav_gchi2;   //!
   TBranch        *b_fwav_gdegfree;   //!
   TBranch        *b_fwav_goffset;   //!
   TBranch        *b_fwav_ggmean;   //!
   TBranch        *b_fwav_ggsigma;   //!
   TBranch        *b_fwav_ggpeak;   //!
   TBranch        *b_fwav_ggmean2;   //!
   TBranch        *b_fwav_ggsigma2;   //!
   TBranch        *b_fwav_ggpeak2;   //!
   TBranch        *b_fwav_ggtime;   //!
   TBranch        *b_fwav_ggtime2;   //!
   TBranch        *b_fwav_ggrisetime;   //!
   TBranch        *b_fwav_ggchi2;   //!
   TBranch        *b_fwav_ggdegfree;   //!
   TBranch        *b_fwav_ggoffset;   //!
   TBranch        *b_fwav_ttime;   //!
   TBranch        *b_fwav_tscale;   //!
   TBranch        *b_fwav_tamp;   //!
   TBranch        *b_fwav_tcharge;   //!
   TBranch        *b_fwav_tchi2;   //!
   TBranch        *b_fwav_tdegfree;   //!
   TBranch        *b_fwav_tNS;   //!
   TBranch        *b_fwav_WavName;   //!
   TBranch        *b_fwav_SplineName;   //!
   TBranch        *b_fwav_SplineTitle;   //!
   TBranch        *b_fwav_PointsPerSpline;   //!
   TBranch        *b_fwav_CFDName;   //!
   TBranch        *b_fwav_CFDTitle;   //!
   TBranch        *b_fwav_CFDSplineName;   //!
   TBranch        *b_fwav_CFDSplineTitle;   //!
   TBranch        *b_fwav_TimeStamp;   //!
   TBranch        *b_fwav_vol;   //!
   TBranch        *b_fwav_vol_fft;   //!
   TBranch        *b_fwav_re_fft;   //!
   TBranch        *b_fwav_im_fft;   //!
   TBranch        *b_fwav_IfDynamicWindow;   //!
   TBranch        *b_fwav_FitWindow_min;   //!
   TBranch        *b_fwav_FitWindow_max;   //!
   TBranch        *b_fwav_DisWindow_min;   //!
   TBranch        *b_fwav_DisWindow_max;   //!
   TBranch        *b_fwav_GaussRange_min;   //!
   TBranch        *b_fwav_GaussRange_max;   //!
   TBranch        *b_fwav_CutoffFrequency;   //!
   TBranch        *b_fwav_AmpThreshold;   //!
   TBranch        *b_fwav_Fraction_CFD;   //!
   TBranch        *b_fwav_Delay_CFD;   //!
   TBranch        *b_fwav_AnalysisOption;   //!
   TBranch        *b_event_fnwav;   //!
   TBranch        *b_event_ClusterSize;   //!
   TBranch        *b_event_cutflag;   //!
   TBranch        *b_event_t;   //!
   TBranch        *b_event_fz;   //!
   TBranch        *b_event_ttrans;   //!
   TBranch        *b_event_tdiff;   //!
   TBranch        *b_event_x;   //!
   TBranch        *b_event_y;   //!

   MyTR(TTree *tree=0);
   virtual ~MyTR();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyTR_cxx
MyTR::MyTR(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("blabla.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("blabla.root");
      }
      f->GetObject("lappd",tree);

   }
   Init(tree);
}

MyTR::~MyTR()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyTR::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyTR::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyTR::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_event_fBits);
   fChain->SetBranchAddress("WavSampleRate", &WavSampleRate, &b_event_WavSampleRate);
   fChain->SetBranchAddress("WavDimSize", &WavDimSize, &b_event_WavDimSize);
   fChain->SetBranchAddress("WavDeltat", &WavDeltat, &b_event_WavDeltat);
   fChain->SetBranchAddress("evtno", &evtno, &b_event_evtno);
   fChain->SetBranchAddress("fwav", &fwav_, &b_event_fwav_);
   fChain->SetBranchAddress("fwav.fUniqueID", fwav_fUniqueID, &b_fwav_fUniqueID);
   fChain->SetBranchAddress("fwav.fBits", fwav_fBits, &b_fwav_fBits);
   fChain->SetBranchAddress("fwav.DoFFT", fwav_DoFFT, &b_fwav_DoFFT);
   fChain->SetBranchAddress("fwav.Samplerate", fwav_Samplerate, &b_fwav_Samplerate);
   fChain->SetBranchAddress("fwav.Deltat", fwav_Deltat, &b_fwav_Deltat);
   fChain->SetBranchAddress("fwav.DimSize", fwav_DimSize, &b_fwav_DimSize);
   fChain->SetBranchAddress("fwav.baseline", fwav_baseline, &b_fwav_baseline);
   fChain->SetBranchAddress("fwav.bnoise", fwav_bnoise, &b_fwav_bnoise);
   fChain->SetBranchAddress("fwav.qfast", fwav_qfast, &b_fwav_qfast);
   fChain->SetBranchAddress("fwav.qtot", fwav_qtot, &b_fwav_qtot);
   fChain->SetBranchAddress("fwav.gain", fwav_gain, &b_fwav_gain);
   fChain->SetBranchAddress("fwav.amp", fwav_amp, &b_fwav_amp);
   fChain->SetBranchAddress("fwav.FWHM", fwav_FWHM, &b_fwav_FWHM);
   fChain->SetBranchAddress("fwav.FW20pM", fwav_FW20pM, &b_fwav_FW20pM);
   fChain->SetBranchAddress("fwav.time", fwav_time, &b_fwav_time);
   fChain->SetBranchAddress("fwav.risingtime", fwav_risingtime, &b_fwav_risingtime);
   fChain->SetBranchAddress("fwav.MinimumTot", fwav_MinimumTot, &b_fwav_MinimumTot);
   fChain->SetBranchAddress("fwav.TotThreshold", fwav_TotThreshold, &b_fwav_TotThreshold);
   fChain->SetBranchAddress("fwav.npeaks", fwav_npeaks, &b_fwav_npeaks);
   fChain->SetBranchAddress("fwav.evtno", fwav_evtno, &b_fwav_evtno);
   fChain->SetBranchAddress("fwav.chno", fwav_chno, &b_fwav_chno);
   fChain->SetBranchAddress("fwav.WavID", fwav_WavID, &b_fwav_WavID);
   fChain->SetBranchAddress("fwav.gmean", fwav_gmean, &b_fwav_gmean);
   fChain->SetBranchAddress("fwav.gsigma", fwav_gsigma, &b_fwav_gsigma);
   fChain->SetBranchAddress("fwav.gpeak", fwav_gpeak, &b_fwav_gpeak);
   fChain->SetBranchAddress("fwav.gtime", fwav_gtime, &b_fwav_gtime);
   fChain->SetBranchAddress("fwav.grisetime", fwav_grisetime, &b_fwav_grisetime);
   fChain->SetBranchAddress("fwav.gchi2", fwav_gchi2, &b_fwav_gchi2);
   fChain->SetBranchAddress("fwav.gdegfree", fwav_gdegfree, &b_fwav_gdegfree);
   fChain->SetBranchAddress("fwav.goffset", fwav_goffset, &b_fwav_goffset);
   fChain->SetBranchAddress("fwav.ggmean", fwav_ggmean, &b_fwav_ggmean);
   fChain->SetBranchAddress("fwav.ggsigma", fwav_ggsigma, &b_fwav_ggsigma);
   fChain->SetBranchAddress("fwav.ggpeak", fwav_ggpeak, &b_fwav_ggpeak);
   fChain->SetBranchAddress("fwav.ggmean2", fwav_ggmean2, &b_fwav_ggmean2);
   fChain->SetBranchAddress("fwav.ggsigma2", fwav_ggsigma2, &b_fwav_ggsigma2);
   fChain->SetBranchAddress("fwav.ggpeak2", fwav_ggpeak2, &b_fwav_ggpeak2);
   fChain->SetBranchAddress("fwav.ggtime", fwav_ggtime, &b_fwav_ggtime);
   fChain->SetBranchAddress("fwav.ggtime2", fwav_ggtime2, &b_fwav_ggtime2);
   fChain->SetBranchAddress("fwav.ggrisetime", fwav_ggrisetime, &b_fwav_ggrisetime);
   fChain->SetBranchAddress("fwav.ggchi2", fwav_ggchi2, &b_fwav_ggchi2);
   fChain->SetBranchAddress("fwav.ggdegfree", fwav_ggdegfree, &b_fwav_ggdegfree);
   fChain->SetBranchAddress("fwav.ggoffset", fwav_ggoffset, &b_fwav_ggoffset);
   fChain->SetBranchAddress("fwav.ttime", fwav_ttime, &b_fwav_ttime);
   fChain->SetBranchAddress("fwav.tscale", fwav_tscale, &b_fwav_tscale);
   fChain->SetBranchAddress("fwav.tamp", fwav_tamp, &b_fwav_tamp);
   fChain->SetBranchAddress("fwav.tcharge", fwav_tcharge, &b_fwav_tcharge);
   fChain->SetBranchAddress("fwav.tchi2", fwav_tchi2, &b_fwav_tchi2);
   fChain->SetBranchAddress("fwav.tdegfree", fwav_tdegfree, &b_fwav_tdegfree);
   fChain->SetBranchAddress("fwav.tNS", fwav_tNS, &b_fwav_tNS);
   fChain->SetBranchAddress("fwav.WavName[100]", fwav_WavName, &b_fwav_WavName);
   fChain->SetBranchAddress("fwav.SplineName[100]", fwav_SplineName, &b_fwav_SplineName);
   fChain->SetBranchAddress("fwav.SplineTitle[100]", fwav_SplineTitle, &b_fwav_SplineTitle);
   fChain->SetBranchAddress("fwav.PointsPerSpline", fwav_PointsPerSpline, &b_fwav_PointsPerSpline);
   fChain->SetBranchAddress("fwav.CFDName[100]", fwav_CFDName, &b_fwav_CFDName);
   fChain->SetBranchAddress("fwav.CFDTitle[100]", fwav_CFDTitle, &b_fwav_CFDTitle);
   fChain->SetBranchAddress("fwav.CFDSplineName[100]", fwav_CFDSplineName, &b_fwav_CFDSplineName);
   fChain->SetBranchAddress("fwav.CFDSplineTitle[100]", fwav_CFDSplineTitle, &b_fwav_CFDSplineTitle);
   fChain->SetBranchAddress("fwav.TimeStamp", fwav_TimeStamp, &b_fwav_TimeStamp);
   fChain->SetBranchAddress("fwav.vol", fwav_vol, &b_fwav_vol);
   fChain->SetBranchAddress("fwav.vol_fft", fwav_vol_fft, &b_fwav_vol_fft);
   fChain->SetBranchAddress("fwav.re_fft", fwav_re_fft, &b_fwav_re_fft);
   fChain->SetBranchAddress("fwav.im_fft", fwav_im_fft, &b_fwav_im_fft);
   fChain->SetBranchAddress("fwav.IfDynamicWindow", fwav_IfDynamicWindow, &b_fwav_IfDynamicWindow);
   fChain->SetBranchAddress("fwav.FitWindow_min", fwav_FitWindow_min, &b_fwav_FitWindow_min);
   fChain->SetBranchAddress("fwav.FitWindow_max", fwav_FitWindow_max, &b_fwav_FitWindow_max);
   fChain->SetBranchAddress("fwav.DisWindow_min", fwav_DisWindow_min, &b_fwav_DisWindow_min);
   fChain->SetBranchAddress("fwav.DisWindow_max", fwav_DisWindow_max, &b_fwav_DisWindow_max);
   fChain->SetBranchAddress("fwav.GaussRange_min", fwav_GaussRange_min, &b_fwav_GaussRange_min);
   fChain->SetBranchAddress("fwav.GaussRange_max", fwav_GaussRange_max, &b_fwav_GaussRange_max);
   fChain->SetBranchAddress("fwav.CutoffFrequency", fwav_CutoffFrequency, &b_fwav_CutoffFrequency);
   fChain->SetBranchAddress("fwav.AmpThreshold", fwav_AmpThreshold, &b_fwav_AmpThreshold);
   fChain->SetBranchAddress("fwav.Fraction_CFD", fwav_Fraction_CFD, &b_fwav_Fraction_CFD);
   fChain->SetBranchAddress("fwav.Delay_CFD", fwav_Delay_CFD, &b_fwav_Delay_CFD);
   fChain->SetBranchAddress("fwav.AnalysisOption", fwav_AnalysisOption, &b_fwav_AnalysisOption);
   fChain->SetBranchAddress("fnwav", &fnwav, &b_event_fnwav);
   fChain->SetBranchAddress("ClusterSize", &ClusterSize, &b_event_ClusterSize);
   fChain->SetBranchAddress("cutflag", &cutflag, &b_event_cutflag);
   fChain->SetBranchAddress("t", &t, &b_event_t);
   fChain->SetBranchAddress("fz", &fz, &b_event_fz);
   fChain->SetBranchAddress("ttrans", &ttrans, &b_event_ttrans);
   fChain->SetBranchAddress("tdiff", &tdiff, &b_event_tdiff);
   fChain->SetBranchAddress("x", &x, &b_event_x);
   fChain->SetBranchAddress("y", &y, &b_event_y);
   Notify();
}

Bool_t MyTR::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyTR::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyTR::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyTR_cxx
