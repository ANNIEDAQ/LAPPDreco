#ifndef WAVEFORM
#define WAVEFORM

#include "TemplateFit.h"
#include "TSplineFit.h"

using namespace std;

class Waveform : public TObject {
private:
  float Calculate_baseline(float *fvol); //
  float Calculate_baseline2(float *fvol, int nottrig);
  float Calculate_baseline3();
  float Waveform_Filter1(float CutoffFrequency,
                         float finput); // ideal low-pass filter
  float Waveform_Filter2(float CutoffFrequency, int T,
                         float finput); // Butterworth low-pass filter
  float Waveform_Filter3(float CutoffFrequency_low, float CutoffFrequency_high,
                         float finput);   // ideal low-pass filter
  float Filter_KURE(int T, float finput); // ideal low-pass filter
  void Calculate_fitrange();

  void Waveform_FFT();
  void Waveform_Fit();
  float CFD_Discriminator1();
  float CFD_Discriminator2();
  float GaussFit();
  float DoubleGaussFit(bool floatbaseline);
  float DoubleDoubleGaussFit();
  float DoTemplateFit(TemplateFit *mtf);
  float Calculate_amp();
  float Calculate_Charge();
  double Calculate_RisingTime();
  int Calculate_Peaks();
  int Calculate_Peaks_nnls();
  float FWPM(float frac);

public:
  int DoFFT;
  double Samplerate; // sampling rate
  float Deltat;   // time per sample
  int DimSize;    // Number of samples
  float baseline; // Base line of this waveform
  float goodbl;
  // float           istrig;
  float bnoise; // noise on the baseline
  float qfast;  // fast charge
  float qtot;   // total charge
  float gain;
  float amp;    // amp of the pulse
  float FWHM;   // full width at half max
  float FW20pM; // full width at 20 percent max
  double time;  // leading edge
  double risingtime;
  int MinimumTot;
  float TotThreshold;
  int npeaks;
  int evtno; // UVevent number
  int chno;
  int WavID; // each waveform must have an ID, ID = 100*evtno + Channel_Number
  double gmean;
  double gsigma;
  double gpeak;
  double gtime;
  double grisetime;
  double gchi2;
  double gdegfree;
  double goffset; // Gaussian parameters

  double ggmean;
  double ggsigma;
  double ggpeak;
  double ggmean2;
  double ggsigma2;
  double ggpeak2;
  double ggtime;
  double ggtime2;
  double ggrisetime;
  double ggchi2;
  double ggdegfree;
  double ggoffset; // Double-Gaussian parameters

  double ttime;
  double tscale;
  double tamp;
  double tcharge;
  double tchi2;
  double tdegfree;
  double tNS; // Template Fit Parameters

  // NNLS parameters
  int nnpeaks;
  double vtime[10];   // array of the times for each pulse
  double vamp[10];    // array of the amplitudes for each pulse
  double vcharge[10]; // array of charges for each pulse
  double vFWPM[10];   // array of FW20pMs for each pulse
  double vchi2[10];
  double vLowB[10];
  double vHiB[10];

  char WavName[100];
  char WavNameRaw[100];
  char xvectorname[100];
  char nnlsoutputname[100];
  char BgName[100];
  char PEDName[100];
  char SplineName[100];
  char SplineTitle[100];
  int PointsPerSpline; // Number of measured points for one spline
  char CFDName[100];
  char CFDTitle[100];
  char CFDSplineName[100];
  char CFDSplineTitle[100];
  float TimeStamp; // Time stamp of this waveform

  vector<float> vol;
  vector<float> vol_raw;
  vector<float> bg;
  vector<float> vol_fft;
  vector<float> re_fft;
  vector<float> im_fft;

  //		float		vol[256];		//Y value of the samples
  //		float		vol_fft[256];	//reconstruced Y value
  //		float		re_fft[256];
  //		float		im_fft[256];
  TH1D *hwav; // Waveform plot
  TH1D *hwav_raw;
  TH1D *xvector; // nnls output vector (weights)
  TH1D *nnlsoutput;
  int steps;
  int maxiter;
  // double          matrixA[632];
  size_t m;
  size_t n;
  int trigno;
  TH1D *hbg;
  TH1D *hcfd; // histogram for CFD
  TH1D *hpedhist;
  TSplineFit *fsplinewav; // SplineFit function
  TSplineFit *fsplinecfd;
  int IfDynamicWindow;
  double FitWindow_min; // Xmin for SplineFit
  double FitWindow_max; // Xmax for SplineFit
  double DisWindow_min;
  double DisWindow_max;
  double GaussRange_min; // Xmin for DoubleGaussFit
  double GaussRange_max; // Xmax for DoubleGaussFit

  float CutoffFrequency;
  double AmpThreshold;
  double Fraction_CFD; // fraction of Attenuation
  int Delay_CFD;       //
  const char *AnalysisOption;

  Waveform();
  virtual ~Waveform();
  void Clear(Option_t *option = "");
  float SetSamplingrate(float frate) {
    Samplerate = frate;
    return Samplerate;
  }
  float SetDeltat(float fdeltat) {
    Deltat = fdeltat;
    return Deltat;
  }
  int SetPointsPerSpline(int fm) {
    PointsPerSpline = fm;
    return PointsPerSpline;
  };
  int Setevtno(int fevtno) {
    evtno = fevtno;
    return evtno;
  }
  int Setchno(int fchno) {
    chno = fchno;
    return fchno;
  }
  int SetDimSize(int fDimSize) {
    DimSize = fDimSize;
    return DimSize;
  }
  float SetTimeStamp(float ft) {
    TimeStamp = ft;
    return TimeStamp;
  }
  double SetAmpThreshold(double fth) {
    AmpThreshold = fth;
    return AmpThreshold;
  }
  double SetFraction_CFD(double fraction) {
    Fraction_CFD = fraction;
    return Fraction_CFD;
  }
  int SetDelay_CFD(double delay) {
    Delay_CFD = delay;
    return Delay_CFD;
  }
  int EnableFFT() {
    DoFFT = 1;
    return DoFFT;
  }
  int EnableDynamicWindow() {
    IfDynamicWindow = 1;
    return IfDynamicWindow;
  }
  float SetCutoffFrequency(float f) {
    CutoffFrequency = f;
    return CutoffFrequency;
  };
  int SetMinimumTot(int fdt) {
    MinimumTot = fdt;
    return MinimumTot;
  }
  float SetTotThreshold(float thr) {
    TotThreshold = thr;
    return TotThreshold;
  }

  void Initialize(float *par, int nottrig); // call it with Set(&par[0])
  void Setup_nnls();
  void Calculate_Variables_nnls(int Npulses);
  void Analyze();
  void Analyze(TemplateFit *mtf);
  void SetGaussRange(double gausslow, double gausshi) {
    GaussRange_min = gausslow;
    GaussRange_max = gausshi;
  }
  void SetFitWindow(double fxmin, double fxmax) {
    FitWindow_min = fxmin;
    FitWindow_max = fxmax;
  }
  void SetDisWindow(double fxmin, double fxmax) {
    DisWindow_min = fxmin;
    DisWindow_max = fxmax;
  }
  const char *SetAnalysisOption(const char *foption) {
    AnalysisOption = foption;
    return AnalysisOption;
  }
  int Getevtno() const { return evtno; }
  float GetTimeStamp() const { return TimeStamp; }
  int GetDimSize() const { return DimSize; }
  float Getbaseline() const { return baseline; }
  float Getqfast() const { return qfast; }
  float Getqtot() const { return qtot; }
  float Getamp() const { return amp; }
  TH1D *GetWavHist() const { return hwav; }
  TSplineFit *GetWavSplineFit() const { return fsplinewav; }
  ClassDef(Waveform, 2)

  float slopeL;
  float slopeR;
  float interceptL;
  float interceptR;
  float timeL;
  float timeR;

  vector<int> HighBound;
  vector<int> LowBound;
};

#endif
