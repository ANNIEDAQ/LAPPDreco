#include "Waveform.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "nnls/nnls.h"
#include "TSpectrum.h"

using namespace std;

Waveform::Waveform() {
  //	cout<<"Class Waveform"<<endl;
  //	hwav = new TH1F("wav","wav",DimSize,0,103.2);
  fsplinewav = 0;
  fsplinecfd = 0;
  hwav = 0;
  hwav_raw = 0;
  xvector = 0;
  nnlsoutput = 0;
  hbg = 0;
  hcfd = 0;
  hpedhist = 0;
  Samplerate = 1e10; // Hz
  Deltat = 100;
  DimSize = 256;
  baseline = 0;
  CutoffFrequency = 5e8;
  PointsPerSpline = 5;
  Fraction_CFD = 0.5;
  Delay_CFD = 100;
  DisWindow_min = 0;
  DisWindow_max = DimSize * Deltat;
  FitWindow_min = 0;
  FitWindow_max = DimSize * Deltat;
  GaussRange_min = 0;
  GaussRange_max = DimSize * Deltat;

  npeaks = 0;
  MinimumTot = 100;
  DoFFT = 0;
  IfDynamicWindow = 0;

  for (int i = 0; i < 10; i++) {
    vtime[i] = 0;
    vamp[i] = 0;
    vcharge[i] = 0;
    vFWPM[i] = 0;
    vchi2[i] = 0;
    vLowB[i] = 0;
    vHiB[i] = 0;
  }
}

Waveform::~Waveform() {
  Clear();
  //	delete hwav;
  //	delete froot;
  //	delete fsplinewav;
}

void Waveform::Clear(Option_t * /*option*/) {
  // Note that we intend on using TClonesArray::ConstructedAt, so we do not
  // need to delete any of the arrays.
  TObject::Clear();
  delete hwav;
  hwav = 0;
  delete hbg;
  hbg = 0;
  delete hwav_raw;
  hwav_raw = 0;
  delete xvector;
  xvector = 0;
  delete nnlsoutput;
  nnlsoutput = 0;
  delete hcfd;
  hcfd = 0;
  delete hpedhist;
  hpedhist = 0;
  delete fsplinewav;
  fsplinewav = 0;
  delete fsplinecfd;
  fsplinecfd = 0;
  bg.clear();
  vol.clear();
  vol_raw.clear();
  vol_fft.clear();
  re_fft.clear();
  im_fft.clear();
}

void Waveform::Initialize(float *fvol,
                          int nottrig) // call it with Set(&fvol[0])
{
  // Set values for the Waveform members
  trigno = nottrig;
  cout << "nottrig " << nottrig << endl;

  steps = 5;    // how finely the template is binned
  maxiter = 10; // the maximum number of iterations

  sprintf(WavName, "wav_ch%d_%d", chno, evtno);
  sprintf(WavNameRaw, "wavraw_ch%d_%d", chno, evtno);
  sprintf(CFDName, "CFD_ch%d_%d", chno, evtno);
  sprintf(xvectorname, "xvec_ch%d_%d", chno, evtno);
  sprintf(nnlsoutputname, "nnlso_ch%d_%d", chno, evtno);
  sprintf(PEDName, "PED_ch%d_%d", chno, evtno);
  sprintf(BgName, "Bg_ch%d_%d", chno, evtno);
  //  hwav = new TH1D(WavName, WavName, DimSize, 0, 10);
  hwav = new TH1D(WavName, WavName, DimSize, 0, DimSize * Deltat);
  hwav_raw = new TH1D(WavNameRaw, WavNameRaw, DimSize, 0, DimSize * Deltat);
  xvector = new TH1D(xvectorname, xvectorname, (DimSize + 60) * steps, 0,
                     (DimSize + 60) * steps);
  nnlsoutput = new TH1D(nnlsoutputname, nnlsoutputname, DimSize * steps, 0,
                        DimSize * Deltat);
  hbg = new TH1D(BgName, BgName, DimSize, 0, DimSize * Deltat);
  hcfd = new TH1D(CFDName, CFDName, DimSize, 0, DimSize * Deltat);
  hpedhist = new TH1D(PEDName, PEDName, 100, -5.0, 5.0);

  for (int i = 0; i < DimSize; i++) {
    vol_raw.push_back(*(fvol + i));
    hwav_raw->SetBinContent(i + 1, vol_raw[i]);
    hwav_raw->SetBinError(i + 1, 0.3);
  }

  /*
  if(nottrig){
    Calculate_baseline3();
  }
  */

  TSpectrum *s = new TSpectrum(1);

  //	Calculate_baseline(&fvol[0]);
  Calculate_baseline2(&fvol[0], nottrig);
  //	TF1 *fsin = new TF1("fsin", "40*sin(2*3.14*1e9*x)", 0,1);
  double Fs = 1e10;       // sampling rate = 10GHz
  double Ts = 1.e12 / Fs; // 100ps/point

  for (int i = 0; i < DimSize; i++) {
    vol.push_back(*(fvol + i));
    // bg.push_back(-(*(fvol+i)-baseline));
  }

  //  s->TSpectrum::Background(&(bg[0]),DimSize,8,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE,
  //  TSpectrum::kBackSmoothing3,kFALSE);
  //  s->TSpectrum::Background(&(bg[0]),DimSize,4,TSpectrum::kBackIncreasingWindow,TSpectrum::kBackOrder2,kTRUE,
  //  TSpectrum::kBackSmoothing3,kFALSE);

  //  TF1* sinit0 = new
  //  TF1("sinit0","([0]*sin(.00066*x+[1]))+[2]",0,DimSize*Deltat);
  //  TF1* sinit1 = new
  //  TF1("sinit1","([0]*sin(.00055*x+[1]))+[2]",0,DimSize*Deltat);
  TF1 *sinit0 = new TF1("sinit0", "([0]*sin([3]*x+[1]))+[2]", 0,
                        DimSize * Deltat); // 0.000557
  TF1 *sinit1 = new TF1("sinit1", "([0]*sin([3]*x+[1]))+[2]", 0,
                        DimSize * Deltat); // 0.000556

  sinit0->SetParameter(3, 0.00055);
  sinit1->SetParameter(3, 0.00055);
  sinit0->SetParLimits(3, 0.0005, 0.0006);
  sinit1->SetParLimits(3, 0.0005, 0.0006);

  bool fitsinewave = false; // CONTROLS WHETHER OR NOT THE SINE WAVE IS TAKEN
                            // OUT

  //  if(nottrig==0) hwav_raw->Fit("sinit0","q","",40000,80000);
  //  if(nottrig==1) hwav_raw->Fit("sinit1","q","",0,DimSize*Deltat);
  //  if(nottrig==2) hwav_raw->Fit("sinit1","q","",0,DimSize*Deltat);
  if (fitsinewave == true) {
    if (nottrig == 0)
      hwav_raw->Fit("sinit0", "q", "", 32000, 50000);
    if (nottrig == 1)
      hwav_raw->Fit("sinit1", "q", "", 32000, 50000);
    if (nottrig == 2)
      hwav_raw->Fit("sinit1", "q", "", 32000, 50000);
  }

  bool clipit = true;
  // double cliplow=56000;
  // double cliphi=65000;
  // double cliplow=56000;
  // double cliphi=66000;
  double cliplow = 50000;
  double cliphi = 55000;

  if (fitsinewave == true) {
    for (int i = 0; i < DimSize; i++) {
      // vol.push_back(*(fvol+i)-baseline);
      //		vol[i] = *(fvol+i)-baseline;
      // hwav->SetBinContent(i+1,*(fvol+i)-baseline);
      if (nottrig == 0)
        hwav->SetBinContent(i + 1,
                            vol[i] - (sinit0->Eval(hwav->GetBinCenter(i + 1))));
      if (nottrig == 1)
        hwav->SetBinContent(i + 1,
                            vol[i] - (sinit1->Eval(hwav->GetBinCenter(i + 1))));
      if (nottrig == 2)
        hwav->SetBinContent(i + 1,
                            vol[i] - (sinit1->Eval(hwav->GetBinCenter(i + 1))));
      if (nottrig < 0)
        hwav->SetBinContent(i + 1, vol[i]);
      // hwav->SetBinError(i+1,0.3);
      // hbg->SetBinContent(i+1,-bg[i]);
      if (clipit) {
        if ((hwav->GetBinCenter(i + 1)) < cliplow)
          hwav->SetBinContent(i + 1, 0.0);
        if ((hwav->GetBinCenter(i + 1)) > cliphi)
          hwav->SetBinContent(i + 1, 0.0);
      }
    }
  } else {
    for (int i = 0; i < DimSize; i++) {
      if (nottrig == 0)
        hwav->SetBinContent(i + 1, vol[i]);
      if (nottrig == 1)
        hwav->SetBinContent(i + 1, vol[i]);
      if (nottrig == 2)
        hwav->SetBinContent(i + 1, vol[i]);
      if (nottrig < 0)
        hwav->SetBinContent(i + 1, vol[i]);
    }
  }
  AmpThreshold = 60;
  Fraction_CFD = 0.5;
  Delay_CFD = 10;
  AnalysisOption = "normal";
}

void Waveform::Setup_nnls() {

  // creates matrixA array and xvector histogram

  cout << "here0" << endl;
  m = (DimSize + 60) * steps;
  n = (DimSize + 60) * steps;
  double *matrixA2 = new double[m * n];

  cout << "here1" << endl;
  // fill matrixA (array version of A)
  TFile *tt = new TFile("pulsecharacteristics.root");
  TH1D *tempp = (TH1D *)tt->Get("templatepulse;1");
  int bincount = 1;
  for (int i = 0; i < (int) m; i++) {
    bincount = 1;
    for (int j = 0; j < (int) m; j++) {
      if (bincount > 0 && bincount <= 30 * steps && j >= i &&
          j <= (i + 30 * steps)) {
        if (tempp->GetBinContent((100 / steps) * (bincount)) < 0) {
          matrixA2[m * i + j] =
              tempp->GetBinContent((100 / steps) * (bincount));
        } else {
          matrixA2[m * i + j] = 0;
        }
        bincount++;
      } else {
        matrixA2[m * i + j] = 0;
      }
    }
  }

  nsNNLS::denseMatrix *matA = new nsNNLS::denseMatrix(m, n, matrixA2);

  double *vv = new double[m];

  // fill vv (array version of b)
  for (int i = 0; i < (int)m; i++) {
    if (i < 30 * steps || i > ((int)m - 30 * steps)) {
      vv[i] = 0;
    } else {
      vv[i] = hwav_raw->GetBinContent((i - 30 * steps) / steps);
    }
  }
  nsNNLS::vector *vecB = new nsNNLS::vector(m, vv);
  nsNNLS::vector *vecX;

  nsNNLS::nnls *solver = new nsNNLS::nnls(matA, vecB, maxiter);

  int flag;
  flag = solver->optimize();

  vecX = solver->getSolution();

  // fill xvector histogram using vecX (vector version of x)
  for (int i = 0; i < (int) m; i++) {
    if (i + 1 - 28 * steps > 0 && i > 30 * steps && i < ((int) m - 30 * steps)) {
      xvector->SetBinContent(i + 1 - 28 * steps, vecX->getData()[i]);
    }
  }

  // create output based on templates and xvector
  double outputa[m];
  for (int i = 0; i < (int) m; i++)
    outputa[i] = 0;

  for (int i = 0; i < matA->nrows(); i++) {
    for (int j = 0; j < matA->ncols(); j++) {
      outputa[j] += (xvector->GetBinContent(i + 1)) * (matA->get(i, j));
    }
  }

  for (int i = 0; i < matA->ncols(); i++) {
    if (i > 30 * steps && i < ((int) m - 30 * steps))
      nnlsoutput->SetBinContent(i + 1 - 30 * steps, outputa[i]);
  }

  delete tt;
  delete solver;
  delete vecX;
  delete vecB;
  delete matA;
  delete vv;
  tempp = 0;
  delete matrixA2;

  cout << "made it" << endl;
}

float Waveform::Calculate_baseline(float *fvol) {
  TH1D *hwavbaseline =
      new TH1D("hwavbaseline", "hwavbaseline", DimSize, 0, DimSize * Deltat);
  for (int i = 0; i < 200; i++) {
    //		SetVol(i, *(fvol+i));
    hwavbaseline->SetBinContent(i + 1, *(fvol + i));
  }
  TF1 *f1 = new TF1("f1", "[0]");
  hwavbaseline->Fit("f1", "NQR");
  baseline = f1->GetParameter(0);
  delete hwavbaseline;
  hwavbaseline = NULL;
  delete f1;
  f1 = NULL;
  return baseline;
}

//============== CALCULATING THE BASELINE: SECOND METHOD ==============//

float Waveform::Calculate_baseline2(float *fvol, int nottrig) {
  double acc = 0;
  TH1D *hnoise = new TH1D("hnoise", "hnoise", 100, -5., 5.);
  TH1D *hbline = new TH1D("hbline", "hbline", 100, -5., 5.);

  bool istrigchan = false;

  for (int i = 0; i < 1000; i++) {
    if (i < 200)
      hnoise->Fill(*(fvol + i));
    hbline->Fill(*(fvol + i));
    if (fabs(*(fvol + i)) > 100.)
      istrigchan = true;

    if ((i >= 30) && (i < 190)) {
      //		SetVol(i, *(fvol+i));
      acc += *(fvol + i);
    }
  }

  double bgmean = 0;
  double meanped = 0;
  baseline = 0;

  TF1 *mgf = new TF1("mgf", "gaus");
  if (nottrig >= 0) {
    hbline->Fit(mgf, "q");
    baseline = mgf->GetParameter(1);
    baseline = hbline->GetMean();
    bnoise = mgf->GetParameter(2);
  } else {
    bnoise = hnoise->GetRMS();
    baseline = acc / 160.;
  }

  for (int i = 0; i < 1000; i++) { // subtracting the baseline
    hpedhist->Fill(*(fvol + i));
  }

  delete mgf;
  delete hnoise;
  delete hbline;
  return baseline;
}

//============== CALCULATING THE BASELINE: THIRD METHOD ==============//
// fit with sine

float Waveform::Calculate_baseline3() {

  TF1 *sinebline = new TF1("sinebline", "[0]+[1]*cos([2]*x+[3])", 0., 100000.);

  sinebline->SetParameter(0, 0.0);
  sinebline->SetParameter(1, 1.0);
  sinebline->SetParameter(2, 0.00055);
  sinebline->SetParameter(3, 0.);
  hwav_raw->Fit("sinebline", "q");
  return hwav_raw //TODO what is the proper return value
}

int Waveform::Calculate_Peaks() {
  npeaks = 0;
  double threshold = TotThreshold * 2;
  double pvol = 0, vollast = 0;
  int nbin = hwav->GetNbinsX();
  int length = 0;
  int MinimumTotBin = (int)(MinimumTot / Deltat);
  for (int i = 0; i < nbin; i++) {
    pvol = TMath::Abs(hwav->GetBinContent(i + 1));
    if (pvol > threshold) {
      length++;
    } else {
      if (length < MinimumTotBin) {
        length = 0;
      } else {
        npeaks++;
        length = 0;
      }
    }
  }
  return npeaks;
}

float Waveform::FWPM(float frac) {

  // cout<<"before "<<endl;

  bool lside = false;
  bool rside = false;
  float fh = frac * amp;
  int mbin = hwav->GetMinimumBin();

  //  cout<<"OK, I'm here..."<<mbin<<" "<<amp<<" "<<fh<<endl;

  int i = 1;
  while (!lside && !rside) {
    if (fabs(hwav->GetBinContent(mbin - i)) < fh)
      lside = true;
    if (fabs(hwav->GetBinContent(mbin + i)) < fh)
      lside = true;
    if (i > 200)
      break;
    //    cout<<i<<endl;
    i++;
  }

  float thewidth = ((float)i) * Deltat;

  //  cout<<"DONE "<<thewidth<<endl;

  return thewidth;
}

float Waveform::Calculate_amp() {
  //	amp = TMath::Abs(fsplinewav->V(hwav->GetMaximumBin()));
  amp = TMath::Abs(hwav->GetMinimum());
  return amp;
}

float Waveform::Calculate_Charge() {
  qtot = TMath::Abs(hwav->Integral() * Deltat / 50); // fC = 1e-15C
  qfast = TMath::Abs(
      hwav->Integral(FitWindow_min / Deltat, FitWindow_max / Deltat) * Deltat /
      50);
  gain = qfast / 1.6e-4; // number of electrons

  // TODO what is the charge
  return 0;
}

float Waveform::Waveform_Filter1(float CutoffFrequency, float finput) {
  double f = 0;
  if (finput < CutoffFrequency)
    f = 1;
  else
    f = 0;
  return f;
}

double Waveform::Calculate_RisingTime() {
  if (amp <= AmpThreshold || npeaks == 0)
    return 0;
  int bin = hwav->GetMinimumBin();
  int peakbin = bin;
  double x0, x1, y0, y1, k;
  double peaktime = hwav->GetBinCenter(peakbin);
  while (hwav->GetBinContent(bin) < -amp * 0.9) {
    bin--;
  }
  x1 = hwav->GetBinCenter(bin);
  x0 = hwav->GetBinCenter(bin - 1);
  y1 = hwav->GetBinContent(bin);
  y0 = hwav->GetBinContent(bin - 1);
  k = (y1 - y0) / (x1 - x0);
  double hightime = (-amp * 0.9 - y0 + k * x0) / k;
  while (hwav->GetBinContent(bin) < -amp * 0.1) {
    bin--;
  }
  x1 = hwav->GetBinCenter(bin);
  x0 = hwav->GetBinCenter(bin - 1);
  y1 = hwav->GetBinContent(bin);
  y0 = hwav->GetBinContent(bin - 1);
  k = (y1 - y0) / (x1 - x0);
  double lowtime = (-amp * 0.1 - y0 + k * x0) / k;
  risingtime = hightime - lowtime;
  return risingtime;
}

float Waveform::Waveform_Filter2(float CutoffFrequency, int T, float finput) {
  double f = 1.0 / (1 + TMath::Power(finput / CutoffFrequency, 2 * T));
  //	TF1 *f1 = new TF1("f1","1/(1+(x/[0])^[1])",0,1e10);
  return f;
}

float Waveform::Filter_KURE(int T, float finput) {

  double f1 = 1.0 / (1 + TMath::Power(finput / 500.0e6, 2 * T));
  double f2 = 1.0 - (1.0 / (1 + TMath::Power(finput / 10.0e6, 2 * T)));

  double f3 = 1.0 / (1 + TMath::Power(finput / 90e6, 2 * T)); // originially
                                                              // 87.0
  double f4 =
      1.0 - (1.0 / (1 + TMath::Power(finput / 80e6, 2 * T))); // originally 90.0

  // return (f1+f2-1);
  return (f1 + f2 - f3 - f4);
}

void Waveform::Waveform_FFT() {
  // Look at the real part of the output
  TH1 *hr = 0;
  hr = hwav->FFT(hr, "RE");
  hr->SetName("hr");
  hr->SetTitle("Real part of the 1st transform");
  hr->GetXaxis()->Set(DimSize, 0, 1e10); // Hz

  /*	//Look at the imaginary part of the output
           TH1 *him = 0;
          him = hwav->FFT(him, "IM");
          him->SetName("him");
          him->SetTitle("Imaginary part of the 1st transform");
          him->GetXaxis()->Set(DimSize,0,1e10);*/

  // Look at the DC component and the Nyquist harmonic:
  double re, im;
  // That's the way to get the current transform object:
  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  // Use the following method to get the full output:
  double *re_full = new double[DimSize];
  double *im_full = new double[DimSize];
  fft->GetPointsComplex(re_full, im_full);
  // filter
  TFile *ff = new TFile("frequencies.root", "RECREATE");
  TH1D *freqq = new TH1D("frequency", "frequency", DimSize, 0, DimSize);
  TH1D *freqq2 = new TH1D("frequency2", "frequency2", DimSize, 0, DimSize);
  for (int i = 0; i < DimSize; i++) {
    //		double f = Waveform_Filter1(5e8, i*1.0e10/DimSize);
    double f = Waveform_Filter2(CutoffFrequency, 4, i * 1.0e10 / DimSize);
    // Filter out 88.5 KURE Ames Alternative
    double f2 = Filter_KURE(200, i * 1.0e10 / DimSize);

    /*
    if( ((i*1.0e10/DimSize)>80.0e6) && ((i*1.0e10/DimSize)<95.0e6) ){
      f=0.0;
    }
    */
    // cout<<f2<<endl;

    re_full[i] = re_full[i] * f2;
    im_full[i] = im_full[i] * f2;

    freqq->SetBinContent(i + 1, re_full[i]);
  }

  // Now let's make a backward transform:
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &DimSize, "C2R M K");
  fft_back->SetPointsComplex(re_full, im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  // Let's look at the output
  hb = TH1::TransformHisto(fft_back, hb, "Re");
  hb->SetTitle("The backward transform result");
  hb->Scale(1.0 / DimSize);
  hb->GetXaxis()->Set(DimSize, 0, DimSize * Deltat); // ps
  for (int i = 0; i < DimSize; i++) {
    hwav->SetBinContent(i + 1, hb->GetBinContent(i + 1));
    vol_fft.push_back(hb->GetBinContent(i + 1));
    re_fft.push_back(hr->GetBinContent(i + 1));
    //	im_fft.push_back(him->GetBinContent(i+1));
    //	vol_fft[i] = hb->GetBinContent(i+1);
    //	re_fft[i] = hr->GetBinContent(i+1);
    //	im_fft[i] = him->GetBinContent(i+1);
    freqq2->SetBinContent(i + 1, abs(re_fft[i]));
  }
  ff->cd();
  freqq->Write();
  freqq2->Write();
  hb->Write();
  hr->Write();

  delete hr;
  //   delete him;
  delete hb;
  delete re_full;
  delete im_full;
  delete fft;
  delete fft_back;
  hr = 0;
  //   him=0;
  hb = 0;
  re_full = 0;
  im_full = 0;
  fft = 0;
  fft_back = 0;
}

void Waveform::Calculate_fitrange() {
  int bin = hwav->GetMinimumBin();
  int peakbin = bin;
  double peaktime = hwav->GetBinCenter(peakbin);
  if (npeaks == 0)
    return;
  while (hwav->GetBinContent(bin) < -amp / 2) {
    bin--;
  }
  double lowtime = hwav->GetBinCenter(bin);
  double fwhm = 2 * (peaktime - lowtime);
  double xmin = peaktime - fwhm;
  double xmax = peaktime + fwhm;
  if (xmin < 0 || xmax > DimSize * Deltat)
    return;
  //	if(3*fwhm<700) return;
  FitWindow_min = peaktime - fwhm;
  FitWindow_max = peaktime + 2 * fwhm;
  //	FitWindow_min = peaktime - 1500;
  //	FitWindow_max = peaktime + 1500;
}

void Waveform::Waveform_Fit() {
  sprintf(SplineName, "SplineFit_ch%d_%d", chno, evtno);
  sprintf(SplineTitle, "SplineFit | SplineFit_ch%d_%d", chno, evtno);
  if (amp <= AmpThreshold || npeaks == 0)
    return;
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();
    fsplinewav = new TSplineFit(SplineName, SplineTitle, 20, PointsPerSpline,
                                hwav, FitWindow_min, FitWindow_max);
    fsplinewav->UpdateFile(
        true); // SplineFit database:  Skip this line for fast analysis
    fsplinewav->ReduceMemory(); // Store only X and Y. Apply this line for fast
                                // analysis
  }
}

float Waveform::CFD_Discriminator1() { // Appoximate CFD
  if (amp <= AmpThreshold) {
    time = 0;
  } // TDC threshold
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();
    Waveform_Fit();
    SetFraction_CFD(0.5);
    double th = -1.0 * Fraction_CFD * amp;
    //		double th = -40;
    double eps = 1e-4;
    int bin = hwav->GetMinimumBin();
    //		while(fsplinewav->V(hwav->GetBinCenter(bin))<=th) {bin--;}
    double xlow = FitWindow_min;
    double xhigh = hwav->GetBinCenter(bin);

    double xmid = (xlow + xhigh) / 2;
    if (fsplinewav->V(xmid) - th == 0)
      time = xmid;

    while ((xhigh - xlow) >= eps) {
      xmid = (xlow + xhigh) / 2;
      if (fsplinewav->V(xmid) - th == 0)
        time = xmid;
      if ((fsplinewav->V(xlow) - th) * (fsplinewav->V(xmid) - th) < 0)
        xhigh = xmid;
      else
        xlow = xmid;
    }
    time = xlow;
  }
  return time;
}

float Waveform::CFD_Discriminator2() { // Exact CFD
  sprintf(CFDSplineName, "CFDSpline_ch%d_%d", chno, evtno);
  sprintf(CFDSplineTitle, "CFDSpline | CFDSpline_ch%d_%d", chno, evtno);
  if (amp <= AmpThreshold || npeaks == 0) {
    time = 0;
  } // TDC threshold
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();
    int Delay_CFDBin = (int)(Delay_CFD / Deltat);
    for (int i = 0; i < DimSize; i++) {
      hcfd->SetBinContent(i + 1, hwav->GetBinContent(i + 1));
    } // copy
    TH1D *hinv = new TH1D("hinv", "hinv", DimSize, 0, DimSize * Deltat);
    for (int i = 0; i < DimSize; i++) {
      if (i + Delay_CFDBin < DimSize)
        hinv->SetBinContent(
            i + 1, hwav->GetBinContent(i + 1 + Delay_CFDBin)); // copy and delay
      else
        hinv->SetBinContent(i + 1, 0);
    }
    hinv->Scale(-1. * Fraction_CFD); // inverse and attenuate
    hcfd->Add(hinv);
    fsplinecfd =
        new TSplineFit(CFDSplineName, CFDSplineTitle, 50, PointsPerSpline, hcfd,
                       FitWindow_min, FitWindow_max);
    //		fsplinecfd->UpdateFile(true); //SplineFit database:  Skip this line
    //for fast analysis
    fsplinecfd->ReduceMemory();

    // solve equation
    double eps = 1e-4;
    double xhigh = hcfd->GetMinimumBin() * Deltat;
    //		double xlow = hcfd->GetMaximumBin()*Deltat;
    double xlow = FitWindow_min;
    //		cout<<xlow<<"\t"<<xhigh<<endl;
    double xmid = (xlow + xhigh) / 2;

    if (fsplinecfd->V(xmid) == 0)
      time = xmid;

    while ((xhigh - xlow) >= eps) {
      xmid = (xlow + xhigh) / 2;
      if (fsplinecfd->V(xmid) == 0)
        time = xmid;
      if (fsplinecfd->V(xlow) * fsplinecfd->V(xmid) < 0)
        xhigh = xmid;
      else
        xlow = xmid;
    }
    time = xlow;
    delete hinv;
    hinv = NULL;
  }
  return time;
}

float Waveform::GaussFit() {
  char gausname[100];
  sprintf(gausname, "gausfit_ch%d_%d", chno, evtno);
  if (amp <= AmpThreshold || npeaks == 0) {
    time = 0;
  } // TDC threshold
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();
    TF1 *mgf = new TF1(gausname, "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", 0,
                       DimSize * Deltat);
    int bin = hwav->GetMinimumBin();
    mgf->SetParameter(0, -amp);
    mgf->SetParameter(1, Deltat * hwav->GetMinimumBin());
    mgf->SetParameter(2, 400);
    mgf->SetParameter(3, 0);
    hwav->Fit(gausname, "QR", "", FitWindow_min, FitWindow_max);
    double par[4];
    mgf->GetParameters(&par[0]);
    double y = -0.5 * amp;
    time = par[2] * sqrt(-2 * log((y - par[3]) / par[0])) + par[1];
    delete mgf;
    mgf = 0;
  }
  return time;
}

// A gaussian fit to the peak, followed by a gaussian fit over the whole
// specified range

Float_t Waveform::DoubleGaussFit(bool floatbaseline) {

  gmean = 0;
  gpeak = 0;
  gsigma = 0;
  gtime = 0;
  gchi2 = 0;

  ggmean = 0;
  ggpeak = 0;
  ggsigma = 0;
  ggmean2 = 0.;
  ggpeak2 = 0;
  ggsigma2 = 0;
  ggtime = 0;
  ggtime2 = 0.;
  ggoffset = 0;
  ggchi2 = 0;

  if (amp <= AmpThreshold || npeaks == 0) {
    time = 0;
  } // TDC threshold
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();

    Double_t amp = TMath::Abs(hwav->GetMinimum());
    int maxbin = hwav->GetMinimumBin();
    double lowrange = hwav->GetBinCenter(maxbin - 4);
    double hirange = hwav->GetBinCenter(maxbin + 4);

    if ((hwav->GetBinCenter(maxbin) < FitWindow_max) &&
        (hwav->GetBinCenter(maxbin) > FitWindow_min)) {
      TF1 *mgf = new TF1("mgf", "gaus", lowrange, hirange);
      mgf->SetParameter(0, -amp);
      mgf->SetParameter(1, Deltat * hwav->GetMinimumBin());
      mgf->SetParameter(2, 400);

      hwav->Fit("mgf", "QR");

      double peaktime = mgf->GetParameter(1);

      double flowrange = peaktime - GaussRange_min;
      double fhirange = peaktime + GaussRange_max;

      if (peaktime < FitWindow_max && peaktime > FitWindow_min) {

        TF1 *mgf2;
        if (!floatbaseline)
          mgf2 = new TF1("mgf2", "gaus", flowrange, fhirange);
        else
          mgf2 = new TF1("mgf2", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", flowrange,
                         fhirange);

        double nparam = 4;
        if (!floatbaseline)
          nparam = 3;

        mgf2->SetParameter(0, -amp);
        mgf2->SetParameter(1, Deltat * hwav->GetMinimumBin());
        mgf2->SetParameter(2, 400);
        if (floatbaseline)
          mgf2->SetParameter(3, 0);

        hwav->Fit("mgf2", "QR");
        double mgpeak = mgf2->GetParameter(0);
        double mgmean = mgf2->GetParameter(1);
        double mgsigma = mgf2->GetParameter(2);
        double mgoffset;
        if (floatbaseline)
          mgoffset = mgf2->GetParameter(3);
        else
          mgoffset = 0;
        double mgchi2 = mgf2->GetChisquare();

        gmean = mgmean;
        gpeak = mgpeak;
        gsigma = mgsigma;
        //		    gtime=mgmean-sqrt(-mgsigma*mgsigma*log(0.37));
        gtime = mgmean - sqrt(-mgsigma * mgsigma * log(0.1));
        goffset = mgoffset;
        gchi2 = mgchi2;
        gdegfree =
            (hwav->FindBin(fhirange)) - (hwav->FindBin(flowrange)) - nparam;

        delete mgf;
        delete mgf2;
      }
    }
  }

  return time;
}

// Not really used - a gaussian fit to the peak, followed by a double-gaussian
// fit over the whole specified range

Float_t Waveform::DoubleDoubleGaussFit() {

  ggmean = 0;
  ggpeak = 0;
  ggsigma = 0;
  ggmean2 = 0.;
  ggpeak2 = 0;
  ggsigma2 = 0;
  ggtime = 0;
  ggtime2 = 0.;
  ggoffset = 0;
  ggchi2 = 0;

  char gausname[100];
  sprintf(gausname, "gausfit_ch%d_%d", chno, evtno);
  if (amp <= AmpThreshold || npeaks == 0) {
    time = 0;
  } // TDC threshold
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();

    Double_t amp = TMath::Abs(hwav->GetMinimum());
    int maxbin = hwav->GetMinimumBin();
    double lowrange = hwav->GetBinCenter(maxbin - 6);
    double hirange = hwav->GetBinCenter(maxbin + 6);

    if ((hwav->GetBinCenter(maxbin) < FitWindow_max) &&
        (hwav->GetBinCenter(maxbin) > FitWindow_min)) {
      TF1 *mgf = new TF1("mgf", "gaus", lowrange, hirange);
      mgf->SetParameter(0, -amp);
      mgf->SetParameter(1, Deltat * hwav->GetMinimumBin());
      mgf->SetParameter(2, 400);

      hwav->Fit("mgf", "QR");

      double peaktime = mgf->GetParameter(1);

      double flowrange = peaktime - 2650.;
      double fhirange = peaktime + 400.;

      if (peaktime < FitWindow_max && peaktime > FitWindow_min) {

        TF1 *mgf2 = new TF1(
            "mgf2",
            "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+[6]",
            flowrange, fhirange);

        mgf2->SetParameter(0, -amp);
        mgf2->SetParameter(1, Deltat * hwav->GetMinimumBin());
        mgf2->SetParameter(2, 400);
        mgf2->SetParameter(3, -amp);
        mgf2->SetParameter(4, Deltat * hwav->GetMinimumBin() + 100);
        mgf2->SetParameter(5, 400);
        mgf2->SetParameter(6, 0);

        hwav->Fit("mgf2", "QR");
        double mgpeak = mgf2->GetParameter(0);
        double mgmean = mgf2->GetParameter(1);
        double mgsigma = mgf2->GetParameter(2);

        double mgpeak2 = mgf2->GetParameter(3);
        double mgmean2 = mgf2->GetParameter(4);
        double mgsigma2 = mgf2->GetParameter(5);
        double mgoffset = mgf2->GetParameter(6);
        double mgchi2 = mgf2->GetChisquare();
        double nparam = 7;

        ggmean = mgmean;
        ggpeak = mgpeak;
        ggsigma = mgsigma;
        ggmean2 = mgmean2;
        ggpeak2 = mgpeak2;
        ggsigma2 = mgsigma2;
        ggtime = mgmean - sqrt(-mgsigma * mgsigma * log(0.03));
        ggtime2 = mgmean2 - sqrt(-mgsigma2 * mgsigma2 * log(0.03));
        ggoffset = mgoffset;
        ggchi2 = mgchi2;
        ggdegfree = ((double)hwav->FindBin(fhirange)) -
                    ((double)hwav->FindBin(flowrange)) - nparam;

        delete mgf;
        delete mgf2;
      } // if seed fit in final time window
        //	Cout<<" "<<Mgmean<<Endl=gtime;
    }   // if raw peak is in time window
  }     // else (more than zero pulse candidates
  return time;
}

Float_t Waveform::DoTemplateFit(TemplateFit *mtf) {

  ttime = 0;
  tscale = 0;
  tchi2 = 0;

  if (amp <= AmpThreshold || npeaks == 0) {
    time = 0;
  } // TDC threshold
  //	if(amp <= AmpThreshold) { time = 0; } //TDC threshold
  //	if(npeaks==0) { time = 0; } //TDC threshold
  else {
    if (IfDynamicWindow == 1)
      Calculate_fitrange();

    Double_t amp = TMath::Abs(hwav->GetMinimum());
    int maxbin = hwav->GetMinimumBin();

    // Attempt to follow the Gaussfit and set range dynamically (not currently
    // used?)
    /*
    double lowrange = hwav->GetBinCenter(maxbin-6);
    double hirange = hwav->GetBinCenter(maxbin+6);

    TF1* mgf = new TF1("mgf","gaus",lowrange,hirange);
    mgf->SetParameter(0, -amp);
    mgf->SetParameter(1, Deltat*hwav->GetMinimumBin());
    mgf->SetParameter(2, 400);

    hwav->Fit("mgf","QR");

    double peaktime = mgf->GetParameter(1);

    double flowrange = peaktime-1100.;
    double fhirange = peaktime+250.;

    mtf->SetTempRange(flowrange,fhirange);
    mtf->SetFitRange(flowrange,fhirange);
    */

    /////////////////////////////////////////////////////////

    if ((hwav->GetBinCenter(maxbin) < FitWindow_max) &&
        (hwav->GetBinCenter(maxbin) > FitWindow_min)) {

      TH1D *gthetemplate = mtf->ReturnTemplate(true);
      mtf->FitWithTemplate(hwav, ttime, tscale, tchi2);
      tamp = tscale * (fabs(gthetemplate->GetMinimum()));

      double binwidth;
      binwidth = gthetemplate->GetBinWidth(1);

      tcharge =
          ((-tscale) * binwidth * (gthetemplate->Integral()) * 6241.5) / 50.;
      tNS = 0;
      if (tamp > 0)
        tNS = bnoise / tamp;
      tdegfree = 0.;
    }
    time = ttime;
  }
  return time;
}

int Waveform::Calculate_Peaks_nnls() {
  nnpeaks = 0;

  HighBound.clear();
  LowBound.clear();

  double threshold = 2; // based on best guess after looking at xvector
  double vol = 0;
  int nbin = (int)m; // bins in xvector
  int length = 0;
  int MinimumTotBin = 10; // minimum length of peak to characterize a peak
                          // (again based on guess)
  for (int i = 30 * steps; i < nbin - 30 * steps; i++) {
    vol = TMath::Abs(xvector->GetBinContent(i + 1));
    if (vol > threshold) {
      length++;
    } else {
      if (length < MinimumTotBin) {
        length = 0;
      } else {
        nnpeaks++;
        TH1D *NEWhwav = (TH1D *)xvector->Clone();
        NEWhwav->GetXaxis()->SetRange((i + 1 - length), (i + 1));
        Int_t MaxBin = NEWhwav->GetMaximumBin() -
                       30 * steps; // subtracting the extra 30 buffer bins (we
                                   // started the loop at 30*steps)
        MaxBin = MaxBin + 12 * steps; // difference in bins between the
                                      // beginning of the template and the peak
                                      // on the template
        cout << "peaktime: " << MaxBin * Deltat / steps << endl;
        LowBound.push_back(i + 1 - length - 30 * steps + 10 * steps);
        HighBound.push_back(i + 1 - 17 * steps);

        if (nnpeaks <= 10)
          vtime[nnpeaks - 1] = MaxBin * Deltat / steps;
        // old way of setting boundaries
        /*
        if (MaxBin+8*steps<nbin)
          {
            HighBound.push_back(MaxBin+8*steps);
            vHiB[nnpeaks]=MaxBin+8;
          }
        else
          {
            HighBound.push_back(nbin);
            vHiB[nnpeaks]=nbin;
          }
        if (MaxBin-8*steps>0)
          {
            LowBound.push_back(MaxBin-8*steps);
            vLowB[nnpeaks]=MaxBin-8*steps;
          }
        else
          {
            LowBound.push_back(0);
            vLowB[nnpeaks]=0;
          }
        */
        delete NEWhwav;
        length = 0;
      }
    }
  }

  for (int i = 0; i < (int) HighBound.size(); i++) {
    cout << "Bounds: " << LowBound[i] << " " << HighBound[i] << endl;
  }

  return nnpeaks;
}

void Waveform::Calculate_Variables_nnls(int Npulses) {

  int pulses;
  if (Npulses <= 10) {
    pulses = Npulses;
  } else {
    pulses = 10;
  }

  // BEGIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TH1D *REBINNEDhwav =
      new TH1D("rebinned", "rebinned", DimSize * steps, 0, DimSize * Deltat);
  for (int i = 0; i < DimSize * steps; i++) {
    REBINNEDhwav->SetBinContent(i + 1,
      hwav_raw->GetBinContent((i + 1) / (steps - 0.01)));
  }

  for (int p = 0; p < pulses; p++) {
    // CALCULATING CHI2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vchi2[p] = 0;
    if (trigno >= 0) {
      for (int i = 0; i < (int) m - 60 * steps; i++) {
        if (LowBound[p] < i && HighBound[p] > i) {
          if (REBINNEDhwav->GetBinContent(i + 1) != 0) {
            vchi2[p] += TMath::Abs(((nnlsoutput->GetBinContent(i + 1) -
                                     (REBINNEDhwav->GetBinContent(i + 1))) *
                                    (nnlsoutput->GetBinContent(i + 1) -
                                     (REBINNEDhwav->GetBinContent(i + 1))) /
                                    (REBINNEDhwav->GetBinContent(i + 1))));
          }
        }
      }
      cout << "chi2 for pulse " << p + 1 << ": " << vchi2[p] << endl;
    }

    // CALCULATING TIMING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // done in calculate_peaks_nnls

    // CALCULATING CHARGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vcharge[p] = 0;
    if (trigno >= 0) {
      vcharge[p] = TMath::Abs(nnlsoutput->Integral(LowBound[p], HighBound[p]));
      cout << "Charge for pulse " << p + 1 << ": " << vcharge[p] << endl;
    }

    // CALCULATING AMP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vamp[p] = 0;
    if (trigno >= 0) {
      double min = 0;
      for (int i = LowBound[p]; i < HighBound[p]; i++) {
        if (nnlsoutput->GetBinContent(i) < min)
          min = nnlsoutput->GetBinContent(i);
      }
      vamp[p] = TMath::Abs(min);
      cout << vamp[p] << endl;
    }

    // CALCULATING FWPM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (trigno >= 0) {
      double frac = 0.2; // choose what fraction you want FWPM to be at

      int MinBin = steps * vtime[p] / Deltat;
      bool lside = false;
      bool rside = false;
      float fh = frac * vamp[p];

      int i = 0;
      int j = 0;

      while (!lside) {
        i++;
        if (fabs(nnlsoutput->GetBinContent(MinBin - i)) < fh)
          lside = true;
      }

      while (!rside) {
        j++;
        if (fabs(nnlsoutput->GetBinContent(MinBin + j)) < fh)
          rside = true;
      }

      float slopeL =
          (nnlsoutput->GetBinContent(MinBin - i + 1) -
           nnlsoutput->GetBinContent(MinBin - i)) /
          ((MinBin - i + 1) * Deltat / steps - (MinBin - i) * Deltat / steps);
      float slopeR =
          (nnlsoutput->GetBinContent(MinBin + j) -
           nnlsoutput->GetBinContent(MinBin + j - 1)) /
          ((MinBin + j) * Deltat / steps - (MinBin + j - 1) * Deltat / steps);

      float interceptL = nnlsoutput->GetBinContent(MinBin - i) -
                         slopeL * (MinBin - i) * Deltat / steps;
      float interceptR = nnlsoutput->GetBinContent(MinBin + j) -
                         slopeR * (MinBin + j) * Deltat / steps;

      float timeL = (fh - interceptL) / slopeL;
      float timeR = (fh - interceptR) / slopeR;

      float thewidth = timeR - timeL;

      vFWPM[p] = thewidth;
      cout << "fwpm for pulse " << p + 1 << ": " << vFWPM[p] << endl;
    }
  }
  delete REBINNEDhwav;
}

void Waveform::Analyze() {
  DisWindow_min = 0;
  DisWindow_max = DimSize * Deltat;
  hwav->GetXaxis()->SetRange(DisWindow_min / Deltat, DisWindow_max / Deltat);
  hcfd->GetXaxis()->SetRange(DisWindow_min / Deltat, DisWindow_max / Deltat);

  bool nnls = true;
  if (nnls) {
    cout << "using nnls algorithm" << endl;
    if (trigno) {
      Setup_nnls();
      int pulses = Calculate_Peaks_nnls();
      Calculate_Variables_nnls(pulses);
    }
  } else {
    if (DoFFT == 1)
      Waveform_FFT();
    Calculate_Peaks();
    Calculate_amp();
    FWHM = FWPM(0.5);
    FW20pM = FWPM(0.2);
    Calculate_Charge();
    //	Calculate_RisingTime();
    //	CFD_Discriminator1();	//approximate CFD
    CFD_Discriminator2(); // exact CFD
    DoubleGaussFit(0);
    //      DoubleDoubleGaussFit();
  }
}

// Overloaded version of the Analyze function that includes the template fit

void Waveform::Analyze(TemplateFit *mtf) {
  DisWindow_min = 0;
  DisWindow_max = DimSize * Deltat;
  hwav->GetXaxis()->SetRange(DisWindow_min / Deltat, DisWindow_max / Deltat);
  hcfd->GetXaxis()->SetRange(DisWindow_min / Deltat, DisWindow_max / Deltat);

  if (DoFFT == 1)
    Waveform_FFT();
  Calculate_Peaks();
  Calculate_amp();
  FWHM = FWPM(0.5);
  FW20pM = FWPM(0.2);
  Calculate_Charge();
  //	Calculate_RisingTime();
  //	CFD_Discriminator1();	//approximate CFD
  CFD_Discriminator2(); // exact CFD
  DoubleGaussFit(0);
  DoTemplateFit(mtf);
  //      DoubleDoubleGaussFit();
}
