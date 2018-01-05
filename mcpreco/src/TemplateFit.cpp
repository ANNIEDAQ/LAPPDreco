#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"
#include "TMath.h"
#include "TemplateFit.h"
#include "TMinuit.h"
#include "TGraph.h"

using namespace std;

ClassImp(TemplateFit)
//ClassImp(Reco)

//TemplateFit fgTemplateFit("extern");

TH1D* ghist;
TH1D* gtemplate;
bool gusefitrange;
bool gusetemprange;
double gtlrange;
double gthrange;
double glrange;
double ghrange;
double gnoise;
double gtstart;
double gscalestart;
double gchi2;
double gdegfree;
int gdegfree_fixed;
TSplineFit *tshist;


extern TH1D* XShiftAndScale(TH1D* ohist, double mshift, double mscale){

  TH1D* rhist = (TH1D*)(ohist->Clone());
  // Double_t origNORM = rhist->Integral();

  rhist->Reset();

  for(int i=1; i<((ohist->GetNbinsX())+1); i++){

    double bcent = ohist->GetBinCenter(i);
    double bcont = ohist->Interpolate(bcent-mshift);

//    cout<<bcont<<" "<<mscale<<" "<<(bcent-mshift)<<" "<<gtlrange<<" "<<gthrange<<endl;
    rhist->SetBinContent(i,bcont*mscale);


    if( gusetemprange && (((bcent-mshift)<gtlrange) || ((bcent-mshift)>gthrange)) ){
        rhist->SetBinContent(i,0);
     }
    }

  // Double_t newNORM = rhist->Integral();

  //  cout<<"NORMS "<<origNORM<<" "<<newNORM<<endl;
  //  if(newNORM>0) rhist->Scale(origNORM/newNORM);

  return rhist;
}


extern TH1D* XShiftAndScale(TH1D* ohist, double mshift, double mscale, bool dorange){

  TH1D* rhist = (TH1D*)(ohist->Clone());
  // Double_t origNORM = rhist->Integral(50,950);

  rhist->Reset();

  for(int i=1; i<((ohist->GetNbinsX())+1); i++){

    double bcent = ohist->GetBinCenter(i);
    double bcont = ohist->Interpolate(bcent-mshift);
    // double bcont = ghist->Eval((bcent-mshift),0,"S");


    rhist->SetBinContent(i,bcont*mscale);
  }

  Double_t newNORM = rhist->Integral(50,950);

  //  if(fabs(mshift)<0.00001) cout<<"WWWW  "<<origNORM<<" "<<newNORM<<endl;;
  //  cout<<"POEU "<<fabs(mshift)<<" "<<origNORM-newNORM<<" "<<(ohist->GetNbinsX())<<endl;

  //  cout<<"NORMS "<<origNORM<<" "<<newNORM<<endl;
  //  if(newNORM>0) rhist->Scale(origNORM/newNORM);

  return rhist;
}

extern double Chi2WithTemplate(TH1D* fhist,TH1D* ftemp, double noise){

  int nbfh = fhist->GetNbinsX();
  int nbt = ftemp->GetNbinsX();

  double chi2=0.0;
  int degfree=0;

  double chi2sum=0;
  double mNorm = fhist->Integral();
  double tNorm = ftemp->Integral();

  if(tNorm!=1) // cout<<"TNORM NOT = TO 1.0!! "<<tNorm<<endl;
  if(tNorm==0) tNorm=1;

  for(int i=1; i<nbfh+1; i++){

      // content of data hist
      double fbc = fhist->GetBinContent(i);
      // bin center of data hist
      double bcent = fhist->GetBinCenter(i);
      if(noise<1.0) noise=1.0;

      //      double tbc = (mNorm/tNorm)*(ftemp->GetBinContent(i));
      //double tbc = (ftemp->GetBinContent(i));

      // find which bin in the template corresponds to the bin center in the data
      int wbin = ftemp->FindBin(bcent);
      double tempval=0;

      // as long as the template value is not 0 or next to a point where the template is cutoff, interpolate
      if( (ftemp->GetBinContent(wbin)!=0) && (ftemp->GetBinContent(wbin-1)!=0) && (ftemp->GetBinContent(wbin+1)!=0) ){
     	 tempval = ftemp->Interpolate(bcent);
      }
      // if the template is after the cutoff, set tempval to 0
      if(ftemp->GetBinContent(wbin)==0) tempval=0;
      // if the template is next to an edge, just use the bin contents and DON'T interpolate
      if( (ftemp->GetBinContent(wbin-1)==0) || (ftemp->GetBinContent(wbin+1)==0) ) tempval=ftemp->GetBinContent(wbin);

      double mchi = (fbc-tempval)/noise;

      if( (!gusefitrange) || (tempval!=0) ){
	chi2sum+= (mchi*mchi);
	degfree++;
      }
  }


    //    chi2 = sqrt(chi2sum/((double)nbt));
  chi2 = chi2sum;

  // if(degfree!=15) cout<<"deg free not equal to 15: "<<degfree<<endl;

  gchi2=chi2;
  gdegfree = (double)degfree;

  // cout<<"THE CHI "<<chi2<<" "<<gdegfree<<" "<<(chi2/gdegfree)<<" "<<gusefitrange<<endl;

  return chi2/(gdegfree+3);
}


// NEW VERSION OF THIS...RATHER THAN SHIFT THE TEMPLATE, WE JUST EVALUATE IT AT THE SHIFTED VALUE ///////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

extern double Chi2WithTemplate(TH1D* fhist,TH1D* ftemp, double mshift, double mscale, double noise){

  int nbfh = fhist->GetNbinsX();
  int nbt = ftemp->GetNbinsX();

  double chi2=0.0;
  int degfree=0;

  double chi2sum=0;
  double mNorm = fhist->Integral();
  double tNorm = ftemp->Integral();

  if(tNorm!=1) // cout<<"TNORM NOT = TO 1.0!! "<<tNorm<<endl;
  if(tNorm==0) tNorm=1;

  for(int i=1; i<nbfh+1; i++){

      // content of data hist
      double fbc = fhist->GetBinContent(i);
      // bin center of data hist
      double bcent = fhist->GetBinCenter(i);
      if(noise<1.0) noise=1.0;

      //      double tbc = (mNorm/tNorm)*(ftemp->GetBinContent(i));
      //double tbc = (ftemp->GetBinContent(i));

      // find which bin in the template corresponds to the bin center in the data
      int wbin = ftemp->FindBin(bcent-mshift);
      double tempval=0;

      // as long as the template value is not 0 or next to a point where the template is cutoff, interpolate

      double shiftx = bcent-mshift;


      bool inrange=false;
      if( (shiftx>gtlrange) &&  (shiftx<=gthrange) ){
	 tempval = mscale*(ftemp->Interpolate(shiftx));
	 inrange=true;
      }

      double mchi = (fbc-tempval)/noise;

      if( (!gusefitrange) || (inrange) ){
	chi2sum+= (mchi*mchi);
	degfree++;
      }
  }

/*
      bool isclose=false;
      if(fabs(mshift)>40000) isclose=true;
*/
  //    chi2 = sqrt(chi2sum/((double)nbt));
  chi2 = chi2sum;

  // cout<<degfree<<" "<<gdegfree_fixed<<" "<<nbfh<<endl;

 // if((degfree!=0) && (degfree!=15) ) cout<<"wtf "<<degfree<<" "<<mshift<<" "<<mscale<<" "<<endl;
 // if(degfree!=50) {chi2=555; degfree=1;}
 //  if( fabs(fabs(mshift)-50)<1) cout<<" tt "<<chi2<<" "<<degfree<<" "<<mshift<<" "<<mscale<<" "<<endl;


  gchi2=chi2;
  gdegfree = (double)degfree;

  // cout<<"THE CHI "<<chi2<<" "<<gdegfree<<" "<<(chi2/gdegfree)<<" "<<gusefitrange<<endl;

//  if((chi2/(gdegfree+3))<1.5) cout<<mshift<<" "<<mscale<<endl;

  return chi2/(gdegfree+3);
}





extern void fit_template_fcn(Int_t&, Double_t*, Double_t &f, Double_t* par, Int_t)
{

  double tfoffset = par[0];
  double tfscale = par[1];
  //  double tfscale = 0.001;
  double tfnoise = gnoise;

  double maxtemp = gtemplate->GetMaximum();
  double maxhist = ghist->GetMaximum();

  if(maxhist>0){
    //    tfscale = maxtemp/maxhist;
    //   cout<<"SCALE "<<tfscale<<endl;
  }

  //  TH1D* gtempadj = XShiftAndScale(gtemplate,tfoffset,tfscale);
  //  f = Chi2WithTemplate(ghistadj,gtempadj,tfnoise);
  //  f = Chi2WithTemplate(ghist,gtempadj,tfnoise);
  f = Chi2WithTemplate(ghist,gtemplate,tfoffset,tfscale,tfnoise);



//  delete gtempadj;

}


Double_t TemplateFit::CalculateFCN(Double_t &f, Double_t *par){

  int myi;
  Double_t *t = new Double_t[2];
  t[0]=1.0;
  t[1]=1.0;
  Double_t mfcn;

  fit_template_fcn(myi,t,mfcn,par,3);

  delete t;

  f=mfcn;
  return mfcn;
}


TemplateFit::TemplateFit()
{
	TemplateFit("derp");
}

TemplateFit::TemplateFit(TString tname)
{
  // Create an Event object.
   // When the constructor is invoked for the first time, the class static
   // variable fgTracks is 0 and the TClonesArray fgTracks is created.

	/*
	if (!fgRecos) fgRecos = new TClonesArray("Reco", 1000);
	fRecos = fgRecos;
	*/
  TString fgname;
  fgname+=tname;
  fgname+="fixed";

  gtstart=0.;
  gscalestart=400.;
  gnoise=10.;

  /*
  cout<<fgTemplateFit<<endl;
  if(!fgTemplateFit){
    cout<<"derp"<<endl;
    fgTemplateFit = new TemplateFit(fgname);
  }
  */


  tempcount=0;

  int _nb=1032;
  int _lb=0.0;
  int _hb=1032*100;

  first_template=true;
  TString fulltname;
  fulltname+="template_";
  fulltname+=tname;
  _ptemplate = new TH1D(fulltname,fulltname,_nb,_lb,_hb);

  gdegfree_fixed=0;

  fMinuitScaleAndOffset = new TMinuit(1);
  fMinuitScaleAndOffset->SetPrintLevel(-1);
  fMinuitScaleAndOffset->SetMaxIterations(500);
}


TemplateFit::~TemplateFit() {
	Clear();
}

bool TemplateFit::IsFirstTemp()
{
  return first_template;
}

void TemplateFit::SetFitRange(Double_t lr, Double_t hr){

  gusefitrange = true;
  glrange = lr;
  ghrange = hr;

}

void TemplateFit::CalcFixedDegFree(TH1D* ttemp)
{
	gdegfree_fixed = 0;
	for(int i=0; i<ttemp->GetNbinsX(); i++){
		double bcent = ttemp->GetBinCenter(i);
		if( (bcent>glrange) && (bcent<ghrange) ){
			gdegfree_fixed++;
		}
	}
//	cout<<"NUM FIXED DEGREES FREE! "<<gdegfree_fixed<<" "<<ttemp->GetNbinsX()<<endl;
}

void TemplateFit::SetTempRange(Double_t lr, Double_t hr){

  gusetemprange = true;
  gtlrange = lr;
  gthrange = hr;

}

void TemplateFit::AddTemplate(TH1D* ntemp){

  if(first_template){
    cout<<"first template"<<endl;
    _ptemplate = (TH1D*) (ntemp->Clone());
    first_template=false;
  } else{
    //    cout<<"not first template"<<endl;
    _ptemplate->Add(ntemp);
  }

  tempcount++;
}

void TemplateFit::AddTemplateNoNorm(TH1D* ntemp){

  _ptemplate->Add(ntemp);

  double pNORM = _ptemplate->Integral();
  if(pNORM>0) _ptemplate->Scale(1/pNORM);

  tempcount++;
}


void TemplateFit::AddSmoothedTemplate(TH1D* sthist){

  this->AddTemplate(sthist);

  if(gtemplate!=0) gtemplate->Clear();

  double fwhm=4000;
  int bin = sthist->GetMinimumBin();
  int peakbin = bin;
  double peaktime = sthist->GetBinCenter(peakbin);
  double FitWindow_min = peaktime - fwhm;
  double FitWindow_max = peaktime + fwhm;
  int 	PointsPerSpline = 4;
  TSplineFit* fsplinehist = new TSplineFit("thespline","thespline",50,PointsPerSpline, sthist,FitWindow_min, FitWindow_max);


  //interpolate template to a finer binning
  int Tnbins = sthist->GetNbinsX();
  double minval=sthist->GetXaxis()->GetBinLowEdge(1);
  double maxval=sthist->GetXaxis()->GetBinUpEdge(Tnbins);
  int numbins = (int)(maxval-minval)/10;
  TH1D* newthist = new TH1D("rebinnedtemphist","rebinnedtemphist",numbins,minval,maxval);


  for(int i=1; i<numbins+1; i++){

   double theT = newthist->GetBinCenter(i);
   double bcont;

   if( (theT>FitWindow_min) && (theT<FitWindow_max) ) bcont = fsplinehist->V(theT);
   else bcont=0;

/*
   if( (theT>(sthist->GetBinCenter(1))) && (theT<(sthist->GetBinCenter(Tnbins))))
   {
   bcont = sthist->Interpolate(theT);
   }
   if(theT<(sthist->GetBinCenter(1))) {bcont = (sthist->GetBinContent(1));}
   if(theT>(sthist->GetBinCenter(Tnbins))) {bcont = (sthist->GetBinContent(Tnbins));}
*/

   newthist->SetBinContent(i,bcont);
  }

 // this->CalcFixedDegFree(newthist);

  // use newly derrived template
  gtemplate = newthist;

  // gtemplate = sthist;
}

void TemplateFit::SetTemplate(TH1D* ntemp){

  _ptemplate = ntemp;

}


TH1D* TemplateFit::ShiftAndScale(TH1D* ohist, double mshift, double mscale){

  TH1D* rhist = XShiftAndScale(ohist,mshift,mscale,0);

  return rhist;
}


void TemplateFit::Combine(TH1D* ohist, double mshift, double mscale){

  TH1D* chist = this->ShiftAndScale(ohist,mshift,mscale);
  this->AddTemplateNoNorm(chist);
  delete chist;
}


void TemplateFit::SetFitHists(TH1D* sfhist, TH1D* sthist){

  gtemplate = sthist;
  ghist = sfhist;
}


void TemplateFit::SetNoise(double tnoise){

  gnoise = tnoise;
}

void TemplateFit::SetStartingValues(double tstart, double scalestart){

  gtstart = tstart;
  gscalestart = scalestart;
}


void TemplateFit::FitWithTemplate(TH1D* fhist, double &moffset, double &mscale, double &mchi2){

  /*
  fgTemplateFit.SetTemplate(thist);
  fgTemplateFit.SetPulseHist(fhist);
  */

  Double_t tNORM = gtemplate->Integral();
  Double_t fNORM = fhist->Integral();

  ghist = fhist;

  Double_t seedshift = gtstart;
  Double_t seedscale = gscalestart;
  if(tNORM!=0) seedscale = fNORM/tNORM;
  Double_t shiftstep = 1.0;
  Double_t scalestep = 50.0;

  Double_t *arglist = new Double_t[10];
  arglist[0]=1;  //1: standard minimization
                 //2: try to improve minimum
  Int_t err = 0;
  Int_t flag=0;

  fMinuitScaleAndOffset->SetFCN(fit_template_fcn);
  fMinuitScaleAndOffset->mnexcm("SET STR",arglist,1,err);
  fMinuitScaleAndOffset->mnparm(0,"shift",seedshift,shiftstep,0,0,err);
  fMinuitScaleAndOffset->mnparm(1,"scale",seedscale,scalestep,0,0,err);

  flag = fMinuitScaleAndOffset->Migrad();

  Double_t toffset,toffseterr,tscale,tscaleerr;
  fMinuitScaleAndOffset->GetParameter(0,toffset,toffseterr);
  fMinuitScaleAndOffset->GetParameter(1,tscale,tscaleerr);

  delete [] arglist;

  moffset=toffset;
  mscale=0.001;

  double maxtemp = gtemplate->GetMaximum();
  double maxhist = ghist->GetMaximum();

  if(maxhist>0){
  //  mscale = maxtemp/maxhist;
  }

  //  TH1D* gtempadj = XShiftAndScale(thist,moffset,mscale);
  //  mchi2 = Chi2WithTemplate(fhist,gtempadj,10.0);
  mchi2 = gchi2/(gdegfree+3);
  // delete gtempadj;

  mscale=tscale;
}


void TemplateFit::FitWithTemplate(TH1D* fhist, TH1D* thist, double &moffset, double &mscale, double &mchi2){

  this->SetFitHists(fhist,thist);
  this->FitWithTemplate(fhist, moffset, mscale, mchi2);
}

TH1D* TemplateFit::ReturnTemplate(bool NormIt){

  //  cout<<"number of histos in template: "<<tempcount<<endl;
  _ptemplate->Scale(1/((Double_t) tempcount));

  if(NormIt){
     Double_t tNORM = fabs(_ptemplate->Integral());
     if(tNORM>0) _ptemplate->Scale(1/tNORM);
  }

  return _ptemplate;

}

TH1D* TemplateFit::ReturnGTemplate(){

  return gtemplate;

}


TH1D* TemplateFit::ReturnPulseHist(){

  cout<<"number of histos in template: "<<tempcount<<endl;

  return _phist;
}

void TemplateFit::SetPulseHist(TH1D* phist){

  _phist=phist;

}


void TemplateFit::Clear(Option_t * /*option*/)
{
//   fWaveforms->Clear("C"); //will also call Track::Clear
//	for(Int_t i=0;i<4;i++) {wav[i].Clear();}
//   delete fWaveforms;

}
