//----------APS Analysis
//----------Author's Name:Jingbo WANG
//----------Copyright:Those valid for ANL
//----------Modified:22/07/2014
#include <RVersion.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TProcessID.h>
#include <TMath.h>
#include "UVevent.h"
#include "Waveform.h"
#include <TSpectrum.h>

#include <ctime>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <cstdarg>

#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>

#include "nnls.h"

#include <iostream>

using namespace std;

//using namespace nsNNLS;

nsNNLS::matrix* load (const char* fname, char type);
nsNNLS::vector* load_vector(const char* fil, size_t sz, char dos);
int     write_vector(const char* fil, nsNNLS::vector* v);
int     write_binvector(const char* fil, nsNNLS::vector* v);

ClassImp(UVevent)
ClassImp(Waveform)

TClonesArray *UVevent::fgWaveforms = 0;

// ------------------------------------------------------
// UVevent
// ------------------------------------------------------
UVevent::UVevent()
{
  // Create a UVevent object.
   // When the constructor is invoked for the first time, the class static
   // variable fgTracks is 0 and the TClonesArray fgTracks is created.
	if (!fgWaveforms) fgWaveforms = new TClonesArray("Waveform", 1000);
	fwav = fgWaveforms;
	fnwav = 0;
	WavDimSize = 1032;
	WavSampleRate = 1.e10;
	WavDeltat = 100;
}

UVevent::~UVevent() {
	Clear();
}

void UVevent::Initialize() {
	for(int i=0;i<WavDimSize;i++) {
		WavDeltat = 1.e12/WavSampleRate;
		t.push_back(i*1.e12/WavSampleRate);
		fz.push_back(i*WavSampleRate/WavDimSize); //Hz
	}
}

void UVevent::UVeventAnalysis()
{	}

void UVevent::Clear(Option_t * /*option*/)
{
   fwav->Clear("C"); //will also call Waveform::Clear
   fnwav = 0;
   t.clear();
   fz.clear();
//	for(int i=0;i<60;i++) {wav[i].Clear();}

//   delete fwav;

}

// ------------------------------------------------------
// Waveform
// ------------------------------------------------------
