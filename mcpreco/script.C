
script()
{
    TString filename = "p7_3200_1650";
    TString filename = filename + ".root";
    TString histname = "/peds/PED_ch0_1361";

    cout << filename << histname << endl;


    TFile* problematic_hist = new TFile(filename);
    TH1D* thehist = (TH1D*)problematic_hist->Get(histname);

    TCanvas* canvas = new TCanvas("c1","c1",600,450);
    thehist->Draw();

    //TF1* fit = new TF1("thefit","gaus");
    Deltat = 100;
    DimSize = 256;
    TF1* thefit = new TF1("thefit","[0]*exp(-0.5*((x-[1])/[2])^2)", 0, DimSize*Deltat);
    int mbin = thehist->GetMaximumBin();
    thefit->SetParameter(0,thehist->GetBinContent(mbin));
    thefit->SetParameter(1,thehist->GetMean());
    thefit->SetParameter(2,0.8);

    thehist->Fit(thefit);
    cout<<thefit->GetChisquare()<<endl;
    
}