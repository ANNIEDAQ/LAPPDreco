void plottingvoltages()
{
  TFile* f = new TFile("VoltagePlots.root","RECREATE");

  TCanvas* ct1 = new TCanvas("ct1","ct1",600,600);
  TCanvas* ct2 = new TCanvas("ct2","ct2",600,600);
  TCanvas* ct3 = new TCanvas("ct3","ct3",600,600);

  TCanvas* cm1 = new TCanvas("cm1","cm1",600,600);
  TCanvas* cm2 = new TCanvas("cm2","cm2",600,600);
  TCanvas* cm3 = new TCanvas("cm3","cm3",600,600);

  TCanvas* cb1 = new TCanvas("cb1","cb1",600,600);
  TCanvas* cb2 = new TCanvas("cb2","cb2",600,600);
  TCanvas* cb3 = new TCanvas("cb3","cb3",600,600);

  double at1[100]; double atx1[100];
  double at2[100]; double atx2[100];
  double at3[100]; double atx3[100];

  double am1[100]; double amx1[100];
  double am2[100]; double amx2[100];
  double am3[100]; double amx3[100];

  double ab1[100]; double abx1[100];
  double ab2[100]; double abx2[100];
  double ab3[100]; double abx3[100];

  at1[0]=4.5; atx1[0]=100;
  at1[1]=5.4; atx1[1]=150;
  at1[2]=5.7; atx1[2]=200;
  at1[3]=6.4; atx1[3]=250;

  at2[0]=4.2; atx2[0]=50;
  at2[1]=5.2; atx2[1]=100;
  at2[2]=6.1; atx2[2]=150;
  at2[3]=7.1; atx2[3]=200;
  at2[4]=7.6; atx2[4]=250;
  
  at3[0]=8.9; atx3[0]=50;
  at3[1]=9.4; atx3[1]=100;
  at3[2]=10.2; atx3[2]=150;
  at3[3]=10.4; atx3[3]=200;

  am1[0]=6; amx1[0]=1050;
  am1[1]=8.9; amx1[1]=1100;
  am1[2]=11.2; amx1[2]=1150;
  am1[3]=12.5; amx1[3]=1200;
  am1[4]=16; amx1[4]=1250;

  am2[0]=4.2; amx2[0]=1000;
  am2[1]=7.8; amx2[1]=1050;
  am2[2]=11.1; amx2[2]=1100;
  am2[3]=14; amx2[3]=1150;
  am2[4]=17.5; amx2[4]=1200;

  am3[0]=5.2; amx3[0]=1000;
  am3[1]=8.6; amx3[1]=1050;
  am3[2]=11.7; amx3[2]=1100;
  am3[3]=15; amx3[3]=1150;

  ab1[0]=4.5; abx1[0]=1500;
  ab1[1]=5.2; abx1[1]=1550;
  ab1[2]=7.2; abx1[2]=1600;
  ab1[3]=10; abx1[3]=1650;
  ab1[4]=12.6; abx1[4]=1700;

  ab2[0]=6; abx2[0]=1500;
  ab2[1]=7.8; abx2[1]=1550;
  ab2[2]=9.9; abx2[2]=1600;
  ab2[3]=13; abx2[3]=1650;
  ab2[4]=17.5; abx2[4]=1700;

  ab3[0]=5.2; abx3[0]=1550;
  ab3[1]=7.2; abx3[1]=1600;
  ab3[2]=10; abx3[2]=1650;
  ab3[3]=12.6; abx3[3]=1700;

  TGraph* gt1 = new TGraph(4,atx1,at1);
  gt1->SetTitle("Changing top voltage (X_1000_1500)");
  gt1->GetXaxis()->SetTitle("Delta V1");
  gt1->GetYaxis()->SetTitle("Amplitude");
  ct1->cd();
  gt1->Draw("AC*");
  gt1->Write();

  TGraph* gt2 = new TGraph(5,atx2,at2);
  gt2->SetTitle("Changing top voltage (X_1000_1550)");
  gt2->GetXaxis()->SetTitle("Delta V1");
  gt2->GetYaxis()->SetTitle("Amplitude");
  ct2->cd();
  gt2->Draw("AC*");
  gt2->Write();

  TGraph* gt3 = new TGraph(4,atx3,at3);
  gt3->SetTitle("Changing top voltage (X_1100_1500)");
  gt3->GetXaxis()->SetTitle("Delta V1");
  gt3->GetYaxis()->SetTitle("Amplitude");
  ct3->cd();
  gt3->Draw("AC*");
  gt3->Write();

  TGraph* gm1 = new TGraph(5,amx1,am1);
  gm1->SetTitle("Changing top&middle voltage (50_X_1500)");
  gm1->GetXaxis()->SetTitle("Delta V2");
  gm1->GetYaxis()->SetTitle("Amplitude");
  cm1->cd();
  gm1->Draw("AC*");
  gm1->Write();

  TGraph* gm2 = new TGraph(5,amx2,am2);
  gm2->SetTitle("Changing top&middle voltage (50_X_1550)");
  gm2->GetXaxis()->SetTitle("Delta V2");
  gm2->GetYaxis()->SetTitle("Amplitude");
  cm2->cd();
  gm2->Draw("AC*");
  gm2->Write();

  TGraph* gm3 = new TGraph(4,amx3,am3);
  gm3->SetTitle("Changing top&middle voltage (100_X_1550)");
  gm3->GetXaxis()->SetTitle("Delta V2");
  gm3->GetYaxis()->SetTitle("Amplitude");
  cm3->cd();
  gm3->Draw("AC*");
  gm3->Write();

  TGraph* gb1 = new TGraph(5,abx1,ab1);
  gb1->SetTitle("Changing top&middle&bottom voltage (100_1100_X)");
  gb1->GetXaxis()->SetTitle("Delta V3");
  gb1->GetYaxis()->SetTitle("Amplitude");
  cb1->cd();
  gb1->Draw("AC*");
  gb1->Write();

  TGraph* gb2 = new TGraph(5,abx2,ab2);
  gb2->SetTitle("Changing top&middle&bottom voltage (50_1100_X)");
  gb2->GetXaxis()->SetTitle("Delta V3");
  gb2->GetYaxis()->SetTitle("Amplitude");
  cb2->cd();
  gb2->Draw("AC*");
  gb2->Write();

  TGraph* gb3 = new TGraph(4,abx3,ab3);
  gb3->SetTitle("Changing top&middle&bottom voltage (100_1050_X)");
  gb3->GetXaxis()->SetTitle("Delta V3");
  gb3->GetYaxis()->SetTitle("Amplitude");
  cb3->cd();
  gb3->Draw("AC*");
  gb3->Write();

}
