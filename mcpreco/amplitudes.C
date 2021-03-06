{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Thu May 18 12:51:50 2017) by ROOT version5.34/04
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",822,76,538,323);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   Canvas_1->Range(-0.1478472,-59.14744,1.672083,376.6278);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetLeftMargin(0.1142061);
   Canvas_1->SetRightMargin(0.08356546);
   Canvas_1->SetBottomMargin(0.1324921);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   TH1F *htemp__1 = new TH1F("htemp__1","fwav[0].amp {fwav[0].evtno<8000&&fwav[0].npeaks==1}",100,0.06,1.52);
   htemp__1->SetBinContent(9,2);
   htemp__1->SetBinContent(10,2);
   htemp__1->SetBinContent(11,1);
   htemp__1->SetBinContent(12,3);
   htemp__1->SetBinContent(13,4);
   htemp__1->SetBinContent(14,1);
   htemp__1->SetBinContent(15,5);
   htemp__1->SetBinContent(16,2);
   htemp__1->SetBinContent(17,3);
   htemp__1->SetBinContent(18,2);
   htemp__1->SetBinContent(19,5);
   htemp__1->SetBinContent(20,2);
   htemp__1->SetBinContent(21,3);
   htemp__1->SetBinContent(22,4);
   htemp__1->SetBinContent(23,13);
   htemp__1->SetBinContent(24,22);
   htemp__1->SetBinContent(25,37);
   htemp__1->SetBinContent(26,52);
   htemp__1->SetBinContent(27,46);
   htemp__1->SetBinContent(28,80);
   htemp__1->SetBinContent(29,92);
   htemp__1->SetBinContent(30,119);
   htemp__1->SetBinContent(31,144);
   htemp__1->SetBinContent(32,170);
   htemp__1->SetBinContent(33,167);
   htemp__1->SetBinContent(34,168);
   htemp__1->SetBinContent(35,190);
   htemp__1->SetBinContent(36,199);
   htemp__1->SetBinContent(37,202);
   htemp__1->SetBinContent(38,213);
   htemp__1->SetBinContent(39,191);
   htemp__1->SetBinContent(40,223);
   htemp__1->SetBinContent(41,229);
   htemp__1->SetBinContent(42,211);
   htemp__1->SetBinContent(43,186);
   htemp__1->SetBinContent(44,201);
   htemp__1->SetBinContent(45,186);
   htemp__1->SetBinContent(46,177);
   htemp__1->SetBinContent(47,212);
   htemp__1->SetBinContent(48,195);
   htemp__1->SetBinContent(49,218);
   htemp__1->SetBinContent(50,151);
   htemp__1->SetBinContent(51,156);
   htemp__1->SetBinContent(52,183);
   htemp__1->SetBinContent(53,150);
   htemp__1->SetBinContent(54,131);
   htemp__1->SetBinContent(55,137);
   htemp__1->SetBinContent(56,141);
   htemp__1->SetBinContent(57,131);
   htemp__1->SetBinContent(58,109);
   htemp__1->SetBinContent(59,103);
   htemp__1->SetBinContent(60,95);
   htemp__1->SetBinContent(61,80);
   htemp__1->SetBinContent(62,80);
   htemp__1->SetBinContent(63,73);
   htemp__1->SetBinContent(64,67);
   htemp__1->SetBinContent(65,68);
   htemp__1->SetBinContent(66,50);
   htemp__1->SetBinContent(67,44);
   htemp__1->SetBinContent(68,36);
   htemp__1->SetBinContent(69,32);
   htemp__1->SetBinContent(70,33);
   htemp__1->SetBinContent(71,27);
   htemp__1->SetBinContent(72,26);
   htemp__1->SetBinContent(73,24);
   htemp__1->SetBinContent(74,19);
   htemp__1->SetBinContent(75,9);
   htemp__1->SetBinContent(76,16);
   htemp__1->SetBinContent(77,10);
   htemp__1->SetBinContent(78,9);
   htemp__1->SetBinContent(79,8);
   htemp__1->SetBinContent(80,6);
   htemp__1->SetBinContent(81,5);
   htemp__1->SetBinContent(82,4);
   htemp__1->SetBinContent(83,6);
   htemp__1->SetBinContent(84,5);
   htemp__1->SetBinContent(85,1);
   htemp__1->SetBinContent(86,2);
   htemp__1->SetBinContent(87,2);
   htemp__1->SetBinContent(88,1);
   htemp__1->SetBinContent(89,2);
   htemp__1->SetBinContent(90,2);
   htemp__1->SetBinContent(91,2);
   htemp__1->SetBinContent(93,1);
   htemp__1->SetMinimum(-1.410661);
   htemp__1->SetMaximum(334.0125);
   htemp__1->SetEntries(6419);
   htemp__1->SetDirectory(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   htemp__1->SetLineColor(ci);
   htemp__1->SetLineWidth(2);
   htemp__1->GetXaxis()->SetTitle("amplitude (mV)");
   htemp__1->GetXaxis()->SetRange(1,100);
   htemp__1->GetXaxis()->CenterTitle(true);
   htemp__1->GetXaxis()->SetLabelFont(42);
   htemp__1->GetXaxis()->SetLabelSize(0.06);
   htemp__1->GetXaxis()->SetTitleSize(0.06);
   htemp__1->GetXaxis()->SetTitleFont(42);
   htemp__1->GetYaxis()->SetLabelFont(42);
   htemp__1->GetYaxis()->SetLabelSize(0.06);
   htemp__1->GetYaxis()->SetTitleSize(0.035);
   htemp__1->GetYaxis()->SetTitleFont(42);
   htemp__1->GetZaxis()->SetLabelFont(42);
   htemp__1->GetZaxis()->SetLabelSize(0.035);
   htemp__1->GetZaxis()->SetTitleSize(0.035);
   htemp__1->GetZaxis()->SetTitleFont(42);
   htemp__1->Draw("");
   
   TH1F *htemp__2 = new TH1F("htemp__2","fwav[0].amp {fwav[0].evtno<8000&&fwav[0].npeaks==1}",100,0.06,1.52);
   htemp__2->SetBinContent(12,1);
   htemp__2->SetBinContent(19,1);
   htemp__2->SetBinContent(20,1);
   htemp__2->SetBinContent(21,12);
   htemp__2->SetBinContent(22,18);
   htemp__2->SetBinContent(23,35);
   htemp__2->SetBinContent(24,69);
   htemp__2->SetBinContent(25,101);
   htemp__2->SetBinContent(26,143);
   htemp__2->SetBinContent(27,197);
   htemp__2->SetBinContent(28,230);
   htemp__2->SetBinContent(29,244);
   htemp__2->SetBinContent(30,285);
   htemp__2->SetBinContent(31,311);
   htemp__2->SetBinContent(32,287);
   htemp__2->SetBinContent(33,311);
   htemp__2->SetBinContent(34,292);
   htemp__2->SetBinContent(35,302);
   htemp__2->SetBinContent(36,271);
   htemp__2->SetBinContent(37,257);
   htemp__2->SetBinContent(38,229);
   htemp__2->SetBinContent(39,210);
   htemp__2->SetBinContent(40,241);
   htemp__2->SetBinContent(41,215);
   htemp__2->SetBinContent(42,201);
   htemp__2->SetBinContent(43,171);
   htemp__2->SetBinContent(44,146);
   htemp__2->SetBinContent(45,135);
   htemp__2->SetBinContent(46,117);
   htemp__2->SetBinContent(47,120);
   htemp__2->SetBinContent(48,92);
   htemp__2->SetBinContent(49,70);
   htemp__2->SetBinContent(50,81);
   htemp__2->SetBinContent(51,66);
   htemp__2->SetBinContent(52,48);
   htemp__2->SetBinContent(53,30);
   htemp__2->SetBinContent(54,33);
   htemp__2->SetBinContent(55,41);
   htemp__2->SetBinContent(56,26);
   htemp__2->SetBinContent(57,16);
   htemp__2->SetBinContent(58,19);
   htemp__2->SetBinContent(59,10);
   htemp__2->SetBinContent(60,10);
   htemp__2->SetBinContent(61,13);
   htemp__2->SetBinContent(62,6);
   htemp__2->SetBinContent(63,6);
   htemp__2->SetBinContent(64,4);
   htemp__2->SetBinContent(65,3);
   htemp__2->SetBinContent(66,2);
   htemp__2->SetBinContent(67,2);
   htemp__2->SetBinContent(68,1);
   htemp__2->SetBinContent(70,1);
   htemp__2->SetBinContent(71,1);
   htemp__2->SetBinContent(73,1);
   htemp__2->SetBinContent(74,1);
   htemp__2->SetBinContent(75,1);
   htemp__2->SetMinimum(-1.410661);
   htemp__2->SetMaximum(334.0125);
   htemp__2->SetEntries(5737);
   htemp__2->SetDirectory(0);
   htemp__2->SetFillColor(51);
   htemp__2->SetFillStyle(3005);

   ci = TColor::GetColor("#663366");
   htemp__2->SetLineColor(ci);
   htemp__2->GetXaxis()->SetTitle("fwav[0].amp");
   htemp__2->GetXaxis()->SetLabelFont(42);
   htemp__2->GetXaxis()->SetLabelSize(0.035);
   htemp__2->GetXaxis()->SetTitleSize(0.035);
   htemp__2->GetXaxis()->SetTitleFont(42);
   htemp__2->GetYaxis()->SetLabelFont(42);
   htemp__2->GetYaxis()->SetLabelSize(0.035);
   htemp__2->GetYaxis()->SetTitleSize(0.035);
   htemp__2->GetYaxis()->SetTitleFont(42);
   htemp__2->GetZaxis()->SetLabelFont(42);
   htemp__2->GetZaxis()->SetLabelSize(0.035);
   htemp__2->GetZaxis()->SetTitleSize(0.035);
   htemp__2->GetZaxis()->SetTitleFont(42);
   htemp__2->Draw("same");
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
