#include "include/Filipad.h"
#include "include/rootcommon.h"
#include "include/detector.h" // define detector acceptances

const int NC = 6;
double cent[NC+1] =           { 0, 5,10,20,30,40,50}; 
const int NEnergy = 2; // 2.76 5.02TeV
const int NObs = 4; // v2{2},v2{4},v3{2},v4{2}
TString S_Obs[NObs] = {"v2{2}","v2{4}","v3{2}","v4{2}"}; //https://doi.org/10.17182/hepdata.73940
TString S_Energy[NEnergy] = {"#sqrt{s_{NN}} = 2.76TeV","#sqrt{s_{NN}} = 5.02TeV"};
const int Npt=3;
double startPtbins[3] = {0.,0.2,0.5};
double endPt = 5.0;
TGraphAsymmErrors *gr_v2intForPtbins[Npt];
TGraphErrors *grflow[D_COUNT][NObs];

void LoadData(){

    TFile *finvnptint = TFile::Open("vnPtintegrated.root");
    TFile *finetaint = TFile::Open("vnetaintegrated.root");
    for(int i=0;i<Npt;i++){
        gr_v2intForPtbins[i] = (TGraphAsymmErrors*)finvnptint->Get(Form("gr_v2intForPtbins%02d",i));
    }

  for(int is = 0; is < D_COUNT; is++){
     for(int io=0;io<NObs;io++) {
       grflow[is][io] = (TGraphErrors*)finetaint->Get(Form("grflow_etaintegratedE%02dD%02dO%02d",0,is,io));
     }
  }
}

void compareInte(){
  LoadData();
  Filipad *fpad;
  fpad = new Filipad(1, 0.9, 0.4, 100, 100, 1.1,7);
    fpad->Draw();
    //==== Upper pad
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    double lowx = 0,highx=60.;
    double ly=0,hy=0.2;
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "centrality[%]", "v_{2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.2,0.5,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.04);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    for(int iPT=0;iPT<2;iPT++) {
      gr_v2intForPtbins[iPT]->SetMarkerStyle(20+iPT);
      gr_v2intForPtbins[iPT]->SetMarkerColor(1+iPT);
      gr_v2intForPtbins[iPT]->Draw("psame");
      leg->AddEntry(gr_v2intForPtbins[iPT],Form("%.1f<p_{T}<%.1f",startPtbins[iPT],endPt),"lp");
    }
    
    int io=0;// v2{2}
    for(int is = 0; is < D_COUNT; is++){
      grflow[is][io]->SetMarkerStyle(22+is);
      grflow[is][io]->SetMarkerColor(3+is);
      grflow[is][io]->SetLineColor(3+is);
      grflow[is][io]->Draw("psame");
      //TString label = Form("%2.0f-%2.0f%%",cent[ic],cent[ic+1]);
      leg->AddEntry(grflow[is][io],Form("%s, %.1f<#eta<%.1f",pdetn[is],decAcc[is][0],decAcc[is][1]),"lp");
    }
    leg->Draw();
  
    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0,5);
    hset( *hfr1, "centrality[%]","#frac{Data}{Ref}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    TGraphAsymmErrors *gr_ratio[NC];
    int iref=0;
   //gPad->GetCanvas()->SaveAs("figs/integratedv2_ptdiff.pdf");
}