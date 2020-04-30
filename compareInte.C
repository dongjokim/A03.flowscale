#include "include/Filipad.h"
#include "include/rootcommon.h"
#include "include/detector.h" // define detector acceptances
int qColor[] = { kBlack, kRed, kBlue, kGreen+2, kOrange-3, kMagenta+1, kViolet, kPink, kGray, 40, 48, 9 };
int qMarker[] = { 25, 24, 21, 28, 29, 33 };

const int NC = 6;
double cent[NC+1] =           { 0, 5,10,20,30,40,50}; 
const int NEnergy = 2; // 2.76 5.02TeV
const int NObs = 4; // v2{2},v2{4},v3{2},v4{2}
TString S_Obs[NObs] = {"v2{2}","v2{4}","v3{2}","v4{2}"}; //https://doi.org/10.17182/hepdata.73940
TString S_Energy[NEnergy] = {"#sqrt{s_{NN}} = 2.76TeV","#sqrt{s_{NN}} = 5.02TeV"};
const int Npt=3;
double startPtbins[3] = {0.,0.2,0.5};
double endPt = 5.0;
TGraphAsymmErrors *gr_v2intForPtbins[Npt]; // 5TeV
TGraphErrors *grflow[D_COUNT][NObs];
TGraphErrors *grflow5tev[D_COUNT];

TGraphErrors *grjoakim[3]; // v2 3 4
void compareJokim();

void LoadData(){

    TFile *finvnptint = TFile::Open("vnPtintegrated.root");
    TFile *finetaint = TFile::Open("vnetaintegrated.root");
    TFile *finetaint5tev = TFile::Open("v2etaintegrated_5tev.root");
    for(int i=0;i<Npt;i++){
        gr_v2intForPtbins[i] = (TGraphAsymmErrors*)finvnptint->Get(Form("gr_v2intForPtbins%02d",i));
    }
  // 2.76 published
  for(int is = 0; is < D_COUNT; is++){
     for(int io=0;io<NObs;io++) {
       grflow[is][io] = (TGraphErrors*)finetaint->Get(Form("grflow_etaintegratedE%02dD%02dO%02d",0,is,io));
     }
  }
  // 5.02 TeV Freja
  for(int is = 0; is < D_COUNT; is++){
       grflow5tev[is] = (TGraphErrors*)finetaint5tev->Get(Form("grflow_etaintegratedE%02dD%02dO%02d",1,is,0));
  }

  TFile *fjk = TFile::Open("data/joakim.root");
  for(int i=0;i<3;i++) {
    grjoakim[i] = (TGraphErrors*)fjk->Get(Form("gr_jk_v%d",i+2));
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
      gr_v2intForPtbins[iPT]->SetMarkerStyle(qMarker[iPT]);
      gr_v2intForPtbins[iPT]->SetMarkerColor(qColor[iPT]);
      gr_v2intForPtbins[iPT]->SetLineColor(qColor[iPT]);
      gr_v2intForPtbins[iPT]->Draw("psame");
      if(iPT==0)leg->AddEntry(gr_v2intForPtbins[iPT],Form("%.1f<p_{T}<%.1f, %s",startPtbins[iPT],endPt, S_Energy[1].Data()),"lp");
      if(iPT==1)leg->AddEntry(gr_v2intForPtbins[iPT],Form("%.1f<p_{T}<%.1f, %s, Ref",startPtbins[iPT],endPt, S_Energy[1].Data()),"lp");
    }
    
    
    int iref=1;
    for(int is = 0; is < D_COUNT; is++){
      grflow5tev[is]->SetMarkerStyle(qMarker[is+2]);
      grflow5tev[is]->SetMarkerColor(qColor[is+2]);
      grflow5tev[is]->SetLineColor(qColor[is+2]);
      grflow5tev[is]->Draw("psame");
      //TString label = Form("%2.0f-%2.0f%%",cent[ic],cent[ic+1]);
      leg->AddEntry(grflow5tev[is],Form("%s",grflow5tev[is]->GetTitle()),"lp");
    }
    leg->Draw();
  
    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.3,1.1);
    hset( *hfr1, "centrality[%]","#frac{Data}{Ref}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    TGraphErrors *gr_ratio[D_COUNT];

    for(int is = 0; is < D_COUNT; is++){
      gr_ratio[is] = GetRatio(grflow5tev[is],gr_v2intForPtbins[1]);
      gr_ratio[is]->SetMarkerStyle(qMarker[is+2]);
      gr_ratio[is]->SetMarkerColor(qColor[is+2]);
      gr_ratio[is]->Draw("psame");
    }
    gPad->GetCanvas()->SaveAs("figs/integratedv2_differteta.pdf");
}

void compareJokim(){
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
    hset( *hfr, "centrality[%]", "v_{2,#Delta#eta>1}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.2,0.5,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.04);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    int iPT =1 ;
    gr_v2intForPtbins[iPT]->SetMarkerStyle(qMarker[iPT]);
    gr_v2intForPtbins[iPT]->SetMarkerColor(qColor[iPT]);
    gr_v2intForPtbins[iPT]->SetLineColor(qColor[iPT]);
    gr_v2intForPtbins[iPT]->Draw("psame");
    leg->AddEntry(gr_v2intForPtbins[iPT],Form("%.1f<p_{T}<%.1f, %s",startPtbins[iPT],endPt, S_Energy[1].Data()),"lp");
    
    grjoakim[0]->SetMarkerStyle(qMarker[iPT+1]);
    grjoakim[0]->SetMarkerColor(qColor[iPT+1]);
    grjoakim[0]->Draw("psame");
    leg->AddEntry(grjoakim[0],"Joakim's ana","p");
    leg->Draw();
    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.8,1.2);
    hset( *hfr1, "centrality[%]","#frac{Data}{Ref}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    TGraphErrors *gr_ratio;


    gr_ratio = GetRatio(grjoakim[0],gr_v2intForPtbins[iPT]);
    gr_ratio->SetMarkerStyle(qMarker[iPT+1]);
    gr_ratio->SetMarkerColor(qColor[iPT+1]);
    gr_ratio->Draw("psame");
    
    gPad->GetCanvas()->SaveAs("figs/integratedv2_jkcomp.pdf");
}