#include "include/Filipad.h"
#include "include/rootcommon.h"
#include "include/detector.h" // define detector acceptances

const int NC = 6;
double cent[NC+1] =           { 0, 5,10,20,30,40,50}; 
const int NEnergy = 2; // 2.76 5.02TeV
const int NObs = 4; // v2{2},v2{4},v3{2},v4{2}
TString S_Obs[NObs] = {"v2{2}","v2{4}","v3{2}","v4{2}"}; //https://doi.org/10.17182/hepdata.73940
TString S_Energy[NEnergy] = {"#sqrt{s_{NN}} = 2.76TeV","#sqrt{s_{NN}} = 5.02TeV"};
TGraphAsymmErrors *gr_vneta[NEnergy][NObs][NC];
TF1 *fgaus[NEnergy][NObs][NC];
// 5TeV PbPb in TH1D from Freja
TH1D *h5tevv2[NC];
TH1D *h5tevv2_ampt[NC];

void LoadHEPData();
void DrawData(int io);
void PlotCentralityDep();
void integratedvneta();



void LoadHEPData() {
	TFile *fin = TFile::Open("data/HEPData-ins1456145-v1-2.76TeVVneta.root");  //https://doi.org/10.17182/hepdata.83737
  for(int io=0;io<NObs;io++) {
    for(int ic=0;ic<NC;ic++) {
      TString title = Form("%s, %2.0f-%2.0f%%",S_Obs[io].Data(),cent[ic],cent[ic+1]);
	   gr_vneta[0][io][ic] = (TGraphAsymmErrors*)fin->Get(Form("Table %d/Graph1D_y%d",ic+1,io+1));
     gr_vneta[0][io][ic]->Fit("gaus");
     gr_vneta[0][io][ic]->SetTitle(title);
     fgaus[0][io][ic] = gr_vneta[0][io][ic]->GetFunction("gaus");
     fgaus[0][io][ic]->SetTitle(title);
    }
  }
  TFile *fin_5tev = TFile::Open("data/FT_v2.root");// from Freja
  for(int ic=0;ic<NC;ic++) {
    h5tevv2[ic] = (TH1D*)fin_5tev->Get(Form("v2_cent_%.0f_%.0f",cent[ic],cent[ic+1]));
    h5tevv2_ampt[ic] = (TH1D*)fin_5tev->Get(Form("ampt_v2_cent_%.0f_%.0f",cent[ic],cent[ic+1]));
    h5tevv2[ic]->Fit("gaus");
    fgaus[1][0][ic] = h5tevv2[ic]->GetFunction("gaus");
  }
}

// Plot cenetrality dep
void PlotCentralityDep(){
   LoadHEPData();
   for(int io=0;io<NObs;io++) {
     for(int ic=0;ic<NC;ic++) {
      DrawData(io);
     }
   }
}

// integrate over eta acceptance from the fitted functions
void integratedvneta5TeV(){
  LoadHEPData();
  double vnintegrated[D_COUNT][NC];
  double v2int;
  for(int is = 0; is < D_COUNT; is++){
     for(int ic=0;ic<NC;ic++) {
        // need to normalize by what ? to get mean
        double norm = TMath::Abs(decAcc[is][0]-decAcc[is][1]);
        vnintegrated[is][ic] = fgaus[1][0][ic]->Integral(decAcc[is][0],decAcc[is][1])/norm;
        //cout << pdetn[is] <<"\t"<< S_Obs[io]<<"\t"<< ic<<"\t"<< vnintegrated[is][io][ic] << endl;
      }
  }
  double centmean[NC];
  for(int ic=0;ic<NC;ic++) {
    centmean[ic] = cent[ic]+(cent[ic+1]-cent[ic])/2.;
  }
  // now calculate v2 as a function of cent for each dector and compare..
  TGraphErrors *grflow[D_COUNT];
  for(int is = 0; is < D_COUNT; is++){
       grflow[is] = new TGraphErrors(NC,centmean,vnintegrated[is],0,0);
       grflow[is]->SetName(Form("grflow_etaintegratedE%02dD%02dO%02d",1,is,0));
       grflow[is]->SetTitle(Form("%s, %s, %.1f<#eta<%.1f %s",S_Obs[0].Data(),pdetn[is],decAcc[is][0],decAcc[is][1],S_Energy[1].Data()));
  }
  TFile *fout = new TFile("v2etaintegrated_5tev.root","recreate");
  fout->cd();
  for(int is = 0; is < D_COUNT; is++){
        grflow[is]->Write();
  }
  fout->Close();
}

// integrate over eta acceptance from the fitted functions
void integratedvneta(){
  LoadHEPData();
  double vnintegrated[D_COUNT][NObs][NC];
  double v2int;
  for(int is = 0; is < D_COUNT; is++){
    for(int io=0;io<NObs;io++) {
     for(int ic=0;ic<NC;ic++) {
        // need to normalize by what ? to get mean
        double norm = TMath::Abs(decAcc[is][0]-decAcc[is][1]);
        vnintegrated[is][io][ic] = fgaus[0][io][ic]->Integral(decAcc[is][0],decAcc[is][1])/norm;
        //cout << pdetn[is] <<"\t"<< S_Obs[io]<<"\t"<< ic<<"\t"<< vnintegrated[is][io][ic] << endl;
      }
    }
  }
  double centmean[NC];
  for(int ic=0;ic<NC;ic++) {
    centmean[ic] = cent[ic]+(cent[ic+1]-cent[ic])/2.;
  }
  // now calculate v2 as a function of cent for each dector and compare..
  TGraphErrors *grflow[D_COUNT][NObs];
  for(int is = 0; is < D_COUNT; is++){
     for(int io=0;io<NObs;io++) {
       grflow[is][io] = new TGraphErrors(NC,centmean,vnintegrated[is][io],0,0);
       grflow[is][io]->SetName(Form("grflow_etaintegratedE%02dD%02dO%02d",0,is,io));
       grflow[is][io]->SetTitle(Form("%s, %s, %.1f<#eta<%.1f %s",S_Obs[io].Data(),pdetn[is],decAcc[is][0],decAcc[is][1],S_Energy[0].Data()));
     }
  }
  // Draw here...
  int io = 1;
  Filipad *fpad;
  fpad = new Filipad(io+1, 0.9, 0.4, 100, 100, 1.1,7);
  fpad->Draw();
  //==== Upper pad
  TPad *p = fpad->GetPad(1); //upper pad
  p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
  double lowx = 0,highx=60.;
  double ly=0,hy=0.18;
  TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
  hset( *hfr, "#eta", S_Obs[io],1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
  hfr->Draw();
  //Legend definition
  TLegend *leg = new TLegend(0.2,0.5,0.85,0.78,"","brNDC");
  leg->SetTextSize(0.04);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
  leg->AddEntry(grflow[0][io],Form("%s",S_Obs[io].Data()),"");
  for(int is = 0; is < D_COUNT; is++){
      grflow[is][io]->SetMarkerStyle(20+is);
      grflow[is][io]->SetMarkerColor(1+is);
      grflow[is][io]->SetLineColor(1+is);
      grflow[is][io]->Draw("psame");
      //TString label = Form("%2.0f-%2.0f%%",cent[ic],cent[ic+1]);
      leg->AddEntry(grflow[is][io],Form("%s, %.1f<#eta<%.1f",pdetn[is],decAcc[is][0],decAcc[is][1]),"lp");
  }

  leg->Draw();
  
  //==== Lower pad
  p = fpad->GetPad(2);
  p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
  TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.,1.2);
  hset( *hfr1, "centrality[%]","#frac{Data}{Ref}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
  hfr1->Draw();
    
  TGraphErrors *gr_ratio[D_COUNT];
  int iref=0;
  for(int is = 0; is < D_COUNT; is++){
    gr_ratio[is] = GetRatio(grflow[is][io],grflow[iref][io]);
    gr_ratio[is]->SetMarkerStyle(20+is);
    gr_ratio[is]->SetMarkerColor(1+is);
    gr_ratio[is]->Draw("psame");
  }
  TFile *fout = new TFile("vnetaintegrated.root","recreate");
  fout->cd();
  for(int is = 0; is < D_COUNT; is++){
     for(int io=0;io<NObs;io++) {
        grflow[is][io]->Write();
      }
  }
  fout->Close();
   //gPad->GetCanvas()->SaveAs("figs/integratedv2_ptdiff.pdf");
 
}

// Draw loaded HEP data for one selected observable.
void DrawData(int io){
 
  Filipad *fpad;
  fpad = new Filipad(io+1, 0.9, 0.4, 100, 100, 1.1,7);
    fpad->Draw();
    //==== Upper pad
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
    double lowx = -5.,highx=5.;
    double ly=0,hy=0.18;
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "#eta", "Obs",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.2,0.5,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.04);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    leg->AddEntry(gr_vneta[0][io][0],Form("%s",S_Obs[io].Data()),"");
    for(int ic=0;ic<NC;ic++){
      gr_vneta[0][io][ic]->SetMarkerStyle(20+ic);
      gr_vneta[0][io][ic]->SetMarkerColor(1+ic);
      gr_vneta[0][io][ic]->SetLineColor(1+ic);
      gr_vneta[0][io][ic]->Draw("psame");
      TString label = Form("%2.0f-%2.0f%%",cent[ic],cent[ic+1]);
      leg->AddEntry(gr_vneta[0][io][ic],Form("%s",label.Data()),"lp");
    }

    leg->Draw();
  
    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0,5);
    hset( *hfr1, "centrality[%]","#frac{Ref-Data}{Ref}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    TGraphAsymmErrors *gr_ratio[NC];
    int iref=0;
    for(int i=0;i<NC;i++){
      gr_ratio[i] = GetRatio(gr_vneta[0][io][i],gr_vneta[0][io][iref]);
      gr_ratio[i]->SetMarkerStyle(20+i);
      gr_ratio[i]->SetMarkerColor(1+i);
      gr_ratio[i]->Draw("psame");
  }
   //gPad->GetCanvas()->SaveAs("figs/integratedv2_ptdiff.pdf");
}

