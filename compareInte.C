#include "include/Filipad.h"
#include "include/rootcommon.h"
#include "include/detector.h" // define detector acceptances

const int NC = 6;
double cent[NC+1] =           { 0, 5,10,20,30,40,50}; 
const int NEnergy = 2; // 2.76 5.02TeV
const int NObs = 4; // v2{2},v2{4},v3{2},v4{2}
TString S_Obs[NObs] = {"v2{2}","v2{4}","v3{2}","v4{2}"}; //https://doi.org/10.17182/hepdata.73940
TString S_Energy[NEnergy] = {"#sqrt{s_{NN}} = 2.76TeV","#sqrt{s_{NN}} = 5.02TeV"};

void compareInte(){

    TFile *fout = new TFile("vnPtintegrated.root","recreate");
    fout->cd();
    for(int i=0;i<Npt;i++){
        gr_v2intForPtbins[i]->SetName(Form("gr_v2intForPtbins%02d",i));
        gr_v2intForPtbins[i]->SetTitle(Form("%.1f<p_{T}<%.1f",startPtbins[i],endPt));
        gr_v2intForPtbins[i]->Write();
    }

  for(int is = 0; is < D_COUNT; is++){
     for(int io=0;io<NObs;io++) {
       grflow[is][io] = new TGraphErrors(NC,centmean,vnintegrated[is][io],0,0);
       grflow[is][io]->SetName(Form("grflow_etaintegratedE%02dD%02dO%02d",0,is,io));
       grflow[is][io]->SetTitle(Form("%s, %s, %.1f<#eta<%.1f %s",S_Obs[io].Data(),pdetn[is],decAcc[is][0],decAcc[is][1],S_Energy[0].Data()));
     }
  }
}
