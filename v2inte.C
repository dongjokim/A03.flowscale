#include "include/Filipad.h"
#include "include/rootcommon.h"
double integralv2(TF1 *fflow, TF1 *fspectra, double lpt, double hpt);
double SampingMethod(TF1 *fflow, TF1 *fspectra, double lpt, double hpt);
double FitVN(double *x,double *par);

double LevyTsallisF0(double *x, double *par){
 double mass = par[3];
 return x[0] * par[0] *
         ( ( par[1]-1 )*( par[1]-2 ) ) /
         ( par[1]*par[2] * (par[1]*par[2] + mass*( par[1]-2 ) ) ) *
         pow( ( 1 + ( sqrt( mass*mass + x[0]*x[0] ) - mass ) / (par[1]*par[2]) ), -par[1] );
}


void v2inte(){
	TGraphAsymmErrors *gr_v2pt;
	TGraphAsymmErrors *gr_dndpt;
	TFile *fflow = TFile::Open("data/HEPData-ins1666817-v1-5TeVPbPb_flow.root");
	TFile *fspectra = TFile::Open("data/HEPData-ins1657384-v1-5TeVPbPb_spectra.root");

	// only for 10-20%
	// table 34 : v2{2,|Δη|>1.}
	gr_v2pt = (TGraphAsymmErrors*)fflow->Get("Table 33/Graph1D_y1");
	gr_dndpt = (TGraphAsymmErrors*)fspectra->Get("Table 2/Graph1D_y3");
	// Levy fit
	double norm = 1.49e+04, slope = 6.1;
	double lowpt=0.;
	double highpt=50.0;
	//TF1 *fLevy = new TF1("fLevy","[0]*exp(-[2]/x)/pow(x,[1])",lowpt,highpt);
	//fLevy->SetParameters(norm,slope,4);
	//gr_dndpt->Fit(fLevy,"R");

	TF1 *fLevy = new TF1("fLevy",LevyTsallisF0,0,100,4); // function to describe the pt spectrum
	fLevy->SetParNames("N","n","C","mass");
	fLevy->SetParLimits( 0, 1e2, 1e04);
	fLevy->SetParLimits( 1, 5, 10);
	fLevy->SetParLimits( 2, 0.05, 0.22);
	fLevy->SetParLimits( 3, 0.02, 0.15);
	fLevy->SetParameters(1.2e+03,7,0.02);

	gr_dndpt->Fit(fLevy,"R");
	//TF1 *fLevy = new TF1("fLevy","[0] + -10*x",lowpt,highpt);
	//fLevy->SetParameters(50,1);
	// Fit vn
	TF1 *f_vn = new TF1("f_vn",FitVN,0.,20.,6);
	f_vn->SetParameters(0,1.2,2.2,2,0.1,1.2);
	f_vn->SetParName(0,"a");
	f_vn->SetParName(1,"n");
	f_vn->SetParName(2,"lamda");
	f_vn->SetParName(3,"m");
	f_vn->SetParName(4,"c1");
	f_vn->SetParName(5,"c2");
	f_vn->SetLineColor(2);
	gr_v2pt->Fit(f_vn,"R");

	for(int j=0;j<6;j++) cout << f_vn->GetParameter(j) << endl;
	double v2_0050 = integralv2(f_vn,fLevy,0.,5.0);
	double v2_0250 = integralv2(f_vn,fLevy,0.2,5.0);
	double v2_0550 = integralv2(f_vn,fLevy,0.5,5.0);
	double v2_samp1 = SampingMethod(f_vn,fLevy,0.2,5.0);
    double v2_samp2 = SampingMethod(f_vn,fLevy,0.,5.0);
	cout << "v2_0050 = "<< v2_0050 << endl;
	cout << "v2_0250 = "<< v2_0250 << endl;
	cout << "v2_0550 = "<< v2_0550 << endl;
	cout << "v2_0050/v2_0250 = "<< v2_0050/v2_0250 << endl;
	cout << "v2_0550/v2_0250 = "<< v2_0550/v2_0250 << endl;
	
	double lowx=0., highx=20.;
    double ly=0., hy=0.20;
    int logx = 0;

    // v2
	Filipad *fpad;
	fpad = new Filipad(1, 0.9, 0.4, 100, 100, 1.2,2);
    fpad->Draw();
    //==== Upper pad
    TPad *p = fpad->GetPad(1); //upper pad
    p->SetTickx(); p->SetLogx(logx); p->SetLogy(0); p->cd();
    TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr, "p_{T}(GeV/c)", "v_{2}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr->Draw();
    //Legend definition
    TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;
    gr_v2pt->Draw("same");

    //==== Lower pad
    p = fpad->GetPad(2);
    p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
    TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -1,1);
    hset( *hfr1, "p_{T}(GeV/c)","#frac{Theory-Data}{Theory}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr1->Draw();
    
    TGraphAsymmErrors *gr_ratio1;
    gr_ratio1 = GetDataOverTheory(gr_v2pt,f_vn);
    gr_ratio1->Draw("same");

    // pt spectra
    logx=1;
	Filipad *fpad1;
	fpad1 = new Filipad(2, 0.9, 0.4, 100, 100, 1.2,2);
    fpad1->Draw();
    //==== Upper pad
    TPad *p1 = fpad1->GetPad(1); //upper pad
    p1->SetTickx(); p1->SetLogx(logx); p1->SetLogy(1); p1->cd();
    ly=1.5e-5,hy=9e3;
    lowx=0.08,highx=70.;
    TH2F *hfr01 = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
    hset( *hfr01, "p_{T}(GeV/c)", "#frac{1}{N_{evt}} #frac{d^{2} N}{dp_{T}d#eta}(Gev^{-1}c)",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
    hfr01->Draw();
    //Legend definition
    TLegend *leg1 = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
    leg1->SetTextSize(0.037);leg1->SetBorderSize(0);leg1->SetFillStyle(0);//legend settings;
    gr_dndpt->Draw("same");

    //==== Lower pad
    p1 = fpad1->GetPad(2);
    p1->SetTickx(); p1->SetGridy(1); p1->SetLogx(logx), p1->SetLogy(0); p1->cd();
    TH2F *hfr02 = new TH2F("hfr1"," ", 100, lowx, highx, 10, -1,1);
    hset( *hfr02, "p_{T}(GeV/c)","#frac{Theory-Data}{Theory}",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);
    hfr02->Draw();
    
    TGraphAsymmErrors *gr_ratio2;
    gr_ratio2 = GetDataOverTheory(gr_dndpt,fLevy);
    gr_ratio2->Draw("same");

}


double SampingMethod(TF1 *fflow, TF1 *fspectra, double lpt, double hpt){
	for(int j=0;j<3;j++) cout << fspectra->GetParameter(j) << endl;
	TH1D *hv2 = new TH1D("hv2","",500,0.,1.0);
	Int_t Nevt  = 1e3;
	for(Int_t i = 0; i < Nevt; i++){
		double pt = fspectra->GetRandom(lpt,hpt);
		double v2 = fflow->Eval(pt);
		hv2->Fill(v2);
	}
	cout << "Simpling Done" << endl;
	cout << lpt <<"-"<<hpt <<"="<< hv2->GetMean() << endl;
	return hv2->GetMean();
}

double integralv2(TF1 *fflow, TF1 *fspectra, double lpt, double hpt){

	int nstep=1000;
	double bw = (hpt-lpt)/nstep;
	double norm = fspectra->Integral(lpt,hpt);
	double sum = 0;

	for(int i=0;i<nstep;i++) {
		double lx = lpt+bw*i;
		double hx = lpt+bw*(i+1);
		double v2 = fflow->Eval(lx+bw/2.); // need bin by bin?
		double sub = fspectra->Integral(lx,hx);
		double weight = sub/norm;
		sum = sum + v2*weight;
		//cout << i <<"\t"<< lx <<"\t"<< hx <<"\t"<< v2 <<"\t"<< norm <<"\t"<< sub <<"\t"<< weight <<"\t" <<sum<< endl;
	}		
	return sum;
}

double FitVN(double *x,double *par){
	double a = par[0];
	double n = par[1];
	double lamda = par[2];
	double m = par[3];
	double c1 = par[4];
	double c2 = par[5];
	if(x[0]<4.3) {
		return  (a + 1/TMath::Power(x[0],n))*(TMath::Power(x[0]/lamda,m)/(1+TMath::Power(x[0]/lamda,m)));            }
	else {
		return  (c1*TMath::Power(x[0],c2));
	}
}

