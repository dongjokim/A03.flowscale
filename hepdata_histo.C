void makeHistHEPDATA(TH1F *srch, TH1F *he1, TH1F *he2, TH1F *tarh);

void hepdata_histo(){

	TH1F *hy1,*hy1_e1,*hy1_e2;
	TGraphAsymmErrors *gr_y1;
	TFile *fin = TFile::Open("data/HEPData-ins1666817-v1-Table_33.root");
	hy1 = (TH1F*)fin->Get("Table 33/Hist1D_y1");
	hy1_e1 = (TH1F*)fin->Get("Table 33/Hist1D_y1_e1");
	hy1_e2 = (TH1F*)fin->Get("Table 33/Hist1D_y1_e2");
	gr_y1 = (TGraphAsymmErrors*)fin->Get("Table 33/Graph1D_y1");

	TH1F *hy1_new = (TH1F*)hy1->Clone();

	TCanvas *c1 = new TCanvas();
	gr_y1->SetMarkerStyle(20);
	gr_y1->SetLineStyle(2);
	gr_y1->SetMarkerColor(2);
	gr_y1->SetLineColor(2);
	gr_y1->Draw("ap");

	makeHistHEPDATA(hy1,hy1_e1,hy1_e2,hy1_new);
	//hy1_new->SetLineWidth(.2);
	hy1_new->SetLineStyle(5);
	hy1_new->Draw("same");
}

void makeHistHEPDATA(TH1F *srch, TH1F *he1, TH1F *he2, TH1F *tarh)
{
    const int nb=srch->GetNbinsX();
    for(int i=1;i<=nb;i++) {
    	double rel_e1 = he1->GetBinContent(i);
    	double rel_e2 = he2->GetBinContent(i);
    	double toterr = TMath::Sqrt(rel_e1*rel_e1+rel_e2*rel_e2);
    	//cout << i <<"\t"<< rel_e1<<"\t"<< rel_e2<<"\t"<<toterr<<endl;
    	tarh->SetBinContent(i,srch->GetBinContent(i));
    	tarh->SetBinError(i,toterr);
    }
}