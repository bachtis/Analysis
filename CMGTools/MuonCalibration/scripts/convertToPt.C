TGraphAsymmErrors* convertToPt(TH1D* curv) {
  TGraphAsymmErrors * g = new TGraphAsymmErrors();
  int N=0;
  for(int i=curv->GetNbinsX();i>=1;i=i-1) {
    double bincenter =curv->GetXaxis()->GetBinCenter(i);
    double bincontent =curv->GetBinContent(i);
    double binerror =curv->GetBinError(i);
    double binlo = curv->GetXaxis()->GetBinLowEdge(i);
    double binhi = curv->GetXaxis()->GetBinUpEdge(i);
 
    g->SetPoint(N,1./bincenter,1/bincontent);
    g->SetPointError(N,1./bincenter - 1./binhi,1/binlo - 1./bincenter,binerror/(bincontent*bincontent),binerror/(bincontent*bincontent));
    N=N+1;
  }
  return g;
}




TGraphAsymmErrors* convertToPtS(TH1D* curv,TH1D*curv2) {
  TGraphAsymmErrors * g = new TGraphAsymmErrors();
  int N=0;
  for(int i=curv2->GetNbinsX();i>=1;i=i-1) {
    double bincenter =-curv2->GetXaxis()->GetBinCenter(i);
    double bincontent =curv2->GetBinContent(i);
    double binerror =curv2->GetBinError(i);
    double binhi = -curv2->GetXaxis()->GetBinLowEdge(i);
    double binlo = -curv2->GetXaxis()->GetBinUpEdge(i);
 
    g->SetPoint(N,1./bincenter,1/bincontent);
    g->SetPointError(N,1./bincenter - 1./binhi,1/binlo - 1./bincenter,binerror/(bincontent*bincontent),binerror/(bincontent*bincontent));
    N=N+1;
  }
  for(int i=curv->GetNbinsX();i>=1;i=i-1) {
    double bincenter =curv->GetXaxis()->GetBinCenter(i);
    double bincontent =curv->GetBinContent(i);
    double binerror =curv->GetBinError(i);
    double binlo = curv->GetXaxis()->GetBinLowEdge(i);
    double binhi = curv->GetXaxis()->GetBinUpEdge(i);
 
    g->SetPoint(N,1./bincenter,1/bincontent);
    g->SetPointError(N,1./bincenter - 1./binhi,1/binlo - 1./bincenter,binerror/(bincontent*bincontent),binerror/(bincontent*bincontent));
    N=N+1;
  }


  return g;
}


TGraphAsymmErrors* convertToPtSc(TH1D* curv,TH1D*curv2) {
  TGraphAsymmErrors * g = new TGraphAsymmErrors();
  int N=0;

  for(int i=curv2->GetNbinsX();i>=1;i=i-1) {
    double bincenter =-curv2->GetXaxis()->GetBinCenter(i);
    double bincontent =curv2->GetBinContent(i);
    double binerror =curv2->GetBinError(i);
    double binhi = -curv2->GetXaxis()->GetBinLowEdge(i);
    double binlo = -curv2->GetXaxis()->GetBinUpEdge(i);
 
    g->SetPoint(N,1./bincenter,bincontent);
    g->SetPointError(N,1./bincenter - 1./binhi,1/binlo - 1./bincenter,binerror,binerror);
    N=N+1;
  }


  for(int i=curv->GetNbinsX();i>=1;i=i-1) {

    double bincenter =curv->GetXaxis()->GetBinCenter(i);
    double bincontent =curv->GetBinContent(i);
    double binerror =curv->GetBinError(i);
    double binlo = curv->GetXaxis()->GetBinLowEdge(i);
    double binhi = curv->GetXaxis()->GetBinUpEdge(i);
 
    g->SetPoint(N,1./bincenter,bincontent);
  
    g->SetPointError(N,1./bincenter - 1./binhi,1/binlo - 1./bincenter,binerror,binerror);
    N=N+1;
  }


  return g;
}
