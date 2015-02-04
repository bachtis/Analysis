using namespace RooFit;

void  convertToRoo(const char* muon1,const char* muon2,const char* file,const char *newfile,const char* tree,const char* preselection,bool isGEN = false,bool isDATA=false,int min=0, int max=-1) {

  TFile * f = new TFile(file);
  TTree  *t = f->Get(tree);
  RooWorkspace *w = new RooWorkspace("w","w");
  w->factory("curvRaw1[0,1000]");
  w->factory("curvRaw2[0,1000]");
  w->factory("curvGenRaw1[0,1000]");
  w->factory("curvGenRaw2[0,1000]");
  w->factory("etaRaw1[-2.5,2.5]");
  w->factory("etaRaw2[-2.5,2.5]");
  w->factory("phiRaw1[-4,4]");
  w->factory("phiRaw2[-4,4]");
  w->factory("massErrRaw[1.0,0,30.]");
  w->factory("massErrRaw1[1.0,0,30.]");
  w->factory("massErrRaw2[1.0,0,30.]");
  w->factory("muMass[0.1056583715]");
  w->factory("massRaw[0,1000]");
  w->factory("rapidity[-6,6]");
  w->factory("eta[-6,6]");
  w->factory("pt[0.0,200.0]");


  TFile * cache = new TFile("__cacheR__.root","RECREATE");
  TTree*  newtree = t->CopyTree(preselection);
  RooArgSet set;
  set.add(*w->var("curvRaw1"));
  set.add(*w->var("etaRaw1"));
  set.add(*w->var("phiRaw1"));
  set.add(*w->var("curvRaw2"));
  set.add(*w->var("etaRaw2"));
  set.add(*w->var("phiRaw2"));
  set.add(*w->var("massRaw"));
  set.add(*w->var("massErrRaw1"));
  set.add(*w->var("massErrRaw2"));
  set.add(*w->var("massErrRaw"));
  set.add(*w->var("curvGenRaw1"));
  set.add(*w->var("curvGenRaw2"));
  set.add(*w->var("rapidity"));
  set.add(*w->var("eta"));
  set.add(*w->var("pt"));




  RooDataSet *  data  = new RooDataSet("data","DATA",set);

  Double_t pt1,eta1,phi1,pt2,eta2,phi2,ptgen1,ptgen2;
  Double_t covPlus[9];
  Double_t covMinus[9];

  newtree->SetBranchAddress(TString::Format("%s_pt",muon1),&pt1);
  newtree->SetBranchAddress(TString::Format("%s_pt",muon2),&pt2);
  if ((!isGEN) && (!isDATA)) {
    newtree->SetBranchAddress("MuPosGenStatus1_pt",&ptgen1);
    newtree->SetBranchAddress("MuNegGenStatus1_pt",&ptgen2);
  }
  newtree->SetBranchAddress(TString::Format("%s_eta",muon1),&eta1);
  newtree->SetBranchAddress(TString::Format("%s_eta",muon2),&eta2);
  newtree->SetBranchAddress(TString::Format("%s_phi",muon1),&phi1);
  newtree->SetBranchAddress(TString::Format("%s_phi",muon2),&phi2);
  newtree->SetBranchAddress("MuPosCovMatrix",covPlus);
  newtree->SetBranchAddress("MuNegCovMatrix",covMinus);


  TLorentzVector *v1 = new TLorentzVector ();
  TLorentzVector *v2 = new TLorentzVector ();

  int lim=0;
  if (max<0 || max>newtree->GetEntries())
    lim=newtree->GetEntries();
  else
    lim=max;
  
  for (Int_t i=min;i<lim;++i) {
    newtree->GetEntry(i);

    w->var("curvRaw1")->setVal(1./pt1);
    if (isGEN) {
      w->var("curvGenRaw1")->setVal(1./pt1);
      w->var("curvGenRaw2")->setVal(1./pt2);

    }
    else {
      w->var("curvGenRaw1")->setVal(1./ptgen1);
      w->var("curvGenRaw2")->setVal(1./ptgen2);
    }


    w->var("etaRaw1")->setVal(eta1);
    w->var("phiRaw1")->setVal(phi1);
    w->var("curvRaw2")->setVal(1./pt2);
    w->var("etaRaw2")->setVal(eta2);
    w->var("phiRaw2")->setVal(phi2);

    v1->SetPtEtaPhiM(pt1,eta1,phi1,w->var("muMass")->getVal());
    v2->SetPtEtaPhiM(pt2,eta2,phi2,w->var("muMass")->getVal());

    TLorentzVector *z = new TLorentzVector( (*v1)+(*v2));
 

    w->var("massRaw")->setVal(((*v1)+(*v2)).M());
    w->var("rapidity")->setVal(((*v1)+(*v2)).Rapidity());
    w->var("eta")->setVal(((*v1)+(*v2)).Eta());
    w->var("pt")->setVal(((*v1)+(*v2)).Pt());


    //Event by event error
    TMatrixDSym covMatrix(6);
    TMatrixD jacobian(1,6);
    covMatrix(0,0) = covPlus[0];
    covMatrix(0,1) = covPlus[1];
    covMatrix(0,2) = covPlus[2];
    covMatrix(1,0) = covPlus[3];
    covMatrix(1,1) = covPlus[4];
    covMatrix(1,2) = covPlus[5];
    covMatrix(2,0) = covPlus[6];
    covMatrix(2,1) = covPlus[7];
    covMatrix(2,2) = covPlus[8];

    covMatrix(3,3) = covMinus[0];
    covMatrix(3,4) = covMinus[1];
    covMatrix(3,5) = covMinus[2];
    covMatrix(4,3) = covMinus[3];
    covMatrix(4,4) = covMinus[4];
    covMatrix(4,5) = covMinus[5];
    covMatrix(5,3) = covMinus[6];
    covMatrix(5,4) = covMinus[7];
    covMatrix(5,5) = covMinus[8];

    jacobian(0,0) = (z->Energy()*(v1->Px()/v1->Energy()) - z->Px())/z->M();
    jacobian(0,1) = (z->Energy()*(v1->Py()/v1->Energy()) - z->Py())/z->M();
    jacobian(0,2) = (z->Energy()*(v1->Pz()/v1->Energy()) - z->Pz())/z->M();

    jacobian(0,3) = (z->Energy()*(v2->Px()/v2->Energy()) - z->Px())/z->M();
    jacobian(0,4) = (z->Energy()*(v2->Py()/v2->Energy()) - z->Py())/z->M();
    jacobian(0,5) = (z->Energy()*(v2->Pz()/v2->Energy()) - z->Pz())/z->M();

    /////DO COMPONENTS
    Double_t errors[2];
    
    Int_t offset=0;
      for (unsigned int j=0;j<2;++j) {
	TMatrixDSym bigCovOne(6);
	for(unsigned int ir =0;ir<3;++ir)
	  for(unsigned int ic =0;ic<3;++ic) {
	    bigCovOne[offset+ir][offset+ic]=covMatrix(offset+ir,offset+ic);
	  }
	TMatrixDSym dmOneCov = bigCovOne.Similarity(jacobian);
	float error = dmOneCov(0,0) ;   
	if (error>0.0)
	  errors[j] = (sqrt(error));
	else    
	  errors[j] = 0.0;

	offset=offset+3;     
      }
    
      w->var("massErrRaw1")->setVal(errors[0]);
      w->var("massErrRaw2")->setVal(errors[1]);
      w->var("massErrRaw")->setVal(sqrt(errors[0]*errors[0]+errors[1]*errors[1]));

    ////////////////////
      data->add(set);
      if (max >0 && i==max)
	break;
      
      
  }

  TFile *f2= new TFile(newfile,"RECREATE");
  f2->cd();
  data->Write();
  f2->Close();
  cache->Close();
  f->Close();

}


