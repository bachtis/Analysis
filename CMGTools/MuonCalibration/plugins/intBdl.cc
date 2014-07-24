/** \file
 *  A simple program to print field value.
 *
 *  $Date: 2009/03/19 10:30:20 $
 *  $Revision: 1.2 $
 *  \author N. Amapane - CERN
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

//#include "DataFormats/GeometryVector/interface/Pi.h"
#include "DataFormats/GeometryVector/interface/CoordinateSets.h"

#include "Utilities/Timing/interface/TimingReport.h"

#include <iostream>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
// #include <boost/program_options/options_description.hpp>
// #include <boost/program_options/variables_map.hpp>
// #include <boost/program_options/parsers.hpp>
// #include <boost/program_options/variables_map.hpp>
// #include <boost/tokenizer.hpp>
// #include <boost/token_functions.hpp>
#include "TH2F.h"
#include "TFile.h"

using namespace edm;
using namespace Geom;
using namespace std;

class intBdl : public edm::EDAnalyzer {
public:
  intBdl(const edm::ParameterSet& pset) {    
  }

  ~intBdl(){}

  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup) {

    ESHandle<MagneticField> magfield;
    setup.get<IdealMagneticFieldRecord>().get(magfield);
    field = magfield.product();

    trackerNew();
    //    scanBarrel();
    //    scanEndcaps();
  }



  void trackerNew() {
    TFile *f = new TFile("map.root","RECREATE");
    TH2F * h = new TH2F("map","map",100,-Geom::pi(),Geom::pi(),100,-2.4,2.4);


    //    string file = "tracker.txt";
    //    ofstream of(file.c_str(),ios_base::trunc);

    int np = 5000;  
    float Rmax   = 112;
    float Zmax   = 274;

    for (int i=1;i< h->GetNbinsX()+1;++i) {
      double phi = h->GetXaxis()->GetBinCenter(i);

      for (int j=1;j<h->GetNbinsY()+1;++j) {
	double eta = h->GetYaxis()->GetBinCenter(j);
	double theta = acos(sinh(eta)/cosh(eta));

	//End at either Rmax or Zmax
	if (fabs(eta)>0.00001) {
		  if (Rmax > fabs(Zmax*tan(theta))) Rmax = fabs(Zmax*tan(theta));
		}
	GlobalPoint start(0,0,0);
	GlobalPoint end(GlobalPoint::Cylindrical(Rmax,phi,Rmax/tan(theta)));
	GlobalVector step;
	double Bzdr, Brdz, Bdl;
	integrate(start, end, np, step, Bzdr, Brdz, Bdl);
	int bin = h->GetBin(i,j);
	cout << i << ',' << j << ','<< eta << ',' << phi << ',' << Bdl << endl; 
	h->SetBinContent(bin,Bdl);

      }

    }
    f->cd();
    h->Write();
    f->Close();


  }



  void scanBarrel() {
    float localPhi=0; // note: no check that this is in the correct sector
    float etastart =-1.2;
    float etaend =+1.2;
    float steps = 300;
    int np = 7000;

    // Distance of detector planes from beamline
    float MB[5] = {0., 444.0, 523.8, 630.6, 753.2};
    float Zmax  = 660; // End of sensitive region: 653.25


    for (int i =0;i<4;++i){ // Layer
      string file = "barrel_L" + boost::lexical_cast<string>(i)+".txt";
      ofstream of(file.c_str(),ios_base::trunc);

      for (int sector = 1; sector<=12; ++sector) {
//	string file = "barrel_L" + boost::lexical_cast<string>(i)+"_s"+ boost::lexical_cast<string>(sector)+".txt";
//	ofstream of(file.c_str(),ios_base::trunc);
	of << "# Layer/Sector " << i << " " << sector << endl;
	for (int j=0; j<=steps; ++j){ // loop on eta
	  float eta = etastart+(etaend-etastart)*j/steps;
	  GlobalPoint start = pointOnBarrelPlane(MB[i], localPhi, eta, sector);
	  GlobalPoint end   = pointOnBarrelPlane(MB[i+1], localPhi, eta, sector);
	  if (fabs(start.z()) > Zmax || fabs(end.z()) > Zmax) continue;
	  GlobalVector step;
	  double Bzdr, Brdz, Bdl;
	  integrate(start, end, np, step, Bzdr, Brdz, Bdl);
	  of << start.x() << " " << start.y() << " " << start.z() << " " 
	     << end.x() << " " << end.y() << " " << end.z() << " " 
	     << step.eta() << " " << step.phi() << " "
	     << Bzdr << " " << Brdz << " " << Bdl	
	     << endl;
	}
	of << endl << endl;
      }
      of.close();
    }
  }
  
  void scanEndcaps() {
    float localPhi=0; // note: no check that this is in the correct sector
    float etastart =0.9;
    float etaend   =2.4;
    float steps = 300;
    int np = 7000;


    float Rmax  = 700; 

    float MElow[5] = {0., 603.2, 828.35, 935.65, 1026.55};
    float MEhigh[5] = {0., 696.85, 828.35, 935.65, 1026.55};

    for (int i =0;i<4;++i){ // Layer
      string file = "endcap_L" + boost::lexical_cast<string>(i)+".txt";
      ofstream of(file.c_str(),ios_base::trunc);

      for (int sector = 1; sector<=12; ++sector) {
	of << "# Layer/Sector " << i << " " << sector << endl;
	float phi = localPhi+(sector-1)*Geom::pi()/6.;

	float prevEta = -999;
	for (int j=0; j<=steps; ++j) { // loop on eta
	  float eta = etastart+(etaend-etastart)*j/steps;
	  float Zstart, Zend;
	  if (eta<1.6) {
	    Zstart = MElow[i];
	    Zend   = MElow[i+1];
	    if (i<=1 && prevEta>=1.6) of << endl;
	  } else {
	    Zstart = MEhigh[i];
	    Zend   = MEhigh[i+1];
	    if (i<=1 && prevEta<1.6) of << endl;
	  }

	  prevEta = eta;
	  
	  float theta = acos(sinh(eta)/cosh(eta));
	  float Rstart = fabs(Zstart*tan(theta));
	  float Rend = fabs(Zend*tan(theta));

	  if (Rend>Rmax) continue; // skip this eta;	  
	    
	  GlobalPoint start(GlobalPoint::Cylindrical(Rstart,phi,Zstart));
	  GlobalPoint end(GlobalPoint::Cylindrical(Rend,phi,Zend));	
	  GlobalVector step;
	  double Bzdr, Brdz, Bdl;
	  integrate(start, end, np, step, Bzdr, Brdz, Bdl);
	  of << start.x() << " " << start.y() << " " << start.z() << " " 
	     << end.x() << " " << end.y() << " " << end.z() << " " 
	     << step.eta() << " " << step.phi() << " "
	     << Bzdr << " " << Brdz << " " << Bdl	
	     << endl;
	}
	of << endl << endl;
      }
      of.close();
    }
  }

  GlobalPoint pointOnBarrelPlane(float planeR, float localPhi, float eta, int sector){
    float theta = acos(sinh(eta)/cosh(eta));
    float localX = planeR;
    float localY = localX*tan(localPhi);
    float R      = sqrt(localX*localX+localY*localY);
    float Z      = R/tan(theta);

    double sectorPhi=(sector-1)*Geom::pi()/6.;
    double sphi = sin(sectorPhi);
    double cphi = cos(sectorPhi);
    return GlobalPoint(localX*cphi-localY*sphi,localX*sphi+localY*cphi,Z);
  }
  
  // NOTE: Bdl is the transverse bending (i.e. excluding longitudinal bending due to Bphi)  
  void integrate(GlobalPoint start, GlobalPoint end, int np, 
		 GlobalVector& step, double& Bzdr, double& Brdz, double& Bdl) {
    
    step = (end-start)/np;
//     double Dz = end.z()-start.z();  
//     double Dr = end.perp()-start.perp();

    Bzdr=0.;
    Brdz=0.;
    Bdl =0.;
//     double Bzdr2=0;
//     double Brdz2=0;    

    GlobalVector B0 = field->inTesla(start);
    GlobalVector Zaxis(0,0,1);
    
    for (int i=0; i<=np; ++i) {
      GlobalPoint pos = start+step*i;
      // Radial unit vector in the transverse plane
      GlobalVector ur  = (GlobalVector(pos.x(),pos.y(),0)).unit();

      // Unit vector transverse to the R,Z plane and oriented along positive phi
      GlobalVector uphi = (GlobalVector(-pos.y(),pos.x(),0)).unit();

      GlobalVector B;
      {
	static TimingReport::Item & timer = (*TimingReport::current())["B"];
	TimeMe t(timer,false); 
	B = field->inTesla(pos);
      }

// Simple computation "by hand"
//       Bzdr2 = Bzdr2 + B.z();
//       Brdz2 = Brdz2 - ( B.x()*cos(pos.phi()) + B.y()*sin(pos.phi()));
      if (i!=0) {
	GlobalVector averageB = (B+B0); // Use the Average field at start and end of step
	GlobalVector Bz(0,0,averageB.z());
	GlobalVector Br = averageB.dot(ur)*ur;
	
	Bdl =  Bdl +(averageB.cross(step)).dot(uphi);
	Bzdr = Bzdr+(Bz.cross(step)).dot(uphi);
	Brdz = Brdz+(Br.cross(step)).dot(uphi);
      }
      B0=B;
    }

    // Return values are in T*m.
//     Bzdr2 = Bzdr2*Dr/(np+1.)/100.;
//     Brdz2 = Brdz2*Dz/(np+1.)/100.;

    // 2. factor since we used B average, see above
    Brdz = Brdz/200.; 
    Bzdr = Bzdr/200.; 
    Bdl  = Bdl/200.; 

  }   
private:
  const MagneticField* field;
};


DEFINE_FWK_MODULE(intBdl);

