#include "RooFit.h"
#include "Riostream.h"
#include "RooStreamParser.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"
#include "TString.h" 
#include "RooRealVar.h"
# include "../interface/Roo2DMapVar.h"

ClassImp(Roo2DMapVar)



//_____________________________________________________________________________
  Roo2DMapVar::Roo2DMapVar(const char *name, const char *title, const TH2* map,RooWorkspace * w, const RooArgList& dependents,double initVal,double minv,double maxv) : 
    RooAbsReal(name,title),
    _actualVars("actualVars","Variables used by the map",this)
{  
  // Constructor with formula expression and list of input variables

  _map = map;
  _title= std::string(title);
  _w =w;
  for (int i=1;i<map->GetNbinsX()+1;++i)
    for (int j=1;j<map->GetNbinsY()+1;++j) {
      int bin = map->GetBin(i,j);
      if (w->var(TString::Format("%s_%d",title,bin))== 0)
	w->factory(TString::Format("%s_%d[%f,%f,%f]",title,bin,initVal,minv,maxv));
      _allVars.Add(w->var(TString::Format("%s_%d",title,bin)));
      _actualVars.add(*((RooAbsArg*)_allVars.find(TString::Format("%s_%d",title,bin))));

      }

  TIterator* iter = dependents.createIterator() ;
  RooAbsArg* arg ;
  while ((arg=(RooAbsArg*)iter->Next())) {
    _allVars.Add(arg) ;
    _inputVars.Add(arg) ;
    _actualVars.add(*arg) ;
  }
  delete iter ;

  
}


RooArgSet& Roo2DMapVar::actualDependents() const {

  _actual.removeAll();
  
  int i ;
  for (i=0 ; i<_allVars.GetSize() ; i++) {
    _actual.add((RooAbsArg&)*_allVars.At(i),kTRUE) ;
  }

  return _actual ;
}





//_____________________________________________________________________________
Roo2DMapVar::Roo2DMapVar(const Roo2DMapVar& other, const char* name) : 
  RooAbsReal(other, name),
  _actualVars("actualVars",this,other._actualVars)
{
  _map = other._map;
  _title= other._title;
  _w = other._w;
  TIterator* iter = other._allVars.MakeIterator() ;
  RooAbsArg* arg ;
  while ((arg=(RooAbsArg*)iter->Next())) {
    _allVars.Add(arg) ;
  }
  delete iter ;

  iter = other._inputVars.MakeIterator() ;
  while ((arg=(RooAbsArg*)iter->Next())) {
    _inputVars.Add(arg) ;
  }
  delete iter ;

  _nset = other._nset;


  // Copy constructor
}



//_____________________________________________________________________________
Roo2DMapVar::~Roo2DMapVar() 
{

}







//_____________________________________________________________________________
Double_t Roo2DMapVar::evaluate() const
{

  // Calculate current value of object from internal formula
  double x  = ((RooAbsReal*)_inputVars.At(0))->getVal();
  double y  = ((RooAbsReal*)_inputVars.At(1))->getVal();



  // Given a point P(x,y), Interpolate approximates the value via bilinear
  // interpolation based on the four nearest bin centers
  // see Wikipedia, Bilinear Interpolation
  // Andy Mastbaum 10/8/2008
  // vaguely based on R.Raja 6-Sep-2008

  Double_t f=0;
  Double_t x1=0,x2=0,y1=0,y2=0;
  Double_t dx,dy;
  int bin_x = _map->GetXaxis()->FindBin(x);
  int bin_y = _map->GetYaxis()->FindBin(y);
  if(bin_x<1 || bin_x>_map->GetNbinsX() || bin_y<1 || bin_y>_map->GetNbinsY()) {
    Error("Interpolate","Cannot interpolate outside histogram domain.");
    //    printf("%f %f %d %d  \n",x,y,bin_x,bin_y);
    return 0;
  }
  Int_t quadrant = 0; // CCW from UR 1,2,3,4
  // which quadrant of the bin (bin_P) are we in?
  dx = _map->GetXaxis()->GetBinUpEdge(bin_x)-x;
  dy = _map->GetYaxis()->GetBinUpEdge(bin_y)-y;
  if (dx<=_map->GetXaxis()->GetBinWidth(bin_x)/2 && dy<=_map->GetYaxis()->GetBinWidth(bin_y)/2)
    quadrant = 1; // upper right
  if (dx>_map->GetXaxis()->GetBinWidth(bin_x)/2 && dy<=_map->GetYaxis()->GetBinWidth(bin_y)/2)
    quadrant = 2; // upper left
  if (dx>_map->GetXaxis()->GetBinWidth(bin_x)/2 && dy>_map->GetYaxis()->GetBinWidth(bin_y)/2)
    quadrant = 3; // lower left
  if (dx<=_map->GetXaxis()->GetBinWidth(bin_x)/2 && dy>_map->GetYaxis()->GetBinWidth(bin_y)/2)
    quadrant = 4; // lower right
  switch(quadrant) {
  case 1:
    x1 = _map->GetXaxis()->GetBinCenter(bin_x);
    y1 = _map->GetYaxis()->GetBinCenter(bin_y);
    x2 = _map->GetXaxis()->GetBinCenter(bin_x+1);
    y2 = _map->GetYaxis()->GetBinCenter(bin_y+1);
    break;
  case 2:
    x1 = _map->GetXaxis()->GetBinCenter(bin_x-1);
    y1 = _map->GetYaxis()->GetBinCenter(bin_y);
    x2 = _map->GetXaxis()->GetBinCenter(bin_x);
    y2 = _map->GetYaxis()->GetBinCenter(bin_y+1);
    break;
  case 3:
    x1 = _map->GetXaxis()->GetBinCenter(bin_x-1);
    y1 = _map->GetYaxis()->GetBinCenter(bin_y-1);
    x2 = _map->GetXaxis()->GetBinCenter(bin_x);
    y2 = _map->GetYaxis()->GetBinCenter(bin_y);
    break;
  case 4:
    x1 = _map->GetXaxis()->GetBinCenter(bin_x);
    y1 = _map->GetYaxis()->GetBinCenter(bin_y-1);
    x2 = _map->GetXaxis()->GetBinCenter(bin_x+1);
    y2 = _map->GetYaxis()->GetBinCenter(bin_y);
    break;
  }
  Int_t bin_x1 = _map->GetXaxis()->FindBin(x1);
  if(bin_x1<1) bin_x1=1;
  Int_t bin_x2 = _map->GetXaxis()->FindBin(x2);
  if(bin_x2>_map->GetNbinsX()) bin_x2=_map->GetNbinsX();
  Int_t bin_y1 = _map->GetYaxis()->FindBin(y1);
  if(bin_y1<1) bin_y1=1;
  Int_t bin_y2 = _map->GetYaxis()->FindBin(y2);
  if(bin_y2>_map->GetNbinsY()) bin_y2=_map->GetNbinsY();
  Int_t bin_q22 = _map->GetBin(bin_x2,bin_y2);
  Int_t bin_q12 = _map->GetBin(bin_x1,bin_y2);
  Int_t bin_q11 = _map->GetBin(bin_x1,bin_y1);
  Int_t bin_q21 = _map->GetBin(bin_x2,bin_y1);
  //  printf("Trying to get the numbers %d %d %d %d \n",bin_q11,bin_q12,bin_q21,bin_q22);
  Double_t q11 = _w->var(TString::Format("%s_%d",_title.c_str(),bin_q11))->getVal();
  Double_t q12 = _w->var(TString::Format("%s_%d",_title.c_str(),bin_q12))->getVal();
  Double_t q21 = _w->var(TString::Format("%s_%d",_title.c_str(),bin_q21))->getVal();
  Double_t q22 = _w->var(TString::Format("%s_%d",_title.c_str(),bin_q22))->getVal();
  Double_t d = 1.0*(x2-x1)*(y2-y1);
  f = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  //  printf("DONE\n");

  return f;
}



//_____________________________________________________________________________
Bool_t Roo2DMapVar::isValidReal(Double_t /*value*/, Bool_t /*printError*/) const 
{
  // Check if given value is valid
  return kTRUE ;
}



//_____________________________________________________________________________
Bool_t Roo2DMapVar::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t /*isRecursive*/)
{
  // Propagate server change information to embedded RooFormula object
  //  changeDependents(newServerList,mustReplaceAll,nameChange);
  return kTRUE;
}



//_____________________________________________________________________________
Double_t Roo2DMapVar::defaultErrorLevel() const 
{
  // Return the default error level for MINUIT error analysis
  // If the formula contains one or more RooNLLVars and 
  // no RooChi2Vars, return the defaultErrorLevel() of
  // RooNLLVar. If the addition contains one ore more RooChi2Vars
  // and no RooNLLVars, return the defaultErrorLevel() of
  // RooChi2Var. If the addition contains neither or both
  // issue a warning message and return a value of 1

  RooAbsReal* nllArg(0) ;
  RooAbsReal* chi2Arg(0) ;

  TIterator* iter = _actualVars.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (dynamic_cast<RooNLLVar*>(arg)) {
      nllArg = (RooAbsReal*)arg ;
    }
    if (dynamic_cast<RooChi2Var*>(arg)) {
      chi2Arg = (RooAbsReal*)arg ;
    }
  }
  delete iter ;

  if (nllArg && !chi2Arg) {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() 
			<< ") Formula contains a RooNLLVar, using its error level" << endl ;
    return nllArg->defaultErrorLevel() ;
  } else if (chi2Arg && !nllArg) {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() 
			<< ") Formula contains a RooChi2Var, using its error level" << endl ;
    return chi2Arg->defaultErrorLevel() ;
  } else if (!nllArg && !chi2Arg) {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() << ") WARNING: "
			<< "Formula contains neither RooNLLVar nor RooChi2Var server, using default level of 1.0" << endl ;
  } else {
    coutI(Minimization) << "RooFormulaVar::defaultErrorLevel(" << GetName() << ") WARNING: "
			<< "Formula contains BOTH RooNLLVar and RooChi2Var server, using default level of 1.0" << endl ;
  }

  return 1.0 ;
}



