#include "RooFit.h"
#include "Riostream.h"
#include "RooStreamParser.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"
#include "TString.h" 
#include "RooRealVar.h"
# include "../interface/Roo3DMapVar.h"

ClassImp(Roo3DMapVar)



//_____________________________________________________________________________
  Roo3DMapVar::Roo3DMapVar(const char *name, const char *title, const TH3* map,RooWorkspace * w, const RooArgList& dependents,double initVal,double minv,double maxv) : 
    RooAbsReal(name,title),
    _actualVars("actualVars","Variables used by the map",this)
{  
  // Constructor with formula expression and list of input variables

  _map = map;
  _title= std::string(name);

  for (int i=1;i<map->GetNbinsX()+1;++i)
    for (int j=1;j<map->GetNbinsY()+1;++j)
      for (int k=1;k<map->GetNbinsZ()+1;++k) {
	int bin = map->GetBin(i,j,k);
	w->factory(TString::Format("%s_%d[%f,%f,%f]",name,bin,initVal,minv,maxv));
	_allVars.Add(w->var(TString::Format("%s_%d",name,bin)));
	_actualVars.add(*((RooAbsArg*)_allVars.find(TString::Format("%s_%d",name,bin))));

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


RooArgSet& Roo3DMapVar::actualDependents() const {

  _actual.removeAll();
  
  int i ;
  for (i=0 ; i<_allVars.GetSize() ; i++) {
    _actual.add((RooAbsArg&)*_allVars.At(i),kTRUE) ;
  }

  return _actual ;
}





//_____________________________________________________________________________
Roo3DMapVar::Roo3DMapVar(const Roo3DMapVar& other, const char* name) : 
  RooAbsReal(other, name),
  _actualVars("actualVars",this,other._actualVars)
{
  _map = other._map;
  _title= other._title;

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
Roo3DMapVar::~Roo3DMapVar() 
{

}







//_____________________________________________________________________________
Double_t Roo3DMapVar::evaluate() const
{
  // Calculate current value of object from internal formula
  double x  = ((RooAbsReal*)_inputVars.At(0))->getVal();
  double y  = ((RooAbsReal*)_inputVars.At(1))->getVal();
  double z  = ((RooAbsReal*)_inputVars.At(2))->getVal();

  int binx = _map->GetXaxis()->FindBin(x);
  int biny = _map->GetYaxis()->FindBin(y);
  int binz = _map->GetZaxis()->FindBin(z);


  int bin = _map->GetBin(binx,biny,binz);
  RooRealVar * var = (RooRealVar*)_allVars.find(TString::Format("%s_%d",_title.c_str(),bin));
  return var->getVal();
}



//_____________________________________________________________________________
Bool_t Roo3DMapVar::isValidReal(Double_t /*value*/, Bool_t /*printError*/) const 
{
  // Check if given value is valid
  return kTRUE ;
}



//_____________________________________________________________________________
Bool_t Roo3DMapVar::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t /*isRecursive*/)
{
  // Propagate server change information to embedded RooFormula object
  //  changeDependents(newServerList,mustReplaceAll,nameChange);
  return kTRUE;
}



//_____________________________________________________________________________
Double_t Roo3DMapVar::defaultErrorLevel() const 
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



