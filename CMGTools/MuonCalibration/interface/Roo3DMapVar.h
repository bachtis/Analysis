#ifndef ROO_3DMAP_VAR
#define ROO_3DMAP_VAR

#include "RooAbsReal.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooLinkedList.h"
#include "RooListProxy.h"
#include <TH3.h>

class RooArgSet ;

class Roo3DMapVar : public RooAbsReal {
public:
  // Constructors, assignment etc
  inline Roo3DMapVar()   { }
  Roo3DMapVar(const char *name, const char *title,const TH3* map, RooWorkspace* w, const RooArgList& inputs,double initVal = 1.0,double minv=0.5,double maxv=1.5);
    Roo3DMapVar(const Roo3DMapVar& other, const char* name=0);
    virtual TObject* clone(const char* newname) const { return new Roo3DMapVar(*this,newname); }
    virtual ~Roo3DMapVar();

    virtual Double_t defaultErrorLevel() const ;


    inline RooAbsArg* getParameter(const char* name) const { 
      // Return pointer to parameter with given name
      return (RooAbsArg*)_actualVars.find(name) ; 
    }
    inline RooAbsArg* getParameter(Int_t index) const { 
      // Return pointer to parameter at given index
      return (RooAbsArg*)_actualVars.at(index) ; 
    }

    // Dependent management
    RooArgSet& actualDependents() const ;


protected:

  // Function evaluation
  virtual Double_t evaluate() const ;

  // Post-processing of server redirection
  virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;

  virtual Bool_t isValidReal(Double_t value, Bool_t printError) const ;
  RooListProxy _actualVars;
  RooLinkedList _allVars ;     // Actual parameters used by formula engine
  RooLinkedList _inputVars ;     // Actual parameters used by formula engine
  mutable RooArgSet _actual;
  const TH3* _map;
  std::string _title;
  mutable RooArgSet* _nset ;     //! Normalization set to be passed along to contents

  ClassDef(Roo3DMapVar,1) // Real-valued function of other RooAbsArgs calculated by a TFormula expression
};

#endif
