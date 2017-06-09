/*  fwdmodel_asl_turboquasar.h - Resting state ASL model for TURBO QUASAR acquisitions

    Moss Zhao and Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class TurboQuasarFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
            ColumnVector& result) const;
  static void ModelUsage();
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return (infertiss?2:0) - (singleti?1:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2+(infertau?1:0)+(infert1?1:0)+(usepve?2:0)):0)+ 2 + (inferart?(artdir?3:4):0) + (calibon?1:0);  

    //return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0); 
  } 

  virtual ~TurboQuasarFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  using FwdModel::SetupARD;
  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  TurboQuasarFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters
  int tiss_index() const {return (infertiss?1:0);} //main tissue parameters: ftiss and delttiss alway come first

  int tau_index() const {  return (infertiss?2:0) + (infertiss?(infertau?1:0):0);  }

  int art_index() const {  return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?1:0); }

  int t1_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?1:0); }
  
  //int inveff_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) +(inferinveff?1:0); }

  //int trailing_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertrailing?1:0); }

  //int taub_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0);}

  int taub_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0);}

  //int R_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0) + (inferart?1:0);}
  int wm_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?1:0); }

  int pv_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?1:0); }

  int disp_index() const {return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?2:0) +1; }

  int crush_index() const {return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?2:0) +2 + (inferart?1:0); }

  int calib_index() const {return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?2:0) +2 + (inferart?(artdir?3:4):0) + (calibon?1:0);}

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;

  // simulation mode
  bool simulation;

  // scan parameters
  double seqtau; //bolus length as set by the sequence
  int repeats;
  int n_bolus;
  double delta_bolus;
  int slice_shifting_factor;
  int delta_ti_gap_factor;
  double t1;
  double t1b;
  double t1wm;
  double lambda;
  double slicedt;
  double pretisat;
  //  bool grase; //to indicate data was collected with GRASE-ASL
  float dti; //TI interval
  float FA; //flip angle

  bool infertiss;
  bool singleti; //specifies that only tissue perfusion should be inferred
  bool infertau;
  bool infertaub;
  bool inferart;
  bool infert1;
  bool inferwm;
  bool usepve;
  bool artdir;
  bool calibon;

  bool onephase;

  string disptype;

  // ard flags
  bool doard;
  bool tissard;
  bool artard;
  bool wmard;

  ColumnVector tis;
  Real timax;
  //ColumnVector crushdir;
  Matrix crushdir;

  // Column vector to save the order of each bolus
  // 1 means bolus and 0 means skip current bolus
  ColumnVector bolus_order;

  // Value to indicate lower limit of bolus duration
  float tau_lowest;

  //kinetic curve functions
  ColumnVector kcblood_nodisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, float deltll, float T_1ll,int n_bolus_total, float delta_bolus, const ColumnVector& bolus_order) const;
  ColumnVector kcblood_gammadisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, float s, float p, float deltll,float T_1ll) const;
  ColumnVector kcblood_gvf(const ColumnVector& tis, float deltblood, float taub, float T_1b, float s, float p, float deltll,float T_1ll) const;
  ColumnVector kcblood_gaussdisp(const ColumnVector& tis, float deltblood, float taub, float T_1b, float sig1, float sig2, float deltll,float T_1ll) const;
  //Tissue
  ColumnVector kctissue_nodisp(const ColumnVector& tis, float delttiss, float tau, float T1_b, float T1_app, float deltll,float T_1ll,int n_bolus_total, float delta_bolus, const ColumnVector& bolus_order) const;
  ColumnVector kctissue_gammadisp(const ColumnVector& tis, float delttiss, float tau, float T1_b, float T1_app, float s, float p, float deltll,float T_1ll) const;
  ColumnVector kctissue_gvf(const ColumnVector& tis, float delttiss, float tau, float T1_b, float T1_app, float s, float p, float deltll,float T_1ll) const;
  ColumnVector kctissue_gaussdisp(const ColumnVector& tis, float delttiss, float tau, float T_1b, float T_1app, float sig1, float sig2, float deltll,float T_1ll) const;

  //useful functions
  float icgf(float a, float x) const;
  float gvf(float t, float s, float p) const;

};
