/*  fwdmodel_asl_quasar.h - Resting state ASL model for QUASAR acquisitions

 Michael Chappell, IBME & FMRIB Image Analysis Group

 Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */
#pragma ONCE

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

class QuasarFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    virtual void Initialize(ArgsType &args);
    virtual void Evaluate(const NEWMAT::ColumnVector &params,
        NEWMAT::ColumnVector &result) const;
    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void NameParams(std::vector<std::string> &names) const;
    virtual int NumParams() const
    {
        return (infertiss ? 2 : 0) - (singleti ? 1 : 0)
            + (infertiss ? (infertau ? 1 : 0) : 0) + (inferart ? 2 : 0)
            + (infert1 ? 2 : 0) + (infertaub ? 1 : 0)
            + (inferwm ? (2 + (infertau ? 1 : 0) + (infert1 ? 1 : 0)
                             + (usepve ? 2 : 0))
                       : 0)
            + 2 + (inferart ? (artdir ? 3 : 4) : 0)
            + (calibon ? 1 : 0);
   }

    virtual ~QuasarFwdModel()
    {
        return;
    }

    virtual void HardcodedInitialDists(MVNDist &prior,
        MVNDist &posterior) const;

    using FwdModel::SetupARD;
    virtual void SetupARD(const MVNDist &posterior, MVNDist &prior,
        double &Fard);
    virtual void UpdateARD(const MVNDist &posterior, MVNDist &prior,
        double &Fard) const;

protected:
    // Constants

    // Lookup the starting indices of the parameters
    int tiss_index() const
    {
        return (infertiss ? 1 : 0);
    } //main tissue parameters: ftiss and delttiss alway come first

    int tau_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0);
    }

    int art_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 1 : 0);
    }

    int t1_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 1 : 0);
    }

    int taub_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0);
    }

    int wm_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0)
            + (inferwm ? 1 : 0);
    }

    int pv_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0)
            + (inferwm ? (2 + (infertau ? 1 : 0) + (infert1 ? 1 : 0)) : 0)
            + (usepve ? 1 : 0);
    }

    int disp_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0)
            + (inferwm ? (2 + (infertau ? 1 : 0) + (infert1 ? 1 : 0)) : 0)
            + (usepve ? 2 : 0) + 1;
    }

    int crush_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0)
            + (inferwm ? (2 + (infertau ? 1 : 0) + (infert1 ? 1 : 0)) : 0)
            + (usepve ? 2 : 0) + 2 + (inferart ? 1 : 0);
    }

    int calib_index() const
    {
        return (infertiss ? 2 : 0) + (infertiss ? (infertau ? 1 : 0) : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0)
            + (inferwm ? (2 + (infertau ? 1 : 0) + (infert1 ? 1 : 0)) : 0)
            + (usepve ? 2 : 0) + 2 + (inferart ? (artdir ? 3 : 4) : 0)
            + (calibon ? 1 : 0);
    }

    // vector indices for the parameters to expereicne ARD
    vector<int> ard_index;

    // scan parameters
    double seqtau; //bolus length as set by the sequence
    int repeats;
    double t1;
    double t1b;
    double t1wm;
    double lambda;
    double slicedt;
    double pretisat;
    
    float dti; //TI interval
    float FA;  //flip angle

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

    std::string disptype;

    // ard flags
    bool doard;
    bool tissard;
    bool artard;
    bool wmard;

    NEWMAT::ColumnVector tis;
    double timax;
    NEWMAT::Matrix crushdir;

    //kinetic curve functions
    NEWMAT::ColumnVector kcblood_nodisp(const NEWMAT::ColumnVector &tis, float deltblood,
        float taub, float T_1b, float deltll, float T_1ll) const;
    NEWMAT::ColumnVector kcblood_gammadisp(const NEWMAT::ColumnVector &tis, float deltblood,
        float taub, float T_1b, float s, float p, float deltll,
        float T_1ll) const;
    NEWMAT::ColumnVector kcblood_gvf(const NEWMAT::ColumnVector &tis, float deltblood,
        float taub, float T_1b, float s, float p, float deltll,
        float T_1ll) const;
    NEWMAT::ColumnVector kcblood_gaussdisp(const NEWMAT::ColumnVector &tis, float deltblood,
        float taub, float T_1b, float sig1, float sig2, float deltll,
        float T_1ll) const;
    //Tissue
    NEWMAT::ColumnVector kctissue_nodisp(const NEWMAT::ColumnVector &tis, float delttiss,
        float tau, float T1_b, float T1_app, float deltll,
        float T_1ll) const;
    NEWMAT::ColumnVector kctissue_gammadisp(const NEWMAT::ColumnVector &tis, float delttiss,
        float tau, float T1_b, float T1_app, float s, float p, float deltll,
        float T_1ll) const;
    NEWMAT::ColumnVector kctissue_gvf(const NEWMAT::ColumnVector &tis, float delttiss,
        float tau, float T1_b, float T1_app, float s, float p, float deltll,
        float T_1ll) const;
    NEWMAT::ColumnVector kctissue_gaussdisp(const NEWMAT::ColumnVector &tis, float delttiss,
        float tau, float T_1b, float T_1app, float sig1, float sig2,
        float deltll, float T_1ll) const;

    //useful functions
    float icgf(float a, float x) const;
    float gvf(float t, float s, float p) const;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, QuasarFwdModel> registration;
};
