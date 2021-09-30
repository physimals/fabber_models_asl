/*  fwdmodel_asl_satrecovdualfa.h - Saturation Recovery curve Dual Flip Angle calibration for ASL

 Michael Chappell, IBME & FMRIB Image Analysis Group

 Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"
#include "fabber_core/inference.h"

#include <string>
#include <vector>

class SatrecovDualFAFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    virtual void Initialize(ArgsType &args);
    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void NameParams(std::vector<std::string> &names) const;
    //virtual int NumParams() const { return (LFAon ? 4 : 3); } // This needs to be changed (you need two additional parameters)
    virtual int NumParams() const { return 6; } // You have six parameters to estimate: M0t, T1, A, g, M0_initial_high_FA(fake), M0_initial_low_FA(fake).
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

protected:
    // Constants

    // Lookup the starting indices of the parameters

    // vector indices for the parameters to expereicne ARD
    std::vector<int> ard_index;

    // scan parameters
    int repeats;
    int nphases;
    double t1;
    double slicedt;

    double FAnom;
    double LFA;
    double dti;
    float dg;

    double t_initial_high_FA;
    double t_initial_low_FA;

    bool looklocker;
    bool LFAon;
    bool fixA;

    NEWMAT::ColumnVector tis;
    NEWMAT::Real timax;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, SatrecovDualFAFwdModel> registration;
};
