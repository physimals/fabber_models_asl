#pragma once

/**
 * fwdmodel_asl_velocityselective.h Velocity selective ASL model
 *
 * Moss Zhao <mosszhao@stanford.edu>
 *
 * Copyright (C) 2020 Stanford University
 */

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

class VelocitySelectiveFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    // Virtual function overrides
    virtual void Initialize(ArgsType &args);
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;

    virtual void NameParams(std::vector<std::string> &names) const;

    // For now we only estimate CBF
    virtual int NumParams() const
    {
        return (infertiss ? 1 : 0); // This indicates that we only estimate CBF (one) parameter
    }

    virtual ~VelocitySelectiveFwdModel() {}
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;

protected:
    // Lookup the starting indices of the parameters
    // If more parameters are included, check out turboquasar model for more sophisticated inplementations
    int tiss_index() const
    {
        return (infertiss ? 1 : 0);
    }

    // scan parameters
    double tau_initial; // bolus duration
    int repeats_initial; // Number of repeats. This is no NEX in GE. GE adds all signal in k-space
    double predelay_initial; // User set variatio typically found in one of the CV values in the scanner. Range should be around 2 secs
    double t1b_initial; // We only need T1 of blood in VS ASL model
    double slicedt_initial; // GE uses spiral readout so this value is zero by default
    double t2b_initial;
    double te_initial;

    // Decision to infer tissue component. Default true.
    bool infertiss;

    NEWMAT::ColumnVector tis;  // Currently only support a single TI VS ASL

    // kinetic curve functions
    NEWMAT::ColumnVector kctissue_model(double ftiss, const NEWMAT::ColumnVector &tis, double tau, double T_1b, double predelay, double te, double T_2b) const;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, VelocitySelectiveFwdModel> registration;
};
