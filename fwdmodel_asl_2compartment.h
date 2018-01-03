#pragma once

/*  fwdmodel_asl_2comparment.h - ASL 2-compartment model

    Martin Craig <martin.craig@eng.ox.ac.uk>

    Copyright (C) 2017 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"

#include <string>
#include <vector>

/**
 * ASL 2-compartment model from Parkes and Tofts 2002
 *
 * This model is now incorporated into aslrest, this code remains as a clean-room
 * implementation for comparison purposes only. It may be removed in the future
 */
class ASL2CompartmentModel : public FwdModel
{
public:
    static FwdModel *NewInstance();
    void Initialize(FabberRunData &args);
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string ModelVersion() const;
    std::string GetDescription() const;
    void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key="") const;
    
protected:
    void GetParameterDefaults(std::vector<Parameter> &params) const;

private:
    enum SolutionType {DIST, SLOW, FAST};
    SolutionType m_sol;
    bool m_infer_bat, m_infer_t1, m_infer_t1b, m_infer_t1e;
    double m_f_initial, m_sig0;
    double m_bat, m_t1, m_t1b, m_t1e, m_ie, m_bolusdur, m_ps, m_vb, m_vbw;
    
    /** Time values of the data volumes - must match number of volumes in data */
    std::vector<double> m_times;

    double dM(double t, double f, double bat, double t1, double t1b, double t1e) const;
    double sig(double sig0, double dm) const;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, ASL2CompartmentModel> registration;
};
