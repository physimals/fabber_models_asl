/**
 * fwdmodel_asl_satrecov.h - Saturation Recovery curve calibration for ASL
 *
 * Michael Chappell, IBME & FMRIB Image Analysis Group
 *
 * Copyright (C) 2010 University of Oxford
 */

/*  CCOPYRIGHT */

#include "fabber_core/fwdmodel.h"

#include "armawrap/newmat.h"

#include <string>
#include <vector>

/**
 * Fabber model for a saturation recovery calibration sequence
 *
 * Supports normal saturation recovery and Look-Locker
 * acquisitions
 */
class SatrecovFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;

    void Initialize(FabberRunData &args);
    void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result,  const std::string &key="") const;

protected:
    void GetParameterDefaults(std::vector<Parameter> &params) const;

    /** List of TIs */
    std::vector<double> m_tis;

    /** Separation between TIs (assumed constant) */
    double m_dti;

    /** Number of repeats at each TI */
    vector<int> m_repeats;

    /** Time between slices in ms */
    double m_slicedt;

    /** Number of slices per band */
    int m_sliceband;

    /** Number of phases */
    int m_nphases;

    /** Presumed tissue T1 - this is inferred with limited variance */
    double m_t1;

    /** Whether to fix the A parameter */
    bool m_fix_a;

    /** Whether this was a Look-Locker acquisition */
    bool m_look_locker;

    /** Whether the lower flip angle was used */
    bool m_lfa_on;

    /** Flip angle for Look-locker */
    double m_fa;

    /** Lower flip angle for Look-Locker */
    double m_lfa;

    /** DG parameter - currently fixed as constant value of 0.023 */
    float m_dg;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, SatrecovFwdModel> registration;
};
