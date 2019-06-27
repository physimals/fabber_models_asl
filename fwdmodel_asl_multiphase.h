/**
 * fwdmodel_asl_multiphase.h
 *
 * Michael Chappell, QuBIc (IBME) & FMRIB Image Analysis Group
 *
 * Copyright (C) 2013 University of Oxford  
 */

/*  CCOPYRIGHT */

#pragma once

#include <fabber_core/fwdmodel.h>

#include <newmat.h>

#include <string>
#include <vector>

class MultiPhaseASLFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;
    void Initialize(ArgsType &args);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void InitVoxelPosterior(MVNDist &posterior) const;
    
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

private:
    // Number of repeats in the data. Phases are assumed to be fastest varying
    int m_repeats;

    // Number of TIs in the data. Data for a given TI/PLD (including all phases and 
    // repeats) is assumed to be blocked within the data
    int m_ntis;

    // Number of parameters inferred for each TI - 2 normally (magnitude + offset)
    // but can be 3 if m_multi_phase_offsets is true
    int m_num_ti_params;

    // Number of phases in data. If phase list not provided these are assumed to 
    // be evenly spaced between 0 and 360
    int m_nphases;

    // If true, each TI/PLD has its own phase offset parameter
    bool m_multi_phase_offsets;

    // Explicitly provided list of phases in degrees
    NEWMAT::ColumnVector m_phases_deg;

    // Inference options
    bool m_incvel;
    bool m_infervel;

    // Name of modulation function
    std::string m_modfn;

    // Fermi function variables
    double m_alpha;
    double m_beta;

    // Modulation function data, used when m_modfn='mat'
    NEWMAT::Matrix m_mod_mat;
    NEWMAT::ColumnVector m_mod_phase;
    NEWMAT::ColumnVector m_mod_v;
    int m_mod_nvelpts;
    double m_mod_vmax;
    double m_mod_vmin;

    // Calculate the modulation function, used when m_modfn='mat'
    double mod_fn(const double inphase, const double v) const;

    // Basic linear interpolation function
    double interp(const NEWMAT::ColumnVector &x, const NEWMAT::ColumnVector &y, const double xi) const;

    // Auto-register with forward model factory
    static FactoryRegistration<FwdModelFactory, MultiPhaseASLFwdModel> registration;
};
