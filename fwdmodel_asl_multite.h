/**
 * fwdmodel_asl_multite.h - Implements the MULTITE model
 * 
 * Michael Chappell, FMRIB Image Analysis Group
 * Josepha Hilmer, FME Bremen, MR Group
 * 
 * Copyright (C) 2007 University of Oxford  
 */

/*  CCOPYRIGHT */
#pragma once

#include <fabber_core/fwdmodel.h>
#include <fabber_core/inference.h>

#include <string>
#include <vector>

class multiTEFwdModel : public FwdModel
{
public:
    std::string ModelVersion() const;
    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;

    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result,
        const std::string &key = "") const;

    static FwdModel *NewInstance();
protected: 
    // Lookup the starting indices of the parameters
    int tiss_index() const 
    { 
        return 1;
    }

    int t1_index() const
    {
        return 2 + (m_infert1 ? 1 : 0);
    }

    int t2_index() const
    {
        return 2 + (m_infert1 ? 2 : 0) + (m_infert2 ? 1 : 0);
    }

    int texch_index() const
    {
        return 2 + (m_infert1 ? 2 : 0) + (m_infert2 ? 2 : 0) + (m_infertexch ? 1 : 0);
    }

    // Scan parameters
    int m_repeats;     // Number of repeats
    double m_t1;       // T1_tis relaxation time of tissue
    double m_t1b;      // T1_bl relaxation time of blood
    double m_t2;       // T2_tis relaxation time of tissue
    double m_t2b;      // T2_bl relaxation time of blood
    double m_texch;    // T1_bl_tis transfer time between blood and tissue
    NEWMAT::ColumnVector m_tis; // column vector of t_i
    NEWMAT::ColumnVector m_tes; // column vector of t_e
    NEWMAT::ColumnVector m_taus; // column vector of bolus durations - same length as tis
    NEWMAT::Real m_timax;

    // Inference options
    bool m_infert1;
    bool m_infert2;
    bool m_infertexch;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, multiTEFwdModel> registration;
};
