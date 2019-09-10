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
    int tiss_index() const { return 1 - (parametermapftiss ? 1 : 0); }

    int tau_index() const
    {
        return 2 - (parametermapftiss ? 1 : 0) - (parametermapdelttiss ? 1 : 0)
            + (infertau ? 1 : 0);
    }

    int art_index() const
    {
        return 2 - (parametermapftiss ? 1 : 0) - (parametermapdelttiss ? 1 : 0) + (infertau ? 1 : 0)
            + (inferart ? 1 : 0);
    }

    int t1_index() const
    {
        return 2 - (parametermapftiss ? 1 : 0) - (parametermapdelttiss ? 1 : 0) + (infertau ? 1 : 0)
            + (inferart ? 2 : 0) + (infert1 ? 1 : 0);
    }

    int taub_index() const
    {
        return 2 - (parametermapftiss ? 1 : 0) - (parametermapdelttiss ? 1 : 0) + (infertau ? 1 : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0);
    }

    int t2_index() const
    {
        return 2 - (parametermapftiss ? 1 : 0) - (parametermapdelttiss ? 1 : 0) + (infertau ? 1 : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0) + (infert2 ? 1 : 0);
    }

    int texch_index() const
    {
        return 2 - (parametermapftiss ? 1 : 0) - (parametermapdelttiss ? 1 : 0) + (infertau ? 1 : 0)
            + (inferart ? 2 : 0) + (infert1 ? 2 : 0) + (infertaub ? 1 : 0) + (infert2 ? 2 : 0)
            + (infertexch ? 1 : 0);
    }

    // Scan parameters
    double seqtau;   // Bolus length as set by the sequence
    int repeats;     // Number of repeats
    double t1;       // T1_tis relaxation time of tissue
    double t1b;      // T1_bl relaxation time of blood
    double t2;       // T2_tis relaxation time of tissue
    double t2b;      // T2_bl relaxation time of blood
    double texch;    // T1_bl_tis transfer time between blood and tissue
    double pretisat; // Pre-TI saturation time
    double lambda;   // lambda
    bool casl;       // casl or pasl
    bool grase;      // Data was collected with GRASE-ASL
    NEWMAT::ColumnVector tis; // column vector of t_i
    NEWMAT::ColumnVector tes; // column vector of t_e
    NEWMAT::Real timax;

    // Inference options
    bool infertau;
    bool infertaub;
    bool inferart;
    bool infert1;
    bool infert2;
    bool infertexch;
    bool doard;

    // FIXME not sure what these are for
    bool parametermapftiss;
    bool parametermapdelttiss;

private:
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, multiTEFwdModel> registration;
};
