/**
 * fwdmodel_asl_multite.cc - Implements the MULTITE ASL model
 * 
 * Michael Chappell, FMRIB Image Analysis Group
 * Josepha Hilmer, FME Bremen, MR Group
 * 
 * Copyright (C) 2007 University of Oxford  
 */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_multite.h"

#include <fabber_core/easylog.h>
#include <fabber_core/tools.h>
#include <fabber_core/priors.h>

#include <newmatio.h>

#include <miscmaths/miscprob.h>

#include <string>
#include <vector>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, multiTEFwdModel> multiTEFwdModel::registration("asl_multite");

FwdModel *multiTEFwdModel::NewInstance() { return new multiTEFwdModel(); }

static OptionSpec OPTIONS[] = {
    { "repeats", OPT_BOOL, "Number of repeats in data", OPT_REQ, "" },
    { "ti<n>", OPT_BOOL, "Sequence of inversion times in seconds (e.g. --ti1=1.4 --ti2=1.8, etc)", OPT_REQ, "" },
    { "te<n>", OPT_BOOL, "Sequence of TE times in seconds (e.g. --te1=1.4 --te2=1.8, etc)", OPT_REQ, "" },
    { "tau", OPT_FLOAT, "Temporal bolus lengths (s)", OPT_NONREQ, "1.0" },
    { "tau<n>", OPT_FLOAT, "Temporal bolus lengths (one per TI) (s)", OPT_NONREQ, "1.0" },
    { "t1", OPT_FLOAT, "T1 of tissue(s)", OPT_NONREQ, "1.0" },
    { "t1b", OPT_FLOAT, "T1 of blood (s)", OPT_NONREQ, "1.2" },
    { "t2", OPT_FLOAT, "T2 of tissue(s)", OPT_NONREQ, "0.030" },
    { "t2b", OPT_FLOAT, "T2 of blood (s)", OPT_NONREQ, "0.120" },
    { "exch2", OPT_FLOAT, "Exchange time", OPT_NONREQ, "0.1" },
    { "infert1", OPT_BOOL, "Infer T1 values", OPT_NONREQ, "" },
    { "infert2", OPT_BOOL, "Infer T2 values", OPT_NONREQ, "" },
    { "infertexch", OPT_BOOL, "Infer exchange time", OPT_NONREQ, "" },
    { "tauboff", OPT_BOOL, "Forces the inference of arterial bolus off", OPT_NONREQ, "" },
    { "" },
};

void multiTEFwdModel::GetOptions(std::vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string multiTEFwdModel::GetDescription() const { return "Model for Multi-TE ASL data"; }

std::string multiTEFwdModel::ModelVersion() const
{
    std::string version = "fwdmodel_asl_multite.cc";
#ifdef GIT_SHA1
    version += std::string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += std::string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void multiTEFwdModel::Initialize(FabberRunData &rundata)
{
    // Acquisition data
    m_repeats = rundata.GetIntDefault("repeats", 1);

    // TIs
    std::vector<double> tis_list = rundata.GetDoubleList("ti", 0.0);
    m_tis.ReSize(tis_list.size());
    for (unsigned int i=0; i<tis_list.size(); i++) 
    {
      m_tis(i+1) = tis_list[i];
    }
    // The final TI
    m_timax = m_tis.Maximum();

    std::vector<double> taus_list = rundata.GetDoubleList("tau", 0.0);
    m_taus.ReSize(tis_list.size());
    if (taus_list.size() < 2) {
        double tau = 1.0;
        if (taus_list.size() == 1) tau = taus_list[0];
        m_taus = tau;
    }
    else if (taus_list.size() == tis_list.size()) 
    {
        for (unsigned int i=0; i<taus_list.size(); i++) 
        {
            m_taus(i+1) = taus_list[i];
        }
    }
    else 
    {
        throw InvalidOptionValue("Incorrect number of taus specified - should match the number of TIs", 
                                 stringify(taus_list.size()), stringify(tis_list.size()));
    }

    // TEs FIXME original code had default of 0.03 for first
    std::vector<double> tes_list = rundata.GetDoubleList("te", 0.0);
    m_tes.ReSize(tes_list.size());
    for (unsigned int i=0; i<tes_list.size(); i++) 
    {
      m_tes(i+1) = tes_list[i];
    }

    // Physiological parameters
    // Defaults matched to the ASLREST model which is
    // different from the original model. These are generally
    // overridden when called by OXASL
    m_t1 = rundata.GetDoubleDefault("t1", 1.3);
    m_t1b = rundata.GetDoubleDefault("t1b", 1.65);
    m_t2 = rundata.GetDoubleDefault("t2", 0.050);
    m_t2b = rundata.GetDoubleDefault("t2b", 0.150);
    m_texch = rundata.GetDoubleDefault("exch2", 0.1);
    m_bat = rundata.GetDoubleDefault("bat", 0.7)
    m_batsd = rundata.GetDoubleDefault("batsd", 0.316)

    // Inference options
    m_infert1 = rundata.GetBool("infert1");   // infer on T1 values
    m_infert2 = rundata.GetBool("infert2");   // infer on T2 values
    m_infertexch = rundata.GetBool("infertexch");

    // Write information about the parameters to the log
    LOG << "Exchange Time: " << m_texch << endl; // TH
    LOG << "Inference using Buxton Kinetic Curve model" << endl;
    LOG << "    Data parameters: repeats = " << m_repeats << ", t1 = " << m_t1 << ", t1b = " << m_t1b;
    LOG << ", bolus lengths (tau) = " << m_taus.t() << endl;
    if (m_infert1)
        LOG << "Infering on T1 values " << endl;
    if (m_infert2)
        LOG << "Infering on T2 values " << endl;
      
    LOG << "TIs: ";
    for (int i = 1; i <= m_tis.Nrows(); i++)
        LOG << m_tis(i) << " ";
    LOG << endl;

    LOG << "TE: ";
    for (int i = 1; i <= m_tes.Nrows(); i++)
        LOG << m_tes(i) << " ";
    LOG << endl;
}

void multiTEFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();
    int p = 0;

    params.push_back(Parameter(p++, "ftiss", DistParams(0, 1e12), DistParams(0.1, 1.0)));
    params.push_back(Parameter(p++, "delttiss", DistParams(m_bat, m_batsd**2), DistParams(m_bat, m_batsd**2)));

    if (m_infert1)
    {
        params.push_back(Parameter(p++, "T_1", DistParams(m_t1, 0.1), DistParams(m_t1, 0.1)));
        params.push_back(Parameter(p++, "t1b", DistParams(m_t1b, 0.1), DistParams(m_t1b, 0.1)));
    }
    
    if (m_infert2)
    {
        params.push_back(Parameter(p++, "T_2", DistParams(m_t2, 1.0), DistParams(m_t2, 1.0)));
        params.push_back(Parameter(p++, "T_2b", DistParams(m_t2b, 1.0), DistParams(m_t2b, 1.0)));
    }

    if (m_infertexch)
    {
        params.push_back(Parameter(p++, "T_exch", DistParams(m_texch, 0.1), DistParams(m_texch, 0.1)));
    }
}

void multiTEFwdModel::EvaluateModel(const ColumnVector &params, ColumnVector &result,
    const std::string &key) const // Evaluate the forward model
{
    // Negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= params.Nrows(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }

    // Sensible limits on transit time
    if (params(tiss_index() + 1) > m_timax - 0.2)
    {
        paramcpy(tiss_index() + 1) = m_timax - 0.2;
    }

    // Parameters - extract and give sensible names
    double f = paramcpy(tiss_index());        // perfusion in ml/100g/min
    double att = paramcpy(tiss_index() + 1);  // arterial transit time in s

    // Set defaults to be used if the parameters are fixed (not inferred)
    double t1 = m_t1;                        // longitudinal relaxation time in tissue in s
    double t1b = m_t1b;                      // longitudinal relaxation time in blood in s
    double t2 = m_t2;                        // transversal relaxation time in tissue in s
    double t2b = m_t2b;                      // transversal relaxation time in blood in s
    double texch = m_texch;                  // characteristic transfer time in s

    if (m_infert1)
    {  
        // T1 cannot get too close to zero
        t1 = params(t1_index());
        t1b = params(t1_index() + 1);
        if (t1 < 1e-2)
            t1 = 1e-2;
        if (t1b < 1e-2)
            t1b = 1e-2;
    }
    
    if (m_infert2)
    { 
        // T2 cannot get too close to zero
        t2b = params(t2_index() + 1);
        t2 = params(t2_index());
        if (t2b < 1e-2)
            t2b = 1e-2;
        if (t2 < 1e-2)
            t2 = 1e-2;
    }
   
    if (m_infertexch)
    {
        // Texch cannot get too close to zero
        texch = params(texch_index());
        if (texch < 1e-2)
            texch = 1e-2;
    }

    // Calculation of delta M using the two compartment model (PCASL) with T1 and T2 decay using
    // only the B-part
    // the functions are calculated with Matlab

    // Matrix with #tis rows and #tes columns
    Matrix delta_M_final(m_tis.Nrows(), m_tes.Nrows());

    for (int j = 1; j < m_tis.Nrows() + 1; j++)
    {
        double tau = m_taus(j);
        double ti = m_tis(j);
        if ((0 < ti) && (ti < att))
        {
            for (int k = 1; k < m_tes.Nrows() + 1; k++)
            {
                delta_M_final(j, k) = 0;
            }
        }

        else if ((att <= ti) && (ti <= (att + tau)))
        {
            for (int k = 1; k < m_tes.Nrows() + 1; k++)
            {
                double te = m_tes(k);
                delta_M_final(j, k) = (f
                    * (exp(-te / t2)
                              * (t1 * exp(-ti / t1) * exp(-att / t1b)
                                        * (exp(ti / t1) - exp(att / t1))
                                    - (t1 * texch * exp(-ti / t1)
                                          * exp(-ti / texch) * exp(-att / t1b)
                                          * (exp(ti / t1) * exp(ti / texch)
                                                - exp(att / t1) * exp(att / texch)))
                                        / (t1 + texch))
                          - exp(-te / t2b) * exp(-te / texch)
                              * exp(-ti / t1b) * exp(-ti / texch)
                              * ((t1b * texch * exp(att / texch)) / (t1b + texch)
                                    - (t1b * texch * exp(ti / t1b)
                                          * exp(ti / texch) * exp(-att / t1b))
                                        / (t1b + texch))
                          + exp(-te / t2) * exp(-ti / t1b)
                              * exp(-ti / texch) * (exp(-te / texch) - 1.0)
                              * ((t1b * texch * exp(att / texch)) / (t1b + texch)
                                    - (t1b * texch * exp(ti / t1b)
                                          * exp(ti / texch) * exp(-att / t1b))
                                        / (t1b + texch)))
                    * 2.0);
            }
        }
        else
        {
            for (int k = 1; k < m_tes.Nrows() + 1; k++)
            {
                double te = m_tes(k);
                delta_M_final(j, k) = (f
                    * (exp(-te / t2)
                              * (t1 * exp(-ti / t1) * exp(-att / t1b)
                                        * exp(att / t1) * (exp(tau / t1) - 1.0)
                                    + (t1 * texch * exp(-ti / t1)
                                          * exp(-ti / texch) * exp(-att / t1b)
                                          * (exp(att / t1) * exp(att / texch)
                                                - exp((att + tau) / t1)
                                                    * exp((att + tau) / texch)))
                                        / (t1 + texch))
                          - (t1b * texch * exp(-te / t2)
                                * exp(-ti / t1b) * exp(-ti / texch)
                                * exp(att / texch) * (exp(-te / texch) - 1.0)
                                * (exp(tau / t1b) * exp(tau / texch) - 1.0))
                              / (t1b + texch)
                          + (t1b * texch * exp(-te / t2b)
                                * exp(-te / texch) * exp(-ti / t1b)
                                * exp(-ti / texch) * exp(att / texch)
                                * (exp(tau / t1b) * exp(tau / texch) - 1.0))
                              / (t1b + texch))
                    * 2.0);
            }
        }
    }

    // Generating the output
    
    // Resize the ColumnVector signal to #m_tis * #tes
    ColumnVector signal;
    signal.ReSize(m_tis.Nrows() * m_tes.Nrows());

    // Fill in the results of the matrix delta_M_final
    for (int j = 1; j <= m_tis.Nrows(); j++)
    {
        for (int k = 1; k <= m_tes.Nrows(); k++)
        {
            signal((j - 1) * m_tes.Nrows() + k) = delta_M_final(j, k);
        }
    }

    // Resize the ColumnVector result to #tis * #tes * #repeats
    result.ReSize(m_tis.Nrows() * m_tes.Nrows() * m_repeats); 
    result = signal;

    // Concatenate signal repeats
    for (int j = 2; j < m_repeats + 1; j++)
    {
        result &= signal;
    }
}
