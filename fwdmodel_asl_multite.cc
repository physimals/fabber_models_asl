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
    { "casl", OPT_BOOL, "Use CASL (or pCASL) preparation rather than PASL", OPT_NONREQ, "" },
    { "grase", OPT_BOOL, "*DEPRECATAED* (data collected using GRASE-ASL: same as --pretissat=0.1)", OPT_NONREQ, "" },
    { "pretisat", OPT_FLOAT, "Blood is saturated a specific time (seconds) before TI image acquired", OPT_NONREQ, "0" },
    { "tau", OPT_FLOAT, "Temporal bolus lengths (s)", OPT_NONREQ, "1.0" },
    { "tau<n>", OPT_FLOAT, "Temporal bolus lengths (one per TI) (s)", OPT_NONREQ, "1.0" },
    { "t1", OPT_FLOAT, "T1 of tissue(s)", OPT_NONREQ, "1.0" },
    { "t1b", OPT_FLOAT, "T1 of blood (s)", OPT_NONREQ, "1.2" },
    { "t2", OPT_FLOAT, "T2 of tissue(s)", OPT_NONREQ, "0.030" },
    { "t2b", OPT_FLOAT, "T2 of blood (s)", OPT_NONREQ, "0.120" },
    { "exch2", OPT_FLOAT, "Exchange time", OPT_NONREQ, "0.1" },
    { "infertau", OPT_BOOL, "Infer bolus length", OPT_NONREQ, "" },
    { "inferart", OPT_BOOL, "Infer arterial component", OPT_NONREQ, "" },
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
    repeats = rundata.GetIntDefault("repeats", 1);

    // Saturation of the bolus a fixed time pre TI measurement
    pretisat = rundata.GetDoubleDefault("pretisat", 0.0);

    // DEPRECATED If data has come from the GRASE-ASL sequence apply pretisat of 0.1s
    bool grase = rundata.GetBool("grase");
    if (grase)
        pretisat = 0.1;

    // Set if the data is CASL or PASL (default)
    casl = rundata.GetBool("casl");

    // TIs
    std::vector<double> tis_list = rundata.GetDoubleList("ti", 0.0);
    tis.ReSize(tis_list.size());
    for (unsigned int i=0; i<tis_list.size(); i++) 
    {
      tis(i+1) = tis_list[i];
    }
    // The final TI
    timax = tis.Maximum();

    std::vector<double> taus_list = rundata.GetDoubleList("tau", 0.0);
    taus.ReSize(tis_list.size());
    if (taus_list.size() < 2) {
        double tau = 1.0;
        if (taus_list.size() == 1) tau = taus_list[0];
        taus = tau;
    }
    else if (taus_list.size() == tis_list.size()) 
    {
        for (unsigned int i=0; i<taus_list.size(); i++) 
        {
            taus(i+1) = taus_list[i];
        }
    }
    else 
    {
        throw InvalidOptionValue("Incorrect number of taus specified - should match the number of TIs", 
                                 stringify(taus_list.size()), stringify(tis_list.size()));
    }

    // TEs FIXME original code had default of 0.03 for first
    std::vector<double> tes_list = rundata.GetDoubleList("te", 0.0);
    tes.ReSize(tes_list.size());
    for (unsigned int i=0; i<tes_list.size(); i++) 
    {
      tes(i+1) = tes_list[i];
    }

    // Physiological parameters
    t1 = rundata.GetDoubleDefault("t1", 1.0);
    t1b = rundata.GetDoubleDefault("t1b", 1.2);
    t2 = rundata.GetDoubleDefault("t2", 0.030);
    t2b = rundata.GetDoubleDefault("t2b", 0.120);
    lambda = rundata.GetDoubleDefault("lambda", 0.9); // add 01/17/17 J.H.
    texch = rundata.GetDoubleDefault("exch2", 0.1);

    // parametermaps?
    // parametermaps given for ftiss and delttiss FIXME what is this for?
    parametermapftiss = rundata.ReadBool("parametermap-ftiss");
    parametermapdelttiss = rundata.ReadBool("parametermap-delttiss"); 

    // Inference options
    infertau = rundata.GetBool("infertau"); // infer on bolus length?
    infert1 = rundata.GetBool("infert1");   // infer on T1 values?
    inferart = rundata.GetBool("inferart"); // infer on arterial compartment?
    infert2 = rundata.GetBool("infert2");   // infer on T2 values?
    infertexch = rundata.GetBool("infertexch");

    // Forces the inference of arterial bolus off
    bool tauboff = rundata.GetBool("tauboff");
    infertaub = (inferart && infertau && !tauboff);

    // Write information about the parameters to the log
    LOG << "Exchange Time: " << texch << endl; // TH
    LOG << "Inference using Buxton Kinetic Curve model" << endl;
    if (!casl)
        LOG << "Data being analysed using PASL inversion profile" << endl;
    if (casl)
        LOG << "Data being analysed using CASL inversion profile" << endl;
    if (pretisat > 0)
        LOG << "Saturation of" << pretisat << "s before TI has been specified" << endl;
    if (grase)
        LOG << "Using pre TI saturation of 0.1 for GRASE-ASL sequence" << endl;
    //if (calib)
    //    LOG << "Input data is in physiological units, using estimated CBF in T_1app calculation"
    //        << endl;
    LOG << "    Data parameters: repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
    LOG << ", bolus lengths (tau) = " << taus.t() << endl;
    if (infertau)
        LOG << "Infering on bolus length " << endl;
    if (inferart)
        LOG << "Infering on artertial compartment " << endl;
    if (doard)
        LOG << "ARD has been set on arterial compartment " << endl;
    if (infert1)
        LOG << "Infering on T1 values " << endl;
    if (infert2)
        LOG << "Infering on T2 values " << endl;
      
    LOG << "TIs: ";
    for (int i = 1; i <= tis.Nrows(); i++)
        LOG << tis(i) << " ";
    LOG << endl;

    LOG << "TE: ";
    for (int i = 1; i <= tes.Nrows(); i++)
        LOG << tes(i) << " ";
    LOG << endl;
}

void multiTEFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();
    int p = 0;

    if (!parametermapftiss && !parametermapdelttiss)
    {
        params.push_back(Parameter(p++, "ftiss", DistParams(0, 1e12), DistParams(0.1, 1.0)));
        // changed from 0.7 to 0.6 (rat) 28.09.2015
        params.push_back(Parameter(p++, "delttiss", DistParams(0.6, 10), DistParams(0.6, 10)));
    }
    //if (infertau)
    //{
    //    params.push_back(Parameter(p++, "tautiss", DistParams(seqtau, 100), DistParams(seqtau, 100)));
    //}
    if (inferart)
    {
        if (doard) 
        {
            params.push_back(Parameter(p++, "fblood", DistParams(0, 1e12), DistParams(0.1, 1.0), PRIOR_ARD));
        }
        else 
        {
            params.push_back(Parameter(p++, "fblood", DistParams(0, 1e12), DistParams(0.1, 1.0), PRIOR_NORMAL));
        }
        params.push_back(Parameter(p++, "deltblood", DistParams(0.5, 10), DistParams(0.5, 10)));
    }
    if (infert1)
    {
        params.push_back(Parameter(p++, "T_1", DistParams(t1, 0.1), DistParams(t1, 0.1)));
        params.push_back(Parameter(p++, "T_1b", DistParams(t1b, 0.1), DistParams(t1b, 0.1)));
    }
    //if (infertaub) {
    //    params.push_back(Parameter(p++, "taublood", DistParams(seqtau, 10), DistParams(seqtau, 10)));
    //}
    if (infert2)
    {
        params.push_back(Parameter(p++, "T_2", DistParams(t2, 1.0), DistParams(t2, 1.0)));
        params.push_back(Parameter(p++, "T_2b", DistParams(t2b, 1.0), DistParams(t2b, 1.0)));
    }
    if (infertexch)
    {
        params.push_back(Parameter(p++, "T_exch", DistParams(texch, 0.1), DistParams(texch, 0.1)));
    }
}

void multiTEFwdModel::EvaluateModel(const ColumnVector &params, ColumnVector &result,
    const std::string &key) const // Evaluate the forward model
{
    // Ensure that values are reasonable

    // Negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= NumParams(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }

    // Sensible limits on transit times
    if (!parametermapftiss && !parametermapdelttiss)
    {
        if (params(tiss_index() + 1) > timax - 0.2)
        {
            paramcpy(tiss_index() + 1) = timax - 0.2;
        }
    }
    if (inferart)
    {
        if (params(art_index() + 1) > timax - 0.2)
        {
            paramcpy(art_index() + 1) = timax - 0.2;
        }
    }

    // Parameters - extract and give sensible names

    // Set defaults to be used if the parameters are fixed (not inferred)
    float ftiss;
    float delttiss;
    float fblood = 0;
    float deltblood = 0;
    float T_1 = t1;
    float T_1b = t1b;
    float T_2 = t2;
    float T_2b = t2b;
    float T_exch = texch;

    // exact parameter values
    if (!parametermapftiss && !parametermapdelttiss)
    {
        ftiss = paramcpy(tiss_index());
        delttiss = paramcpy(tiss_index() + 1);
    }
    else
    {
        // FIXME
        // ftiss = parametermapftiss(1);
        // delttiss = parametermapdelttiss(1);
    }

    // FIXME
    //if (infertau)
   // {
   // }

    if (inferart)
    {
        fblood = paramcpy(art_index());
        deltblood = paramcpy(art_index() + 1);
    }

    if (infert1)
    {  
        // T1 cannot be zero!
        // FIXME Why here params and not paramscpy
        T_1 = params(t1_index());
        T_1b = params(t1_index() + 1);
        if (T_1 < 0.01)
            T_1 = 0.01;
        if (T_1b < 0.01)
            T_1b = 0.01;
    }
    
    if (infert2)
    { 
        // T2 cannot be zero!
        T_2b = params(t2_index() + 1);
        T_2 = params(t2_index());
        if (T_2b < 1e-2)
            T_2b = 1e-2;
        if (T_2 < 1e-2)
            T_2 = 1e-2;
    }
   
    if (infertexch)
    {
        T_exch = params(texch_index());
        if (T_exch < 1e-2)
            T_exch = 1e-2;
    }
    
    // October 2016: Josepha Hilmer
    // Compatibility of parameter names
    float M_0_bl = 1;          // equilibrium magnetization blood in A/m
    float M_0_tis = M_0_bl;    // equilibrium magnetization tissue in A/m
    float alpha = 1;

    double f = ftiss;          // perfusion in ml/100g/min
    double t_tra = delttiss;   // arterial transit time in s
    double T1_bl = T_1b;       // longitudinal relaxation time in blood in s
    double T1_tis = T_1;       // longitudinal relaxation time in tissue in s
    double T2_bl = T_2b;       // transversal relaxation time in blood in s
    double T2_tis = T_2;       // transversal relaxation time in tissue in s
    double T1_bl_tis = T_exch; // characteristic transfer time in s

    // Calculation of delta M using the two compartment model (pcasl) with T1 and T2 decay using
    // only the B-part
    // the functions are calculated with Matlab

    // Matrix with #tis rows and #tes columns
    Matrix delta_M_final(tis.Nrows(), tes.Nrows());

    for (int j = 1; j < tis.Nrows() + 1; j++)
    {
        double tau = taus(j);
        if ((0 < tis(j)) && (tis(j) < t_tra))
        {
            for (int k = 1; k < tes.Nrows() + 1; k++)
            {
                delta_M_final(j, k) = 0;
            }
        }

        else if ((t_tra <= tis(j)) && (tis(j) <= (t_tra + tau)))
        {
            for (int k = 1; k < tes.Nrows() + 1; k++)
            {
                delta_M_final(j, k) = (alpha * f
                    * (exp(-tes(k) / T2_tis)
                              * (M_0_tis * T1_tis * exp(-tis(j) / T1_tis) * exp(-t_tra / T1_bl)
                                        * (exp(tis(j) / T1_tis) - exp(t_tra / T1_tis))
                                    - (M_0_tis * T1_tis * T1_bl_tis * exp(-tis(j) / T1_tis)
                                          * exp(-tis(j) / T1_bl_tis) * exp(-t_tra / T1_bl)
                                          * (exp(tis(j) / T1_tis) * exp(tis(j) / T1_bl_tis)
                                                - exp(t_tra / T1_tis) * exp(t_tra / T1_bl_tis)))
                                        / (T1_tis + T1_bl_tis))
                          - M_0_bl * exp(-tes(k) / T2_bl) * exp(-tes(k) / T1_bl_tis)
                              * exp(-tis(j) / T1_bl) * exp(-tis(j) / T1_bl_tis)
                              * ((T1_bl * T1_bl_tis * exp(t_tra / T1_bl_tis)) / (T1_bl + T1_bl_tis)
                                    - (T1_bl * T1_bl_tis * exp(tis(j) / T1_bl)
                                          * exp(tis(j) / T1_bl_tis) * exp(-t_tra / T1_bl))
                                        / (T1_bl + T1_bl_tis))
                          + M_0_bl * exp(-tes(k) / T2_tis) * exp(-tis(j) / T1_bl)
                              * exp(-tis(j) / T1_bl_tis) * (exp(-tes(k) / T1_bl_tis) - 1.0)
                              * ((T1_bl * T1_bl_tis * exp(t_tra / T1_bl_tis)) / (T1_bl + T1_bl_tis)
                                    - (T1_bl * T1_bl_tis * exp(tis(j) / T1_bl)
                                          * exp(tis(j) / T1_bl_tis) * exp(-t_tra / T1_bl))
                                        / (T1_bl + T1_bl_tis)))
                    * 2.0);
            }
        }
        else
        {
            for (int k = 1; k < tes.Nrows() + 1; k++)
            {
                delta_M_final(j, k) = (alpha * f
                    * (exp(-tes(k) / T2_tis)
                              * (M_0_tis * T1_tis * exp(-tis(j) / T1_tis) * exp(-t_tra / T1_bl)
                                        * exp(t_tra / T1_tis) * (exp(tau / T1_tis) - 1.0)
                                    + (M_0_tis * T1_tis * T1_bl_tis * exp(-tis(j) / T1_tis)
                                          * exp(-tis(j) / T1_bl_tis) * exp(-t_tra / T1_bl)
                                          * (exp(t_tra / T1_tis) * exp(t_tra / T1_bl_tis)
                                                - exp((t_tra + tau) / T1_tis)
                                                    * exp((t_tra + tau) / T1_bl_tis)))
                                        / (T1_tis + T1_bl_tis))
                          - (M_0_bl * T1_bl * T1_bl_tis * exp(-tes(k) / T2_tis)
                                * exp(-tis(j) / T1_bl) * exp(-tis(j) / T1_bl_tis)
                                * exp(t_tra / T1_bl_tis) * (exp(-tes(k) / T1_bl_tis) - 1.0)
                                * (exp(tau / T1_bl) * exp(tau / T1_bl_tis) - 1.0))
                              / (T1_bl + T1_bl_tis)
                          + (M_0_bl * T1_bl * T1_bl_tis * exp(-tes(k) / T2_bl)
                                * exp(-tes(k) / T1_bl_tis) * exp(-tis(j) / T1_bl)
                                * exp(-tis(j) / T1_bl_tis) * exp(t_tra / T1_bl_tis)
                                * (exp(tau / T1_bl) * exp(tau / T1_bl_tis) - 1.0))
                              / (T1_bl + T1_bl_tis))
                    * 2.0);
            }
        }
    }

    // Generating the output
    
    // Resize the ColumnVector signal to #tis * #tes
    ColumnVector signal;
    signal.ReSize(tis.Nrows() * tes.Nrows());

    // Fill in the results of the matrix delta_M_final
    for (int j = 1; j < tis.Nrows() + 1; j++)
    {
        /*if (isnan(signal))
        {
            kctissue = 0;
            double ti = tis(j);
            LOG << "Warning NaN in tissue curve at TI:" << ti << " with f:" << ftiss
                << " delt:" << delttiss << " tau:" << tauset << " T1:" << T_1 << " T1b:" << T_1b
                << endl;
        }*/

        for (int k = 1; k < tes.Nrows() + 1; k++)
        {
            signal((j - 1) * tes.Nrows() + k) = delta_M_final(j, k);
        }
    }

    // Resize the ColumnVector result to #tis * #tes * #repeats
    result.ReSize(tis.Nrows() * tes.Nrows() * repeats); 
    result = signal;

    // Concatenate signal repeats
    for (int j = 2; j < repeats + 1; j++)
    {
        result &= signal;
    }
}
