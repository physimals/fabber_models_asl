/**
 * fwdmodel_asl_2compartment.cc
 *
 * Martin Craig <martin.craig@eng.ox.ac.uk>
 *
 * Copyright (C) 2008 University of Oxford
 */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_2compartment.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>
#include <fabber_core/tools.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

// Variance for non-informative priors
#define VAR_NONINFORM 1e8

FactoryRegistration<FwdModelFactory, ASL2CompartmentModel> ASL2CompartmentModel::registration(
    "asl_2comp");

FwdModel *ASL2CompartmentModel::NewInstance() { return new ASL2CompartmentModel(); }
static OptionSpec OPTIONS[] = {
    { "bolusdur", OPT_FLOAT, "Bolus duration", OPT_REQ, "" },
    { "solution", OPT_STR, "Solution type - dist, fast or slow", OPT_NONREQ, "dist" },
    { "ie", OPT_FLOAT, "Inversion efficiency", OPT_NONREQ, "0.7" },
    { "ps", OPT_FLOAT, "Permeability surface area", OPT_NONREQ, "0.8" },
    { "vb", OPT_FLOAT, "Blood volume fraction", OPT_NONREQ, "0.03" },
    { "vbw", OPT_FLOAT, "", OPT_NONREQ, "0.7" },
    { "f", OPT_FLOAT, "Initial CBF value", OPT_NONREQ, "0.1" },
    { "sig0", OPT_FLOAT, "Initial signal scale factor", OPT_NONREQ, "1" },
    { "infer-bat", OPT_BOOL, "Infer bolus arrival time", OPT_NONREQ, "" },
    { "bat", OPT_FLOAT, "Initial bolus arrival time (fixed if not inferred)", OPT_NONREQ, "0.7" },
    { "infer-t1", OPT_BOOL, "Infer tissue T1", OPT_NONREQ, "" },
    { "t1", OPT_FLOAT, "Initial tissue T1 (fixed if not inferred)", OPT_NONREQ, "1.3" },
    { "infer-t1b", OPT_BOOL, "Infer blood T1", OPT_NONREQ, "" },
    { "t1b", OPT_FLOAT, "Initial blood T1 (fixed if not inferred)", OPT_NONREQ, "1.65" },
    { "infer-t1e", OPT_BOOL, "Infer extracellular T1", OPT_NONREQ, "" },
    { "t1e", OPT_FLOAT, "Initial extracellular T1 (fixed if not inferred)", OPT_NONREQ, "0.6" },
    { "" },
};

void ASL2CompartmentModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string ASL2CompartmentModel::GetDescription() const
{
    return "ASL 2-compartment model from Parkes/Tofts 2002";
}

string ASL2CompartmentModel::ModelVersion() const
{
    string version = "fwdmodel_asl_2compartment.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void ASL2CompartmentModel::Initialize(FabberRunData &rundata)
{
    string sol = rundata.GetStringDefault("solution", "dist");
    if (sol == "dist")
    {
        m_sol = DIST;
    }
    else if (sol == "fast")
    {
        m_sol = FAST;
    }
    else if (sol == "slow")
    {
        m_sol = SLOW;
    }
    else
    {
        throw InvalidOptionValue("solution", sol, "Must be 'dist', 'fast' or 'slow'");
    }

    // Fixed parameters
    m_ie = rundata.GetDoubleDefault("ie", 1.0, 0, 1);
    m_bolusdur = rundata.GetDouble("bolusdur", 0);
    m_ps = rundata.GetDoubleDefault("ps", 0.8);
    m_vb = rundata.GetDoubleDefault("vb", 0.03);
    m_vbw = rundata.GetDoubleDefault("vbw", 0.7);

    // Timing parameters (TIs / PLDs)
    vector<double> tis = rundata.GetDoubleList("ti");
    bool ti_set = (tis.size() >= 1);

    vector<double> plds = rundata.GetDoubleList("pld");
    bool pld_set = (plds.size() >= 1);

    if (ti_set && pld_set)
    {
        throw FabberRunDataError("Cannot specify TIs and PLDs at the same time");
    }
    if (pld_set)
    {
        for (int i = 0; i < plds.size(); i++)
        {
            m_times.push_back(plds[i] + m_bolusdur);
        }
    }
    else if (ti_set)
    {
        m_times = tis;
    }
    else
    {
        throw FabberRunDataError("Either TIs or PLDs must be specified");
    }

    // Initial values of inferred parameters
    m_f_initial = rundata.GetDoubleDefault("f", 0.1, 0);
    m_sig0 = rundata.GetDoubleDefault("sig0", 1, 0);

    // Optionally inferred parameters
    m_bat = rundata.GetDoubleDefault("bat", 0.7, 0);
    m_infer_bat = rundata.GetBool("infer-bat");

    m_t1 = rundata.GetDoubleDefault("t1", 1.3, 0);
    m_infer_t1 = rundata.GetBool("infer-t1");

    m_t1b = rundata.GetDoubleDefault("t1b", 1.65, 0);
    m_infer_t1b = rundata.GetBool("infer-t1b");

    m_t1e = rundata.GetDoubleDefault("t1e", 1.3, 0);
    m_infer_t1e = rundata.GetBool("infer-t1e");
}

void ASL2CompartmentModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p = 0;
    params.push_back(Parameter(p++, "f", DistParams(m_f_initial, VAR_NONINFORM),
        DistParams(m_f_initial, 10), PRIOR_NORMAL, TRANSFORM_LOG()));
    // params.push_back(Parameter(p++, "sig0", DistParams(m_sig0_initial, 1e12),
    // DistParams(m_sig0_initial, 10)));
    if (m_infer_bat)
    {
        params.push_back(Parameter(p++, "bat", DistParams(m_bat, VAR_NONINFORM),
            DistParams(m_bat, 10), PRIOR_NORMAL, TRANSFORM_LOG()));
    }
    if (m_infer_t1)
    {
        params.push_back(Parameter(p++, "t1", DistParams(m_t1, VAR_NONINFORM), DistParams(m_t1, 10),
            PRIOR_NORMAL, TRANSFORM_LOG()));
    }
    if (m_infer_t1b)
    {
        params.push_back(Parameter(p++, "t1b", DistParams(m_t1b, VAR_NONINFORM),
            DistParams(m_t1b, 10), PRIOR_NORMAL, TRANSFORM_LOG()));
    }
    if (m_infer_t1e)
    {
        params.push_back(Parameter(p++, "t1e", DistParams(m_t1e, VAR_NONINFORM),
            DistParams(m_t1e, 10), PRIOR_NORMAL, TRANSFORM_LOG()));
    }
}

/**
 * Get change in magnetization
 *
 * Equations [22] and [23] in Parkes and Tofts 2002
 *
 * dM(t)	change in magnetization 							-> signal
 * f		blood perfusion										-> parameter
 * bat		arrival time    		 							-> parameter (optional)
 * t1		T1 in tissue									    -> parameter (optional)
 * t1b		T1 in blood											-> parameter (optional)
 * t1e		T1 in extravascular space							-> parameter (optional)
 * m0a 		equib magnetization of arterial blood				-> parameter (absorbed into sig0?)
 * m_ie   	inversion efficiency							    -> fixed
 * m_bolusdur	labelling time									-> fixed
 * m_ps		Permeability / surface area product				    -> fixed
 * m_vbw	Blood water volume                                  -> fixed
 * m_vb		Blood volume                                        -> fixed
 * tp       t - bat                                             -> derived
 * v_bw     m_vb * m_vbw                                        -> derived
 * A		m_ps / v_bw										    -> derived
 * C		1/t1e											    -> derived
 * D		1/t1b											    -> derived
 * J		A+D for t < MTT									    -> derived
 *   		A+D+f/m_vb for t >= MTT							    -> derived
 */
double ASL2CompartmentModel::dM(
    double t, double f, double bat, double t1, double t1b, double t1e) const
{
    // Whatever... Just a scale factor
    double m0a = 1;

    // Derived variables
    double tp = t - bat;
    double v_bw = m_vb * m_vbw;
    double A = m_ps / v_bw;
    double C = 1 / t1e;
    double D = 1 / t1b;

    // J depends on solution method
    double J;
    double MTT;
    // Flow in correct units to compare with m_vb
    double fs = f / 6000;
    switch (m_sol)
    {
    case FAST:
        J = A + D + fs / m_vb;
        break;
    case SLOW:
        J = A + D;
        break;
    case DIST:
    default:
        MTT = m_vb / fs;
        if (t < MTT)
            J = A + D;
        else
            J = A + D + fs / m_vb;
        break;
    }

    if (tp < 0)
    {
        return 0;
    }
    else if (tp <= m_bolusdur)
    {
        // Equation [22]
        double term1 = (1 - exp(-J * tp)) / J;
        double term2 = A * ((J - C + C * exp(-J * tp) - J * exp(-C * tp)) / (J * C * (J - C)));
        double factor = 2 * f * m0a * m_ie * exp(-D * bat);
        return factor * (term1 + term2);
    }
    else
    {
        // Equation [23]
        double term1 = (1 / J + A / (J * (C - J))) * ((exp(J * m_bolusdur) - 1) * exp(-J * tp));
        double term2 = -A * exp(-C * tp) * (exp(C * m_bolusdur) - 1) / (C * (C - J));
        double factor = 2 * f * m0a * m_ie * exp(-D * bat);
        return factor * (term1 + term2);
    }
}

double ASL2CompartmentModel::sig(double sig0, double dm) const
{
    // double lambda = 0.9; // (vew + vbw)/m_vbw
    return sig0 * dm;
}

void ASL2CompartmentModel::EvaluateModel(
    const ColumnVector &params, ColumnVector &result, const std::string &key) const
{
    int p = 1;
    int nt = data.Nrows();
    if (nt != m_times.size())
    {
        throw FabberRunDataError("The number of data volumes was " + stringify(nt)
            + " but the number of TIs/PLDs specified was " + stringify(m_times.size()));
    }

    // Perfusion parameter
    double f = params(p++);

    // Equilibrium signal
    // double sig0 = params(p++);

    double bat = m_bat;
    if (m_infer_bat)
    {
        bat = params(p++);
    }

    double t1 = m_t1;
    if (m_infer_t1)
    {
        t1 = params(p++);
    }

    double t1b = m_t1b;
    if (m_infer_t1b)
    {
        t1b = params(p++);
    }

    double t1e = m_t1e;
    if (m_infer_t1e)
    {
        t1e = params(p++);
    }

    result.ReSize(nt);
    for (int i = 1; i <= nt; i++)
    {
        result(i) = sig(m_sig0, dM(m_times[i - 1], f, bat, t1, t1b, t1e));
    }
}
