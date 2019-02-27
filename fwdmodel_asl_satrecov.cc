/*  fwdmodel_asl_satrecov.cc - Saturation Recovery curve calibration for ASL

 Michael Chappell, IBME & FMRIB Image Analysis Group

 Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_satrecov.h"

#include "fabber_core/fwdmodel.h"

#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

FactoryRegistration<FwdModelFactory, SatrecovFwdModel> SatrecovFwdModel::registration("satrecov");

FwdModel *SatrecovFwdModel::NewInstance() { return new SatrecovFwdModel(); }

static OptionSpec OPTIONS[] = {
    { "repeats", OPT_INT, "Number of repeats in data", OPT_NONREQ, "1" },
    { "t1", OPT_FLOAT, "T1 value (s)", OPT_NONREQ, "1.3" },
    { "phases", OPT_INT, "Number of phases", OPT_NONREQ, "1" },
    { "slicedt", OPT_FLOAT, "Increase in TI per slice", OPT_NONREQ, "0.0" },
    { "fixa", OPT_BOOL, "Fix the A parameter where it will be ambiguous", OPT_NONREQ, "" },
    { "FA", OPT_FLOAT, "Flip angle in degrees for Look-Locker readout", OPT_NONREQ, "0" },
    { "m_lfa", OPT_FLOAT, "Low flip angle in degrees for Look-Locker readout", OPT_NONREQ, "0" },
    { "ti<n>", OPT_FLOAT, "List of TI values", OPT_NONREQ, "" }, { "" },
};

void SatrecovFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string SatrecovFwdModel::ModelVersion() const
{
    string version = "fwdmodel_asl_satrecov.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

string SatrecovFwdModel::GetDescription() const 
{ 
    return "Saturation recovery ASL model"; 
}

void SatrecovFwdModel::Initialize(ArgsType &args)
{
    // Basic acquisition parameters
    m_tis = args.GetDoubleList("ti", 0);
    m_repeats = args.GetIntDefault("repeats", 1);
    m_t1 = args.GetDoubleDefault("t1", 1.3);
    m_nphases = args.GetIntDefault("phases", 1);
    m_slicedt = args.GetDoubleDefault("slicedt", 0);

    // Assuming even sampling - this only applies to LL acquisitions
    if (m_tis.size() < 2) throw InvalidOptionValue("Number of TIs", stringify(m_tis.size()), "Need at least 2 TIs");
    m_dti = m_tis[1] - m_tis[0]; 

    // To fix the A parameter where it will be ambiguous
    m_fix_a = args.GetBool("fixa");

    // With a look locker readout
    double fa_deg = args.GetDoubleDefault("FA", 0);
    m_look_locker = (fa_deg > 0.1);
    m_fa = fa_deg * M_PI / 180; // convert to radians
        
    double lfa_deg = args.GetDoubleDefault("LFA", 0);
    m_lfa_on = (m_lfa > 0);
    m_lfa = lfa_deg * M_PI / 180; // convert to radians
    m_dg = 0.023;
    if (m_look_locker)
    {
        LOG << "SatrecovFwdModel::Look-Locker mode" << endl;
        LOG << "SatrecovFwdModel::FA (degrees)=" << fa_deg << endl;
        if (m_lfa_on) LOG << "SatrecovFwdModel::LFA (degrees)=" << lfa_deg << endl;
    }
}

void SatrecovFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    params.push_back(Parameter(p++, "M0t", DistParams(0, 1e12), DistParams(10, 10)));
    params.push_back(Parameter(p++, "T1t", DistParams(m_t1, 0.1), DistParams(m_t1, 0.1)));
    double a_var = 0.1;
    if (m_fix_a) a_var = 1e-12;
    params.push_back(Parameter(p++, "A", DistParams(1, a_var), DistParams(1, a_var)));
    if (m_lfa_on) {
        params.push_back(Parameter(p++, "g", DistParams(1, 0.01), DistParams(1, 0.01)));
    }
}

void SatrecovFwdModel::EvaluateModel(const ColumnVector &params, ColumnVector &result, const std::string &key) const
{
    // Ensure that values are reasonable - negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= NumParams(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }

    float M0t = paramcpy(1);
    float T1t = paramcpy(2);
    float A = paramcpy(3);

    float g = 1.0;
    if (m_lfa_on) g = paramcpy(4);

    // Duo flip angle correction
    // Details refer to Moss's Turbo QUASAR methodology paper or Esben's thesis p122
    float FA = g * m_fa;
    float lFA = (g + dg) * m_lfa; // dg only corrects for the low flip angle.

    float T1tp = T1t;
    float M0tp = M0t;

    if (m_look_locker)
    {
        // Note that we do not have sin(FA) here - we actually estiamte the M0
        // at the flip angle used for the readout!
        T1tp = 1 / (1 / T1t - log(cos(FA)) / m_dti); // FA is in radians
        M0tp = M0t * (1 - exp(-m_dti / T1t)) / (1 - cos(FA) * exp(-m_dti / T1t));
    }

    if (m_lfa_on)
        result.ReSize(m_tis.size() * (m_nphases + 1) * m_repeats);
    else
        result.ReSize(m_tis.size() * m_nphases * m_repeats);

    int nti = m_tis.size();
    if (nti != data.Nrows()) throw InvalidOptionValue("ti<n>", stringify(nti), "Number of TIs does not match number of volumes in data");

    for (int ph = 1; ph <= m_nphases; ph++)
    {
        for (int it = 0; it < nti; it++)
        {
            for (int rpt = 1; rpt <= m_repeats; rpt++)
            {
                double ti = m_tis[it] + m_slicedt * coord_z;
                result((ph - 1) * (nti * m_repeats) + it * m_repeats + rpt)
                    = M0tp * (1 - A * exp(-ti / T1tp));
            }
        }
    }

    if (m_lfa_on)
    {
        int ph = m_nphases + 1;
        T1tp = 1 / (1 / T1t - log(cos(lFA)) / m_dti);
        M0tp = M0t * (1 - exp(-m_dti / T1t)) / (1 - cos(lFA) * exp(-m_dti / T1t));
        for (int it = 0; it < nti; it++)
        {
            for (int rpt = 1; rpt <= m_repeats; rpt++)
            {
                // slicedt accounts for increase in delay between slices
                double ti = m_tis[it] + m_slicedt * coord_z;
                // Note the sin(m_lfa)/sin(FA) term since the M0 we estimate is
                // actually MOt*sin(FA)
                result((ph - 1) * (nti * m_repeats) + it * m_repeats + rpt)
                    = M0tp * sin(lFA) / sin(FA) * (1 - A * exp(-m_tis[it] / T1tp));
            }
        }
    }
}
