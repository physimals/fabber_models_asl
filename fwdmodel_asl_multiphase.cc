/**
 * fwdmodel_asl_multiphase.cc
 * 
 * Implements a model for correcting off resonance effect
 * for multiphase pcASL
 *
 * Michael Chappell, QuBIc (IBME) & FMRIB Image Analysis Group
 *
 * Copyright (C) 2013 University of Oxford  
 */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_multiphase.h"

#include <fabber_core/tools.h>
#include <fabber_core/fwdmodel.h>
#include <fabber_core/priors.h>

#include <newmat.h>

#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;

string MultiPhaseASLFwdModel::GetDescription() const 
{ 
    return "ASL multiphase model";
}

static OptionSpec OPTIONS[] = {
    { "repeats", OPT_INT, "Number of repeats in data", OPT_NONREQ, "1" },
    { "modfn", OPT_STR, "Modulation function", OPT_NONREQ, "fermi" },
    { "modmat", OPT_MATRIX, "Modulation function matrix file, used if modfn=mat", OPT_NONREQ, "" },
    { "alpha", OPT_FLOAT, "Shape of the modulation function - alpha", OPT_NONREQ, "66" },
    { "beta", OPT_FLOAT, "Shape of the modulation function - beta", OPT_NONREQ, "12" },
    { "incvel", OPT_BOOL, "Include vel parameter", OPT_NONREQ, "" },
    { "infervel", OPT_BOOL, "Infer value of vel parameter", OPT_NONREQ, "" },
    { "nph", OPT_INT, "Number of evenly-spaced phases between 0 and 360", OPT_NONREQ, "8" },
    { "ph<n>", OPT_FLOAT, "Individually-specified phase angles in degrees", OPT_NONREQ, "" },
    { "ntis", OPT_INT, "Number of delay times (TIs/PLDs)", OPT_NONREQ, "1" },
    { "multi-phase-offsets", OPT_BOOL, "Each TI/PLD has its own independent phase offset parameter", OPT_NONREQ, "" },
    { "" },
};

void MultiPhaseASLFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string MultiPhaseASLFwdModel::ModelVersion() const
{
    string version = "fwdmodel_asl_multiphase.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void MultiPhaseASLFwdModel::Initialize(ArgsType &rundata)
{
    // Number of TIs/PLDs
    m_ntis = rundata.GetIntDefault("ntis", 1);
    LOG << "MultiPhaseASLFwdModel::Data has " << m_ntis << " TIs/PLDs" << endl;

    // Number of repeats in data
    m_repeats = rundata.GetIntDefault("repeats", 1);
    LOG << "MultiPhaseASLFwdModel::Data has " << m_repeats << " repeats" << endl;

    m_multi_phase_offsets = (m_ntis > 1) && rundata.GetBool("multi-phase-offsets");
    if (m_multi_phase_offsets) LOG << "MultiPhaseASLFwdModel::Each TI/PLD has its own phase offset" << endl;
    else LOG << "MultiPhaseASLFwdModel::Using single phase offset for all TIs/PLDs" << endl;

    // The number of parameters inferred for each TI
    m_num_ti_params = m_multi_phase_offsets ? 3 : 2;
    
    // Phases - either specified explicitly or evenly spaced between 0 and 360
    std::vector<double> phases_deg = rundata.GetDoubleList("ph", 0, 360);
    m_nphases = phases_deg.size();
    if (m_nphases > 0) 
    {
        m_phases_deg.ReSize(m_nphases);
        LOG << "MultiPhaseASLFwdModel::Using " << m_nphases << " user-specified phases: ";
        for (int i=0; i<m_nphases; i++) 
        {
            m_phases_deg(i+1) = phases_deg[i];
            LOG << " " << phases_deg[i];
        }
    }
    else
    {
        // Phases have not been specified explicitly - read the number and space them
        // out evenly between 0 and 360
        m_nphases = rundata.GetIntDefault("nph", 8);
        m_phases_deg.ReSize(m_nphases);
        LOG << "MultiPhaseASLFwdModel::Using " << m_nphases << " evenly spaced phases: ";
        for (int i=0; i<m_nphases; i++)
        {
            m_phases_deg(i+1) = i * 360 / m_nphases;
            LOG << " " << m_phases_deg(i+1);
        }
    }
    LOG << endl;

    m_infervel = false;
    m_incvel = false;

    // modulation function
    m_modfn = rundata.GetStringDefault("modfn", "fermi");
    if (m_modfn == "mat")
    {
        // modmat
        string modmatstring;
        Matrix mod_temp;
        modmatstring = rundata.GetString("modmat");
        mod_temp = fabber::read_matrix_file(modmatstring);
        int nphasepts = mod_temp.Nrows() - 1;
        m_mod_nvelpts = mod_temp.Ncols() - 1;
        m_mod_phase = (mod_temp.SubMatrix(2, nphasepts + 1, 1, 1)).AsColumn();
        m_mod_v = (mod_temp.SubMatrix(1, 1, 2, m_mod_nvelpts + 1)).AsColumn();
        m_mod_mat = mod_temp.SubMatrix(2, nphasepts + 1, 2, m_mod_nvelpts + 1);

        m_mod_vmax = m_mod_v(m_mod_nvelpts);
        m_mod_vmin = m_mod_v(1);

        m_infervel = rundata.ReadBool("infervel");
        m_incvel = m_infervel || rundata.ReadBool("incvel");

        LOG << "MultiPhaseASLFwdModel::Inference using numerical modulation function" << endl;
        LOG << "MultiPhaseASLFwdModel::File is: " << modmatstring << endl;
    }
    else if (m_modfn == "fermi")
    {
        // shape of the fermi function
        m_alpha = rundata.GetDoubleDefault("alpha", 55);
        m_beta = rundata.GetDoubleDefault("beta", 12);
        
        LOG << "MultiPhaseASLFwdModel::Inference using Fermi model" << endl;
        LOG << "MultiPhaseASLFwdModel::alpha=" << m_alpha << " ,beta=" << m_beta << endl;
    }
    else
    {
        throw InvalidOptionValue("modfn", m_modfn, "Must be fermi or mat");
    }
}

void MultiPhaseASLFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    if (m_ntis == 1) 
    {
        // For compatibility don't number outputs if there is only one TI/PLD
        params.push_back(Parameter(p++, "mag", DistParams(0, 1e12), DistParams(0, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
        params.push_back(Parameter(p++, "offset", DistParams(0, 1e12), DistParams(0, 1e12), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    }
    else 
    {
        // Separate magnitude/offset for each TI/PLD 
        for (int i=0; i<m_ntis; i++) 
        {
            params.push_back(Parameter(p++, "mag" + stringify(i+1), DistParams(0, 1e12), DistParams(0, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
            params.push_back(Parameter(p++, "offset" + stringify(i+1), DistParams(0, 1e12), DistParams(0, 1e12), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
            if (m_multi_phase_offsets) 
            {
                params.push_back(Parameter(p++, "phase" + stringify(i+1), DistParams(0, 100.0), DistParams(0, 10.0), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
            }
        }
    }

    // Parameters common to all TIs
    if (!m_multi_phase_offsets || (m_ntis == 1))
    {
        params.push_back(Parameter(p++, "phase", DistParams(0, 100.0), DistParams(0, 10.0), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    }

    if (m_incvel)
    {
        // If we are not inferring the velocity set the variance to be very low
        double var = m_infervel ? 0.1 : 1e-12;
        params.push_back(Parameter(p++, "vel", DistParams(0.3, var), DistParams(0.3, var), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    }
}

void MultiPhaseASLFwdModel::InitVoxelPosterior(MVNDist &posterior) const
{
    // Initialize the magntidue and offset parameters
    if (data.Nrows() != m_nphases*m_ntis*m_repeats)
    {
        throw InvalidOptionValue("Num phases * num TIs * num repeats", 
                                 stringify(m_nphases*m_ntis*m_repeats), 
                                 "Should match number of data volumes: " + stringify(data.Nrows()));
    }

    // For each TI/PLD, take mean over the repeats and initialize the offset
    // to the average of the maximum and minimum intensity
    for (int t=0; t<m_ntis; t++) 
    {
        ColumnVector dmean(m_nphases);
        dmean = 0.0;
        for (int r=0; r<m_repeats; r++)
        {
            for (int p=1; p<=m_nphases; p++)
            {
                dmean(p) = dmean(p) + data(t*m_repeats*m_nphases + r*m_nphases + p);
            }
        }
        dmean = dmean / m_repeats;
        double dmax = dmean.Maximum();
        double dmin = dmean.Minimum();

        posterior.means(m_num_ti_params*t + 1) = (dmax - dmin) / 2;
        posterior.means(m_num_ti_params*t + 2) = (dmax + dmin) / 2;
    }

    // Initialize the phase value from the point of max intensity (averaged over
    // all repeats/TIs). 
    // With a phase offset of 0, the maximum will be at 180 and the mininum 
    // at 0. Phase offset of 90 gives max/min at 270/90, etc. In general the minimum
    // occurs at the phase offset value.
    // Note that the phase offset is in radians but the the measured phases are input 
    // in degrees, and we initialize the phase offset in the range -pi to pi.
    ColumnVector dmean(m_nphases);
    dmean = 0.0;
    for (int i = 1; i <= m_nphases; i++)
    {
        for (int j = 0; j < m_repeats*m_ntis; j++)
        {
            dmean(i) = dmean(i) + data(j * m_nphases + i);
        }
    }

    int ind;
    dmean.Minimum1(ind);
    float phase_init = m_phases_deg(ind);
    while (phase_init > 180) phase -= 360; 
    if (m_multi_phase_offsets)
    {
        for (int t=0; t<m_ntis; t++) 
        {
            posterior.means(m_num_ti_params*t + 3) = phase_init * M_PI / 180;
        }
    }
    else 
    {
        posterior.means(m_num_ti_params*m_ntis + 1) = phase_init * M_PI / 180;
    }
}

void MultiPhaseASLFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Extract inferred parameters
    ColumnVector mag(m_ntis);
    ColumnVector offset(m_ntis);
    ColumnVector phase_off_rad(m_ntis);
    for (int t=0; t<m_ntis; t++) 
    {
        mag(t+1) = params(m_num_ti_params*t+1);
        offset(t+1) = params(m_num_ti_params*t+2);
        if (m_multi_phase_offsets)
        {
            phase_off_rad(t+1) = params(m_num_ti_params*t+3);
        }
        else
        {
            phase_off_rad(t+1) = params(m_num_ti_params*m_ntis+1);
        }
    }

    ColumnVector phase_off_deg = phase_off_rad * 180 / M_PI;    // in degrees
    double flowvel = 0.3;
    if (m_incvel)
    {
        flowvel = params(m_num_ti_params*m_ntis + 2);
    }

    // Loop over TIs and phases to create result including repeated measurements
    int nn = m_nphases * m_repeats * m_ntis;
    result.ReSize(nn);
    for (int t=0; t<m_ntis; t++) 
    {
        for (int p=1; p<=m_nphases; p++)
        {
            double evalfunc = 0;

            // Extract the measurement phase from list, apply offset and convert
            // to radians. Force range to -180 -> 180 / -pi -> pi
            double phase_meas_deg = m_phases_deg(p);
            double phase_actual_deg = phase_meas_deg - phase_off_deg(t+1);
            while (phase_actual_deg < -180) phase_actual_deg += 360;
            while (phase_actual_deg > 180) phase_actual_deg -= 360;
            double phase_actual_rad = phase_actual_deg * M_PI / 180;
            if (m_modfn == "fermi")
            {
                // Use the Fermi modulation function
                // Note using the given values requires phases here to be in degrees
                evalfunc = mag(t+1) * (-2 / (1 + exp((abs(phase_actual_deg) - m_alpha) / m_beta)))
                    + offset(t+1);
            }
            else if (m_modfn == "mat")
            {
                // Evaluation modulation function from interpolation of values
                // FIXME range of phase values: 0 to 2pi or -pi to pi???
                evalfunc = mag(t+1) * (mod_fn(phase_actual_rad, flowvel)) + offset(t+1);
            }

            // Write the same output for each of the repeats 
            for (int r=0; r<m_repeats; r++)
            {
                result(t*m_repeats*m_nphases + r*m_nphases + p) = evalfunc;
            }
        }
    }
}

double MultiPhaseASLFwdModel::mod_fn(const double inphase, const double v) const
{
    // the modulation function - evaluated from interpolation of modmat
    double ans;
    double phase = inphase;

    // phase will be normally in range 0 --> 2*pi
    if (phase < 0.0)
        phase = 0.0;
    else if (phase > 2 * M_PI)
        phase = 2 * M_PI;
    // old from veasl model
    // deal with phase outside range -pi --> +pi
    // phase = asin(sin(phase)); //this assumes symmtery of function
    // if (phase>0) phase=std::fmod(phase+M_PI,2*M_PI)-M_PI;
    // else if (phase<0) phase=std::fmod(phase-M_PI,2*M_PI)+M_PI;
    // ** end old

    // bilinear interpolation
    if (v >= m_mod_vmax)
    {
        ColumnVector usecolumn = m_mod_mat.Column(m_mod_nvelpts);
        ans = interp(m_mod_phase, usecolumn, phase);
    }
    else if (v <= m_mod_vmin)
    {
        ColumnVector usecolumn = m_mod_mat.Column(1);
        ans = interp(m_mod_phase, usecolumn, phase);
    }
    else
    {
        int ind = 1;
        while (v >= m_mod_v(ind))
            ind++;

        ColumnVector usecolumn = m_mod_mat.Column(ind - 1);
        double mod_l = interp(m_mod_phase, usecolumn, phase);
        ColumnVector usecolumn2 = m_mod_mat.Column(ind);
        double mod_u = interp(m_mod_phase, usecolumn2, phase);
        ans = mod_l + (v - m_mod_v(ind - 1)) / (m_mod_v(ind) - m_mod_v(ind - 1)) * (mod_u - mod_l);
    }

    return ans;
}

// Look-up function for data table defined by x, y
// Returns the values yi at xi using linear interpolation
// Assumes that x is sorted in ascending order
// ? could be replaced my MISCMATHS:interp1 ?
double MultiPhaseASLFwdModel::interp(
    const ColumnVector &x, const ColumnVector &y, const double xi) const
{
    double ans;
    if (xi >= x.Maximum())
        ans = y(x.Nrows());
    else if (xi <= x.Minimum())
        ans = y(1);
    else
    {
        int ind = 1;
        while (xi >= x(ind))
            ind++;
        double xa = x(ind - 1), xb = x(ind), ya = y(ind - 1), yb = y(ind);
        ans = ya + (xi - xa) / (xb - xa) * (yb - ya);
    }
    return ans;
}

FwdModel *MultiPhaseASLFwdModel::NewInstance() { return new MultiPhaseASLFwdModel(); }

FactoryRegistration<FwdModelFactory, MultiPhaseASLFwdModel> MultiPhaseASLFwdModel::registration(
    "asl_multiphase");
