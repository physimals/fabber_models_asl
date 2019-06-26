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
    // number of repeats in data
    m_repeats = rundata.GetIntDefault("repeats", 1);
    
    // phases
    std::vector<double> phases_deg = rundata.GetDoubleList("ph", 0, 360);
    m_nphases = phases_deg.size();
    if (m_nphases > 0) 
    {
        m_phases_deg.ReSize(m_nphases);
        for (int i=0; i<m_nphases; i++) 
        {
            m_phases_deg(i+1) = phases_deg[i];
        }
    }
    else
    {
        // Phases have not been specified explicitly - read the number and space them
        // out evenly between 0 and 360
        m_nphases = rundata.GetIntDefault("nph", 8);
        m_phases_deg.ReSize(m_nphases);
        for (int i=0; i<m_nphases; i++)
        {
            m_phases_deg(i+1) = i * 360 / m_nphases;
        }
    }

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

        LOG << "Inference using numerical modulation function" << endl;
        LOG << "File is: " << modmatstring << endl;
    }
    else if (m_modfn == "fermi")
    {
        // shape of the fermi function
        m_alpha = rundata.GetDoubleDefault("alpha", 55);
        m_beta = rundata.GetDoubleDefault("beta", 12);
        
        LOG << "Inference using Fermi model" << endl;
        LOG << "alpha=" << m_alpha << " ,beta=" << m_beta << endl;
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
    params.push_back(Parameter(p++, "mag", DistParams(0, 1e12), DistParams(0, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    params.push_back(Parameter(p++, "phase", DistParams(0, 10.0/M_PI), DistParams(0, 10.0/M_PI), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    params.push_back(Parameter(p++, "offset", DistParams(0, 1e12), DistParams(0, 1e12), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
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

    // Take mean over the repeats and initialize the offset
    // to the average of the maximum and minimum intensity
    ColumnVector dmean(8);
    dmean = 0.0;
    for (int i = 1; i <= 8; i++)
    {
        for (int j = 1; j <= m_repeats; j++)
        {
            dmean(i) = dmean(i) + data((j - 1) * m_nphases + i);
        }
    }
    dmean = dmean / m_repeats;
    double dmax = dmean.Maximum();
    double dmin = dmean.Minimum();

    posterior.means(1) = (dmax - dmin) / 2;
    posterior.means(3) = (dmax + dmin) / 2;

    // Initialize the phase value from the point of max intensity. A peak
    // at 180 means a phase of 0.
    int ind;
    dmean.Maximum1(ind);
    float phase = m_phases_deg(ind) - 180;
    posterior.means(2) = phase * M_PI / 180;
}

void MultiPhaseASLFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Extract inferred parameters
    double mag = params(1);
    double phaserad = params(2);           // in radians
    double phase = params(2) * 180 / M_PI; // in degrees
    double offset = params(3);

    double flowvel = 0.3;
    if (m_incvel)
    {
        flowvel = params(4);
    }

    // Loop over phases to create result including repeated measurements
    int nn = m_nphases * m_repeats;
    result.ReSize(nn);
    for (int i = 1; i <= m_nphases; i++)
    {
        double evalfunc = 0;

        // Extract the measurement phase from list and convert to radians
        double ph_deg = m_phases_deg(i);
        if (ph_deg > 179) ph_deg -= 360;
        double ph_rad = ph_deg * M_PI / 180;

        if (m_modfn == "fermi")
        {
            // Use the Fermi modulation function
            // Note using the given values requires phases here to be in degrees
            evalfunc = mag * (-2 / (1 + exp((abs(ph_deg - phase) - m_alpha) / m_beta)))
                + offset;
        }
        else if (m_modfn == "mat")
        {
            // Evaluation modulation function from interpolation of values
            evalfunc = mag * (mod_fn(ph_rad - phaserad, flowvel)) + offset;
        }

        // Write the same output for each of the repeats 
        for (int j = 1; j <= m_repeats; j++)
        {
            result((j - 1) * m_nphases + i) = evalfunc;
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
