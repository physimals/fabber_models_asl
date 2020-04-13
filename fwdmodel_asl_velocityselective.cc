/**
 * fwdmodel_asl_velocityselective.cc
 *
 * Moss Zhao <mosszhao@stanford.edu>
 *
 * Copyright (C) 2020 Stanford University
 */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_velocityselective.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>
#include <fabber_core/inference.h>
#include <fabber_core/tools.h>
#include <newmatio.h>
#include <iostream>
#include <stdexcept>
#include <miscmaths/miscprob.h>

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;

FactoryRegistration<FwdModelFactory, VelocitySelectiveFwdModel> VelocitySelectiveFwdModel::registration(
    "velocity_selective");

FwdModel *VelocitySelectiveFwdModel::NewInstance() { return new VelocitySelectiveFwdModel(); }

string VelocitySelectiveFwdModel::ModelVersion() const
{
    string version = "fwdmodel_asl_velocityselective.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

static OptionSpec OPTIONS[] = {
    { "tau", OPT_FLOAT, "Bolus duration in seconds. GE DATA DICOM field <Private_0043_10a5>", OPT_REQ, "1.59" },
    { "repeats", OPT_INT, "Number of repeats in data.", OPT_REQ, "1" },
    { "t1b", OPT_FLOAT, "T1 blood value in seconds.", OPT_NONREQ, "1.65" },
    { "slicedt", OPT_FLOAT, "Increase in TI per slice.", OPT_NONREQ, "0" },
    { "t2b", OPT_FLOAT, "T2 blood value in seconds.", OPT_NONREQ, "0.116" }, //J. Zhao et al 2007 MRM, with Hct=0.44 and SatO2=0.99, more suitable for the labeling pulse
    //In GE's raw header file:
    //effective_te=header.RawHeader.echotimes(1) - header.RawHeader.te/1e6
    //effective_te=0.0360-0.010704
    { "te", OPT_FLOAT, "Effective TE in seconds. Note:In GE's raw header file, effective_te=header.RawHeader.echotimes(1) - header.RawHeader.te/1e6", OPT_NONREQ, "0.0253" },
    { "ti<n>", OPT_FLOAT, "List of TI values in seconds. Note: this should be the time between labeling and imaging. DICOM field <InversionTime> of GE's VS ASL data.", OPT_NONREQ, "" },
    { "predelay", OPT_FLOAT, "Predelay time in seconds. This should be a user defined variable during scanning.", OPT_NONREQ, "1.9" },
    { "tissoff", OPT_BOOL, "Turn off tissue CPT. Not advised to use", OPT_NONREQ, "" },
    { "" },
};

void VelocitySelectiveFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string VelocitySelectiveFwdModel::GetDescription() const
{
    return "Velocity Selective ASL model";
}



void VelocitySelectiveFwdModel::Initialize(ArgsType &args)
{   
    repeats_initial = convertTo<int>(args.ReadWithDefault("repeats", "1")); // number of repeats in data
    t1b_initial = convertTo<double>(args.ReadWithDefault("t1b", "1.65"));
    predelay_initial = convertTo<double>(args.ReadWithDefault("predelay", "1.9"));
    tau_initial = convertTo<double>(args.ReadWithDefault("tau", "1.59"));
    slicedt_initial = convertTo<double>(args.ReadWithDefault("slicedt", "0"));
    te_initial = convertTo<double>(args.ReadWithDefault("te", "0.0253"));
    t2b_initial = convertTo<double>(args.ReadWithDefault("t2b", "0.116"));

    // special - turn off tissue cpt
    infertiss = true;
    bool tissoff = args.ReadBool("tissoff");
    if (tissoff)
        infertiss = false;

    // Deal with tis
    tis.ReSize(1); // will add extra values onto end as needed
    tis(1) = atof(args.Read("ti1").c_str());
    // get the rest of the tis
    while (true)
    {
        int N = tis.Nrows() + 1;
        string tiString = args.ReadWithDefault("ti" + stringify(N), "stop!");
        if (tiString == "stop!")
            break; // we have run out of tis

        // append the new ti onto the end of the list
        ColumnVector tmp(1);
        tmp = convertTo<double>(tiString);
        tis &= tmp; // vertical concatenation
    }

    // add information about the parameters to the log
    LOG << "Inference using development model" << endl;
    LOG << "    Data parameters: #repeats = " << repeats_initial
        //<< ", t1 = " << t1
        << ", t1b = " << t1b_initial
        << ", predelay = " << predelay_initial
        << ", t2b = " << t2b_initial
        << ", te = " << te_initial;
    LOG << ", bolus length (tau) = " << tau_initial << endl;
    if (infertiss)
    {
        LOG << "Infertting on tissue component " << endl;
    }
    //if (infert1)
    //{
        LOG << "Infering on T1 values " << endl;
    //}
    LOG << "TIs: ";
    for (int i = 1; i <= tis.Nrows(); i++)
        LOG << tis(i) << " ";
    LOG << endl;

}


void VelocitySelectiveFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    // In VS ASL, we don't need to estimate arrival time
    if (infertiss) {
        names.push_back("ftiss");
    }
}


void VelocitySelectiveFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
    if (infertiss)
    {   
        // Index start from zero. We only need CBF parameter for now
        prior.means(tiss_index()) = 0;
        precisions(tiss_index(), tiss_index()) = 1e-12;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital
    // posterior
    // Tissue perfusion
    if (infertiss)
    {
        posterior.means(tiss_index()) = 10;
        precisions(tiss_index(), tiss_index()) = 1;
    }

    posterior.SetPrecisions(precisions);
}


void VelocitySelectiveFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{

    // ensure that values are reasonable
    // negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= NumParams(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }

    // parameters to be inferred - extract and give sensible names
    double ftiss;
    if (infertiss) {
        ftiss = paramcpy(tiss_index());
    }



    // sequence parameter TIs - we copy them from the user proivded values
    ColumnVector the_tis;
    the_tis = tis;
    the_tis += slicedt_initial * coord_z; // account here for an increase in delay between slices

    // sequence parameter tau - we copy them from the user proivded values
    double the_tau = tau_initial;
    // sequence parameter T_1b - we copy them from the user proivded values
    double the_T_1b = t1b_initial;
    // sequence parameter T_1b - we copy them from the user proivded values
    double the_predelay = predelay_initial; 
     // sequence parameter T_1b - we copy them from the user proivded values
    double the_te = te_initial; 
    // sequence parameter T_1b - we copy them from the user proivded values
    double the_T_2b = t2b_initial;    


    /*****************************/
    /* Begin model fitting here  */
    /*****************************/

    // ASL difference signal. Tissue kinetic model
    ColumnVector kctissue(tis.Nrows());
    kctissue = 0.0;

    if (infertiss) {
        kctissue = kctissue_model(ftiss, the_tis, the_tau, the_T_1b, the_predelay, the_te, the_T_2b);
    }

    // Nan catching
    bool cont = true;
    int it = 1;
    while (cont)
    {
        if (isnan(kctissue(it)) | isinf(kctissue(it)))
        {
            LOG << "Warning NaN in kctissue" << endl;
            LOG << "params: " << params.t() << endl;
            LOG << "kctissue: " << kctissue.t() << endl;
            cont = false;
            kctissue = 0.0;
        }
        it++;
        if (it > kctissue.Nrows())
            cont = false;
    }

    result = kctissue;

    /*
    // assemble the result
    int nti = tis.Nrows();
    result.ReSize(tis.Nrows() * repeats);

    for (int it = 1; it <= tis.Nrows(); it++)
    {

        
        // loop over the repeats
        for (int rpt = 1; rpt <= repeats; rpt++)
        {
            int tiref = (it - 1) * repeats + rpt;

            result(tiref) = ftiss * kctissue(it);
        }
        
    }
    */
}

// --- Kinetic curve functions ---
// Reference: Eric Wong, 2006. https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.20906
ColumnVector VelocitySelectiveFwdModel::kctissue_model(double ftiss, const ColumnVector &tis, double tau, double T_1b, double predelay, double te, double T_2b) const {

    ColumnVector kctissue(tis.Nrows());
    kctissue = 0.0;

    // Perfusion quantification equation of Wong's 2006 paper 
    double alpha = exp((-1) * te / T_2b);
    double M_ZB = alpha * (1 - exp((-1) * predelay / T_1b));

    // VS ASL uses saturation to in the labeling pulse, thus it is a 90 degree inversion so we don't need to 'times 2' in the model
    // The predelay component is to account for the time that between the readout and next inversion
    // Ref: Figure 1 in Perfusion imaging using FAIR with a short predelay, Jinyuan Zhou, 1999
    for (int it = 1; it <= tis.Nrows(); it++) {

        kctissue(it) = M_ZB * ftiss * tau * exp((-1) * tis(it)/ T_1b) * exp((-1) * te / T_2b);
    }

    return kctissue;
}






