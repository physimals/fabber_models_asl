/*  fwdmodel_asl_grase.cc - Implements the dual echo GRASE resting state ASL
 model

 Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_grase.h"

#include "fabber_core/inference.h"
#include "fabber_core/tools.h"

#include "miscmaths/miscprob.h"
#include "armawrap/newmat.h"
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;

FactoryRegistration<FwdModelFactory, GraseFwdModel> GraseFwdModel::registration("buxton");

static OptionSpec OPTIONS[] = {
    { "t1", OPT_FLOAT, "T1 value", OPT_NONREQ, "1.3" },
    { "t1b", OPT_FLOAT, "T1b value", OPT_NONREQ, "1.65" },
    { "lambda", OPT_FLOAT, "lambda value", OPT_NONREQ, "0.9" },
    { "ti<n>", OPT_FLOAT, "List of TI values. At least one required", OPT_NONREQ, "" },
    { "repeats", OPT_INT, "Number of repeats in data", OPT_NONREQ, "1" },
    { "pretisat", OPT_FLOAT, "Deal with saturation of the bolus a fixed time pre TI measurement",
        OPT_NONREQ, "0.0" },
    { "tau", OPT_FLOAT, "Bolus duration. Default is effectively infinite", OPT_NONREQ, "1000" },
    { "casl", OPT_BOOL, "Data is CASL (not PASL)", OPT_NONREQ, "PASL" },
    { "infertau", OPT_BOOL, "Infer bolus duration parameter", OPT_NONREQ, "" },
    { "infert1", OPT_BOOL, "Infer T1 parameter", OPT_NONREQ, "" },
    { "inferart", OPT_BOOL, "Infer arterial parameters", OPT_NONREQ, "" },
    { "bat", OPT_FLOAT, "Bolus arrival time", OPT_NONREQ, "0.7" },
    { "batsd", OPT_FLOAT, "Bolus arrival time standard deviation", OPT_NONREQ, "0.316" },
    { "slicedt", OPT_FLOAT, "Increase in TI per slice", OPT_NONREQ, "0.0" },
    { "calib", OPT_BOOL, "Data has already been subjected to calibration", OPT_NONREQ, "" },
    { "ardoff", OPT_BOOL, "Disable ARD", OPT_NONREQ, "" },
    { "tauboff", OPT_BOOL, "Forces the inference of arterial bolus off", OPT_NONREQ, "" }, { "" },
};

void GraseFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void GraseFwdModel::Initialize(ArgsType &args)
{
    repeats = convertTo<int>(args.ReadWithDefault("repeats", "1")); // number of repeats in data
    t1 = convertTo<double>(args.ReadWithDefault("t1", "1.3"));
    t1b = convertTo<double>(args.ReadWithDefault("t1b", "1.65"));
    lambda = convertTo<double>(args.ReadWithDefault("lambda", "0.9"));

    pretisat
        = convertTo<double>(args.ReadWithDefault("pretisat", "0")); // deal with saturation of the
                                                                    // bolus a fixed time pre TI
                                                                    // measurement
    grase = args.ReadBool("grase"); // DEPRECEATED data has come from the
                                    // GRASE-ASL sequence - therefore apply
                                    // pretisat of 0.1s
    if (grase)
        pretisat = 0.1;

    casl = args.ReadBool("casl"); // set if the data is CASL or PASL (default)
    slicedt = convertTo<double>(args.ReadWithDefault("slicedt", "0.0")); // increase in TI per slice

    calib = args.ReadBool("calib");
    infertau = args.ReadBool("infertau"); // infer on bolus length?
    infert1 = args.ReadBool("infert1");   // infer on T1 values?
    inferart = args.ReadBool("inferart"); // infer on arterial compartment?
    // inferinveff = args.ReadBool("inferinveff"); //infer on a linear decrease
    // in inversion efficiency?
    // infertrailing = args.ReadBool("infertrailing"); //infers a trailing edge
    // bolus slope using new model
    seqtau
        = convertTo<double>(args.ReadWithDefault("tau", "1000")); // bolus length as set by sequence
                                                                  // (default of 1000 is effectively
                                                                  // infinite
    setdelt = convertTo<double>(args.ReadWithDefault("bat", "0.7"));
    double deltsd; // std dev for delt prior
    deltsd = convertTo<double>(args.ReadWithDefault("batsd", "0.316"));
    deltprec = 1 / (deltsd * deltsd);

    bool ardoff = false;
    ardoff = args.ReadBool("ardoff");
    bool tauboff = false;
    tauboff = args.ReadBool("tauboff"); // forces the inference of arterial bolus off

    // combination options
    infertaub = false;
    if (inferart && infertau && !tauboff)
        infertaub = true;

    // deal with ARD selection
    doard = false;
    if (inferart == true && ardoff == false)
    {
        doard = true;
    }

    /* if (infertrailing) {
     if (!infertau) {
     // do not permit trailing edge inference without inferring on bolus length
     throw Invalid_option("--infertrailing has been set without setting
     --infertau");
     }
     else if (inferinveff)
     //do not permit trailing edge inference and inversion efficiency inference
     (they are mututally exclusive)
     throw Invalid_option("--infertrailing and --inferinveff may not both be
     set");
     }*/

    // Deal with tis
    tis.ReSize(1); // will add extra values onto end as needed
    tis(1) = atof(args.Read("ti1", "0").c_str());

    while (true) // get the rest of the tis
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
    timax = tis.Maximum(); // dtermine the final TI

    // need to set the voxel coordinates to a deafult of 0 (for the times we
    // call the model before we start handling data)
    coord_x = 0;
    coord_y = 0;
    coord_z = 0;

    singleti = false; // normally we do multi TI ASL
                      /* This option is currently disabled since it is not compatible with basil
                       if (tis.Nrows()==1) {
                       //only one TI therefore only infer on CBF and ignore other inference
                       options
                       LOG << "--Single inversion time mode--" << endl;
                       LOG << "Only a sinlge inversion time has been supplied," << endl;
                       LOG << "Therefore only tissue perfusion will be inferred." << endl;
                       LOG << "-----" << endl;
                       singleti = true;
                       // force other inference options to be false
                       //infertau = false; infert1 = false; inferart = false; //inferinveff =
                       false;
                       }
                       */

    // add information about the parameters to the log
    LOG << "Inference using Buxton Kinetic Curve model" << endl;
    if (!casl)
        LOG << "Data being analysed using PASL inversion profile" << endl;
    if (casl)
        LOG << "Data being analysed using CASL inversion profile" << endl;
    if (pretisat > 0)
        LOG << "Saturation of" << pretisat << "s before TI has been specified" << endl;
    if (grase)
        LOG << "Using pre TI saturation of 0.1 for GRASE-ASL sequence" << endl;
    if (calib)
        LOG << "Input data is in physioligcal units, using estimated CBF in "
               "T_1app calculation"
            << endl;
    LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
    LOG << ", bolus length (tau) = " << seqtau << endl;
    if (infertau)
    {
        LOG << "Infering on bolus length " << endl;
    }
    if (inferart)
    {
        LOG << "Infering on artertial compartment " << endl;
    }
    if (doard)
    {
        LOG << "ARD has been set on arterial compartment " << endl;
    }
    if (infert1)
    {
        LOG << "Infering on T1 values " << endl;
    }
    /*if (inferinveff) {
     LOG << "Infering on Inversion Efficency slope " << endl; }
     if (infertrailing) {
     LOG << "Infering bolus trailing edge period" << endl; }*/
    LOG << "TIs: ";
    for (int i = 1; i <= tis.Nrows(); i++)
        LOG << tis(i) << " ";
    LOG << endl;
}

std::string GraseFwdModel::GetDescription() const
{
    return "Implements the dual echo GRASE resting state ASL model";
}

string GraseFwdModel::ModelVersion() const
{
    string version = "fwdmodel_asl_grase.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void GraseFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
    prior.means(tiss_index()) = 0;
    precisions(tiss_index(), tiss_index()) = 1e-12;

    if (!singleti)
    {
        // Tissue bolus transit delay
        prior.means(tiss_index() + 1) = setdelt;
        precisions(tiss_index() + 1, tiss_index() + 1) = deltprec;
    }

    // Tissue bolus length
    if (infertau)
    {
        prior.means(tau_index()) = seqtau;
        precisions(tau_index(), tau_index()) = 10;
    }

    if (infertaub)
    {
        prior.means(taub_index()) = seqtau;
        precisions(taub_index(), taub_index()) = 10;
    }

    // Arterial Perfusion & bolus delay
    if (inferart)
    {
        int aidx = art_index();
        prior.means(aidx) = 0;
        prior.means(aidx + 1) = 0.5;
        precisions(aidx + 1, aidx + 1) = 10;
        precisions(aidx, aidx) = 1e-12;
    }

    // T1 & T1b
    if (infert1)
    {
        int tidx = t1_index();
        prior.means(tidx) = t1;
        prior.means(tidx + 1) = t1b;
        precisions(tidx, tidx) = 100;
        precisions(tidx + 1, tidx + 1) = 100;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital
    // posterior

    // Tissue perfusion
    posterior.means(tiss_index()) = 10;
    precisions(tiss_index(), tiss_index()) = 1;

    // Arterial perfusion
    if (inferart)
    {
        posterior.means(art_index()) = 10;
        precisions(art_index(), art_index()) = 1;
    }
    posterior.SetPrecisions(precisions);
}

void GraseFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
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

    // sensible limits on transit times
    if (!singleti)
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

    // parameters that are inferred - extract and give sensible names
    float ftiss;
    float delttiss;
    float tauset; // the value of tau set by the sequence (may be effectively
                  // infinite)
    float taubset;
    float fblood;
    float deltblood;
    float T_1;
    float T_1b;

    ftiss = paramcpy(tiss_index());
    if (!singleti)
    {
        delttiss = paramcpy(tiss_index() + 1);
    }
    else
    {
        // only inferring on tissue perfusion, assume fixed value for tissue
        // arrival time
        delttiss = setdelt;
    }

    if (infertau)
    {
        tauset = paramcpy(tau_index());
    }
    else
    {
        tauset = seqtau;
    }

    if (infertaub)
    {
        taubset = paramcpy(taub_index());
    }
    else
    {
        taubset = tauset;
    }

    if (inferart)
    {
        fblood = paramcpy(art_index());
        deltblood = paramcpy(art_index() + 1);
    }
    else
    {
        fblood = 0;
        deltblood = 0;
    }

    if (infert1)
    {
        T_1 = paramcpy(t1_index());
        T_1b = paramcpy(t1_index() + 1);

        // T1 cannot be zero!
        if (T_1 < 1e-12)
            T_1 = 0.01;
        if (T_1b < 1e-12)
            T_1b = 0.01;
    }
    else
    {
        T_1 = t1;
        T_1b = t1b;
    }

    // float lambda = 0.9;

    float f_calib;
    // if we are using calibrated data then we can use ftiss to calculate T_1app
    if (calib)
        f_calib = ftiss;
    else
        f_calib = 0.01; // otherwise assume sensible value (units of s^-1)

    float T_1app = 1 / (1 / T_1 + f_calib / lambda);
    float R = 1 / T_1app - 1 / T_1b;

    float tau;  // bolus length as seen by kintic curve
    float taub; // bolus length of blood as seen in signal

    float F = 0;
    float kctissue;
    float kcblood;

    // loop over tis
    float ti;
    result.ReSize(tis.Nrows() * repeats);

    for (int it = 1; it <= tis.Nrows(); it++)
    {
        ti = tis(it) + slicedt * coord_z; // account here for an increase in the
                                          // TI due to delays between slices
        // cout << coord_z << " : " << ti << endl;
        if (casl)
            F = 2 * ftiss;
        else
            F = 2 * ftiss * exp(-ti / T_1app);

        /* According to EAGLE GRASE sequence bolus length is current TI - 0.1s
         (assuming infite length 'true' bolus)
         However, also allow here bolus length to be finite as recorded in
         tauset
         NB tauset is the 'true' bolus length, tau is what the tissue actually
         sees as a result of the sequence
         25-3-2009 now deal with this scenario via pretisat parameter
         */

        /*	if (grase)
         {
         //GRASE -  deal with bolus length (see above) */

        // Deal with saturation of the bolus before the TI - defined by pretisat
        if (tauset < ti - pretisat)
        {
            tau = tauset;
        }
        else
        {
            tau = ti - pretisat;
        }

        if (taubset < ti - pretisat)
        {
            taub = taubset;
        }
        else
        {
            taub = ti - pretisat;
        }

        /* }
         else {
         tau = tauset;
         taub = taubset;
         }
         */

        /*	if (infertrailing)
         {
         ///////////////// New bolus length & trailing edge model //////////////
         // tau is now defined as the point where inv eff starts to drop
         // bt1 is start of bolus trailing edge, bt2 is end of bolus trailing
         edge

         // constrain trailing epriod to sensible limits
         if (trailingperiod < 1e-6) trailingperiod=1e-6;
         //if (trailingperiod > 2*tau) trailingperiod = 2*tau;


         // define bolus trailgin edge start and stop times
         float bt1 = delttiss + tau ;
         float bt2 = delttiss + tau + trailingperiod;

         // tissue contribution
         if(ti < delttiss)
         { kctissue = 0;}

         // once bolus is arriving, but before the trailing edge has arrived
         else if(ti >= delttiss && ti <= (delttiss + bt1))
         {
         kctissue = F/R * (exp(R*ti)-exp(R*delttiss) );
         }

         // once the trailing edge of the bolus is arriving
         else if(ti >= bt1 && ti <= bt2)
         {
         kctissue = F/R * ( (exp(R*bt1)-exp(R*delttiss)) +
         (exp(R*ti)-exp(R*bt1))*(1+inveffslope*bt1+inveffslope/R) -
         inveffslope*(ti*exp(R*ti) - bt1*exp(R*bt1)) );
         if (kctissue<0) {kctissue = 0; } //dont allow negative values (shoudld
         be redundant)
         }

         // once the trailing edge of the bolus is past
         else //(ti > delttiss + bt2)
         {
         kctissue = F/R * ( (exp(R*bt1)-exp(R*delttiss)) +
         (exp(R*bt2)-exp(R*bt1))*(1+inveffslope*bt1+inveffslope/R) -
         inveffslope*(bt2*exp(R*bt2) - bt1*exp(R*bt1)) );
         if (kctissue<0) { kctissue = 0; } //dont allow negative values
         }


         // arterial contribution
         // calc the correct bt1 for the arterial bolus
         bt1 = deltblood + tau ;
         //bt2 = deltblood + tau + trailingperiod;

         if(ti < deltblood)
         {
         //kcblood = 0;
         // use a artifical lead in period for arterial bolus to improve model
         fitting
         kcblood = fblood * exp(-deltblood/T_1b) * (0.98 * exp(
         (ti-deltblood)/0.1 ) + 0.02 * ti/deltblood );
         }

         // once bolus is arriving, but before the trailing edge has arrived
         else if(ti >= deltblood && ti <= bt1)
         {
         kcblood = fblood * exp(-ti/T_1b);
         }

         // Once we are into the trailing edge of the bolus
         else if(ti >= bt1 && ti <= bt2)
         {
         kcblood = fblood * exp(-ti/T_1b);
         kcblood = kcblood * (1 - 1/trailingperiod*(ti- bt1));
         if (kcblood<0) { kcblood = 0; } //dont allow negative values
         }

         // Once the trailing edge of the bolus has passed
         else //(ti >  bt2)
         {
         kcblood = 0; //end of bolus
         }
         }




         /////////// Older version of model (implements Inversion efficiency
         slope) /////////
         else
         {*/
        // deal with the case where the inveffslope is severe and cuts off the
        // bolus)
        // if (tau>bollen) tau=bollen;
        // if (taub>bollen) taub=bollen;
        // --[tissue contribution]------
        if (ti < delttiss)
        {
            kctissue = 0;
        }
        else if (ti >= delttiss && ti <= (delttiss + tau))
        {
            if (casl)
                kctissue
                    = F * T_1app * exp(-delttiss / T_1b) * (1 - exp(-(ti - delttiss) / T_1app));
            else
                kctissue = F / R * ((exp(R * ti) - exp(R * delttiss)));
        }
        else //(ti > delttiss + tau)
        {
            if (casl)
                kctissue = F * T_1app * exp(-delttiss / T_1b) * exp(-(ti - tau - delttiss) / T_1app)
                    * (1 - exp(-tau / T_1app));
            else
                kctissue = F / R * ((exp(R * (delttiss + tau)) - exp(R * delttiss)));
        }

        // --[arterial contribution]------
        if (ti < deltblood)
        {
            kcblood = fblood * exp(-deltblood / T_1b)
                * (0.98 * exp((ti - deltblood) / 0.05) + 0.02 * ti / deltblood);
            // use a artifical lead in period for arterial bolus to improve
            // model fitting
            // NB same equation for PASL and CASL
        }
        else if (ti >= deltblood && ti <= (deltblood + taub))
        {
            if (casl)
                kcblood = fblood * exp(-deltblood / T_1b);
            else
                kcblood = fblood * exp(-ti / T_1b);
        }
        else //(ti > deltblood + tau)
        {
            kcblood = 0; // end of bolus
            if (casl)
                kcblood = fblood * exp(-deltblood / T_1b);
            else
                kcblood = fblood * exp(-(deltblood + taub) / T_1b);
            kcblood *= (0.98 * exp(-(ti - deltblood - taub) / 0.05)
                + 0.02 * (1 - (ti - deltblood - taub) / 5));
            // artifical lead out period for taub model fitting
            if (kcblood < 0)
                kcblood = 0; // negative values are possible with the lead out
                             // period equation
        }

        if (isnan(kctissue))
        {
            kctissue = 0;
            LOG << "Warning NaN in tissue curve at TI:" << ti << " with f:" << ftiss
                << " delt:" << delttiss << " tau:" << tau << " T1:" << T_1 << " T1b:" << T_1b
                << endl;
        }
        //}

        /* output */
        // loop over the repeats
        for (int rpt = 1; rpt <= repeats; rpt++)
        {
            result((it - 1) * repeats + rpt) = kctissue + kcblood;
        }
    }

    return;
}

FwdModel *GraseFwdModel::NewInstance() { return new GraseFwdModel(); }
void GraseFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("ftiss");
    if (!singleti)
        names.push_back("delttiss");
    if (infertau)
    {
        names.push_back("tautiss");
    }
    if (inferart)
    {
        names.push_back("fblood");
        names.push_back("deltblood");
    }
    if (infert1)
    {
        names.push_back("T_1");
        names.push_back("T_1b");
    }
    /* if (inferinveff) {
     names.push_back("Inveffslope");
     }
     if (infertrailing) {
     names.push_back("trailingperiod");
     }*/
    if (infertaub)
    {
        names.push_back("taublood");
    }
}

void GraseFwdModel::SetupARD(const MVNDist &theta, MVNDist &thetaPrior, double &Fard) const
{
    int ardindex = ard_index();

    if (doard)
    {
        SymmetricMatrix PriorPrec;
        PriorPrec = thetaPrior.GetPrecisions();

        PriorPrec(ardindex, ardindex) = 1e-12; // set prior to be initally non-informative

        thetaPrior.SetPrecisions(PriorPrec);

        thetaPrior.means(ardindex) = 0;

        // set the Free energy contribution from ARD term
        SymmetricMatrix PostCov = theta.GetCovariance();
        double b
            = 2 / (theta.means(ardindex) * theta.means(ardindex) + PostCov(ardindex, ardindex));
        Fard = -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5)
            - 0.5 * log(b); // taking c as 0.5 - which it will be!
    }

    return;
}

void GraseFwdModel::UpdateARD(const MVNDist &theta, MVNDist &thetaPrior, double &Fard) const
{
    int ardindex = ard_index();

    if (doard)
    {
        SymmetricMatrix PriorCov;
        SymmetricMatrix PostCov;
        PriorCov = thetaPrior.GetCovariance();
        PostCov = theta.GetCovariance();

        PriorCov(ardindex, ardindex)
            = theta.means(ardindex) * theta.means(ardindex) + PostCov(ardindex, ardindex);

        thetaPrior.SetCovariance(PriorCov);

        // Calculate the extra terms for the free energy
        double b
            = 2 / (theta.means(ardindex) * theta.means(ardindex) + PostCov(ardindex, ardindex));
        Fard = -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5)
            - 0.5 * log(b); // taking c as 0.5 - which it will be!
    }

    return;
}
