/*  fwdmodel_biexp.cc - Implements a model for correcting off resonance effect for multiphase pcASL

 Michael Chappell, QuBIc (IBME) & FMRIB Image Analysis Group

 Copyright (C) 2013 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_multiphase.h"

#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"

FactoryRegistration<FwdModelFactory, MultiPhaseASLFwdModel> MultiPhaseASLFwdModel::registration(
    "asl_multiphase");

static OptionSpec OPTIONS[] = {
    { "repeats", OPT_INT, "Number of repeats in data", OPT_NONREQ, "1" },
    { "modfn", OPT_STR, "Modulation function", OPT_NONREQ, "fermi" },
    { "modmat", OPT_MATRIX,
        "Modulation function matrix file, used if modfn=mat",
        OPT_NONREQ, "" },
    { "alpha", OPT_FLOAT, "Shape of the modulation function - alpha",
        OPT_NONREQ, "66" },
    { "beta", OPT_FLOAT, "Shape of the modulation function - beta",
        OPT_NONREQ, "12" },
    { "incvel", OPT_BOOL, "Include vel parameter", OPT_NONREQ, "" },
    { "infervel", OPT_BOOL, "Infer value of vel parameter",
        OPT_NONREQ, "" },
    { "" },
};

void MultiPhaseASLFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string MultiPhaseASLFwdModel::GetDescription() const
{
    return "ASL multiphase model";
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

void MultiPhaseASLFwdModel::HardcodedInitialDists(MVNDist &prior,
    MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e12;

    // Set priors

    // magnitude
    prior.means(1) = 0;
    precisions(1, 1) = 1e-12;

    // phase (radians)
    prior.means(2) = 0;
    precisions(2, 2) = M_PI / 10;

    // offset
    prior.means(3) = 0;
    precisions(3, 3) = 1e-12;

    //flow vel
    if (incvel)
    {
        prior.means(4) = 0.3;
        if (infervel)
        {
            precisions(4, 4) = 10;
        }
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;
}

void MultiPhaseASLFwdModel::InitParams(MVNDist &posterior) const
{
    // init the magntidue and offset parameters

    // mean over the repeats
    ColumnVector dmean(8);
    dmean = 0.0;
    for (int i = 1; i <= 8; i++)
    {
        for (int j = 1; j <= repeats; j++)
        {
            dmean(i) = dmean(i) + data((j - 1) * 8 + i);
        }
    }
    dmean = dmean / repeats;

    double dmax = dmean.Maximum();
    double dmin = dmean.Minimum();

    posterior.means(1) = (dmax - dmin) / 2;
    posterior.means(3) = (dmax + dmin) / 2;

    //init the mid phase value - by finding the point where the max intensity is
    int ind;
    float val;
    val = dmean.Maximum1(ind);    // find the max
    val = (ind - 1) * 180 / M_PI; //frequency of the minimum in ppm
    if (val > 179)
        val -= 360;
    val *= M_PI / 180;
    posterior.means(2) = val;

    if (incvel)
    {
        posterior.means(4) = 0.3;
    }
}

void MultiPhaseASLFwdModel::Evaluate(const ColumnVector &params,
    ColumnVector &result) const
{
    // ensure that values are reasonable
    // negative check

  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }

  // parameters that are inferred - extract and give sensible names
   double mag;
   double phase;
   double phaserad;
   double offset;
   double flowvel;

   mag=params(1);
   phaserad = params(2);  //in radians
   phase=params(2)*180/M_PI; //in degrees
   offset=params(3);

   if (incvel) {
     flowvel = params(4);
   }
   else {
     flowvel = 0.3;
   }   

   int nn=nph*repeats;
   result.ReSize(nn);
   // loop to create result
   for (int i=1; i<=nph; i++)
     {
       double evalfunc;

       
       double ph;
       ph = ph_list(i); //extract the measurement phase from list (vector)
       if (ph>179) ph -= 360; 
       double ph_rad = ph*M_PI/180; //in radians

       if (modfn == "fermi") {
	 // use the Fermi modulation function
	 evalfunc = mag*( -2/(1+exp( (abs(ph - phase) -alpha)/beta)) ) + offset; //note using the given values requires phases here to be in degrees
       }
       else if (modfn == "mat") {
	 //evaluation modulation function from interpolation of values
	 evalfunc = mag*(mod_fn(ph_rad - phaserad,flowvel)) + offset;
       }


       for (int j=1; j<=repeats; j++) {
	 result((j-1)*nph+i) = evalfunc;
       }
 
      }
    //cout << result.t();

  return;

}

FwdModel *MultiPhaseASLFwdModel::NewInstance()
{
    return new MultiPhaseASLFwdModel();
}

void MultiPhaseASLFwdModel::Initialize(ArgsType &args)
{

      // specify command line parameters here
      repeats = convertTo<int>(args.ReadWithDefault("repeats","1")); // number of repeats in data

      // phases
      string ph_temp;
      ph_temp = args.ReadWithDefault("ph1","none");
	if (ph_temp != "none" ) {
	  // a list of phases
	  ph_list.ReSize(1); //will add extra values onto end as needed
	  ph_list(1) = convertTo<double>(ph_temp);

	  while (true) //get the rest of the phases
	    {
	      int N = ph_list.Nrows()+1;
	      ph_temp = args.ReadWithDefault("ph"+stringify(N), "stop!");
	      if (ph_temp == "stop!") break; //we have run out of phases
	 
	      // append the new phase onto the end of the list
	      ColumnVector tmp(1);
	      tmp = convertTo<double>(ph_temp);
	      ph_list &= tmp; //vertical concatenation
	    }
	}
	else {
	  // phases have not been specified on command line  - use defaults
	  nph = convertTo<int>(args.ReadWithDefault("nph","8")); //number of phases
	  for (int i=1; i<=nph; i++)
	    {
	      //evenly spaced phases
	      double ph = (360/nph*(i-1) ); //in degrees
	      ColumnVector tmp(1);
	      tmp = ph;
	      ph_list &= tmp; //vertical concatention
	    }
	}

      // modulation function
      modfn = args.ReadWithDefault("modfn","fermi");      

      //modmat
      string modmatstring;
      Matrix mod_temp;
      modmatstring = args.ReadWithDefault("modmat","none");
      if (modmatstring != "none") {
	mod_temp = read_ascii_matrix(modmatstring);
      }

      // shape of the fermi function
      alpha = convertTo<double>(args.ReadWithDefault("alpha","55"));
      beta = convertTo<double>(args.ReadWithDefault("beta","12"));
      

      // deal with ARD selection
      //doard=false;
      //if (inferart==true && ardoff==false) { doard=true; }

      infervel=false;
      incvel=false;
      if ( modfn == "mat" ) {
	assert(mod_temp(1,1)==99);
	int nphasepts = mod_temp.Nrows()-1;
	nvelpts = mod_temp.Ncols()-1;
	
	mod_phase = ( mod_temp.SubMatrix(2,nphasepts+1,1,1) ).AsColumn();
	mod_v = ( mod_temp.SubMatrix(1,1,2,nvelpts+1) ).AsColumn();
	mod_mat = mod_temp.SubMatrix(2,nphasepts+1,2,nvelpts+1);
	
	vmax = mod_v(nvelpts);
	vmin = mod_v(1);
	
	
	infervel = args.ReadBool("infervel");
	if (infervel) {
	  incvel=true;
	}
	else {
	  incvel = args.ReadBool("incvel");
	}
	
      }

     
      // add information about the parameters to the log
      // test correctness of specified modulation function
      if (modfn == "fermi") {
	LOG << "Inference using Fermi model" << endl;
	LOG << "alpha=" << alpha << " ,beta=" << beta << endl;
      }
      else if (modfn == "mat") {
	LOG << "Inference using numerical modulation function" << endl;
	LOG << "File is: " << modmatstring << endl;
      }
      else {
	throw invalid_argument("Unrecognised modulation function");
      }
      
}

void MultiPhaseASLFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("mag");
    names.push_back("phase");
    names.push_back("offset");
    if (incvel)
    {
        names.push_back("vel");
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
    //deal with phase outside range -pi --> +pi
    //phase = asin(sin(phase)); //this assumes symmtery of function
    //if (phase>0) phase=std::fmod(phase+M_PI,2*M_PI)-M_PI;
    //else if (phase<0) phase=std::fmod(phase-M_PI,2*M_PI)+M_PI;
    // ** end old

    //bilinear interpolation
    if (v >= vmax)
    {
        ColumnVector usecolumn = mod_mat.Column(nvelpts);
        ans = interp(mod_phase, usecolumn, phase);
    }
    else if (v <= vmin)
    {
        ColumnVector usecolumn = mod_mat.Column(1);
        ans = interp(mod_phase, usecolumn, phase);
    }
    else
    {
        int ind = 1;
        while (v >= mod_v(ind))
            ind++;

        ColumnVector usecolumn = mod_mat.Column(ind - 1);
        double mod_l = interp(mod_phase, usecolumn, phase);
        ColumnVector usecolumn2 = mod_mat.Column(ind);
        double mod_u = interp(mod_phase, usecolumn2, phase);
        ans = mod_l
            + (v - mod_v(ind - 1)) / (mod_v(ind) - mod_v(ind - 1))
                * (mod_u - mod_l);
    }

    return ans;
}

double MultiPhaseASLFwdModel::interp(const ColumnVector &x,
    const ColumnVector &y, const double xi) const
// Look-up function for data table defined by x, y
// Returns the values yi at xi using linear interpolation
// Assumes that x is sorted in ascending order
// ? could be replaced my MISCMATHS:interp1 ?
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
