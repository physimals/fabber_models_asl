/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Michael Chappell, FMRIB Image Analysis & IBME QuBIc Groups

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include <fabber_core/fabber_core.h>

#include "fwdmodel_asl_2compartment.h"
#include "fwdmodel_asl_grase.h"
#include "fwdmodel_asl_multiphase.h"
#include "fwdmodel_asl_quasar.h"
#include "fwdmodel_asl_rest.h"
#include "fwdmodel_asl_satrecov.h"
#include "fwdmodel_asl_satrecovdualfa.h"
#include "fwdmodel_asl_turboquasar.h"
#include "fwdmodel_asl_multite.h"
#include "fwdmodel_asl_velocityselective.h"

int main(int argc, char **argv)
{
    // This should not be necessary! They are supposed to auto-register
    // but this is not working in current build against FSL 6.0.6
    // Possibly side-effect of move to shared libraries vs static
    ASLFwdModel::NewInstance();
    GraseFwdModel::NewInstance();
    MultiPhaseASLFwdModel::NewInstance();
    QuasarFwdModel::NewInstance();
    TurboQuasarFwdModel::NewInstance();
    ASL2CompartmentModel::NewInstance();
    SatrecovFwdModel::NewInstance();
    multiTEFwdModel::NewInstance();
    VelocitySelectiveFwdModel::NewInstance();

    return execute(argc, argv);
}
