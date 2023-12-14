#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void ML_BSSN_bench_ParamCompat(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(my_initial_boundary_condition, "extrapolate-gammas")) {
    // Override initial_boundary_condition
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing "
                               "ML_BSSN_bench::initial_boundary_condition="
                               "\"extrapolate-gammas\" because "
                               "ML_BSSN_bench::my_initial_boundary_condition="
                               "\"extrapolate-gammas\"");
    int ierr = CCTK_ParameterSet("initial_boundary_condition", "ML_BSSN_bench",
                                 "extrapolate-gammas");
    if (ierr) {
      CCTK_ERROR("Could not set "
                 "ML_BSSN_bench::initial_boundary_condition=\"extrapolate-gammas\"");
    }
  }

  if (CCTK_EQUALS(my_rhs_boundary_condition, "NewRad")) {
    // Override rhs_boundary_condition
    CCTK_WARN(CCTK_WARN_ALERT,
              "Forcing ML_BSSN_bench::rhs_boundary_condition=\"NewRad\" because "
              "ML_BSSN_bench::my_rhs_boundary_condition=\"NewRad\"");
    int ierr = CCTK_ParameterSet("rhs_boundary_condition", "ML_BSSN_bench", "NewRad");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::rhs_boundary_condition=\"NewRad\"");
    }
  } else if (CCTK_EQUALS(my_rhs_boundary_condition, "static")) {
    // Override rhs_boundary_condition
    CCTK_WARN(CCTK_WARN_ALERT,
              "Forcing ML_BSSN_bench::rhs_boundary_condition=\"scalar\" because "
              "ML_BSSN_bench::my_rhs_boundary_condition=\"static\" (sic!)");
    int ierr = CCTK_ParameterSet("rhs_boundary_condition", "ML_BSSN_bench", "scalar");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::rhs_boundary_condition=\"scalar\"");
    }
  }

  if (CCTK_EQUALS(apply_dissipation, "never")) {
    // Override epsDiss to disable dissipation
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::epsDiss=0.0 because "
                               "ML_BSSN_bench::apply_dissipation=\"never\"");
    int ierr = CCTK_ParameterSet("epsDiss", "ML_BSSN_bench", "0.0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::epsDiss=0.0");
    }
  }

  if (LapseACoeff == 0.0) {
    // Override evolveA
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::evolveA=0 because "
                               "ML_BSSN_bench::LapseACoeff=0.0");
    int ierr = CCTK_ParameterSet("evolveA", "ML_BSSN_bench", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::evolveA=0");
    }
  } else if (LapseACoeff == 1.0) {
    // Override evolveA
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::evolveA=1 because "
                               "ML_BSSN_bench::LapseACoeff=1.0");
    int ierr = CCTK_ParameterSet("evolveA", "ML_BSSN_bench", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::evolveA=1");
    }
  }

  // Make a modifiable copy
  CCTK_INT my_evolveB = evolveB;
  if (ShiftBCoeff == 0.0) {
    // Override evolveB
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::evolveB=0 because "
                               "ML_BSSN_bench::ShiftBCoeff=0.0");
    int ierr = CCTK_ParameterSet("evolveB", "ML_BSSN_bench", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::evolveB=0");
    }
    my_evolveB = 0;
  } else if (ShiftBCoeff == 1.0) {
    // Override evolveB
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::evolveB=1 because "
                               "ML_BSSN_bench::ShiftBCoeff=1.0");
    int ierr = CCTK_ParameterSet("evolveB", "ML_BSSN_bench", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::evolveB=1");
    }
    my_evolveB = 1;
  }

  if (my_evolveB && shiftGammaCoeff == 0.0) {
    // Override evolveB
    CCTK_WARN(
        CCTK_WARN_ALERT,
        "Forcing ML_BSSN_bench::evolveB=0 because ML_BSSN_bench::shiftGammaCoeff=0.0");
    int ierr = CCTK_ParameterSet("evolveB", "ML_BSSN_bench", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::evolveB=0");
    }
  }

  if (LapseAdvectionCoeff == 0.0) {
    // Override advectLapse
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::advectLapse=0 because "
                               "ML_BSSN_bench::LapseAdvectionCoeff=0.0");
    int ierr = CCTK_ParameterSet("advectLapse", "ML_BSSN_bench", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::advectLapse=0");
    }
  } else if (LapseAdvectionCoeff == 1.0) {
    // Override advectLapse
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::advectLapse=1 because "
                               "ML_BSSN_bench::LapseAdvectionCoeff=1.0");
    int ierr = CCTK_ParameterSet("advectLapse", "ML_BSSN_bench", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::advectLapse=1");
    }
  }

  if (ShiftAdvectionCoeff == 0.0) {
    // Override advectShift
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::advectShift=0 because "
                               "ML_BSSN_bench::ShiftAdvectionCoeff=0.0");
    int ierr = CCTK_ParameterSet("advectShift", "ML_BSSN_bench", "0");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::advectShift=0");
    }
  } else if (ShiftAdvectionCoeff == 1.0) {
    // Override advectShift
    CCTK_WARN(CCTK_WARN_ALERT, "Forcing ML_BSSN_bench::advectShift=1 because "
                               "ML_BSSN_bench::ShiftAdvectionCoeff=1.0");
    int ierr = CCTK_ParameterSet("advectShift", "ML_BSSN_bench", "1");
    if (ierr) {
      CCTK_ERROR("Could not set ML_BSSN_bench::advectShift=1");
    }
  }
}

void ML_BSSN_bench_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (!CCTK_EQUALS(my_initial_data, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::my_initial_data is "
                               "outdated; please update the parameter file. Do "
                               "not use this parameter, and set up initial "
                               "conditions via ADMBase as usual.");
  }

  if (!CCTK_EQUALS(my_initial_boundary_condition, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT,
              "Parameter ML_BSSN_bench::my_initial_boundary_condition is outdated; "
              "please update the parameter file. Do not use this parameter, "
              "and set up initial boundary conditions as usual.");
  }

  if (!CCTK_EQUALS(my_rhs_boundary_condition, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::my_rhs_boundary_condition "
                               "is outdated; please update the parameter file. "
                               "Do not use this parameter, and set up RHS "
                               "boundary conditions as usual.");
  }

  if (!CCTK_EQUALS(my_boundary_condition, "default")) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::my_boundary_condition is "
                               "outdated; please update the parameter file. Do "
                               "not use this parameter, and set up RHS "
                               "boundary conditions as usual.");
  }

  if (evolveB && shiftGammaCoeff == 0.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::shiftGammaCoeff should not "
                               "be set to 0.0 any more to disable evolving "
                               "B^i. Instead, set ML_BSSN_bench::evolveB=0.");
  }

  if (LapseACoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::LapseACoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_BSSN_bench::evolveA.");
  }

  if (ShiftBCoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::ShiftBCoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_BSSN_bench::evolveB.");
  }

  if (LapseAdvectionCoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::LapseAdvectionCoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_BSSN_bench::advectLapse.");
  }
  if (ShiftAdvectionCoeff != -1.0) {
    CCTK_WARN(CCTK_WARN_ALERT, "Parameter ML_BSSN_bench::ShiftAdvectionCoeff is "
                               "outdated; please update the parameter file. "
                               "Instead of using this parameter, you should "
                               "set ML_BSSN_bench::advectShift.");
  }
}
