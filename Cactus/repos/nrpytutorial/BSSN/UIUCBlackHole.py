# This module sets up UIUC Black Hole initial data in terms of
# the variables used in BSSN_RHSs.py

# Authors: Zachariah B. Etienne, zachetie **at** gmail **dot** com
#          Terrence Pierre Jacques, terrencepierrej **at** gmail **dot** com
#          Ian Ruchlin

# ## This module sets up initial data for a merging black hole system in spherical coordinates. We can convert from spherical to any coordinate system defined in [reference_metric.py](../edit/reference_metric.py) (e.g., SinhSpherical, Cylindrical, Cartesian, etc.) using the [Exact ADM Spherical-to-BSSNCurvilinear converter module](Tutorial-ADM_Initial_Data-Converting_Exact_ADM_Spherical_or_Cartesian_to_BSSNCurvilinear.ipynb)
#
# ### Here we set up UIUC Black Hole initial data ([Liu, Etienne, & Shapiro, PRD 80 121503, 2009](https://arxiv.org/abs/1001.4077)):
#
# UIUC black holes have the advantage of finite coordinate radius in the maximal spin limit. It is therefore excellent for studying very highly spinning black holes. This module sets the UIUC black hole at the origin.
#
# **Inputs for initial data**:
#
# * The black hole mass $M$.
# * The dimensionless spin parameter $\chi = a/M$
#
# **Additional variables needed for spacetime evolution**:
#
# * Desired coordinate system
# * Desired initial lapse $\alpha$ and shift $\beta^i$. We will choose our gauge conditions as $\alpha=1$ and $\beta^i=B^i=0$. $\alpha = \psi^{-2}$ will yield much better behavior, but the conformal factor $\psi$ depends on the desired *destination* coordinate system (which may not be spherical coordinates).

# Step P0: Load needed modules
import sympy as sp              # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par  # NRPy+: Parameter interface
import indexedexp as ixp        # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm  # NRPy+: Reference metric support
import sys                      # Standard Python module for multiplatform OS-level functions
import BSSN.BSSN_in_terms_of_ADM as BitoA  # Needed to compute conformal factor
import BSSN.BSSN_quantities as Bq  # Sets default for EvolvedConformalFactor_cf

thismodule = __name__

# The UIUC initial data represent a Kerr black hole with mass M
#  and dimensionless spin chi in UIUC quasi-isotropic coordinates,
#   see https://arxiv.org/abs/1001.4077
# Input parameters:
M, chi = par.Cparameters("REAL", thismodule, ["M", "chi"], [1.0, 0.99])


def UIUCBlackHole():
    global r,th,ph, gammaDD, KDD, alpha, betaU, BU

    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 1: Set psi, the conformal factor:
    # Spin per unit mass
    a = M*chi

    # Defined under equation 1 in Liu, Etienne, & Shapiro (2009) https://arxiv.org/pdf/1001.4077.pdf
    # Boyer - Lindquist outer horizon
    rp = M + sp.sqrt(M**2 - a**2)
    # Boyer - Lindquist inner horizon
    rm = M - sp.sqrt(M**2 - a**2)

    # Boyer - Lindquist radius in terms of UIUC radius
    # Eq. 11
    # r_{BL} = r * ( 1 + r_+ / 4r )^2
    rBL = r*(1 + rp / (4*r))**2

    # Expressions found below Eq. 2
    # Sigma = r_{BL}^2 + a^2 cos^2 theta
    SIG = rBL**2 + a**2*sp.cos(th)**2

    # Delta = r_{BL}^2 - 2Mr_{BL} + a^2
    DEL = rBL**2 - 2*M*rBL + a**2

    # A = (r_{BL}^2 + a^2)^2 - Delta a^2 sin^2 theta
    AA = (rBL**2 + a**2)**2 - DEL*a**2*sp.sin(th)**2

    # *** The ADM 3-metric in spherical basis ***
    gammaDD = ixp.zerorank2()
    # Declare the nonzero components of the 3-metric
    # (Eq. 13 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf):

    # ds^2 = Sigma (r + r_+/4)^2 / ( r^3 (r_{BL} - r_- ) * dr^2 +
                # Sigma d theta^2  +  (A sin^2 theta) / Sigma  *  d\phi^2

    gammaDD[0][0] = ((SIG*(r + rp/4)**2)/(r**3*(rBL - rm)))
    gammaDD[1][1] = SIG
    gammaDD[2][2] = AA/SIG*sp.sin(th)**2

    # *** The physical trace-free extrinsic curvature in spherical basis ***
    # Nonzero components of the extrinsic curvature K, given by
    # Eq. 14 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf:
    KDD     = ixp.zerorank2()


    # K_{r phi} = K_{phi r} = (Ma sin^2 theta) / (Sigma sqrt{A Sigma}) *
    #     [3r^4_{BL} + 2a^2 r^2_{BL} - a^4 - a^2 (r^2_{BL} - a^2) sin^2 theta] *
    #     (1 + r_+ / 4r) (1 / sqrt{r(r_{BL} - r_-)})

    KDD[0][2] = KDD[2][0] = (M*a*sp.sin(th)**2)/(SIG*sp.sqrt(AA*SIG))*\
                    (3*rBL**4 + 2*a**2*rBL**2 - a**4- a**2*(rBL**2 - a**2)*\
                     sp.sin(th)**2)*(1 + rp/(4*r))*1/sp.sqrt(r*(rBL - rm))

    # Components of the extrinsic curvature K, given by
    # Eq. 15 of Liu, Etienne, & Shapiro, https://arxiv.org/pdf/1001.4077.pdf:

    # K_{theta phi} = K_{phi theta} = -(2a^3 Mr_{BL} cos theta sin^3 theta) /
    #         (Sigma sqrt{A Sigma}) x (r - r_+ / 4) sqrt{(r_{BL} - r_-) / r }

    KDD[1][2] = KDD[2][1] = -((2*a**3*M*rBL*sp.cos(th)*sp.sin(th)**3)/ \
                    (SIG*sp.sqrt(AA*SIG)))*(r - rp/4)*sp.sqrt((rBL - rm)/r)

    betaU = ixp.zerorank1() # We generally choose \beta^i = 0 for these initial data
    BU    = ixp.zerorank1() # We generally choose B^i = 0 for these initial data

    # Validated against original SENR: KDD[0][2], KDD[1][2], gammaDD[2][2], gammaDD[0][0], gammaDD[1][1]
    #print(sp.mathematica_code(gammaDD[1][1]))

    # Finally set alpha. We generally choose alpha = 1/psi**2 (psi = BSSN conformal factor)
    #                    for these initial data
    try:
        cf_type = par.parval_from_str("EvolvedConformalFactor_cf")
    except:
        print("UIUCBlackHole Error: Must set BSSN_quantities::EvolvedConformalFactor_cf;")
        print("                     the lapse is set in terms of the BSSN conformal factor")
        sys.exit(1)

    try:
        cf_type = par.parval_from_str("EvolvedConformalFactor_cf")
    except:
        print("UIUCBlackHole Error: Must set BSSN_quantities::EvolvedConformalFactor_cf;")
        print("                     the lapse is set in terms of the BSSN conformal factor")
        sys.exit(1)

    rfm.reference_metric()  # BitoA.cf_from_gammaDD requires reference_metric() first be called.
    BitoA.cf_from_gammaDD(gammaDD)
    cf = BitoA.cf

    # Let's choose alpha = 1/psi**2 (psi = BSSN conformal factor) for these initial data,
    # where psi = exp(phi); chi = 1/psi**4; W = 1/psi**2
    if cf_type == "phi":
        alpha = sp.exp(-2*cf)
    elif cf_type == "chi":
        alpha = sp.sqrt(cf)
    elif cf_type == "W":
        alpha = cf
    else:
        print("Error EvolvedConformalFactor_cf type = \""+cf_type+"\" unknown.")
        sys.exit(1)
