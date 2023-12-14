# This module sets up Static Trumpet initial data in terms of
# the variables used in BSSN_RHSs.py

# Authors: Terrence Pierre Jacques, terrencepierrej **at** gmail **dot** com
#          Zachariah B. Etienne, zachetie **at** gmail **dot** com
#          Ian Ruchlin

# ## This module sets up initial data for a static trumpet geometry in spherical coordinates. We can convert from spherical to any coordinate system defined in [reference_metric.py](../edit/reference_metric.py) (e.g., SinhSpherical, Cylindrical, Cartesian, etc.) using the [Exact ADM Spherical-to-BSSNCurvilinear converter module](Tutorial-ADM_Initial_Data-Converting_Exact_ADM_Spherical_or_Cartesian_to_BSSNCurvilinear.ipynb)
#
# ### NRPy+ Source Code for this module: [BSSN/BrillLindquist.py](../edit/BSSN/BrillLindquist.py)
#
# <font color='green'>**All quantities have been validated against the [original SENR code](https://bitbucket.org/zach_etienne/nrpy).**</font>

# ### Here we set up Static Trumpet initial data ([Dennison and Baumgarte, 2014](https://arxiv.org/abs/1403.5484)):
#
# Description of Static Trumpet geometry.
#
# **Inputs for initial data**:
#
# * The black hole mass $M$.
#
# **Additional variables needed for spacetime evolution**:
#
# * Desired coordinate system
# * Desired initial lapse $\alpha$ and shift $\beta^i$. We will choose our gauge conditions as $\alpha=1$ and $\beta^i=B^i=0$. $\alpha = \psi^{-2}$ will yield much better behavior, but the conformal factor $\psi$ depends on the desired *destination* coordinate system (which may not be spherical coordinates).

# Step 1: Initialize core Python/NRPy+ modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

thismodule = __name__

# Input parameters:
M = par.Cparameters("REAL", thismodule, ["M"], [1.0])


def StaticTrumpet():
    global r,th,ph, gammaDD, KDD, alpha, betaU, BU

    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 1: Set psi, the StaticTrumpet conformal factor
    # Dennison and Baumgarte (2014) Eq. 13
    # https://arxiv.org/pdf/1403.5484.pdf

    # psi = sqrt{1 + M/r }
    psi0 = sp.sqrt(1 + M/r)

    # *** The physical spatial metric in spherical basis ***
    # Set the upper-triangle of the matrix...
    # Eq. 15
    # gamma_{ij} = psi^4 * eta_{ij}
    # eta_00 = 1, eta_11 = r^2, eta_22 = r^2 * sin^2 (theta)
    gammaDD = ixp.zerorank2()
    gammaDD[0][0] = psi0**4
    gammaDD[1][1] = psi0**4 * r**2
    gammaDD[2][2] = psi0**4 * r**2*sp.sin(th)**2
    # ... then apply symmetries to get the other components

    # *** The physical trace-free extrinsic curvature in spherical basis ***
    # Set the upper-triangle of the matrix...

    # Eq.19 and 20
    KDD = ixp.zerorank2()

    # K_{rr} = M / r^2
    KDD[0][0] = -M / r**2

    # K_{theta theta} = K_{phi phi} / sin^2 theta = M
    KDD[1][1] = M

    KDD[2][2] = M * sp.sin(th)**2
    # ... then apply symmetries to get the other components

    # Lapse function and shift vector
    # Eq. 15
    # alpha = r / (r+M)
    alpha = r / (r + M)

    betaU = ixp.zerorank1()
    # beta^r = Mr / (r + M)^2
    betaU[0] = M*r / (r + M)**2

    BU    = ixp.zerorank1()
