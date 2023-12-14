Brief description of included thorns

**NSNS_parameter_files**
Contains parameter files for magnetized and unmagnetized BNS evolutions.

**Seed_Magnetic_Fields_BNS**
Extended Seed_Magnetic_Fields thorn for binary neutron stars.

**VolumeIntegrals_GRMHD**
Nice GRMHD volume integration thorn, depends only on ADMBase & HydroBase. Performs volume integrals on arbitrary "Swiss-cheese"-like topologies, and even interoperates with Carpet to track NS centers of mass.

**VolumeIntegrals_vacuum**
Nice GRMHD volume integration thorn, currently depends on ML_BSSN. There is a bit of code duplication and duplicated functionality between VI_GRMHD and VI_vacuum, to ensure that VI_vacuum can be used without enabling a GRMHD code. There is probably a better way of doing this, but I haven't had the time to think deeply about this.

**particle_tracerET**
Solves the ODE
\partial_t x^i = v^i
for typically thousands of tracer particles, using an RK4 integration atop the current timestepping. E.g., one RK4 substep in the particle integration might occur every 16 RK4 substeps in the GRMHD evolution. These tracer particle positions are quite useful for visualizing magnetic field lines in a consistent way from frame-to-frame in a movie (recall that in the GRMHD approximation, the magnetic field lines stay attached the fluid elements they thread ["Flux Freezing"]). Note that the velocity must be consistent with the velocity appearing in the GRMHD induction equation. This thorn reads in the HydroBase vel[] vector gridfunction, which assumes the Valencia formalism, and converts it into the induction equation velocity.

**smallbPoynET**
Computes b^i, b^2, and three spatial components of Poynting flux. It also computes (-1-u_0), which is useful for tracking unbound matter.