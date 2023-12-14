#include <array>
#include <vector>
#include <string>
#include <iostream>
#include "importer.h"
#include "mpi.h"

int main(int argc, char **argv) {
  // initialize MPI
  int rc = MPI_Init(&argc, &argv);

  if (rc != MPI_SUCCESS) {
    std::cerr << "Error starting MPI" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  int const npoints = 1000;
  std::vector<double> xx(npoints), yy(npoints), zz(npoints);

	double plim = 30.;
  double mlim = -30.;

  double dx = (plim-mlim)/npoints;
	double dy = (plim-mlim)/npoints;
	double dz = (plim-mlim)/npoints;

  for (int i = 0; i < npoints; ++i) {
    xx[i] = mlim + i*dx;
    yy[i] = mlim + i*dy;
    zz[i] = mlim + i*dz;
    //std::cout << xx[i] << std::endl;
  }
  KadathImporter(npoints, "BNS", xx, yy, zz);
  MPI_Finalize();
  return 0;
}
