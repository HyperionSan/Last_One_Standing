      module null_coortranvars

      implicit none

	double precision, allocatable, dimension(:,:), save :: xnn, ynn
        double precision,allocatable, dimension(:,:), save :: xno, yno
                                                                
	double precision, allocatable, dimension(:,:), save :: xsn, ysn
        double precision,allocatable, dimension(:,:), save :: xso, yso

            	
	double precision, allocatable, dimension(:,:), save :: intu1n, intu2n
	double precision, allocatable, dimension(:,:), save :: intu1s, intu2s

	double complex, allocatable, dimension(:,:), save :: newsout
	integer , allocatable, dimension(:,:), save :: maskscri

	double complex, allocatable, dimension(:,:), save :: newsout_s
	integer , allocatable, dimension(:,:), save :: maskscri_s


      double precision, allocatable, dimension(:,:), save:: utimebon, asfac
      double precision, allocatable, dimension(:,:), save:: utimebon_s, asfac_s

contains
  subroutine null_coordtran_allocate

    use null_grid, only : nn, nx


  allocate (xnn(nn,nn), ynn(nn,nn), xno(nn,nn), yno(nn,nn), &
       &   xsn(nn,nn), ysn(nn,nn), xso(nn,nn), yso(nn,nn), &
       &   intu1n(nn,nn), intu2n(nn,nn),intu1s(nn,nn), intu2s(nn,nn), &
       &   newsout(nn,nn), maskscri(nn,nn),newsout_s(nn,nn),maskscri_s(nn,nn) &
       &   ,utimebon(nn,nn), asfac(nn,nn),utimebon_s(nn,nn), asfac_s(nn,nn) )


  end subroutine null_coordtran_allocate

      end module null_coortranvars
