      module null_newsvarnorth
	use null_grid

  double complex, allocatable,dimension(:,:), save :: j_lnn, j_lon, j_lom1n, jom1n
  double precision, allocatable,dimension(:,:), save :: alphaom1n, alphaon, nwo, nwom1
  double precision, allocatable,dimension(:,:), save ::  alphann, nwn 
  double complex, allocatable,dimension(:,:,:), save :: jjn

contains
  subroutine null_newsnorth_allocate

    use null_grid, only : nn, nx


  allocate ( alphaon(nn,nn),   alphann(nn,nn), alphaom1n(nn,nn),&
       &           jjn(nn,nn,nx), nwo(nn,nn),   nwom1(nn,nn), nwn(nn,nn),&
       &            j_lnn(nn,nn),  j_lon(nn,nn), j_lom1n(nn,nn),jom1n(nn,nn))

  end subroutine null_newsnorth_allocate

      end module null_newsvarnorth
