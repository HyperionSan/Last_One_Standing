      module null_newsvarsouth
	use null_grid

  double complex, allocatable,dimension(:,:), save :: j_lns, j_los, j_lom1s, jom1s
  double precision, allocatable,dimension(:,:), save :: alphaom1s, alphaos, swo, swom1
  double precision, allocatable,dimension(:,:), save ::  alphans, swn
  double complex, allocatable,dimension(:,:,:), save :: jjs

contains
  subroutine null_newssouth_allocate

    use null_grid, only : nn, nx

  allocate ( alphaos(nn,nn),   alphans(nn,nn), alphaom1s(nn,nn),&
       &       jjs(nn,nn,nx), swo(nn,nn), swom1(nn,nn), swn(nn,nn),&
       &       j_lns(nn,nn),  j_los(nn,nn), j_lom1s(nn,nn),jom1s(nn,nn))

  end subroutine null_newssouth_allocate

      end module null_newsvarsouth
