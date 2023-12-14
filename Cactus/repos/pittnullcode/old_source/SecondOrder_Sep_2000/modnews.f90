	module null_news
	implicit none
	contains

	subroutine news( newsn, nun, nuo, betann, betaon,   & 
     &                           newss, sun, suo, betans, betaos   )

      use null_grid
      use null_eth
      use null_newsvarsouth
      use null_newsvarnorth
      use null_vars, only : jos, jon
      use null_interp
      use calc_facnorth
      use calc_facsouth
      use calc_newsnorth
      use calc_newssouth

! news <-- ::wn, wo , wom1, alphan, alphao, alphaom1
! jscri <-- :: jn, jo, jom1, j_ln, j_lo, j_lom1

! actually alpha stands for delta= alpha - gamma

      implicit none

	double complex, dimension(nn,nn) :: newsn, newss
	double complex, dimension(nn,nn) :: nun, nuo, sun, suo
	double precision, dimension(nn,nn) :: betann, betaon, betans, betaos


	call calcconfactnorth(nun, nuo, betann, betaon)
	call calcconfactsouth(sun, suo, betans, betaos)

	call null_rnsint(nwn, swn)
	call null_rnsint(swn, nwn)
        call null_rnsint(alphann, alphans)
        call null_rnsint(alphans, alphann)


	call calcnewssouth(newss, sun, suo, betans, betaos)
	call calcnewsnorth(newsn, nun, nuo, betann, betaon)

       call null_cnsint(newss, newsn, 2)
       call null_cnsint(newsn, newss, 2)

!	call wt_csmooth(nn,newss)
!	call wt_csmooth(nn,newsn)

	swom1 =swo
	swo = swn
	alphaom1s = alphaos
	alphaos = alphans

	nwom1 =nwo
	nwo = nwn
	alphaom1n = alphaon
	alphaon = alphann



!update *****************************************
	jom1n = jon(:,:,nx)
	j_lom1n = j_lon
	jjn = jon(:,:,:)

	jom1s = jos(:,:,nx)
	j_lom1s = j_los
	jjs = jos(:,:,:)

	
!end update *****************************************



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!_

	return
	end subroutine news

	end module null_news

