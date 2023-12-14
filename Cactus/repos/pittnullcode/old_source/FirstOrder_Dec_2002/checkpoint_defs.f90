module checkpoint_defs
  implicit none
 
  logical, save :: checkpointrecover=.false.
  integer, save :: checkpointcount = -1
  logical, save :: alreadydumped = .false.
  character*(13), save :: folder0 = 'checkpoint/0/'
  character*(13), save :: folder1 = 'checkpoint/1/'
  character*(13), save :: folderu
  character*(8), save :: folderr= 'recover/'
  double precision, dimension(:), allocatable :: checkTimeLvls

  integer :: NJFile         = 400
  integer :: NJlFile        = 401
  integer :: NcBFile        = 402
  integer :: NUFile         = 403
  integer :: NcomegaFile    = 404
  integer :: NDeltaFile     = 405
  integer :: NbetaFile      = 406
  integer :: NomegaFile     = 407
  integer :: NuBFile        = 408
  integer :: NzBondiFile    = 409
  integer :: NUyFile        = 410
  integer :: NmaskFile      = 411

  integer :: MJnFile  = 400 
  integer :: MJsFile  = 401 
  integer :: MnunFile = 402 
  integer :: MnusFile = 403 
  integer :: McknFile = 404 
  integer :: McksFile = 405 
  integer :: MbnFile  = 406 
  integer :: MbsFile  = 407 
  integer :: McbnFile = 408 
  integer :: McbsFile = 409 
  integer :: MunFile  = 410 
  integer :: MusFile  = 411
  integer :: MwnFile  = 412 
  integer :: MwsFile  = 413
  integer :: MVarsFile = 414 
  integer :: MMlossFile = 415 
  integer :: MVarsAsciiFile = 416 
  integer :: MTimeLevels = 417 

  integer :: HJFile      = 400
  integer :: HrhoFile    = 401
  integer :: HrhodotFile = 402
  integer :: HomegaFile  = 403
  integer :: HrholFile   = 404
  integer :: HJlFile     = 405
end module checkpoint_defs
