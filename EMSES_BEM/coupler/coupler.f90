program coupler
  use mpi
  use ctca
  implicit none

  integer(kind=4) :: ierr, myrank, nprocs, areaid(3), fromrank, progid
  integer(kind=4) :: reqinfo(4)
  integer(kind=4) :: dataint(2)
  integer(kind=4) :: world_myrank
  logical :: CTCAC_verbose = .false.

  call CTCAC_init()   ! couplerの初期化

  call MPI_Comm_size(CTCA_subcomm, nprocs, ierr)
  call MPI_Comm_rank(CTCA_subcomm, myrank, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, world_myrank, ierr)

  ! if (myrank.eq.0) print *, world_myrank, "coupler init done"

  call CTCAC_regarea_real8(areaid(1))
  call CTCAC_regarea_real8(areaid(2))
  call CTCAC_regarea_real8(areaid(3))

  ! メインループ
  do while (.true.)
    call CTCAC_pollreq(reqinfo, fromrank, dataint, 2) ! リクエスト問い合わせ
    ! print *, "coupler poll req done"

    ! リクエスターでCoToCoA 終了か確認
    if (CTCAC_isfin()) then
      ! print *, "coupler got fin"
      exit
    end if

    if (fromrank >= 0) then
      ! if(CTCAC_verbose) print *, "coupler got req"

      progid = 0

      ! enqueue request with data
      call CTCAC_enqreq(reqinfo, progid, dataint, 2)  ! リクエストをキューに追加
      ! if(CTCAC_verbose) print *, "coupler enqreq done"
    end if
  end do

  call CTCAC_finalize()

end program coupler

