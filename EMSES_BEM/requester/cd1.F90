  subroutine wrt3r_mp(f_id,dsname,myid,nodes,dsize,tab,data,istat,status,odoms)
    use HDF5
    implicit none
    integer(kind=4),parameter :: rank=3
!
    integer(kind=HID_T),intent(in) :: f_id
    character(len=*),intent(in) :: dsname
    integer(kind=4),intent(in) :: myid
    integer(kind=4),intent(in) :: nodes(rank)
    integer(kind=HSIZE_T),intent(in) :: dsize(rank)
    integer(kind=HSIZE_T),intent(in) :: tab
    real(kind=4),intent(in) :: data(-tab:dsize(1)-1+tab,-tab:dsize(2)-1+tab,-tab:dsize(3)-1+tab)
    integer(kind=4),intent(out) :: istat, status
    integer(kind=HSIZE_T),intent(in),optional :: odoms(3,rank)
!
    integer(kind=HID_T) :: s_id, d_id, filespace, memspace, p_id
    integer(kind=HSIZE_T) :: all_data_size(rank)
    integer(kind=HSIZE_T) :: local_data_size(rank)
    integer(kind=HSIZE_T) :: lb(rank)
    integer(kind=HSIZE_T) :: ub(rank)
    integer(kind=HSIZE_T) :: lbg(rank)
    integer(kind=HSIZE_T) :: ubg(rank)
    integer(kind=HSIZE_T) :: stride(rank)
    integer(kind=HSIZE_T) :: offset(rank)
    integer :: coord(rank)
    integer :: i
!
    istat = -1

    coord(1) = mod(myid,nodes(1))
    coord(2) = mod(myid/nodes(1),nodes(2))
    coord(3) = myid/(nodes(1)*nodes(2))

    do i=1,rank
      lbg(i) = 0
      ubg(i) = (dsize(i)-1)*nodes(i)
      stride(i) = 1
      if(present(odoms)) lbg(i) = odoms(1,i)
      if(present(odoms)) ubg(i) = odoms(2,i)
      if(present(odoms)) stride(i) = odoms(3,i)
      all_data_size(i) = ubg(i) - lbg(i) + 1
      all_data_size(i) = all_data_size(i) &
     &                 + modulo(-all_data_size(i),stride(i))
      all_data_size(i) = all_data_size(i)/stride(i)
!
      if(coord(i).eq.0.and.coord(i).eq.nodes(i)-1) then
        offset(i) = 0
        lb(i) = lbg(i)
        ub(i) = ubg(i)
      else if(coord(i).eq.0) then
        offset(i) = 0
        lb(i) = lbg(i)
        ub(i) = min(dsize(i)-2,ubg(i))
      else if(coord(i).eq.nodes(i)-1) then
        offset(i) = max((dsize(i)-1)*coord(i)-lbg(i),0)
        lb(i) = max(modulo(-offset(i),stride(i)),lbg(i)-(dsize(i)-1)*coord(i))
        ub(i) = ubg(i) - (dsize(i) - 1)*coord(i)
      else
        offset(i) = max((dsize(i)-1)*coord(i)-lbg(i),0)
        lb(i) = max(modulo(-offset(i),stride(i)),lbg(i)-(dsize(i)-1)*coord(i))
        ub(i) = min(dsize(i)-2,ubg(i)-(dsize(i)-1)*coord(i))
      end if
      offset(i) = ((dsize(i)-1)*coord(i) - lbg(i) + lb(i))/stride(i)
    end do
    local_data_size(1:rank) = shape(data(lb(1):ub(1):stride(1), &
   &                                     lb(2):ub(2):stride(2), &
   &                                     lb(3):ub(3):stride(3)))

    call h5screate_simple_f(rank,all_data_size,filespace,status)
    call h5dcreate_f(f_id,dsname,H5T_NATIVE_REAL,filespace,d_id,status)
    call h5sclose_f(filespace,status)

    call h5screate_simple_f(rank, local_data_size, memspace, status)
    call h5dget_space_f(d_id, filespace, status)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, local_data_size, status)

    call h5pcreate_f(H5P_DATASET_XFER_F, p_id, status)
    call h5pset_dxpl_mpio_f(p_id, H5FD_MPIO_COLLECTIVE_F, status)

    call h5dwrite_f(d_id, H5T_NATIVE_REAL, &
   &                data(lb(1):ub(1):stride(1),lb(2):ub(2):stride(2),lb(3):ub(3):stride(3)), &
   &                all_data_size, istat, &
   &                file_space_id=filespace, mem_space_id=memspace, xfer_prp=p_id)

    call h5sclose_f(filespace, status)
    call h5sclose_f(memspace, status)
    call h5dclose_f(d_id, status)
    call h5pclose_f(p_id, status)
!
    return
  end subroutine wrt3r_mp
