module procedures
  implicit none
contains

  subroutine linked_list(pos, ngrid, gridmin, gridmax, ll, lirst, rgrid)
    implicit none
    integer*8 :: i, ng, ipx, ipy, ipz
    integer*4, intent(in) :: ngrid
    real*4, intent(out) :: rgrid
    real*4, intent(in) :: gridmin, gridmax
    real*4, dimension(:,:), intent(in) :: pos
    integer*8, dimension(:,:,:), intent(out) :: lirst
    integer*8, dimension(:), intent(out) :: ll

    rgrid = (gridmax - gridmin) / ngrid

    ng = size(pos, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i

    end do

    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
      endif
    end do

  end subroutine linked_list

  subroutine read_unformatted(input_filename, data, weight, np)
    implicit none
    integer*8 :: nrows, ncols
    character(len=500), intent(in) :: input_filename
    integer*8, intent(out) :: np
    real*4, allocatable, dimension(:,:), intent(out) :: data
    real*4, allocatable, dimension(:), intent(out) :: weight

    open(20, file=input_filename, status='old', form='unformatted')
    read(20) nrows
    read(20) ncols
    allocate(data(ncols, nrows))
    allocate(weight(nrows))
    read(20) data
    close(20)
    np = nrows
    if (ncols .ge. 4) then
      weight = data(4, :)
    else
      weight = 1.0
    end if
  end subroutine read_unformatted

  subroutine binning(rmin, rmax, nbin, bin, bin_edges, rwidth)
    implicit none

    integer*4 :: i
    integer*4, intent(in) :: nbin
    real*4, intent(in) :: rmin, rmax
    real*4, intent(out) :: rwidth
    real*4, allocatable, dimension(:), intent(out) :: bin, bin_edges

    allocate(bin(nbin))
    allocate(bin_edges(nbin + 1))

    rwidth = (rmax - rmin) / nbin
    do i = 1, nbin + 1
      bin_edges(i) = rmin + (i - 1) * rwidth
    end do
    do i = 1, nbin
      bin(i) = bin_edges(i + 1) - rwidth / 2.
    end do

  end subroutine binning

end module procedures