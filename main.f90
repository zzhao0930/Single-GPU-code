program solution
    implicit none

    integer(kind=4), parameter :: ny=81, n=5
    integer(kind=4) :: i, j, k
    real(kind=8) :: dy(ny+1)
    real(kind=8) :: pi=4.0d0*atan(1.0d0)
    real(kind=8), parameter :: constA=3.2d0
    real(kind=8) :: dt0, dx0, tn
    real(kind=8) :: term_sec
    real(kind=8) :: yGrid(0:ny+1)
    real(kind=8) :: u(ny)
    real(kind=8) :: H
    real(kind=8) :: viscosity
    real(kind=8), parameter :: Reynolds=100.0d0
    real(kind=8), parameter :: u_max=0.1d0
    integer(kind=4) :: itc(n)
    data itc/20788, 41562, 83124, 166262, 332511/
    character(len=100) :: filename

    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    dy(1:ny+1) = yGrid(1:ny+1) - yGrid(0:ny)
    dx0 = dy(1)
    dt0 = dx0

    yGrid = yGrid-dx0/2.d0

    H = (yGrid(ny+1)-yGrid(1))/2.0d0

    viscosity=2.0d0*u_max*H/Reynolds

    do j = 1, n
        tn = itc(j)*dt0*viscosity/H**2.0d0
        write(*,*)"tn=",tn
        write(filename,*) itc(j)
        filename = adjustl(filename)
        open(unit=j,file='Theoretical solution-'//trim(filename)//'.dat',status="unknown")
        do i=1,ny
            term_sec=0.0d0
            do k=0,300
                term_sec=term_sec+((-1)**k*4.0d0)*cos((k+0.5d0)*pi*(yGrid(i)/H-1.0d0))&
                        *exp(-(k+0.5d0)**2*pi**2*tn)/((k+0.5d0)*pi)**3
            end do
            u(i)=(1.0d0-(yGrid(i)/H-1.0d0)**2)-term_sec
            write(j,*)  u(i), yGrid(i)
        enddo
        close(j)
    end do

end program

