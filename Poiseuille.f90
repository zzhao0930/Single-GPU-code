module commondata
        implicit none
        ! Grid parameters
        integer(kind=4), parameter :: ny=81, nx=ny*2+1
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Reynolds=100.0d0
        integer(kind=4), parameter :: B = 2          !aspect ratio
        real(kind=8), parameter :: U0 = 0.1d0        !characteristic velocity

        ! Parameters in 0 system
        real(kind=8), parameter :: rho0 = 1.0d0
        real(kind=8) :: dt0                       ! Time step
        real(kind=8) :: dx0                       ! Grid spacing
        real(kind=8) :: viscosity0

        ! Parameters in LB system
        real(kind=8) :: length_LB                 ! Characteristic length
        real(kind=8) :: lengthScale
        real(kind=8) :: viscosity_LB              ! Kinematic viscosity
        real(kind=8) :: tauf                ! Relaxation time
        real(kind=8) :: dt                        !Time step
        real(kind=8) :: ax
        ! Iteration control
        integer(kind=4) :: itc = 0            ! Current iteration count
        integer(kind=4) :: itc_max             ! Maximum iterations

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), yGrid(0:ny+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: up(:,:), vp(:,:)
        real(kind=8), allocatable :: rho(:,:)            ! Density
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)    ! force components

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)  ! Current and post-collision distributions

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:8)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:8), ey(0:8)  ! Lattice velocity directions
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq

        integer(kind=4) :: inter_x(nx,3), inter_y(ny,3)

        real(kind=8) :: errorU = 100.0d0
        real(kind=8), parameter :: epsU=1e-13
        real(kind=8) :: time(5)
        data time/0.1546, 0.3091, 0.6182, 1.2365, 2.4729/
        integer(kind=4) :: itcU(5)
end module commondata

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    real(kind=8) :: dy(ny+1)
    real(kind=8) :: constA

    ! Allocate flow variables
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (rho(nx,ny))
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))

    allocate (f(nx,ny,0:8))
    allocate (f_post(nx,ny,0:8))

    ! Compute grid coordinates
    constA = 3.2d0
    do i = 1, nx
        xGrid(i) = (dble(i) - 1.0d0) / (dble(nx) - 1.0d0)
    end do
    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    xGrid=xGrid*(yGrid(ny+1)-yGrid(1))*dble(B)

    ! Compute grid spacing using array slicing
    dy(1:ny+1) = yGrid(1:ny+1) - yGrid(0:ny)
    dx0 = dy(1)
    dt0 = dx0
    write(*,*) "---in 0 system---"
    write(*,*) "deltaX = ", dx0
    write(*,*) "deltaT = ", dt0

    ! Compute grid spacing in system LB
    length_LB = 1.0d0 / dy(1)
    dt = dt0 * length_LB
    lengthScale = length_LB-1.0d0

    ! Initialize flow variables
    rho = rho0
    u = 0.0d0
    v = 0.0d0

    omega_U(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega_U(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega_U(alpha) = 1.0d0/36.0d0
    enddo

    ! Output initial parameters for verification
    write(*,*) "nx =", nx, ", ny =", ny

    ! Initialize iteration count and error
    itc = 0
    ! Calculate viscosity based on LB unit
    viscosity_LB = U0*lengthScale/Reynolds

    ax = 8.0d0*Reynolds*rho0*viscosity_LB**2.0d0/lengthScale**3.0d0

    itcU = nint(time*(lengthScale/2.0d0)**2.0d0/viscosity_LB)
    itc_max = maxval(itcU)+1
    write(*,*) "---in LB unit---"
    write(*,*) "characteristic length   =", real(length_LB), "l.u."
    write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
    write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
    write(*,*) "    "

    ! Calculate relaxation time
    tauf = viscosity_LB * 3.0d0 + 0.5d0
    write(*,*) "tauf =", real(tauf)

    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf
    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)

    do j = 1, ny
        do i = 1, nx
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(i, j, alpha) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo
        enddo
    enddo

    do i = 1, nx
        if(i == 1)then
            inter_x(i,:) = (/i+1, i, i+2/)

        elseif(i == nx)then
            inter_x(i,:) = (/i-1, i, i-2/)

        else
            inter_x(i,:) = (/i-1, i, i+1/)
        end if
    enddo

    do j = 1, ny
        if(j == 1)then
            inter_y(j,:) = (/j+1, j, j+2/)

        elseif(j == ny)then
            inter_y(j,:) = (/j-1, j, j-2/)

        else
            inter_y(j,:) = (/j-1, j, j+1/)
        end if
    enddo

    return
end subroutine initial

program main
    use commondata
    implicit none
    real(8) :: timestart, timeEnd
    call cpu_time(timestart)

    call initial()

    !$acc data copy(u,v,rho) copyin(xGrid,yGrid,ex,ey,f,inter_x,inter_y,up,vp)&
    !$acc create(f_post,Fx,Fy)
    do while((errorU > epsU).AND.(itc < itc_max))
        itc = itc+1

        call collision_U()

        call interpolate()

        call bounceback_u()

        call macro_u()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

        if(any(itc == itcU)) then
            call getU()
        endif
    enddo
    !$acc end data

    call cpu_time(timeEnd)
    write(*,*)"Time=", timeEnd-timestart, "s"
    write(*,*) "MLUPS = ", real(dble(nx*ny)/1e6*dble(itc)/(timeEnd-timeStart))
    call output_ASCII()

    deallocate(u)
    deallocate(v)
    deallocate(up)
    deallocate(vp)
    deallocate(rho)
    deallocate(f)
    deallocate(f_post)
    stop
end program main

subroutine collision_U()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: alpha
    real(kind=8) :: s(0:8)
    real(kind=8) :: m(0:8)
    real(kind=8) :: m_post(0:8)
    real(kind=8) :: meq(0:8)
    real(kind=8) :: fSource(0:8)

!$acc parallel loop private(m,m_post,s,meq,fSource) gang vector collapse(2)
    do j=1,ny
        do i=1,nx

            m(0) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(1) = -4.0d0*f(i,j,0)-f(i,j,1)-f(i,j,2)-f(i,j,3)-f(i,j,4)+2.0d0*(f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8))
            m(2) = 4.0d0*f(i,j,0)-2.0d0*(f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4))+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            m(3) = f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(4) = -2.0d0*f(i,j,1)+2.0d0*f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)
            m(5) = f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(6) = -2.0d0*f(i,j,2)+2.0d0*f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)
            m(7) = f(i,j,1)-f(i,j,2)+f(i,j,3)-f(i,j,4)
            m(8) = f(i,j,5)-f(i,j,6)+f(i,j,7)-f(i,j,8)

            meq(0) = rho(i,j)
            meq(1) = rho(i,j)*( -2.0d0+3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(2) = rho(i,j)*( 1.0d0-3.0d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
            meq(3) = rho(i,j)*u(i,j)
            meq(4) = -rho(i,j)*u(i,j)
            meq(5) = rho(i,j)*v(i,j)
            meq(6) = -rho(i,j)*v(i,j)
            meq(7) = rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
            meq(8) = rho(i,j)*(u(i,j)*v(i,j))

            s(0) = 0.0d0
            s(1) = Snu
            s(2) = Snu
            s(3) = 0.0d0
            s(4) = Sq
            s(5) = 0.0d0
            s(6) = Sq
            s(7) = Snu
            s(8) = Snu

            Fx(i,j) = ax*rho(i,j)
            Fy(i,j) = 0.0d0

            fSource(0) = 0.0d0
            fSource(1) = (6.0d0-3.0d0*s(1))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(2) = -(6.0d0-3.0d0*s(2))*(u(i,j)*Fx(i,j)+v(i,j)*Fy(i,j))
            fSource(3) = Fx(i,j)
            fSource(4) = -(1.0d0-0.5d0*s(4))*Fx(i,j)
            fSource(5) = Fy(i,j)
            fSource(6) = -(1.0d0-0.5d0*s(6))*Fy(i,j)
            fSource(7) = (2.0d0-s(7))*(u(i,j)*Fx(i,j)-v(i,j)*Fy(i,j))
            fSource(8) = (1-0.5d0*s(8))*(v(i,j)*Fx(i,j)+u(i,j)*Fy(i,j))

            do alpha=0,8
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)*dt
            enddo

            f_post(i,j,0) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
            f_post(i,j,1) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                            +m_post(7)*0.25d0
            f_post(i,j,2) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0
            f_post(i,j,3) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                            +m_post(7)*0.25d0
            f_post(i,j,4) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
            f_post(i,j,5) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
            f_post(i,j,6) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
            f_post(i,j,7) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
            f_post(i,j,8) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0

            f(i,j,0) = f_post(i,j,0)
        enddo
    enddo
!$acc end parallel loop

    return
end subroutine collision_U

subroutine interpolate()
    use commondata
    implicit none
    real(kind=8) :: interpolateF, delta_x, delta_y
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: f0, f1, f2
!$acc routine (interpolateF) seq
!$acc parallel loop present(f,f_post,ex,ey,xGrid,yGrid) gang vector collapse(2)
        do j = 1, ny
            do i = 1, nx
                do alpha = 1, 8
                    delta_x=dble(ex(alpha))*dt0
                    delta_y=dble(ey(alpha))*dt0

            f0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(inter_x(i,1), inter_y(j,1), alpha), f_post(inter_x(i,1), inter_y(j,2)&
                , alpha), f_post(inter_x(i,1), inter_y(j,3), alpha))

            f1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(inter_x(i,2), inter_y(j,1), alpha), f_post(inter_x(i,2), inter_y(j,2)&
                , alpha), f_post(inter_x(i,2), inter_y(j,3), alpha))

            f2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                , yGrid(inter_y(j,2)), f_post(inter_x(i,3), inter_y(j,1), alpha), f_post(inter_x(i,3), inter_y(j,2)&
                , alpha), f_post(inter_x(i,3), inter_y(j,3), alpha))

            f(i, j, alpha) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                            xGrid(inter_x(i,3))+delta_x, xGrid(inter_x(i,2)), f0, f1, f2)

                end do
            enddo
        enddo
!$acc end parallel loop
end subroutine

!!NOTE: consider using compiler-specific directives to suggest inlining if necessary.
pure function interpolateF(x0, x1, x2, x, f0, f1, f2) result(f_interp)
    implicit none
    !$acc routine (interpolateF) seq
    real(kind=8), intent(in) :: x0, x1, x2, x, f0, f1, f2
    real(kind=8) :: f_interp

    ! Interpolation formula
    f_interp = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) * f0 + &
               ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) * f1 + &
               ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1)) * f2

    return
end function interpolateF

subroutine bounceback_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !left
    !$acc parallel loop present(f)
    do j=1,ny
        f(1,j,1) = f(nx,j,1)
        f(1,j,5) = f(nx,j,5)
        f(1,j,8) = f(nx,j,8)

    !right
        f(nx,j,3) = f(1,j,3)
        f(nx,j,6) = f(1,j,6)
        f(nx,j,7) = f(1,j,7)
    end do
    !$acc end parallel

    !$acc parallel loop  present(f,f_post)
    do i=1, nx
        !Bottom side (j=1)
        f(i,1,2) = f_post(i,1,4)
        f(i,1,5) = f_post(i,1,7)
        f(i,1,6) = f_post(i,1,8)

        !Top side
        f(i,ny,4) = f_post(i,ny,2)
        f(i,ny,7) = f_post(i,ny,5)
        f(i,ny,8) = f_post(i,ny,6)
    enddo
    !$acc end parallel
return
end subroutine bounceback_u

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(rho,u,v,f,Fx,Fy) gang vector collapse(2)
    do j=1, ny
        do i=1, nx
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0*dt*Fx(i,j))/rho(i,j)
            v(i,j) = (f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0*dt*Fy(i,j))/rho(i,j)
        enddo
    enddo
!$acc end parallel loop
    return
end subroutine macro_u

subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(8) :: error1, error2

    error1 = 0.0d0
    error2 = 0.0d0
!$acc parallel loop default(none) reduction(+:error1,error2) present(u,v,up,vp)
    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))*(u(i,j)-up(i,j))+(v(i,j)-vp(i,j))*(v(i,j)-vp(i,j))
            error2 = error2+u(i,j)*u(i,j)+v(i,j)*v(i,j)

            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo
!$acc end parallel
    errorU = sqrt(error1)/sqrt(error2)

    write(*,*) itc,' ',errorU

    return
end subroutine check

subroutine output_ASCII()
    use commondata
    implicit none
    integer :: i, j
    character(len=100) :: filename

    write(filename,*) ny
    filename = adjustl(filename)

    !!$acc update self(u, v, temp, rho)
    open(unit=01,file='MRT-'//trim(filename)//'.dat',status='unknown')

    write(01,*) 'TITLE="thermal convective flows"'
    write(01,*) 'VARIABLES="X" "Y" "U" "V" "rho"'
    write(01,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(01,100) xGrid(i), yGrid(j), u(i,j), v(i,j), rho(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(01)

    return
end subroutine output_ASCII

subroutine getU()
    use commondata
    implicit none
    integer :: j
    character(len=100) :: filename

     write(filename,*) itc
    filename = adjustl(filename)
    !$acc update self(u)
    open(unit=03,file='u-'//trim(filename)//'.dat',status="unknown")
    do j=1,ny
        write(03,*)  u(nxHalf,j)/u0, yGrid(j)-dx0/2.0d0
    enddo
    close(03)

    return
end subroutine getU
