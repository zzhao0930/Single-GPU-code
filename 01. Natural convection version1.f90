module commondata
        implicit none

        ! Grid parameters
        integer(kind=4), parameter :: nx=513, ny=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Ra=1e8          !Rayleigh number
        real(kind=8), parameter :: Ma=0.1d0         !Mach number
        real(kind=8), parameter :: Pr=0.71d0        !Prandtl number
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0
        real(kind=8) :: kappa                        !Thermal expansion coefficient
        real(kind=8) :: gbeta                        !Volume expansion coefficient * Gravitational acceleration

        ! Parameters in 0 system
        real(kind=8), parameter :: rho0 = 1.0d0
        real(kind=8) :: length0                   ! Characteristic length
        real(kind=8) :: viscosity0                ! Kinematic viscosity
        real(kind=8) :: dt0                       ! Time step
        real(kind=8) :: dx0                       ! Grid spacing

        ! Parameters in LB system
        real(kind=8) :: length_LB                 ! Characteristic length
        real(kind=8) :: viscosity_LB              ! Kinematic viscosity
        real(kind=8) :: tauf, taut                ! Relaxation time
        real(kind=8) :: dt                       ! Time step

        ! Iteration control
        integer(kind=4) :: itc = 0            ! Current iteration count
        integer(kind=4), parameter :: itc_max = 20000000  ! Maximum iterations

        ! Convergence criteria
        real(kind=8) :: errorU, errorT              ! Current error
        real(kind=8), parameter :: epsU=1e-6, epsT=1e-6  ! Convergence threshold

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), yGrid(0:ny+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:), v(:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:)            ! Density
        real(kind=8), allocatable :: up(:,:), vp(:,:)    ! Previous velocity components for error checking
        real(kind=8), allocatable :: utemp(:,:)          ! Previous temperature for error checking
        real(kind=8), allocatable :: Fx(:,:), Fy(:,:)    ! force components

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:), f_post(:,:,:)  ! Current and post-collision distributions
        real(kind=8),allocatable :: g(:,:,:), g_post(:,:,:)

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:8), omega_T(0:4)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:8), ey(0:8)  ! Lattice velocity directions
        data ex/0, 1, 0, -1,  0, 1, -1, -1,  1/
        data ey/0, 0, 1,  0, -1, 1,  1, -1, -1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq, sig_k

end module commondata

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    real(kind=8) :: dx(nx+1), dy(ny+1)
    real(kind=8) :: constA

    ! Output initial parameters for verification
    write(*,*) "nx =", nx, ", ny =", ny
    write(*,*) "Ra =", real(Ra)

    ! Initialize iteration count and error
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0

    ! Compute grid coordinates
    constA = 3.2d0
    do i = 0, nx+1
        xGrid(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do
    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    write(*,*) "xGrid(0) = ", xGrid(0)
    write(*,*) "xGrid(nx+1) = ", xGrid(nx+1)
    write(*,*) "xGrid(1) = ", xGrid(1)
    write(*,*) "xGrid(nx) = ", xGrid(nx)
    write(*,*) "xGrid(1)+xGrid(nx) = ", xGrid(1)+xGrid(nx)
    write(*,*) "    "

    ! Compute grid spacing using array slicing
    dx(1:nx+1) = xGrid(1:nx+1) - xGrid(0:nx)
    dy(1:ny+1) = yGrid(1:ny+1) - yGrid(0:ny)
    dx0 = dx(1)
    dt0 = dx0
    write(*,*) "---in 0 system---"
    write(*,*) "deltaX = ", dx0
    write(*,*) "deltaT = ", dt0

    ! Calculate characteristic length
    length0 = xGrid(nx+1)-xGrid(1)
    write(*,*) "length0 = ", real(length0)

    ! Calculate viscosity based on Reynolds number in 0 system
    viscosity0 = Ma*length0*dsqrt(Pr)/dsqrt(3.0d0*Ra)

    write(*,*) "viscosity0 = ", real(viscosity0)
    write(*,*) "    "

    ! Calculate viscosity based on LB unit
    length_LB = 1.0d0 / dx(1)
    dt = dt0 * length_LB
    viscosity_LB = (Ma*(length_LB-1.0d0)*dsqrt(Pr))/dsqrt(3.0d0*Ra)
    kappa = viscosity_LB/Pr
    gbeta = Ra*viscosity_LB*kappa/((length_LB-1.0d0)**3)
    write(*,*) "---in LB unit---"
    write(*,*) "characteristic length   =", real(length_LB), "l.u."
    write(*,*) "viscosity_LB =", real(viscosity_LB), "l.u.^2/t.s."
    write(*,*) "timeStep ratio for (uniform) / (non-uniform) : ", real(length_LB / dble(nx))
    write(*,*)  gbeta
    write(*,*) "    "

    write(*,*)"dt=", dt
    write(*,*)"xGrid(1)=", xGrid(1)
    ! Calculate relaxation time
    tauf = viscosity_LB * 3.0d0 + 0.5d0
    write(*,*) "tauf =", real(tauf)

    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf
    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)
    sig_k = 1.0d0/(0.5d0+4.0d0*(1.0d0/Snu-0.5d0)/(3.0d0*Pr))

    ! Allocate flow variables
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (rho(nx,ny))
    allocate (up(nx,ny))
    allocate (vp(nx,ny))
    allocate (temp(nx,ny))
    allocate (utemp(nx,ny))
    allocate (Fx(nx,ny))
    allocate (Fy(nx,ny))

    allocate (f(nx,ny,0:8))
    allocate (f_post(nx,ny,0:8))
    allocate (g(nx,ny,0:4))
    allocate (g_post(nx,ny,0:4))

    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do j=1,ny
        do i=1,nx
            temp(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
        enddo
    enddo


    omega_U(0) = 4.0d0/9.0d0
    do alpha=1,4
        omega_U(alpha) = 1.0d0/9.0d0
    enddo
    do alpha=5,8
        omega_U(alpha) = 1.0d0/36.0d0
    enddo

    omega_T(0) = 1.0d0/2.0d0
    do alpha=1,4
        omega_T(alpha) = 1.0d0/8.0d0
    enddo

    do j = 1, ny
        do i = 1, nx
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(i, j, alpha) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo

            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(i,j,alpha)=omega_T(alpha)*temp(i,j)*(1.0d0+4.0d0*un(alpha))
            end do
        enddo
    enddo

    return
end subroutine initial

program main
    use commondata
    implicit none
    real(8) :: timestart, timeEnd

    call cpu_time(timestart)

    call initial()

!$acc data copy(u,v,rho,temp) copyin(xGrid,yGrid,ex,ey,omega_U,omega_T,f,g) create(g_post,f_post,up,vp,utemp,Fx,Fy)
    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call interpolate()

        call interpolate_boundary()

        call bounceback_u()

        call bounceback_T()

        call macro_uv()

        call macro_t()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif
    enddo
    !$acc end data

    call cpu_time(timeEnd)
    write(*,*)"Time=", timeEnd-timestart, "s"
    write(*,*) "MLUPS = ", real(dble(nx*ny)/1e6*dble(itc)/(timeEnd-timeStart))
    call output_ASCII()

    deallocate(u)
    deallocate(v)
    deallocate(rho)
    deallocate(up)
    deallocate(vp)
    deallocate(f)
    deallocate(f_post)
    deallocate(g)
    deallocate(g_post)
    deallocate(temp)
    deallocate(utemp)

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

!$acc parallel loop default(none) present(f,f_post,u,v,rho,Fx,Fy,temp) private(m,m_post,s,meq,fSource) independent gang vector collapse(2)
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

            Fx(i,j) = 0.0d0
            Fy(i,j) = gbeta*temp(i,j)*rho(i,j)

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
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)
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
        enddo
    enddo
!$acc end parallel loop

!$acc parallel loop default(none) present(f,f_post) independent gang vector collapse(2)
    do j = 1, ny
        do i = 1, nx
            f(i,j,0) = f_post(i,j,0)
        enddo
    enddo
!$acc end parallel loop

    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: Q(0:4)

!$acc parallel loop default(none) present(g,g_post,u,v,temp) private(n,n_post,Q,neq) independent gang vector collapse(2)
    do j=1,ny
        do i=1,nx
            n(0) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(1) = g(i,j,1)-g(i,j,3)
            n(2) = g(i,j,2)-g(i,j,4)
            n(3) = g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
            n(4) = g(i,j,1)-g(i,j,2)+g(i,j,3)-g(i,j,4)

            Q(0) = 1.0d0
            Q(1) = sig_k
            Q(2) = sig_k
            Q(3) = 1.9d0
            Q(4) = 1.9d0

            neq(0) = temp(i,j)
            neq(1) = temp(i,j)*u(i,j)
            neq(2) = temp(i,j)*v(i,j)
            neq(3) = temp(i,j)*0.5d0
            neq(4) = 0.0d0

            do alpha=0,4
                n_post(alpha)=n(alpha)-Q(alpha)*(n(alpha)-neq(alpha))
            enddo

            g_post(i,j,0) = n_post(0)-n_post(3)
            g_post(i,j,1) = n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(i,j,2) = n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0
            g_post(i,j,3) = -n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(i,j,4) = -n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0

        enddo
    enddo
!$acc end parallel loop

!$acc parallel loop default(none) present(g,g_post) independent gang vector collapse(2)
    do j = 1, ny
        do i = 1, nx
            g(i,j,0) = g_post(i,j,0)
        enddo
    enddo
!$acc end parallel loop

    return
end subroutine collision_T

subroutine interpolate()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: interpolateF, f0, f1, f2, g0, g1, g2, delta_x, delta_y

!$acc routine (interpolateF)

!$acc parallel loop default(none) present(f,f_post,ex,ey,xGrid,yGrid) independent gang vector collapse(2)
    do j = 2, ny-1
        do i = 2, nx-1
            do alpha = 1, 8
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                f0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                f_post(i-1, j-1, alpha), f_post(i-1, j, alpha), f_post(i-1, j+1, alpha))

                f1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                f_post(i, j-1, alpha), f_post(i, j, alpha), f_post(i, j+1, alpha))

                f2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                f_post(i+1, j-1, alpha), f_post(i+1, j, alpha), f_post(i+1, j+1, alpha))

                f(i, j, alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                f0, f1, f2)

            end do
        end do
    end do
!$acc end parallel loop
!$acc parallel loop default(none) present(g,g_post,ex,ey,xGrid,yGrid) independent gang vector collapse(2)
    do j = 2, ny-1
        do i = 2, nx-1
            do alpha = 1, 4
                delta_x=dble(ex(alpha)) * dt0
                delta_y=dble(ey(alpha)) * dt0

                g0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                g_post(i-1, j-1, alpha), g_post(i-1, j, alpha), g_post(i-1, j+1, alpha))

                g1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                g_post(i, j-1, alpha), g_post(i, j, alpha), g_post(i, j+1, alpha))

                g2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                g_post(i+1, j-1, alpha), g_post(i+1, j, alpha), g_post(i+1, j+1, alpha))

                g(i, j, alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                g0, g1, g2)
            end do
        end do
    end do
!$acc end parallel loop
return
end subroutine interpolate

subroutine interpolate_upwind()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: interpolateF, f0, f1, f2, g0, g1, g2, delta_x, delta_y
    integer(kind=4) :: x1, x2, y1, y2

!$acc routine (interpolateF)
!$acc parallel loop default(none) present(f, f_post, g, g_post, ex, ey, xGrid, yGrid) independent gang vector collapse(2)
    do j = 3, ny-2
        do i = 3, nx-2
            do alpha = 1, 8
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                x1 = i - ex(alpha)
                y1 = j - ey(alpha)
                x2 = i - 2 * ex(alpha)
                y2 = j - 2 * ey(alpha)

                select case (alpha)
                case (2, 4)
                    f(i, j, alpha) = interpolateF(yGrid(j) + delta_y, yGrid(y1) + delta_y, yGrid(y2) + delta_y, yGrid(j), &
                        f_post(i, j, alpha), f_post(i, y1, alpha), f_post(i, y2, alpha))

                case (1, 3)
                    f(i, j, alpha) = interpolateF(xGrid(i) + delta_x, xGrid(x1) + delta_x, xGrid(x2) + delta_x, xGrid(i), &
                        f_post(i, j, alpha), f_post(x1, j, alpha), f_post(x2, j, alpha))

                case default
                    f0 = interpolateF(yGrid(j) + delta_y, yGrid(y1) + delta_y, yGrid(y2) + delta_y, yGrid(j), &
                        f_post(i, j, alpha), f_post(i, y1, alpha), f_post(i, y2, alpha))

                    f1 = interpolateF(yGrid(j) + delta_y, yGrid(y1) + delta_y, yGrid(y2) + delta_y, yGrid(j), &
                        f_post(x1, j, alpha), f_post(x1, y1, alpha), f_post(x1, y2, alpha))

                    f2 = interpolateF(yGrid(j) + delta_y, yGrid(y1) + delta_y, yGrid(y2) + delta_y, yGrid(j), &
                        f_post(x2, j, alpha), f_post(x2, y1, alpha), f_post(x2, y2, alpha))

                    f(i, j, alpha) = interpolateF(xGrid(i) + delta_x, xGrid(x1) + delta_x, xGrid(x2) + delta_x, xGrid(i), &
                        f0, f1, f2)
                end select
            end do

            do alpha = 1, 4
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                x1 = i - ex(alpha)
                y1 = j - ey(alpha)
                x2 = i - 2 * ex(alpha)
                y2 = j - 2 * ey(alpha)

                select case (alpha)
                case (2, 4)
                    g(i, j, alpha) = interpolateF(yGrid(j) + delta_y, yGrid(y1) + delta_y, yGrid(y2) + delta_y, yGrid(j), &
                        g_post(i, j, alpha), g_post(i, y1, alpha), g_post(i, y2, alpha))

                case (1, 3)
                    g(i, j, alpha) = interpolateF(xGrid(i) + delta_x, xGrid(x1) + delta_x, xGrid(x2) + delta_x, xGrid(i), &
                        g_post(i, j, alpha), g_post(x1, j, alpha), g_post(x2, j, alpha))
                end select
            enddo
        end do
    end do
    !$acc end parallel

!$acc parallel loop default(none) present(f, f_post, g, g_post, ex, ey, xGrid, yGrid)
    do j = 2, ny-1
        ! Loop over the i index with a step of nx-3, targeting near-Boundary-like points (i=2 and i=nx-1).
        do i = 2, nx-1, nx-3
            do alpha = 1, 8
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                f0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            f_post(i-1,j-1,alpha), f_post(i-1,j,alpha), f_post(i-1,j+1,alpha))

                f1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            f_post(i,j-1,alpha), f_post(i,j,alpha), f_post(i,j+1,alpha))

                f2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            f_post(i+1,j-1,alpha),f_post(i+1,j,alpha), f_post(i+1,j+1,alpha))

                f(i,j,alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                            f0, f1, f2)

            end do

            do alpha = 1, 4
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                g0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            g_post(i-1,j-1,alpha), g_post(i-1,j,alpha), g_post(i-1,j+1,alpha))

                g1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            g_post(i,j-1,alpha), g_post(i,j,alpha), g_post(i,j+1,alpha))

                g2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            g_post(i+1,j-1,alpha),g_post(i+1,j,alpha), g_post(i+1,j+1,alpha))

                g(i,j,alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                            g0, g1, g2)

            end do
        end do
    end do
!$acc end parallel loop

!$acc parallel loop default(none) present(f,f_post,g,g_post,ex,ey,yGrid,xGrid)
    do i = 2, nx-1
        ! Loop over the j index with a step of ny-3, targeting near-Boundary-like points (j=2 and j=ny-1).
        do j = 2, ny-1 ,ny-3
            do alpha=1, 8
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                f0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                        f_post(i-1,j-1,alpha), f_post(i-1,j,alpha), f_post(i-1,j+1,alpha))

                f1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                        f_post(i,j-1,alpha), f_post(i,j,alpha), f_post(i,j+1,alpha))

                f2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                        f_post(i+1,j-1,alpha), f_post(i+1,j,alpha), f_post(i+1,j+1,alpha))

                f(i,j,alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                        f0, f1, f2)

            end do

            do alpha = 1, 4
                delta_x = dble(ex(alpha)) * dt0
                delta_y = dble(ey(alpha)) * dt0

                g0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            g_post(i-1,j-1,alpha), g_post(i-1,j,alpha), g_post(i-1,j+1,alpha))

                g1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            g_post(i,j-1,alpha), g_post(i,j,alpha), g_post(i,j+1,alpha))

                g2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                                            g_post(i+1,j-1,alpha),g_post(i+1,j,alpha), g_post(i+1,j+1,alpha))

                g(i,j,alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                            g0, g1, g2)

            end do
        end do
    end do
!$acc end parallel loop
return
end subroutine interpolate_upwind

subroutine interpolate_boundary()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: interpolateF, f0, f1, f2, g0, g1, g2, delta_x, delta_y

!$acc routine (interpolateF)
!-------------------------------------boundary-----------------------------------------------
!$acc parallel loop default(none) present(f, f_post, g, g_post, ex, ey, yGrid, xGrid)
    do j = 2, ny-1
        do alpha = 1, 8

            if(alpha==2 .or. alpha==3 .or. alpha==4 .or. alpha==6 .or. alpha==7) then
            !~ i = 1
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

f0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                            f_post(1,j-1,alpha), f_post(1,j,alpha), f_post(1,j+1,alpha))

f1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                            f_post(2,j-1,alpha), f_post(2,j,alpha), f_post(2,j+1,alpha))

f2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                            f_post(3,j-1,alpha), f_post(3,j,alpha), f_post(3,j+1,alpha))

f(1,j,alpha) = interpolateF(xGrid(1)+delta_x, xGrid(2)+delta_x, xGrid(3)+delta_x, xGrid(1), &
                            f0, f1, f2)
            endif

            if(alpha==1 .or. alpha==2 .or. alpha==4 .or. alpha==5 .or. alpha==8) then
            !~ i = nx
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

f0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                            f_post(nx,j-1,alpha), f_post(nx,j,alpha), f_post(nx,j+1,alpha))

f1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                            f_post(nx-1,j-1,alpha), f_post(nx-1,j,alpha), f_post(nx-1,j+1,alpha))

f2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                            f_post(nx-2,j-1,alpha), f_post(nx-2,j,alpha), f_post(nx-2,j+1,alpha))

f(nx,j,alpha) = interpolateF(xGrid(nx)+delta_x, xGrid(nx-1)+delta_x, xGrid(nx-2)+delta_x, xGrid(nx), &
                            f0, f1, f2)
            endif
        end do
!------------------------Temperature distribution function--------------------------
        do alpha=1,4

            if(alpha==2 .or. alpha==3 .or. alpha==4)then
                !~ i = 1
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

g0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                    g_post(1, j-1, alpha), g_post(1, j, alpha), g_post(1, j+1, alpha))

g1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                    g_post(2, j-1, alpha), g_post(2, j, alpha), g_post(2, j+1, alpha))

g2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                    g_post(3, j-1, alpha), g_post(3, j, alpha),g_post(3, j+1, alpha))

g(1, j, alpha) = interpolateF(xGrid(1)+delta_x, xGrid(2)+delta_x, xGrid(3)+delta_x, xGrid(1), &
                                g0, g1, g2)
            end if

            if(alpha==1 .or.alpha==2 .or. alpha==4)then
                !~ i = nx
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

g0 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                    g_post(nx, j-1, alpha),g_post(nx, j, alpha),g_post(nx, j+1, alpha))

g1 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                    g_post(nx-1, j-1, alpha), g_post(nx-1, j, alpha), g_post(nx-1, j+1, alpha))

g2 = interpolateF(yGrid(j-1)+delta_y, yGrid(j)+delta_y, yGrid(j+1)+delta_y, yGrid(j), &
                    g_post(nx-2, j-1, alpha), g_post(nx-2, j, alpha), g_post(nx-2, j+1, alpha))

g(nx,j,alpha) = interpolateF(xGrid(nx)+delta_x, xGrid(nx-1)+delta_x, xGrid(nx-2)+delta_x, xGrid(nx), &
                            g0, g1, g2)
            end if
        end do
    end do
!$acc end parallel

!$acc parallel loop default(none) present(f, f_post, g, g_post, ex, ey, yGrid, xGrid)
    do i = 2, nx-1
        do alpha = 1, 8

            if(alpha==1 .or. alpha==3 .or. alpha==4 .or. alpha==7 .or. alpha==8) then
            !~ j = 1
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

f0 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                            f_post(i-1,1,alpha), f_post(i-1,2,alpha), f_post(i-1,3,alpha))

f1 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                            f_post(i,1,alpha), f_post(i,2,alpha), f_post(i,3,alpha))

f2 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                            f_post(i+1,1,alpha), f_post(i+1,2,alpha), f_post(i+1,3,alpha))

f(i,1,alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                            f0, f1, f2)
            endif

            if(alpha==1 .or. alpha==2 .or. alpha==3 .or. alpha==5 .or. alpha==6) then
            !~ j = ny
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

f0 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                            f_post(i-1,ny,alpha), f_post(i-1,ny-1,alpha), f_post(i-1,ny-2,alpha))

f1 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                            f_post(i,ny,alpha), f_post(i,ny-1,alpha), f_post(i,ny-2,alpha))

f2 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                            f_post(i+1,ny,alpha), f_post(i+1,ny-1,alpha), f_post(i+1,ny-2,alpha))

f(i,ny,alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                            f0, f1, f2)
            endif
        end do
!---------------------------Temperature distribution function-----------------------------------
        do alpha=1,4

            if(alpha==1 .or. alpha==3 .or.alpha==4)then
                !~ j = 1
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

g0 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                    g_post(i-1, 1, alpha), g_post(i-1, 2, alpha), g_post(i-1, 3, alpha))

g1 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                    g_post(i, 1, alpha), g_post(i, 2, alpha), g_post(i, 3, alpha))

g2 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                    g_post(i+1, 1, alpha), g_post(i+1, 2, alpha), g_post(i+1, 3, alpha))

g(i, 1, alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                    g0, g1, g2)
            end if

            if(alpha==1 .or.alpha==2 .or. alpha==3)then
                !~ j = ny
                delta_x = dble(ex(alpha))*dt0
                delta_y = dble(ey(alpha))*dt0

g0 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                    g_post(i-1, ny, alpha), g_post(i-1, ny-1, alpha),g_post(i-1, ny-2, alpha))

g1 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                    g_post(i, ny, alpha), g_post(i, ny-1, alpha), g_post(i, ny-2, alpha))

g2 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                    g_post(i+1, ny, alpha), g_post(i+1, ny-1, alpha), g_post(i+1, ny-2, alpha))

g(i, ny, alpha) = interpolateF(xGrid(i-1)+delta_x, xGrid(i)+delta_x, xGrid(i+1)+delta_x, xGrid(i), &
                                g0, g1, g2)
            end if
        end do
    end do
!$acc end parallel

!--------------------------------corner-------------------------------------------------
!$acc parallel loop default(none) present(f,f_post,ex,ey,yGrid,xGrid)
do alpha = 1, 8

    if(alpha==3 .or.alpha==4 .or. alpha==7)then
        delta_x = dble(ex(alpha))*dt0
        delta_y = dble(ey(alpha))*dt0

        f0 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                                    f_post(1,1,alpha), f_post(1,2,alpha), f_post(1,3,alpha))

        f1 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                                    f_post(2,1,alpha), f_post(2,2,alpha), f_post(2,3,alpha))

        f2 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                                    f_post(3,1,alpha), f_post(3,2,alpha), f_post(3,3,alpha))

        f(1,1,alpha) = interpolateF(xGrid(1)+delta_x, xGrid(2)+delta_x, xGrid(3)+delta_x, xGrid(1), &
                                    f0, f1, f2)
    end if

    if(alpha==2 .or.alpha==3 .or.alpha==6)then
        delta_x = dble(ex(alpha))*dt0
        delta_y = dble(ey(alpha))*dt0

        f0 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                                    f_post(1,ny,alpha), f_post(1,ny-1,alpha), f_post(1,ny-2,alpha))

        f1 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                                    f_post(2,ny,alpha), f_post(2,ny-1,alpha), f_post(2,ny-2,alpha))

        f2 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                                    f_post(3,ny,alpha), f_post(3,ny-1,alpha), f_post(3,ny-2,alpha))

        f(1,ny,alpha) = interpolateF(xGrid(1)+delta_x, xGrid(2)+delta_x, xGrid(3)+delta_x, xGrid(1), &
                                    f0, f1, f2)
    end if

    if(alpha==1 .or.alpha==2 .or.alpha==5)then
        delta_x = dble(ex(alpha))*dt0
        delta_y = dble(ey(alpha))*dt0

        f0 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                                    f_post(nx,ny,alpha), f_post(nx,ny-1,alpha), f_post(nx,ny-2,alpha))

        f1 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                                    f_post(nx-1,ny,alpha), f_post(nx-1,ny-1,alpha), f_post(nx-1,ny-2,alpha))

        f2 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                                    f_post(nx-2,ny,alpha), f_post(nx-2,ny-1,alpha), f_post(nx-2,ny-2,alpha))

        f(nx,ny,alpha) = interpolateF(xGrid(nx)+delta_x, xGrid(nx-1)+delta_x, xGrid(nx-2)+delta_x, xGrid(nx), &
                                    f0,f1,f2)
    end if

    if (alpha==1 .or.alpha==4 .or.alpha==8)then
        delta_x = dble(ex(alpha))*dt0
        delta_y = dble(ey(alpha))*dt0
        f0 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                                    f_post(nx,1,alpha), f_post(nx,2,alpha), f_post(nx,3,alpha))

        f1 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                                    f_post(nx-1,1,alpha), f_post(nx-1,2,alpha), f_post(nx-1,3,alpha))

        f2 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                                    f_post(nx-2,1,alpha), f_post(nx-2,2,alpha), f_post(nx-2,3,alpha))

        f(nx,1,alpha) = interpolateF(xGrid(nx)+delta_x, xGrid(nx-1)+delta_x, xGrid(nx-2)+delta_x, xGrid(nx), &
                                    f0, f1, f2)
    end if
end do
!$acc end parallel

!------------------------------Temperature distribution function----------------------------
!$acc parallel loop default(none) present(g, g_post, ex, ey, xGrid, yGrid)
    do alpha=1,4

        if(alpha==3 .or.alpha==4)then
            delta_x = dble(ex(alpha))*dt0
            delta_y = dble(ey(alpha))*dt0

        g0=interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                        g_post(1, 1, alpha), g_post(1, 2, alpha), g_post(1, 3, alpha))

        g1=interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                        g_post(2, 1, alpha), g_post(2, 2, alpha), g_post(2, 3, alpha))

        g2=interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                        g_post(3, 1, alpha), g_post(3, 2, alpha), g_post(3, 3, alpha))

        g(1, 1, alpha)=interpolateF(xGrid(1)+delta_x, xGrid(2)+delta_x, xGrid(3)+delta_x, xGrid(1), &
                                    g0, g1, g2)
        end if

        if(alpha==2 .or.alpha==3)then
            delta_x = dble(ex(alpha))*dt0
            delta_y = dble(ey(alpha))*dt0

        g0=interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                        g_post(1, ny, alpha), g_post(1, ny-1, alpha), g_post(1, ny-2, alpha))

        g1=interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                        g_post(2, ny, alpha), g_post(2, ny-1, alpha), g_post(2, ny-2, alpha))

        g2=interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                        g_post(3, ny, alpha), g_post(3, ny-1, alpha), g_post(3, ny-2, alpha))

        g(1, ny, alpha)=interpolateF(xGrid(1)+delta_x, xGrid(2)+delta_x, xGrid(3)+delta_x, xGrid(1), &
                                    g0, g1, g2)
        end if

        if(alpha==1 .or.alpha==2)then
            delta_x = dble(ex(alpha))*dt0
            delta_y = dble(ey(alpha))*dt0

        g0 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                            g_post(nx, ny, alpha), g_post(nx, ny-1, alpha), g_post(nx, ny-2, alpha))

        g1 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                            g_post(nx-1, ny, alpha), g_post(nx-1, ny-1, alpha), g_post(nx-1, ny-2, alpha))

        g2 = interpolateF(yGrid(ny)+delta_y, yGrid(ny-1)+delta_y, yGrid(ny-2)+delta_y, yGrid(ny), &
                            g_post(nx-2, ny, alpha), g_post(nx-2, ny-1, alpha), g_post(nx-2, ny-2, alpha))

        g(nx, ny, alpha) = interpolateF(xGrid(nx)+delta_x, xGrid(nx-1)+delta_x, xGrid(nx-2)+delta_x, xGrid(nx), &
                                    g0, g1, g2)
        end if

        if(alpha==1 .or.alpha==4)then
            delta_x = dble(ex(alpha))*dt0
            delta_y = dble(ey(alpha))*dt0

        g0 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                            g_post(nx, 1, alpha), g_post(nx, 2, alpha), g_post(nx, 3, alpha))

        g1 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                            g_post(nx-1, 1, alpha), g_post(nx-1, 2, alpha),g_post(nx-1, 3, alpha))

        g2 = interpolateF(yGrid(1)+delta_y, yGrid(2)+delta_y, yGrid(3)+delta_y, yGrid(1), &
                            g_post(nx-2, 1, alpha), g_post(nx-2, 2, alpha),g_post(nx-2, 3, alpha))

        g(nx, 1, alpha) =interpolateF(xGrid(nx)+delta_x, xGrid(nx-1)+delta_x, xGrid(nx-2)+delta_x, xGrid(nx), &
                                        g0, g1, g2)
        end if
    end do
!$acc end parallel
return
end subroutine

!!NOTE: consider using compiler-specific directives to suggest inlining if necessary.
pure function interpolateF(x0, x1, x2, x, f0, f1, f2) result(f_interp)
    implicit none
    !$acc routine (interpolateF)
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

    !$acc parallel loop default(none) present(f,f_post)
    do j=1, ny
        !Left side
        f(1,j,1) = f_post(1,j,3)
        f(1,j,5) = f_post(1,j,7)
        f(1,j,8) = f_post(1,j,6)

        !Right side
        f(nx,j,3) = f_post(nx,j,1)
        f(nx,j,6) = f_post(nx,j,8)
        f(nx,j,7) = f_post(nx,j,5)
    enddo
    !$acc end parallel
    !$acc parallel loop default(none) present(f,f_post)
    do i=1, nx
        !Bottom side
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

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(g,g_post,omega_T)
    do j=1, ny
        !left
        g(1,j,1) = -g_post(1,j,3) + 2.0d0 * omega_T(1) * Thot
        !!right
        g(nx,j,3) = -g_post(nx,j,1) + 2.0d0 * omega_T(3) * Tcold
    end do
!$acc end parallel
!$acc parallel loop default(none) present(g,g_post)
    do i=1, nx
        !top
        g(i,ny,4) = g_post(i,ny,2)
        !bottom
        g(i,1,2) = g_post(i,1,4)
    end do
!$acc end parallel
    return
end subroutine

subroutine macro_uv()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(rho,u,v,f,Fx,Fy)independent gang vector collapse(2)
    do j=1, ny
        do i=1, nx
            rho(i,j) = f(i,j,0)+f(i,j,1)+f(i,j,2)+f(i,j,3)+f(i,j,4)+f(i,j,5)+f(i,j,6)+f(i,j,7)+f(i,j,8)
            u(i,j) = (f(i,j,1)-f(i,j,3)+f(i,j,5)-f(i,j,6)-f(i,j,7)+f(i,j,8)+0.5d0  * Fx(i,j) )/rho(i,j)
            v(i,j) = (f(i,j,2)-f(i,j,4)+f(i,j,5)+f(i,j,6)-f(i,j,7)-f(i,j,8)+0.5d0  * Fy(i,j) )/rho(i,j)
        enddo
    enddo
!$acc end parallel loop
    return
end subroutine macro_uv

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j

!$acc parallel loop default(none) present(temp,g)independent gang vector collapse(2)
    do j=1, ny
        do i=1, nx
            temp(i,j) = g(i,j,0)+g(i,j,1)+g(i,j,2)+g(i,j,3)+g(i,j,4)
        end do
    end do
!$acc end parallel loop
    return
end subroutine macro_t

subroutine check()
    use commondata
    implicit none
    integer :: i, j
    real(kind=8) :: error1, error2,error3, error4
    error1 = 0.0d0
    error2 = 0.0d0
    error3 = 0.0d0
    error4 = 0.0d0

 !$acc parallel loop default(none) reduction(+:error1,error2) present(u,v,up,vp)
    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))**2.0d0+(v(i,j)-vp(i,j))**2.0d0
            error2 = error2+u(i,j)**2.0d0+v(i,j)**2.0d0
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo
!$acc end parallel loop

    errorU = dsqrt(error1)/dsqrt(error2)

!$acc parallel loop default(none) reduction(+:error3,error4) present(temp,utemp)
    do j=1,ny
        do i=1,nx
            error3 = error3+(temp(i,j)-utemp(i,j))**2
            error4 = error4+temp(i,j)*temp(i,j)

            utemp(i,j) = temp(i,j)
        end do
    end do
!$acc end parallel

    errorT = dsqrt(error3)/dsqrt(error4)

    write(*,*) itc,' ',errorU,' ',errorT

    open(unit=10,file='errorU.dat',status='unknown',position='append')
        write(10,*) itc, errorU
    close(10)

    open(unit=11,file='errorT.dat',status='unknown',position='append')
        write(11,*) itc, errorT
    close(11)

    return
end subroutine check

subroutine output_ASCII()
    use commondata
    implicit none
    integer :: i, j
    character(len=100) :: filename

    write(filename,*) itc
    filename = adjustl(filename)


    open(unit=02,file='MRTcavity-'//trim(filename)//'.dat',status='unknown')

    write(02,*) 'TITLE="thermal convective flows"'
    write(02,*) 'VARIABLES="X" "Y" "U" "V" "T" "rho"'
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xGrid(i), yGrid(j), u(i,j), v(i,j), temp(i,j), rho(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
end subroutine output_ASCII
