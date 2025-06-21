module commondata
        implicit none

        ! Grid parameters
        integer(kind=4), parameter :: nx=257, ny=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8), parameter :: Ra=1e8           !Rayleigh number
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
        integer(kind=4), parameter :: itc_max = 2000000000 ! Maximum iterations

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

        integer(kind=4) :: inter_x(nx,3), inter_y(ny,3)

end module commondata

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: un(0:8)
    real(kind=8) :: velocitySquare
    real(kind=8) :: dx(nx+1), dy(ny+1)
    real(kind=8), parameter :: constA= 3.2d0

    ! Output initial parameters for verification
    write(*,*) "nx =", nx, ", ny =", ny
    write(*,*) "Ra =", real(Ra)

    ! Initialize iteration count and error
    itc = 0
    errorU = 100.0d0
    errorT = 100.0d0

    ! Compute grid coordinates
    !$omp parallel do default(none) shared(xGrid) private(i)
    do i = 0, nx+1
        xGrid(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do
    !$omp end parallel do

    !$omp parallel do default(none) shared(yGrid) private(j)
    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do
    !$omp end parallel do
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
    xGrid(1:nx) = xGrid(1:nx)-dx0/2.0d0
    yGrid(1:ny) = yGrid(1:ny)-dx0/2.0d0
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

    allocate (f(0:8,nx,ny))
    allocate (f_post(0:8,nx,ny))
    allocate (g(0:4,nx,ny))
    allocate (g_post(0:4,nx,ny))

    ! Initialize flow variables
    rho = rho0
    u = 0.0d0
    v = 0.0d0
    temp = 0.0d0
    utemp=0.0d0
    up = 0.0d0
    vp = 0.0d0

    !$omp parallel do default(none) shared(f,g,temp) private(i,j)
    do j=1,ny
        do i=1,nx
            temp(i,j) = dble(i-1)/dble(nx-1)*(Tcold-Thot)+Thot
        enddo
    enddo
    !$omp end parallel do
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
    !$omp parallel do default(none) shared(f,g,u,v,rho,temp,ex,ey,omega_U,omega_T) private(i,j,alpha,velocitySquare,un)
    do j = 1, ny
        do i = 1, nx
            velocitySquare = u(i,j)*u(i,j)+v(i,j)*v(i,j)
            do alpha = 0,8
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                f(alpha,i,j) = rho(i,j)*omega_U(alpha)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*velocitySquare)
            enddo

            do alpha=0,4
                un(alpha) = u(i,j)*ex(alpha)+v(i,j)*ey(alpha)
                g(alpha,i,j)=omega_T(alpha)*temp(i,j)*(1.0d0+4.0d0*un(alpha))
            end do
        enddo
    enddo
    !$omp end parallel do
    !$omp parallel do default(none) shared(inter_x) private(i)
    do i = 1, nx
        if(i == 1)then
            inter_x(i,:) = (/i, i+1, i+2/)

        elseif(i == nx)then
            inter_x(i,:) = (/i-2, i-1, i/)

        else
            inter_x(i,:) = (/i-1, i, i+1/)
        end if
    enddo
    !$omp end parallel do
    !$omp parallel do default(none) shared(inter_y) private(j)
    do j = 1, ny
        if(j == 1)then
            inter_y(j,:) = (/j, j+1, j+2/)

        elseif(j == ny)then
            inter_y(j,:) = (/j-2, j-1, j/)

        else
            inter_y(j,:) = (/j-1, j, j+1/)
        end if
    enddo
    !$omp end parallel do
    return
end subroutine initial

program main
    use omp_lib
    use commondata
    implicit none
    real(kind=8) :: timestart, timeEnd
    integer(kind=4) :: myMaxThreads

    call OMP_set_num_threads(12)
    myMaxThreads = OMP_get_max_threads()

    write(*,*) "Max Running threads=",myMaxThreads

    timestart = OMP_get_wtime()

    call initial()

    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call interpolate()

        call bounceback_u()

        call bounceback_T()

        call macro_u()

        call macro_t()

        if(MOD(itc,2000).EQ.0) then
            call check()
        endif

    enddo

    timeEnd = OMP_get_wtime()

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
    !$omp parallel do default(none) shared(f,f_post,rho,u,v,Fx,Fy,Snu,Sq,temp,gbeta,dt)private(i,j,alpha,s,m,m_post,meq,fSource)
    do j=1,ny
        do i=1,nx

            m(0) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            m(1) = -4.0d0*f(0,i,j)-f(1,i,j)-f(2,i,j)-f(3,i,j)-f(4,i,j)+2.0d0*(f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j))
            m(2) = 4.0d0*f(0,i,j)-2.0d0*(f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j))+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            m(3) = f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
            m(4) = -2.0d0*f(1,i,j)+2.0d0*f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)
            m(5) = f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
            m(6) = -2.0d0*f(2,i,j)+2.0d0*f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)
            m(7) = f(1,i,j)-f(2,i,j)+f(3,i,j)-f(4,i,j)
            m(8) = f(5,i,j)-f(6,i,j)+f(7,i,j)-f(8,i,j)

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
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)*dt
            enddo

            f_post(0,i,j) = ( m_post(0)-m_post(1)+m_post(2) )/9.0d0
            f_post(1,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0+m_post(3)/6.0d0-m_post(4)/6.0d0 &
                            +m_post(7)*0.25d0
            f_post(2,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            +m_post(5)/6.0d0-m_post(6)/6.0d0-m_post(7)*0.25d0
            f_post(3,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0-m_post(3)/6.0d0+m_post(4)/6.0d0 &
                            +m_post(7)*0.25d0
            f_post(4,i,j) = m_post(0)/9.0d0-m_post(1)/36.0d0-m_post(2)/18.0d0 &
                            -m_post(5)/6.0d0+m_post(6)/6.0d0-m_post(7)*0.25d0
            f_post(5,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0+m_post(8)*0.25d0
            f_post(6,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            +m_post(5)/6.0d0+m_post(6)/12.0d0-m_post(8)*0.25d0
            f_post(7,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0-m_post(3)/6.0d0-m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0+m_post(8)*0.25d0
            f_post(8,i,j) = m_post(0)/9.0d0+m_post(1)/18.0d0+m_post(2)/36.0d0+m_post(3)/6.0d0+m_post(4)/12.0d0 &
                            -m_post(5)/6.0d0-m_post(6)/12.0d0-m_post(8)*0.25d0

            f(0,i,j) = f_post(0,i,j)
        enddo
    enddo
    !$omp end parallel do
    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: n(0:4), n_post(0:4), neq(0:4)
    real(kind=8) :: Q(0:4)
    !$omp parallel do default(none) shared(g,g_post,temp,u,v,sig_k) private(i,j,alpha,Q,n,n_post,neq)
    do j=1,ny
        do i=1,nx
            n(0) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(1) = g(1,i,j)-g(3,i,j)
            n(2) = g(2,i,j)-g(4,i,j)
            n(3) = g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
            n(4) = g(1,i,j)-g(2,i,j)+g(3,i,j)-g(4,i,j)

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

            g_post(0,i,j) = n_post(0)-n_post(3)
            g_post(1,i,j) = n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(2,i,j) = n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0
            g_post(3,i,j) = -n_post(1)/2.0d0+n_post(3)/4.0d0+n_post(4)/4.0d0
            g_post(4,i,j) = -n_post(2)/2.0d0+n_post(3)/4.0d0-n_post(4)/4.0d0

            g(0,i,j) = g_post(0,i,j)
        enddo
    enddo
    !$omp end parallel do
    return
end subroutine collision_T

subroutine interpolate()
    use commondata
    implicit none
    real(kind=8) :: interpolateF2, delta_x, delta_y
    integer(kind=4) :: i, j, alpha
    real(kind=8) :: f0, f1, f2, g0, g1, g2
    !$omp parallel do default(none)shared(f,f_post,inter_x,inter_y,ex,ey,dt0,xGrid,yGrid)private(i,j,alpha,f0,f1,f2,delta_x,delta_y)
    do j = 1, ny
        do i = 1, nx
            do alpha = 1, 8
                delta_x=dble(ex(alpha))*dt0
                delta_y=dble(ey(alpha))*dt0

        f0 = interpolateF2(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
            , f_post(alpha, inter_x(i,1), inter_y(j,1)), f_post(alpha, inter_x(i,1), inter_y(j,2))&
            , f_post(alpha, inter_x(i,1), inter_y(j,3)), yGrid(j))

        f1 = interpolateF2(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
            , f_post(alpha, inter_x(i,2), inter_y(j,1)), f_post(alpha, inter_x(i,2), inter_y(j,2))&
            , f_post(alpha, inter_x(i,2), inter_y(j,3)), yGrid(j))

        f2 = interpolateF2(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
            , f_post(alpha, inter_x(i,3), inter_y(j,1)), f_post(alpha, inter_x(i,3), inter_y(j,2))&
            , f_post(alpha, inter_x(i,3), inter_y(j,3)), yGrid(j))

        f(alpha, i, j) = interpolateF2(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                        xGrid(inter_x(i,3))+delta_x, f0, f1, f2, xGrid(i))

            end do
        enddo
    enddo
    !$omp end parallel do
    !$omp parallel do default(none)shared(g,g_post,inter_x,inter_y,ex,ey,dt0,xGrid,yGrid)private(i,j,alpha,g0,g1,g2,delta_x,delta_y)
    do j = 1, ny
        do i = 1, nx
            do alpha = 1, 4
                delta_x=dble(ex(alpha))*dt0
                delta_y=dble(ey(alpha))*dt0

        g0 = interpolateF2(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
            , g_post(alpha, inter_x(i,1), inter_y(j,1)), g_post(alpha, inter_x(i,1), inter_y(j,2))&
            , g_post(alpha, inter_x(i,1), inter_y(j,3)), yGrid(j))

        g1 = interpolateF2(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
            , g_post(alpha, inter_x(i,2), inter_y(j,1)), g_post(alpha, inter_x(i,2), inter_y(j,2))&
            , g_post(alpha, inter_x(i,2), inter_y(j,3)), yGrid(j))

        g2 = interpolateF2(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
            , g_post(alpha, inter_x(i,3), inter_y(j,1)), g_post(alpha, inter_x(i,3), inter_y(j,2))&
            , g_post(alpha, inter_x(i,3), inter_y(j,3)), yGrid(j))

        g(alpha, i, j) = interpolateF2(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, &
                        xGrid(inter_x(i,3))+delta_x, g0, g1, g2 ,xGrid(i))
            enddo
        end do
    end do
    !$omp end parallel do
end subroutine

pure function interpolateF2(x0, x1, x2, f0, f1, f2, x) result(f_interp)
    implicit none
    real(kind=8), intent(in) :: x0, x1, x2, x, f0, f1, f2
    real(kind=8) :: f_interp

    ! Interpolation formula
    f_interp = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) * f0 + &
               ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) * f1 + &
               ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1)) * f2

    return
end function interpolateF2

subroutine bounceback_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !$omp parallel do default(none) shared(f,f_post) private(j)
    do j=1, ny
        !left side
        f(1,1,j) = f_post(3,1,j)
        f(5,1,j) = f_post(7,1,j)
        f(8,1,j) = f_post(6,1,j)

        !right side
        f(3,nx,j) = f_post(1,nx,j)
        f(6,nx,j) = f_post(8,nx,j)
        f(7,nx,j) = f_post(5,nx,j)
    enddo
    !$omp end parallel do
    !$omp parallel do default(none) shared(f,f_post) private(i)
    do i=1, nx
        !Top side
        f(4,i,ny) = f_post(2,i,ny)
        f(7,i,ny) = f_post(5,i,ny)
        f(8,i,ny) = f_post(6,i,ny)

        !Bottom side
        f(2,i,1) = f_post(4,i,1)
        f(5,i,1) = f_post(7,i,1)
        f(6,i,1) = f_post(8,i,1)
    enddo
    !$omp end parallel do
return
end subroutine bounceback_u

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j

    !$omp parallel do default(none) shared(g,g_post,omega_T) private(j)
    do j=1, ny
        !left
        g(1,1,j) = -g_post(3,1,j) + 2.0d0 * omega_T(3) * Thot
        !!right
        g(3,nx,j) = -g_post(1,nx,j) + 2.0d0 * omega_T(1) * Tcold
    end do
    !$omp end parallel do
    !$omp parallel do default(none) shared(g,g_post) private(j)
    do i=1, nx
        !top
        g(4,i,ny) = g_post(2,i,ny)
        !bottom
        g(2,i,1) = g_post(4,i,1)
    end do
    !$omp end parallel do
    return
end subroutine

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !$omp parallel do default(none) shared(f,rho,u,v,Fx,Fy,dt) private(i,j)
    do j=1, ny
        do i=1, nx
            rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
            u(i,j) = (f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j)+0.5d0*dt*Fx(i,j) )/rho(i,j)
            v(i,j) = (f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j)+0.5d0*dt*Fy(i,j) )/rho(i,j)
        enddo
    enddo
    !$omp end parallel do
    return
end subroutine macro_u

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    !$omp parallel do default(none) shared(g,temp) private(i,j)
    do j=1, ny
        do i=1, nx
            temp(i,j) = g(0,i,j)+g(1,i,j)+g(2,i,j)+g(3,i,j)+g(4,i,j)
        end do
    end do
    !$omp end parallel do
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
    !$omp parallel do default(none) shared(u,up,v,vp) private(i,j) reduction(+:error1,error2)
    do j=1,ny
        do i=1,nx
            error1 = error1+(u(i,j)-up(i,j))**2.0d0+(v(i,j)-vp(i,j))**2.0d0
            error2 = error2+u(i,j)**2.0d0+v(i,j)**2.0d0
            up(i,j) = u(i,j)
            vp(i,j) = v(i,j)
        enddo
    enddo
    !$omp end parallel do
    errorU = dsqrt(error1)/dsqrt(error2)
    !$omp parallel do default(none) shared(temp,utemp) private(i,j) reduction(+:error3,error4)
    do j=1,ny
        do i=1,nx
            error3 = error3+(temp(i,j)-utemp(i,j))**2
            error4 = error4+temp(i,j)*temp(i,j)

            utemp(i,j) = temp(i,j)
        end do
    end do
    !$omp end parallel do
    errorT = dsqrt(error3)/dsqrt(error4)

    write(*,*) itc,' ',errorU,' ',errorT

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
