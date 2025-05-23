module commondata
        implicit none

        ! Grid parameters
        integer(kind=4), parameter :: nx=257, ny=nx, nz=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1, nzHalf=(nz-1)/2+1

        real(kind=8), parameter :: Ra=1e6          !Rayleigh number
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
        integer(kind=4), parameter :: itc_max = 20000 ! Maximum iterations

        ! Convergence criteria
        real(kind=8) :: errorU, errorT              ! Current error
        real(kind=8), parameter :: epsU=1e-6, epsT=1e-6  ! Convergence threshold

        ! Grid coordinates
        real(kind=8) :: xGrid(0:nx+1), yGrid(0:ny+1), zGrid(0:nz+1)

        ! Flow variables
        real(kind=8), allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)      ! Velocity components
        real(kind=8), allocatable :: temp(:,:,:)           ! temperature
        real(kind=8), allocatable :: rho(:,:,:)            ! Density
        real(kind=8), allocatable :: up(:,:,:), vp(:,:,:), wp(:,:,:)    ! Previous velocity components for error checking
        real(kind=8), allocatable :: utemp(:,:,:)          ! Previous temperature for error checking
        real(kind=8), allocatable :: Fz(:,:,:)    ! force components

        ! Distribution functions
        real(kind=8), allocatable :: f(:,:,:,:), f_post(:,:,:,:)  ! Current and post-collision distributions
        real(kind=8),allocatable :: g(:,:,:,:), g_post(:,:,:,:)

        ! MRT relaxation parameters
        real(kind=8) :: omega_U(0:18), omega_T(0:6)  ! Relaxation rates for MRT

        ! Lattice directions
        integer(kind=4) :: ex(0:18), ey(0:18), ez(0:18)  ! Lattice velocity directions
        data ex/ 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0/
        data ey/ 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1/
        data ez/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1, -1,-1, 1, 1,-1,-1/

        ! Additional MRT parameters
        real(kind=8) :: Snu, Sq, Se, sig_k

        integer(kind=4) :: inter_x(nx,3), inter_y(ny,3), inter_z(nz,3)

end module commondata
subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j, K, alpha
    real(kind=8) :: un(0:18)
    real(kind=8) :: us2
    real(kind=8) :: dx(nx+1)
    real(kind=8) :: constA

    ! Output initial parameters for verification
    write(*,*) "nx =", nx, ", ny =", ny,", nz =", nz
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
    do k = 0, nz+1
        zGrid(k) = 0.5d0 * (erf(constA  * (dble(k) / dble(nz+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dx(1:nx+1) = xGrid(1:nx+1) - xGrid(0:nx)
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

    ! Compute grid spacing in system LB
    xGrid(1:nx) = xGrid(1:nx)-dx0/2.0d0
    yGrid(1:ny) = yGrid(1:ny)-dx0/2.0d0
    zGrid(1:nz) = zGrid(1:nz)-dx0/2.0d0
    xGrid=xGrid*length_LB
    yGrid=yGrid*length_LB
    zGrid=zGrid*length_LB
    write(*,*)"dt=", dt
    write(*,*)"xGrid(1)=", xGrid(1)
    ! Calculate relaxation time
    tauf = viscosity_LB * 3.0d0 + 0.5d0
    write(*,*) "tauf =", real(tauf)

    ! Calculate MRT relaxation parameters
    Snu = 1.0d0/tauf

    Sq = 8.0d0*(2.0d0-Snu)/(8.0d0-Snu)

    Se = Snu

    sig_k = 1.0d0/(0.5d0+(1.0d0/Snu-0.5d0)/(0.5d0*Pr))

    ! Allocate flow variables
    allocate (u(nx,ny,nz))
    allocate (v(nx,ny,nz))
    allocate (w(nx,ny,nz))
    allocate (rho(nx,ny,nz))
    allocate (up(nx,ny,nz))
    allocate (vp(nx,ny,nz))
    allocate (wp(nx,ny,nz))
    allocate (temp(nx,ny,nz))
    allocate (utemp(nx,ny,nz))
    allocate (Fz(nx,ny,nz))

    allocate (f(nx,ny,nz,0:18))
    allocate (f_post(nx,ny,nz,0:18))
    allocate (g(nx,ny,nz,0:6))
    allocate (g_post(nx,ny,nz,0:6))

    ! Initialize flow variables
    rho = rho0
    temp = 0.0d0
    utemp=0.0d0
    u = 0.0d0
    v = 0.0d0
    up = 0.0d0
    vp = 0.0d0

    do k=1,nz
        do j=1,ny
            do i=1,nx
                temp(i,j,k) = dble(j-1)/dble(ny-1)*(Tcold-Thot)+Thot
            enddo
        enddo
    end do

    omega_U(0)=1.0d0/3.0d0
    do alpha=1,6
        omega_U(alpha)=1.0d0/18.0d0
    end do
    do alpha=7,18
        omega_U(alpha)=1.0d0/36.0d0
    end do

    omega_T(0) = 1.0d0/2.0d0
    do alpha=1,6
        omega_T(alpha) = 1.0d0/12.0d0
    enddo

    do k=1,nz
        do j=1,ny
            do i=1,nx
                us2 = u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)
                do alpha=0,18
                    un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                    f(i,j,k,alpha) = omega_U(alpha)*rho(i,j,k)*(1.0d0+3.0d0*un(alpha)+4.5d0*un(alpha)*un(alpha)-1.5d0*us2)
                enddo

                do alpha=0,6
                    un(alpha) = u(i,j,k)*ex(alpha)+v(i,j,k)*ey(alpha)+w(i,j,k)*ez(alpha)
                    g(i,j,k,alpha)=omega_T(alpha)*temp(i,j,k)*(1.0d0+6.0d0*un(alpha))
                end do
            enddo
        enddo
    end do

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

    do k = 1, nz
        if(k == 1)then
            inter_z(k,:) = (/k+1, k, k+2/)
        elseif(k == nz)then
            inter_z(k,:) = (/k-1, k, k-2/)
        else
            inter_z(k,:) = (/k-1, k, k+1/)
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
!$acc data copy(u,v,w,rho,temp) copyin(xGrid,yGrid,zGrid,ex,ey,ez,omega_U,omega_T,f,g,inter) &
!$acc create(g_post,f_post,up,vp,wp,utemp,Fz)
    do while(((errorU > epsU).or.(errorT > epsT)).AND.(itc < itc_max))

        itc = itc+1

        call collision_U()

        call collision_T()

        call interploate()

        call bounceback_u()

        call bounceback_T()

        call macro_u()

        call macro_t()
        !$acc wait(1,2)
        if(MOD(itc,2000).EQ.0) then
            call check()
        endif
    enddo
    !$acc end data
    call cpu_time(timeEnd)
    write(*,*)"Time=", timeEnd-timestart, "s"
    write(*,*) "MLUPS = ", real(dble(nx*ny*nz)/1e6*dble(itc)/(timeEnd-timeStart))
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
    deallocate(Fz)
    stop
end program main

subroutine collision_U()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
    integer(kind=4) :: alpha
    real(kind=8) :: m(0:18), m_post(0:18), meq(0:18)
    real(kind=8) :: s(0:18)
    real(kind=8) :: fSource(0:18)
    real(kind=8) :: Vsquare

!$acc parallel loop private(m,m_post,s,meq,fSource) gang vector collapse(3) async(1)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                Vsquare=u(i,j,k)**2.0d0+v(i,j,k)**2.0d0+w(i,j,k)**2.0d0

                m(0)=f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)&
                    &+f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)&
                    &+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)

                m(1)=-30.0d0*f(i,j,k,0)-11.0d0*(f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+&
                    &f(i,j,k,5)+f(i,j,k,6))+8.0d0*(f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+&
                    &f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+&
                    &f(i,j,k,17)+f(i,j,k,18))

                m(2)=12.0d0*f(i,j,k,0)-4.0d0*(f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+&
                    &f(i,j,k,5)+f(i,j,k,6))+(f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+&
                    &f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+&
                    &f(i,j,k,17)+f(i,j,k,18))

                m(3)=f(i,j,k,1)-f(i,j,k,2)+f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,11)&
                    &-f(i,j,k,12)+f(i,j,k,13)-f(i,j,k,14)

                m(4)=4.0d0*(-f(i,j,k,1)+f(i,j,k,2))+f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)&
                    &+f(i,j,k,11)-f(i,j,k,12)+f(i,j,k,13)-f(i,j,k,14)

                m(5)=f(i,j,k,3)-f(i,j,k,4)+f(i,j,k,7)+f(i,j,k,8)-f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,15)&
                    &-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)

                m(6)=-4.0d0*(f(i,j,k,3)-f(i,j,k,4))+f(i,j,k,7)+f(i,j,k,8)-f(i,j,k,9)-f(i,j,k,10)+f(i,j,k,15)&
                    &-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)

                m(7)=f(i,j,k,5)-f(i,j,k,6)+f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)&
                    &+f(i,j,k,15)+f(i,j,k,16)-f(i,j,k,17)-f(i,j,k,18)

                m(8)=-4.0d0*(f(i,j,k,5)-f(i,j,k,6))+f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)&
                    &+f(i,j,k,15)+f(i,j,k,16)-f(i,j,k,17)-f(i,j,k,18)

                m(9)=2.0d0*(f(i,j,k,1)+f(i,j,k,2))-(f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6))&
                    &+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)&
                    &+f(i,j,k,14)-2.0d0*(f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18))

                m(10)=-4.0d0*(f(i,j,k,1)+f(i,j,k,2))+2.0d0*(f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)+f(i,j,k,6))&
                    &+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)+f(i,j,k,12)+f(i,j,k,13)&
                    &+f(i,j,k,14)-2.0d0*(f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18))

                m(11)=f(i,j,k,3)+f(i,j,k,4)-f(i,j,k,5)-f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)&
                    &+f(i,j,k,10)-f(i,j,k,11)-f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)

                m(12)=-2.0d0*(f(i,j,k,3)+f(i,j,k,4)-f(i,j,k,5)-f(i,j,k,6))+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)&
                     &+f(i,j,k,10)-f(i,j,k,11)-f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)

                m(13)=f(i,j,k,7)-f(i,j,k,8)-f(i,j,k,9)+f(i,j,k,10)

                m(14)=f(i,j,k,15)-f(i,j,k,16)-f(i,j,k,17)+f(i,j,k,18)

                m(15)=f(i,j,k,11)-f(i,j,k,12)-f(i,j,k,13)+f(i,j,k,14)

                m(16)=f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)&
                    &-f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)+f(i,j,k,14)

                m(17)=-f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)&
                    &+f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18)

                m(18)=f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)&
                    &-f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)

                meq(0)=rho(i,j,k)
                meq(1)=-11.0d0*rho(i,j,k)+19.0d0*rho(i,j,k)*Vsquare
                meq(2)=3.0d0*rho(i,j,k)-5.5d0*rho(i,j,k)*Vsquare
                meq(3)=rho(i,j,k)*u(i,j,k)
                meq(4)=-2.0d0*rho(i,j,k)*u(i,j,k)/3.0d0
                meq(5)=rho(i,j,k)*v(i,j,k)
                meq(6)=-2.0d0*rho(i,j,k)*v(i,j,k)/3.0d0
                meq(7)=rho(i,j,k)*w(i,j,k)
                meq(8)=-2.0d0*rho(i,j,k)*w(i,j,k)/3.0d0
                meq(9)=3.0d0*rho(i,j,k)*u(i,j,k)**2.0d0-rho(i,j,k)*Vsquare
                meq(10)=-0.5d0*(3.0d0*rho(i,j,k)*u(i,j,k)**2.0d0-rho(i,j,k)*Vsquare)
                meq(11)=rho(i,j,k)*(v(i,j,k)**2.0d0-w(i,j,k)**2.0d0)
                meq(12)=-0.5d0*rho(i,j,k)*(v(i,j,k)**2.0d0-w(i,j,k)**2.0d0)
                meq(13)=rho(i,j,k)*u(i,j,k)*v(i,j,k)
                meq(14)=rho(i,j,k)*v(i,j,k)*w(i,j,k)
                meq(15)=rho(i,j,k)*u(i,j,k)*w(i,j,k)
                meq(16)=0.0d0
                meq(17)=0.0d0
                meq(18)=0.0d0

                s(0)=0.0d0 !-----------S_rho
                s(1)=Snu !------------S_e
                s(2)=Snu !------------S_vapersilon
                s(3)=0.0d0 !------------S_j
                s(4)=Sq !---------------S_q
                s(5)=0.0d0 !----------------S_j
                s(6)=Sq !-----------------S_q
                s(7)=0.0d0 !----------------S_j
                s(8)=Sq !-----------------S_q
                s(9)=Snu !---------------S_nu
                s(10)=Snu !--------------S_Pi
                s(11)=Snu !--------------S_nu
                s(12)=Snu !--------------S_Pi
                s(13)=Snu !--------------S_nu
                s(14)=Snu !--------------S_nu
                s(15)=Snu !--------------S_nu
                s(16)=Sq !----------------S_m
                s(17)=Sq !----------------S_m
                s(18)=Sq !----------------S_m


            Fz(i,j,k)=gbeta*temp(i,j,k)*rho(i,j,k)

            fSource(0) = 0.0d0
            fSource(1) = (38.0d0-19.0d0*s(1))*w(i,j,k)*Fz(i,j,k)
            fSource(2) = -(11.0d0-5.5d0*s(2))*w(i,j,k)*Fz(i,j,k)
            fSource(3) = 0.0d0
            fSource(4) = 0.0d0
            fSource(5) = 0.0d0
            fSource(6) = 0.0d0
            fSource(7) = (1.0d0-0.5d0*s(7))*Fz(i,j,k)
            fSource(8) = -2.0d0*(1.0d0-0.5d0*s(8))*Fz(i,j,k)/3.0d0
            fSource(9) = (-2.0d0+s(9))*w(i,j,k)*Fz(i,j,k)
            fSource(10) = (1.0d0-0.5d0*s(10))*w(i,j,k)*Fz(i,j,k)
            fSource(11) = -(2.0d0-s(11))*w(i,j,k)*Fz(i,j,k)
            fSource(12) = (1.0d0-0.5d0*s(12))*w(i,j,k)*Fz(i,j,k)
            fSource(13) = 0.0d0
            fSource(14) = (1.0d0-0.5d0*s(14))*v(i,j,k)*Fz(i,j,k)
            fSource(15) = (1.0d0-0.5d0*s(15))*u(i,j,k)*Fz(i,j,k)
            fSource(16) = 0.0d0
            fSource(17) = 0.0d0
            fSource(18) = 0.0d0

            do alpha=0,18
                m_post(alpha) = m(alpha)-s(alpha)*(m(alpha)-meq(alpha))+fSource(alpha)*dt
            enddo

            f_post(i,j,k,0)=m_post(0)/19.0d0-5.0d0*m_post(1)/399.0d0+m_post(2)/21.0d0

                f_post(i,j,k,1)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0+m_post(3)/10.0d0&
                    &-m_post(4)/10.0d0+m_post(9)/18.0d0-m_post(10)/18.0d0

                f_post(i,j,k,2)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0-m_post(3)/10.0d0&
                    &+m_post(4)/10.0d0+m_post(9)/18.0d0-m_post(10)/18.0d0

                f_post(i,j,k,3)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0+m_post(5)/10.0d0&
                    &-m_post(6)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0+m_post(11)/12.0d0-m_post(12)/12.0d0

                f_post(i,j,k,4)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0-m_post(5)/10.0d0&
                    &+m_post(6)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0+m_post(11)/12.0d0-m_post(12)/12.0d0

                f_post(i,j,k,5)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0+m_post(7)/10.0d0&
                    &-m_post(8)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0-m_post(11)/12.0d0+m_post(12)/12.0d0

                f_post(i,j,k,6)=m_post(0)/19.0d0-11.0d0*m_post(1)/2394.0d0-m_post(2)/63.0d0-m_post(7)/10.0d0&
                    &+m_post(8)/10.0d0-m_post(9)/36.0d0+m_post(10)/36.0d0-m_post(11)/12.0d0+m_post(12)/12.0d0

                f_post(i,j,k,7)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0+m_post(5)/10.0d0+m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0+m_post(13)/4.0d0+m_post(16)/8.0d0-m_post(17)/8.0d0

                f_post(i,j,k,8)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0+m_post(5)/10.0d0+m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0-m_post(13)/4.0d0-m_post(16)/8.0d0-m_post(17)/8.0d0

                f_post(i,j,k,9)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0-m_post(5)/10.0d0-m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0-m_post(13)/4.0d0+m_post(16)/8.0d0+m_post(17)/8.0d0

                f_post(i,j,k,10)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0-m_post(5)/10.0d0-m_post(6)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &+m_post(11)/12.0d0+m_post(12)/24.0d0+m_post(13)/4.0d0-m_post(16)/8.0d0+m_post(17)/8.0d0

                f_post(i,j,k,11)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0+m_post(15)/4.0d0-m_post(16)/8.0d0+m_post(18)/8.0d0

                f_post(i,j,k,12)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0-m_post(15)/4.0d0+m_post(16)/8.0d0+m_post(18)/8.0d0

                f_post(i,j,k,13)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(3)/10.0d0&
                    &+m_post(4)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0-m_post(15)/4.0d0-m_post(16)/8.0d0-m_post(18)/8.0d0

                f_post(i,j,k,14)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(3)/10.0d0&
                    &-m_post(4)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0+m_post(9)/36.0d0+m_post(10)/72.0d0&
                    &-m_post(11)/12.0d0-m_post(12)/24.0d0+m_post(15)/4.0d0+m_post(16)/8.0d0-m_post(18)/8.0d0

                f_post(i,j,k,15)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(5)/10.0d0&
                    &+m_post(6)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &+m_post(14)/4.0d0+m_post(17)/8.0d0-m_post(18)/8.0d0

                f_post(i,j,k,16)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(5)/10.0d0&
                    &-m_post(6)/40.0d0+m_post(7)/10.0d0+m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &-m_post(14)/4.0d0-m_post(17)/8.0d0-m_post(18)/8.0d0

                f_post(i,j,k,17)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0+m_post(5)/10.0d0&
                    &+m_post(6)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &-m_post(14)/4.0d0+m_post(17)/8.0d0+m_post(18)/8.0d0

                f_post(i,j,k,18)=m_post(0)/19.0d0+4.0d0*m_post(1)/1197.0d0+m_post(2)/252.0d0-m_post(5)/10.0d0&
                    &-m_post(6)/40.0d0-m_post(7)/10.0d0-m_post(8)/40.0d0-m_post(9)/18.0d0-m_post(10)/36.0d0&
                    &+m_post(14)/4.0d0-m_post(17)/8.0d0+m_post(18)/8.0d0
            end do
        end do
    end do
!$acc end parallel loop

!$acc parallel loop gang vector collapse(3) async(1)
    do k = 1, nz
       do j = 1, ny
            do i = 1, nx
                f(i,j,k,0)=f_post(i,j,k,0)
            enddo
        enddo
    end do
!$acc end parallel loop

    return
end subroutine collision_U

subroutine collision_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, alpha
    real(kind=8) :: n(0:6), n_post(0:6), neq(0:6)
    real(kind=8) :: Q(0:6)

!$acc parallel loop private(n,n_post,Q,neq) gang vector collapse(3) async(2)
     do k=1,nz
        do j=1,ny
            do i=1,nx
                n(0) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
                n(1) = g(i,j,k,1)-g(i,j,k,2)
                n(2) = g(i,j,k,3)-g(i,j,k,4)
                n(3) = g(i,j,k,5)-g(i,j,k,6)
                n(4) = g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
                n(5) = g(i,j,k,1)+g(i,j,k,2)-g(i,j,k,3)-g(i,j,k,4)
                n(6) = g(i,j,k,1)+g(i,j,k,2)-g(i,j,k,5)-g(i,j,k,6)

                Q(0) = 1.0d0
                Q(1) = sig_k
                Q(2) = sig_k
                Q(3) = sig_k
                Q(4) = 1.2d0
                Q(5) = 1.2d0
                Q(6) = 1.2d0

                neq(0) = temp(i,j,k)
                neq(1) = temp(i,j,k)*u(i,j,k)
                neq(2) = temp(i,j,k)*v(i,j,k)
                neq(3) = temp(i,j,k)*w(i,j,k)
                neq(4) = temp(i,j,k)*0.5d0
                neq(5) = 0.0d0
                neq(6) = 0.0d0

                do alpha=0,6
                    n_post(alpha)=n(alpha)-Q(alpha)*(n(alpha)-neq(alpha))
                enddo

                g_post(i,j,k,0) = n_post(0)-n_post(4)
                g_post(i,j,k,1) = n_post(1)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0+n_post(6)/6.0d0
                g_post(i,j,k,2) = -n_post(1)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0+n_post(6)/6.0d0
                g_post(i,j,k,3) = n_post(2)/2.0d0+n_post(4)/6.0d0-n_post(5)/3.0d0+n_post(6)/6.0d0
                g_post(i,j,k,4) = -n_post(2)/2.0d0+n_post(4)/6.0d0-n_post(5)/3.0d0+n_post(6)/6.0d0
                g_post(i,j,k,5) = n_post(3)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0-n_post(6)/3.0d0
                g_post(i,j,k,6) = -n_post(3)/2.0d0+n_post(4)/6.0d0+n_post(5)/6.0d0-n_post(6)/3.0d0

            enddo
        enddo
    end do

!$acc end parallel loop

!$acc parallel loop  gang vector collapse(3) async(2)
    do k = 1, nz
       do j = 1, ny
            do i = 1, nx
                g(i,j,k,0)=g_post(i,j,k,0)
            enddo
        enddo
    end do
!$acc end parallel loop

    return
end subroutine collision_T

subroutine interploate()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k, alpha
    real(kind=8) :: delta_x, delta_y, delta_z
    real(kind=8) :: interpolateF
    real(kind=8) :: f0,f1,f2

!$acc routine (interpolateF)
!--------------------------------------------------------------------------------------------------------------------
!$acc parallel loop present(f,f_post,ex,ey,ez,zGrid,yGrid,xGrid) gang vector collapse(3)async(1)
    do k=1, nz
        do j=1, ny
            do i=1, nx
                do alpha=1, 18
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt
                    delta_z=dble(ez(alpha))*dt
!------------------------------------------yoz------------------------------------------------------------------------
    if(alpha==5 .or. alpha==6 .or. alpha==15 .or. alpha==16 .or. alpha==17 .or. alpha==18)then
        f0 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                    , f_post(inter_x(i,2),inter_y(j,1),inter_z(k,1),alpha)&
                    , f_post(inter_x(i,2), inter_y(j,1),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,2),inter_y(j,1),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

        f1 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                    , f_post(inter_x(i,2),inter_y(j,2),inter_z(k,1),alpha)&
                    , f_post(inter_x(i,2), inter_y(j,2),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,2),inter_y(j,2),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

        f2 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                    , f_post(inter_x(i,2),inter_y(j,3),inter_z(k,1),alpha)&
                    , f_post(inter_x(i,2), inter_y(j,3),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,2),inter_y(j,3),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

        f(i,j,k,alpha) = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))&
                    +delta_y, f0, f1, f2, yGrid(inter_y(j,2)))

!-----------------------------------------xoz---------------------------------------------------------------------------
    else if(alpha==1 .or. alpha==2 .or. alpha==11 .or. alpha==12 .or. alpha==13 .or. alpha==14)then
        f0 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                    , f_post(inter_x(i,1),inter_y(j,2),inter_z(k,1),alpha)&
                    , f_post(inter_x(i,1), inter_y(j,2),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,1),inter_y(j,2),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

        f1 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                    , f_post(inter_x(i,2),inter_y(j,2),inter_z(k,1),alpha)&
                    , f_post(inter_x(i,2),inter_y(j,2),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,2),inter_y(j,2),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

        f2 = interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                    , f_post(inter_x(i,3),inter_y(j,2),inter_z(k,1),alpha)&
                    , f_post(inter_x(i,3),inter_y(j,2),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,3),inter_y(j,2),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

        f(i,j,k,alpha) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, xGrid(inter_x(i,3))&
                    +delta_x, f0, f1, f2, xGrid(inter_x(i,2)))

!-----------------------------------------xoy---------------------------------------------------------------------------
    else if(alpha==3 .or. alpha==4 .or. alpha==7 .or. alpha==8 .or. alpha==9 .or. alpha==10)then
        f0 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                    , f_post(inter_x(i,1),inter_y(j,1),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,1),inter_y(j,2),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,1),inter_y(j,3),inter_z(k,2),alpha), yGrid(inter_y(j,2)))

        f1 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                   , f_post(inter_x(i,2),inter_y(j,1),inter_z(k,2),alpha)&
                   , f_post(inter_x(i,2),inter_y(j,2),inter_z(k,2),alpha)&
                   , f_post(inter_x(i,2),inter_y(j,3),inter_z(k,2),alpha), yGrid(inter_y(j,2)))

        f2 = interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                    , f_post(inter_x(i,3),inter_y(j,1),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,3),inter_y(j,2),inter_z(k,2),alpha)&
                    , f_post(inter_x(i,3),inter_y(j,3),inter_z(k,2),alpha), yGrid(inter_y(j,2)))

        f(i,j,k,alpha) = interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, xGrid(inter_x(i,3)) &
                    +delta_x, f0, f1, f2, xGrid(inter_x(i,2)))
    end if
                end do
            end do
        end do
    end do
!$acc end parallel loop

!$acc parallel loop present(g,g_post,ex,ey,ez,zGrid,yGrid,xGrid) gang vector collapse(3)async(2)
    do k=1, nz
        do j=1, ny
            do i=1, nx
                do alpha=1, 6
                    delta_x=dble(ex(alpha))*dt
                    delta_y=dble(ey(alpha))*dt
                    delta_z=dble(ez(alpha))*dt
    if(alpha==3 .or. alpha==4)then
        g(i,j,k,alpha)=interpolateF(yGrid(inter_y(j,1))+delta_y, yGrid(inter_y(j,2))+delta_y, yGrid(inter_y(j,3))+delta_y&
                                    , g_post(inter_x(i,2),inter_y(j,1),inter_z(k,2),alpha)&
                                    , g_post(inter_x(i,2),inter_y(j,2),inter_z(k,2),alpha)&
                                    , g_post(inter_x(i,2),inter_y(j,3),inter_z(k,2),alpha), yGrid(inter_y(j,2)))

    elseif(alpha==5 .or. alpha==6)then
        g(i,j,k,alpha)=interpolateF(zGrid(inter_z(k,1))+delta_z, zGrid(inter_z(k,2))+delta_z, zGrid(inter_z(k,3))+delta_z&
                                    , g_post(inter_x(i,2),inter_y(j,2),inter_z(k,1),alpha)&
                                    , g_post(inter_x(i,2),inter_y(j,2),inter_z(k,2),alpha)&
                                    , g_post(inter_x(i,2),inter_y(j,2),inter_z(k,3),alpha), zGrid(inter_z(k,2)))

    elseif(alpha==1 .or. alpha==2)then
        g(i,j,k,alpha)=interpolateF(xGrid(inter_x(i,1))+delta_x, xGrid(inter_x(i,2))+delta_x, xGrid(inter_x(i,3))+delta_x&
                                    , g_post(inter_x(i,1),inter_y(j,2),inter_z(k,2),alpha)&
                                    , g_post(inter_x(i,2),inter_y(j,2),inter_z(k,2),alpha)&
                                    , g_post(inter_x(i,3),inter_y(j,2),inter_z(k,2),alpha), xGrid(inter_x(i,2)))
    end if
                end do
            end do
        end do
    end do
!$acc end parallel loop
return
end subroutine

pure function interpolateF(x0, x1, x2, f0, f1, f2, x) result(f_interp)
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
    integer(kind=4) :: i, j, k
    !left
    !$acc parallel loop default(none) present(f,f_post) async(1)
    do k=1, nz
        do i=1, nx
            f(i,1,k,3) = f_post(i,1,k,4)
            f(i,1,k,7) = f_post(i,1,k,10)
            f(i,1,k,8) = f_post(i,1,k,9)
            f(i,1,k,15) = f_post(i,1,k,18)
            f(i,1,k,17) = f_post(i,1,k,16)

            f(i,ny,k,4) = f_post(i,ny,k,3)
            f(i,ny,k,9) = f_post(i,ny,k,8)
            f(i,ny,k,10) = f_post(i,ny,k,7)
            f(i,ny,k,16) = f_post(i,ny,k,17)
            f(i,ny,k,18) = f_post(i,ny,k,15)
        end do
    end do
    !$acc end parallel

    !$acc parallel loop default(none) present(f,f_post) async(1)
    do k=1, nz
        do j=1, ny
            !front side
            f(nx,j,k,2) = f_post(nx,j,k,1)
            f(nx,j,k,8) = f_post(nx,j,k,9)
            f(nx,j,k,10) = f_post(nx,j,k,7)
            f(nx,j,k,12) = f_post(nx,j,k,13)
            f(nx,j,k,14) = f_post(nx,j,k,11)

            !back side
            f(1,j,k,1) = f_post(1,j,k,2)
            f(1,j,k,9) = f_post(1,j,k,8)
            f(1,j,k,7) = f_post(1,j,k,10)
            f(1,j,k,13) = f_post(1,j,k,12)
            f(1,j,k,11) = f_post(1,j,k,14)
        end do
    end do
    !$acc end parallel
    !$acc parallel loop default(none) present(f,f_post) async(1)
    do j=1, ny
        do i=1, nx
            !Bottom side (j=0)
            f(i,j,1,5) = f_post(i,j,1,6)
            f(i,j,1,11) = f_post(i,j,1,14)
            f(i,j,1,12) = f_post(i,j,1,13)
            f(i,j,1,15) = f_post(i,j,1,18)
            f(i,j,1,16) = f_post(i,j,1,17)

            !Top side (j=ny)
            f(i,j,nz,6) = f_post(i,j,nz,5)
            f(i,j,nz,13) = f_post(i,j,nz,12)
            f(i,j,nz,14) = f_post(i,j,nz,11)
            f(i,j,nz,17) = f_post(i,j,nz,16)
            f(i,j,nz,18) = f_post(i,j,nz,15)
        end do
    enddo
    !$acc end parallel
    return
end subroutine bounceback_u

subroutine bounceback_T()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k

    !left
    !$acc parallel loop default(none) present(g,g_post,omega_T) async(2)
    do k=1, nz
        do i=1, nx
            g(i,1,k,3) = -g_post(i,1,k,4) + 2.0d0*omega_T(3)*Thot

            g(i,ny,k,4) = -g_post(i,ny,k,3) + 2.0d0*omega_T(4)*Tcold
        end do
    end do
    !$acc end parallel

    !$acc parallel loop default(none) present(g,g_post) async(2)
    do k=1, nz
        do j=1, ny
            !front side
            g(nx,j,k,2) = g_post(nx,j,k,1)

            !back side
            g(1,j,k,1) = g_post(1,j,k,2)
        end do
    end do
    !$acc end parallel
    !$acc parallel loop default(none) present(g,g_post) async(2)
    do j=1, ny
        do i=1, nx
            !Bottom side (j=0)
            g(i,j,1,5) = g_post(i,j,1,6)

            !Top side (j=ny)
            g(i,j,nz,6) = g_post(i,j,nz,5)
        end do
    enddo
    !$acc end parallel
    return
end subroutine

subroutine macro_u()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
!$acc parallel loop default(none) present(rho,u,v,w,f,Fz)  gang vector collapse(3) async(1)
    do k=1,nz
        do j=1,ny
            do i=1,nx
               rho(i,j,k)=f(i,j,k,0)+f(i,j,k,1)+f(i,j,k,2)+f(i,j,k,3)+f(i,j,k,4)+f(i,j,k,5)&
                    &+f(i,j,k,6)+f(i,j,k,7)+f(i,j,k,8)+f(i,j,k,9)+f(i,j,k,10)+f(i,j,k,11)&
                    &+f(i,j,k,12)+f(i,j,k,13)+f(i,j,k,14)+f(i,j,k,15)+f(i,j,k,16)+f(i,j,k,17)+f(i,j,k,18)

                u(i,j,k)=(f(i,j,k,1)-f(i,j,k,2)+f(i,j,k,7)-f(i,j,k,8)+f(i,j,k,9)-f(i,j,k,10)&
                &+f(i,j,k,11)-f(i,j,k,12)+f(i,j,k,13)-f(i,j,k,14))/rho(i,j,k)

                v(i,j,k)=(f(i,j,k,3)-f(i,j,k,4)+f(i,j,k,7)+f(i,j,k,8)-f(i,j,k,9)-f(i,j,k,10)&
                &+f(i,j,k,15)-f(i,j,k,16)+f(i,j,k,17)-f(i,j,k,18))/rho(i,j,k)

                w(i,j,k)=(f(i,j,k,5)-f(i,j,k,6)+f(i,j,k,11)+f(i,j,k,12)-f(i,j,k,13)-f(i,j,k,14)&
                &+f(i,j,k,15)+f(i,j,k,16)-f(i,j,k,17)-f(i,j,k,18)+0.5d0*dt*Fz(i,j,k))/rho(i,j,k)
            end do
        end do
    end do
!$acc end parallel
    return
end subroutine macro_u

subroutine macro_t()
    use commondata
    implicit none
    integer(kind=4) :: i, j, k
!$acc parallel loop default(none) present(temp,g) gang vector collapse(3) async(2)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                temp(i,j,k) = g(i,j,k,0)+g(i,j,k,1)+g(i,j,k,2)+g(i,j,k,3)+g(i,j,k,4)+g(i,j,k,5)+g(i,j,k,6)
            end do
        end do
    end do
!$acc end parallel loop
    return
end subroutine macro_t

subroutine check()
    use commondata
    implicit none
    integer(kind=4) :: i , j, k
    real(kind=8) :: error1, error2,error3, error4
    error1 = 0.0d0
    error2 = 0.0d0
    error3 = 0.0d0
    error4 = 0.0d0

!$acc parallel loop default(none) reduction(+:error1,error2) present(u,v,w,up,vp,wp)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                error1 = error1+(u(i,j,k)-up(i,j,k))**2+(v(i,j,k)-vp(i,j,k))**2+(w(i,j,k)-wp(i,j,k))**2
                error2 = error2+u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k)

                up(i,j,k) = u(i,j,k)
                vp(i,j,k) = v(i,j,k)
                wp(i,j,k) = w(i,j,k)
            end do
        end do
    end do
!$acc end parallel
    errorU=dsqrt(error1)/dsqrt(error2)

!$acc parallel loop default(none) reduction(+:error3,error4) present(temp,utemp)
    do k=1,nz
        do j=1,ny
            do i=1,nx
                error3 = error3+(temp(i,j,k)-utemp(i,j,k))**2
                error4 = error4+temp(i,j,k)*temp(i,j,k)

                utemp(i,j,k) = temp(i,j,k)
            end do
        end do
    end do
!$acc end parallel
    errorT = dsqrt(error3)/dsqrt(error4)

    write(*,*) itc,' ',errorU,' ',errorT

    return
end subroutine check

subroutine output_ASCII
    use commondata
    implicit none
    integer :: i, j, k
    character(len=100) :: filename

    write(filename,*) int(nx)
    filename = adjustl(filename)

    open(unit=02,file="side-heated convection"//trim(filename)//'.dat',status='unknown')
    write(02,*) 'TITLE="side-heated convection "'
    write(02,*) 'VARIABLES="X" "Y" "Z" "U" "V" "W" "RHO" "T" '
    write(02,101) nx, ny ,nz
    do k=1,nz
        do j=1,ny
            do i=1,nx
                write(02,100) xGrid(i), yGrid(j), zGrid(k),u(i,j,k), v(i,j,k), w(i,j,k), rho(i,j,k), temp(i,j,k)
            enddo
        enddo
    end do
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'K=',1x,i5,1x,'F=POINT')
    close(02)
    return
end subroutine output_ASCII
