module commondata
        implicit none

        integer(kind=4), parameter :: nx=513, ny=nx
        integer(kind=4), parameter :: nxHalf=(nx-1)/2+1, nyHalf=(ny-1)/2+1

        real(kind=8) :: dt                                                       ! Time step
        real(kind=8) :: length_LB                                                ! System conversion coefficient
         real(kind=8), parameter :: Ra=1e6                                      ! Rayleigh number
        real(kind=8), parameter :: Pr=0.71d0                                    ! Prandtl number
        real(kind=8), parameter :: Ma=0.1d0                                     ! Mach number
        real(kind=8) :: kappa                                                    ! Thermal expansion coefficient
        real(kind=8) :: gbeta                                    !Volume expansion coefficient * Gravitational acceleration
        real(kind=8) :: viscosity_LB                                             ! Kinematic viscosity
        real(kind=8), parameter :: Thot=1.0d0, Tcold=0.0d0, Tref=0.0d0

        real(kind=8) :: xGrid(0:nx+1), yGrid(0:ny+1)                             ! Grid parameters
        real(kind=8) :: lengthScale                                              ! Characteristic Length

        real(kind=8), allocatable :: u(:,:), v(:,:)                             ! Velocity components
        real(kind=8), allocatable :: T(:,:)                                     ! temperature
        real(kind=8), allocatable :: P(:,:)                                     ! pressure

        real(kind=8) :: coeff_integral(ny)                                       ! Integral coefficient
end module commondata

program main
    use commondata
    implicit none

    call initial()
    call loadData()
    call dimensionless()
    call getspeed_extreme()
    call calNu()
    call streaming_function()

    deallocate (u)
    deallocate (v)
    deallocate (T)
    deallocate (P)

    write(*,*)"Calculation completed"
    return
end program main

subroutine initial()
    use commondata
    implicit none
    integer(kind=4) :: i, j
    real(kind=8) :: constA
    real(kind=8) :: dy(ny+1)

    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (T(nx,ny))
    allocate (P(nx,ny))

    constA = 3.2d0
    do i = 0, nx+1
        xGrid(i) = 0.5d0 * (erf(constA  * (dble(i) / dble(nx+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do
    do j = 0, ny+1
        yGrid(j) = 0.5d0 * (erf(constA  * (dble(j) / dble(ny+1) - 0.5d0)) / erf(0.5d0 * constA ) + 1.0d0)
    end do

    ! Compute grid spacing using array slicing
    dy(1:ny+1) = yGrid(1:ny+1) - yGrid(0:ny)

    ! Compute grid spacing in system LB
    length_LB = 1.0d0 / dy(1)
    yGrid = yGrid-dy(1)/2.0d0
    xGrid = xGrid-dy(1)/2.0d0

    xGrid=xGrid*length_LB
    yGrid=yGrid*length_LB

    lengthScale = xGrid(nx)+xGrid(1)

    viscosity_LB = (Ma*lengthScale*dsqrt(Pr))/dsqrt(3.0d0*Ra)
    kappa = viscosity_LB/Pr
    gbeta = Ra*viscosity_LB*kappa/lengthScale**3

    return
end subroutine

subroutine loadData()
     use commondata
    implicit none
    integer(kind=4) :: i, j
    integer(kind=4) :: ios
    real(kind=8) :: value1, value2, value3, value4, value5, value6
    character(len=100) :: line

    open(unit=01,file='MRTcavity-7312000.dat',form="formatted", access="sequential", status='old', iostat=ios)

    !Check if the file has been successfully opened
    if (ios /= 0) then
        print *, 'Error opening file!'
        stop
    end if

    do i = 1, 3
        read(01, '(A)') line
    end do

    do j = 1, ny
        do i = 1, nx
            read(01,*) value1, value2, value3, value4, value5, value6
            u(i, j) = value3
            v(i, j) = value4
            T(i, j) = value5
        end do
    end do

    close(01)

    return
end subroutine

subroutine dimensionless()
    use commondata
    implicit none

    u = u/dsqrt(gbeta*lengthScale)
    v = v/dsqrt(gbeta*lengthScale)
    T = (T-Tref)/(Thot-Tcold)

    xGrid = xGrid/lengthScale
    yGrid = yGrid/lengthScale
    lengthScale = xGrid(nx)+xGrid(1)

    return
end subroutine

subroutine getspeed_extreme()
    use commondata
    implicit none

    real(kind=8) :: u_max, v_max

    u_max = maxval(u(nxHalf,:))
    v_max = maxval(v(:,nyHalf))

     open(unit=01,file='Calculation results.dat',position='append')
        write(01,*) "Maximum horizontal velocity along the vertical centerline", u_max
        write(01,*) "Maximum vertical velocity along the horizontal centerline", v_max
    close(01)

    return
end subroutine

subroutine calNu()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: heat_flux(nx,ny)
    real(kind=8) :: integral_x(ny), integral_y
    real(kind=8) :: delta_x1, delta_x2
    real(kind=8) :: Nu_left(ny), Nu_left_max, Nu_left_min, Nu_left_ave
    real(kind=8) :: Nu_mid(ny), Nu_mid_ave
    real(kind=8) :: Nu_vol

!-----------------------------Initialize integral coefficient------------------------------------------------
    do j = 1, ny
        coeff_integral(j) = (yGrid(j+1)-yGrid(j-1))/2.0d0
    end do

!-----------------------------Thermal diffusivity coefficient------------------------------------------------
    do j = 1, ny
        do i = 1, nx
            heat_flux(i,j) = u(i,j)*T(i,j)
        end do
    end do

!--------------------------------Nusselt number on the left wall---------------------------------------------------------
    delta_x1 = xGrid(1)
    delta_x2 = xGrid(2)-xGrid(1)
    do j = 1, ny
        Nu_left(j) = delta_x2*(Thot-delta_x1**2/delta_x2**2*T(2,j)+(delta_x1**2/delta_x2**2-1.0d0)*T(1,j))/&
                    (delta_x1*delta_x2+delta_x1**2)
    end do
    Nu_left_max = maxval(Nu_left)
    Nu_left_min = minval(Nu_left)
    Nu_left_ave = 0.0d0
    do j = 1, ny
        Nu_left_ave = Nu_left_ave+coeff_integral(j)*Nu_left(j)
    end do

!--------------------------------Midline Nusselt number--------------------------------------------------------------------
    delta_x1 = xGrid(nxHalf)-xGrid(nxHalf-1)
    delta_x2 = xGrid(nxHalf+1)-xGrid(nxHalf)
    Nu_mid_ave = 0.0d0
    do j = 1, ny
        Nu_mid(j) = dsqrt(Ra*Pr)*u(nxHalf,j)*T(nxHalf,j)+delta_x2*(T(nxHalf-1,j)-delta_x1**2/delta_x2**2*T(nxHalf+1,j)&
                    +(delta_x1**2/delta_x2**2-1.0d0)*T(nxHalf,j))/(delta_x1*delta_x2+delta_x1**2)
    end do

    do j = 1, ny
        Nu_mid_ave = Nu_mid_ave+coeff_integral(j)*Nu_mid(j)

    end do

!-------------------------------Volume averaged Nusselt number-------------------------------------------------------------
    integral_x = 0.0d0
    integral_y = 0.0d0
    do j = 1, ny
        do i = 1, nx
            integral_x(j) = integral_x(j)+coeff_integral(i)*heat_flux(i,j)
        end do
    end do
    do j = 1, ny
        integral_y = integral_y+coeff_integral(j)*integral_x(j)
    end do
    Nu_vol = dsqrt(Ra*Pr)*integral_y+1.0d0

    open(unit=01,file='Calculation results.dat',position='append')
        write(01,*) "Nu_left_max", Nu_left_max
        write(01,*) "Nu_left_min", Nu_left_min
        write(01,*) "Nu_left_ave", Nu_left_ave
        write(01,*) "Nu_mid_ave", Nu_mid_ave
        write(01,*) "Nu_vol", Nu_vol
    close(01)

    return
end subroutine

subroutine streaming_function()
    use commondata
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: phi(0:nx,ny), phi_max, phi_mid
    real(kind=8) :: phi_max_posx, phi_max_posy

    phi = 0.0d0
    do j = 1, ny
        do i = 1, nx
           phi(i,j) = -v(i,j)*coeff_integral(i)+phi(i-1,j)
        end do
    end do

    phi_max = maxval(dabs(phi))
    do j = 1, ny
        do i = 1, nx
            if(dabs(phi(i,j)) == phi_max)then
                phi_max_posx = xGrid(i)
                phi_max_posy = yGrid(j)
            end if
        end do
    end do

    phi_mid = dabs(phi(nxHalf,nyHalf))

    open(unit=01,file='Calculation results.dat',position='append')
        write(01,*) "Absolute value of flow function at the center", phi_mid
    close(01)
    return
end subroutine
