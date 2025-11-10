module global
    implicit none
    ! 時間
    real(8), parameter ::dt = 0.01d0
    integer, parameter :: steps = 10000

    ! 格子数
    !integer, parameter :: Nx = 100
    !integer, parameter :: Ny = 100

    ! 計算領域サイズ
    real(8), parameter :: Lx = 2.0d0 * acos(-1.0d0)
    real(8), parameter :: Ly = 2.0d0 * acos(-1.0d0)
    real(8) :: dx 
    real(8) :: dy 
    ! 粘性係数μ
    real(8), parameter :: mu = 1.0d0

end module global

module subprogs
    use global
    implicit none
    contains

    subroutine viscosity_cal(u, v, B, Nx, Ny) ! 粘性項　←　格子の取り方合ってる？
        integer, intent(in) :: Nx, Ny
        real(8), intent(in) :: u(Nx+2, Ny+2), v(Nx+2, Ny+2)
        real(8), intent(out) :: B(2, Nx+2, Ny+2)
        real(8) :: dudxx, dudyy, dvdxx, dvdyy
        real(8) :: re_dx, re_dy
        integer :: i, j
        re_dx = 1.0d0 / dx
        re_dy = 1.0d0 / dy
        do j = 1, Ny
            do i = 1, Nx

                ! dudxx = (-(-u(i-1, j)+u(i, j))*re_dx + (-u(i, j)+u(i+1, j))*re_dx) * re_dx
                ! dudyy = (-(-u(i, j-1)+u(i, j))*re_dy + (-u(i, j)+u(i, j+1))*re_dy) * re_dy
                ! dvdxx = (-(-v(i-1, j)+v(i, j))*re_dx + (-v(i, j)+v(i+1, j))*re_dx) * re_dx
                ! dvdyy = (-(-v(i, j-1)+v(i, j))*re_dy + (-v(i, j)+v(i, j+1))*re_dy) * re_dy
                dudxx = (u(i-1, j) - 2*u(i, j) + u(i+1, j)) / dx**2
                dudyy = (u(i, j-1) - 2*u(i, j) + u(i, j+1)) / dy**2
                dvdxx = (v(i-1, j) - 2*v(i, j) + v(i+1, j)) / dx**2
                dvdyy = (v(i, j-1) - 2*v(i, j) + v(i, j+1)) / dy**2

                B(1, i, j) = mu * (dudxx + dudyy)
                B(2, i, j) = mu * (dvdxx + dvdyy)

            enddo
        enddo

    end subroutine viscosity_cal
    
end module subprogs

program main
    use global
    use subprogs
    implicit none
    real(8) :: t, error_x, error_y
    real(8), allocatable :: x(:, :, :)
    real(8), allocatable :: u(:, :), v(:, :)
    real(8), allocatable :: B(:, :, :), B_exact(:, :, :)
    integer :: exp, Nx, Ny, N, i, j

    open(50, file = '/home/shishida/SMAC/.d/err_viscosity.d') ! 出力ファイルを開く

    do exp = 3, 7 ! Nx, Ny = 8 ~ 128
        Nx = 2 ** exp + 1
        Ny = 2 ** exp + 1
        N = (Nx-1) * (Ny-1)

        if (allocated(x)) deallocate(x)
        if (allocated(u)) deallocate(u)
        if (allocated(v)) deallocate(v)
        if (allocated(B)) deallocate(B)
        if (allocated(B_exact)) deallocate(B_exact)

        allocate(x(2, Nx+2, Ny+2))
        allocate(u(Nx+2, Ny+2))
        allocate(v(Nx+2, Ny+2))
        allocate(B(2, Nx+2, Ny+2))
        allocate(B_exact(2, Nx+2, Ny+2))

        ! 初期化
        x(:, :, :) = 0.0d0
        u(:, :) = 0.0d0
        v(:, :) = 0.0d0
        B(:, :, :) = 0.0d0
        B_exact(:, :, :) = 0.0d0

        dx = Lx / real(Nx-1,kind=8)
        dy = Ly / real(Ny-1,kind=8)


        ! 格子生成
        do j = 1, Ny
            do i = 1, Nx
                x(1, i, j) = (i-1) * dx
                x(2, i, j) = (j-1) * dy
            enddo
        enddo

        ! 与える速度(Taylor-green)
        do j = 1, Ny
            do i = 1, Nx
                u(i, j) =   sin(real(i-1+0.5, kind=8) * dx) * cos(real(j-1, kind=8) * dy)
                v(i, j) = - cos(real(i-1, kind=8) * dx) * sin(real(j-1+0.5, kind=8) * dy)
            enddo
        enddo

        ! 数値計算
        call viscosity_cal(u, v, B, Nx, Ny)

        ! 理論解
        do j = 1, Ny
            do i = 1, Nx
                B_exact(1, i, j) = -2.0d0 * mu * sin(real(i-1+0.5, kind=8) * dx) * cos(real(j-1, kind=8) * dy)
                B_exact(2, i, j) =  2.0d0 * mu * cos(real(i-1, kind=8) * dx) * sin(real(j-1+0.5, kind=8) * dy)
            enddo
        enddo

        ! RSMEを計算
        error_x = 0.0d0
        error_y = 0.0d0
        do j = 2, Ny-1
            do i = 2, Nx-1
                error_x = error_x + (B(1,i,j)-B_exact(1,i,j))**2
                error_y = error_y + (B(2,i,j)-B_exact(2,i,j))**2
                write(*, *) B(1, i, j), B_exact(1, i, j)
            enddo
        enddo
        error_x = sqrt(error_x / real((Nx-2)*(Ny-2),kind=8))
        error_y = sqrt(error_y / real((Nx-2)*(Ny-2),kind=8))

        write(50, *)  Nx-1, error_x, error_y


    enddo


end program main
