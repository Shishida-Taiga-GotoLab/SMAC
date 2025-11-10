module global
    implicit none
    ! 時間
    real(8), parameter ::dt = 0.01d0
    integer, parameter :: steps = 10000

    ! 格子数
    !integer :: Nx
    !integer :: Ny

    ! 計算領域サイズ
    real(8), parameter :: Lx = 1.0d0
    real(8), parameter :: Ly = 1.0d0
    real(8) :: dx 
    real(8) :: dy
end module global


module subprogs
    implicit none
    contains

    ! subroutine gauss_seidel(BXP, BXM, BYP, BYM, B0, u_pre, phi, x)
    !     use global
    !     real(8), intent(in) :: BXP, BXM, BYP, BYM, B0, u_pre(Nx+1, Ny+1), x(2, Nx+1, Ny+1)
    !     real(8), intent(inout) :: phi(Nx+1, Ny+1)
    !     real(8) :: phi_new(Nx+1, Ny+1), er, er0, RHS, phi_pre
    !     integer :: i, j, itr, itr_max

    !     itr_max = 10000
    !     phi_new(:, :) = 0.0d0
    !     er = 0.0d0
    !     er0 = 10e-6

    !     do itr = 1, itr_max
    !         er = 0.0d0
    !         do j = 2, Ny-1
    !             do i = 2, Nx-1

    !                 RHS = -2 * (x(1, i, j) * (1-x(1, i, j)) + x(2, i, j) * (1-x(2, i, j)))
    !                 !RHS = ((u_pre(i+1,j)-u_pre(i-1,j))/(2*dx)) / dt →いったん簡単な右辺を考える(20251002)
    !                 !phi_new(i, j) = phi(i, j) + &
    !                 !& ( BYM*phi(i, j-1) + BXM*phi(i-1, j) - B0*phi(i, j) &
    !                 !& + BXP*phi(i+1, j) + BYP*phi(i, j+1) - RHS) / B0 →長いので一旦ボツ(20251002)
    !                 phi_pre = phi(i,j)
    !                 phi(i,j) = (BYM*phi(i,j-1) + BXM*phi(i-1,j) + BXP*phi(i+1,j) + BYP*phi(i,j+1) - RHS ) / B0
    !                 er = er + abs(phi(i, j) - phi_pre)
    !             enddo
    !         enddo
    !         !phi(:, :) = phi_new(:, :)
    !         write(*, *) 'itr = ', itr, 'err = ', er
    !         if (er < er0) then
    !             write(*, *) "exit itr = ", itr
    !             exit
    !         endif
            
    !     enddo
    !     end subroutine gauss_seidel

        subroutine SOR(BXP, BXM, BYP, BYM, B0, u_pre, phi, x, Nx, Ny)
        use global
        integer, intent(in) :: Nx, Ny
        real(8), intent(in) :: BXP, BXM, BYP, BYM, B0, u_pre(Nx+1, Ny+1), x(2, Nx+1, Ny+1)
        real(8), intent(inout) :: phi(Nx+1, Ny+1)
        real(8) ::  er, er0, RHS(Nx+1, Ny+1), RHS_sum, RHS_norm, omg, phi_pre, E(Nx+1, Ny+1) ! RHS : 右辺　r : 離散残差
        integer :: i, j, itr, itr_max

        itr_max = 10000
        er = 0.0d0
        er0 = 10e-3
        omg = 1.7d0
        !RHS(i, j) = 0.0d0
        RHS_sum = 0.0d0


        ! 右辺を計算
        do j = 2, Ny-1
            do i = 2, Nx-1
                !RHS(i, j) = - 2.0d0 * (x(1, i, j) * (Lx - x(1, i, j)) + x(2, i, j) * (Ly - x(2, i, j))) ! 解φが簡単な形で表現され、ディリクレ境界条件はすべて0でok
                RHS(i, j) = - 2.0d0 *(acos(-1.0d0) / Lx)**2 * sin(acos(-1.0d0)*x(1, i, j) / Lx) * sin(acos(-1.0d0)*x(2, i, j) / Ly)
                RHS_sum = RHS_sum + RHS(i, j) ** 2
            enddo
        enddo
        RHS_norm = sqrt(RHS_sum / real((Nx-1)*(Ny-1), kind = 8))
        !write(*,*) RHS_norm

        do itr = 1, itr_max
            er = 0.0d0 ! 誤差をリセット

            do j = 2, Ny-1
                do i = 2, Nx-1

                    !RHS = - 2.0d0 * (x(1, i, j) * (Lx - x(1, i, j)) + x(2, i, j) * (Ly - x(2, i, j))) ! 解φが簡単な形で表現され、ディリクレ境界条件はすべて0でok
                    !RHS = ((u_pre(i+1,j)-u_pre(i-1,j))/(2*dx)) / dt →いったん簡単な右辺を考える(20251002)
                    !RHS_sum = RHS_sum + RHS(i, j) ** 2
                    !phi_pre = phi(i, j) ! 更新前のphiの値を格納
                    E(i, j) = BYM*phi(i,j-1) + BXM*phi(i-1,j) - B0*phi(i, j)  + BXP*phi(i+1,j) + BYP*phi(i,j+1) - RHS(i, j) ! E : ∇φ-(RHS)
                    phi(i, j) = phi(i, j) + omg * E(i, j) / B0
                    !er = er + (phi(i, j) - phi_pre) ** 2
                    er = er + (E(i, j)) ** 2

                    !write(*, *) i, j, RHS, E(i, j)
                enddo
            enddo
            !er = sqrt(er) / abs(RHS)
            er = sqrt(er) / RHS_norm ! er : 相対誤差
            
            !write(*, *) 'itr = ', itr, 'err = ', er
            if (er < er0) then
                 write(*, *) "exit  at itr = " , itr , 'Nx = ', Nx-1
                exit
            endif
            
        enddo

    end subroutine SOR

end module subprogs

program main
    use global
    use subprogs
    implicit none
    real(8) :: t
    real(8) :: error, diff
    real(8), allocatable :: x(:, :, :)
    real(8) :: BXP, BXM, BYP, BYM, B0
    real(8), allocatable :: u_pre(:, :)
    real(8), allocatable :: phi(:, :), phi_exact(:, :), phi_num(:, :)
    integer :: exp, Nx, Ny, N, i, j, k

    open(50, file = '/home/shishida/SMAC/.d/err.d') ! 出力ファイルを開く


    do exp = 3, 7 ! Nx, Ny = 8 ~ 128
        Nx = 2 ** exp + 1
        Ny = 2 ** exp + 1
        N = (Nx-1) * (Ny-1)

        if (allocated(x)) deallocate(x)
        if (allocated(u_pre)) deallocate(u_pre)
        if (allocated(phi)) deallocate(phi)
        if (allocated(phi_exact)) deallocate(phi_exact)
        if (allocated(phi_num)) deallocate(phi_num)

        allocate(x(2, Nx+1, Ny+1))
        allocate(u_pre(Nx+1, Ny+1))
        allocate(phi(Nx+1, Ny+1))
        allocate(phi_exact(Nx+1, Ny+1))
        allocate(phi_num(Nx+1, Ny+1))

        u_pre(:, :) = 0.0d0
        phi(:, :) = 0.0d0
        phi_exact(:, :) = 0.0d0
        phi_num(:, :) = 0.0d0

        dx = Lx / real(Nx-1,kind=8)
        dy = Ly / real(Ny-1,kind=8)

        ! ディリクレ境界条件
        phi(:, :) = 0.0d0
        phi(:, :) = 0.0d0

        ! 格子生成
        do j = 1, Ny
            do i = 1, Nx
                x(1, i, j) = (i-1) * dx
                x(2, i, j) = (j-1) * dy
            enddo
        enddo


        error = sqrt(error / real(Nx * Ny, kind=8))


        ! 圧力ポアソン方程式　定数項
        BXP = 1.0d0 / (dx ** 2)
        BXM = 1.0d0 / (dx ** 2)
        BYP = 1.0d0 / (dy ** 2)
        BYM = 1.0d0 / (dy ** 2)
        B0 = 2.0d0 / (dx ** 2) + 2.0d0 / (dy ** 2)

        !call gauss_seidel(BXP, BXM, BYP, BYM, B0, u_pre, phi, x)
        call SOR(BXP, BXM, BYP, BYM, B0, u_pre, phi, x, Nx, Ny)

        ! 二条平均平方根誤差を計算
        error = 0.0d0
        do j = 1, Ny
            do i = 1, Nx
                phi_exact(i, j) = sin(acos(-1.0d0)*x(1, i, j) / Lx) * sin(acos(-1.0d0)*x(2, i, j) / Ly)
                diff = phi(i, j) - phi_exact(i, j)
                error = error + diff ** 2
            enddo
        enddo
        error = sqrt(error / real(N, kind = 8))
        write(50, '(I8, 1X, E15.8)') Nx-1, error
    enddo

    close(50)



end program main