! correl.f90
! Compilation (requires BLAS & LAPACK): gfortran correl.f90 -O3 -lblas -llapack -o correl
! Use: ./correl N M
! Chain length = N
! Chain filling = M
! Creates/overwrites 2 files:
! - correl_data/correl.dat: symmetric C(N, N) matrix (upper triangular part)
! - correl_data/eigvals.dat: nu(N, N) eigenvalues of truncated C, (eigenvector, L) indices
program correl 
    use iso_fortran_env, only: RP=>real64
    implicit none
    integer, parameter :: MAXLEN = 100000
    integer :: N, M
    real(RP), allocatable :: Phi(:, :), C(:, :), nu(:, :), CL(:, :), eig(:)
    real(RP), allocatable :: WORK(:)
    integer :: L, INFO

    ! Initialisation
    call read_data(N, M, Phi)
    allocate(C(N, N), nu(N, N))

    ! Correlation matrix
    call dsyrk('U', 'N', N, M, 1._RP, Phi, N, 0._RP, C, N)

    ! Eigenvalues
    nu = 0._RP
    do L = 1, N
        write(*,'(A,A,I0,A,I0)', advance='no') char(13), 'Diagonalising L = ', L, ' / ', N
        allocate(CL(L, L), eig(L), WORK(3*L-1))
        CL = C(1:L, 1:L)
        call dsyev('N', 'U', L, CL, L, eig, WORK, 3*L-1, INFO)
        if (INFO /= 0) stop "Failed to diagonalise"
        nu(1:L, L) = eig
        deallocate(CL, eig, WORK)
    end do

    ! Save data
    call write_arrays(N, C, nu)

    ! Cleanup
    deallocate(Phi, C, nu)


    contains


    subroutine read_data(N, M, Phi)
        integer, intent(out) :: N, M
        real(RP), allocatable, intent(out) :: Phi(:, :)
        character(len=256) :: arg
        character(len=MAXLEN) :: line
        integer :: i, j, ios

        ! Get input file and N from command line
        call get_command_argument(1, arg)
        read(arg, *) N
        call get_command_argument(2, arg)
        read(arg, *) M

        allocate(Phi(N, M))

        ! Read Phi matrix from file
        open(unit=10, file='spin_chain_data/eigenvectors.dat', status='old')
        do i = 1, N
            read(10, *, iostat=ios) (Phi(i, j), j = 1, M)
            if (ios /= 0) then
                print *, 'Error reading row', i
                stop
            end if
        end do

        close(10)

    end subroutine read_data


    subroutine write_arrays(N, C, nu)
        integer, intent(in) :: N
        real(RP), dimension(N, N), intent(in) :: C, nu
        integer :: i, j, ios

        ! Open files
        open(unit=20, file='correl_data/correl.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) stop "Error opening correl_data/poly_data.dat file"

        open(unit=30, file='correl_data/eigvals.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) stop "Error opening correl_data/eigvals.dat file"

        ! Write correlation matrix, correlation eigenvalue matrix
        do i = 1, N
            write(20, '(*(G0.6,:," "))') (C(i,j), j=1,N)
            write(30, '(*(G0.6,:," "))') (nu(i,j), j=1,N)
        end do

        ! Close files
        close(20)
        close(30)

    end subroutine write_arrays


end program correl