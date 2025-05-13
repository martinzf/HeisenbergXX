! spin_chain.f90
! Compilation (requires LAPACK): gfortran spin_chain.f90 -O3 -llapack -o spin_chain
! Use: ./spin_chain [input.dat] N
! input.dat file:
!   0.0 1.0 2.0 ...   J(1) to J(N-1)
!   1.0 1.0 1.0 ...   b(1) to b(N)
! Creates/overwrites 3 files:
! - spin_chain_data/poly_data.dat: E(N), w(N), g(N) rows
! - spin_chain_data/polynomials.dat: P(N+1, N+1) polynomial coefficient matrix, (power, poly index) indices
! - spin_chain_data/eigenvectors.dat: Phi(N, N) eigenvector matrix, (position, energy) indices
program orth_polys
    use iso_fortran_env, only: RP=>real64
    implicit none
    integer, parameter :: MAXLEN = 100000
    integer :: N
    real(RP), allocatable :: a(:), J(:), b(:)
    real(RP), allocatable :: P(:, :), Phi(:, :), E(:), w(:), log_g(:)
    real(RP), allocatable :: dPN1(:), log_norm(:)
    integer :: i
    integer :: MM, INFO
    integer, allocatable :: ISUPPZ(:), IWORK(:)
    real(RP), allocatable :: WORK(:)

    ! Initialisation
    call read_data(N, J, b)
    allocate(a(N))
    a = J**2
    allocate(P(N+1, N+1), Phi(N, N), E(N), w(N), log_g(N))
    allocate(dPN1(N+1), log_norm(N))
    allocate(ISUPPZ(2*N), WORK(18*N), IWORK(10*N))

    ! Recursive definition: P(n, i) = n-1th power coef. of ith polynomial
    ! P_{n+1} = (x - b_n) P_n - a_{n-1} P_{n-1}
    P = 0._RP
    P(1, 1) = 1._RP
    P(1:2, 2) = (/-b(1), 1._RP/)
    do i = 2, N
        P(:, i+1) = cshift(P(:, i), -1) - b(i) * P(:, i) - a(i-1) * P(:, i-1)
    end do

    ! Normalisation coefficients
    log_g(1) = 0._RP
    do i = 2, N 
        log_g(i) = log(a(i-1)) + log_g(i-1)
    end do

    ! Roots
    call dstegr('V', 'A', N, b, J, 0._RP, 0._RP, 0, 0, 0.d0, MM, E, Phi, N, ISUPPZ, WORK, 18*N, IWORK, 10*N, INFO)
    if (INFO /= 0) stop "Error diagonalising matrix"

    ! Weights
    dPN1(N+1) = 0._RP
    dPN1(1:N) = [(i * P(i+1,N+1), i = 1, N)]
    w = exp(log_g(N)) / polyval(P(:, N), E) / polyval(dPN1, E)

    ! Save data
    call write_arrays(N, E, w, log_g, P, Phi)
    
    ! Cleanup
    deallocate(a, J, b)
    deallocate(P, Phi, E, w, log_g)
    deallocate(dPN1, log_norm)
    deallocate(ISUPPZ, WORK, IWORK)


    contains


    subroutine read_data(N, J, b)
        integer, intent(out) :: N
        real(RP), allocatable, intent(out) :: J(:), b(:)
        character(len=256) :: input_file, arg
        character(len=MAXLEN) :: line
        integer :: i, ios

        ! Get input file and N from command line
        call get_command_argument(1, input_file)
        call get_command_argument(2, arg)
        read(arg, *) N

        allocate(J(N), b(N))

        open(unit=10, file=trim(input_file), status='old')

        read(10, '(A)', iostat=ios) line
        if (ios /= 0) stop "Error reading first line"
        read(line, *) (J(i), i = 1, N-1)
        J(N) = 0.

        read(10, '(A)', iostat=ios) line
        if (ios /= 0) stop "Error reading second line"
        read(line, *) (b(i), i = 1, N)

        close(10)

    end subroutine read_data


    function polyval(coeffs, x) result(P)
        real(RP), intent(in) :: coeffs(:), x(:)
        real(RP) :: P(size(x))
        integer :: i

        P = 0.
        do i = 1, size(coeffs)
            P = P + coeffs(i) * x**(i-1)
        end do

    end function polyval


    subroutine write_arrays(N, E, w, log_g, P, Phi)
        integer, intent(in) :: N
        real(RP), dimension(N+1, N+1), intent(in) :: P
        real(RP), intent(in) :: Phi(N, N)
        real(RP), dimension(N), intent(in) :: E, w, log_g
        integer :: i, j, ios

        ! Open files
        open(unit=20, file='spin_chain_data/poly_data.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) stop "Error opening spin_chain_data/poly_data.dat file"

        open(unit=30, file='spin_chain_data/polynomials.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) stop "Error opening spin_chain_data/polynomials.dat file"

        open(unit=40, file='spin_chain_data/eigenvectors.dat', status='replace', action='write', iostat=ios)
        if (ios /= 0) stop "Error opening spin_chain_data/eigenvectors.dat file"

        ! Write three arrays as rows
        write(20, '(*(G0.6,:," "))') (- E(i), i = 1, N)
        write(20, '(*(G0.6,:," "))') (w(i), i = 1, N)
        write(20, '(*(G0.6,:," "))') (exp(log_g(i)), i = 1, N)

        ! Write polynomial, eigenvector and correlation matrices
        do i = 1, N
            write(30, '(*(G0.6,:," "))') (P(i,j), j=1,N+1)
            write(40, '(*(G0.6,:," "))') (Phi(i,j), j=1,N)
        end do
        write(30, '(*(G0.6,:," "))') (P(N+1,j), j=1,N+1)

        ! Close files
        close(20)
        close(30)
        close(40)

    end subroutine write_arrays


end program orth_polys