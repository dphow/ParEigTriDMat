! ============================================================================
! Name        : ParEigTriDMatScaLA.f90
! Author      : Daniel Howard
! Version     : 0.1 beta
! Copyright   : GBU General Public License 2018
! Description : Calculate Eigenvalues and optionally Eigenvectors of a symmetric tridiagonal matrix
! ============================================================================
! Usage Description
! TODO

program ParEigTriDMatScaLA
    use mpi
    use EigenSolver
    use SecularSolver
    implicit none

    integer, parameter :: DIMENSIONS = 32       ! matrix size
    logical, parameter :: TEST = .true.           ! Run with test values or load from file input

    integer            :: ierror                  ! error code
    integer            :: my_rank                 ! rank of process
    integer            :: num_nodes               ! number of processes
    integer            :: source                  ! rank of sender
    integer            :: dest                    ! rank of receiver
    integer            :: status(MPI_STATUS_SIZE) ! return status for receive

    integer            :: i, j, k, ll, m, iStart, iEnd, treeDepth, currentDepth, cnt, N, N1, N2

    real (KIND (0.0D0))                    :: norm
    real (KIND (0.0D0))                    :: bn
    real (KIND (0.0D0)), allocatable       :: a(:)          ! diagonal of TriDiag matrix
    real (KIND (0.0D0)), allocatable       :: b(:)          ! off diagonal of TriDiag matrix
    real (KIND (0.0D0)), allocatable       :: Q(:,:)        ! Q as determined from Q1 and Q2
    real (KIND (0.0D0)), allocatable       :: Q1(:,:)       ! first block orthogonal matrix
    real (KIND (0.0D0)), allocatable       :: Q2(:,:)       ! second block orthogonal matrix
    real (KIND (0.0D0)), allocatable       :: lambda(:)     ! will store eigenvalues
    real (KIND (0.0D0)), allocatable       :: xi(:)
    real (KIND (0.0D0)), allocatable       :: d(:)          ! will store computed eigenvalues
    real (KIND (0.0D0)), allocatable       :: tmpd(:)
    real (KIND (0.0D0)), allocatable       :: tmpdd(:,:)

    integer, allocatable    :: probSizes(:,:)
    real, allocatable       :: adjust(:,:)

    dest = 0

    if (TEST) then
        N = DIMENSIONS
        allocate(a(N))
        allocate(b(N)) ! treat b(end) as 0 but instantiate full vector for convenience
        do i = 1, N
            a(i) = i-0.50d0
            b(i) = 0.30d0
        end do
    else
        ! TODO: Load from command line arguments and/or local files
        N = 64
    endif
    allocate(Q(N,N))
    allocate(Q1(N,N))
    allocate(Q2(N,N))
    allocate(xi(N))
    allocate(d(N))
    allocate(lambda(N))
    allocate(tmpd(N))
    allocate(tmpdd(N,N))

    ! TODO: Initialize and implement Scalapack routines

    ! start up MPI
    call MPI_Init(ierror)

    ! find process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)

    ! find number of processes and allocate indexing and adjustment matrices
    call MPI_Comm_size(MPI_COMM_WORLD, num_nodes, ierror)

    ! Create problem tree of depth relative to number of compute procresses available
    treeDepth = floor(log(real(num_nodes)) / log(2.)) + 1
    currentDepth = treeDepth

    ! Determine problem sizes for breaking up global problem into sub-problem tree
    ! Record b_k*zz' adjustments for use later
    allocate(probSizes(treeDepth,num_nodes))
    allocate(adjust(treeDepth,num_nodes))
    probSizes = 0
    adjust = 0

    probSizes(1,1) = N
    do i = 2, treeDepth
        iStart = 1
        do j = 1, 2**(i-2)
            probSizes(i,2*j-1) = probSizes(i-1,j) / 2
            probSizes(i,2*j) = probSizes(i-1,j) - probSizes(i,2*j-1)
            iStart = iStart + probSizes(i,2*j-1)
            adjust(i-1,j) = b(iStart-1)
            a(iStart-1) = a(iStart-1) - b(iStart-1)
            a(iStart) = a(iStart) - b(iStart-1)
            iStart = iStart + probSizes(i,2*j)
        end do
    end do

    if (my_rank == 0) then
        print *, "Problem Size Tree"
        do i = 1, treeDepth
            print *, probSizes(i,:)
        end do
    end if

    ! Compute start and end probSizes for every process' subproblem
    iStart = 1
    k = 1
    do while (k <= my_rank)
        iStart = iStart + probSizes(treeDepth, k)
        k = k + 1
    end do
    iEnd = iStart + probSizes(treeDepth, my_rank+1) - 1

    ! Solve Eigenproblem at bottom-most depth of problem tree
    call TDQREigensolver(probSizes(treeDepth,my_rank+1), a(iStart:iEnd), b(iStart:iEnd), d, Q1, my_rank)
    if (my_rank == 0) then
    print *, "Process rank:", my_rank, " D:", d(1:8), "a", a
    end if
    !print *, "Process rank:", my_rank, " Q1:", Q1(1:8,1:8)
    ! Fan-in algorithm to finish solving system for each level up on tree
    do i = 0, treeDepth-2
        iStart = 1
        cnt = 1
        do j = 0, num_nodes-1, 2**i
            if (my_rank == j) then
                N1 = probSizes(currentDepth,iStart)
                N2 = probSizes(currentDepth,iStart+1)
                if (MOD(iStart,2) == 1)  then
                    source = j+2**i
                    call MPI_Recv(d(N1+1:N1+N2), N2, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, status, ierror)
                    do k = 1, probSizes(currentDepth,iStart+1)
                        call MPI_Recv(Q2(k,1:N2), N2, MPI_DOUBLE_PRECISION, source, k, MPI_COMM_WORLD, status, ierror)
                    end do
                    bn = adjust(treeDepth-1,cnt)
                    cnt = cnt + 1
                else
                    dest = j-2**i
                    call MPI_Send(d(1:N1), N1, MPI_DOUBLE_PRECISION, dest, 0, MPI_COMM_WORLD, ierror)
                    do k = 1, N1
                        call MPI_Send(Q1(k,1:N1), N1, MPI_DOUBLE_PRECISION, dest, k, MPI_COMM_WORLD, ierror)
                    end do
                end if
            end if
            iStart = iStart + 1
        end do

        do j = 0, num_nodes-1, 2**(i+1)
            if (my_rank == j) then
                cnt = 1
                do k = 1, N1
                    xi(cnt) = Q1(N1,k)
                    cnt = cnt + 1
                end do
                do k = 1, N2
                    xi(cnt) = Q2(1,k)
                    cnt = cnt + 1
                end do
                ! Solve for the secular equation to
                ! obtain eigenvalues by bisection method
                call SolveSecularEq(bn,N1+N2,d,xi,lambda, my_rank);
                ! Form the Q matrix from Q1 and Q2
                do k = 1, N1
                    do ll = 1, N1+N2
                        Q(k,ll) = 0.0
                        do m = 1, N1
                            Q(k,ll) = Q(k,ll) + (Q1(k,m)*xi(m) / (d(m)-lambda(ll)))
                        end do
                    end do
                end do
                do k = 1, N2
                    do ll = 1, N1+N2
                        Q(N1+k,ll) = 0.0;
                        do m = 1, N2
                            Q(k+N1,ll) = Q(k+N1,ll) + (Q2(k,m)*xi(N1+m) / (d(N1+m)-lambda(ll)))
                        end do
                    end do
                end do

                ! Normalize the Q matrix so that each eigenvector
                ! has length one
                do k = 1, N1+N2
                    norm = 0.0
                    do ll = 1, N1+N2
                        norm = norm + Q(ll,k)*Q(ll,k)
                    end do
                    norm = sqrt(norm)
                    do ll = 1, N1+N2
                        Q(ll,k) = Q(ll,k)/norm
                    end do
                end do
                ! Swap d and lambda arrays for use in the
                ! next part of the fan-in algorithm
                tmpd = d
                d = lambda
                lambda = tmpd

                ! Swap Q and Q1 for use in the
                ! next part of the fan-in algorithm
                tmpdd = Q1
                Q1 = Q
                Q = tmpdd
            end if
        end do
        currentDepth = currentDepth - 1
    end do

    if (my_rank == 0) then
        print *, "The eigenvalues are: "
        print *, d
    end if

    deallocate(a, b, Q, Q1, Q2, xi, lambda, d, tmpd, tmpdd)
    call MPI_Finalize(ierror)
    stop
end program
