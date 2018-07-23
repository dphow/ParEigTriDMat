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
    integer            :: tag                     ! tag for messages

    integer            :: i, j, k, ll, m, iStart, iEnd, treeDepth, cnt, N, N1, N2
    real               :: norm

    real                    :: bn
    real, allocatable       :: a(:)          ! diagonal of TriDiag matrix
    real, allocatable       :: b(:)          ! off diagonal of TriDiag matrix
    real, allocatable       :: Q(:,:)        ! Q as determined from Q1 and Q2
    real, allocatable       :: Q1(:,:)       ! first block orthogonal matrix
    real, allocatable       :: Q2(:,:)       ! second block orthogonal matrix
    real, allocatable       :: lambda(:)     ! will store eigenvalues
    real, allocatable       :: xi(:)
    real, allocatable       :: d(:)          ! will store computed eigenvalues
    real, allocatable       :: tmpd(:)
    real, allocatable       :: tmpdd(:,:)

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
        N = 64
    endif
    allocate(Q(N,N))
    allocate(Q1(N,N))
    allocate(Q2(N,N))
    allocate(xi(N))
    allocate(d(N))
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

    ! Determine problem sizes for breaking up global problem into sub-problem tree
    ! Record b_k*zz' adjustments for use later
    allocate(probSizes(treeDepth,num_nodes))
    allocate(adjust(treeDepth,num_nodes))
    probSizes = 0
    adjust = 0

    probSizes(1,1) = N
    do i = 1, treeDepth-1
        iStart = 1
        do j = 1, 2**(i-1)
            probSizes(i+1,2*j-1) = probSizes(i,j) / 2
            probSizes(i+1,2*j) = probSizes(i,j) - probSizes(i+1,2*j-1)
            iStart = iStart + probSizes(i+1,2*j-1)
            adjust(i,j) = b(iStart-1)
            a(iStart-1) = a(iStart-1) - b(iStart-1)
            a(iStart) = a(iStart) - b(iStart-1)
            iStart = iStart + probSizes(i+1,2*j)
        end do
    end do

    !if (my_rank == 0) then
    !    print *, "Problem Size Tree"
    !    do i = 1, treeDepth
    !        print *, probSizes(i,:)
    !    end do
    !end if

    ! Compute start and end probSizes for every process' subproblem
    iStart = 1
    k = 1
    do while (k <= my_rank)
        iStart = iStart + probSizes(treeDepth, k)
        k = k + 1
    end do
    iEnd = iStart + probSizes(treeDepth, my_rank+1) - 1

    ! Debug
    !print *, my_rank, "pSize", probSizes(treeDepth,my_rank+1)
    !print *, my_rank, "a", a(iStart:iEnd)
    !print *, my_rank, "b", b(iStart:iEnd)

    ! Solve Eigenproblem at bottom-most depth of problem tree
    call TDQREigensolver(probSizes(treeDepth,my_rank+1), a(iStart:iEnd), b(iStart:iEnd), d ,Q1, my_rank)
    print *, my_rank, d
    call exit(-1)
    ! Fan-in algorithm to finish solving system for each level up on tree
    tag = 1
    do i = 1, treeDepth-1
        iStart = 1
        cnt = 1
        do j = 1, num_nodes, 2**(i-1)
            if (my_rank+1 == j) then
                if (MOD(iStart,2) == 1)  then
                    source = j+2**(i-1)
                    call MPI_Recv(d+probSizes(treeDepth,iStart), probSizes(treeDepth,iStart+1), &
                        MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierror)
                    do k = 1, probSizes(treeDepth,iStart+1)
                        call MPI_Recv(Q2(k,:), probSizes(treeDepth,iStart+1), &
                            MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierror)
                        ! TODO: double check MPI Recv/Send of submatrix vector data
                    end do
                    N1 = probSizes(treeDepth,iStart)
                    N2 = probSizes(treeDepth,iStart+1)
                    bn = adjust(treeDepth-1,cnt)
                    cnt = cnt + 1
                else
                    dest = j-2**(i-1)
                    call MPI_Send(d, probSizes(treeDepth,iStart), &
                        MPI_REAL, dest, tag, MPI_COMM_WORLD, ierror)
                    do k = 1, probSizes(treeDepth,iStart)
                        call MPI_Send(Q1(k,:), probSizes(treeDepth,iStart), &
                            MPI_REAL, dest, tag, MPI_COMM_WORLD, ierror)
                    end do
                end if
            end if
            iStart = iStart + 1
        end do

        do j = 1, num_nodes, 2**i
            if (my_rank+1 == j) then
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
                call SolveSecularEq(bn,N1+N2,d,xi,lambda);

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
        treeDepth = treeDepth - 1
    end do

    if (my_rank == 0) then
        print *, "The eigenvalues are: "
        print *, d
    end if

!    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    deallocate(a, b, Q, Q1, Q2, xi, d, tmpd, tmpdd)
    call MPI_Finalize(ierror)
    stop
end program
