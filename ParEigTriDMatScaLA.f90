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
    implicit none

    integer, parameter :: SIZE = 50               ! matrix size
    logical, parameter :: TEST = .true.           ! Run with test values or load from file input

    integer            :: ierror                  ! error code
    integer            :: my_rank                 ! rank of process
    integer            :: num_nodes               ! number of processes
    integer            :: source                  ! rank of sender
    integer            :: dest                    ! rank of receiver
    integer            :: status(MPI_STATUS_SIZE) ! return status for receive
    integer            :: tag                     ! tag for messages

    integer            :: i, j, k, ll, m, isum, ioffset, cnt, N, N1, N2
    real               :: norm

    real                    :: bn
    real, allocatable       :: a(:)          ! diagonal of TriDiag matrix
    real, allocatable       :: b(:)          ! off diagonal of TriDiag matrix
    real, allocatable       :: Q(:,:)          ! Q as determined from Q1 and Q2
    real, allocatable       :: Q1(:,:)         ! first block orthogonal matrix
    real, allocatable       :: Q2(:,:)         ! second block orthogonal matrix
    real, allocatable       :: lambda(:)     ! will store eigenvectors
    real, allocatable       :: xi(:)
    real, allocatable       :: d(:)          ! will store computed eigenvalues
    real, allocatable       :: tmpd(:)
    real, allocatable       :: tmpdd(:,:)

    integer, allocatable    :: indices(:,:)
    real, allocatable       :: adjust(:,:)

    dest = 0

    if (TEST) then
        N = SIZE
        allocate(a(N))
        allocate(b(N)) ! treat b(end) as 0 but instantiate full vector for convenience
        do i = 1, N
            a(i) = i+0.5;
            b(i) = 0.3;
        end do
    else
        ! TODO: Load vector data files
        N = 50
    endif
    allocate(Q(N,N))
    allocate(Q1(N,N))
    allocate(Q2(N,N))
    allocate(xi(N))
    allocate(d(N))
    allocate(tmpd(N))
    allocate(tmpdd(N,N))

    !TODO: Initialize and implement Scalapack routines

    ! start up MPI
    call MPI_Init(ierror)

    ! find process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)

    ! find number of processes and allocate indexing and adjustment matrices
    call MPI_Comm_size(MPI_COMM_WORLD, num_nodes, ierror)
    allocate(indices(num_nodes,num_nodes))
    allocate(adjust(num_nodes,num_nodes))

    ! Determine indices for breaking up global problem into sub-problem tree
    indices(1,1) = N
    ioffset = ceiling(log(real(num_nodes)) / log(2.))
    do i = 1, ioffset
        isum = 1
        do j = 1, 2**i
            indices(i+1,2*j-1) = indices(i,j) / 2
            indices(i+1,2*j) = indices(i,j) - indices(i+1,2*j-1)
            isum = isum + indices(i+1,2*j-1)
            adjust(i,j) = b(isum-1)
            a(isum-1) = a(isum-1) - b(isum-1)
            a(isum) = a(isum) - b(isum-1)
            isum = isum + indices(i+1,2*j)
        end do
    end do

    isum = 1
    k = 0
    do while (k .LT. my_rank)
        isum = isum + indices(ioffset, k+1)
        k = k + 1
    end do
    ! TODO: configure ScaLAPACK function here
    call TDQREigensolver(indices(ioffset,my_rank+1), a(isum), b(isum),d ,Q1)

    ! Fan-in algorithm to finish solving system
    ! TODO: Check for OpenMP !omp paralellization
    tag = 1
    do i = 1, ioffset
        isum = 1
        cnt = 1
        do j = 1, num_nodes, 2**(i-1)
            if (my_rank+1 == j) then
                if (MOD(isum,2) == 1)  then
                    source = j+2**(i-1)
                    call MPI_Recv(d+indices(ioffset,isum), indices(ioffset,isum+1), &
                        MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierror)
                    do k = 1, indices(ioffset,isum+1)
                        call MPI_Recv(Q2(k,:), indices(ioffset,isum+1), &
                            MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierror)
                        ! TODO: double check MPI Recv/Send of submatrix vector data
                    end do
                    N1 = indices(ioffset,isum)
                    N2 = indices(ioffset,isum+1)
                    bn = adjust(ioffset-1,cnt)
                    cnt = cnt + 1
                else
                    dest = j-2**(i-1)
                    call MPI_Send(d, indices(ioffset,isum), &
                        MPI_REAL, dest, tag, MPI_COMM_WORLD, ierror)
                    do k = 1, indices(ioffset,isum)
                        call MPI_Send(Q1(k,:), indices(ioffset,isum), &
                            MPI_REAL, dest, tag, MPI_COMM_WORLD, ierror)
                    end do
                end if
            end if
            isum = isum + 1
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
                ! obtain eigenvalues
                ! TODO: Bisection method?
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
        ioffset = ioffset - 1
    end do

!    if (mynode == 0) then
!        cout << "The eigenvalues are: " << endl;
!        for(k=0; k<size; k++)
!            cout << d[k] << endl;
!    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    ! TODO: deallocate memory
    call MPI_Finalize(ierror)
    stop
end program
