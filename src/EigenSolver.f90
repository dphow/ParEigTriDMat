module EigenSolver
    implicit none

contains
    recursive subroutine TDQREigensolver(N, a, b, lambda, Q, my_rank)
        use SecularSolver
        implicit none
        integer, intent(in)         :: N, my_rank
        real (KIND (0.0D0)), intent(inout)            :: a(:), b(:)
        real (KIND (0.0D0)), intent(out)           :: lambda(:)
        real (KIND (0.0D0)), intent(inout)         :: Q(:,:)

        integer                     :: i, j, k, N1, N2, cnt
        real (KIND (0.0D0))                        :: d(N), total
        real (KIND (0.0D0)), allocatable           :: xi(:), Q1(:,:), Q2(:,:)

        if (N == 1) then
            Q(1,1) = 1.0
            lambda(1) = a(1)
        else
            N1 = N / 2
            N2 = N - N1
            allocate(Q1(N1,N1))
            allocate(Q2(N2,N2))
            Q1 = 0.0
            Q2 = 0.0

            a(N1) = a(N1) - b(N1)
            a(N1+1)   = a(N1+1) - b(N1)

            ! Recursively solve each subproblem
            call TDQREigensolver(N1, a(1:N1), b(1:N1), d(1:N1), Q1, my_rank)
            call TDQREigensolver(N2, a(N1+1:N), b(N1+1:N), d(N1+1:N), Q2, my_rank)

            ! Store total eigenvector from rank 1 modification
            cnt = 1
            allocate(xi(N))
            do i = 1, N1
                xi(cnt) = Q1(N1,i)
                cnt = cnt + 1
            end do
            do i = 1, N2
                xi(cnt) = Q2(1,i)
                cnt = cnt + 1
            end do

            call SolveSecularEq(b(N1),N,d,xi,lambda, my_rank)

            ! Compute the eigenvectors
            do i = 1, N1
                do j = 1, N
                    Q(i,j) = 0.0
                    do k = 1, N1
                        Q(i,j) = Q(i,j) + Q1(i,k)*xi(k)/(d(k)-lambda(j))
                    end do
                end do
            end do
            do i = 1, N2
                do j = 1, N
                    Q(N1+i,j) = 0.0
                    do k = 1, N2
                        Q(N1+i,j) = Q(N1+i,j) + Q2(i,k)*xi(N1+k)/(d(N1+k)-lambda(j))
                    end do
                end do
            end do

            ! Normalize each eigenvector in Q
            do j = 1, N
                total = 0.0
                do i = 1, N
                    total = total + Q(i,j)**2
                end do
                total = sqrt(total)
                do i = 1, N
                    Q(i,j) = Q(i,j)/total
                end do
            end do

            deallocate(Q1, Q2, xi)
        end if
    end
end module EigenSolver
