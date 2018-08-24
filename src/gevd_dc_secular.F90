!> @brief    Compute the distances from the old eigenvalues (poles) to new eigenvlaues
!!           of the deflated system \bar{D} - \rho z z' + \lambda (I - z z').
!! @note     We check the distances from the nearest poles  for all new !eigenvalues.
!!
!! @param n                    Matrix dimension.
!! @param d                    d(:) must be stored in ascending order. dimension(n)
!! @param z                    dimension(n)
!! @param rho                  Scalar rho.
!! @param y                    Distance between new eigenvalues and their nearest poles. dimension(n).
!!                             New eigenvalue lambda(j) = d(nearest_pole_index(j)) - y(j).
!! @param nearest_pole_index   See the description of y. dimension(n).
!!                             y(n) := d(nearest_root_index(j)) - lambda_i
!! @param rwork                Workspace. dimension(6 * n).
!! @param info                 Return 0 if success.
subroutine gevd_dc_secular(n, d, z, rho, y, nearest_root_index, rwork, info)
  ! Compute eigenvalues of \bar{D} - \rho z z' + \lambda (I - z z')
  implicit none
  integer, intent(in) :: n
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho
  real(kind(0.0d0)), intent(out) :: y(*)
  integer, intent(out) :: nearest_root_index(*)
  real(kind(0.0d0)), intent(out) :: rwork(6 * n)
  integer, intent(out) :: info
  !
  integer :: i
  
!$omp parallel
!$omp do private(i, info, rwork) schedule(guided, 2)
  do i = 1, n
    call gevd_dc_secular_elsner(n, i, d, z, rho, y(i), nearest_root_index(i), rwork, info)
    if (info /= 0)  then
      ! write (*, *) '[Warning] gevd_dc_secular_elsner failed at i = ', i, '. Trying bisection.'
      call gevd_dc_secular_bisection(n, i, d, z, rho, y(i), nearest_root_index(i), info)
      if (info /= 0)  then
        write (*, *) '[Error] gevd_dc_secular failed at i = ', i
        stop
      endif
    endif
  enddo
!$omp end parallel

  info = 0
  return
end

