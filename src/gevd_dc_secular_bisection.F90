!> @brief    Evaluate a secular equation f(x) accurately  where
!!           f(x) := \\sum_{k = 1}^{n} \\frac{z_k^2 (d_k - \\rho)}{x - d_k} - 1 + \\bm{z}^\top \\bm{z}.
!! @note    x is implicitly given as y := x - d(i).
!! @note    This subroutine does not check argument validity.
!!
!! @param n        The dimension of the vectors d and z.
!! @param d        The vector d. dimension(n).
!! @param z        The vector z. dimension(n).
!! @param rho      The scalar \\rho.
!! @param y        The scalar argument of the function. x \\ne d_i (i = 1, ..., n) is required.
!! @param i        The nearest poll from the eigenvalue to compute.
!! @param f        The function value f(x).
!! @param info     Error code. This never fail (always returns info = 0).
subroutine gevd_dc_secular_evaluation(n, d, z, rho, y, i, f, info)
  implicit none
  integer, intent(in) :: n
  real(kind(0.0d0)), intent(in) :: d(*), z(*), rho, y
  integer, intent(in) :: i
  real(kind(0.0d0)), intent(out) :: f
  integer, intent(out) :: info
  !
  integer :: k, m
  real(kind(0.0d0)) :: plus, minus

  !> Find m such that d(m) < rho < d(m + 1).
  if (rho < d(1)) then
    m = 0
  else
    m = n
    do k = 1, n - 1
      if (rho < d(k + 1)) then
        m = k
        exit
      endif
    enddo
  endif
  
  !> The main computation (include summation loop).
  !! The positive and negative terms are separately accumulated and finally compute the sum of them
  plus = 0.0d0
  minus = 0.0d0
  if (i .le. m) then
    ! The i-th smallest eigenvalue is at the left side of rho.
    do k = n, m + 1, -1
      minus = minus + (z(k) ** 2.0d0) * (d(k) - rho) / ((d(i) - d(k)) + y)
    enddo
    do k = 1, i - 1
      minus = minus + (z(k) ** 2.0d0) * (d(k) - rho) / ((d(i) - d(k)) + y)
    enddo
    do k = m, i, -1
      plus = plus + (z(k) ** 2.0d0) * (d(k) - rho) / ((d(i) - d(k)) + y)
    enddo
  else
    ! The i-th smallest eigenvalue is at the right side of rho.
    do k = 1, m
      minus = minus + (z(k) ** 2.0d0) * (d(k) - rho) / ((d(i) - d(k)) + y)
    enddo
    do k = m + 1, i
      plus = plus + (z(k) ** 2.0d0) * (d(k) - rho) / ((d(i) - d(k)) + y)
    enddo
    do k = n, i + 1, -1
      minus = minus + (z(k) ** 2.0d0) * (d(k) - rho) / ((d(i) - d(k)) + y)
    enddo
  endif
  !
  do k = 1, n
    plus = plus + z(k) ** 2.0d0
  enddo
  minus = minus - 1.0d0
  !
  f = plus + minus

  info = 0
  return
end

!> @brief    Compute i-th smallest eigenvalues of
!!           D - \\rho \\bm{z} \\bm{z}^\\top - \\lambda (I - \\bm{z} \\bm{z}^\\top)
!!           accurately by bisection method.
!!
!! @param n                    The dimension of the matrix.
!! @param i                    The eigenvalue number to compute.
!! @param d                    The diagonal entries of D (is diagonal). d(1) < ... < d(n). dimension(n).
!! @param z                    The vector \bm{z}. z_i /= 0 (i = 1, ..., n) is required.
!! @param rho                  The scalar \rho.
!! @param y                    Distance between new eigenvalue and from the newrest pole.
!!                             y := d(nearest_root_index) - \lambda (= new eigenvalue).
!! @param nearest_pole_index   See the description of y. dimension(n).
!!                             y(n) := d(nearest_root_index(j)) - lambda_i
!! @param info                 Error code. This subroutine always returns info = 0.
subroutine gevd_dc_secular_bisection(n, i, d, z, rho, y, nearest_root_index, info)
  implicit none
  integer, intent(in) :: n, i
  real(kind(0.0d0)), intent(in) :: d(n), z(n), rho
  real(kind(0.0d0)), intent(out) :: y 
  integer, intent(out) :: nearest_root_index
  integer, intent(out) :: info
  !
  integer :: k
  integer, parameter :: maxit = 81
  real(kind(0.0d0)) :: my ! my := - y
  real(kind(0.0d0)) :: lower, upper, one_minus_zzt, bound_diff, f
  logical :: lambda_is_left ! true when lambda < rho, else rho > lambda

  !> Find the nearest poll from the eigenvalue to compute.
  one_minus_zzt = 1.0d0
  do k = 1, n
    one_minus_zzt = one_minus_zzt - z(k) ** 2.0d0
  enddo
  bound_diff = 0.0d0
  do k = 1, n
    bound_diff = bound_diff + (z(k) ** 2.0d0) * dabs(d(k) - rho)
  enddo
  bound_diff = dble(n) * bound_diff / one_minus_zzt
  !
  if (d(i) < rho) then
    ! d(i - 1) < lambda < d(i) < rho
    if (i .ne. 1) then
      my = 0.5d0 * (d(i - 1) - d(i))
      call gevd_dc_secular_evaluation(n, d, z, rho, my, i, f, info)
      if (f .lt. 0.0d0) then
        nearest_root_index = i
      else
        nearest_root_index = i - 1
      endif
    else
      nearest_root_index = i
    endif
  else
    ! rho < d(i) < lambda < d(i + 1)
    if (i .ne. n) then
      my = 0.5d0 * (d(i + 1) - d(i))
      call gevd_dc_secular_evaluation(n, d, z, rho, my, i, f, info)
      if (f .lt. 0.0d0) then
        nearest_root_index = i
      else
        nearest_root_index = i + 1
      endif
    else
      nearest_root_index = i
    endif
  endif
  
  !> Compute the upper and lower bounds of the eigenvalue and
  !! iteratively compute the eigenvalue.
  if (nearest_root_index .eq. i) then
    if (d(i) < rho) then
      ! d(i - 1) < lambda < d(i) < rho
      lambda_is_left = .true.
      if (i .ne. 1) then
        lower = d(i - 1) - d(i)
      else
        lower = - bound_diff
      endif
      upper = 0.0d0
    else
      ! rho < d(i) < lambda < d(i + 1)
      lambda_is_left = .false.
      lower = 0.0d0
      if (i .ne. n) then
        upper = d(i + 1) - d(i)
      else
        upper = bound_diff
      endif
    endif
    do k = 1, maxit
      my = 0.5d0 * (lower + upper)
      call gevd_dc_secular_evaluation(n, d, z, rho, my, i, f, info)
      if ((f < 0.0d0 .and. lambda_is_left) .or. &
          (f > 0.0d0 .and. (.not. lambda_is_left))) then
        lower = my
      else
        upper = my
      endif
    enddo
    my = 0.5d0 * (lower + upper)
  else
    if (d(i) < rho) then
      ! d(i - 1) < lambda < d(i) < rho
      lambda_is_left = .true.
      lower = 0.0d0
      upper = d(i) - d(i - 1) ! i .ne. 1 assumed
    else
      ! rho < d(i) < lambda < d(i + 1)
      lambda_is_left = .false.
      lower = d(i) - d(i + 1)
      upper = 0.0d0
    endif
    do k = 1, maxit
      my = 0.5d0 * (lower + upper)
      if (lambda_is_left) then
        call gevd_dc_secular_evaluation(n, d, z, rho, my, i - 1, f, info)
      else
        call gevd_dc_secular_evaluation(n, d, z, rho, my, i + 1, f, info)
      endif
      if ((f < 0.0d0 .and. lambda_is_left) .or. &
          (f > 0.0d0 .and. (.not. lambda_is_left))) then
        lower = my
      else
        upper = my
      endif
    enddo
    my = 0.5d0 * (lower + upper)
  endif
  y = - my

  info = 0
  return
end

