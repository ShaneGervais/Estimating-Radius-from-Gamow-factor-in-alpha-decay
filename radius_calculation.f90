program gamow_radius

  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  ! physical constants
  real(dp), parameter :: hbarc = 197.3269804_dp ! Mev*fm [reduced plank * speed of light]
  real(dp), parameter :: e2 = 1.44_dp ! MeV*fm [e^2 / 4 pi epsi_0]
  real(dp), parameter :: amu = 931.494_dp ! MeV/c^2 [amu to MeV]
  real(dp), parameter :: me = 0.511_dp ! MeV/c^2 [mass of electron]
  real(dp), parameter :: r0 = 1.2_dp ! fm [Wood-Saxon approx -> R = r0 A^1/3]

  ! inputs
  real(dp), parameter :: tau = 7.3e10_dp ! s lifetime of 226Ra
  real(dp), parameter :: tau0 = 6.3e-23_dp ! s from exp
  integer(dp), parameter :: Z_D = 86
  integer(dp), parameter :: A_D = 222
  integer(dp), parameter :: A_alpha = 4
  integer(dp), parameter :: Z_alpha = 2
  integer(dp), parameter :: Z_P = 88
  integer(dp), parameter :: A_P = 226
  
  ! atomic masses
  real(dp), parameter :: m_alpha = 4.002603_dp*amu ! - Z_alpha*me
  real(dp), parameter :: m_D     = 222.0175_dp*amu ! - Z_D*me
  real(dp), parameter :: m_P     = 226.0254_dp*amu ! - Z_P*me
  
  ! Netwon's method settings
  real(dp), parameter :: tol = 1.0e-12_dp
  integer, parameter :: max_iter = 1000

  real(dp) :: Q, m_red, G, pref, C, b, x, x_new, g_val, a, R_alpha, R_D, r0_fit
  real(dp) :: pi
  integer(dp) :: n
  logical :: converged

  Q = m_P - m_D - m_alpha
  m_red = m_alpha * m_D / (m_D + m_alpha)

  pi = acos(-1.0_dp)

  ! Gamow factor from lifetime
  G = log(tau/tau0)

  ! prefactor and target value -> C
  pref = 4.0_dp*real(Z_D, dp)*(e2/hbarc)*sqrt(2.0_dp*m_red/Q)
  C = G/pref

  ! coulomb -> classical
  b = 2.0_dp*real(Z_D, dp)*e2/Q

  ! AI-gen formatting
  write(*, '(a)')
  write(*, '(a, f12.4,a)')
  write(*, '(a,f12.4)')
  write(*, '(a,f12.4)')
  write(*, '(a,f12.6)')
  write(*,'(a,f12.4,a)')

  ! thick barrier starting guess
  x = ((pi/2.0_dp - C)/2.0_dp)**2

  ! itterate
  write(*,'(/,a)') ' -- Netwon iteractions -- '
  write(*,'(a)')    'n        x_n        g(x_n)'
  converged = .false.

  do n = 0, max_iter
    g_val = f(x) - C
    write(*, '(i4,2es20.10)') n, x, g_val

    x_new = x - g_val/fprime(x) ! see derivative function below
    
    ! boundaries 0 - 1
    if (x_new <= 0.0_dp) x_new = 0.5_dp*x
    if (x_new >= 1.0_dp) x_new = 0.5_dp*(x + 1.0_dp)

    if (abs(x_new - x) < tol) then
      x = x_new
      converged = .true.
      exit
    end if
    x = x_new
  end do

  if (.not. converged) then
    write(*, '(a)') 'WARNING: Newton did not converge'
  end if

  ! estimate radius
  a = x*b
  R_alpha = r0*real(A_alpha, dp)**(1.0_dp/3.0_dp)
  R_D = a - R_alpha
  r0_fit = a/real(A_alpha, dp)**(1.0_dp/3.0_dp) + real(A_D, dp)**(1.0_dp/3.0_dp)

  ! display results
  write(*, '(/,a)')       '-- Results --'
  write(*,'(a,f14.10)')   ' x = a/b                = ', x
  write(*,'(a,es14.4)')   ' residual f(x) - C      = ', f(x) - C
  write(*,'(a,f10.4,a)')  ' touching radius a      = ', a, ' fm'
  write(*,'(a,f10.4,a)')  ' R_alpha (r0 = 1.2 fm)  = ', R_alpha, ' fm'
  write(*,'(a,f10.4,a)')  ' R(222Rn) = a - R_alpha = ', R_D, ' fm'
  write(*,'(a,f10.4,a)')  ' fitted r0 from a       = ', r0_fit, ' fm'
  write(*,'(a,f10.4,a)')  ' compare 1.2*222^(1/3)  = ', r0*real(A_D, dp)**(1.0_dp/3.0_dp), ' fm'

  contains

  pure function f(xx) result(val)
    real(dp), intent(in) :: xx
    real(dp) :: val
    val = acos(sqrt(xx)) - sqrt(xx*(1.0_dp - xx))
  end function f

  pure function fprime(xx) result(val)
    real(dp), intent(in) :: xx
    real(dp) :: val

    val = -sqrt((1.0_dp - xx)/xx)
  end function fprime

end program gamow_radius
