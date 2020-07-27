program EOB
implicit none

real * 16 , parameter :: Pi=3.1415926535897932384626433832795028841971693993751_16
COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
real * 16 :: coordinates(3),Momentum(3),u,iota,p,e,theta_min,theta_max,r_min,r_max,R,theta,phi
real * 16 :: m_1,m_2,S_1,S_2,S_Kerr,alpha,beta,gamma,g(5)
real * 16 :: Derivative_of_parameters_(10)
real * 16 :: The_guy_with_alpha,H_hat,H_eff_over_mu,Partial_H_hat_partial_H_eff_over_mu,Energy
real * 16 :: Carter
integer :: i

real * 16 :: g_1(5),g_2(5),a1,a2,b1,b2,c1,c2


!——————————————————————以下参数为采用Runge-Kutta方法求解微分方程而定义————————————————————————————————

real * 16 :: K1(5),K2(5),K3(5),K4(5),h
real * 16 :: Temp_X(3),Temp_P(3)
real * 16 :: t,t_f

!——————————————————————————————————————————————
open (10,file='./1.txt',status='unknown')
open (20,file='./2.txt',status='unknown')
write(20,*)"t ","Energy ","Carter "

100 format(6f25.10)
200 format(1f10.4,2f25.15)
!————————————————————————————————————————————————————————————————————————————————————————————
!——————————————————————————————————————————————质量——————————————————————————————————————————————
M = 1
mu = 0
nu = 0
!——————————————————————————————————————————————自旋——————————————————————————————————————————————
!S_1 = 0.5 * m_1 **2
!S_2 = 0.5 * m_2 **2
!S_Kerr = S_1 + S_2
S_Kerr = 0.5*M**2
a = S_Kerr/M
chi_Kerr = a/M
!——————————————————————————————————————————————K参数——————————————————————————————————————————————
K = 1.447 - 1.715*nu - 3.246*nu**2
omega_1 = 0._16
omega_2 = 0._16
!——————————————————————————————————————————————Delta参数——————————————————————————————————————————————
Delta_0 = K*(nu*K-2)

Delta_1 = -2*(nu*K-1)*(K+Delta_0)

Delta_2 = 1/2.*Delta_1*(-4 * nu * K + Delta_1 + 4) - chi_Kerr**2*(nu*K - 1)**2*Delta_0

Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
                (-nu*K+Delta_2 + 1)-3*chi_Kerr**2*(nu*K-1)**2*Delta_1)


Delta_4 = 1/12.*(6*chi_Kerr**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2+3*Delta_1**4-8*(nu*K-1)*&
              Delta_1**3 - 12*Delta_2*Delta_1**2+12*(2*(nu*K-1)*Delta_2+Delta_3)*Delta_1 +&
              12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2+6*(Delta_2**2-4*Delta_3*(nu*K-1)))
!——————————————————————————————————————————————轨道参数——————————————————————————————————————————————
p = 7 * M 
e = 0.3
!R = p/(1+e*Cos(theta))
r_min = p/(1+e)
r_max = p/(1-e)
u = M/r_min

theta_min = Pi/4.
theta_max = Pi-theta_min
iota = theta_min

coordinates(:) = (/r_min,theta_min,0._16/)
Momentum(:) = (/0._16,0._16,0._16/)

h = 0.03 * M
t = 0
t_f = 10000 * M


!——————————————————————————————————————————————求解初始条件————————————————————————————————————————————————
g_1 = Metric((/r_min,theta_min,0._16/)) 																
g_2 = Metric((/r_max,theta_max,Pi/))																	

b1 = g_1(5)/g_1(1)
b2 = g_2(5)/g_2(1)

a1 = 1/Sqrt(-g_1(1))
a2 = 1/Sqrt(-g_2(1))

c1 = g_1(4) - g_1(5)**2/g_1(1)
c2 = g_2(4) - g_2(5)**2/g_2(1)

Momentum(1) = 0
Momentum(2) = 0
Momentum(3) = Sqrt((b1**2*(a1**2+a2**2)-2*b1*b2*(a1**2+a2**2)+b2**2*(a1**2+a2**2)-a1**4*c1+a1**2*a2**2*c1-&
              2*Sqrt((b1-b2)**2*a1**2*a2**2*((b1-b2)**2-(a1**2-a2**2)*(c1-c2)))+a1**2*a2**2*c2-a2**4*c2)/&
              (b1**4-4*b1**3*b2+b2**4+(a1**2*c1-a2**2*c2)**2-2*b2**2*(a1**2*c1+a2**2*c2)+4*b1*b2*(-b2**2+&
                a1**2*c1+a2**2*c2)+b1**2*(6*b2**2-2*(a1**2*c1+a2**2*c2))))

Energy = b1*Momentum(3) + a1*Sqrt(1+c1*Momentum(3)**2)

Carter = Cos(theta_min)**2*((1-Energy**2)*a**2+Momentum(3)**2/Sin(theta_min)**2)

write(*,*)"The Angular Momentum L equals to:",Momentum(3)
write(*,*)"The Effective Energy E equals to:",Energy
write(*,*)"The Carter Constant Q equals to:",Carter
!——————————————————————————————————————————————初始条件求解完毕————————————————————————————————————————————————

!————————————————————————————————————————————————-更新参数——————————————————————————————————————————————-————-
!1--->P_R_Dot,2--->P_Theta_Dot,3--->R_Dot,4--->Theta_Dot,5--->Phi_Dot 									    ｜
!   Derivative_of_parameters_ = D_Parameter(coordinates,Momentum) 											｜
!   coordinates(1) = coordinates(1) + R_Dot(coordinates,Momentum)*h/10										｜
!   coordinates(2) = coordinates(2) + Theta_Dot(coordinates,Momentum)*h/10									｜
!   coordinates(3) = coordinates(3) + Phi_Dot(coordinates,Momentum)*h/10									｜
!   Momentum(1) = Momentum(1) + P_R_Dot(coordinates,Momentum,Derivative_of_parameters_)*h/10				｜
!   Momentum(2) = Momentum(2) + P_Theta_Dot(coordinates,Momentum,Derivative_of_parameters_)*h/10			｜
!————————————————————————————————————————————————-更新参数——————————————————————————————————————————————-————-

do while (t .le. t_f)

  	Derivative_of_parameters_ = D_Parameter(coordinates,Momentum)

	do i = 1,5
	    K1(i) = X_P_Dot(i,coordinates,Momentum,Derivative_of_parameters_)
	 end do

	Temp_X(:) = coordinates(:) + h/2.*(/K1(3),K1(4),K1(5)/)
	Temp_P(:) = Momentum(:) + h/2.*(/K1(1),K1(2),0._16/)
    Derivative_of_parameters_ = D_Parameter(Temp_X,Temp_P)
	do i = 1,5
	   K2(i) = X_P_Dot(i,Temp_X,Temp_P,Derivative_of_parameters_)
	end do

	Temp_X(:) = coordinates(:) + h/2.*(/K2(3),K2(4),K2(5)/)
	Temp_P(:) = Momentum(:) + h/2.*(/K2(1),K2(2),0._16/)
	Derivative_of_parameters_ = D_Parameter(Temp_X,Temp_P)
	do i = 1,5
	   K3(i) = X_P_Dot(i,Temp_X,Temp_P,Derivative_of_parameters_)
	end do

	Temp_X(:) = coordinates(:) + h*(/K3(3),K3(4),K3(5)/)
	Temp_P(:) = Momentum(:) + h*(/K3(1),K3(2),0._16/)
	Derivative_of_parameters_ = D_Parameter(Temp_X,Temp_P)
	do i = 1,5
	   K4(i) = X_P_Dot(i,Temp_X,Temp_P,Derivative_of_parameters_)
	end do

	do i = 1,2
	   Momentum(i) = Momentum(i) + h/6.*(K1(i) + 2*K2(i) + 2*K3(i) + K4(i))
	end do

	do i = 1,3
	   coordinates(i) = coordinates(i) + h/6.*(K1(i+2) + 2*K2(i+2) + 2*K3(i+2) + K4(i+2))
	end do

	g = Metric(coordinates)

	alpha = 1/Sqrt(-g(1))

	beta = g(5)/g(1)

	gamma = g(4)-g(5)**2/g(1)

	Energy = beta*Momentum(3) + alpha*Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2)

	Carter = Momentum(2)**2 + Cos(coordinates(2))**2*((1-Energy**2)*a**2 + Momentum(3)**2/Sin(coordinates(2))**2)
!————————————————输出结果——————————————————————

	t = t + h
	write(10,100)coordinates,Momentum
	write(20,200)t,Energy,Carter

end do

write(*,*) "程序运行成功，可以查看分析数据了哦～"
close(10)
close(20)
contains

!————————————METRIC_TENSOR——————————————
!tt,rr,thetatheta,phiphi,tphi

function Metric(coordinates)
    implicit none
	COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
	real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 ,INTENT(IN):: coordinates(3)
    real * 16 :: m_1,m_2,u,Delta_Bar_u,Delta_u,R,theta
    real * 16 :: Delta_t,varpi_square,D_Minus_u,Delta_r,omega_tilde,Sigma,Lambda_t
    real * 16 :: g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi,Metric(5)

    R = coordinates(1)

    theta = coordinates(2)

    u = M/R

    Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2

    Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))

    Delta_t = R**2*Delta_u

    varpi_square = R**2 + a**2

    D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)

    Delta_r = Delta_t * D_Minus_u

    omega_tilde = 2*a*M*R + omega_1*nu*a*M**3/R + omega_2*nu*M*a**3/R

    Sigma = R**2+a**2*Cos(theta)**2

    Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

    g_tt = -Lambda_t/(Delta_t*Sigma)
    g_rr = Delta_r/Sigma
    g_thetatheta = 1/Sigma
    g_phiphi = 1/Lambda_t*(-omega_tilde**2/(Delta_t*Sigma)+Sigma/Sin(theta)**2)
    g_tphi = -omega_tilde/(Delta_t*Sigma)
    Metric(:)= (/g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi/)

    return
    end function
!——————————————----Derivative_Of_Metric————————————————————---

function D_Parameter(coordinates,Momentum)
  implicit none
  !D_parameter--->Dalpha,Dbeta,dgamma     D_Metric--->Dg_xx
  !dgtt_dr,dgtt_dth,dgrr_dr,dgrr_dth,dgthth_dr,dgthth_dth,dgphph_dr,dgphph_dth,dgtph_dr,dgtph_dth
  !dalpha_dr,dalpha_dth,dbeta_dr,dbeta_dth,dgammarr_dr,dgammarr_dth,dgammatt_dr,dgammatt_dth,dgammapp_dr,dgammapp_dth
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: u,Delta_u,Delta_t,dDeltau_du,varpi_square,dDeltat_dr,dDeltar_dr,dD_Minus_u_du
  real * 16 :: dLambdat_dr,Sigma,Lambda_t,omega_tilde,D_Minus_u,Delta_r,Delta_Bar_u
  real * 16 :: m_1,m_2,gtt,gtphi,R,theta
  real * 16 :: D_Metric(10),coordinates(3),Momentum(3),g(5)
  real * 16 :: D_Parameter(10)

  R = coordinates(1)

  theta = coordinates(2)

  u = M/R

  omega_tilde = 2*a*M*R + omega_1*nu*a*M**3/R + omega_2*nu*M*a**3/R

  varpi_square = R**2 + a**2

  Sigma = R**2 + a**2*Cos(theta)**2

  Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2

  Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))

  Delta_t = R**2*Delta_u

  D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)

  Delta_r = Delta_t * D_Minus_u

  Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

  gtt = -Lambda_t/(Delta_t*Sigma)

  gtphi = -omega_tilde/(Delta_t*Sigma)

  dDeltau_du =(nu*((-1+K*nu)**(-2)+(2*u)/(-1+K*nu)+chi_Kerr**2*u**2)*(Delta_1+u*(2*Delta_2+u*(3*Delta_3+4*Delta_4*u))))/ &
              (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u))))+2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
              (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))

  dDeltat_dr = 2*R*Delta_u - dDeltau_du*M

  dD_Minus_u_du = (12*u*nu+6*u**2*(26-3*nu)*nu)/(1+6*u**2*nu+2*u**3*(26-3*nu)*nu)

  dDeltar_dr = dDeltat_dr * D_Minus_u - Delta_t * dD_Minus_u_du * M /R**2

  dLambdat_dr = 4*R*varpi_square - a**2*Sin(theta)**2*dDeltat_dr
  !dg_tt/dr
  D_Metric(1) = (Lambda_t*(Sigma*dDeltat_dr + 2*R*Delta_t) - dLambdat_dr * Delta_t * Sigma)/(Delta_t*Sigma)**2
  !dg_rr/dr
  D_Metric(3) = (dDeltar_dr*Sigma - 2*R*Delta_r)/Sigma**2
  !dg_thth/dr
  D_Metric(5) = -2*R/Sigma**2
  !dt_phph/dr
  D_Metric(7) = (2*R*Lambda_t-dLambdat_dr*Sigma)/(Lambda_t**2*Sin(theta)**2)+ &
    (-8*a**2*M**2*R*Lambda_t*Delta_t*Sigma+omega_tilde**2*(dLambdat_dr*Delta_t*Sigma+ &
    Lambda_t*dDeltat_dr*Sigma+2*R*Lambda_t*Delta_t))/(Lambda_t*Sigma*Delta_t)**2
  !dg_tph/dr
  D_Metric(9) = (omega_tilde*(dDeltat_dr * Sigma + 2 * R * Delta_t)-2*a*M*Sigma*Delta_t)/(Delta_t*Sigma)**2

  !dg_tt/dth
  D_Metric(2) =  a**2*Sin(2*theta)*varpi_square*(Delta_t-varpi_square)/(Delta_t*Sigma**2)
  !dg_rr/dth
  D_Metric(4) = Delta_r*a**2*Sin(2*theta)/Sigma**2
  !dg_thth/dth
  D_Metric(6) = a**2*Sin(2*theta)/Sigma**2
  !dg_phph/dth
  D_Metric(8) = a**2*Sin(2*theta)/Lambda_t**2*(-Lambda_t/Sin(theta)**2+&
                Sigma*Delta_t/Sin(theta)**2-Lambda_t*Sigma/(a**2*Sin(theta)**4)-&
                omega_tilde**2*(Sigma*Delta_t+Lambda_t)/(Sigma**2*Delta_t))

  D_Metric(10) = -a**2*omega_tilde*Sin(2*theta)/(Sigma**2*Delta_t)


  D_Parameter(1) = 1/(2.*(-gtt)**1.5d0)*D_Metric(1)

  D_Parameter(2) = 1/(2.*(-gtt)**1.5d0)*D_Metric(2)

  D_Parameter(3) = (2*a*M*Lambda_t-dLambdat_dr*2*a*M*R)/Lambda_t**2

  D_Parameter(4) = omega_tilde*a**2*Sin(2*theta)*Delta_t/Lambda_t**2

  D_Parameter(5) = D_Metric(3)

  D_Parameter(6) = D_Metric(4)

  D_Parameter(7) = D_Metric(5)

  D_Parameter(8) = D_Metric(6)

  D_Parameter(9) = D_Metric(7) - (2*gtphi*gtt*D_Metric(9)-gtphi**2*D_Metric(1))/gtt**2

  D_Parameter(10) = D_Metric(8) - (2*gtphi*gtt*D_Metric(10)-gtphi**2*D_Metric(2))/gtt**2

  return
  end function
!——————————————————————P_Dot————————————————————————————————

function X_P_Dot(choice,coordinates,Momentum,Derivative_of_parameters)
  implicit none
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16,INTENT(IN) :: coordinates(3),Momentum(3),Derivative_of_parameters(10)
  integer :: choice
  real * 16 ::X_P_Dot
  if (choice .eq. 1) then
    X_P_Dot = P_R_Dot(coordinates,Momentum,Derivative_of_parameters)
  else if (choice .eq. 2) then
    X_P_Dot = P_Theta_Dot(coordinates,Momentum,Derivative_of_parameters)
  else if (choice .eq. 3) then
    X_P_Dot = R_Dot(coordinates,Momentum)
  else if (choice .eq. 4) then
    X_P_Dot = Theta_Dot(coordinates,Momentum)
  else if (choice .eq. 5)then
    X_P_Dot = Phi_Dot(coordinates,Momentum)
  end if  
  return
  end function
!——————————————————————P_R_Dot————————————————————————————————

function P_R_Dot(coordinates,Momentum,Derivative_of_parameters)
  implicit none
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16, INTENT(IN) :: coordinates(3),Momentum(3),Derivative_of_parameters(10)
  real * 16 :: The_guy_with_alpha,Energy,H_hat,dH_dHeff
  real * 16 :: alpha,beta,gamma
  real * 16 :: g(5),dbeta_dr,dalpha_dr,dgammarr_dr,dgammatt_dr,dgammapp_dr
  real * 16 :: P_R_Dot

  g = Metric(coordinates)

  alpha = 1/Sqrt(-g(1))

  beta = g(5)/g(1)

  gamma = g(4)-g(5)**2/g(1)

  The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2)

  Energy = beta*Momentum(3) + alpha*The_guy_with_alpha

  if (abs(nu) .ge. 1E-16) then
    H_hat = M/mu*(Sqrt(1+2*nu*(E-1))-1)
    dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
  elseif (abs(nu) .le. 1E-16) then
    dH_dHeff = M
  endif

  dalpha_dr = Derivative_of_parameters(1)

  dbeta_dr = Derivative_of_parameters(3)

  dgammarr_dr = Derivative_of_parameters(5)

  dgammatt_dr = Derivative_of_parameters(7)

  dgammapp_dr = Derivative_of_parameters(9)


  P_R_Dot = -dH_dHeff*(dbeta_dr*Momentum(3)+dalpha_dr*The_guy_with_alpha+alpha*(dgammarr_dr*Momentum(1)**2+&
            dgammatt_dr*Momentum(2)**2 + dgammapp_dr*Momentum(3)**2)/(2*The_guy_with_alpha))
  return
  end function
!——————————————————————P_Theta_Dot————————————————————————————————

function P_Theta_Dot(coordinates,Momentum,Derivative_of_parameters)
  implicit none
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16, INTENT(IN) :: coordinates(3),Momentum(3),Derivative_of_parameters(10)
  real * 16 :: Energy,H_hat,dH_dHeff
  real * 16 :: g(5),alpha,beta,gamma
  real * 16 :: dbeta_dth,dalpha_dth,dgammarr_dth,dgammatt_dth,dgammapp_dth,dgtt_dth,The_guy_with_alpha
  real * 16 :: P_Theta_Dot

  g = Metric(coordinates)

  alpha = 1/Sqrt(-g(1))

  beta = g(5)/g(1)

  gamma = g(4)-g(5)**2/g(1)

  The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2)

  Energy = beta*Momentum(3) + alpha*The_guy_with_alpha

  if (abs(nu) .ge. 1E-16) then
    H_hat = M/mu*(Sqrt(1+2*nu*(E-1))-1)
    dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
  elseif (abs(nu) .le. 1E-16) then
    dH_dHeff = M
  endif

  dalpha_dth = Derivative_of_parameters(2)

  dbeta_dth = Derivative_of_parameters(4)

  dgammarr_dth = Derivative_of_parameters(6)

  dgammatt_dth = Derivative_of_parameters(8)

  dgammapp_dth = Derivative_of_parameters(10)

  P_Theta_Dot = -dH_dHeff*(dbeta_dth*Momentum(3)+dalpha_dth*The_guy_with_alpha+alpha*(dgammarr_dth*Momentum(1)**2+&
              dgammatt_dth*Momentum(2)**2+dgammapp_dth*Momentum(3)**2)/(2*The_guy_with_alpha))
  return
  end function
!——————————————————————R_Dot————————————————————————————————

function R_Dot(coordinates,Momentum)
  implicit none
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16, INTENT(IN) :: coordinates(3),Momentum(3)
  real * 16 :: Energy,H_hat,dH_dHeff
  real * 16 :: g(5),alpha,beta,gamma,i,The_guy_with_alpha
  real * 16 :: R_Dot

  g = Metric(coordinates)

  alpha = 1/Sqrt(-g(1))

  beta = g(5)/g(1)

  gamma = g(4)-g(5)**2/g(1)

  The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2)

  Energy = beta*Momentum(3) + alpha*The_guy_with_alpha

  if (abs(nu) .ge. 1E-16) then
    H_hat = M/mu*(Sqrt(1+2*nu*(E-1))-1)
    dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
  elseif (abs(nu) .le. 1E-16) then
    dH_dHeff = M
  endif

  R_Dot = dH_dHeff*(alpha * g(2) * Momentum(1) / The_guy_with_alpha)

  return
  end function
!——————————————————————Theta_Dot————————————————————————————————

function Theta_Dot(coordinates,Momentum) 
  implicit none
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16, INTENT(IN) :: coordinates(3),Momentum(3)
  real * 16 :: H_hat,Energy,dH_dHeff
  real * 16 :: g(5),alpha,beta,gamma,The_guy_with_alpha
  real * 16 :: Theta_Dot

  g = Metric(coordinates)

  alpha = 1/Sqrt(-g(1))

  beta = g(5)/g(1)

  gamma = g(4)-g(5)**2/g(1)

  The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2)

  Energy = beta*Momentum(3) + alpha*The_guy_with_alpha

  if (abs(nu) .ge. 1E-16) then
    H_hat = M/mu*(Sqrt(1+2*nu*(E-1))-1)
    dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
  elseif (abs(nu) .le. 1E-16) then
    dH_dHeff = M
  endif

  Theta_Dot = dH_dHeff*(alpha*g(3)*Momentum(2)/The_guy_with_alpha)

  return
  end function
!——————————————————————Phi_Dot————————————————————————————————

function Phi_Dot(coordinates,Momentum)
  implicit none
  COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16, INTENT(IN) :: coordinates(3),Momentum(3)
  real * 16 :: Energy,H_hat,dH_dHeff
  real * 16 :: g(5),alpha,beta,gamma,The_guy_with_alpha
  real * 16 :: Phi_Dot

  g = Metric(coordinates)

  alpha = 1/Sqrt(-g(1))

  beta = g(5)/g(1)

  gamma = g(4)-g(5)**2/g(1)

  The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2)

  Energy = beta*Momentum(3) + alpha*The_guy_with_alpha

  if (abs(nu) .ge. 1E-16) then
    H_hat = M/mu*(Sqrt(1+2*nu*(E-1))-1)
    dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
  elseif (abs(nu) .le. 1E-16) then
    dH_dHeff = M
  endif

  Phi_Dot = dH_dHeff*(beta + alpha*gamma*Momentum(3) / The_guy_with_alpha)

  return
  end function

end program