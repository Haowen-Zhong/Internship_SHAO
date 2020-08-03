!做这个修改的目的是将RK4替换成RK78以求得到更高的精度

program EOB
implicit none
!定义物理量以及求解初始条件
  real * 16 , parameter :: Pi=3.1415926535897932384626433832795028841971693993751_16
  !COMMON M,mu,nu,a,chi_Kerr,K,omega_1,omega_2,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4

  COMMON /group1/ M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  COMMON /group3/ Lz
  real * 16 :: M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: Lz

  real * 16 :: coordinates(3),Momentum(3)
  real * 16 :: u,p,e,theta_min,theta_max,r_min,r_max
  real * 16 :: Energy,Carter,g_1(8),g_2(8),parameters(9),a1,a2,b1,b2,c1,c2
  integer :: i
  !——————————————————————以下参数为采用Runge-Kutta方法求解微分方程而定义——————————————————————————————
  real * 16 :: h,t,t_f

  !————————————————————————————————————————————————————————————————————————————————————————————

  open (10,file='./1.txt',status='unknown')
  open (20,file='./2.txt',status='unknown')
  write(20,*)"t ","Energy ","Carter "

  !————————————————————————————————————————————————————————————————————————————————————————————
  !——————————————————————————————————————————————质量——————————————————————————————————————————————
  M = 1
  nu = 1E-5_16
  mu = nu * M
  !——————————————————————————————————————————————自旋——————————————————————————————————————————————
  !S_1 = 0.5 * m_1 **2
  !S_2 = 0.5 * m_2 **2
  !S_Kerr = S_1 + S_2
  !S_Kerr = 0.5*M**2
  a = 0.5*M
  chi_Kerr = a/M

  !——————————————————————————————————————————————K参数——————————————————————————————————————————————
  K0 = 1.4467
  K = K0*(1-4*nu)**2 + 4*(1-2*nu)*nu             !文献45   式(6.11)
  omega_1 = 0._16
  omega_2 = 0._16
  !——————————————————————————————————————————————Delta参数——————————————————————————————————————————————
  Delta_0 = K*(nu*K-2)

  Delta_1 = -2*(nu*K-1)*(K+Delta_0)

  Delta_2 = 1/2.*Delta_1*(-4*nu*K+Delta_1+4)-chi_Kerr**2*(nu*K-1)**2*Delta_0

  Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
              (-nu*K+Delta_2+1)-3*chi_Kerr**2*(nu*K-1)**2*Delta_1)


  Delta_4 = 1/12.*(6*chi_Kerr**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2 + 3*Delta_1**4 - 8*(nu*K-1)*&
            Delta_1**3 - 12*Delta_2*Delta_1**2 + 12*(2*(nu*K-1)*Delta_2+Delta_3) * Delta_1 + &
            12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2 + 6*(Delta_2**2-4*Delta_3*(nu*K-1)))
  !——————————————————————————————————————————————轨道参数——————————————————————————————————————————————
  p = 5 * M 
  e = 0.2
  ! R = p/(1+e*Cos(theta))
  r_min = p/(1+e)
  r_max = p/(1-e)

  theta_min = Pi/4.
  theta_max = Pi-theta_min

  !——————————————————————————————————————————————求解初始条件————————————————————————————————————————————————
  parameters(:) = Parameter_Calculator((/r_min,theta_min,0._16/))
  g_1(:) = Metric(parameters,theta_min) 			
  parameters(:) = Parameter_Calculator((/r_max,theta_max,Pi/))													
  g_2(:) = Metric(parameters,theta_max)																	

  a1 = g_1(6)
  a2 = g_2(6)

  b1 = g_1(7)
  b2 = g_2(7)

  c1 = g_1(8)
  c2 = g_2(8)


  Momentum(1) = 0._16
  Momentum(2) = 0._16
  Momentum(3) = Sqrt((b1**2*(a1**2+a2**2)-2*b1*b2*(a1**2+a2**2)+b2**2*(a1**2+a2**2)-a1**4*c1+a1**2*a2**2*c1-&
                2*Sqrt((b1-b2)**2*a1**2*a2**2*((b1-b2)**2-(a1**2-a2**2)*(c1-c2)))+a1**2*a2**2*c2-a2**4*c2)/&
                (b1**4-4*b1**3*b2+b2**4+(a1**2*c1-a2**2*c2)**2-2*b2**2*(a1**2*c1+a2**2*c2)+4*b1*b2*(-b2**2+&
                a1**2*c1+a2**2*c2)+b1**2*(6*b2**2-2*(a1**2*c1+a2**2*c2))))

  Energy = b1*Momentum(3) + a1*Sqrt(1+c1*Momentum(3)**2)
  Carter = Cos(theta_min)**2*((1-Energy**2)*a**2+Momentum(3)**2/Sin(theta_min)**2)
  ! write(*,*)"The Angular Momentum L equals to:",Momentum(3)
  ! write(*,*)"The Effective Energy E equals to:",Energy
  Lz = Momentum(3)
  write(*,*)"The Angular Momentum L equals to:",Lz
  write(*,*)"The Effective Energy E equals to:",Energy
  write(*,*)"The Carter Constant Q equals to:",Carter
!求解演化

  h = 0.005 * M
  t = 0_16
  t_f = 20000 * M
  coordinates(:) = (/r_min,theta_min,0._16/)
  Momentum(:) = (/0._16,0._16,Lz/)
  call Evolution_RK78(h,coordinates,Momentum,t,t_f)

close(10)
close(20)
contains

!——————————————————————————————————————————————————————————————————————————————————————————————————————————
!——————————————————————————————————————————————————————————————————————————————————————————————————————————
!——————————————————————————————————————————————————————————————————————————————————————————————————————————

function Parameter_Calculator(coordinates)

  implicit None
  COMMON /group1/ M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4

  real * 16 ,  INTENT(IN)::coordinates(3)
  real * 16 :: R,theta,u,omega_tilde,Delta_Bar_u,Delta_u,Delta_t,varpi_square
  real * 16 :: D_Minus_u,Delta_r,Sigma,Lambda_t,gtt,gtphi,dDeltau_du,dDeltat_dr,dD_Minus_u_du
  real * 16 :: dDeltar_dr,dLambdat_dr
  real * 16 :: Parameter_Calculator(9)

    R = coordinates(1) 

    theta = coordinates(2) 

    u = M/R 

    varpi_square = R**2 + a**2 

    omega_tilde = 2*a*M*R + omega_1*nu*a*M**3/R + omega_2*nu*M*a**3/R 

    Sigma = R**2 + a**2*Cos(theta)**2 

    Delta_Bar_u = chi_Kerr**2*u**2 + 2*u/(nu*K-1) + 1/(nu*K-1)**2 

    Delta_u = Delta_Bar_u*(1+ nu*Delta_0 + nu*Log(1+u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + u*Delta_4))))) 

    Delta_t = R**2*Delta_u

    D_Minus_u = 1 + Log(1 + 6*nu*u**2 + 2*(26-3*nu)*nu*u**3) 

    Delta_r = Delta_t * D_Minus_u 

    Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2 

    dDeltau_du =(nu*((-1 + K*nu)**(-2) + (2*u)/(-1 + K*nu) + chi_Kerr**2*u**2)* &
              (Delta_1 + u*(2*Delta_2 + u*(3*Delta_3 + 4*Delta_4*u))))/ &
              (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))) +  &
              2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
              (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u))))) 

    dDeltat_dr = 2*R*Delta_u - dDeltau_du*M 

    dD_Minus_u_du = (12*u*nu + 6*nu*u**2*(26-3*nu))/(1 + 6*nu*u**2 + 2*nu*u**3*(26-3*nu)) 

    dDeltar_dr = dDeltat_dr * D_Minus_u - Delta_t * dD_Minus_u_du * M /R**2 

    dLambdat_dr = 4*R*varpi_square - a**2*Sin(theta)**2*dDeltat_dr 

    !                         1       2           3             4           5       6       7             8             9          
    Parameter_Calculator = (/Sigma, Lambda_t, varpi_square, omega_tilde, Delta_t, Delta_r, dDeltat_dr, dDeltar_dr, dLambdat_dr/)

  return
  end function
!————————————METRiC_TENSOR——————————————
!tt,rr,thetatheta,phiphi,tphi

function Metric(parameters,theta)
    implicit none
    real * 16 ,iNTENT(iN):: parameters(9),theta
    real * 16 :: g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi,alpha,beta,gamma
    real * 16 :: Metric(8)

    g_tt = - parameters(2)/(parameters(5)*parameters(1)) 
    g_rr = parameters(6)/parameters(1)
    g_thetatheta = 1/parameters(1)
    g_phiphi = 1/parameters(2)*(-parameters(4)**2/(parameters(5)*parameters(1)) + parameters(1)/Sin(theta)**2) 
    g_tphi = - parameters(4)/(parameters(5)*parameters(1)) 

    alpha = 1/Sqrt(-g_tt)

    beta = g_tphi/g_tt

    gamma = g_phiphi-g_tphi**2/g_tt

    Metric(:)= (/g_tt, g_rr, g_thetatheta, g_phiphi, g_tphi, alpha, beta, gamma/) 

    return
    end function
!——————————————----Derivative_Of_Metric————————————————————---

function D_Parameter(parameters,coordinates)
  implicit none

  COMMON /group1/ M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4

  real * 16 ,  INTENT(IN)::parameters(9),coordinates(3)
  real * 16 :: R,theta,omega_tilde,Delta_t,varpi_square
  real * 16 :: Delta_r,Sigma,Lambda_t,gtt,gtphi,dDeltau_du,dDeltat_dr,dDeltar_dr,dLambdat_dr
  real * 16 :: D_Metric(10),D_Parameter(10)

  R = coordinates(1)
  theta = coordinates(2)

  Sigma = parameters(1)
  Lambda_t = parameters(2)
  varpi_square = parameters(3)
  omega_tilde = parameters(4)
  Delta_t = parameters(5)
  Delta_r = parameters(6)
  dDeltat_dr = parameters(7)
  dDeltar_dr = parameters(8)
  dLambdat_dr = parameters(9)

  gtt = - parameters(2)/(parameters(5)*parameters(1))
  gtphi =  -parameters(4)/(parameters(5)*parameters(1)) 
  !dg_tt/dr  
  D_Metric(1) = (Lambda_t*(Sigma*dDeltat_dr + 2*R*Delta_t) - dLambdat_dr * Delta_t * Sigma)/(Delta_t*Sigma)**2
  !dg_rr/dr  
  D_Metric(3) = (dDeltar_dr*Sigma - 2*R*Delta_r)/Sigma**2
  !dg_thth/dr  
  D_Metric(5) = - 2*R/Sigma**2
  !dg_phph/dr 
  D_Metric(7) = (2*R*Lambda_t-dLambdat_dr*Sigma)/(Lambda_t**2*Sin(theta)**2)+ (-8*a**2*M**2*R*Lambda_t*Delta_t*&
                Sigma+omega_tilde**2*(dLambdat_dr*Delta_t*Sigma+Lambda_t*dDeltat_dr*Sigma&
                +2*R*Lambda_t*Delta_t))/(Lambda_t*Sigma*Delta_t)**2
  !dg_tph/dr 
  D_Metric(9) = 2*a*M*((dDeltat_dr*Sigma+2*R*Delta_t)*R-Delta_t*Sigma)/(Delta_t*Sigma)**2

  !dg_tt/dth 
  D_Metric(2) = a**2*Sin(2*theta)*(Delta_t*Sigma-Lambda_t)/(Delta_t*Sigma**2)
  !dg_rr/dth  
  D_Metric(4) = Delta_r*a**2*Sin(2*theta)/Sigma**2
  !dg_thth/dth   
  D_Metric(6) = a**2*Sin(2*theta)/Sigma**2
  !dg_phph/dth 
  D_Metric(8) = a**2*Sin(2*theta)/Lambda_t**2*&
                (-Lambda_t/Sin(theta)**2 + Sigma*Delta_t/Sin(theta)**2 - &
                Lambda_t*Sigma/(a**2*Sin(theta)**4)-&
                omega_tilde**2*(Sigma*Delta_t+Lambda_t)/(Sigma**2*Delta_t))
  !dg_tph/dth 
  D_Metric(10) = -a**2*omega_tilde*Sin(2*theta)/(Sigma**2*Delta_t)

  !dalpha_dr
  D_Parameter(1) = 1/(2.*(-gtt)**1.5_16)*D_Metric(1)
  !dalpha_dth
  D_Parameter(2) = 1/(2.*(-gtt)**1.5_16)*D_Metric(2)
  !dbeta_dr
  D_Parameter(3) = 2*a*M*(Lambda_t-dLambdat_dr*R)/Lambda_t**2
  !dbeta_dth
  D_Parameter(4) = 2*a**3*M*R*Sin(2*theta)*Delta_t/Lambda_t**2

  D_Parameter(5) = D_Metric(3)

  D_Parameter(6) = D_Metric(4)

  D_Parameter(7) = D_Metric(5)

  D_Parameter(8) = D_Metric(6)

  D_Parameter(9) = D_Metric(7) - (2*gtphi*gtt*D_Metric(9)-gtphi**2*D_Metric(1))/gtt**2

  D_Parameter(10) = D_Metric(8) - (2*gtphi*gtt*D_Metric(10)-gtphi**2*D_Metric(2))/gtt**2

  return
  end function

!—————————————————————— P_Dot & X_Dot ————————————————————————————————

function P_Dot(choice,g,Momentum,Derivative_of_parameters)
  implicit none

  COMMON /group1/ M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  real * 16 ::    M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2

  real * 16 ,  INTENT(IN) :: g(8),Momentum(3),Derivative_of_parameters(10)
  real * 16 :: The_guy_with_alpha,Energy,dH_dHeff
  real * 16 :: dalpha_dr,dbeta_dr,dgammarr_dr,dgammatt_dr,dgammapp_dr
  real * 16 :: dalpha_dth,dbeta_dth,dgammarr_dth,dgammatt_dth,dgammapp_dth
  real * 16 :: P_Dot,P_R_Dot,P_Theta_Dot
  integer :: choice

  The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+g(8)*Momentum(3)**2)

  Energy = g(7)*Momentum(3) + g(6)*The_guy_with_alpha

  dH_dHeff = 1/Sqrt(1+2*nu*(Energy-1))

  dalpha_dr = Derivative_of_parameters(1)

  dbeta_dr = Derivative_of_parameters(3)

  dgammarr_dr = Derivative_of_parameters(5)

  dgammatt_dr = Derivative_of_parameters(7)

  dgammapp_dr = Derivative_of_parameters(9)

  dalpha_dth = Derivative_of_parameters(2)

  dbeta_dth = Derivative_of_parameters(4)

  dgammarr_dth = Derivative_of_parameters(6)

  dgammatt_dth = Derivative_of_parameters(8)

  dgammapp_dth = Derivative_of_parameters(10)

  P_Theta_Dot = -dH_dHeff*(dbeta_dth*Momentum(3)+dalpha_dth*The_guy_with_alpha+g(6)*(dgammarr_dth*Momentum(1)**2+&
              dgammatt_dth*Momentum(2)**2 + dgammapp_dth*Momentum(3)**2)/(2*The_guy_with_alpha))

  P_R_Dot = -dH_dHeff*(dbeta_dr*Momentum(3)+dalpha_dr*The_guy_with_alpha+g(6)*(dgammarr_dr*Momentum(1)**2+&
            dgammatt_dr*Momentum(2)**2 + dgammapp_dr*Momentum(3)**2)/(2*The_guy_with_alpha))
  

  if (choice .eq. 1) then
    P_Dot = P_R_Dot
  else if (choice .eq. 2) then
    P_Dot = P_Theta_Dot
  end if 

  return
  end function  

function X_Dot(choice,g,Momentum)
    implicit none
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
    real * 16 :: M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2

    real * 16 ,  INTENT(IN) :: g(8),Momentum(3)
    real * 16 :: The_guy_with_alpha,Energy,dH_dHeff
    real * 16 :: X_Dot
    integer :: choice

    The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+g(8)*Momentum(3)**2)

    Energy = g(7)*Momentum(3) + g(6)*The_guy_with_alpha

    dH_dHeff = 1/Sqrt(1+2*nu*(Energy-1))

    if (choice .eq. 1) then
      !R_Dot
      X_Dot = dH_dHeff*(g(6) * g(2) * Momentum(1) / The_guy_with_alpha)
    else if (choice .eq. 2)then
      !Theta_Dot
      X_Dot = dH_dHeff*(g(6) * g(3) * Momentum(2) / The_guy_with_alpha)
    else if (choice .eq. 3)then
      !Phi_Dot
      X_Dot = dH_dHeff*(g(7) + g(6)*g(8)*Momentum(3) / The_guy_with_alpha)
    end if
 
    return
    end function
!—————————————————————— P_Dot & X_Dot ————————————————————————————————

Subroutine Evolution_RK78(h,coordinates,Momentum,t,t_f)
  implicit none
  COMMON /group1/ M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  real * 16 :: M,mu,nu,a,chi_Kerr,K0,K,omega_1,omega_2
  real * 16 :: h,coordinates(3),Momentum(3),t,t_f
  real * 16 :: X(5),X1(5),g(8)
  real * 16 :: Y0(5),Y1(5),Y2(5),Y3(5),Y4(5),Y5(5),Y6(5)
  real * 16 :: Y7(5),Y8(5),Y9(5),Y10(5),Y11(5),Y12(5)
  real * 16 :: Energy,Carter,parameters(9)
  integer :: i

  100 format(6f25.10)
  200 format(1f10.4,2f25.15)

  do while( t .le. t_f)

    X(1) = coordinates(1)
    X(2) = coordinates(2)
    X(3) = coordinates(3)
    X(4) = Momentum(1)
    X(5) = Momentum(2)
    do i = 1,5
      X1(i) = X(i)
    end do
    call Differential_Eq(X1,Y0)

    do i = 1,5
      X1(i) = X(i) + h*2_16/27_16*Y0(i)
    end do
    call Differential_Eq(X1,Y1)

    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)+3_16*Y1(i))/36_16
    end do 
    call Differential_Eq(X1,Y2)

    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)+3_16*Y2(i))/24_16
    end do
    call Differential_Eq(X1,Y3)

    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*20_16+(-Y2(i)+Y3(i))*75_16)/48_16
    end do 
    call Differential_Eq(X1,Y4)

    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)+Y3(i)*5_16+Y4(i)*4_16)/20_16   
    end do    
    call Differential_Eq(X1,Y5)

    do i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*25_16+Y3(i)*125_16-Y4(i)*260_16+Y5(i)*250_16)/108_16
    end do 
    call Differential_Eq(X1,Y6)


    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*93_16+Y4(i)*244_16-Y5(i)*200_16+Y6(i)*13_16)/900_16
    end do 
    call Differential_Eq(X1,Y7)


    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*180_16-Y3(i)*795_16+Y4(i)*1408_16-Y5(i)*1070_16 + &      
                  Y6(i)*67_16+Y7(i)*270_16)/90_16
    end do 
    call Differential_Eq(X1,Y8)

    do i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*455_16+Y3(i)*115_16-Y4(i)*3904_16+Y5(i)*3110_16 - &
                  Y6(i)*171_16+Y7(i)*1530_16-Y8(i)*45_16)/540_16
    end do
    call Differential_Eq(X1,Y9)


    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*2383_16-Y3(i)*8525_16+Y4(i)*17984_16-Y5(i)*15050_16+ &
                  Y6(i)*2133_16+Y7(i)*2250_16+Y8(i)*1125_16+Y9(i)*1800_16)/4100_16
    end do
    call Differential_Eq(X1,Y10)


    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*60_16-Y5(i)*600_16-Y6(i)*60_16+(Y8(i)-Y7(i)+&
                  2_16*Y9(i))*300_16)/4100_16
    end do
    call Differential_Eq(X1,Y11)


    do i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*1777_16-Y3(i)*8525_16+Y4(i)*17984_16-&
                  Y5(i)*14450_16+Y6(i)*2193_16+Y7(i)*2550_16 + &
                  Y8(i)*825_16+Y9(i)*1200_16+Y11(i)*4100_16)/4100_16
    end do
    call Differential_Eq(X1,Y12)

    do i = 1,5
      X(i) = X(i) + h*(Y5(i)*272_16+(Y6(i)+Y7(i))*216_16+(Y8(i)+Y9(i))*27_16 + &
                    (Y11(i)+Y12(i))*41_16)/840_16
    end do

    t = t + h
    coordinates(1) = X(1)
    coordinates(2) = X(2)
    coordinates(3) = X(3)
    Momentum(1) = X(4)
    Momentum(2) = X(5)

    parameters(:) = Parameter_Calculator(coordinates)
    g(:) = Metric(parameters,coordinates(2))

    Energy = g(7)*Momentum(3) + g(6)*Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+g(8)*Momentum(3)**2)
                  
    Carter = Momentum(2)**2 + Cos(coordinates(2))**2*((1-Energy**2)*a**2 + Momentum(3)**2/Sin(coordinates(2))**2)

    write(10,100)coordinates,Momentum
    write(20,200)t,Energy,Carter

  end do

  return
  end Subroutine


Subroutine Differential_Eq(X1,Y)
  implicit none
  COMMON /group3/ Lz
  real * 16 :: Lz
  real * 16 :: X1(5),Y(5),Temp_X(3),Temp_P(3)
  real * 16 :: g(8),Derivative_of_parameters(10)
  real * 16 :: parameters(9)

  Temp_X(:) = X1(1:3)
  Temp_P(:) = (/X1(4),X1(5),Lz/)
  parameters(:) = Parameter_Calculator(Temp_X)
  g(:) = Metric(parameters,Temp_X(2))
  Derivative_of_parameters = D_Parameter(parameters,Temp_X)

  Y(1) = X_Dot(1,g,Temp_P)                          !dot(R)
  Y(2) = X_Dot(2,g,Temp_P)                          !dot(theta)
  Y(3) = X_Dot(3,g,Temp_P)                          !dot(phi)
  Y(4) = P_Dot(1,g,Temp_P,Derivative_of_parameters) !dot(pr)
  Y(5) = P_Dot(2,g,Temp_P,Derivative_of_parameters) !dot(ptheta)

  return 

  end Subroutine

end program