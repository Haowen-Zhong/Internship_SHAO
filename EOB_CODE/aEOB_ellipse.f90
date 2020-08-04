!————————————————————————————————————————————说明—————————————————————————————————————————————————————————————————————————————
!                                                                                                                           _|
!说明:本程序中所有的公式全部源自于文献:PHYSICAL REVIEW D 86,024011 (2012)以及其所引用的参考文献12,34,35                             |_                          
!本程序分为几大部分  第一部分定义相关物理量:两天体质量,自旋,轨道倾角.并且利用近心点和远心点处哈密顿量所满足关系计算出初始条件                _|
!本程序中所有数组代入计算时都都使用(:)表示这是数组的计算                                                                           |_
!                                                                                                                            ｜
!————————————————————————————————————————————说明—————————————————————————————————————————————————————————————————————————————

program EOB
  implicit none

  real * 16 , parameter :: Pi=3.1415926535897932384626433832795028841971693993751_16


  COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
  COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

  real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
  real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS

  real * 16 :: coordinates(3),Momentum(3),u,p,e,theta_min,theta_max,R_min,R_max,R,theta,phi,psi,zeta
  real * 16 :: Rotation(3,3),S_1(3),S_2(3),S(3),norm_S_Star
  real * 16 :: alpha,beta,gamma,g(5),Derivative_of_parameters_(10),Q_4
  real * 16 :: The_guy_with_alpha,H_hat,Energy,Carter
  real * 16 :: S_Star_(3),S_Star_1(3),S_Star_2(3),parameters(9),Difference=1._16,epsilon = 1E-20_16
  real * 16 :: g_1(5),g_2(5),a1,a2,b1,b2,c1,c2,Z1,Z2,parameters1(9),parameters2(9),old_L
  integer :: i = 1

  !——————————————————————以下参数为采用Runge-Kutta方法求解微分方程而定义————————————————————————————————

  real * 16 :: K1(5),K2(5),K3(5),K4(5),h
  real * 16 :: Temp_X(3),Temp_P(3)
  real * 16 :: t,t_f

  !————————————————————————————————————————————————————————————————————————————————————————————————

  open (10,file='./coordinates.txt',status='unknown')
  open (20,file='./E&Q.txt',status='unknown')
  write(20,*)"t ","Energy ","Carter "

  100 format(6f25.10)
  200 format(1f10.4,2f25.15)
  !——————————————————————————————————————————————质量————————————————————————————————————————————————
  m_1 = 1_16
  m_2 = 1_16
  M = m_1 + m_2
  mu = m_1*m_2/M
  nu = mu/M
  !——————————————————————————————————————————自旋以及求解旋转矩阵——————————————————————————————————————
  d_SO = -69.5_16                                                 !文献式(39)
  d_SS = 2.75_16                                                  !文献式(39)
  S_1(:) = (/0.0_16,0.0_16,0.5_16/) * M **2
  S_2(:) = (/0._16,0._16,0._16/) * m_2 **2
  S_Kerr = S_1 + S_2                                              !最原始的坐标系中的自旋分量
  sigma_Star = m_1/m_2*S_2 + m_2/m_1*S_1                          !最原始的坐标系中的自旋分量
  norm_S_Kerr = Sqrt(S_Kerr(1)**2+S_Kerr(2)**2+S_Kerr(3)**2)
  if ((Abs(S_Kerr(1)) .ge. 1E-10) .or. (Abs(S_Kerr(2)) .ge. 1E-10)) then
    psi  = Acos(S_Kerr(3)/norm_S_Kerr)
    zeta = Atan(S_Kerr(2)/S_Kerr(1))
    Rotation(1,:) = (/-Sin(zeta),Cos(zeta),0._16/)
    Rotation(2,:) = (/-Cos(psi)*Cos(zeta),-Cos(psi)*Sin(zeta),Sin(psi)/)
    Rotation(3,:) = (/Sin(psi)*Cos(zeta),Sin(psi)*Sin(zeta),Cos(psi)/)
    S_Kerr(:) = matmul(Rotation,S_Kerr)                                    !旋转坐标系后的自旋分量
    sigma_Star(:) = matmul(Rotation,sigma_Star)                            !旋转坐标系后的自旋分量
  end if

  a = norm_S_Kerr/M
  chi_Kerr = a/M
  !——————————————————————————————————————————————K参数——————————————————————————————————————————————
  K = 1.447 - 1.715*nu - 3.246*nu**2                              !文献式(37)
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

  p = 5 * M 
  e = 0.2
  !R = p/(1+e*Cos(theta))
  R_min = p/(1+e)
  R_max = p/(1-e)
  u = M/R_min

  theta_min = Pi/4.
  theta_max = Pi-theta_min

  !从近心点，同时也是theta最小值处开始运动
  coordinates(:) = (/R_min,theta_min,0._16/)
  Momentum(:) = (/0._16,0._16,0._16/)

  !步长0.03M     模拟时长为10000M
  h = 0.05 * M
  t = 0
  t_f = 5000 * M


  !————————————————————————————————————————求解不考虑自旋效应时的初始条件——————————————————————————————————————————_
  !                                                                                                            |_
  !为了计算度规  首先调用函数计算所需要的参数                                                                        _|
  !下标1表示近心点       下标2表示远心点                                                                           |_
  !事实上由于轨道有倾角一定会产生进动 远心点的phi不可能是Pi 但是由于 Kerr度规是轴对称度规与phi无关 这里的取值不会影响结果     _|
  !我只要人为令初始状态时在近心点 并且phi为0即可 后续演化依旧正确                                                      |_
  !这里利用近心点和远心点处的哈密顿量相等来求解轨道角动量(H_eff本身守恒，这两个点处pr_dot,p_theta_dot=0更容易计算)          _|
  !                                                                                                            |
  !————————————————————————————————————————求解不考虑自旋效应时的初始条件——————————————————————————————————————————

  parameters(:) = Parameter_Calculator((/R_min,theta_min,0._16/))
  g_1 = Metric(parameters,(/R_min,theta_min,0._16/))

  parameters(:) = Parameter_Calculator((/R_max,theta_max,Pi/))											
  g_2 = Metric(parameters,(/R_max,theta_max,Pi/))																	

  b1 = g_1(5)/g_1(1)
  b2 = g_2(5)/g_2(1)

  a1 = 1/Sqrt(-g_1(1))
  a2 = 1/Sqrt(-g_2(1))

  c1 = g_1(4) - g_1(5)**2/g_1(1)
  c2 = g_2(4) - g_2(5)**2/g_2(1)

  Momentum(1) = 0._16
  Momentum(2) = 0._16
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
  parameters1 = Parameter_Calculator((/R_min,theta_min,0._16/))
  parameters2 = Parameter_Calculator((/R_max,theta_max,Pi/))
  write(*,*)"Energy caused by Spin",d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4+&
  H_SO_SS(parameters,coordinates,Momentum)+Spin_to_itself(parameters,coordinates,Momentum)
  do while((Difference .ge. epsilon) .and. (i .le. 200))
    old_L = Momentum(3)
    Z1 = H_SO_SS(parameters1,(/R_min,theta_min,0._16/),Momentum)/mu + &
         d_SS*nu*dot_product(S_Kerr,sigma_Star)/R_min**4+&
         Spin_to_itself(parameters1,(/R_min,theta_min,0_16/),Momentum)
    Z2 = H_SO_SS(parameters2,(/R_max,theta_max,Pi/),Momentum)/mu + &
         d_SS*nu*dot_product(S_Kerr,sigma_Star)/R_max**4+&
         Spin_to_itself(parameters2,(/R_max,theta_max,Pi/),Momentum)

    Momentum(3) = Sqrt(((((b2-b1)*old_L+Z2-Z1+a2*Sqrt(1+c2*old_L**2))/a1)**2-1)/c1)
    Difference = abs(Momentum(3)-old_L)
    i = i + 1
    write(*,*)"i=",i,"Difference=",Difference,"Momentum=",Momentum(3)
    write(*,*)"Eq=",(b2*Momentum(3) + a2*Sqrt(1+c2*Momentum(3)**2)+Z2)-&
                      (b1*Momentum(3) + a1*Sqrt(1+c1*Momentum(3)**2)+Z1)
  end do





  !————————————————————————————————————————————————-更新参数——————————————————————————————————————————————-————-
  !1--->P_R_Dot,2--->P_Theta_Dot,3--->R_Dot,4--->Theta_Dot,5--->Phi_Dot 									                    ｜
  !   Euler Method:可以用欧拉法判断是否存在系统误差                                                                |
  !   Derivative_of_parameters_ = D_Parameter(coordinates,Momentum) 											                    ｜
  !   coordinates(1) = coordinates(1) + R_Dot(coordinates,Momentum)*h/10										                  ｜
  !   coordinates(2) = coordinates(2) + Theta_Dot(coordinates,Momentum)*h/10									                ｜
  !   coordinates(3) = coordinates(3) + Phi_Dot(coordinates,Momentum)*h/10									                  ｜
  !   Momentum(1) = Momentum(1) + P_R_Dot(coordinates,Momentum,Derivative_of_parameters_)*h/10				        ｜
  !   Momentum(2) = Momentum(2) + P_Theta_Dot(coordinates,Momentum,Derivative_of_parameters_)*h/10		      	｜
  !————————————————————————————————————————————————-更新参数——————————————————————————————————————————————-————-

  i = 1
  coordinates(:) = (/R_min,theta_min,0._16/)
  do while (t .le. t_f)

      parameters = Parameter_Calculator(coordinates)
      g(:) = Metric(parameters,coordinates)

      alpha = 1/Sqrt(-g(1))

      beta = g(5)/g(1)

      gamma = g(4)-g(5)**2/g(1)

      Q_4 = 2*(4-3*nu)*nu*M**2/parameters(1)*Momentum(1)**4

      The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)

      Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + &
               H_SO_SS(parameters,coordinates,Momentum)/mu + &
               d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4 + &
               Spin_to_itself(parameters,coordinates,Momentum)
             
      Carter = Momentum(2)**2 + Cos(coordinates(2))**2*&
               ((1-Energy**2)*a**2+Momentum(3)**2/Sin(coordinates(2))**2)

      write(10,100)coordinates,Momentum
      write(20,200)t,Energy,Carter

    Derivative_of_parameters_ = D_Parameter(coordinates,Momentum)
  	do i = 1,5
  	    K1(i) = X_P_Dot(i,parameters,coordinates,Momentum,Derivative_of_parameters_)
  	end do

  	Temp_X(:) = coordinates(:) + h/2.*(/K1(3),K1(4),K1(5)/)
  	Temp_P(:) = Momentum(:) + h/2.*(/K1(1),K1(2),0._16/)
    parameters = Parameter_Calculator(Temp_X)
    Derivative_of_parameters_ = D_Parameter(Temp_X,Temp_P)
  	do i = 1,5
  	   K2(i) = X_P_Dot(i,parameters,Temp_X,Temp_P,Derivative_of_parameters_)
  	end do

  	Temp_X(:) = coordinates(:) + h/2.*(/K2(3),K2(4),K2(5)/)
  	Temp_P(:) = Momentum(:) + h/2.*(/K2(1),K2(2),0._16/)
  	Derivative_of_parameters_ = D_Parameter(Temp_X,Temp_P)
    parameters = Parameter_Calculator(Temp_X)
  	do i = 1,5
  	   K3(i) = X_P_Dot(i,parameters,Temp_X,Temp_P,Derivative_of_parameters_)
  	end do

  	Temp_X(:) = coordinates(:) + h*(/K3(3),K3(4),K3(5)/)
  	Temp_P(:) = Momentum(:) + h*(/K3(1),K3(2),0._16/)
  	Derivative_of_parameters_ = D_Parameter(Temp_X,Temp_P)
    parameters = Parameter_Calculator(Temp_X)
  	do i = 1,5
  	   K4(i) = X_P_Dot(i,parameters,Temp_X,Temp_P,Derivative_of_parameters_)
  	end do

  	do i = 1,2
  	   Momentum(i) = Momentum(i) + h/6.*(K1(i) + 2*K2(i) + 2*K3(i) + K4(i))
  	end do

  	do i = 1,3
  	   coordinates(i) = coordinates(i) + h/6.*(K1(i+2) + 2*K2(i+2) + 2*K3(i+2) + K4(i+2))
  	end do
  !————————————————输出结果——————————————————————

  	t = t + h

  end do

  write(*,*) "程序运行成功，可以查看分析数据了哦～"
  close(10)
  close(20)
  contains

  !————————————————————————————————-PARAMETER_CALCULATOR————————————————————————————————————————————————————
  !通过这个函数直接计算出后续大量使用到的参数，防止参数被重复计算，浪费计算资源                                       ｜
  !1--->Sigma      2--->varpi_square    3--->omega_tilde   4---> Delta_t    5--->Delta_r                   ｜
  !6--->Lambda_t   7--->Delta_t'        8--->Lambda_t'     9--->Delta_r'                                    |
  !————————————————————————————————-PARAMETER_CALCULATOR————————————————————————————————————————————————————

  function Parameter_Calculator(coordinates)
      implicit none
      COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
      COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4

      real * 16 ,INTENT(IN):: coordinates(3)
      real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
      real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
      real * 16 :: u,Delta_Bar_u,Delta_u,R,theta
      real * 16 :: Delta_t,varpi_square,D_Minus_u,Delta_r,omega_tilde,Sigma,Lambda_t
      real * 16 :: dDeltau_du,dDeltar_dr,dD_Minus_u_du,dLambdat_dr,dDeltat_dr
      real * 16 :: Parameter_Calculator(9)

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

      dDeltau_du =(nu*((-1+K*nu)**(-2)+(2*u)/(-1+K*nu)+chi_Kerr**2*u**2)*(Delta_1+u*(2*Delta_2+u*(3*Delta_3+4*Delta_4*u))))/ &
                (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u))))+2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
                (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))

      dDeltat_dr = 2*R*Delta_u - dDeltau_du*M

      dD_Minus_u_du = (12*u*nu+6*u**2*(26-3*nu)*nu)/(1 + 6*u**2*nu + 2*nu*u**3*(26-3*nu))

      dDeltar_dr = dDeltat_dr * D_Minus_u - Delta_t * dD_Minus_u_du * M /R**2

      dLambdat_dr = 4*R*varpi_square - a**2*Sin(theta)**2*dDeltat_dr

      Parameter_Calculator(:) = (/Sigma,varpi_square,omega_tilde,Delta_t,Delta_r,Lambda_t,dDeltat_dr,dLambdat_dr,dDeltar_dr/)

      return
      end function
  !——————————————————————————————————————METRIC_TENSOR——————————————————————————————————————————————————————
  !tt,rr,thetatheta,phiphi,tphi

  function Metric(parameters,coordinates)
      implicit none
      real * 16 ,INTENT(IN):: parameters(9),coordinates(3)
      real * 16 :: Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta
      real * 16 :: g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi,Metric(5)

      theta = coordinates(2)

      Sigma = parameters(1)
      omega_tilde = parameters(3)
      Delta_t = parameters(4)
      Delta_r = parameters(5)
      Lambda_t = parameters(6)

      g_tt = -Lambda_t/(Delta_t*Sigma)
      g_rr = Delta_r/Sigma
      g_thetatheta = 1/Sigma
      g_phiphi = 1/Lambda_t*(-omega_tilde**2/(Delta_t*Sigma)+Sigma/Sin(theta)**2)
      g_tphi = -omega_tilde/(Delta_t*Sigma)
      Metric(:)= (/g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi/)

      return
      end function
  !——————————————----Derivative_Of_Metric———————————————————————

  function D_Parameter(coordinates,Momentum)
    implicit none
    !D_parameter--->Dalpha,Dbeta,dgamma     D_Metric--->Dg_xx
    !dgtt_dr,dgtt_dth,dgrr_dr,dgrr_dth,dgthth_dr,dgthth_dth,dgphph_dr,dgphph_dth,dgtph_dr,dgtph_dth
    !dalpha_dr,dalpha_dth,dbeta_dr,dbeta_dth,dgammarr_dr,dgammarr_dth,dgammatt_dr,dgammatt_dth,dgammapp_dr,dgammapp_dth

    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4

    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4

    real * 16 :: Delta_t,dDeltau_du,varpi_square,dDeltat_dr,dDeltar_dr,dD_Minus_u_du
    real * 16 :: dLambdat_dr,Sigma,Lambda_t,omega_tilde,D_Minus_u,Delta_r,Delta_u
    real * 16 :: gtt,gtphi,R,theta
    real * 16 :: D_Metric(10),coordinates(3),Momentum(3),g(5),parameters(9)
    real * 16 :: D_Parameter(10)

    R = coordinates(1)

    theta = coordinates(2)

    parameters(:) = Parameter_Calculator(coordinates)
    Sigma = parameters(1)
    varpi_square = parameters(2)
    omega_tilde = parameters(3)
    Delta_t = parameters(4)
    Delta_r = parameters(5)
    Lambda_t = parameters(6)
    dDeltat_dr = parameters(7)
    dLambdat_dr = parameters(8)
    dDeltar_dr = parameters(9)
    gtt = -Lambda_t/(Delta_t*Sigma)
    gtphi = -omega_tilde/(Delta_t*Sigma)

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
  !——————————————————————P_Dot——————————————————————————————————

  function X_P_Dot(choice,parameters,coordinates,Momentum,Derivative_of_parameters)
    implicit none
    real * 16,INTENT(IN) :: parameters(9),coordinates(3),Momentum(3),Derivative_of_parameters(10)
    real * 16 ::X_P_Dot
    integer :: choice
    if (choice .eq. 1) then
      X_P_Dot = P_R_Dot(parameters,coordinates,Momentum,Derivative_of_parameters)
    else if (choice .eq. 2) then
      X_P_Dot = P_Theta_Dot(parameters,coordinates,Momentum,Derivative_of_parameters)
    else if (choice .eq. 3) then
      X_P_Dot = R_Dot(parameters,coordinates,Momentum)
    else if (choice .eq. 4) then
      X_P_Dot = Theta_Dot(parameters,coordinates,Momentum)
    else if (choice .eq. 5)then
      X_P_Dot = Phi_Dot(parameters,coordinates,Momentum)
    end if  
    return
    end function
  !——————————————————————P_R_Dot————————————————————————————————

  function P_R_Dot(parameters,coordinates,Momentum,Derivative_of_parameters)
    implicit none
   
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

    real * 16, INTENT(IN) :: parameters(9),coordinates(3),Momentum(3),Derivative_of_parameters(10)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS

    real * 16 :: The_guy_with_alpha,Energy,H_hat,dH_dHeff,Q_4
    real * 16 :: alpha,beta,gamma,dQ_dr,dH_Spin_dr,dSpin_to_itself_dr
    real * 16 :: g(5),dbeta_dr,dalpha_dr,dgammarr_dr,dgammatt_dr,dgammapp_dr
    real * 16 :: P_R_Dot


    g(:) = Metric(parameters,coordinates)

    alpha = 1/Sqrt(-g(1))

    beta = g(5)/g(1)

    gamma = g(4)-g(5)**2/g(1)

    Q_4 = 2*(4-3*nu)*nu*M**2/parameters(1)*Momentum(1)**4
    

    The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)

    dalpha_dr = Derivative_of_parameters(1)

    dbeta_dr = Derivative_of_parameters(3)

    dgammarr_dr = Derivative_of_parameters(5)

    dgammatt_dr = Derivative_of_parameters(7)

    dgammapp_dr = Derivative_of_parameters(9)



    dQ_dr = -4*(4-3*nu)*nu*M**2*coordinates(1)*Momentum(1)**4/parameters(1)**2

    dH_Spin_dr = numerical_differentiation(1,coordinates,Momentum)

    dSpin_to_itself_dr = numerical_differentiation(2,coordinates,Momentum)

    Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + &
             H_SO_SS(parameters,coordinates,Momentum)/mu + &
             d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4 + &
             Spin_to_itself(parameters,coordinates,Momentum)

    if (abs(nu) .ge. 1E-16) then
      H_hat = M/mu*(Sqrt(1+2*nu*(Energy-1))-1)
      dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
    elseif (abs(nu) .le. 1E-16) then
      dH_dHeff = M
    endif

    P_R_Dot = -dH_dHeff*(dbeta_dr*Momentum(3)+dalpha_dr*The_guy_with_alpha+alpha*(dgammarr_dr*Momentum(1)**2+&
              dgammatt_dr*Momentum(2)**2 + dgammapp_dr*Momentum(3)**2 + dQ_dr)/(2*The_guy_with_alpha) + &
              dH_Spin_dr/mu - 4*d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**5  + dSpin_to_itself_dr)
    return
     end function
  !——————————————————————P_Theta_Dot————————————————————————————

  function P_Theta_Dot(parameters,coordinates,Momentum,Derivative_of_parameters)
    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

    real * 16, INTENT(IN) :: parameters(9),coordinates(3),Momentum(3),Derivative_of_parameters(10)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS
    real * 16 :: Energy,H_hat,dH_dHeff,Q_4
    real * 16 :: g(5),alpha,beta,gamma,dQ_dth,dH_Spin_dth,dSpin_to_itself_dth
    real * 16 :: dbeta_dth,dalpha_dth,dgammarr_dth,dgammatt_dth,dgammapp_dth,dgtt_dth,The_guy_with_alpha
    real * 16 :: P_Theta_Dot

    g(:) = Metric(parameters,coordinates)

    alpha = 1/Sqrt(-g(1))

    beta = g(5)/g(1)

    gamma = g(4)-g(5)**2/g(1)

    Q_4 = 2*(4-3*nu)*nu*M**2/parameters(1)*Momentum(1)**4
    
    The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)

    Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + &
             H_SO_SS(parameters,coordinates,Momentum)/mu + &
             d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4 + &
             Spin_to_itself(parameters,coordinates,Momentum)

    if (abs(nu) .ge. 1E-16) then
      H_hat = M/mu*(Sqrt(1+2*nu*(Energy-1))-1)
      dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
    elseif (abs(nu) .le. 1E-16) then
      dH_dHeff = M
    endif

    dalpha_dth = Derivative_of_parameters(2)

    dbeta_dth = Derivative_of_parameters(4)

    dgammarr_dth = Derivative_of_parameters(6)

    dgammatt_dth = Derivative_of_parameters(8)

    dgammapp_dth = Derivative_of_parameters(10)

    dQ_dth = 2*(4-3*nu)*nu*M**2*a**2*Sin(2*coordinates(2))*Momentum(1)**4/parameters(1)**2


    dH_Spin_dth = numerical_differentiation(5,coordinates,Momentum)

    dSpin_to_itself_dth = numerical_differentiation(6,coordinates,Momentum)


    P_Theta_Dot = -dH_dHeff*(dbeta_dth*Momentum(3)+dalpha_dth*The_guy_with_alpha+alpha*(dgammarr_dth*Momentum(1)**2+&
                dgammatt_dth*Momentum(2)**2+dgammapp_dth*Momentum(3)**2+dQ_dth)/(2*The_guy_with_alpha) + &
                dH_Spin_dth/mu + dSpin_to_itself_dth)
    return
    end function
   !————————————————————————R_Dot————————————————————————————————

  function R_Dot(parameters,coordinates,Momentum)
    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS
    real * 16, INTENT(IN) :: parameters(9),coordinates(3),Momentum(3)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS

    real * 16 :: g(5),alpha,beta,gamma,The_guy_with_alpha,dQ_dpr,dH_Spin_dpr,dSpin_to_itself_dpr
    real * 16 :: Energy,H_hat,dH_dHeff,Q_4
    real * 16 :: R_Dot

    g(:) = Metric(parameters,coordinates)

    alpha = 1/Sqrt(-g(1))

    beta = g(5)/g(1)

    gamma = g(4)-g(5)**2/g(1)

    Q_4 = 2*(4-3*nu)*nu*M**2/parameters(1)*Momentum(1)**4
    
    The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2+Q_4) 

    Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + &
             H_SO_SS(parameters,coordinates,Momentum)/mu + &
             d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4 + &
             Spin_to_itself(parameters,coordinates,Momentum)

    if (abs(nu) .ge. 1E-16) then
      H_hat = M/mu*(Sqrt(1+2*nu*(Energy-1))-1)
      dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
    elseif (abs(nu) .le. 1E-16) then
      dH_dHeff = M
    endif

    dQ_dpr = 8*(4-3*nu)*nu*M**2*Momentum(1)**3/parameters(1)

    dH_Spin_dpr = numerical_differentiation(3,coordinates,Momentum)

    dSpin_to_itself_dpr = numerical_differentiation(4,coordinates,Momentum)

    R_Dot = dH_dHeff*(alpha * (2*g(2)*Momentum(1)+dQ_dpr)/(2*The_guy_with_alpha)+dH_Spin_dpr/mu + &
            dSpin_to_itself_dpr)

    return
    end function
   !————————————————————Theta_Dot————————————————————————————————

  function Theta_Dot(parameters,coordinates,Momentum)
    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

    real * 16, INTENT(IN) :: parameters(9),coordinates(3),Momentum(3)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS
    real * 16 :: g(5),alpha,beta,gamma,The_guy_with_alpha,dH_Spin_dpth,dSpin_to_itself_dpth
    real * 16 :: H_hat,Energy,dH_dHeff,Q_4
    real * 16 :: Theta_Dot


    g(:) = Metric(parameters,coordinates)

    alpha = 1/Sqrt(-g(1))

    beta = g(5)/g(1)

    gamma = g(4)-g(5)**2/g(1)

    Q_4 = 2*(4-3*nu)*nu*M**2/parameters(1)*Momentum(1)**4
    
    The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)


    Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + &
             H_SO_SS(parameters,coordinates,Momentum)/mu + &
             d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4 + &
             Spin_to_itself(parameters,coordinates,Momentum)

    if (abs(nu) .ge. 1E-16) then
      H_hat = M/mu*(Sqrt(1+2*nu*(Energy-1))-1)
      dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
    elseif (abs(nu) .le. 1E-16) then
      dH_dHeff = M
    endif

    dH_Spin_dpth = numerical_differentiation(7,coordinates,Momentum)

    dSpin_to_itself_dpth = numerical_differentiation(8,coordinates,Momentum)

    Theta_Dot = dH_dHeff*(alpha*g(3)*Momentum(2)/The_guy_with_alpha+dH_Spin_dpth/mu + dSpin_to_itself_dpth)

    return
    end function
   !——————————————————————Phi_Dot————————————————————————————————

  function Phi_Dot(parameters,coordinates,Momentum)
    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

    real * 16, INTENT(IN) :: parameters(9),coordinates(3),Momentum(3)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS
    real * 16 :: g(5),alpha,beta,gamma,The_guy_with_alpha,dH_Spin_dpph,dSpin_to_itself_dpph
    real * 16 :: Energy,H_hat,dH_dHeff,Q_4
    real * 16 :: Phi_Dot

    g(:) = Metric(parameters,coordinates)

    alpha = 1/Sqrt(-g(1))

    beta = g(5)/g(1)

    gamma = g(4)-g(5)**2/g(1)

    Q_4 = 2*(4-3*nu)*nu*M**2/parameters(1)*Momentum(1)**4
    
    The_guy_with_alpha = Sqrt(1+g(2)*Momentum(1)**2+g(3)*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)

    Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + &
             H_SO_SS(parameters,coordinates,Momentum)/mu + &
             d_SS*nu*dot_product(S_Kerr,sigma_Star)/coordinates(1)**4 + &
             Spin_to_itself(parameters,coordinates,Momentum)

    if (abs(nu) .ge. 1E-16) then
      H_hat = M/mu*(Sqrt(1+2*nu*(Energy-1))-1)
      dH_dHeff = M**2*nu/(mu*(mu*H_hat+M))
    elseif (abs(nu) .le. 1E-16) then
      dH_dHeff = M
    endif

    dH_Spin_dpph = numerical_differentiation(9,coordinates,Momentum)

    dSpin_to_itself_dpph = numerical_differentiation(10,coordinates,Momentum)

    Phi_Dot = dH_dHeff*(beta + alpha*gamma*Momentum(3)/The_guy_with_alpha+dH_Spin_dpph/mu + dSpin_to_itself_dpph)

    return
    end function

  !————————————————————————S_star————————————————————————————————

  function S_Star(parameters,coordinates,Momentum)

    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

    real * 16 ,INTENT(IN) :: parameters(9),coordinates(3),Momentum(3)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS
    real * 16 :: Delta_sigma_1(3),Delta_sigma_2(3),Delta_sigma_3(3),S_Star(3)
    real * 16 :: u,R,theta,p_R,p_theta,p_phi,Q,g(5)

    g(:) = Metric(parameters,coordinates)
    R = coordinates(1)
    u = M/R
    theta = coordinates(2)

    p_R = Momentum(1)
    p_theta = Momentum(2)
    p_phi = Momentum(3)

    Q = 1 + g(2)*p_R**2 + g(3)*p_theta**2 + (g(4)-g(5)**2/g(1))*p_phi**2


    Delta_sigma_1(:) = sigma_Star(:) * (7/6.*nu*u + nu/3.*(Q-1)-5/2.*nu*g(2)*p_R**2)+&
                       S_Kerr(:) * (-2/3.*nu*u+1/4.*nu*(Q-1)-3*nu*g(2)*p_R**2)

    Delta_sigma_2(:) = sigma_Star(:) * (1/36.*(353*nu-27*nu**2)*u**2+5*nu**2*g(2)**2*p_R**4-&
                       1/72.*(23*nu+3*nu**2)*(Q-1)**2+1/36.*(-103*nu+60*nu**2)*u*(Q-1)+&
                       1/12.*(16*nu-21*nu**2)*g(2)*p_R**2*(Q-1)+1/12.*(47*nu-54*nu**2)*u*g(2)*p_R**2)+&
                       S_Kerr(:) * (1/9.*(-56*nu-21*nu**2)*u**2+45/8.*nu**2*g(2)**2*p_R**4-&
                        5/16.*nu*(Q-1)**2+1/36.*(-109*nu+51*nu**2)*u*(Q-1)+&
                        1/8.*(2*nu-13*nu**2)*g(2)*p_R**2*(Q-1)-1/24.*(16*nu+147*nu**2)*u*g(2)*p_R**2)

    Delta_sigma_3(:) = d_SO*nu/R**3*sigma_Star(:)

    S_Star(:) = sigma_Star(:) + Delta_sigma_1(:) + Delta_sigma_2(:) + Delta_sigma_3(:)
    S_Star(:) = (/0._16,0._16,0._16/)
    return 
    end function

  !————————————————————————H_SO————————————————————————————————

  function H_SO_SS(parameters,coordinates,Momentum)
    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS

    real * 16 ,INTENT(IN) :: parameters(9),coordinates(3),Momentum(3)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
    real * 16 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS
    real * 16 :: S(3),xi(3),v(3),n(3)
    real * 16 :: R,theta,phi,p_R,p_theta,p_phi,Q,u,varpi_square,Delta_t,Delta_r,omega_cos,omega_tilde
    real * 16 :: Sigma,Lambda_t,Lambda_t_prime,omega,omega_r,Delta_t_prime
    real * 16 :: e_2nu,e_2mutilde,B_rtilde,B_tilde,J_tilde,mu_r,mu_cos,nu_r,nu_cos,H_SO,H_SS
    real * 16 :: H_SO_SS

    R = coordinates(1)
    theta = coordinates(2)
    phi = coordinates(3)
    u = M/R
    n(:) = (/Sin(theta)*Cos(phi),Sin(theta)*Sin(phi),Cos(theta)/)

    p_R = Momentum(1)
    p_theta = Momentum(2)
    p_phi = Momentum(3)
    !************************************************************************************************************
    !1--->Sigma      2--->varpi_square    3--->omega_tilde   4---> Delta_t    5--->Delta_r                     !*       
    !6--->Lambda_t   7--->Delta_t'        8--->Lambda_t'     9--->Delta_u                                      !*
    Sigma = parameters(1)                                                                                      !*
    varpi_square = parameters(2)                                                                               !*
    omega_tilde = parameters(3)                                                                                !*
    Delta_t = parameters(4)                                                                                    !*
    Delta_r = parameters(5)                                                                                    !*
    Lambda_t = parameters(6)                                                                                   !*
    Delta_t_prime = parameters(7)                                                                              !*
    Lambda_t_prime = parameters(8)                                                                             !*
    !************************************************************************************************************
    omega = omega_tilde / Lambda_t                                                                             !*
    omega_r = (-Lambda_t_prime*omega_tilde+Lambda_t*2*a*M)/Lambda_t**2                                         !* 
    omega_cos = -2*a**2*Cos(theta)*Delta_t*omega_tilde/Lambda_t**2                                             !* 
    e_2nu = Delta_t*Sigma/Lambda_t                                                                             !*
    e_2mutilde = Sigma                                                                                         !*
    B_tilde = Sqrt(Delta_t)                                                                                    !*
    B_rtilde = (Sqrt(Delta_r)*Delta_t_prime-2*Delta_t)/(2*Sqrt(Delta_t*Delta_r))                               !*
    J_tilde = Sqrt(Delta_r)                                                                                    !*
    mu_r = R/Sigma - 1/Sqrt(Delta_r)                                                                           !*
    mu_cos = a**2*Cos(theta)/Sigma                                                                             !*
    nu_r = R/Sigma + varpi_square*(varpi_square*Delta_t_prime-4*R*Delta_t)/(2*Lambda_t*Delta_t)                !*
    nu_cos = a**2*varpi_square*Cos(theta)*(varpi_square-Delta_t)/(Lambda_t*Sigma)                              !*
    !************************************************************************************************************
    !************************************************************************************************************
    Q  = 1 + Delta_r*p_R**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma                 !*
    S(:) = nu * S_Star(parameters,coordinates,Momentum)                                                        !*
    xi(:) = Sin(theta)*(/-Sin(phi),Cos(phi),0._16/)                                                            !*
    v(:) = -Sin(theta)*(/Cos(theta)*Cos(phi),Cos(theta)*Sin(phi),-Sin(theta)/)                                 !*  
    !************************************************************************************************************
    H_SO = e_2nu/Sqrt(e_2mutilde)*(Sqrt(e_2mutilde*e_2nu)-B_tilde)*p_phi*&
           dot_product(S,S_Kerr/norm_S_Kerr)/(B_tilde**2*Sqrt(Q)*Sin(theta)**2)+&
           Sqrt(e_2nu)/(e_2mutilde*B_tilde**2*(Sqrt(Q)+1)*Sqrt(Q)*Sin(theta)**2)*&
           (dot_product(S,xi)*J_tilde*(-mu_r*Sin(theta)*p_theta*(Sqrt(Q)+1)-mu_cos*p_R*Sin(theta)**2-&
           Sqrt(Q)*(-nu_r*Sin(theta)*p_theta+(mu_cos-nu_cos)*p_R*Sin(theta)**2))*B_tilde**2 +&
           Sqrt(e_2mutilde*e_2nu)*p_phi*(2*Sqrt(Q)+1)*B_tilde*(J_tilde*nu_r*dot_product(S,v)-nu_cos*dot_product(S,n)*&
           Sin(theta)**2)-J_tilde*B_rtilde*Sqrt(e_2mutilde*e_2nu)*p_phi*(Sqrt(Q)+1)*dot_product(S,v))
    !************************************************************************************************************
    H_SS = omega*dot_product(S,S_Kerr/norm_S_Kerr) + J_tilde*omega_r/(2*B_tilde*(Sqrt(Q)+1)*Sqrt(Q)*Sin(theta)**2*&
           (Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*(Sqrt(e_2mutilde*e_2nu)*Sin(theta)*p_theta*p_phi*dot_product(S,xi)*B_tilde+&
           e_2nu*e_2mutilde*p_phi**2*dot_product(S,v)+e_2mutilde*(1+Sqrt(Q))*Sqrt(Q)*dot_product(S,v)*Sin(theta)**2*B_tilde**2+&
           J_tilde*p_R*(-Sin(theta)*p_theta*dot_product(S,n)-J_tilde*p_R*dot_product(S,v))*Sin(theta)**2*B_tilde**2)+&
           omega_cos/(2*B_tilde*(Sqrt(Q)+1)*Sqrt(Q)*(Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*&
           (-p_phi**2*dot_product(S,n)*e_2mutilde*e_2nu+Sqrt(e_2nu*e_2mutilde)*J_tilde*p_R*p_phi*dot_product(S,xi)*B_tilde+&
           (dot_product(S,n)*(p_theta*Sin(theta))**2+J_tilde*p_R*dot_product(S,v)*Sin(theta)*p_theta-&
           e_2mutilde*(1+Sqrt(Q))*Sqrt(Q)*dot_product(S,n)*Sin(theta)**2)*B_tilde**2)
           
    !************************************************************************************************************
    H_SO_SS = H_SO + H_SS

    return
    end function

  !————————————————————————Spin_to_itself————————————————————————————————

  function Spin_to_itself(parameters,coordinates,Momentum)
    implicit none
    COMMON /group1/ m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2

    real * 16, INTENT(IN) :: parameters(9),coordinates(3),Momentum(3)
    real * 16 :: m_1,m_2,M,mu,nu,a,chi_Kerr,K,omega_1,omega_2
    real * 16 :: S_Star_(3),X(3),R,theta,phi,norm_S_Star
    real * 16 :: Spin_to_itself

    R = coordinates(1)
    theta = coordinates(2)
    phi = coordinates(3)

    X(:) = R*(/Sin(theta)*Cos(phi),Sin(theta)*Sin(phi),Cos(theta)/)
    S_Star_(:) = S_Star(parameters,coordinates,Momentum)
    norm_S_Star = Sqrt(S_Star_(1)**2+S_Star_(2)**2+S_Star_(3)**2)
    Spin_to_itself = -(R**2*norm_S_Star**2 - 3*dot_product(X,S_Star_)**2)/(2*M*R**5)

    return
    end function

  function numerical_differentiation(choice,coordinates,Momentum)

!     用差分代替偏微分
!     choice == 1  ===> (H_SO_SS)'_R                choice ==2   ===> (Spin_to_itself)'_R
!     choice == 3  ===> (H_SO_SS)'_PR               choice ==4   ===> (Spin_to_itself)'_PR
!     choice == 5  ===> (H_SO_SS)'_theta            choice ==6   ===> (Spin_to_itself)'_theta
!     choice == 7  ===> (H_SO_SS)'_ptheta           choice ==8   ===> (Spin_to_itself)'_ptheta
!     choice == 9  ===> (H_SO_SS)'_pphi             choice ==10  ===> (Spin_to_itself)'_pphi

    implicit none
    real * 16 , INTENT(IN)  ::  coordinates(3),Momentum(3)
    integer   , INTENT(IN)  ::  choice
    real * 16 :: f1,f2,f3,f4,h,Temp_X(3),Temp_P(3),parameters(9)
    real * 16 :: numerical_differentiation

    h = 1E-4_16 * M
    !choice == 1  ===> (H_SO_SS)'_R 
    if (choice .eq. 1) then

      Temp_X(:) = coordinates(:)

      Temp_X(1) = coordinates(1) + 2*h
      parameters(:) = Parameter_Calculator(Temp_X)
      f1 = H_SO_SS(parameters,Temp_X,Momentum)

      Temp_X(1) = coordinates(1) + h
      parameters(:) = Parameter_Calculator(Temp_X)
      f2 = H_SO_SS(parameters,Temp_X,Momentum)

      Temp_X(1) = coordinates(1) - h
      parameters(:) = Parameter_Calculator(Temp_X)
      f3 = H_SO_SS(parameters,Temp_X,Momentum)

      Temp_X(1) = coordinates(1) - 2*h 
      parameters(:) = Parameter_Calculator(Temp_X)
      f4 = H_SO_SS(parameters,Temp_X,Momentum)

    !choice ==2   ===> (Spin_to_itself)'_R
    else if (choice .eq. 2) then

      Temp_X(:) = coordinates(:)

      Temp_X(1) = coordinates(1) + 2*h
      parameters(:) = Parameter_Calculator(Temp_X)
      f1 = Spin_to_itself(parameters,Temp_X,Momentum)

      Temp_X(1) = coordinates(1) + h
      parameters(:) = Parameter_Calculator(Temp_X)
      f2 = Spin_to_itself(parameters,Temp_X,Momentum)

      Temp_X(1) = coordinates(1) - h
      parameters(:) = Parameter_Calculator(Temp_X)
      f3 = Spin_to_itself(parameters,Temp_X,Momentum)

      Temp_X(1) = coordinates(1) - 2*h
      parameters(:) = Parameter_Calculator(Temp_X)
      f4 = Spin_to_itself(parameters,Temp_X,Momentum)

    !choice == 3  ===> (H_SO_SS)'_PR 
    else if (choice .eq. 3) then

      Temp_P(:) = Momentum(:)
      parameters(:) = Parameter_Calculator(coordinates)

      Temp_P(1) = Momentum(1) + 2*h
      f1 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(1) = Momentum(1) + h
      f2 = H_SO_SS(parameters,coordinates,Temp_P)


      Temp_P(1) = Momentum(1) - h
      f3 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(1) = Momentum(1) - 2*h 
      f4 = H_SO_SS(parameters,coordinates,Temp_P)

    !choice ==4   ===> (Spin_to_itself)'_PR
    else if (choice .eq. 4 ) then

      Temp_P(:) = Momentum(:)
      parameters(:) = Parameter_Calculator(coordinates)

      Temp_P(1) = Momentum(1) + 2*h
      f1 = Spin_to_itself(parameters,coordinates,Temp_P)

      Temp_P(1) = Momentum(1) + h
      f2 = Spin_to_itself(parameters,coordinates,Temp_P)

      Temp_P(1) = Momentum(1) - h
      f3 = Spin_to_itself(parameters,coordinates,Temp_P)
      
      Temp_P(1) = Momentum(1) - 2*h
      f4 = Spin_to_itself(parameters,coordinates,Temp_P)

    !choice == 5  ===> (H_SO_SS)'_theta 
    else if (choice .eq. 5 ) then

      Temp_X(:) = coordinates(:)

      Temp_X(2) = coordinates(2) + 2*h
      parameters(:) = Parameter_Calculator(Temp_X)
      f1 = H_SO_SS(parameters,Temp_X,Momentum)

      Temp_X(2) = coordinates(2) + h
      parameters(:) = Parameter_Calculator(Temp_X)
      f2 = H_SO_SS(parameters,Temp_X,Momentum)

      Temp_X(2) = coordinates(2) - h
      parameters(:) = Parameter_Calculator(Temp_X)
      f3 = H_SO_SS(parameters,Temp_X,Momentum)

      Temp_X(2) = coordinates(2) - 2*h 
      parameters(:) = Parameter_Calculator(Temp_X)
      f4 = H_SO_SS(parameters,Temp_X,Momentum)

    !choice ==6   ===> (Spin_to_itself)'_theta
    else if (choice .eq. 6) then

      Temp_X(:) = coordinates(:)
      parameters(:) = Parameter_Calculator(Temp_X)

      Temp_X(2) = coordinates(2) + 2*h
      parameters(:) = Parameter_Calculator(Temp_X)
      f1 = Spin_to_itself(parameters,Temp_X,Momentum)

      Temp_X(2) = coordinates(2) + h
      parameters(:) = Parameter_Calculator(Temp_X)
      f2 =  Spin_to_itself(parameters,Temp_X,Momentum)

      Temp_X(2) = coordinates(2) - h
      parameters(:) = Parameter_Calculator(Temp_X)
      f3 =  Spin_to_itself(parameters,Temp_X,Momentum)

      Temp_X(2) = coordinates(2) - 2*h 
      parameters(:) = Parameter_Calculator(Temp_X)
      f4 =  Spin_to_itself(parameters,Temp_X,Momentum)

    !choice == 7  ===> (H_SO_SS)'_ptheta  
    else if (choice .eq. 7) then

      Temp_P(:) = Momentum(:)
      parameters(:) = Parameter_Calculator(coordinates)

      Temp_P(2) = Momentum(2) + 2*h
      f1 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(2) = Momentum(2) + h
      f2 = H_SO_SS(parameters,coordinates,Temp_P)


      Temp_P(2) = Momentum(2) - h
      f3 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(2) = Momentum(2) - 2*h 
      f4 = H_SO_SS(parameters,coordinates,Temp_P)

    !choice ==8   ===> (Spin_to_itself)'_ptheta
    else if (choice .eq. 8) then

      Temp_P(:) = Momentum(:)
      parameters(:) = Parameter_Calculator(coordinates)

      Temp_P(2) = Momentum(2) + 2*h
      f1 = Spin_to_itself(parameters,coordinates,Temp_P)


      Temp_P(2) = Momentum(2) + h
      f2 = Spin_to_itself(parameters,coordinates,Temp_P)


      Temp_P(2) = Momentum(2) - h
      f3 = Spin_to_itself(parameters,coordinates,Temp_P)

      Temp_P(2) = Momentum(2) - 2*h
      f4 = Spin_to_itself(parameters,coordinates,Temp_P)


    !choice == 9  ===> (H_SO_SS)'_pphi 
    else if (choice .eq. 9) then

      Temp_P(:) = Momentum(:)
      parameters(:) = Parameter_Calculator(coordinates)

      Temp_P(3) = Momentum(3) + 2*h
      f1 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(3) = Momentum(3) + h
      f2 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(3) = Momentum(3) - h
      f3 = H_SO_SS(parameters,coordinates,Temp_P)

      Temp_P(3) = Momentum(3) - 2*h 
      f4 = H_SO_SS(parameters,coordinates,Temp_P)

    !choice ==10  ===> (Spin_to_itself)'_pphi
    else if (choice .eq. 10) then

      Temp_P(:) = Momentum(:)
      parameters(:) = Parameter_Calculator(coordinates)

      Temp_P(3) = Momentum(3) + 2*h
      f1 = Spin_to_itself(parameters,coordinates,Temp_P)


      Temp_P(3) = Momentum(3) + h
      f2 = Spin_to_itself(parameters,coordinates,Temp_P)


      Temp_P(3) = Momentum(3) - h
      f3 = Spin_to_itself(parameters,coordinates,Temp_P)

      Temp_P(3) = Momentum(3) - 2*h
      f4 = Spin_to_itself(parameters,coordinates,Temp_P)


    end if

    numerical_differentiation = (-f1+8*f2-8*f3+f4)/(12*h)

    return
    end function

end program