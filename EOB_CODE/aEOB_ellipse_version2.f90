program EOB
  implicit none
  real * 8 , parameter :: Pi=3.1415926535897932384626433832795028841971693993751D0
  COMMON /group1/ M,mu,nu,a,chi_Kerr,K
  COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star
  COMMON /group4/ Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
  COMMON /group5/ gtt,grr,gthth,gphph,gtph,alpha,beta,gamma
  COMMON /group6/ dDeltar_dr,dDeltat_dr,dLambdat_dr
  COMMON /group7/ D_Parameter
  COMMON /group8/ Energy,dH_dHeff,Carter,Lz 
  real * 8 :: M,mu,nu,a,chi_Kerr,K
  real * 8 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  real * 8 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3)
  real * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
  real * 8 :: gtt,grr,gthth,gphph,gtph,alpha,beta,gamma
  real * 8 :: dDeltar_dr,dDeltat_dr,dLambdat_dr
  real * 8 :: D_Parameter(10)
  real * 8 :: Energy,dH_dHeff,Carter,Lz

  real * 8 :: m_1,m_2,K0,coordinates(3),Momentum(3),u,p,e,theta_min,theta_max,R_min,R_max,R,theta,phi,psi,zeta
  real * 8 :: Rotation(3,3),S_1(3),S_2(3),S(3),norm_S_Star
  real * 8 :: Q_4,The_guy_with_alpha
  real * 8 :: Difference=1.D0,epsilon = 1E-16
  real * 8 :: a1,a2,b1,b2,c1,c2,Z1,Z2
  real * 8 :: Time_Begin,Time_End,Time
  integer :: i = 1
  real * 8 :: K1(5),K2(5),K3(5),K4(5),h
  real * 8 :: Temp_X(3),Temp_P(3)
  real * 8 :: t,t_f
  Call CPU_TIME(Time_Begin)
  open (10,file='./1.txt',status='unknown')
  open (20,file='./2.txt',status='unknown')
  write(20,*)"t ","Energy ","Carter "

  100 format(6f25.10)
  200 format(1f10.4,2f25.15)
  m_1 = 1D0
  m_2 = 1D0
  M = m_1 + m_2
  mu = m_1*m_2/M
  nu = mu/M                                                 !文献式(39)
  S_1(:) = (/0.0D0,0.0D0,0.5D0/) * m_1 **2
  S_2(:) = (/0.D0,0.D0,0.D0/) * m_2 **2
  S_Kerr = S_1 + S_2                                              !最原始的坐标系中的自旋分量
  sigma_Star = m_1/m_2*S_2 + m_2/m_1*S_1                          !最原始的坐标系中的自旋分量
  norm_S_Kerr = Sqrt(S_Kerr(1)**2+S_Kerr(2)**2+S_Kerr(3)**2)
  if ((Abs(S_Kerr(1)) .ge. 1E-10) .or. (Abs(S_Kerr(2)) .ge. 1E-10)) then
    psi  = Acos(S_Kerr(3)/norm_S_Kerr)
    zeta = Atan(S_Kerr(2)/S_Kerr(1))
    Rotation(1,:) = (/-Sin(zeta),Cos(zeta),0.D0/)
    Rotation(2,:) = (/-Cos(psi)*Cos(zeta),-Cos(psi)*Sin(zeta),Sin(psi)/)
    Rotation(3,:) = (/Sin(psi)*Cos(zeta),Sin(psi)*Sin(zeta),Cos(psi)/)
    S_Kerr(:) = matmul(Rotation,S_Kerr)                                    !旋转坐标系后的自旋分量
    sigma_Star(:) = matmul(Rotation,sigma_Star)                            !旋转坐标系后的自旋分量
  end if
  a = norm_S_Kerr/M
  chi_Kerr = a/M
  K0 = 1.4467
  K = K0*(1-4*nu)**2 + 4*(1-2*nu)*nu                            !文献式(37)
  Delta_0 = K*(nu*K-2)
  Delta_1 = -2*(nu*K-1)*(K+Delta_0)
  Delta_2 = 1/2.*Delta_1*(-4 * nu * K + Delta_1 + 4) - chi_Kerr**2*(nu*K - 1)**2*Delta_0
  Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
                  (-nu*K+Delta_2 + 1)-3*chi_Kerr**2*(nu*K-1)**2*Delta_1)
  Delta_4 = 1/12.*(6*chi_Kerr**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2+3*Delta_1**4-8*(nu*K-1)*&
                Delta_1**3 - 12*Delta_2*Delta_1**2+12*(2*(nu*K-1)*Delta_2+Delta_3)*Delta_1 +&
                12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2+6*(Delta_2**2-4*Delta_3*(nu*K-1)))
  p = 5 * M 
  e = 0.2
  R_min = p/(1+e)
  R_max = p/(1-e)
  u = M/R_min
  theta_min = Pi/4.
  theta_max = Pi-theta_min
  coordinates(:) = (/R_min,theta_min,0.D0/)
  Momentum(:) = (/0.D0,0.D0,0.D0/)
  call Parameter_Calculator((/R_min,theta_min,0D0/))
  call Metric(theta_min)
  a1 = alpha
  b1 = beta
  c1 = gamma
  call Parameter_Calculator((/R_max,theta_max,Pi/))
  call Metric(theta_min)
  a2 = alpha
  b2 = beta
  c2 = gamma

  Momentum(1) = 0.D0
  Momentum(2) = 0.D0
  Momentum(3) = Sqrt((b1**2*(a1**2+a2**2)-2*b1*b2*(a1**2+a2**2)+b2**2*(a1**2+a2**2)-a1**4*c1+a1**2*a2**2*c1-&
                2*Sqrt((b1-b2)**2*a1**2*a2**2*((b1-b2)**2-(a1**2-a2**2)*(c1-c2)))+a1**2*a2**2*c2-a2**4*c2)/&
                (b1**4-4*b1**3*b2+b2**4+(a1**2*c1-a2**2*c2)**2-2*b2**2*(a1**2*c1+a2**2*c2)+4*b1*b2*(-b2**2+&
                  a1**2*c1+a2**2*c2)+b1**2*(6*b2**2-2*(a1**2*c1+a2**2*c2))))
  Energy = b1*Momentum(3) + a1*Sqrt(1+c1*Momentum(3)**2)
  Carter = Cos(theta_min)**2*((1-Energy**2)*a**2+Momentum(3)**2/Sin(theta_min)**2)
  write(*,*)"The Angular Momentum L equals to:",Momentum(3)
  write(*,*)"The Effective Energy E equals to:",Energy
  write(*,*)"The Carter Constant Q equals to:",Carter

!   !——————————————————————————————————————————————初始条件求解完毕————————————————————————————————————————————————
  do while(Difference .ge. epsilon)
    Lz = Momentum(3)
    call Parameter_Calculator((/R_min,theta_min,0.D0/))
    call Metric(theta_min)
    Z1 = H_S((/R_min,theta_min,0D0/),Momentum)
    call Parameter_Calculator((/R_max,theta_max,Pi/))
    call Metric(theta_max)
    Z2 = H_S((/R_max,theta_max,Pi/),Momentum)
    Momentum(3) = Sqrt(((((b2-b1)*Lz+Z2-Z1+a2*Sqrt(1+c2*Lz**2))/a1)**2-1)/c1)
    Difference = abs(Momentum(3)-Lz)
    i = i + 1
    write(*,*)"i=",i,"Difference=",Difference,"Momentum=",Momentum(3)
  end do
  call CPU_TIME(Time_End)
  Time = Time_End - Time_Begin
  write(*,*) Time,"s"
!   h = 0.01 * M
!   t = 0
!   t_f = 10000 * M
!   coordinates(:) = (/R_min,theta_min,0.D0/)
!   call Evolution_RK4(h,coordinates,Momentum,t,t_f)
!   write(*,*) "程序运行成功，可以查看分析数据了哦～"
!   close(10)
!   close(20)
!   Call CPU_TIME(Time_End)
!   Time = Time_End - Time_Begin
!   write(*,*)Time

  contains
  Subroutine Parameter_Calculator(coordinates)
      implicit none
      COMMON /group1/ M,mu,nu,a,chi_Kerr,K
      COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
      COMMON /group4/ Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
      COMMON /group6/ dDeltar_dr,dDeltat_dr,dLambdat_dr
      real * 8 :: M,mu,nu,a,chi_Kerr,K
      real * 8 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
      real * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
      real * 8 :: dDeltar_dr,dDeltat_dr,dLambdat_dr
      real * 8 :: coordinates(3)
      real * 8 :: R,theta,u,Delta_Bar_u,Delta_u,D_Minus_u,dDeltau_du,dD_Minus_u_du
      R = coordinates(1)
      theta = coordinates(2)
      u = M/R
      Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2
      Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))
      Delta_t = R**2*Delta_u
      varpi_square = R**2 + a**2
      D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)
      Delta_r = Delta_t * D_Minus_u
      omega_tilde = 2*a*M*R
      Sigma = R**2+a**2*Cos(theta)**2
      Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2
      dDeltau_du =(nu*((-1+K*nu)**(-2)+(2*u)/(-1+K*nu)+chi_Kerr**2*u**2)*(Delta_1+u*(2*Delta_2+u*(3*Delta_3+4*Delta_4*u))))/ &
                (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u))))+2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
                (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))
      dDeltat_dr = 2*R*Delta_u - dDeltau_du*M
      dD_Minus_u_du = (12*u*nu+6*u**2*(26-3*nu)*nu)/(1 + 6*u**2*nu + 2*nu*u**3*(26-3*nu))
      dDeltar_dr = dDeltat_dr * D_Minus_u - Delta_t * dD_Minus_u_du * M /R**2
      dLambdat_dr = 4*R*varpi_square - a**2*Sin(theta)**2*dDeltat_dr
      return
      end Subroutine
  Subroutine Metric(theta)
      implicit none
      COMMON /group4/ Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
      COMMON /group5/ gtt,grr,gthth,gphph,gtph,alpha,beta,gamma
      real * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
      real * 8 :: gtt,grr,gthth,gphph,gtph,alpha,beta,gamma
      real * 8 :: theta

      gtt = -Lambda_t/(Delta_t*Sigma)
      grr = Delta_r/Sigma
      gthth = 1/Sigma
      gphph = 1/Lambda_t*(-omega_tilde**2/(Delta_t*Sigma)+Sigma/Sin(theta)**2)
      gtph = -omega_tilde/(Delta_t*Sigma)
      alpha = 1/Sqrt(-gtt)
      beta = gtph/gtt
      gamma = gphph - gtph**2/gtt
      return
      end Subroutine
!   Subroutine Derivative_Of_Parameter(coordinates)
!     implicit none
!     !dgtt_dr,dgtt_dth,dgrr_dr,dgrr_dth,dgthth_dr,dgthth_dth,dgphph_dr,dgphph_dth,dgtph_dr,dgtph_dth
!     !dalpha_dr,dalpha_dth,dbeta_dr,dbeta_dth,dgammarr_dr,dgammarr_dth,dgammatt_dr,dgammatt_dth,dgammapp_dr,dgammapp_dth
!     COMMON /group1/ M,mu,nu,a
!     COMMON /group4/ Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
!     COMMON /group5/ gtt,grr,gthth,gphph,gtph,alpha,beta,gamma
!     COMMON /group6/ dDeltar_dr,dDeltat_dr,dLambdat_dr
!     COMMON /group7/ D_Parameter
!     real * 8 :: M,mu,nu,a
!     real * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
!     real * 8 :: gtt,grr,gthth,gphph,gtph,alpha,beta,gamma
!     real * 8 :: dDeltar_dr,dDeltat_dr,dLambdat_dr
!     real * 8 :: D_Parameter(10)
!     real * 8 :: R,theta,D_Metric(10),coordinates(3)
!     R = coordinates(1)
!     theta = coordinates(2)
!     D_Metric(1) = (Lambda_t*(Sigma*dDeltat_dr + 2*R*Delta_t) - dLambdat_dr * Delta_t * Sigma)/(Delta_t*Sigma)**2
!     D_Metric(3) = (dDeltar_dr*Sigma - 2*R*Delta_r)/Sigma**2
!     D_Metric(5) = -2*R/Sigma**2
!     D_Metric(7) = (2*R*Lambda_t-dLambdat_dr*Sigma)/(Lambda_t**2*Sin(theta)**2)+ &
!                   (-8*a**2*M**2*R*Lambda_t*Delta_t*Sigma+omega_tilde**2*(dLambdat_dr*Delta_t*Sigma+ &
!                   Lambda_t*dDeltat_dr*Sigma+2*R*Lambda_t*Delta_t))/(Lambda_t*Sigma*Delta_t)**2
!     D_Metric(9) = (omega_tilde*(dDeltat_dr * Sigma + 2 * R * Delta_t)-2*a*M*Sigma*Delta_t)/(Delta_t*Sigma)**2
!     D_Metric(2) =  a**2*Sin(2*theta)*varpi_square*(Delta_t-varpi_square)/(Delta_t*Sigma**2)
!     D_Metric(4) = Delta_r*a**2*Sin(2*theta)/Sigma**2
!     D_Metric(6) = a**2*Sin(2*theta)/Sigma**2
!     D_Metric(8) = a**2*Sin(2*theta)/Lambda_t**2*(-Lambda_t/Sin(theta)**2+&
!                   Sigma*Delta_t/Sin(theta)**2-Lambda_t*Sigma/(a**2*Sin(theta)**4)-&
!                   omega_tilde**2*(Sigma*Delta_t+Lambda_t)/(Sigma**2*Delta_t))
!     D_Metric(10) = -a**2*omega_tilde*Sin(2*theta)/(Sigma**2*Delta_t)
!     D_Parameter(1) = 1/(2.*(-gtt)**1.5d0)*D_Metric(1)
!     D_Parameter(2) = 1/(2.*(-gtt)**1.5d0)*D_Metric(2)
!     D_Parameter(3) = (2*a*M*Lambda_t-dLambdat_dr*2*a*M*R)/Lambda_t**2
!     D_Parameter(4) = omega_tilde*a**2*Sin(2*theta)*Delta_t/Lambda_t**2
!     D_Parameter(5) = D_Metric(3)
!     D_Parameter(6) = D_Metric(4)
!     D_Parameter(7) = D_Metric(5)
!     D_Parameter(8) = D_Metric(6)
!     D_Parameter(9) = D_Metric(7) - (2*gtph*gtt*D_Metric(9)-gtph**2*D_Metric(1))/gtt**2
!     D_Parameter(10) = D_Metric(8) - (2*gtph*gtt*D_Metric(10)-gtph**2*D_Metric(2))/gtt**2
!     return
!     end Subroutine

  function H_S(coordinates,Momentum)
    Implicit None
    COMMON /group1/ M,mu,nu,a
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star
    COMMON /group4/ Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    COMMON /group6/ dDeltar_dr,dDeltat_dr,dLambdat_dr
    real * 8 :: M,mu,nu,a
    real * 8 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3)
    real * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    real * 8 :: dDeltar_dr,dDeltat_dr,dLambdat_dr
    real * 8 :: coordinates(3),Momentum(3)
    real * 8 :: d_SO,d_SS,R,theta,phi,p_R,p_theta,p_phi,u,n(3),grr
    real * 8 :: omega,omega_r,omega_cos,e_2nu,e_2mutilde,B_rtilde,B_tilde,J_tilde,mu_r,mu_cos,nu_r,nu_cos,Q
    real * 8 :: Delta_sigma_1(3),Delta_sigma_2(3),Delta_sigma_3(3),S_Star(3),S(3),xi(3),v(3),X(3),norm_S_Star
    real * 8 :: H_SO,H_SS,Spin_to_itself,H_S
    d_SO = -69.5D0                                              
    d_SS = 2.75D0 
    R = coordinates(1)
    u = M/R
    theta = coordinates(2)
    phi = coordinates(3)
    p_R = Momentum(1)
    p_theta = Momentum(2)
    p_phi = Momentum(3)
    grr = Delta_r/Sigma
    Q  = 1 + Delta_r*p_R**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma 
    Delta_sigma_1(:) = sigma_Star(:) * (7/6.*nu*u + nu/3.*(Q-1)-5/2.*nu*grr*p_R**2)+&
                       S_Kerr(:) * (-2/3.*nu*u+1/4.*nu*(Q-1)-3*nu*grr*p_R**2)
    Delta_sigma_2(:) = sigma_Star(:) * (1/36.*(353*nu-27*nu**2)*u**2+5*nu**2*grr**2*p_R**4-&
                       1/72.*(23*nu+3*nu**2)*(Q-1)**2+1/36.*(-103*nu+60*nu**2)*u*(Q-1)+&
                       1/12.*(16*nu-21*nu**2)*grr*p_R**2*(Q-1)+1/12.*(47*nu-54*nu**2)*u*grr*p_R**2)+&
                       S_Kerr(:) * (1/9.*(-56*nu-21*nu**2)*u**2+45/8.*nu**2*grr**2*p_R**4-&
                        5/16.*nu*(Q-1)**2+1/36.*(-109*nu+51*nu**2)*u*(Q-1)+&
                        1/8.*(2*nu-13*nu**2)*grr*p_R**2*(Q-1)-1/24.*(16*nu+147*nu**2)*u*grr*p_R**2)
    Delta_sigma_3(:) = d_SO*nu/R**3*sigma_Star(:)
    S_Star(:) = sigma_Star(:) + Delta_sigma_1(:) + Delta_sigma_2(:) + Delta_sigma_3(:)     
    X(:) = R*(/Sin(theta)*Cos(phi),Sin(theta)*Sin(phi),Cos(theta)/)               
    norm_S_Star = Sqrt(S_Star(1)**2+S_Star(2)**2+S_Star(3)**2)      
    Spin_to_itself = -(R**2*norm_S_Star**2 - 3*dot_product(X,S_Star)**2)/(2*M*R**5)                                                                     
    S(:) = nu * S_Star(:) 
    n(:) = (/Sin(theta)*Cos(phi),Sin(theta)*Sin(phi),Cos(theta)/) 
    omega = omega_tilde / Lambda_t                                                                             
    omega_r = (-dLambdat_dr*omega_tilde+Lambda_t*2*a*M)/Lambda_t**2                                          
    omega_cos = -2*a**2*Cos(theta)*Delta_t*omega_tilde/Lambda_t**2                                              
    e_2nu = Delta_t*Sigma/Lambda_t                                                                             
    e_2mutilde = Sigma                                                                                         
    B_tilde = Sqrt(Delta_t)                                                                                    
    B_rtilde = (Sqrt(Delta_r)*dDeltat_dr-2*Delta_t)/(2*Sqrt(Delta_t*Delta_r))                               
    J_tilde = Sqrt(Delta_r)                                                                                    
    mu_r = R/Sigma - 1/Sqrt(Delta_r)                                                                           
    mu_cos = a**2*Cos(theta)/Sigma                                                                             
    nu_r = R/Sigma + varpi_square*(varpi_square*dDeltat_dr-4*R*Delta_t)/(2*Lambda_t*Delta_t)                
    nu_cos = a**2*varpi_square*Cos(theta)*(varpi_square-Delta_t)/(Lambda_t*Sigma)                                                                 
    xi(:) = Sin(theta)*(/-Sin(phi),Cos(phi),0.D0/)                                                            
    v(:) = -Sin(theta)*(/Cos(theta)*Cos(phi),Cos(theta)*Sin(phi),-Sin(theta)/)          
    H_SO = e_2nu/Sqrt(e_2mutilde)*(Sqrt(e_2mutilde*e_2nu)-B_tilde)*p_phi*&
           dot_product(S,S_Kerr/norm_S_Kerr)/(B_tilde**2*Sqrt(Q)*Sin(theta)**2)+&
           Sqrt(e_2nu)/(e_2mutilde*B_tilde**2*(Sqrt(Q)+1)*Sqrt(Q)*Sin(theta)**2)*&
           (dot_product(S,xi)*J_tilde*(-mu_r*Sin(theta)*p_theta*(Sqrt(Q)+1)-mu_cos*p_R*Sin(theta)**2-&
           Sqrt(Q)*(-nu_r*Sin(theta)*p_theta+(mu_cos-nu_cos)*p_R*Sin(theta)**2))*B_tilde**2 +&
           Sqrt(e_2mutilde*e_2nu)*p_phi*(2*Sqrt(Q)+1)*B_tilde*(J_tilde*nu_r*dot_product(S,v)-nu_cos*dot_product(S,n)*&
           Sin(theta)**2)-J_tilde*B_rtilde*Sqrt(e_2mutilde*e_2nu)*p_phi*(Sqrt(Q)+1)*dot_product(S,v))
    H_SS = omega*dot_product(S,S_Kerr/norm_S_Kerr) + J_tilde*omega_r/(2*B_tilde*(Sqrt(Q)+1)*Sqrt(Q)*Sin(theta)**2*&
           (Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*(Sqrt(e_2mutilde*e_2nu)*Sin(theta)*p_theta*p_phi*dot_product(S,xi)*B_tilde+&
           e_2nu*e_2mutilde*p_phi**2*dot_product(S,v)+e_2mutilde*(1+Sqrt(Q))*Sqrt(Q)*dot_product(S,v)*Sin(theta)**2*B_tilde**2+&
           J_tilde*p_R*(-Sin(theta)*p_theta*dot_product(S,n)-J_tilde*p_R*dot_product(S,v))*Sin(theta)**2*B_tilde**2)+&
           omega_cos/(2*B_tilde*(Sqrt(Q)+1)*Sqrt(Q)*(Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*&
           (-p_phi**2*dot_product(S,n)*e_2mutilde*e_2nu+Sqrt(e_2nu*e_2mutilde)*J_tilde*p_R*p_phi*dot_product(S,xi)*B_tilde+&
           (dot_product(S,n)*(p_theta*Sin(theta))**2+J_tilde*p_R*dot_product(S,v)*Sin(theta)*p_theta-&
           e_2mutilde*(1+Sqrt(Q))*Sqrt(Q)*dot_product(S,n)*Sin(theta)**2)*B_tilde**2)
    H_S = (H_SS + H_SO)/mu + Spin_to_itself + d_SS*nu*dot_product(S_Kerr,sigma_Star)/R**4
    return
    end function
!   function P_Dot(choice,coordinates,Momentum)
!     implicit none
!     COMMON /group1/ M,mu,nu,a
!     COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star,d_SO,d_SS
!     COMMON /group4/ Sigma,Lambda_t,omega_tilde,Delta_t,Delta_r
!     COMMON /group5/ grr,gthth,gphph,gtph,alpha,beta,gamma
!     COMMON /group7/ D_Parameter
!     COMMON /group8/ The_guy_with_alpha,Energy,dH_dHeff
!     real * 8 :: M,mu,nu,a
!     real * 8 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3),d_SO,d_SS
!     real * 8 :: Sigma,Lambda_t,omega_tilde,Delta_t,Delta_r
!     real * 8 :: grr,gthth,gphph,gtph,alpha,beta,gamma
!     real * 8 :: D_Parameter(10)
!     real * 8 :: Energy,dH_dHeff
!     real * 8 :: coordinates(3),Momentum(3)
!     real * 8 :: Q_4,The_guy_with_alpha,dQ_dr,dQ_dth,h
!     real * 8 :: Temp_X(3),f1,f2,dHS_dr,dHS_dth
!     real * 8 :: P_Dot
!     integer :: choice
!     h = 1E-4 * M
!     if (choice .eq. 1)then
!       dQ_dr = -4*(4-3*nu)*nu*M**2*coordinates(1)*Momentum(1)**4/Sigma**2
!       Temp_X(1) = coordinates(1) + h/2D0
!       Temp_X(2:3) = coordinates(2:3)
!       call Parameter_Calculator(Temp_X)
!       f1 = H_S(Temp_X,Momentum)
!       Temp_X(1) = coordinates(1) - h/2D0
!       call Parameter_Calculator(Temp_X)
!       f2 = H_S(Temp_X,Momentum)
!       dHS_dr = (f1-f2)/h
!       P_Dot = -dH_dHeff*(D_parameter(3)*Momentum(3)+D_parameter(1)*The_guy_with_alpha+alpha*(D_parameter(5)*Momentum(1)**2+&
!                 D_parameter(7)*Momentum(2)**2+D_parameter(9)*Momentum(3)**2 + dQ_dr)/(2*The_guy_with_alpha) + dHS_dr)
!     else if (choice .eq. 2)then
!       dQ_dth = 2*(4-3*nu)*nu*M**2*a**2*Sin(2*coordinates(2))*Momentum(1)**4/Sigma**2
!       Temp_X(:) = coordinates(:)
!       Temp_X(2) = coordinates(2) + h/2D0
!       call Parameter_Calculator(Temp_X)
!       f1 = H_S(Temp_X,Momentum)
!       Temp_X(2) = coordinates(2) - h/2D0
!       call Parameter_Calculator(Temp_X)
!       f2 = H_S(Temp_X,Momentum)
!       dHS_dth = (f1-f2)/h
!       P_Dot = -dH_dHeff*(D_parameter(4)*Momentum(3)+D_parameter(2)*The_guy_with_alpha+alpha*(D_parameter(6)*Momentum(1)**2+&
!                   D_parameter(8)*Momentum(2)**2+D_parameter(10)*Momentum(3)**2+dQ_dth)/(2*The_guy_with_alpha) + dHS_dth)
!     end if
!     end function
!   function X_Dot(choice,coordinates,Momentum)
!     implicit none
!     COMMON /group1/ M,mu,nu
!     COMMON /group4/ Sigma,Lambda_t,omega_tilde,Delta_t,Delta_r
!     COMMON /group5/ grr,gthth,gphph,gtph,alpha,beta,gamma
!     COMMON /group7/ D_Parameter
!     COMMON /group8/ Energy,dH_dHeff
!     real * 8 :: M,mu,nu,a
!     real * 8 :: Sigma,Lambda_t,omega_tilde,Delta_t,Delta_r
!     real * 8 :: grr,gthth,gphph,gtph,alpha,beta,gamma
!     real * 8 :: D_parameter(10)
!     real * 8 :: Energy,dH_dHeff
!     real * 8 :: coordinates(3),Momentum(3)
!     real * 8 :: Q_4,The_guy_with_alpha,dQ_dpr
!     real * 8 :: h,Temp_P(3),f1,f2,dHS_dpr,dHS_dpth,dHS_dpph
!     real * 8 :: X_Dot
!     integer :: choice
!     h = 1E-4 * M
!     if (choice .eq. 1) then
!       dQ_dpr = 8*(4-3*nu)*nu*M**2*Momentum(1)**3/Sigma
!       Temp_P(:) = Momentum(:)
!       call Parameter_Calculator(coordinates)
!       Temp_P(1) = Momentum(1) + h/2D0
!       f1 = H_S(coordinates,Temp_P)
!       Temp_P(1) = Momentum(1) - h/2D0 
!       f2 = H_S(coordinates,Temp_P)
!       dHS_dpr = (f1-f2)/h
!       X_Dot = dH_dHeff*(alpha*(2*grr*Momentum(1)+dQ_dpr)/(2*The_guy_with_alpha)+dH_dpr)
!     else if (choice .eq. 2)then
!       Temp_P(:) = Momentum(:)
!       call Parameter_Calculator(coordinates)
!       Temp_P(2) = Momentum(2) + h/2D0
!       f1 = H_S(coordinates,Temp_P)
!       Temp_P(2) = Momentum(2) - h/2D0 
!       f2 = H_S(coordinates,Temp_P)
!       dHS_dpth = (f1-f2)/h
!       X_Dot = dH_dHeff*(alpha*gthth*Momentum(2)/The_guy_with_alpha+dHS_dpth)
!     else if (choice .eq. 3)then
!       Temp_P(:) = Momentum(:)
!       call Parameter_Calculator(coordinates)
!       Temp_P(3) = Momentum(3) + h/2D0
!       f1 = H_S(coordinates,Temp_P)
!       Temp_P(3) = Momentum(3) - h/2D0 
!       f2 = H_S(coordinates,Temp_P)
!       dHS_dpph = (f1-f2)/h
!       X_Dot = dH_dHeff*(beta + alpha*gamma*Momentum(3)/The_guy_with_alpha+dHS_dpph)
!     end if 
!     return
!     end function
!   Subroutine Evolution_RK4(h,coordinates,Momentum,t,t_f)
!     implicit none
!     COMMON /group1/ M,mu,nu,a
!     COMMON /group8/ Q4,The_guy_with_alpha,Energy,dH_dHeff,Carter,Lz
!     real * 8 :: M,mu,nu,a
!     real * 8 :: h,coordinates(3),Momentum(3),t,t_f
!     real * 8 :: X(5),X1(5),g(8)
!     real * 8 :: Y0(5),Y1(5),Y2(5),Y3(5)
!     real * 8 :: Energy,Carter,parameters(9)
!     integer :: i
!     100 format(6f25.10)
!     200 format(1f10.4,2f25.15)

!     do while (t .le. t_f)
!       X(1:3) = coordinates(:)
!       X(4:5) = Momentum(1:2)
!       X1(:) = X(:)
!       call Differential_Eq(X1,Y0)
!       X1(:) = X(:) + h*Y0(:)/2D0
!       call Differential_Eq(X1,Y1)
!       X1(:) = X(:) + h*Y1(:)/2D0
!       call Differential_Eq(X1,Y2)
!       X1(:) = X(:) + h*Y2(:)
!       call Differential_Eq(X1,Y3)
!       X(:) = X(:) + h*(Y0(:) + 2D0*Y1(:) + 2D0*Y2(:) + Y3(:))/6D0
!       t = t + h
!       coordinates(:) = X(:)
!       Momentum(1:2) = X(4:5)
!       call Parameter_Calculator(coordinates)
!       call Derivative_of_parameter(Temp_X)
!       call Metric(Temp_X(2))
!       g(:) = Metric(coordinates(2))
!       Q_4 = 2*(4-3*nu)*nu*M**2/Sigma*Momentum(1)**4
!       The_guy_with_alpha = Sqrt(1+grr*Momentum(1)**2+gthth*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)
!       Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + H_S(coordinates,Momentum)
!       dH_dHeff = 1/Sqrt(1+2*nu*(Energy-1))
!       Carter = Momentum(2)**2 + Cos(coordinates(2))**2*((1-Energy**2)*a**2 + Momentum(3)**2/Sin(coordinates(2))**2)
!       write(10,100)coordinates,Momentum
!       write(20,200)t,Energy,Carter
!     end do
!     return
!     end Subroutine
!   Subroutine Differential_Eq(X1,Y)
!     implicit none
!     COMMON /group1/ M,mu,nu,a
!     COMMON /group5/ grr,gthth,gphph,gtph,alpha,beta,gamma
!     COMMON /group8/ Q4,The_guy_with_alpha,Energy,dH_dHeff,Carter,Lz
!     real * 8 :: M,mu,nu,a
!     real * 8 :: grr,gthth,gphph,gtph,alpha,beta,gamma
!     real * 8 :: Q4,The_guy_with_alpha,Energy,dH_dHeff,Carter,Lz 
!     real * 8 :: X1(5),Y(5),Temp_X(3),Temp_P(3)
!     Temp_X(:) = X1(1:3)
!     Temp_P(:) = (/X1(4),X1(5),Lz/)
!     call Parameter_Calculator(Temp_X)
!     call Derivative_of_parameter(Temp_X)
!     call Metric(Temp_X(2))
!     Q_4 = 2*(4-3*nu)*nu*M**2/Sigma*Momentum(1)**4
!     The_guy_with_alpha = Sqrt(1+grr*Momentum(1)**2+gthth*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)
!     Energy = beta*Momentum(3) + alpha*The_guy_with_alpha + H_S(coordinates,Momentum)
!     dH_dHeff = 1/Sqrt(1+2*nu*(Energy-1))
!     Y(1) = X_Dot(1,Temp_X,Temp_P)                    
!     Y(2) = X_Dot(2,Temp_X,Temp_P)                       
!     Y(3) = X_Dot(3,Temp_X,Temp_P)                          
!     Y(4) = P_Dot(1,Temp_X,Temp_P) 
!     Y(5) = P_Dot(2,Temp_X,Temp_P) 
!     return 
!     end Subroutine

end program