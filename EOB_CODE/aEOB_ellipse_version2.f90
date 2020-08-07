program EOB
  IMPLICIT NONE
  REAL * 8 , parameter :: Pi=3.141592653589793D0
  COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
  COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star
  COMMON /group4/ Lz
  REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
  REAL * 8 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
  REAL * 8 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3)
  REAL * 8 :: Lz
  REAL * 8 :: m_1,m_2,K0,coordinates(3),Momentum(3),u,p,e,theta_min,theta_max,R_min,R_max,R,theta,phi,psi,zeta
  REAL * 8 :: Rotation(3,3),S_1(3),S_2(3),S(3),norm_S_Star,Q_4
  REAL * 8 :: Difference=1.D0,epsilon = 1E-16_8
  REAL * 8 :: Para1(9),Para2(9),a1,a2,b1,b2,c1,c2,Z1,Z2,square_Q,Energy,dH_dHeff,Carter
  REAL * 8 :: Time_Begin,Time_END,Time
  REAL * 8 :: h,t,t_f
  INTEGER :: i = 1
  OPEN (10,file='./1.txt',status='unknown')
  OPEN (20,file='./2.txt',status='unknown')
  WRITE(20,*)"t ","Energy ","Carter "
  m_1 = 1.D0
  m_2 = 1E-10
  M = m_1 + m_2
  mu = m_1*m_2/M
  nu = mu/M                                                 
  S_1(:) = (/0D0,0D0,0.5D0/) * m_1 **2
  S_2(:) = (/0D0,0D0,0D0/) * m_2 **2
  S_Kerr = S_1 + S_2                                             
  sigma_Star = m_1/m_2*S_2 + m_2/m_1*S_1                         
  norm_S_Kerr = Sqrt(S_Kerr(1)**2+S_Kerr(2)**2+S_Kerr(3)**2)
  IF ((ABS(S_Kerr(1)) .GE. 1E-4_8) .OR. (ABS(S_Kerr(2)) .GE. 1E-4_8)) THEN
    psi  = ACOS(S_Kerr(3)/norm_S_Kerr)
    zeta = ATAN(S_Kerr(2)/S_Kerr(1))
    Rotation(1,:) = (/-SIN(zeta),COS(zeta),0D0/)
    Rotation(2,:) = (/-COS(psi)*COS(zeta),-COS(psi)*SIN(zeta),SIN(psi)/)
    Rotation(3,:) = (/SIN(psi)*COS(zeta),SIN(psi)*SIN(zeta),COS(psi)/)
    S_Kerr(:) = matmul(Rotation,S_Kerr)                                    !旋转坐标系后的自旋分量
    sigma_Star(:) = matmul(Rotation,sigma_Star)                            !旋转坐标系后的自旋分量
  END IF
  a = norm_S_Kerr/M
  chi_Kerr = a/M
  K0 = 1.4467D0
  K = K0*(1D0-4D0*nu)**2+4D0*(1D0-2D0*nu)*nu 
  Delta_0 = K*(nu*K-2D0)
  Delta_1 = -2D0*(nu*K-1D0)*(K+Delta_0)
  Delta_2 = 1D0/2D0*Delta_1*(-4D0*nu*K+Delta_1+4D0)-chi_Kerr**2D0*(nu*K-1D0)**2D0*Delta_0
  Delta_3 = 1D0/3D0*(-Delta_1**3D0+3D0*(nu*K-1D0)*Delta_1**2D0+3D0*Delta_1*Delta_2-6D0*(nu*K-1D0)*&
            (-nu*K+Delta_2+1D0)-3D0*chi_Kerr**2D0*(nu*K-1D0)**2D0*Delta_1)
  Delta_4 = 1D0/12D0*(6D0*chi_Kerr**2D0*(Delta_1**2D0-2D0*Delta_2)*(nu*K-1D0)**2+3D0*Delta_1**4-8D0*(nu*K-1D0)*&
            Delta_1**3-12D0*Delta_2*Delta_1**2+12D0*(2D0*(nu*K-1D0)*Delta_2+Delta_3)*Delta_1 +&
            12D0*(94D0/3D0-41D0/32D0*Pi**2D0)*(nu*K-1D0)**2+6D0*(Delta_2**2-4D0*Delta_3*(nu*K-1D0)))
  nuK_1 = 1/(nu*K-1) 
  p = 5.D0 * M 
  e = 0.2D0
  R_min = p/(1D0+e)
  R_max = p/(1D0-e)
  u = M/R_min
  theta_min = Pi/4D0
  theta_max = Pi-theta_min
  coordinates(:) = (/R_min,theta_min,0D0/)
  write(*,*)coordinates
  Momentum(:) = (/0D0,0D0,0D0/)
  Para1 =  Parameter_Calculator((/R_min,theta_min,0D0/))
  a1 = 1D0/SQRT(Para1(2)/(Para1(1)*Para1(5)))
  b1 = Para1(4)/Para1(2)
  c1 = 1D0/Para1(2)*(-Para1(4)**2D0/(Para1(5)*Para1(1))+Para1(1)/SIN(theta_min)**2D0)+Para1(4)**2D0/(Para1(1)*Para1(2)*Para1(5))
  Para2 = Parameter_Calculator((/R_max,theta_max,Pi/))
  a2 = 1D0/SQRT(Para2(2)/(Para2(1)*Para2(5)))
  b2 = Para2(4)/Para2(2)
  c2 = 1D0/Para2(2)*(-Para2(4)**2D0/(Para2(5)*Para2(1))+Para2(1)/SIN(theta_max)**2D0)+Para2(4)**2D0/(Para2(1)*Para2(2)*Para2(5))
  Momentum(1) = 0D0
  Momentum(2) = 0D0
  Momentum(3) = SQRT((b1**2D0*(a1**2D0+a2**2D0)-2D0*b1*b2*(a1**2D0+a2**2D0)+b2**2D0*(a1**2D0+a2**2D0)-&
                a1**4D0*c1+a1**2D0*a2**2D0*c1-2D0*Sqrt((b1-b2)**2D0*a1**2D0*a2**2D0*((b1-b2)**2D0-(a1**2D0-&
                a2**2D0)*(c1-c2)))+a1**2D0*a2**2D0*c2-a2**4D0*c2)/(b1**4D0-4D0*b1**3D0*b2+b2**4D0+(a1**2D0*c1-&
                a2**2D0*c2)**2D0-2D0*b2**2D0*(a1**2D0*c1+a2**2D0*c2)+4D0*b1*b2*(-b2**2D0+&
                a1**2D0*c1+a2**2D0*c2)+b1**2D0*(6D0*b2**2D0-2D0*(a1**2D0*c1+a2**2D0*c2))))
  Energy = b1*Momentum(3) + a1*SQRT(1+c1*Momentum(3)**2)
  Carter = COS(theta_min)**2*((1-Energy**2)*a**2+Momentum(3)**2/SIN(theta_min)**2)
  WRITE(*,*)"The Angular Momentum L equals to:",Momentum(3)
  WRITE(*,*)"The Effective Energy E equals to:",Energy
  WRITE(*,*)"The Carter Constant Q equals to:",Carter
  DO WHILE(Difference .GE. epsilon)
    Lz = Momentum(3)
    Z1 = H_S((/R_min,theta_min,0D0/),Momentum,Para1)
    Z2 = H_S((/R_max,theta_max,Pi/),Momentum,Para2)
    Momentum(3) = Sqrt(((((b2-b1)*Lz+Z2-Z1+a2*Sqrt(1+c2*Lz**2))/a1)**2-1)/c1)
    Difference = ABS(Momentum(3)-Lz)
    WRITE(*,*)"i=",i,"Difference=",Difference,"Momentum=",Momentum(3)
    i = i + 1
  END DO
  h = 0.1D0 * M
  t = 0D0
  t_f = 20000D0 * M
  CALL Evolution_RK78(h,coordinates,Momentum,t,t_f)
  WRITE(*,*) "程序运行成功"
  close(10)
  close(20)
  CALL CPU_TIME(Time_END)
  Time = Time_END - Time_Begin
  WRITE(*,*)Time
CONTAINS
  FUNCTION Parameter_Calculator(coordinates)
      IMPLICIT NONE
      COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
      COMMON /group2/ Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
      REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
      REAL * 8 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
      REAL * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
      REAL * 8 :: dDeltar_dr,dDeltat_dr,dLambdat_dr
      REAL * 8 :: coordinates(3)
      REAL * 8 :: R,theta,u,Delta_Bar_u,Delta_u,D_Minus_u,dDeltau_du,dD_Minus_u_du
      REAL * 8 :: Long_Delta1,Long_Delta2
      REAL * 8 :: Parameter_Calculator(9)
      R = coordinates(1)
      theta = coordinates(2)
      u = M/R
      Delta_Bar_u = chi_Kerr**2*u**2+2D0*u*nuK_1+nuK_1**2
      Long_Delta1 = 1D0+u*(Delta_1+u*(Delta_2+u*(Delta_3+u*Delta_4)))
      Long_Delta2 = 1D0+nu*Delta_0+nu*Log(Long_Delta1)
      Delta_u = Delta_Bar_u*Long_Delta2
      Delta_t = R**2*Delta_u
      varpi_square = R**2 + a**2
      D_Minus_u = 1D0 + Log(1D0+6D0*nu*u**2+2D0*(26D0-3D0*nu)*nu*u**3)
      Delta_r = Delta_t*D_Minus_u
      omega_tilde = 2D0*a*M*R
      Sigma = R**2+a**2*COS(theta)**2
      Lambda_t = varpi_square**2 - a**2*Delta_t*SIN(theta)**2
      dDeltau_du =(nu*(nuK_1**2+2D0*u*nuK_1+chi_Kerr**2*u**2)*(Delta_1+u*(2D0*Delta_2+u* &
                  (3D0*Delta_3+4D0*Delta_4*u))))/Long_Delta1+2D0*(nuK_1+chi_Kerr**2*u)*Long_Delta2
      dDeltat_dr = 2D0*R*Delta_u-dDeltau_du*M
      dD_Minus_u_du = (12D0*u*nu+6D0*u**2*(26D0-3D0*nu)*nu)/(1D0+6D0*u**2*nu+2D0*nu*u**3*(26D0-3D0*nu))
      dDeltar_dr = dDeltat_dr*D_Minus_u-Delta_t*dD_Minus_u_du*M/R**2
      dLambdat_dr = 4D0*R*varpi_square-a**2*SIN(theta)**2*dDeltat_dr
      Parameter_Calculator(:)=(/Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r,dDeltar_dr,dDeltat_dr,dLambdat_dr/)
      RETURN
      END FUNCTION
  SUBROUTINE Derivative_Of_Parameter(coordinates,Parameters,D_Parameter)
    IMPLICIT NONE
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: R,theta,Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    REAL * 8 :: dDeltar_dr,dDeltat_dr,dLambdat_dr
    REAL * 8 :: gtt,gtf,dgtt_dr,dgff_dr,dgtf_dr,dgtt_dth,dgff_dth,dgtf_dth
    REAL * 8 :: coeff1,coeff2,coeff3
    REAL * 8 :: Lambda_t2,Sigma_2,theta_2
    REAL * 8 :: coordinates(3),Parameters(9),D_Parameter(10)
    R = coordinates(1)
    theta = coordinates(2)
    Sigma = Parameters(1)
    Lambda_t = Parameters(2)
    varpi_square = Parameters(3)
    omega_tilde = Parameters(4)
    Delta_t = Parameters(5)
    Delta_r = Parameters(6)
    dDeltar_dr = Parameters(7)
    dDeltat_dr = Parameters(8)
    dLambdat_dr = Parameters(9) 
    gtt = -Lambda_t/(Delta_t*Sigma)
    gtf = -omega_tilde/(Delta_t*Sigma)
    dgtt_dr = (Lambda_t*(Sigma*dDeltat_dr+2*R*Delta_t)-dLambdat_dr*Delta_t*Sigma)/(Delta_t*Sigma)**2
    dgtt_dth =  a**2*Sin(2*theta)*varpi_square*(Delta_t-varpi_square)/(Delta_t*Sigma**2)
    dgff_dr = (2*R*Lambda_t-dLambdat_dr*Sigma)/(Lambda_t**2*Sin(theta)**2)+ &
              (-8*a**2*M**2*R*Lambda_t*Delta_t*Sigma+omega_tilde**2*(dLambdat_dr*Delta_t*Sigma+ &
              Lambda_t*dDeltat_dr*Sigma+2*R*Lambda_t*Delta_t))/(Lambda_t*Sigma*Delta_t)**2
    dgff_dth = a**2*Sin(2*theta)/Lambda_t**2*(-Lambda_t/Sin(theta)**2+&
               Sigma*Delta_t/Sin(theta)**2-Lambda_t*Sigma/(a**2*Sin(theta)**4)-&
               omega_tilde**2*(Sigma*Delta_t+Lambda_t)/(Sigma**2*Delta_t))
    dgtf_dr = (omega_tilde*(dDeltat_dr*Sigma+2*R*Delta_t)-2*a*M*Sigma*Delta_t)/(Delta_t*Sigma)**2
    dgtf_dth = -a**2*omega_tilde*Sin(2*theta)/(Sigma**2*Delta_t)
    coeff1 = 1/(2D0*(-gtt)**1.5D0)
    coeff2 =  2D0*gtf/gtt
    coeff3 = (gtf/gtt)**2
    Lambda_t2 =1D0/Lambda_t**2
    Sigma_2 = 1D0/Sigma**2
    theta_2 = 2D0*theta
    D_Parameter(1) = coeff1*dgtt_dr
    D_Parameter(2) = coeff1*dgtt_dth
    D_Parameter(3) = omega_tilde*(Lambda_t/R-dLambdat_dr)*Lambda_t2
    D_Parameter(4) = omega_tilde*a**2*Sin(theta_2)*Delta_t*Lambda_t2
    D_Parameter(5) = (dDeltar_dr*Sigma - 2*R*Delta_r)*Sigma_2
    D_Parameter(6) = Delta_r*a**2*Sin(theta_2)*Sigma_2
    D_Parameter(7) =  -2*R*Sigma_2
    D_Parameter(8) = a**2*Sin(theta_2)*Sigma_2
    D_Parameter(9) = dgff_dr - (coeff2*dgtf_dr-coeff3*dgtt_dr)
    D_Parameter(10) = dgff_dth - (coeff2*dgtf_dth-coeff3*dgtt_dth)
    RETURN
    END SUBROUTINE
  FUNCTION   H_S(coordinates,Momentum,Parameters)
    IMPLICIT NONE
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
    COMMON /group3/ S_Kerr,norm_S_Kerr,sigma_Star
    REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: S_Kerr(3),norm_S_Kerr,sigma_Star(3)
    REAL * 8 :: coordinates(3),Momentum(3),Parameters(9)
    REAL * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r,dDeltar_dr,dDeltat_dr,dLambdat_dr
    REAL * 8 :: d_SO,d_SS,R,theta,phi,p_R,p_theta,p_phi,u,n(3),grr
    REAL * 8 :: omega,omega_r,omega_COS,e_2nu,e_2mutilde,B_rtilde,B_tilde,J_tilde,mu_r,mu_COS,nu_r,nu_COS,Q
    REAL * 8 :: e_nu,e_mutilde,Sqrt_Q,S_v,S_n,S_xi,S_SKerr
    REAL * 8 :: Delta_sigma_1(3),Delta_sigma_2(3),Delta_sigma_3(3),S_Star(3),S(3),xi(3),v(3),X(3),norm_S_Star
    REAL * 8 :: H_SO,H_SS,Spin_to_itself,H_S
    d_SO = -69.5D0                                              
    d_SS = 2.75D0 
    R = coordinates(1)
    u = M/R
    theta = coordinates(2)
    phi = coordinates(3)
    p_R = Momentum(1)
    p_theta = Momentum(2)
    p_phi = Momentum(3)
    Sigma = Parameters(1)
    Lambda_t = Parameters(2)
    varpi_square = Parameters(3)
    omega_tilde = Parameters(4)
    Delta_t = Parameters(5)
    Delta_r = Parameters(6)
    dDeltat_dr = Parameters(8)
    dLambdat_dr = Parameters(9) 
    grr = Delta_r/Sigma
    Q  = 1 + grr*p_R**2D0+p_phi**2D0*Sigma/(Lambda_t*SIN(theta)**2D0)+p_theta**2D0/Sigma 
    Sqrt_Q = SQRT(Q)
    Delta_sigma_1(:) = sigma_Star(:) * (7D0/6D0*nu*u + nu/3D0*(Q-1)-5D0/2D0*nu*grr*p_R**2D0)+&
                       S_Kerr(:) * (-2D0/3D0*nu*u+1D0/4D0*nu*(Q-1)-3D0*nu*grr*p_R**2D0)
    Delta_sigma_2(:) = sigma_Star(:)*(1D0/36D0*(353D0*nu-27D0*nu**2D0)*u**2D0+5D0*nu**2D0*grr**2D0*p_R**4D0-&
                       1D0/72D0*(23D0*nu+3D0*nu**2D0)*(Q-1D0)**2D0+1D0/36D0*(-103D0*nu+60D0*nu**2D0)*u*(Q-1D0)+&
                       1D0/12D0*(16D0*nu-21D0*nu**2D0)*grr*p_R**2D0*(Q-1D0)+1D0/12D0*(47D0*nu-54D0*nu**2D0)*u*grr*p_R**2D0)+&
                       S_Kerr(:)*(1D0/9D0*(-56D0*nu-21D0*nu**2D0)*u**2D0+45D0/8D0*nu**2D0*grr**2D0*p_R**4D0-&
                        5D0/16D0*nu*(Q-1D0)**2D0+1D0/36D0*(-109D0*nu+51D0*nu**2D0)*u*(Q-1D0)+&
                        1D0/8D0*(2D0*nu-13D0*nu**2D0)*grr*p_R**2D0*(Q-1D0)-1D0/24D0*(16D0*nu+147D0*nu**2D0)*u*grr*p_R**2D0)
    Delta_sigma_3(:) = d_SO*nu/R**3D0*sigma_Star(:)
    S_Star(:) = sigma_Star(:) + Delta_sigma_1(:) + Delta_sigma_2(:) + Delta_sigma_3(:)     
    X(:) = R*(/SIN(theta)*COS(phi),SIN(theta)*SIN(phi),COS(theta)/)               
    norm_S_Star = Sqrt(S_Star(1)**2D0+S_Star(2)**2D0+S_Star(3)**2D0)      
    Spin_to_itself = -(R**2D0*norm_S_Star**2D0 - 3D0*DOT_PRODUCT(X,S_Star)**2D0)/(2D0*M*R**5D0)                                                                     
    S(:) = nu * S_Star(:) 
    n(:) = (/SIN(theta)*COS(phi),SIN(theta)*SIN(phi),COS(theta)/) 
    omega = omega_tilde / Lambda_t                                                                             
    omega_r = omega_tilde*(-dLambdat_dr+Lambda_t/R)/Lambda_t**2                                          
    omega_cos = -2*a**2*COS(theta)*Delta_t*omega_tilde/Lambda_t**2                                              
    e_2nu = Delta_t*Sigma/Lambda_t 
    e_nu = Sqrt(e_2nu)                                                                            
    e_2mutilde = Sigma
    e_mutilde = Sqrt(e_2mutilde)                                                                                        
    B_tilde = Sqrt(Delta_t)                                                                                    
    B_rtilde = (Sqrt(Delta_r)*dDeltat_dr-2*Delta_t)/(2*Sqrt(Delta_t*Delta_r))                               
    J_tilde = Sqrt(Delta_r)                                                                                    
    mu_r = R/Sigma - 1D0/Sqrt(Delta_r)                                                                           
    mu_cos = a**2*COS(theta)/Sigma                                                                             
    nu_r = R/Sigma + varpi_square*(varpi_square*dDeltat_dr-4D0*R*Delta_t)/(2*Lambda_t*Delta_t)                
    nu_cos = a**2*varpi_square*COS(theta)*(varpi_square-Delta_t)/(Lambda_t*Sigma)                                                                 
    xi(:) = SIN(theta)*(/-SIN(phi),COS(phi),0D0/)                                                            
    v(:) = -SIN(theta)*(/COS(theta)*COS(phi),COS(theta)*SIN(phi),-SIN(theta)/) 
    S_v = DOT_PRODUCT(S,v)
    S_n = DOT_PRODUCT(S,n)
    S_xi = DOT_PRODUCT(S,xi)
    S_SKerr = DOT_PRODUCT(S,S_Kerr/norm_S_Kerr)   
    H_SO = e_2nu/e_mutilde*(e_mutilde*e_nu-B_tilde)*p_phi*&
           S_SKerr/(B_tilde**2*Sqrt_Q*Sin(theta)**2)+&
           e_nu/(e_2mutilde*B_tilde**2*(Sqrt_Q+1)*Sqrt_Q*Sin(theta)**2)*&
           (S_xi*J_tilde*(-mu_r*Sin(theta)*p_theta*(Sqrt_Q+1)-mu_cos*p_R*Sin(theta)**2-&
           Sqrt_Q*(-nu_r*Sin(theta)*p_theta+(mu_cos-nu_cos)*p_R*Sin(theta)**2))*B_tilde**2 +&
           Sqrt(e_2mutilde*e_2nu)*p_phi*(2*Sqrt_Q+1)*B_tilde*(J_tilde*nu_r*S_v-nu_cos*S_n*&
           Sin(theta)**2)-J_tilde*B_rtilde*Sqrt(e_2mutilde*e_2nu)*p_phi*(Sqrt_Q+1)*S_v)
    H_SS = omega*S_SKerr + J_tilde*omega_r/(2*B_tilde*(Sqrt_Q+1)*Sqrt_Q*Sin(theta)**2*&
           (Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*(Sqrt(e_2mutilde*e_2nu)*Sin(theta)*p_theta*p_phi*S_xi*B_tilde+&
           e_2nu*e_2mutilde*p_phi**2*S_v+e_2mutilde*(1+Sqrt_Q)*Sqrt_Q*S_v*Sin(theta)**2*B_tilde**2+&
           J_tilde*p_R*(-Sin(theta)*p_theta*S_n-J_tilde*p_R*S_v)*Sin(theta)**2*B_tilde**2)+&
           omega_cos/(2*B_tilde*(Sqrt_Q+1)*Sqrt_Q*(Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*&
           (-p_phi**2*S_n*e_2mutilde*e_2nu+e_nu*e_mutilde*J_tilde*p_R*p_phi*S_xi*B_tilde+&
           (S_n*(p_theta*Sin(theta))**2+J_tilde*p_R*S_v*Sin(theta)*p_theta-&
           e_2mutilde*(1+Sqrt_Q)*Sqrt_Q*S_n*Sin(theta)**2)*B_tilde**2)
    H_S = (H_SS+H_SO)/mu + Spin_to_itself + d_SS*nu*DOT_PRODUCT(S_Kerr,sigma_Star)/R**4
    return
    END FUNCTION
  SUBROUTINE Evolution_RK78(h,coordinates,Momentum,t,t_f)
    IMPLICIT NONE
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
    COMMON /group4/ Lz
    REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: Lz
    REAL * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    REAL * 8 :: Parameters(9),D_Parameter(10)
    REAL * 8 :: grr,gthth,alpha,beta,gamma
    REAL * 8 :: Q_4,square_Q,Energy,dH_dHeff,Carter,h,delta,coordinates(3),Momentum(3),t,t_f
    REAL * 8 :: X(5),X1(5),Y0(5),Y1(5),Y2(5),Y3(5),Y4(5),Y5(5)
    REAL * 8 :: Y6(5),Y7(5),Y8(5),Y9(5),Y10(5),Y11(5),Y12(5)
    100 format(6f25.10)
    200 format(1f10.4,2f25.15)
    delta = 1E-5*M
    DO while( t .LE. t_f)
      X(1:3) = coordinates(:)
      X(4:5) = Momentum(1:2)
      X1(:) = X(:)
    DO i = 1,5
      X1(i) = X(i)
    END DO
    CALL Differential_Eq(delta,X1,Y0)
    DO i = 1,5
      X1(i) = X(i) + h*2D0/27D0*Y0(i)
    END DO
    CALL Differential_Eq(delta,X1,Y1)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)+3D0*Y1(i))/36D0
    END DO 
    CALL Differential_Eq(delta,X1,Y2)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)+3D0*Y2(i))/24D0
    END DO
    CALL Differential_Eq(delta,X1,Y3)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)*20D0+(-Y2(i)+Y3(i))*75D0)/48D0
    END DO 
    CALL Differential_Eq(delta,X1,Y4)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)+Y3(i)*5D0+Y4(i)*4D0)/20D0   
    END DO   
    CALL Differential_Eq(delta,X1,Y5)
    DO i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*25D0+Y3(i)*125D0-Y4(i)*260D0+Y5(i)*250D0)/108D0
    END DO 
    CALL Differential_Eq(delta,X1,Y6)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)*93D0+Y4(i)*244D0-Y5(i)*200D0+Y6(i)*13D0)/900D0
    END DO 
    CALL Differential_Eq(delta,X1,Y7)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)*180D0-Y3(i)*795D0+Y4(i)*1408D0-Y5(i)*1070D0 + &      
                  Y6(i)*67D0+Y7(i)*270D0)/90D0
    END DO 
    CALL Differential_Eq(delta,X1,Y8)
    DO i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*455D0+Y3(i)*115D0-Y4(i)*3904D0+Y5(i)*3110D0 - &
                  Y6(i)*171D0+Y7(i)*1530D0-Y8(i)*45D0)/540D0
    END DO
    CALL Differential_Eq(delta,X1,Y9)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)*2383D0-Y3(i)*8525D0+Y4(i)*17984D0-Y5(i)*15050D0+ &
                  Y6(i)*2133D0+Y7(i)*2250D0+Y8(i)*1125D0+Y9(i)*1800D0)/4100D0
    END DO
    CALL Differential_Eq(delta,X1,Y10)
    DO i = 1,5
      X1(i) = X(i) + h*(Y0(i)*60D0-Y5(i)*600D0-Y6(i)*60D0+(Y8(i)-Y7(i)+&
                  2D0*Y9(i))*300D0)/4100D0
    END DO
    CALL Differential_Eq(delta,X1,Y11)
    DO i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*1777D0-Y3(i)*8525D0+Y4(i)*17984D0-&
                  Y5(i)*14450D0+Y6(i)*2193D0+Y7(i)*2550D0 + &
                  Y8(i)*825D0+Y9(i)*1200D0+Y11(i)*4100D0)/4100D0
    END DO
    CALL Differential_Eq(delta,X1,Y12)
    DO i = 1,5
      X(i) = X(i) + h*(Y5(i)*272D0+(Y6(i)+Y7(i))*216D0+(Y8(i)+Y9(i))*27D0 + &
                    (Y11(i)+Y12(i))*41D0)/840D0
    END DO
      t = t + h
      coordinates(:) = X(1:3)
      Momentum(1:2) = X(4:5)
      t = t + h
      Parameters =  Parameter_Calculator(coordinates)
      Sigma = Parameters(1)
      Lambda_t = Parameters(2)
      varpi_square = Parameters(3)
      omega_tilde = Parameters(4)
      Delta_t = Parameters(5)
      Delta_r = Parameters(6)
      grr = Delta_r/Sigma
      gthth = 1/Sigma
      alpha = 1/SQRT(Lambda_t/(Delta_t*Sigma))
      beta = omega_tilde/Lambda_t
      gamma = 1/Lambda_t*(-omega_tilde**2/(Delta_t*Sigma)+Sigma/SIN(coordinates(2))**2)+omega_tilde**2/(Lambda_t*Delta_t*Sigma)
      Q_4 = 2*(4-3*nu)*nu*M**2/Sigma*Momentum(1)**4
      square_Q = SQRT(1+grr*Momentum(1)**2+gthth*Momentum(2)**2+gamma*Momentum(3)**2+Q_4)
      Energy = beta*Lz+ alpha*square_Q + H_S(coordinates,Momentum,Parameters)
      dH_dHeff = 1/SQRT(1+2*nu*(Energy-1))
      Carter = Momentum(2)**2 + COS(coordinates(2))**2*((1-Energy**2)*a**2+Lz**2/SIN(coordinates(2))**2)
      WRITE(10,100)coordinates,Momentum
      WRITE(20,200)t/M,Energy,Carter
    END DO
    RETURN
    END SUBROUTINE
  SUBROUTINE Differential_Eq(delta,X1,Y)
    IMPLICIT NONE
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
    COMMON /group4/ Lz
    REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: Lz
    REAL * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    REAL * 8 :: grr,gthth,alpha,beta,gamma
    REAL * 8 :: Parameters(9),D_Parameter(10)
    REAL * 8 :: Q_4,square_Q,Energy,dH_dHeff,Carter
    REAL * 8 :: delta,X1(5),Y(5)
    REAL * 8 :: Temp_X(3),Temp_P(3),f1,f2

    Parameters =  Parameter_Calculator(X1(1:3))
    CALL Derivative_Of_Parameter(X1(1:3),Parameters,D_Parameter)
    Temp_X(:) = X1(1:3)
    Temp_P(1:2) = X1(4:5)
    Temp_P(3) = Lz 
    Sigma = Parameters(1)
    Lambda_t = Parameters(2)
    varpi_square = Parameters(3)
    omega_tilde = Parameters(4)
    Delta_t = Parameters(5)
    Delta_r = Parameters(6)
    grr = Delta_r/Sigma
    gthth = 1/Sigma
    alpha = 1/SQRT(Lambda_t/(Delta_t*Sigma))
    beta = omega_tilde/Lambda_t
    gamma = 1/Lambda_t*(-omega_tilde**2/(Delta_t*Sigma)+Sigma/SIN(Temp_X(2))**2)+omega_tilde**2/(Lambda_t*Delta_t*Sigma)
    Q_4 = 2*(4-3*nu)*nu*M**2/Sigma*Temp_P(1)**4
    square_Q = SQRT(1+grr*X1(4)**2+gthth*X1(5)**2+gamma*Lz**2+Q_4)
    Energy = beta*Lz+ alpha*square_Q + H_S(Temp_X,Temp_P,Parameters)
    dH_dHeff = 1/SQRT(1+2*nu*(Energy-1))

    Temp_P(1) = X1(4) + delta/2D0
    f1 = H_S(Temp_X,Temp_P,Parameters)
    Temp_P(1) = X1(4) - delta/2D0 
    f2 = H_S(Temp_X,Temp_P,Parameters)
    Y(1) = dH_dHeff*(alpha*(2*grr*X1(4)+8*(4-3*nu)*nu*M**2*X1(4)**3/Sigma)/(2*square_Q)+(f1-f2)/delta) !R_Dot
    Temp_P(1) = X1(4)
    Temp_P(2) = X1(5) + delta/2D0
    f1 = H_S(Temp_X,Temp_P,Parameters)
    Temp_P(2) = X1(5) - delta/2D0 
    f2 = H_S(Temp_X,Temp_P,Parameters)
    Y(2) = dH_dHeff*(alpha*gthth*X1(5)/square_Q+(f1-f2)/delta) !Theta_Dot
    Temp_P(2) = X1(5)
    Temp_P(3) = Lz + delta/2D0
    f1 = H_S(Temp_X,Temp_P,Parameters)
    Temp_P(3) = Lz - delta/2D0 
    f2 = H_S(Temp_X,Temp_P,Parameters)   
    Y(3) = dH_dHeff*(beta + alpha*gamma*Lz/square_Q+(f1-f2)/delta) !Phi_Dot
    Temp_P(3) = Lz
    Temp_X(1) = X1(1) + delta/2D0
    Parameters =  Parameter_Calculator(Temp_X)
    f1 = H_S(Temp_X,Temp_P,Parameters)
    Temp_X(1) = X1(1) - delta/2D0
    Parameters =  Parameter_Calculator(Temp_X)
    f2 = H_S(Temp_X,Temp_P,Parameters)
    Y(4) = -dH_dHeff*(D_Parameter(3)*Lz+D_Parameter(1)*square_Q+alpha*(D_Parameter(5)*X1(4)**2+&
            D_Parameter(7)*X1(5)**2+D_Parameter(9)*Lz**2-4*(4-3*nu)*nu*M**2*X1(1)*&
            X1(4)**4/Sigma**2)/(2*square_Q)+(f1-f2)/delta) !P_R_Dot
    Temp_X(1) = X1(1)
    Temp_X(2) = X1(2) + delta/2D0
    Parameters =  Parameter_Calculator(Temp_X)   
    f1 = H_S(Temp_X,Temp_P,Parameters)
    Temp_X(2) = X1(2) - delta/2D0
    Parameters =  Parameter_Calculator(Temp_X)
    f2 = H_S(Temp_X,Temp_P,Parameters)
    Y(5) = -dH_dHeff*(D_Parameter(4)*Lz+D_Parameter(2)*square_Q+alpha*(D_Parameter(6)*X1(4)**2+&
            D_Parameter(8)*X1(5)**2+D_Parameter(10)*Lz**2+2*(4-3*nu)*nu*M**2*a**2*&
            SIN(2*X1(2))*X1(4)**4/Sigma**2)/(2*square_Q)+(f1-f2)/delta) !P_Theta_Dot
    return 
    END Subroutine
END program