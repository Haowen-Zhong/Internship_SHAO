program EOB
implicit none
real * 16, parameter :: Pi=3.1415926535879792328_16
real * 16 :: m_1,m_2,M,mu,S,S_1,S_2,S_Kerr,S_star_,a,chi_Kerr,nu,K,omega_1,omega_2,alpha,beta,g(5)
real * 16 :: R,theta,phi,p_r,p_rstar,p_theta,p_phi,p_phi_old,u,theta_min,theta_max
real * 16 :: Sigma,varpi_square,Lambda_t,omega_tilde,Delta_t,Delta_r,D_Minus_u,Delta_Bar_u,Delta_u
real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
real * 16 :: partial_alpha_partial_r,partial_alpha_partial_theta,partial_beta_phi_partial_theta,partial_gamma_thetatheta_partial_r
real * 16 :: partial_gamma_thetatheta_partial_theta,partial_gamma_phiphi_partial_r,partial_gamma_phiphi_partial_theta
real * 16 :: The_guy_with_alpha,H_hat,H_eff_over_mu,Partial_H_hat_partial_H_eff_over_mu,Energy
real * 16 :: theta_dot_,phi_dot_,p_theta_dot_
real * 16 :: omega,e_2nu,omega_r,nu_r,mu_r,omega_cos,nu_cos,mu_cos,B_tilde,B_rtilde,e_2mutilde,J_tilde,Q
real * 16 :: Q4,S_Kerr_hat,SS,Spin_to_self,Carter
real * 16 :: lower_sigma,sigma_star,Delta_sigma1,Delta_sigma2,Delta_sigma3,d_SO,d_SS
!***********************************************************************************************
!以下参数为根据守恒关系求出初始的p_\theta和p_\phi，E而定义
!w := \parial_r\beta^\phi; b := \partial_r\alpha; d := \partial_r\gamma^{\phi\phi}
!gtt_prime := \partial_r g^{tt}; gtphi_prime := \partial_r g^{t\phi}; gphiphi_prime := \partial_r g^{\phiphi}
!Lambda_t_prime := \partial_r \Lambda_t; Delta_t_prime := \partial_r \Delta_t
!Delta_u_prime := \partial_u\Delta_u
!***********************************************************************************************

real * 16 :: w,b,c,d,e,f,h,j,gtt_prime,gtphi_prime,gphiphi_prime,Lambda_t_prime,Delta_t_prime,Delta_u_prime

!***********************************************************************************************
!以下参数为采用Runge-Kutta方法求解微分方程以及Jacobi Iteration求解角动量而定义

real * 16 :: K11,K12,K13,K14,K21,K22,K23,K24,K31,K32,K33,K34,epsilon=1E-15_16,error=1
integer :: i = 0
!***********************************************************************************************
!输出数据以及定义输出格式

open (10,file='data_aEOB_RK4_2.txt',status='unknown')
100 format(4 f15.5)
!***********************************************************************************************
!——————————————————人为给定两黑洞的参数——————————————————————————
m_1 = 1.d0
m_2 = 0.000001d0

M = m_1 + m_2


mu = (m_1 * m_2) / M
nu = mu / M


S_1 = 0.5 * m_1 **2
S_2 = 0.5 * m_2 **2
S_Kerr = S_1 + S_2

S_Kerr_hat = Sign(1._16,S_Kerr)!取1的绝对值，S_kerr的符号
lower_sigma = S_1 + S_2
sigma_star = m_2/m_1*S_1 +m_1/m_2*S_2 

a = S_Kerr/M
chi_Kerr = a/M

d_SS = 2.75
d_SO = -69.5


K = 1.447 - 1.715*nu - 3.246*nu**2
omega_1 = 0
omega_2 = 0

!theta极值
theta_min = Pi/4.
theta_max = 3*Pi/4.

write(*,*)theta_min

!开始的时候从极值点开始积分，如果用守恒条件的话
!会非常难以处理，这样去计算的话对于Pi/2情况也适用
p_theta = 0
p_r = 0
!圆轨道初值
R = 8. * M 

u = M/R
theta = theta_min
phi = 0



h = 0.05 * M

!********************常数不需要调整********************

Delta_0 = K*(nu*K-2)

Delta_1 = -2*(nu*K-1)*(K+Delta_0)

Delta_2 = 1/2.*Delta_1*(-4*nu*K+Delta_1+4)-a**2/M**2*(nu*K-1)**2*Delta_0

Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
            (-nu*K+Delta_2+1)-3*a**2/M**2*(nu*K-1)**2*Delta_1)


Delta_4 = 1/12.*(6*a**2/M**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2+3*Delta_1**4-8*(nu*K-1)*&
          Delta_1**3 - 12*Delta_2*Delta_1**2+12*(2*(nu*K-1)*Delta_2+Delta_3)*Delta_1 +&
          12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2+6*(Delta_2**2-4*Delta_3*(nu*K-1)))

Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2


Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))

!******************************************************

Delta_t = R**2*Delta_u

varpi_square = R**2 + a**2

D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)

Delta_r = Delta_t * D_Minus_u

omega_tilde = 2*a*M*R + omega_1*nu*a*M**3/R + omega_2*nu*M*a**3/R

!********************常数不需要调整********************

Sigma = R**2+a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

!*****************************************根据守恒量计算p_phi,"E"****************************************
!这里的p_phi是在不考虑其他项的存在的前提下求出的，我将作为零级近似利用雅可比迭代法去更加逼近真值
Delta_u_prime =(nu*((-1 + K*nu)**(-2) + (2*u)/(-1 + K*nu) + chi_Kerr**2*u**2)* &
          (Delta_1 + u*(2*Delta_2 + u*(3*Delta_3 + 4*Delta_4*u))))/ &
        (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))) +  &
       2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
        (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))

Delta_t_prime = 2*R*Delta_u - Delta_u_prime*M

Lambda_t_prime = 4*R*varpi_square - a**2*Sin(theta)**2*Delta_t_prime

gtt_prime = (Lambda_t*(Sigma*Delta_t_prime + 2*R*Delta_t) - Lambda_t_prime * Delta_t * Sigma)/(Delta_t*Sigma)**2

gtphi_prime =(omega_tilde*(Delta_t_prime * Sigma + 2 * R * Delta_t)-2*a*M*Sigma*Delta_t)/(Delta_t*Sigma)**2

gphiphi_prime = (2*R*Lambda_t-Lambda_t_prime*Sigma)/(Lambda_t**2*Sin(theta)**2)+ (-8*a**2*M**2*R*Lambda_t*Delta_t*&
    Sigma+omega_tilde**2*(Lambda_t_prime*Delta_t*Sigma+Lambda_t*Delta_t_prime*Sigma&
    +2*R*Lambda_t*Delta_t))/(Lambda_t*Sigma*Delta_t)**2



!w = partial_r\beta^\phi 

w = (2*a*M*Lambda_t-Lambda_t_prime*2*a*M*R)/Lambda_t**2

!b = paritla_r\alpha 

b = 1/(2.*(-g(1))**1.5d0)*(Lambda_t*(Sigma*Delta_t_prime+2*R*Delta_t)-Lambda_t_prime*Delta_t*Sigma)/(Delta_t*Sigma)**2

!c = \gamma^{\phi\phi}
c = g(4) - g(5)**2/g(1)

!d = \partial_r\gamma^{\phi\phi}
d = gphiphi_prime - (2*g(5)*g(1)*gtphi_prime-g(5)**2*gtt_prime)/g(1)**2

alpha = 1/Sqrt(-g(1))

p_phi = Sqrt(2*b**2/(w**2-2*b**2*c-b*d*alpha+Sqrt(w**4-2*b*d*alpha*w**2)))
p_phi_old = p_phi

beta = g(5)/g(1)

!********************计算自旋带来的修正*******************************************

Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma

S_star_ = S_star(theta,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)


S = nu*S_star_

!原文式(6),根据量纲分析，原文错误，分母应当是R

SS = d_SS*nu*sigma_star*Sigma/R**4

!原文式(5)最后一项

Spin_to_self = -1/(2*M*R**3)*(1-3*Cos(theta)**2)*S_star_**2

!********************计算自旋带来的修正*******************************************


!************************Jacobi iteration****************************
!利用不考虑自旋项时候计算得到的角动量作为迭代的初值开始雅可比迭代得到修正的角动量值

do while (error .ge. epsilon)
! e =Partial_R(H_SO/mu+H_SS/mu)
  e = (numerical_differentiation(1,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)+&
    numerical_differentiation(2,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma))/mu
!f = 4d_SS\nu/R^5\sigma^*
  f = 4*d_SS*nu*sigma_star*lower_sigma/R**5
!h=3/(2MR^4)(1-3cos^2\theta)S_Star^2
  h =3*(1-3*Cos(theta)**2)*S_star_**2/(2*M*R**4)

  j = 2*S_star_*(1-3*Cos(theta)**2)/(M*R**3)*&
  numerical_differentiation(3,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)

  !p_phi = -1/w*(b*sqrt(1+c*p_phi**2)+alpha*d*p_phi**2/(2*sqrt(1+c*p_phi**2))-f+h-j)
  p_phi = Sqrt(-(w*p_phi+b*Sqrt(1+c*p_phi**2)+e-f+h-j)*2*Sqrt(1+c*p_phi**2)/(alpha*d))
  error = Abs(p_phi_old-p_phi)
  p_phi_old = p_phi
  i = i+1
write(*,*) "error =",error,"p_phi=",p_phi,"times=",i
end do
i = 1
!************************Jacobi iteration****************************

!Q_4用🐢坐标来提高程序的稳定性。好像是在计算波形的时候会用到，现在可以不做变换

!p_rstar = p_r*Sqrt(Delta_t*Delta_r)/varpi_square
!这也是无量纲化的，理当正确
!如果不做变换：Q4=2*(4-3*nu)*nu*M**2/Sigma*p_r**4
Q4=2*(4-3*nu)*nu*M**2/Sigma*p_r**4
!Q4 = 2*(4-3*nu)*nu*M**2/Sigma*p_rstar**4*varpi_square**4/(Delta_r*Delta_t)**2

Energy = beta*p_phi + alpha*Sqrt(1+c*p_phi**2)+H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi,K)/mu+&
H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi,K)/mu + SS +Spin_to_self

write(*,*)"The angular momentum L equals to:",p_phi
write(*,*)"The Effective Energy E equals to:",Energy
!\hat{H}_{real}
H_hat = M/mu*(Sqrt(1+2*nu*(Energy-1))-1)
!\frac{\partial\hat{H}_{real}}{H_eff/\mu}
Partial_H_hat_partial_H_eff_over_mu = M**2*nu/(mu*(mu*H_hat+M))
!——————————————————计算完成————————————————————————————————————

do while (i <= 50000)

h = 0.05 * M

S_star_ = S_star(theta,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
!————————————————微分方程————————————————————

p_theta_dot_ = P_Theta_Dot(p_r,p_theta,p_phi,H_hat,R,theta,a,M,mu,varpi_square,&
              Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde,S_star_,sigma_star,lower_sigma)

theta_dot_ = Theta_dot(p_r,p_theta,p_phi,theta,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
              Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)

phi_dot_ = Phi_dot(p_r,p_theta,p_phi,theta,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
          Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)
!————————————————微分方程——————————————————————

!————————————————更新参数——————————————————————

K11 = p_theta_dot_
K21 = theta_dot_
K31 = phi_dot_

K12 = P_Theta_Dot(p_r,p_theta+h/2.*K11,p_phi,H_hat,R,theta+h/2.*K21,a,M,mu,varpi_square,&
      Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde,S_star_,sigma_star,lower_sigma)
K22 = Theta_dot(p_r,p_theta+h/2.*K11,p_phi,theta+h/2.*K21,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
      Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)
K32 = Phi_dot(p_r,p_theta+h/2.*K11,p_phi,theta+h/2.*K21,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
      Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)

K13 = P_Theta_Dot(p_r,p_theta+h/2.*K12,p_phi,H_hat,R,theta+h/2.* K22,a,M,mu,varpi_square,&
      Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde,S_star_,sigma_star,lower_sigma)
K23 = Theta_dot(p_r,p_theta+h/2.*K12,p_phi,theta+h/2.*K22,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
      Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)
K33 = Phi_dot(p_r,p_theta+h/2.*K12,p_phi,theta+h/2.*K22,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
      Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)


K14 = P_Theta_Dot(p_r,p_theta+h*K13,p_phi,H_hat,R,theta+h* K23,a,M,mu,varpi_square,&
      Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde,S_star_,sigma_star,lower_sigma)
K24 = Theta_dot(p_r,p_theta+h*K13,p_phi,theta+h*K23,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
      Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)
K34 = Phi_dot(p_r,p_theta+h*K13,p_phi,theta+h*K23,R,a,M,mu,varpi_square,H_hat,Lambda_t,&
      Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)

p_theta = p_theta + h/6.*(K11+2*K12+2*K13+K14)
theta = theta + h/6.*(K21+2*K22+2*K23+K24)
phi = phi + h/6.*(K31+2*K32+2*K33+K34)


Sigma = R**2+a**2*Cos(theta)**2
Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2
g = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)
alpha = 1/Sqrt(-g(1))
beta = g(5)/g(1)
c = g(4) - g(5)**2/g(1)
SS = d_SS*nu*sigma_star*Sigma/R**4
Spin_to_self = -1/(2*M*R**3)*(1-3*Cos(theta)**2)*S_star_**2

Energy = beta*p_phi + alpha*Sqrt(1+c*p_phi**2)+H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi,K)/mu+&
H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi,K)/mu + SS +Spin_to_self

i = i+1
!————————————————输出结果——————————————————————

if (MOD(i,500) .eq. 0) then
  write(*,*)"Energy=",Energy
end if 
!write(10,100) R,theta,phi,theta_dot_


! Carter = (p_theta*M)**2+Cos(theta)**2*(1-Energy**2)*a**2+(p_phi*M)**2/Sin(theta)**2)
! write(*,*)Carter

end do

write(*,*) "程序运行成功，可以查看分析数据了哦～"
close(10)

contains

!————————————METRIC_TENSOR——————————————
!tt,rr,thetatheta,phiphi,tphi

function Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)
    real* 16 :: Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta
    real* 16 :: g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi,Metric(5)
    g_tt = -Lambda_t/(Delta_t*Sigma)
    g_rr = Delta_r/Sigma
    g_thetatheta = 1/Sigma
    g_phiphi = 1/Lambda_t*(-omega_tilde**2/(Delta_t*Sigma)+Sigma/Sin(theta)**2)
    g_tphi = -omega_tilde/(Delta_t*Sigma)
    Metric(:)= (/g_tt,g_rr,g_thetatheta,g_phiphi,g_tphi/)
    return
end function

!——————————————————————P_Theta_Dot————————————————————————————————

function P_Theta_Dot(p_r,p_theta,p_phi,H_hat,R,theta,a,M,mu,varpi_square,Partial_H_hat_partial_H_eff_over_mu,&
                     Delta_t,Delta_r,omega_tilde,S_star_,sigma_star,lower_sigma)
  real * 16 :: p_r,p_theta,p_phi,H_hat,R,theta,a,M,mu,nu,Delta_t,Delta_r,omega_tilde,u,Q
  real * 16 :: Partial_H_hat_partial_H_eff_over_mu,partial_beta_phi_partial_theta,partial_alpha_partial_theta
  real * 16 :: partial_g_tt_partial_theta,partial_g_tphi_partial_theta
  real * 16 :: partial_g_thetatheta_partial_theta,partial_g_phiphi_partial_theta
  real * 16 :: The_guy_with_alpha,alpha,partial_gamma_thetatheta_partial_theta,partial_gamma_phiphi_partial_theta
  real * 16 :: g(5),varpi_square,S_star_,sigma_star,lower_sigma,P_Theta_Dot

  u = M/R

  Sigma = R**2+a**2*Cos(theta)**2

  Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

  g(:) = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

  alpha = 1/Sqrt(-g(1))

  Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma
  !———————为计算alpha,beta,gamma导数先计算度规的导数—————————————————
  partial_g_tt_partial_theta = a**2*Sin(2*theta)*(a**2*Delta_t+R**2*Delta_t-varpi_square**2)/(Delta_t*Sigma**2)

  partial_g_tphi_partial_theta= -1/Delta_t*(a**2*Sin(2*theta)*omega_tilde)/Sigma**2

  partial_g_thetatheta_partial_theta = a**2*Sin(2*theta)/Sigma**2

  partial_g_phiphi_partial_theta = a**2*Sin(2*theta)/Lambda_t**2*((2*Sigma*Delta_t-Lambda_t)/Sin(theta)**2- &
  omega_tilde**2*(Delta_t*Sigma+Lambda_t)/(Delta_t*Sigma**2))

  !————————————alpha,beta,gamma的导数——————————————————————————————
  partial_alpha_partial_theta = 1./2*(-g(1))**(1.5)*partial_g_tt_partial_theta

  partial_beta_phi_partial_theta =(partial_g_tphi_partial_theta*g(1)-partial_g_tt_partial_theta*g(5))/g(1)**2

  partial_gamma_thetatheta_partial_theta = partial_g_thetatheta_partial_theta

  partial_gamma_phiphi_partial_theta = partial_g_phiphi_partial_theta-(2*g(5)*g(1)*partial_g_tphi_partial_theta &
  -g(5)**2*partial_g_tt_partial_theta)/g(1)**2

  !The_guy_with_alpha 用来指代根号下的一长串东西
  !由于圆轨道的便利不存在Q4的影响
  The_guy_with_alpha = Sqrt(1+g(3)*p_theta**2+(g(4)-g(5)**2/g(1))*p_phi**2)

  p_theta_dot = - Partial_H_hat_partial_H_eff_over_mu*(partial_beta_phi_partial_theta*p_phi+ &
          partial_alpha_partial_theta*The_guy_with_alpha+alpha*(partial_gamma_thetatheta_partial_theta*&
              p_theta**2+partial_gamma_phiphi_partial_theta*p_phi**2)/The_guy_with_alpha+&
          numerical_differentiation(13,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)/mu+&
          numerical_differentiation(14,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)/mu+&
          3*Sin(2*theta)/(2*M*R**3)*S_star_**2-(1-3*Cos(theta)**2)*S_star_/(M*R**3)*&
          numerical_differentiation(15,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma))
          

  return 
  end function

!——————————————————————Theta_Dot————————————————————————————————

function Theta_dot(p_r,p_theta,p_phi,theta,R,a,M,mu,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,&
                   Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)
real * 16 :: p_r,p_theta,p_phi,theta,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,u,Q,M,mu
real * 16 :: Partial_H_hat_partial_H_eff_over_mu,alpha,The_guy_with_alpha,S_star_,sigma_star,lower_sigma
real * 16 :: Sigma,g(5)
real * 16 :: Theta_dot

u = M/R

Sigma = R**2 + a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g(:) = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

alpha = 1/Sqrt(-g(1))

Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma
!The_guy_with_alpha 用来指代根号下的一长串东西
The_guy_with_alpha = Sqrt(1+g(3)*p_theta**2+(g(4)-g(5)**2/g(1))*p_phi**2)

Theta_dot = Partial_H_hat_partial_H_eff_over_mu*(alpha*g(3)*p_theta/The_guy_with_alpha+&
   numerical_differentiation(7,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)/mu+&
   numerical_differentiation(8,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)/mu-&
   (1-3*Cos(theta)**2)*S_star_/(M*R**3)* &
   numerical_differentiation(9,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma))

return 
end function
!——————————————————————Phi_Dot————————————————————————————————

function Phi_dot(p_r,p_theta,p_phi,theta,R,a,M,mu,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,&
                 Partial_H_hat_partial_H_eff_over_mu,S_star_,sigma_star,lower_sigma)
real * 16 :: p_r,p_theta,p_phi,theta,r,a,H_hat,Lambda_t,Delta_t,Delta_r,M,u,Q,mu
real * 16 :: Partial_H_hat_partial_H_eff_over_mu ,alpha,beta,The_guy_with_alpha,S_star_,sigma_star,lower_sigma
real * 16 :: Sigma,varpi_square,g(5)
real * 16 :: Phi_dot

u = M/R

Sigma = r**2 + a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g(:) = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

alpha = 1/Sqrt(-g(1))

beta = g(5)/g(1)

Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma
!The_guy_with_alpha 用来指代根号下的一长串东西

The_guy_with_alpha = Sqrt(1+g(3)*p_theta**2+(g(4)-g(5)**2/g(1))*p_phi**2)

Phi_dot = Partial_H_hat_partial_H_eff_over_mu*(beta + alpha*(g(4)-g(5)**2/g(1))*p_phi/The_guy_with_alpha+&
   numerical_differentiation(10,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)/mu+&
   numerical_differentiation(11,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)/mu-&
   (1-3*Cos(theta)**2)*S_star_/(M*R**3)* &
   numerical_differentiation(12,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma))

return 
end function 

function H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi,K)
real * 16, parameter :: Pi=3.1415926535d0
real * 16 :: omega_tilde,Lambda_t,Sigma,Delta_t,Delta_r,Delta_t_prime,varpi_square,theta,R,a,u,Q,M,chi_Kerr,K
real * 16 :: e_2nu,e_2mutilde,B_tilde,B_rtilde,S,S_Kerr_hat,J_tilde,nu_r,nu_cos,p_r,p_theta,p_phi,omega,nu
real * 16 :: D_Minus_u,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4,Delta_Bar_u,Delta_u,Delta_u_prime,Lambda_t_prime
real * 16 :: mu_r,mu_cos
real * 16 :: H_SO

u = M/R
chi_Kerr = a/M
varpi_square = R**2 + a**2
D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)
omega_tilde = 2*a*M*R
Delta_0 = K*(nu*K-2)
Delta_1 = -2*(nu*K-1)*(K+Delta_0)
Delta_2 = 1/2.*Delta_1*(-4*nu*K+Delta_1+4)-a**2/M**2*(nu*K-1)**2*Delta_0
Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
            (-nu*K+Delta_2+1)-3*a**2/M**2*(nu*K-1)**2*Delta_1)
Delta_4 = 1/12.*(6*a**2/M**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2+3*Delta_1**4-8*(nu*K-1)*&
          Delta_1**3 - 12*Delta_2*Delta_1**2+12*(2*(nu*K-1)*Delta_2+Delta_3)*Delta_1 +&
          12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2+6*(Delta_2**2-4*Delta_3*(nu*K-1)))
Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2
Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))
Delta_t = R**2*Delta_u
Delta_r = Delta_t * D_Minus_u

Delta_u_prime =(nu*((-1 + K*nu)**(-2) + (2*u)/(-1 + K*nu) + chi_Kerr**2*u**2)* &
          (Delta_1 + u*(2*Delta_2 + u*(3*Delta_3 + 4*Delta_4*u))))/ &
        (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))) +  &
       2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
        (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))
Delta_t_prime = 2*R*Delta_u - Delta_u_prime*M


Sigma = R**2 + a**2*Cos(theta)**2
Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2
Lambda_t_prime = 4*R*varpi_square - a**2*Sin(theta)**2*Delta_t_prime
omega = omega_tilde /Lambda_t
e_2nu = Delta_t*Sigma/Lambda_t
e_2mutilde = Sigma
B_tilde = Sqrt(Delta_t)
B_rtilde = (Sqrt(Delta_r)*Delta_t_prime-2*Delta_t)/(2*Sqrt(Delta_t*Delta_r))
J_tilde = Sqrt(Delta_r)
mu_r = R/Sigma - 1/Sqrt(Delta_r)
mu_cos = a**2*Cos(theta)/Sigma
nu_r = R/Sigma + varpi_square*(varpi_square*Delta_t_prime-4*R*Delta_t)/(2*Lambda_t*Delta_t)
nu_cos = a**2*varpi_square*Cos(varpi_square-Delta_t)/(Lambda_t*Sigma)
Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma


H_SO = e_2nu/Sqrt(e_2mutilde)*(Sqrt(e_2mutilde*e_2nu)-B_tilde)*p_phi*S*S_Kerr_hat/(B_tilde**2*Sqrt(Q)*&
Sin(theta)**2)+Sqrt(e_2nu)/e_2mutilde/(B_tilde**2*(Sqrt(Q+1)*Sqrt(Q)*Sin(theta)**2))*&
(Sqrt(e_2mutilde*e_2nu)*p_phi*(2*Sqrt(Q)+1)*B_tilde*(J_tilde*nu_r*(S*Sin(theta)**2)-nu_cos*(S*Cos(theta)*Sin(theta)**2)-&
J_tilde*B_rtilde*Sqrt(e_2mutilde*e_2nu)*p_phi*(Sqrt(Q)+1)*S*Sin(theta)**2))


return
end function

function H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi,K)
real * 16, parameter :: Pi=3.1415926535d0
real * 16 :: omega_tilde,Lambda_t,Sigma,Delta_t,Delta_r,Delta_t_prime,varpi_square,theta,R,a,u,Q,M,chi_Kerr,K
real * 16 :: e_2nu,e_2mutilde,B_tilde,B_rtilde,S,S_Kerr_hat,J_tilde,nu_r,nu_cos,p_r,p_theta,p_phi,omega,omega_r
real * 16 :: D_Minus_u,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4,Delta_Bar_u,Delta_u,Delta_u_prime,Lambda_t_prime
real * 16 :: mu_r,mu_cos,nu
real * 16 :: H_SS

u = M/R
chi_Kerr = a/M
varpi_square = R**2 + a**2
D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)
omega_tilde = 2*a*M*R
Delta_0 = K*(nu*K-2)
Delta_1 = -2*(nu*K-1)*(K+Delta_0)
Delta_2 = 1/2.*Delta_1*(-4*nu*K+Delta_1+4)-a**2/M**2*(nu*K-1)**2*Delta_0
Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
            (-nu*K+Delta_2+1)-3*a**2/M**2*(nu*K-1)**2*Delta_1)
Delta_4 = 1/12.*(6*a**2/M**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2+3*Delta_1**4-8*(nu*K-1)*&
          Delta_1**3 - 12*Delta_2*Delta_1**2+12*(2*(nu*K-1)*Delta_2+Delta_3)*Delta_1 +&
          12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2+6*(Delta_2**2-4*Delta_3*(nu*K-1)))
Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2
Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))
Delta_t = R**2*Delta_u
Delta_r = Delta_t * D_Minus_u

Delta_u_prime =(nu*((-1 + K*nu)**(-2) + (2*u)/(-1 + K*nu) + chi_Kerr**2*u**2)* &
          (Delta_1 + u*(2*Delta_2 + u*(3*Delta_3 + 4*Delta_4*u))))/ &
        (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))) +  &
       2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
        (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))
Delta_t_prime = 2*R*Delta_u - Delta_u_prime*M


Sigma = R**2 + a**2*Cos(theta)**2
Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2
Lambda_t_prime = 4*R*varpi_square - a**2*Sin(theta)**2*Delta_t_prime
omega = omega_tilde /Lambda_t
omega_r = (-Lambda_t_prime*omega_tilde+Lambda_t*2*a*M)/Lambda_t**2
e_2nu = Delta_t*Sigma/Lambda_t
e_2mutilde = Sigma
B_tilde = Sqrt(Delta_t)
B_rtilde = (Sqrt(Delta_r)*Delta_t_prime-2*Delta_t)/(2*Sqrt(Delta_t*Delta_r))
J_tilde = Sqrt(Delta_r)
mu_r = R/Sigma - 1/Sqrt(Delta_r)
mu_cos = a**2*Cos(theta)/Sigma
nu_r = R/Sigma + varpi_square*(varpi_square*Delta_t_prime-4*R*Delta_t)/(2*Lambda_t*Delta_t)
nu_cos = a**2*varpi_square*Cos(varpi_square-Delta_t)/(Lambda_t*Sigma)
Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma

H_SS = omega*(S*S_Kerr_hat) + J_tilde*omega_r/(2*B_tilde*(Sqrt(Q+1))*Sqrt(Q)*Sin(theta)**2*(Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*&
(e_2nu*e_2mutilde*p_phi**2*S*Sin(theta)**2+e_2mutilde*(1+Sqrt(Q))*Sqrt(Q)*S*Sin(theta)**4*B_tilde**2+&
J_tilde*p_r*(-Sin(theta)*p_theta*Cos(theta)-J_tilde*p_r*S*Sin(theta)**2)*Sin(theta)**2*B_tilde**2)+&
omega_cos/(2*B_tilde*(Sqrt(Q+1))*Sqrt(Q)*(Sqrt(e_2mutilde)**3*Sqrt(e_2nu)))*&
(p_phi**2*Cos(theta)/(e_2mutilde*e_2nu)+(Cos(theta)*(p_theta*Sin(theta))**2+J_tilde*p_r*S*Sin(theta)**3*p_theta-&
  e_2mutilde*(1+Sqrt(Q))*Sqrt(Q)*Cos(theta)*Sin(theta)**2)*B_tilde**2)
return
end function

function S_star(theta,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
real * 16 ,  parameter :: Pi=3.1415926535d0
real * 16 :: omega_tilde,Lambda_t,Sigma,Delta_t,Delta_r,Delta_t_prime,varpi_square,theta,R,a,u,Q,M,chi_Kerr,K
real * 16 :: e_2nu,e_2mutilde,B_tilde,B_rtilde,S,S_Kerr_hat,J_tilde,nu_r,nu_cos,p_r,p_theta,p_phi,omega,nu
real * 16 :: D_Minus_u,Delta_0,Delta_1,Delta_2,Delta_3,Delta_4,Delta_Bar_u,Delta_u,Delta_u_prime,Lambda_t_prime
real * 16 :: mu_r,mu_cos,d_SO,sigma_star,lower_sigma,Delta_sigma1,Delta_sigma2,Delta_sigma3
real * 16 :: S_star

u = M/R
D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)
chi_Kerr = a/M
varpi_square = R**2 + a**2
Sigma = R**2 + a**2*Cos(theta)**2
Delta_0 = K*(nu*K-2)
Delta_1 = -2*(nu*K-1)*(K+Delta_0)
Delta_2 = 1/2.*Delta_1*(-4*nu*K+Delta_1+4)-a**2/M**2*(nu*K-1)**2*Delta_0
Delta_3 = 1/3.*(-Delta_1**3+3*(nu*K-1)*Delta_1**2+3*Delta_1*Delta_2-6*(nu*K-1)*&
            (-nu*K+Delta_2+1)-3*a**2/M**2*(nu*K-1)**2*Delta_1)
Delta_4 = 1/12.*(6*a**2/M**2*(Delta_1**2-2*Delta_2)*(nu*K-1)**2+3*Delta_1**4-8*(nu*K-1)*&
          Delta_1**3 - 12*Delta_2*Delta_1**2+12*(2*(nu*K-1)*Delta_2+Delta_3)*Delta_1 +&
          12*(94/3.-41/32.*Pi**2)*(nu*K-1)**2+6*(Delta_2**2-4*Delta_3*(nu*K-1)))
Delta_Bar_u = chi_Kerr**2*u**2+2*u/(nu*K-1)+1/(nu*K-1)**2
Delta_u = Delta_Bar_u*(1+nu*Delta_0+nu*Log(1+u*Delta_1+u**2*Delta_2+u**3*Delta_3+u**4*Delta_4))
Delta_t = R**2*Delta_u
Delta_r = Delta_t * D_Minus_u
Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

Q  = 1 + Delta_r*p_r**2/Sigma + p_phi**2*Sigma/(Lambda_t*Sin(theta)**2) + p_theta**2/Sigma

Delta_sigma1 = sigma_star*(7./6*nu*u+nu/3.*(Q-1)-5/2.*nu*Delta_r/Sigma*p_r**2)+&
              lower_sigma*(-2./3*nu*u+1./4*nu*(Q-1)-3*nu*Delta_r/Sigma*p_r**2)

Delta_sigma2 = sigma_star*(1./36*(353*nu-27*nu**2)*u**2+5*nu**2*(Delta_r/Sigma)**2*p_r**4-&
              1./72*(23*nu+3*nu**2)*(Q-1)**2+1./36*(-103*nu+60*nu**2)*u*(Q-1)+&
              1./12*(16*nu-21*nu**2)*Delta_r/Sigma*p_r**2*(Q-1)+1./12*(47*nu-54*nu**2)*u*Delta_r/Sigma*p_r**2)+&
              lower_sigma*(1./9*(-56*nu-21*nu**2)*u**2+45./8*nu**2*(Delta_r/Sigma)**2*p_r**4-&
              5./16*nu*(Q-1)**2+1./36*(-109*nu+51*nu**2)*u*(Q-1)+&
              1./8*(2*nu-13*nu**2)*Delta_r/Sigma*p_r**2*(Q-1)+&
              1./24*(-16*nu-147*nu**2)*u*Delta_r/Sigma*p_r**2)

Delta_sigma3 = d_SO*nu*sigma_star/(R/M)**3

S_star = sigma_star + Delta_sigma1 + Delta_sigma2 + Delta_sigma3




return
end function

function numerical_differentiation(choice,theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma)

!用差分代替偏微分
!choice == 1  ===> (H_SO)'_R         choice ==2  ===> (H_SS)'_R         choice ==3  ===> (S)'_R
!choice == 4  ===> (H_SO)'_PR        choice ==5  ===> (H_SS)'_PR        choice ==6  ===> (S)'_PR
!choice == 7  ===> (H_SO)'_ptheta    choice ==8  ===> (H_SS)'_ptheta    choice ==9  ===> (S)'_ptheta
!choice == 10 ===> (H_SO)'_pphi      choice ==11 ===> (H_SS)'_pphi      choice ==12 ===> (S)'_pphi
!choice == 13 ===> (H_SO)'_theta     choice ==14 ===> (H_SS)'_theta     choice ==15 ===> (S)'_theta

real * 16 :: theta,M,R,a,u,Q,p_r,p_theta,p_phi,sigma_star,lower_sigma
real * 16 :: f1,f2,f3,f4,h
integer :: choice
real * 16 :: numerical_differentiation
h = 0.00000001_16

if (choice .eq. 1) then
  h  = h * M
  f1 = H_SO(theta,R+2*h,a,nu,M,p_r,p_theta,p_phi,K)
  f2 = H_SO(theta,R+h,a,nu,M,p_r,p_theta,p_phi,K)
  f3 = H_SO(theta,R-h,a,nu,M,p_r,p_theta,p_phi,K)
  f4 = H_SO(theta,R-2*h,a,nu,M,p_r,p_theta,p_phi,K)
else if (choice .eq. 2) then
  h  = h * M 
  f1 = H_SS(theta,R+2*h,a,nu,M,p_r,p_theta,p_phi,K)
  f2 = H_SS(theta,R+h,a,nu,M,p_r,p_theta,p_phi,K)
  f3 = H_SS(theta,R-h,a,nu,M,p_r,p_theta,p_phi,K)
  f4 = H_SS(theta,R+2*h,a,nu,M,p_r,p_theta,p_phi,K)
else if (choice .eq. 3) then
  h  = h * M
  f1 = S_star(theta,R+2*h,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f2 = S_star(theta,R+h,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f3 = S_star(theta,R-h,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f4 = S_star(theta,R-2*h,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
else if (choice .eq. 4 ) then
  f1 = H_SO(theta,R,a,nu,M,p_r+2*h,p_theta,p_phi,K)
  f2 = H_SO(theta,R,a,nu,M,p_r+h,p_theta,p_phi,K)
  f3 = H_SO(theta,R,a,nu,M,p_r-h,p_theta,p_phi,K)
  f4 = H_SO(theta,R,a,nu,M,p_r-2*h,p_theta,p_phi,K)
else if (choice .eq. 5 ) then
  f1 = H_SS(theta,R,a,nu,M,p_r+2*h,p_theta,p_phi,K)
  f2 = H_SS(theta,R,a,nu,M,p_r+h,p_theta,p_phi,K)
  f3 = H_SS(theta,R,a,nu,M,p_r-h,p_theta,p_phi,K)
  f4 = H_SS(theta,R,a,nu,M,p_r-2*h,p_theta,p_phi,K)
else if (choice .eq. 6) then
  f1 = S_star(theta,R,a,nu,M,p_r+2*h,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f2 = S_star(theta,R,a,nu,M,p_r+h,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f3 = S_star(theta,R,a,nu,M,p_r-h,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f4 = S_star(theta,R,a,nu,M,p_r-2*h,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
else if (choice .eq. 7) then
  f1 = H_SO(theta,R,a,nu,M,p_r,p_theta+2*h,p_phi,K)
  f2 = H_SO(theta,R,a,nu,M,p_r,p_theta+h,p_phi,K)
  f3 = H_SO(theta,R,a,nu,M,p_r,p_theta-h,p_phi,K)
  f4 = H_SO(theta,R,a,nu,M,p_r,p_theta-2*h,p_phi,K)
else if (choice .eq. 8) then
  f1 = H_SS(theta,R,a,nu,M,p_r,p_theta+2*h,p_phi,K)
  f2 = H_SS(theta,R,a,nu,M,p_r,p_theta+h,p_phi,K)
  f3 = H_SS(theta,R,a,nu,M,p_r,p_theta-h,p_phi,K)
  f4 = H_SS(theta,R,a,nu,M,p_r,p_theta-2*h,p_phi,K)
else if (choice .eq. 9) then
  f1 = S_star(theta,R,a,nu,M,p_r,p_theta+2*h,p_phi,K,d_SO,sigma_star,lower_sigma)
  f2 = S_star(theta,R,a,nu,M,p_r,p_theta+h,p_phi,K,d_SO,sigma_star,lower_sigma)
  f3 = S_star(theta,R,a,nu,M,p_r,p_theta-h,p_phi,K,d_SO,sigma_star,lower_sigma)
  f4 = S_star(theta,R,a,nu,M,p_r,p_theta-2*h,p_phi,K,d_SO,sigma_star,lower_sigma)
else if (choice .eq. 10) then
  f1 = H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi+2*h,K)
  f2 = H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi+h,K)
  f3 = H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi-h,K)
  f4 = H_SO(theta,R,a,nu,M,p_r,p_theta,p_phi-2*h,K)
else if (choice .eq. 11) then
  f1 = H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi+2*h,K)
  f2 = H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi+h,K)
  f3 = H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi-h,K)
  f4 = H_SS(theta,R,a,nu,M,p_r,p_theta,p_phi-2*h,K)
else if (choice .eq. 12) then
  f1 = S_star(theta,R,a,nu,M,p_r,p_theta+2*h,p_phi,K,d_SO,sigma_star,lower_sigma)
  f2 = S_star(theta,R,a,nu,M,p_r,p_theta+h,p_phi,K,d_SO,sigma_star,lower_sigma)
  f3 = S_star(theta,R,a,nu,M,p_r,p_theta-h,p_phi,K,d_SO,sigma_star,lower_sigma)
  f4 = S_star(theta,R,a,nu,M,p_r,p_theta-2*h,p_phi,K,d_SO,sigma_star,lower_sigma)
else if (choice .eq. 13) then
  h  = h * M
  f1 = H_SO(theta+2*h,R,a,nu,M,p_r,p_theta,p_phi+2*h,K)
  f2 = H_SO(theta,R+h,a,nu,M,p_r,p_theta,p_phi+h,K)
  f3 = H_SO(theta,R-h,a,nu,M,p_r,p_theta,p_phi-h,K)
  f4 = H_SO(theta-2*h,R,a,nu,M,p_r,p_theta,p_phi-2*h,K)
else if (choice .eq. 14) then
  h = h * M
  f1 = H_SS(theta+2*h,R,a,nu,M,p_r,p_theta,p_phi,K)
  f2 = H_SS(theta+h,R,a,nu,M,p_r,p_theta,p_phi,K)
  f3 = H_SS(theta-h,R,a,nu,M,p_r,p_theta,p_phi,K)
  f4 = H_SS(theta-2*h,R,a,nu,M,p_r,p_theta,p_phi,K)
else if (choice .eq. 15) then
  h = h * M
  f1 = S_star(theta+2*h,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f2 = S_star(theta+h,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f3 = S_star(theta-h,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
  f4 = S_star(theta-2*h,R,a,nu,M,p_r,p_theta,p_phi,K,d_SO,sigma_star,lower_sigma)
end if

numerical_differentiation = (-f1+8*f2-8*f3+f4)/(12*h)

return
end function


end program
