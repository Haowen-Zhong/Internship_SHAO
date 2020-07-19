program EOB
implicit none
real * 16, parameter :: Pi=3.1415926535d0
real * 16 :: m_1,m_2,M,mu,S_1,S_2,S_Kerr,S_star,a,chi_Kerr,nu,K,omega_1,omega_2,alpha,beta,g(5)
real * 16 :: r,theta,phi,p_theta,p_phi,u,theta_min,theta_max
real * 16 :: Sigma,varpi_square,Lambda_t,omega_tilde,Delta_t,Delta_r,D_Minus_u,Delta_Bar_u,Delta_u
real * 16 :: Delta_0,Delta_1,Delta_2,Delta_3,Delta_4
real * 16 :: partial_alpha_partial_r,partial_alpha_partial_theta,partial_beta_phi_partial_theta,partial_gamma_thetatheta_partial_r
real * 16 :: partial_gamma_thetatheta_partial_theta,partial_gamma_phiphi_partial_r,partial_gamma_phiphi_partial_theta
real * 16 :: The_guy_with_alpha,H_hat,H_eff_over_mu,Partial_H_hat_partial_H_eff_over_mu,E
real * 16 :: theta_dot_,phi_dot_,p_theta_dot_

!***********************************************************************************************
!以下参数为根据守恒关系求出初始的p_\theta和p_\phi，E而定义
!q := \parial_r\beta^\phi; b := \partial_r\alpha; d := \partial_r\gamma^{\phi\phi}
!gtt_prime := \partial_r g^{tt}; gtphi_prime := \partial_r g^{t\phi}; gphiphi_prime := \partial_r g^{\phiphi}
!Lambda_t_prime := \partial_r \Lambda_t; Delta_t_prime := \partial_r \Delta_t
!Delta_u_prime := \partial_u\Delta_u
!***********************************************************************************************

real * 16 :: q,b,c,d,gtt_prime,gtphi_prime,gphiphi_prime,Lambda_t_prime,Delta_t_prime,Delta_u_prime
!***********************************************************************************************
!以下参数为采用Runge-Kutta方法求解微分方程而定义

real * 16 :: K11,K12,K13,K14,K15,K16
real * 16 :: K21,K22,K23,K24,K25,K26
real * 16 :: K31,K32,K33,K34,K35,K36
real * 16 :: t_f,t,t_0,h,e_max,e_min,h_max,h_min,n_max
real * 16 :: p_theta_1,p_theta_2,theta_1,theta_2,phi_1,phi_2,R_1,R_2,R_3,Diffrence

integer :: i

!***********************************************************************************************
!输出数据以及定义输出格式

open (10,file='data_EOB_RKF45.txt',status='unknown')
100 format(5f15.5)
!***********************************************************************************************
!——————————————————人为给定两黑洞的参数——————————————————————————
!这一部分经过检查没有问题

m_1 = 1d0
m_2 = 0.1d0

M = m_1 + m_2

mu = (m_1 * m_2) / M
nu = mu / M

S_1 = 0.5 * m_1 **2
S_2 = 0.5 * m_2 **2
S_Kerr = S_1 + S_2
a = S_Kerr/M
chi_Kerr = a/M

K = 1.447 - 1.715*nu - 3.246*nu**2
omega_1 = 0
omega_2 = 0

!theta极值
theta_min = Pi/2.
theta_max = Pi/2.
!开始的时候从极值点开始积分，如果用守恒条件的话
!会非常难以处理，这样去计算的话对于Pi/2情况也适用
p_theta = 0

!圆轨道初值
r = 8. * M 
u = M/r
theta = theta_min
phi = 0

t_0 = 0
t_f = 10000
t = t_0
h_min = 0.005*M
h_max = 0.05*M
e_min = 1E-7
e_max = 1E-5
n_max = t_f/h_min
h = h_max
!********************常数不需要调整********************
!经过检查没有问题

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

Delta_t = Delta_u*r**2

varpi_square = r**2 + a**2

D_Minus_u = 1 + Log(1+6*nu*u**2+2*(26-3*nu)*nu*u**3)

Delta_r = Delta_t * D_Minus_u

omega_tilde = 2*a*M*r + omega_1*nu*a*M**3*u + omega_2*nu*M*a**3*u

!********************常数不需要调整********************

Sigma = r**2+a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

!*****************************************根据守恒量计算p_phi,"E"****************************************
!由于Delta_t_prime的结果正确，Delta_u_prime理应正确
Delta_u_prime =(nu*((-1 + K*nu)**(-2) + (2*u)/(-1 + K*nu) + chi_Kerr**2*u**2)* &
          (Delta_1 + u*(2*Delta_2 + u*(3*Delta_3 + 4*Delta_4*u))))/ &
        (1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))) +  &
       2*(1/(-1 + K*nu) + chi_Kerr**2*u)* &
        (1 + Delta_0*nu + nu*Log(1 + u*(Delta_1 + u*(Delta_2 + u*(Delta_3 + Delta_4*u)))))

!和Mathematica的结果相同
Delta_t_prime = 2*r*Delta_u - Delta_u_prime*M

!和Mathematica的结果相同
Lambda_t_prime = 4*r*varpi_square - a**2*Sin(theta)**2*Delta_t_prime

!和Mathematica的结果相同
gtt_prime = (Lambda_t*(Sigma*Delta_t_prime + 2*r*Delta_t) - Lambda_t_prime * Delta_t * Sigma)/(Delta_t*Sigma)**2

!和Mathematica的结果相同
gtphi_prime =(omega_tilde*(Delta_t_prime * Sigma + 2 * r * Delta_t)-2*a*M*Sigma*Delta_t)/(Delta_t*Sigma)**2

!和Mathematica的结果相同
gphiphi_prime = (2*r*Lambda_t-Lambda_t_prime*Sigma)/(Lambda_t**2*Sin(theta)**2)+ (-8*a**2*M**2*r*Lambda_t*Delta_t*&
    Sigma+omega_tilde**2*(Lambda_t_prime*Delta_t*Sigma+Lambda_t*Delta_t_prime*Sigma&
    +2*r*Lambda_t*Delta_t))/(Lambda_t*Sigma*Delta_t)**2



!q = partial_r\beta^\phi 和Mathematica的结果相同

q = (2*a*M*Lambda_t-Lambda_t_prime*2*a*M*r)/Lambda_t**2

!b = paritla_r\alpha 和Mathematica的结果相同

b = 1/(2.*(-g(1))**1.5d0)*(Lambda_t*(Sigma*Delta_t_prime+2*r*Delta_t)-Lambda_t_prime*Delta_t*Sigma)/(Delta_t*Sigma)**2

!c = \gamma^{\phi\phi}
c = g(4) - g(5)**2/g(1)

!d = \partial_r\gamma^{\phi\phi}
d = gphiphi_prime - (2*g(5)*g(1)*gtphi_prime-g(5)**2*gtt_prime)/g(1)**2

alpha = 1/Sqrt(-g(1))

!对得上！是正确的！
p_phi = Sqrt(2*b**2/(q**2-2*b**2*c-b*d*alpha+Sqrt(q**4-2*b*d*alpha*q**2)))


!可以对得上！
beta = g(5)/g(1)
E = beta*p_phi + alpha*Sqrt(1+c*p_phi**2)

write(*,*)"The angular momentum L equals to:",p_phi
write(*,*)"The Effective Energe E equals to:",E
!\hat{H}_{real}
H_hat = M/mu*(Sqrt(1+2*nu*(E-1))-1)
!\frac{\partial\hat{H}_{real}}{H_eff/\mu}
Partial_H_hat_partial_H_eff_over_mu = M**2*nu/(mu*(mu*H_hat+M))
!——————————————————计算完成————————————————————————————————————

do while (i .le. n_max .and. t .le. t_f )
!————————————————微分方程————————————————————

K11 = h * P_Theta_Dot(p_theta,p_phi,H_hat,r,theta,a,varpi_square,&
                           Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)

K21 = h * Theta_dot(p_theta,p_phi,theta,r,a,varpi_square,H_hat,Lambda_t,&
                       Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)

K31 = h * Phi_dot(p_theta,p_phi,theta,r,a,varpi_square,H_hat,Lambda_t,&
                   Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)

K12 = h * P_Theta_Dot(p_theta + 1/4.*K11,p_phi,H_hat,r,theta + 1/4.*K21,a,varpi_square,&
                           Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)
K22 = h * Theta_dot(p_theta + 1/4.*K11,p_phi,theta + 1/4.*K21,r,a,varpi_square,H_hat,Lambda_t,&
                       Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
K32 = h * Phi_dot(p_theta + 1/4.*K11,p_phi,theta + 1/4.*K21,r,a,varpi_square,H_hat,Lambda_t,&
                   Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)

K13 = h * P_Theta_Dot(p_theta +3/32.*K11+9/32.*K12,p_phi,H_hat,r,theta+3/32.*K21+9/32.*K22,a,varpi_square,&
                           Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)
K23 = h * Theta_dot(p_theta+3/32.*K11+9/32.*K12,p_phi,theta+3/32.*K21+9/32.*K22,r,a,varpi_square,H_hat,Lambda_t,&
                       Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
K33 = h * Phi_dot(p_theta+3/32.*K11+9/32.*K12,p_phi,theta+3/32.*K21+9/32.*K22,r,a,varpi_square,H_hat,Lambda_t,&
                   Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)


K14 = h * P_Theta_Dot(p_theta+1932/2197.*K11-7200/2197.*K12+7296/2197.*K13,p_phi,H_hat,r,theta+1932/2197.*&
  K21-7200/2197.*K22+7296/2197.*K23,a,varpi_square,Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)
K24 = h * Theta_dot(p_theta+1932/2197.*K11-7200/2197.*K12+7296/2197.*K13,p_phi,theta+1932/2197.*&
  K21-7200/2197.*K22+7296/2197.*K23,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
K34 = h * Phi_dot(p_theta+1932/2197.*K11-7200/2197.*K12+7296/2197.*K13,p_phi,theta+1932/2197.*&
  K21-7200/2197.*K22+7296/2197.*K23,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)

K15 = h * P_Theta_Dot(p_theta+439/216.*K11-8.*K12+3680/513.*K13-845/4104.*K14,p_phi,H_hat,r,theta+439/216.*&
  K21-8.*K22+3680/513.*K23-845/4104.*K24,a,varpi_square,Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)
K25 = h * Theta_dot(p_theta+439/216.*K11-8.*K12+3680/513.*K13-845/4104.*K14,p_phi,theta+439/216.*&
  K21-8.*K22+3680/513.*K23-845/4104.*K24,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
K35 = h * Phi_dot(p_theta+439/216.*K11-8.*K12+3680/513.*K13-845/4104.*K14,p_phi,theta+439/216.*&
  K21-8.*K22+3680/513.*K23-845/4104.*K24,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)


K16 = h * P_Theta_Dot(p_theta-8/27.*K11+2.*K12-3544/2565.*K13+1859/4104.*K14-11/40.*K15,p_phi,H_hat,r,theta-8/27.*&
  K21+2.*K22-3544/2565.*K23+1859/4104.*K24-11/40.*K25,&
  a,varpi_square,Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)
K26 = h * Theta_dot(p_theta-8/27.*K11+2.*K12-3544/2565.*K13+1859/4104.*K14-11/40.*K15,p_phi,theta-8/27.*&
  K21+2.*K22-3544/2565.*K23+1859/4104.*K24-11/40.*K25,&
  r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
K36 = h * Phi_dot(p_theta-8/27.*K11+2.*K12-3544/2565.*K13+1859/4104.*K14-11/40.*K15,p_phi,theta-8/27.*&
  K21+2.*K22-3544/2565.*K23+1859/4104.*K24-11/40.*K25,&
  r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)

p_theta_1 = p_theta + 25/256.*K11 + 1408/2565.*K13 + 2197/4101.*K14 - 1/5.*K15
p_theta_2 = p_theta + 16/135.*K11 + 6656/12858.*K13 + 28561/56430.*K14 -9/50.*K15 + 2/55.*K16
R_1 = Abs(p_theta_1-p_theta_2)/h

theta_1 = theta + 25/256.*K21 + 1408/2565.*K23 + 2197/4101.*K24 - 1/5.*K25
theta_2 = theta + 16/135.*K21 + 6656/12858.*K23 + 28561/56430.*K24 -9/50.*K25 + 2/55.*K26
R_2 = Abs(theta_1-theta_2)/h

phi_1 = phi + 25/256.*K31 + 1408/2565.*K33 + 2197/4101.*K34 - 1/5.*K35
phi_2 = phi + 16/135.*K31 + 6656/12858.*K33 + 28561/56430.*K34 -9/50.*K35 + 2/55.*K36
R_3 = Abs(phi_1-phi_2)/h

Diffrence = Max(R_1,R_2,R_3)
if (Diffrence .ge. e_max .and. h .ge. h_min) then
  h = h/2
Else 
  p_theta = p_theta_2
  theta = theta_2
  phi = phi_2
  t = t + h
  i = i + 1
end if
if (Diffrence < e_min) then
  h = 2.*h
end if  
!————————————————输出结果——————————————————————

write(10,100) r*M,theta,phi
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

function P_Theta_Dot(p_theta,p_phi,H_hat,r,theta,a,varpi_square,Partial_H_hat_partial_H_eff_over_mu,Delta_t,Delta_r,omega_tilde)
real * 16 :: p_theta,p_phi,H_hat,r,theta,a,M,mu,nu,Delta_t,Delta_r,omega_tilde
real * 16 :: Partial_H_hat_partial_H_eff_over_mu,partial_beta_phi_partial_theta,partial_alpha_partial_theta
real * 16 :: partial_g_tt_partial_theta,partial_g_tphi_partial_theta
real * 16 :: partial_g_thetatheta_partial_theta,partial_g_phiphi_partial_theta
real * 16 :: The_guy_with_alpha,alpha,partial_gamma_thetatheta_partial_theta,partial_gamma_phiphi_partial_theta
real * 16 :: g(5),varpi_square
real * 16 :: p_theta_dot

Sigma = r**2+a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g(:) = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

alpha = 1/Sqrt(-g(1))

!———————为计算alpha,beta,gamma导数先计算度规的导数—————————————————
partial_g_tt_partial_theta = a**2*Sin(2*theta)*(a**2*Delta_t+r**2*Delta_t-varpi_square**2)/(Delta_t*Sigma**2)

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
The_guy_with_alpha = Sqrt(1+g(3)*p_theta**2+(g(4)-g(5)**2/g(1))*p_phi**2)

p_theta_dot = - Partial_H_hat_partial_H_eff_over_mu*(partial_beta_phi_partial_theta*p_phi+ &
        partial_alpha_partial_theta*The_guy_with_alpha+alpha*(partial_gamma_thetatheta_partial_theta*&
            p_theta**2+partial_gamma_phiphi_partial_theta*p_phi**2)/The_guy_with_alpha)


return 
end function

!——————————————————————Theta_Dot————————————————————————————————

function Theta_dot(p_theta,p_phi,theta,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
real * 16 :: p_theta,p_phi,theta,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r
real * 16 :: Partial_H_hat_partial_H_eff_over_mu,alpha,The_guy_with_alpha
real * 16 :: Sigma,g(5)
real * 16 :: Theta_dot


Sigma = r**2 + a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g(:) = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

alpha = 1/Sqrt(-g(1))

!The_guy_with_alpha 用来指代根号下的一长串东西
The_guy_with_alpha = Sqrt(1+g(3)*p_theta**2+(g(4)-g(5)**2/g(1))*p_phi**2)

Theta_dot = Partial_H_hat_partial_H_eff_over_mu*(alpha*g(3)*p_theta/The_guy_with_alpha)

return 
end function
!——————————————————————Phi_Dot————————————————————————————————

function Phi_dot(p_theta,p_phi,theta,r,a,varpi_square,H_hat,Lambda_t,Delta_t,Delta_r,Partial_H_hat_partial_H_eff_over_mu)
real * 16 :: p_theta,p_phi,theta,r,a,H_hat,Lambda_t,Delta_t,Delta_r
real * 16 :: Partial_H_hat_partial_H_eff_over_mu ,alpha,The_guy_with_alpha
real * 16 :: Sigma,varpi_square,g(5)
real * 16 :: Phi_dot


Sigma = r**2 + a**2*Cos(theta)**2

Lambda_t = varpi_square**2 - a**2*Delta_t*Sin(theta)**2

g(:) = Metric(Lambda_t,Delta_t,Delta_r,Sigma,omega_tilde,theta)

alpha = 1/Sqrt(-g(1))

!The_guy_with_alpha 用来指代根号下的一长串东西

The_guy_with_alpha = Sqrt(1+g(3)*p_theta**2+(g(4)-g(5)**2/g(1))*p_phi**2)

Phi_dot = Partial_H_hat_partial_H_eff_over_mu*(alpha*(g(4)-g(5)**2/g(1))*p_phi/The_guy_with_alpha)

return 
end function

end program
