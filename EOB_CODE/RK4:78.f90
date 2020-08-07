  SUBROUTINE Evolution_RK4(h,coordinates,Momentum,t,t_f)
    IMPLICIT NONE
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
    COMMON /group4/ Lz
    REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: Lz
    REAL * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    REAL * 8 :: Parameters(9),D_Parameter(10)
    REAL * 8 :: grr,gthth,alpha,beta,gamma
    REAL * 8 :: Q_4,square_Q,Energy,dH_dHeff,Carter,h,delta,coordinates(3),Momentum(3),t,t_f
    REAL * 8 :: X(5),X1(5),Y0(5),Y1(5),Y2(5),Y3(5)
    100 format(6f25.10)
    200 format(1f10.4,2f25.15)
    delta = 1E-5 * M
    DO WHILE (t .LE. t_f)
      X(1:3) = coordinates(:)
      X(4:5) = Momentum(1:2)
      X1(:) = X(:)
      CALL Differential_Eq(delta,X1,Y0)
      X1(:) = X(:) + h*Y0(:)/2D0
      CALL Differential_Eq(delta,X1,Y1)
      X1(:) = X(:) + h*Y1(:)/2D0
      CALL Differential_Eq(delta,X1,Y2)
      X1(:) = X(:) + h*Y2(:)
      CALL Differential_Eq(delta,X1,Y3)
      X(:) = X(:) + h*(Y0(:)+2D0*Y1(:)+2D0*Y2(:)+Y3(:))/6D0
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
  Subroutine Evolution_RK78(h,coordinates,Momentum,t,t_f)
    Implicit None
    COMMON /group1/ M,mu,nu,a,chi_Kerr,K,nuK_1
    COMMON /group4/ Lz
    REAL * 8 :: M,mu,nu,a,chi_Kerr,K,nuK_1
    REAL * 8 :: Lz
    REAL * 8 :: Sigma,Lambda_t,varpi_square,omega_tilde,Delta_t,Delta_r
    REAL * 8 :: Parameters(9),D_Parameter(10)
    REAL * 8 :: grr,gthth,alpha,beta,gamma
    REAL * 8 :: Q_4,square_Q,Energy,dH_dHeff,Carter,h,delta,coordinates(3),Momentum(3),t,t_f
    real * 8 :: X(5),X1(5),Y0(5),Y1(5),Y2(5),Y3(5),Y4(5),Y5(5)
    real * 8 :: Y6(5),Y7(5),Y8(5),Y9(5),Y10(5),Y11(5),Y12(5)
    100 format(6f25.10)
    200 format(1f10.4,2f25.15)
    delta = 1E-5*M
    do while( t .le. t_f)
      X(1:3) = coordinates(:)
      X(4:5) = Momentum(1:2)
      X1(:) = X(:)
    do i = 1,5
      X1(i) = X(i)
    end do
    call Differential_Eq(delta,X1,Y0)
    do i = 1,5
      X1(i) = X(i) + h*2D0/27D0*Y0(i)
    end do
    call Differential_Eq(delta,X1,Y1)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)+3D0*Y1(i))/36D0
    end do 
    call Differential_Eq(delta,X1,Y2)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)+3D0*Y2(i))/24D0
    end do
    call Differential_Eq(delta,X1,Y3)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*20D0+(-Y2(i)+Y3(i))*75D0)/48D0
    end do 
    call Differential_Eq(delta,X1,Y4)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)+Y3(i)*5D0+Y4(i)*4D0)/20D0   
    end do    
    call Differential_Eq(delta,X1,Y5)
    do i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*25D0+Y3(i)*125D0-Y4(i)*260D0+Y5(i)*250D0)/108D0
    end do 
    call Differential_Eq(delta,X1,Y6)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*93D0+Y4(i)*244D0-Y5(i)*200D0+Y6(i)*13D0)/900D0
    end do 
    call Differential_Eq(delta,X1,Y7)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*180D0-Y3(i)*795D0+Y4(i)*1408D0-Y5(i)*1070D0 + &      
                  Y6(i)*67D0+Y7(i)*270D0)/90D0
    end do 
    call Differential_Eq(delta,X1,Y8)
    do i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*455D0+Y3(i)*115D0-Y4(i)*3904D0+Y5(i)*3110D0 - &
                  Y6(i)*171D0+Y7(i)*1530D0-Y8(i)*45D0)/540D0
    end do
    call Differential_Eq(delta,X1,Y9)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*2383D0-Y3(i)*8525D0+Y4(i)*17984D0-Y5(i)*15050D0+ &
                  Y6(i)*2133D0+Y7(i)*2250D0+Y8(i)*1125D0+Y9(i)*1800D0)/4100D0
    end do
    call Differential_Eq(delta,X1,Y10)
    do i = 1,5
      X1(i) = X(i) + h*(Y0(i)*60D0-Y5(i)*600D0-Y6(i)*60D0+(Y8(i)-Y7(i)+&
                  2D0*Y9(i))*300D0)/4100D0
    end do
    call Differential_Eq(delta,X1,Y11)
    do i = 1,5
      X1(i) = X(i) + h*(-Y0(i)*1777D0-Y3(i)*8525D0+Y4(i)*17984D0-&
                  Y5(i)*14450D0+Y6(i)*2193D0+Y7(i)*2550D0 + &
                  Y8(i)*825D0+Y9(i)*1200D0+Y11(i)*4100D0)/4100D0
    end do
    call Differential_Eq(delta,X1,Y12)
    do i = 1,5
      X(i) = X(i) + h*(Y5(i)*272D0+(Y6(i)+Y7(i))*216D0+(Y8(i)+Y9(i))*27D0 + &
                    (Y11(i)+Y12(i))*41D0)/840D0
    end do
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
    end do
    return
    end Subroutine