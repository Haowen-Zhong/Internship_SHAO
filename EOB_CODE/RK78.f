  
      Subroutine Evolution_fit(h,r,phi,pr,Lz,tob,EE)
      implicit real*8 (a-h,o-z)
      parameter (N=4)
      real*8 Lz
      DIMENSION X(N),X1(N),Y0(N),Y1(N),Y2(N),Y3(N),Y4(N),Y5(N),
     &          Y6(N),Y7(N),Y8(N),Y9(N),Y10(N),Y11(N),Y12(N)
      common /group1/ Ub,Up,a,Pi
      common /group4/ Utot,Ueff,rat

c     H=0.01d0*Utot
      X(1)=r
      X(2)=phi
      X(3)=pr
      X(4)=Lz      
      T=tob
      
      DO 16 I=1,N
   16 X1(I)=X(I)
      T1=T
      CALL Evolution_YHC(X1,Y0,N)
      DO 1 I=1,N
    1 X1(I)=X(I)+H*2D0/27D0*Y0(I)
      T1=T+H*2D0/27D0
      CALL Evolution_YHC(X1,Y1,N)
      DO 2 I=1,N
    2 X1(I)=X(I)+H*(Y0(I)+3D0*Y1(I))/36D0
      T1=T+H*1D0/9D0
      CALL Evolution_YHC(X1,Y2,N)
      DO 3 I=1,N
    3 X1(I)=X(I)+H*(Y0(I)+3D0*Y2(I))/24D0
      T1=T+H*1D0/6D0
      CALL Evolution_YHC(X1,Y3,N)
      DO 4 I=1,N
    4 X1(I)=X(I)+H*(Y0(I)*20D0+(-Y2(I)+Y3(I))*75D0)/48D0
      T1=T+H*5D0/12D0
      CALL Evolution_YHC(X1,Y4,N)
      DO 5 I=1,N
    5 X1(I)=X(I)+H*(Y0(I)+Y3(I)*5D0+Y4(I)*4D0)/20D0
      T1=T+H*1D0/2D0
      CALL Evolution_YHC(X1,Y5,N)
      DO 6 I=1,N
    6 X1(I)=X(I)+H*(-Y0(I)*25D0+Y3(I)*125D0-Y4(I)*260D0+Y5(I)*250D0)
     &      /108D0
      T1=T+H*5D0/6D0
      CALL Evolution_YHC(X1,Y6,N)
      DO 7 I=1,N
    7 X1(I)=X(I)+H*(Y0(I)*93D0+Y4(I)*244D0-Y5(I)*200D0+Y6(I)*13D0)/900D0
      T1=T+H*1D0/6D0
      CALL Evolution_YHC(X1,Y7,N)
      DO 8 I=1,N
    8 X1(I)=X(I)+H*(Y0(I)*180D0-Y3(I)*795D0+Y4(I)*1408D0-Y5(I)*1070D0
     &      +Y6(I)*67D0+Y7(I)*270D0)/90D0
      T1=T+H*2D0/3D0
      CALL Evolution_YHC(X1,Y8,N)
      DO 9 I=1,N
    9 X1(I)=X(I)+H*(-Y0(I)*455D0+Y3(I)*115D0-Y4(I)*3904D0+Y5(I)*3110D0
     &      -Y6(I)*171D0+Y7(I)*1530D0-Y8(I)*45D0)/540D0
      T1=T+H*1D0/3D0
      CALL Evolution_YHC(X1,Y9,N)
      DO 10 I=1,N
   10 X1(I)=X(I)+H*(Y0(I)*2383D0-Y3(I)*8525D0+Y4(I)*17984D0-Y5(I)*15050
     &   D0+Y6(I)*2133D0+Y7(I)*2250D0+Y8(I)*1125D0+Y9(I)*1800D0)/4100D0
      T1=T+H
      CALL Evolution_YHC(X1,Y10,N)
      DO 11 I=1,N
   11 X1(I)=X(I)+H*(Y0(I)*60D0-Y5(I)*600D0-Y6(I)*60D0+(Y8(I)-Y7(I)
     &      +2D0*Y9(I))*300D0)/4100D0
      T1=T
      CALL Evolution_YHC(X1,Y11,N)
      DO 12 I=1,N
   12 X1(I)=X(I)+H*(-Y0(I)*1777D0-Y3(I)*8525D0+Y4(I)*17984D0
     &      -Y5(I)*14450D0+Y6(I)*2193D0+Y7(I)*2550D0
     &      +Y8(I)*825D0+Y9(I)*1200D0+Y11(I)*4100D0)/4100D0
      T1=T+H
      CALL Evolution_YHC(X1,Y12,N)
      EE=0D0
      DO 14 I=1,N
       X(I)=X(I)+H*(Y5(I)*272D0+(Y6(I)+Y7(I))*216D0+(Y8(I)+Y9(I))*27D0
     &     +(Y11(I)+Y12(I))*41D0)/840D0
   14 EE=EE+DABS((Y0(I)+Y10(I)-Y11(I)-Y12(I))*H*41D0/840D0)
      T=T+H
      r=X(1)
      phi=X(2)
      pr=X(3)
      Lz=X(4)
      tob=T

      RETURN
      END
C********************************************************************
      SUBROUTINE Evolution_YHC(X,Y,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Lz, a8(0:10),ah(0:10)
      DIMENSION X(N),Y(N)
      common /group1/ Ub,Up,a,Pi
      common /group4/ Utot,Ueff,rat,s2
      common /e/ dEdt,dLdt,w
      
      
      Y(1)=(delt*delr*pr/((tlb*E-wfd*Lz))+Hs_pr)/Eob !prs*At*Ar/E*dett
      Y(2)=(wns+Hs_pf)/Eob
      Y(3)=-(funp(r,Lz,pr)+Hs_r+3*s2**2/(2*Utot*r**4))/Eob+Fr 
      Y(4)=dLdt
      !write(*,*) dLdt,Fr
      !pause
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc