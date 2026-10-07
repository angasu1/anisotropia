 module potparam                  
   use types
       save                             
       integer(ik)                     ::  indicador               
       
 end module potparam 

Module Potencial_HeOCS
   use types  
   implicit none
   real(rk)                            :: vx(3),ax(3),alpha(19,13),x(19)
   integer(ik)                         :: n_lam, np, npts
   contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function PHEOCS(r,cth)
   real(rk)                            :: PHEOCS
   real(rk), intent(in)                :: r, cth
   real (rk)                           :: Re!, costh  
   integer(ik)                         :: indicador
!        program potocs
   !    use potparam
   !    implicit double precision (a-h,o-z)

        Re=r*0.52917706_rk !debe estar en bhor, porque el factor este convierte a Angstrom 
        !costh = dcos(th)
       
       !if (indicador.ne.1) then
       open(unit=35,file='../tendencia_aniso/ocs_datos.dat',status="old",action="read")
       call reading_OCS 
       ! call EXTINT
       !read *
       close(35)
       !indicador=1
       !endif 

        PHEOCS=Extpot(Re,cth)
        !write(*,*) PHEOCS
        !return     
   end function PHEOCS

   FUNCTION EXTPOT(R,COSTH)
!
!     Program to project out in legendre polys
!     in the angular co-ordinate.
   real(rk)                            :: EXPOT
   real (rk),intent(in)                :: R, costh  
   !integer(ik)                         :: MXSITE,mxparm,mxdata       
   !integer,PARAMETER                   :: MXSITE=10, mxparm=500,mxdata=500
   !integer(ik)                         :: mxparm, mxdata, mxsite
   !real(rk)                            :: parm0(mxparm), rdata(mxdata),p0(mxparm)
   !real(rk)                            :: vx(3),ax(3),alpha(19,13),x(19) 
   real(rk)                            :: RLAM(50),BLAM(50)!,vx(10),ax(10)!,sumleg
   real(rk)                            :: rnw,cthnw !alpha(19,13),x(19)
   real(rk)                            :: x_k,pot,vpp,ra,ro,rs,rc,rb
   real(rk)                            :: eps_scal, r_scal, extpot, extint,ang,aur,aue
   integer(ik)                         :: n_r,k
   !integer(ik)                         :: n_lam, np, npts
   integer(ik)                         :: i,j,n!,npts
   
!      common npts
!     
      DATA RLAM/50*0._rk/, BLAM/50*0._rk/
      DATA AUE/219474.6354_rk/
      DATA AUR/0.52917706_rk/

      !n_lam = 13_ik !nlam
      !np = 3_ik     !no se aun
      !npts = 18_ik !npoint = 19
      ro = 1.68029_rk
      rs = 1.03731_rk

      ra = (r*r+rs*rs+2.0_rk*r*rs*costh)**0.5_rk
      rb = (r*r+ro*ro-2.0_rk*r*ro*costh)**0.5_rk

      rnw = (ra + rb) /(ro+rs)
      cthnw = (ra - rb) /(ro+rs)
     
      eps_scal = sumleg(vx,np,cthnw)
      r_scal = rnw*sumleg(ax,np,cthnw)
      
      pot=0.0_rk

      do j = 1,n_lam  
         vpp = 0.0_rk
         do i=1,npts 
            x_k=x(i)
            vpp=vpp+alpha(i,j)*rad_kernel(x_k,r_scal) !vpp parte radial 
         enddo
         pot = pot+ vpp*plm(j-1,0,cthnw)
      enddo

      extpot=eps_scal*pot
      !write(*,*) extpot, eps_scal/10000

      return
   end FUNCTION EXTPOT!(R,COSTH)

!-------------------------------------------------------------------
   subroutine reading_OCS
   integer(ik), parameter              :: maxLinea=10000
   real(rk),allocatable                :: matA(:,:), matB(:), matC(:,:)
   integer(ik)                         :: i,j,k, nLinea,ioerror, linepot
   character(4)                        :: temp

   nLinea = 0_ik
       do i = 1, maxLinea
           read(35,*,iostat=ioerror) temp
           if(ioerror /= 0) exit
           if (nLinea > maxLinea) stop "Archivo demasiado grande"
           nLinea = nLinea + 1   
       end do! i = 1, maxLinea

       rewind(35) ! (linea del caset)  

       allocate(matA(nLinea, 2))        

       ! Leer la matriz, de los coeficientes                        
       do i =1, nLinea                      
           read(35,*) (matA(i,j), j = 1, 2)
       end do! i =1, nLinea       
       
!     read number of scaling paras
       np = matA(1,1)
     
!      read scaling paras
       do i = 1,np  
           vx(i) = matA(i+1,1)
           ax(i) = matA(i+1,2)
       end do !i = 2,np

       linepot = np + 2_ik              
       allocate(matB(nLinea-linepot)) 
!     read number of alphas(=no. of rvalues and theta values)
       npts = matA(linepot,1)
       n_lam = matA(linepot,2) 
       
!     NB:x(i)=radial pt at which have abinitio pt
!     NB:y(i)=(1-cos(theta))/2.0
!     where theta is angular pt at which have abinitio pt
!     read x,y and alphas for  RK
       
       do i = 1, nLinea-linepot
             matB(i) = matA(i+linepot,2)
       end do !i = 1, nLinea-linepot

!     NB:x(i)=radial pt at which have abinitio pt
       do i = 1, npts
           x(i) = matA(i+linepot,1) 
           !write(*,*) "xi",x(i)
           !read *
       end do !i = 1, npts
       
!      alphas for RK
       alpha = reshape(matB,(/19,13/))

   end subroutine reading_OCS
!-------------------------------------------------------------------


!===================================================================
 !      ENTRY EXTINT                      
 !     read number of scaling paras      
 !      read(35,*) np                      
 !     read scaling paras                
 !      do i=1,np                         
 !         read(35,*) vx(i), ax(i)         
 !      enddo                             
 !     read number of alphas(=no. of rvalues and theta values)
 !      read(35,*) npts, n_lam             
 !     NB:x(i)=radial pt at which have abinitio pt
 !     NB:y(i)=(1-cos(theta))/2.0        
 !     where theta is angular pt at which have abinitio pt
 !     read x,y and alphas for  RK       
                                         
 !      do j = 1,n_lam                    
 !         do i=1, npts                   
 !            read(35,*) x(i),alpha(i,j)   
 !         enddo                          
 !      enddo                             
                                         
 !      return                            
 !      end  
!===================================================================

 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function rad_kernel(x,r)
      !implicit none
      real(rk),intent(in)              :: x,r
      real(rk)                         :: rad_kernel,x_l,x_s,R0
      !real*8                      rad_kernel,x,r,x_l,x_s,R0
      
      x_l=r
      x_s=x
      if(r.lt.x) then
         x_s=r
         x_l=x
      endif
!     kernel for protonated systems (long range prop to R-4)
!     rad_kernel=2.d0/(15.d0*x_l**5)*(1.d0-5.d0*x_s/(7.d0*x_l))
!     kernel for protonated systems (long range prop to R-6)m=5,n=2
!     rad_kernel=2.d0/(21.d0*x_l**6.0d0)*(1.d0-3.d0*x_s/(4.d0*x_l))
!     m=6, n=2
      rad_kernel=1._rk/(14._rk*x_l**7.0_rk)*(1._rk-7._rk*x_s/(9._rk*x_l))
!     Increase order of smoothness.n=3
!     rad_kernel=1.0d0/(28.0d0*x_l**7.0d0)*(1.d0-7.d0*x_s/(5.d0*x_l)
!     . + 28.d0*x_s*x_s/(55.d0*x_l*x_l))
!     n=3 m=5
!     rad_kernel=3.d0/(56.0d0*x_l**6)*(1.d0-4.d0*x_s/(3.d0*x_l)
!     . + 7.d0*x_s*x_s/(15.d0*x_l*x_l))
      
!     rad_kernel=1.d0/(14.0d0*(R0+x_l)**7)*
!     .    (1.d0-7.d0*(x_s+R0)/(9.d0*(x_l+R0)))
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




   FUNCTION PLM(L,MM,X)

!     COMPUTES NORMALIZED ASSOC. LEGENDRE POLYNOMIALS BY RECURSION.
!C     THE VALUES RETURNED ARE NORMALIZED FOR INTEGRATION OVER X
!C     (I.E. INTEGRATION OVER COS THETA BUT NOT PHI).
!C     NOTE THAT THE NORMALIZATION GIVES 
!C           PLM(L,0,1)=SQRT(L+0.5)
!C           PLM(L,0,X)=SQRT(L+0.5) P(L,X)
!C     FOR M.NE.0, THE VALUE RETURNED DIFFERS FROM THE USUAL
!C           DEFINITION OF THE ASSOCIATED LEGENDRE POLYNOMIAL
!C           (E.G. EDMONDS PAGES 23-24) 
!C           BY A FACTOR OF (-1)**M*SQRT(L+0.5)*SQRT((L-M)!/(L+M)!)
!C     THUS THE SPHERICAL HARMONICS ARE 
!C          CLM = PLM * EXP(I*M*PHI) / SQRT(L+0.5)
!C          YLM = PLM * EXP(I*M*PHI) / SQRT(2*PI)
!C
!C     PROGRAM OF R. NERF MODIFED BY S. GREEN.
!C       MOD FEB. 82 BY S.G. ACCORDING TO R.T PACK'S SUGGESTION
!C       FOR IMPROVED ACCURACY LARGE ABS(X)
!C     MODIFIED AUG 88 BY S.G. TO KEEP RAT FROM BLOWING UP AT LARGE L
!C       BY INCLUDING FACTOR (-Z)**L IN PRECEDING LOOP
!C     ERROR IN STMT. NO. 211 CORRECTED MAY 93 BY SG
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   real(rk)                            :: PLM, X, P1,P2,P3,RAT
   integer (ik)                        :: L, MM, LM1, M, MU  
   integer(ik)                         :: MTEST, I, XI,XL,XLP,XM,XMU
   real(rk)                            :: XNORM, XTEST, Z, XLM
  
      M=IABS(MM)
      IF ( M.GT.L .OR. L.LT.0)  GO TO 9999
      XL=L
      XM=M
      P1=1.D0
      P2=X
      IF(L.EQ.0) GO TO 20
!C *** USE ALTERNATE RECURSION FOR LARGE M OR LARGE ABS(X)
      IF (M.EQ.0)  GO TO 10
      XTEST=0.49999D0*(1.D0+1.D0/XM)
      MTEST=L/3
      IF (M.GT.MTEST .OR. ABS(X).GT.XTEST)  GO TO 210
!C ***
   10 DO 100 I=1,L
      XI=I
      P3=((2.D0*XI+1.D0)*X*P2-XI*P1)/(XI+1.D0)
      P1=P2
  100 P2=P3
      IF (M.EQ.0) GO TO 20
!C AT END OF LOOP P1=P(L,0,X)
      IF (ABS(X).GT.1.D0)  GO TO 9999
      Z=SQRT(1.D0-X*X)
      IF (Z.LE.1.D-10) GO TO 999
  201 P2=(XL+1.D0)*(P2-X*P1)/Z
      DO 200 I=1,M
      XI=I
      P3=-2.D0*X/Z*P2*XI-(XL+XI)*(XL-XI+1.D0)*P1
      P1=P2
  200 P2=P3
      GO TO 20
!C ***
!C     BELOW RECURS DOWN IN M, AS SUGGESTED BY R. T PACK
  210 IF (ABS(X).GT.1.D0)  GO TO 9999
      Z=SQRT(1.D0-X*X)
      IF (Z.LE.1.D-10)  GO TO 999
!C     CALCULATE RATIO OF FACTORIALS FOR PLL
      RAT=1.D0
      XI=0.D0
      DO 211 I=1,L
      XI=XI+1.D0
  211 RAT=-0.5D0*Z*(XL+XI)*RAT
!C     CALCULATE PLL   ----  N.B. ABOVE INCL (-Z)**L COMPUTATION
!C     P1=RAT*(-Z)**L
      P1=RAT
      IF (M.EQ.L)  GO TO 20
      P2=P1
      P1=-X*P2/Z
      IF (M.EQ.L-1)  GO TO 20
      LM1=L-M-1
!C     RECUR DOWNWARD IN M
      DO 213 I=1,LM1
      MU=L-I-1
      XMU=MU
      P3=P2
      P2=P1
      P1=2.D0*(XMU+1.D0)*X*P2/Z+P3
  213 P1=-P1/((XL-XMU)*(XL+XMU+1.D0))
      GO TO 20
!C ***
!C     NORMALIZATION . . .
   20 XNORM=(2.D0*XL+1.D0)/2.D0
      IF (M.LE.0)  GO TO 1000
      XLM=XL+1.D0
      XLP=XL
      DO 1100 I=1,M
      XLM=XLM-1.D0
      XLP=XLP+1.D0
 1100 XNORM=XNORM/(XLM*XLP)
 1000 PLM=P1*SQRT(XNORM)
      RETURN
 9999 WRITE(6,699)  L,MM,X
  699 FORMAT('0 * * * ERROR.  ARGUMENT OUT OF RANGE FOR PLM(',2I6,D16.8, ' ).')
    
!C     IF Z=0, THEN X=1 AND PLM(1.)=0 FOR M.GT.0.
!C       IN THAT CASE, WE HAVE BRANCHED TO 999
  999 PLM=0.D0
      RETURN
      END


      FUNCTION SUMLEG(COEFF,NP,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!C   
!C     SUMLEG EVALUATES A LEGENDRE SERIES AT A GIVEN ANGLE THETA, USING
!C     THE RECURSION RELATIONSHIP FOR LEGENDRE POLYNOMIALS.
!C
!C     ON INPUT, COEFF IS THE LEGENDRE SERIES, STARTING AT P0
!C               NP    IS THE ORDER OF THE LEGENDRE SERIES
!C               X     IS COS(THETA)
!C
      DIMENSION COEFF(NP)
      !integer(ik)                      :: NP, k
      integer                      :: NP, k
      IF(NP.GE.1) GOTO 1
      WRITE(6,601)NP
  601 FORMAT(/' **** ERROR IN SUMLEG: NP =',I5)
      STOP
    1 SUMLEG=COEFF(1)
      IF(NP.EQ.1) RETURN
      SUMLEG=SUMLEG+X*COEFF(2)
      IF(NP.EQ.2) RETURN
      P0=1.D0
      P1=X
      DO 10 K=3,NP
      TEMP=(DBLE(K+K-3)*X*P1 - DBLE(K-2)*P0) / DBLE(K-1)
      P0=P1
      P1=TEMP
   10 SUMLEG=SUMLEG+P1*COEFF(K)
      RETURN
      END
end Module Potencial_HeOCS
