program anisotropia_calculation                 
use types                               
use leg_integrator                      
use pothx                               
use math                                
implicit none                           
integer(ik)                            :: ll
integer(ik), parameter                 :: maxLinea=10000, nColumna = 22
integer(ik)                            :: nLinea, ioerror
real(rk)                               :: x1,x2, vmedia, vpromi, vmedmin, temporal
integer(ik)                            :: i,j,k,nlam,npoin, npointh
real(rk), allocatable                  :: vlam(:), pl(:)
real(rk)                               :: rr,th,cth, vserie,vv
real(rk)                               :: rmin,rmax,thmin,thmax, sumvlam, no, s_alpha
real(rk), dimension(:,:), allocatable  :: vlamda, plam
character(4)                           :: temp
real(rk), allocatable                  :: vi(:), vprom(:), vicuad(:), vpromcuad(:)
                                        
                                        
!====== abrir archivo para calculo de pot_interpol=================================                         
open(unit=300,file='polinomios.dat',status='old',action='read',iostat=ioerror)
open(unit=200,file='potlamHCl.dat',status='old',action='read')
     if (ioerror /= 0) stop "No puede abrir el archivo"
open (unit=100,file='datos_HCl_prueba.dat')    
open (unit=500,file='cme_HCl_prueba.dat')    
open (unit=700,file='sep_HCl_prueba.dat')    
                                       
! determinar el número de datos en filas
nLinea = 0                              
do i = 1, maxLinea                      
   read(300,*,iostat=ioerror) temp      
   if (ioerror /= 0) exit               
   if (nLinea > maxLinea) stop "Archivo demasiado grande"
   nLinea = nLinea + 1                  
end do ! i = 1, maxLinea                
write(*,*) "linea", nLinea
rewind(300) ! (linea del caset)         
                                        
! Colocar memoria para la matriz        
allocate(plam(nLinea, nColumna))        
allocate(vlamda(nLinea, nColumna))      
                                       
! Leer la matriz                        
do i = 1, nLinea                        
   read(300,1000) (plam(i,j),j = 1, nColumna)
   read(200,1000) (vlamda(i,j),j = 1, nColumna)
end do ! i = 1, nLinea   
!==================================================================================
x1 = -1.0_rk !lim de integración                             
x2 = 1.0_rk  !lim de integración                             
ind= 4 ! 1 HF, 2 HCl, 3 HBr, 4 HCN, 5 OCS                                 
nlam = 20 !no. de lamdas para el potencial (el máximo son 20)                       
!allocate(vprom(nLinea))
allocate(vpromcuad(nLinea))
call pars_assignment(ind)               

rmin = 4._rk                          
rmax = 12._rk                          
thmin = 0._rk                           
thmax = pi                              
npoin= 100_ik                          

allocate(vprom(npoin))
!Medida de anisotropía
vmedia = 0.0_rk
s_alpha = 0._rk
allocate(vi(npoin))
allocate(vicuad(npoin))
!s = sqrt(sum(Vi - Vmed)^2)/abs(Vmed)*npoint


!=======Se obtiene vmedia en el camino de min energia y el cme ============
do i = 0, npoin                         
   th = thmin+(thmax-thmin)*i/npoin 
   cth= dcos(th)                      
   do j = 0, npoin                      
       rr = rmin+(rmax-rmin)*j/npoin        
       if(ind.eq.5)then
          vv = PHEOCS(rr,cth)
       else 
          vv = V(rr,cth,ind)
       end if!(ind.eq.5)then
       call interpolacion(rr,cth,vlamda,plam,nlam,sumvlam)
       vi(j) = sumvlam  
       !vicuad(j) = sumvlam**2
       !vi(j) = vv  
       !vicuad(j) = vv**2
   end do !j = 0, npoin                 
   write(500,*)  th,th*180._rk/pi, minval(vi)
   vmedia = vmedia + minval(vi) 
end do !i= 1, npoin                     
vmedia = (vmedia/npoin)
write(*,*) "media", vmedia
write(100,*) "media", vmedia

!====== Se obtiene valor de anisotropia ====================   
do i = 0,  npoin                         
   th = thmin+(thmax-thmin)*i/npoin 
   cth= dcos(th)                      
   do j = 0, npoin                      
       rr = rmin+(rmax-rmin)*j/npoin        
       if(ind.eq.5)then
          vv = PHEOCS(rr,cth)
       else 
          vv = V(rr,cth,ind)
       end if!(ind.eq.5)then
       call interpolacion(rr,cth,vlamda,plam,nlam,sumvlam)
       vi(j) = sumvlam 
       !vi(j) = vv 
   end do !j = 0, npoin                 
   s_alpha = s_alpha + (minval(vi) - vmedia)**2!/abs(vmedia) 
end do !i= 1, npoin                     
s_alpha = (1._rk/npoin)*sqrt(s_alpha)/abs(vmedia) 
write(*,*) "anisotropia", s_alpha
write(100,*) "anisotropia", s_alpha


!=============Escritura de contornos========================       
do i = 0, npoin                         
   rr = rmin+(rmax-rmin)*i/npoin        
   do j = 0, npoin                      
       th = thmin+(thmax-thmin)*j/npoin 
       cth= dcos(th)                      
       if(ind.eq.5)then
          vv = PHEOCS(rr,cth)
       else 
          vv = V(rr,cth,ind)
       end if!(ind.eq.5)then
       call interpolacion(rr,cth,vlamda,plam,nlam,sumvlam)
       write(700,*) th*180._rk/pi,rr, vv,sumvlam 
   end do !j = 0, npoin                 
end do !i= 1, npoin                     
   
                                        
1000 format(*(f28.16,2x))               
close(100)                              
end program anisotropia_calculation                 
