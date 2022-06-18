module math                             
use types                            
                                        
implicit none                           
                                        
Contains 


  
Subroutine interpolacion(r0,th0,vlam, plam,nlam, vsum)                

real(rk), intent(in)                   :: th0, r0 
real(rk), allocatable,intent(in)       :: vlam(:,:), plam(:,:)
integer(ik),intent(in)                 :: nlam
real(rk), intent(out)                  :: vsum 
real(rk)                               :: pi=dacos(-1.0_rk) 
integer(ik)                            :: i, j, nColumna,lineath, linear!, lineacth
character(4)                           :: temp
real(rk)                               :: th1, th2, r1, r2!,cth !, th0, r0, vsum
real(rk), allocatable                  :: pl(:), vl(:) 

nColumna = nlam + 2 !2 significa que aun tengo otros 2 términos de la matriz que leer, el valor de R o th y el 0
allocate(pl(nColumna-1))        
allocate(vl(nColumna-1))        
                                        
!interpolación, encontrar los vlam, plam y sacar el potencial, para nuevos r0 y th0
!encontrar linea en donde está el valor theta y r
linear = 1 + nint((1000._rk*(r0-2._rk))/33._rk)
!lineath = 1 + nint((1000._rk*(th0))/pi)
lineath = 1 - nint(500._rk*(th0-1._rk))

!write(*,*) linear, lineath
   
   th1 = plam(lineath,1)
   th2 = plam(lineath+1,1)
   r1 = vlam(linear,1)
   r2 = vlam(linear+1,1)
   !write(*,*) "puntos de interpol",th0, th1, th2
   !write(*,*) "puntos de interpol", r0,r1,r2
   
   ! Asegura estar en la linea correcta
   if(abs(th0-th1)>(1._rk/1000._rk).or.abs(r0-r1)>(33._rk/1000._rk))then
   write(*,*) "no estas en la línea correcta"
   write(*,*) "th",lineath, th0, th1, th2
   write(*,*) "radio", linear,r0, r1, r2
   stop
   end if
   ! encuentra valores nuevos de vlam y plam
   do j = 2, nColumna
       pl(j-1) = plam(lineath,j) + ((th0-th1)/(th2-th1))*(plam(lineath+1,j) - plam(lineath,j))
       !write(*,*) pl(j-1)
   end do !j = 2, nColumna

   do j = 2, nColumna
       vl(j-1) = vlam(linear,j) + ((r0-r1)/(r2-r1))*(vlam(linear+1,j) - vlam(linear,j))
       !Multiplica los terminos Vlam por individual, empezando desde la columna 2 = lam0, j=2=lam0 
       !if( j == 4) then
       !    vl(j-1) = vl(j-1)/2._rk
       !end if
       !write(*,*) vl(j-1)
   end do !j = 2, nColumna
       vsum = sum(vl*pl)
       !write(*,*) "vsum", vsum
       !read *

1000 format(*(f28.16,2x))
end subroutine interpolacion                 
end module math                             
