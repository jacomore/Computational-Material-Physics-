program Schro_1D_Shoot
implicit none

! Allocating variables
integer :: N                                 ! number of subdivisions
integer :: i,j,k                             ! loop integers
 

real(kind = 8) , parameter :: eps = 1.d-5    ! machine epsilon
real(kind = 8) :: Dx                         ! width of the x partition
real(kind = 8) :: De                         ! Energy increase
real(kind = 8) :: E                          ! energy of the well level
real(kind = 8) :: xmin, xmax                 ! edges of the well
real(kind = 8) :: psi_old, psi_new           ! the two values of the WF in N+1 positions
real(kind = 8) :: sign_prod, A

real(kind = 8), allocatable :: psi(:)        ! wavefunction
real(kind = 8), allocatable :: phi(:)        ! phi array defined by the Numerov's algorithm



! I/O unit

write(*,*) "Number of subdivisions:"
read(*,*) N
write(*,*) "Energy increase:"
read(*,*) De
write(*,*) "left bound of the partition:"
read(*,*) xmin
write(*,*) "right bound of the partition:"
read(*,*) xmax
write(*,*) "Initial energy of the system:"
read(*,*) E

! allocate the wavefunction psi and the auxiliary function phi
allocate(psi(N+1), phi(N+1))

! define the width of the partition intervals
Dx = (xmax-xmin)/N

! calculate all the elements of phi and only the last point of psi
call Numerov(Dx, E, N, phi, psi_old)

! perform the do loop while the sign of psi_old*psi_new > 0
print*, "Course search"
sign_prod = 1.d0
do while (sign_prod > 0.d0)
        E = E + De
        print*, "Energy:", E
        call Numerov(Dx, E, N, phi, psi_new)
        sign_prod = psi_old*psi_new
        psi_old = psi_new
end do 

!The correct value of energy is contained within the interval [E, E-De]

print*, "Start bisection:"
print*, "The wavefunction is store in the file psi.dat"
call bisection(E-De, E, Dx, N, eps, E, psi)

A = norm(psi,Dx,N)


open(unit=10, file = "psi.dat", action = "WRITE")
do i=1,N+1
        write(10,*) xmin +  (i-1)*Dx, A*psi(i)
enddo 
close(unit=10)
              
! Functions and subroutines

contains 
!----------------------------------------------------------
        function psi_to_phi(psi,E,Dx)
               implicit none 

               real(kind=8) :: psi, E, Dx
               real(kind=8) :: psi_to_phi
              
               psi_to_phi = psi*(1 - Dx**2*(-2*E)/12)
      
               return
       end function psi_to_phi
!----------------------------------------------------------

!----------------------------------------------------------
         function phi_to_psi(phi,E, Dx)
                 implicit none 

                 real(kind=8) :: phi, E, Dx
                 real(kind=8) :: phi_to_psi

                 phi_to_psi = phi/(1- Dx**2*(-2*E)/12)
                
                 return 
         end function phi_to_psi
!-----------------------------------------------------------

!-----------------------------------------------------------
          subroutine Numerov(Dx,E, N, phi, psi_last)
                 implicit none 

                 real(kind=8), intent(in) :: Dx, E
                 integer, intent(in):: N

                 real(kind=8), dimension(N+1), intent(out) :: phi
                 real(kind=8), intent(out) :: psi_last

                 phi(1) = psi_to_phi(0.d0,E,Dx)
                 phi(2) = psi_to_phi(1.d0,E,Dx)

                 do i=2,N
                        phi(i+1) = 2*phi(i) - phi(i-1) + Dx**2*(-2*E)*phi_to_psi(phi(i),E,Dx)
                 enddo

                 psi_last = phi_to_psi(phi(N+1),E,Dx)

                 return 
           end subroutine
!----------------------------------------------------------

!----------------------------------------------------------
          subroutine bisection(E_min, E_max, Dx, N, eps, E_best, psi)
                  implicit none 
                  
                  real(kind=8), intent(in) :: Dx, eps
                  real(kind=8), intent(in) :: E_max, E_min
                  integer, intent(in) :: N

                  real(kind=8), intent(out) :: E_best
                  real(kind=8), dimension(N+1), intent(out) :: psi

                  real(kind=8), dimension(N+1) :: phi
                  real(kind=8) :: psi_last_avg, psi_last_min, psi_last_max
                  real(kind=8) :: E_avg, E_max_tmp, E_min_tmp

                  E_max_tmp = E_max
                  E_min_tmp = E_min

                  do    
                        E_avg = (E_min_tmp + E_max_tmp)/2.d0

                        call Numerov(Dx, E_min_tmp, N, phi, psi_last_min)
                        call Numerov(Dx, E_max_tmp, N, phi, psi_last_max)
                        call Numerov(Dx, E_avg, N, phi, psi_last_avg)
                        

                        print*, "Average energy:",E_avg

                        if (abs(psi_last_avg) < eps) then 
                                do i=1,N+1
                                        psi(i) = phi_to_psi(phi(i), E_avg, Dx)
                                enddo 
                                E_best = E_avg
                                return

                        elseif (psi_last_min*psi_last_avg < 0) then 
                               E_max_tmp = E_avg
                               
                       
                       elseif (psi_last_max*psi_last_avg < 0) then 
                               E_min_tmp = E_avg
                              
                     
                       end if
                  enddo
          end subroutine bisection

!------------------------------------------------------------------------------
          function norm(psi,Dx,N)
                implicit none 
                 
                integer :: N 
                real(kind=8) :: Dx
                real(kind=8), dimension(N+1) :: psi, psi_square

                real(kind=8) :: norm,integral
                
                ! #1: calculate the integral of the absolute value of the wavefunction by means of the midpoint rule for quadrature

                integral = 0
                psi_square(:) = psi(:)**2   ! coherently with the assumption that the WF is real. In general psi should be complex
                
                do i=1,N
                        integral = integral + Dx*(psi_square(i) + psi_square(i+1))/2
                enddo
                norm = 1/sqrt(integral)
                
                return
          end function norm           



                               





               
                  















            
            



end program 
