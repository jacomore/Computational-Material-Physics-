program var_1d
implicit none 
!---------------------------------------------------------------------------
! var_1d solves the 1D Schrodinger equation for an infinite barrier well
! in the region [-1,1] (a.u), with a soft potential V0 in the region [a,b], 
! where a > -1 and b < 1. Eigenvectors and eigenvalues are calculated by 
! applying the variational principle. In particular, the WF is expanded 
! in a series of cos(kpix/2L (k odd) and sin(kpix/2L) (k even) and the 
! expansion coefficients are evaluated through Lapack subroutine "dsyev".           
!---------------------------------------------------------------------------

real(kind=8) :: V0,a,b          ! V0 is the inner potential in the region between [a,b]
integer :: N                    ! N is the number of basis functions          
integer ::  n_eig               ! number of eigenvalues to be computed
character(len=1) :: job         ! if "N" --> eigenvalues, if "V" --> eigenvalues + eigenvectors
 
real(kind=8), dimension(100) :: xdata, psi      ! Wavefunction visualization  
real(kind=8), parameter :: L = 1                ! width of either the sides of the well (-1,1)
real(kind=8), allocatable :: ham(:,:)           ! Hamiltonian
real(kind=8), allocatable :: W(:)               ! contains the eigenvalues in ascending order

integer :: i,j,ios, info        ! service variables


! Input unit
write(*,*) "Reading input from INPUT.dat"
open(unit=10, file="INPUT.dat", iostat = ios)
if (ios == 0) read(10,*,iostat = ios) N
if (ios == 0) read(10,*,iostat = ios) V0
if (ios == 0) read(10,*,iostat = ios) a
if (ios == 0) read(10,*,iostat = ios) b
if (ios == 0) read(10,*, iostat = ios) n_eig
if (ios == 0) read(10,*, iostat = ios) job
close(10)

! Allocate the Hamiltonian and the eigenvalues vector
allocate(ham(N,N), W(N))

! Calculating the hamiltonian of the system by calling the function H
write(*,*) "Calculating the hamiltonian of the system and storing it in ham.dat"
ham = H(V0,L,a,b,N)

! Writing the hamiltonian in "ham.dat"
open(unit=11, file = "ham.dat")
do i = 1,N
        write(11,*) (ham(i,j), j = 1,N)
enddo
close(11)

write(*,*) "Calculating eigenvalues and, optionally, eigenvectors"
! Calculating eigevalues and, optionally (if job = V), eigenvectors
call diagonal(ham,job,N,W,info)

write(*,*) "Storing eigenvectors in 'evalues.dat' file"
! Writing eigenvalues in the output file "evalues.txt"
open(unit=11, file = "evalues.dat")
do i = 1,n_eig
write(11,*) i,W(i)
enddo
close(11)


! Writing the first wavefunction, if eigenvectors are calculated
if (job == "V") then
        write(*,*) "Storing the square value of the wavefunction in 'wave1square.dat'"
        open(unit=11, file = "wave1square.dat")
        call wavefunc(V0,L,a,b,N,1,xdata,psi)
        do i = 1,100
                  write(11,*) xdata(i), psi(i)**2
                  enddo 
        close(unit=11)
endif

contains 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine linspace(xmin,xmax,n_points,xdata)
        implicit none 
        
        real(kind=8) :: xmin, xmax
        integer :: n_points,i

        real(kind=8), dimension(n_points) :: xdata
        real(kind=8) :: h         
        h = (xmax-xmin)/(n_points-1)
        do i = 0,n_points-1
                xdata(i+1) = xmin + h*i  
        enddo
        return
end subroutine 

!-------------------------------------------------------------------------
subroutine Romberg(a,b,k1,k2,L,Ival)
           implicit none 
           real(kind=8), intent(in) :: a,b,L
           integer, intent(in) ::  k1, k2

           integer :: row, col, k
           real(kind=8), dimension(0:n,0:n) :: R
           real(kind=8) :: h, sum_f

           real(kind=8), intent(out) :: Ival
           integer, parameter :: n = 5   ! number of romberg iterations

           Ival = 0
                
           h = 0.5*(b-a)
           R(0,0) = h*(basis_k12(b,k1,k2,L)+basis_k12(a,k1,k2,L))

           open(1,file='Romberg.txt')
           write(1,*) R(0,0)

           do row=1,n
           h = (b-a)/(2**row)        
           sum_f = 0
                do k=1,2**(row-1) 
                        sum_f = sum_f + basis_k12(a+(2*k-1)*h,k1,k2,L)
                enddo
                       
                R(row,0) = 0.5*R(row-1,0) + h*sum_f
                do col = 1,row  
                        R(row,col) = (4**col*R(row,col-1) - R(row-1,col-1))/(4.d0**col -1)             
                enddo
                write(1,*) R(row,:)
           enddo   
           write(1,*) k1,k2, R(n,n)
           Ival = R(n,n)
           return 
end subroutine Romberg
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
function V(V0,L,a,b,k1,k2)
        implicit none 
        real(kind=8) :: V0,a,b,L
        integer :: k1,k2
       
        real(kind=8) :: V, V_val

        call Romberg(a,b,k1,k2,L,V_val) 
        V = V0*V_val
        return 
end function V
!-----------------------------------------------------------------------------------------   


!-----------------------------------------------------------------------------------------
function basis_k12(x,k1,k2,L)
        implicit none 
        real(kind=8), intent(in) :: x
        integer, intent(in) :: k1, k2
        real(kind=8), intent(in) :: L 

        real(kind=8), parameter :: pi = acos(-1.0)
        real(kind=8) :: basis_k12

        if (mod(k1,2).ne. mod(k2,2)) then
                basis_k12 = 0.d0
        else
                basis_k12 = basis_k(x,k1,L)*basis_k(x,k2,L)
        endif
        return
end function basis_k12
!----------------------------------------------------------------------------------------


!----------------------------------------------------------------------------------------
function basis_k(x,k,L)
        implicit none 
        real(kind=8), intent(in) :: x
        integer, intent(in) :: k
        real(kind=8), intent(in) :: L 

        real(kind=8), parameter :: pi = acos(-1.0)
        real(kind=8) :: A
        real(kind=8) :: basis_k
        
        A = 1.d0/sqrt(L)
        if (mod(k,2) == 0) then 
                basis_k = A*sin(pi*k*x/(2.*L))
        elseif (mod(k,2) == 1) then
                basis_k = A*cos(pi*k*x/(2.*L))
        endif 
        return
end function basis_k

!----------------------------------------------------------------------------------------
function H(V0,L,a,b,N)
        implicit none
        integer :: N 
        real(kind=8) :: V0,L,a,b

        real(kind=8), dimension(N,N) :: H

        integer :: i,j 
        real(kind=8), parameter :: pi = acos(-1.0) 

        do i =1,N
                H(i,i) = (dble(i)*pi)**2/8 + V(V0,L,a,b,i,i)
                do j = 1,i-1
                        H(i,j) = V(V0,L,a,b,i,j)
                        H(j,i) = H(i,j)
                enddo
        enddo
        return 
end function H

!---------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------
subroutine diagonal(A,job,N,W,info)
        implicit none 
        integer, intent(in) :: N                ! Leading dimension of the matrix
        character(len=1), intent(in) :: job     ! If job = N --> eigenvalue only, if job=V --> eigen_value/vectors

        real(kind=8), dimension(N,N), intent(inout) :: A        ! Input: SYM real matrix, output: Eigenvalues
        real(kind=8), dimension(N) :: W         ! gives the eigenvalues in output
        integer, intent(out) :: info            ! inform whether the task was completed successfully or not
        integer :: LDA
        integer:: LWORK   
        real(kind=8), dimension(max(1,3*N-1)) :: WORK
        LDA = N 
        LWORK = max(1,3*N-1)
        call dsyev(job,"U",N,A,LDA,W,WORK,LWORK,info)
        return
end subroutine 
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------

subroutine wavefunc(V0,L,a,b,N,state,xdata,psi)
        implicit none
        integer, intent(in) :: N, state
        real(kind=8), dimension(100), intent(out) :: psi,xdata
        real(kind=8), intent(in) :: L,V0,a,b
        real(kind=8), dimension(N,N) :: A_matrix
        real(kind=8), dimension(N) :: W
        integer :: i,j

        call linspace(-1.d0,1.d0,100, xdata)
        
        A_matrix = H(V0,L,a,b,N)
        
        call diagonal(A_matrix,"V",N,W,info)
        
        psi = 0
        do i= 1,N
          do j = 1,100
             psi(j) = psi(j)+ basis_k(xdata(j),i,L)*A_matrix(i,state)
          enddo
        ! Writing the square value of the first wavefunction (charge density of the ground state)
        enddo
        return 
end subroutine wavefunc

end program



