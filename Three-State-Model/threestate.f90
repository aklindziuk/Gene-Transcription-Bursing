!---------------------------------
! Alena Klindziuk
! 06/01/2018
!
! Stochastic simulation based on Gillespie algorithm 
! Simulates the dynamics of a gene transcription system with three states
! Records the probability of visiting any allowed state in the system
!
! Input:
! mRNA generation rate a, degeneration rate b, transition rates kf and kn
!
! Output: 
! - command line: mean, variance, fano factor from simulation vs analytically calc
! -"timeseries. dat" (time, m, n) records every state visited and the time when it !  happened
! -"dist.dat" (n, m=1 bin, m=2 bin, m=3 bin) records how many times each n, m was !  visited by the system
! -"dist_norm.dat" (n, normalized bins) records probibility distribution  
!  of system visiting each state
! -"histogram.dat" (n, probibility dens) probibility density for each n, summing 
!  over m=1, m=2, m=3 for each n
!
! To run: gfortran threestate.f90 -o threestate.exe
! To execute: ./threestate.exe   
!---------------------------------


program threestate_ssa

implicit none

integer :: n, m, i, j, k, dim_x, dim_y
real :: a1, a2, a3, b, kn2, kn3, kf1, kf2, c0
real :: time, tau, t_max,rand1, rand2, mean, mean2, var, var2, sum
real :: frac_1, frac_2, frac_3, dnt
real :: analytic_mean, analytic_fano, analytic_var, intermediate 
real :: i1, i2, i3, i4, i5, i6, i7, i8
real, dimension(1:3, 0:1000) :: dist
real, dimension(1:3, 0:1000) :: dist_norm
real, dimension(0:1000) :: histogram
real :: start, finish
real :: x1, x2, x3, y1, y2

call cpu_time(start)

!---------select system parameters-----------------
a1 = 100
a2 = 100
a3 = 100
b = 1
kn2 = .10
kn3 = .10
kf1 = 10
kf2 = 10
!-----------------------------------------------------

! run simulation for this long
t_max = 10000

call random_seed

dim_x = 1000
dim_y = 3

dist_norm = 0
dist = 0
histogram = 0

!set initial state to n=0, m=1
n = 0
m = 1

!start simulation at time = 0
time = 0.0


	open(2, file = "timeseries.dat", status = "replace")

	do while (time < t_max) 
		
		call random_number (rand1)
		call random_number (rand2)
		
		if (rand1 == 0.0) then
		call random_number(rand1)
		end if 
		
		if ( n == 0) then
			if ( m == 1) then

				!!! case 1 ---------->
				c0 = kf1 + a1
				tau = (1/c0)*log(1/rand1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if ( rand2 < a1/c0) then
					n = n + 1
				else
					m = m + 1
				end if

			else if ( m == 2 ) then

				!!! case 2 ------------>
				c0 = a2 + kf2 + kn2
				tau = (1/c0)*log(1/rand1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if ( rand2 < kf2/c0 ) then
					m = m - 1
				else if ( (rand2 > kf2/c0) .and. (rand2 < (kf2+a2)/c0) ) then
					n = n + 1
				else 
					m = m + 1
				end if
			
			else 

				!!! case 3 -------------->
				c0 = kn3 + a3
				tau = (1/c0)*log(1/rand1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if ( rand2 < kn3/c0 ) then
					m = m - 1
				else 
					n = n + 1
				end if

			end if

		else

			if ( m == 1) then

				!!! case 4 ---------->
				c0 = n*b + a1 + kf1
				tau = (1/c0)*log(1/rand1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if ( rand2 < n*b/c0 ) then
					n = n - 1
				else if ( (rand2 > n*b/c0) .and. (rand2 < (n*b + kf1)/c0) ) then
					m = m + 1
				else
					n = n + 1
				end if

			else if ( m == 2 ) then

				!!! case 5 --------------->
				c0 = kn2 + n*b + a2 + kf2
				tau = (1/c0)*log(1/rand1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if ( rand2 < kn2/c0 ) then
					m = m - 1
				else if ( (rand2 > (kn2)/c0) .and. (rand2<(kn2+n*b)/c0) ) then
					n = n - 1
				else if ( (rand2 > (kn2+n*b)/c0) .and. (rand2<(kn2+n*b+a2)/c0)) then
					n = n + 1
				else
					m = m + 1
				end if

			else

				!!! case 6 -------------->
				c0 = kn3 + n*b + a3
				tau = (1/c0)*log(1/rand1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if ( rand2 < kn3/c0 ) then
					m = m - 1
				else if ( (rand2>kn3/c0) .and. (rand2<(kn3+n*b)/c0) ) then
					n = n - 1
				else
					n = n + 1
				end if

			end if
		end if

	write(2,*) time, m, n

	end do
	close(2)

sum= 0
mean = 0
mean2 = 0
var2 = 0
frac_1 = 0
dnt = 0

! total time
do i = 1, dim_y
 do j = 0, dim_x
 	sum = sum + dist(i,j)
 end do
end do

! normalize distribution, find mean
do  i = 1, dim_y
 do j = 0, dim_x
 	dist_norm(i,j) = real(dist(i,j))/sum
 	mean = mean + real(j)*dist_norm(i,j)
 end do
end do

! find variance, check sum of normalized bins = 1
do i = 1, dim_y
	do j = 0, dim_x
		var2 = var2 + dist_norm(i,j)*(real(j)-mean)**2
		dnt = dnt + dist_norm(i,j)
	end do 
end do

do i = 1, dim_y
 histogram(:) = histogram(:) + dist_norm(i,:)
end do

!find mean from hist, find fraction in each state
do i = 0, dim_x
	mean2 = mean2 + histogram(i)*real(i)
	frac_1 = frac_1 + dist_norm(1, i)
	frac_2 = frac_2 + dist_norm(2, i)
	frac_3 = frac_3 + dist_norm(3, i)
end do 

! write all findings to file
open(1, file = "dist_norm.dat", status = "replace")
open(3, file = "dist.dat", status = "replace")
open(4, file = "histogram.dat", status = "replace")
do i = 0,dim_x
		write(1, *) i, dist_norm(1,i), dist_norm(2,i), dist_norm(3,i)
		write(3, *) i, dist(1,i), dist(2,i), dist(3,i)
		write(4, *) i, histogram(i)
end do
close(4)
close(1)
close(3)


! ---- calculate mean, var, and fano analytically-----

x1 = a1/b
x2 = a2/b
x3 = a3/b

y1 = kf1/kn2
y2 = kf2/kn3

analytic_mean = ( x1 + x2*y1 + x3*y1*y2 )/(1 + y1 + y2*y1 )
i1 = x1*x1*(b**2 + b*kf2 + b*kn2 + b*kn3 + kf2*kn3) 
i2 = 2*x1*x2*(b*kf1 + kf1*kn3 ) 
i3 = 2*x1*x3*(kf1*kf2)     
i4 = 2*x2*x3*y1*(b*kf2 + kf1*kf2) 
i5 = x2*x2*y1*(b**2 + b*kf1 + b*kn3 + kf1*kn3) 
i6 = x3*x3*y2*(b**2*y1 + b*kf1*y1 + b*kf2*y1 + kf1*kf2*y1 + b*kf1)
i7 = (1 + y1 + y1*y2) 
i8 = (b**2 + b*kf1 + b*kf2 + kf1*kf2 + b*kn2 + b*kn3 + kf1*kn3 + kn2*kn3)
intermediate = (i1+ i2 + i3 + i4 + i5 + i6)/(i8*i7)
analytic_fano = 1 + intermediate/analytic_mean - analytic_mean
analytic_var = analytic_mean + intermediate - analytic_mean**2

call cpu_time(finish)

print*, '-------------------'
print*, "results:   ","simulaton       ", 'analytic        ','error'
print*, "mean =  ", mean, analytic_mean, (mean - analytic_mean)
print*, "fano =  ", var2/mean, analytic_fano, (var2/mean - analytic_fano)
print*, "variance", var2, analytic_var, (var2 - analytic_var)
print*, "-------------------"
print*, "frac_1 =  ", frac_1
print*, "frac_2 =  ", frac_2
print*, "frac_3 =  ", frac_3
print*, 'total time', sum
print*, 'check norm time', dnt
!print *,"simul time =", finish-start

end program threestate_ssa


