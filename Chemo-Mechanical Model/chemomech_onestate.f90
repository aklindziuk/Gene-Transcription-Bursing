program twostate_ssa

implicit none

integer :: n, m, i, j, dim
real :: time, t_max, random1, random2, a, b, kn, kf, y,alpha_0, alpha_1, alpha_2
real :: dist_sum, norm_sum, tau,  mean, m2, var, p_on, p_off, pf_check
real :: t
real, dimension(1:2, 0:500) :: dist
real, dimension(1:2, 0:500) :: dist_norm
real, dimension(0:500) :: histogram
real :: start, finish


call cpu_time(start)

! ---------------select parameters----------------

a = 10
b = 2
kn = 100000
kf = 10
y = 1

!--------------------------------------------------

t_max = 10000
dim = 500

call random_seed

dist = 0
histogram = 0
dist_norm = 0
t = 0

	n = 0
	m = 1
	j=0
	time = 0.0

	open(2, file = "timeseries.dat", status = "replace")

	do while (time < t_max) 
		
		call random_number (random1)
		!because 0 random will make tau = Infinity
		if (random1 == 0.0) then
			call random_number (random1)
		end if
		call random_number (random2)
		
		if ( m == 1) then
			if ( n == 0) then

				!!! case 1 ----------> ON (n=0)

				alpha_0 = kf + a 
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
					t = t + tau
					write(2,*) time, 'case 1', m, n, j, tau
				end if
				time = time + tau
				
				if ( random2 < a/alpha_0) then
					n = n + 1
				else
					m = m + 1
				end if

			else

				!!! case 2 ------------> ON (n>0)

				alpha_0 = kf + a + n*b
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
					t = t + tau
					write(2,*) time, 'case 2', m, n, j, tau, random1
				end if
				time = time + tau

				alpha_1 = kf/alpha_0
				alpha_2 = (kf + a)/alpha_0
				if (random2 < alpha_1) then !(if kf)
					m = m + 1 
				else if (alpha_1 < random2 .and. random2 < alpha_2) then !(if a)
					n = n + 1
				else !(if nb)
					n = n - 1 
				end if

			end if
		else 
			if ( n == 0) then

				!!! case 3 ----------> OFF (n=0)

				alpha_0 = kn + a/(y**(j+1))
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
					t = t + tau
					write(2,*) time,'case 3', m, n, j, tau, random1
				end if
				time = time + tau

				alpha_1 = kn/alpha_0
				if (random2 < alpha_1) then !(if kn)
					m = m - 1
					j = 0 
				else 
					n = n + 1
					j = j + 1
				end if
			else

				!!! case 4 ---------------> OFF (n>0)

				alpha_0 = kn + n*b + a/(y**(j+1))
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
					t = t + tau
					write(2,*) time,'case 4', m, n, j, tau, random1
				end if
				time = time + tau

				alpha_1 = kn/alpha_0
				alpha_2 = (kn + n*b)/alpha_0
				if ( random2 < alpha_1) then !(if kn)
					m = m - 1
					j = 0
				else if (random2>alpha_1 .and. random2<alpha_2) then !(if nb)
					n = n - 1
				else !(if a/y**j)
					n = n + 1
					j = j + 1
				end if 

			end if
		end if
end do

!--------------------note -------------------------------
!   if time is set to a large value, rounding error seems to occur, 
!   which makes dist_sum >> (time - time*.2)... strange
!--------------------------------------------------------

dist_sum = 0
norm_sum = 0
p_on = 0
p_off = 0
pf_check = 0
m2 = 0
var = 0
mean = 0

!------------find sum of times in the distribution ----------------

do j = 1, 2
 do i = 0, dim
 	dist_sum = dist_sum + dist(j,i)
 end do
end do
print*, " "
print*, "dist sum", dist_sum

!------------normalize distribution, check sum = 1 -------------------

do j = 1,2
	do i = 0, dim
		dist_norm(j,i) = dist(j,i)/dist_sum
		norm_sum = norm_sum + dist_norm(j,i)
	end do
end do
print*, "norm sum", norm_sum

!------------------------find the histogram -------------------------

histogram(:) = dist_norm(1,:) + dist_norm(2,:)

!-------------------find moments, variance-------------------------------

do i = 0,dim
	mean = mean + i*histogram(i)
	m2 = m2 + i**2*histogram(i)
	var = var + (i - mean)**2*histogram(i)
end do 
print*, " "
print*, 'mean ', mean
print*, '<n^2>', m2
print*, 'var  ', var

!------------------find prob in each state-----------------------------

do i = 0,dim
	p_on = p_on + dist_norm(1,i)
	p_off = p_off + dist_norm(2,i)
end do

print*, " "
print*, 'P_on', p_on
print*, 'P_off', p_off



!------------------write all the distributions----------------------

open(3, file = "dist.dat", status = "replace")
open(4, file = 'dist_norm.dat', status = 'replace')
open(5, file = 'histogram.dat', status = 'replace')

do i = 0, dim
 	write(3, *)  i, dist(1,i), dist(2,i)
 	write(4, *)  i, dist_norm(1,i), dist_norm(2,i)
 	write(5, *)  histogram(i)
end do
close(3)
close(4)
close(5)


end program twostate_ssa
