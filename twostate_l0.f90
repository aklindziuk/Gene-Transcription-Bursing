program twostate_ssa

implicit none

integer :: n, m, i, j, k, dim
real :: time, t_max, random1, random2, a, b, kn, kf, alpha_0, alpha_1, alpha_2, mean, mean2, var, var2, sum, tau, fon, dnt
real :: nn, ff, on
real, dimension(1:2, 0:250) :: dist
real, dimension(1:2, 0:250) :: dist_norm
real, dimension(0:250) :: histogram
real :: start, finish


call cpu_time(start)

! I
a = 100
b = 16
kn = 1
kf = 3

t_max = 10000
dim = 250

call random_seed

dist_norm = 0
dist = 0

	n = 0
	m = 1
	time = 0.0

	open(2, file = "timeseries.dat", status = "replace")

	do while (time < t_max) 
		
		call random_number (random1)
		call random_number (random2)
		
		if ( n == 0) then
			if ( m == 1) then

				!!! case 1 ---------->

				alpha_0 = kf + a
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				
				if ( random2 < a/alpha_0) then
					n = n + 1
				else
					m = m + 1
				end if
			else

				!!! case 2 ------------>

				tau = (1/kn)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
		
				m = m - 1
			end if
		else
			if ( m == 1) then

				!!! case 3 ---------->

				alpha_0 = a + n*b
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				alpha_1 = a/alpha_0
				if (random2 < alpha_1) then
					n = n + 1
				else 
					n = n - 1 
				end if
			end if
		end if

	write(2,*) time, m, n

	end do

sum= 0
mean = 0
mean2 = 0
var2 = 0
fon = 0
dnt = 0

! total time
do j = 1, 2
 do i = 0, dim
 	sum = sum + dist(j,i)
 end do
end do

! normalize distribution, find mean
do j = 1,2
 do i = 0, dim
 	dist_norm(j,i) = real(dist(j,i))/sum
 	mean = mean + i*dist_norm(j,i)
 end do
end do

! find variance, check normalized time sum
do j = 1, 2
	do i = 0, dim
		var2 = var2 + dist_norm(j,i)*(real(i)-mean)**2
		dnt = dnt + dist_norm(j, i)
	end do 
end do

histogram(:) = dist_norm(1,:) + dist_norm(2,:)

!find mean, find fraction on
do i = 1, dim
	mean2 = mean2 + histogram(i)*i
	fon = fon + dist_norm(1, i)
end do 

open(1, file = "dist_norm.dat", status = "replace")
open(3, file = "dist.dat", status = "replace")
open(4, file = "histogram.dat", status = "replace")
do i = 0,dim
		write(1, *) i, dist_norm(1,i), dist_norm(2,i)
		write(3, *) i, dist(1,i), dist(2,i)
		write(4, *) i, histogram(i)
end do
close(4)
close(1)
close(3)

nn = (a/b)/(1+(kf/kn)*exp(-a/b))
ff = 1 + ((a/b)*(kf/kn)*exp(-a/b))/(1+(kf/kn)*exp(-a/b))
on = kn/(kn + kf)

call cpu_time(finish)

print*, "results:   ","simulaton       ", 'analytical      ','error'
print*, '-------------------'
print*, "mean =  ", mean, nn, (mean-nn)
print*, "fano2 =", var2/mean, ff, (var2/mean - ff)
!print*, "fon =  ", fon, on, (fon- on)
print*, '-------------------'
print*, "variance2", var2
print*, 'total time', sum
print*, 'check norm time', dnt
!print *,"simul time =", finish-start



end program twostate_ssa


