program twostate_ssa

implicit none

integer :: n, m, i, j, k, dim
real :: time, t_max, g, e, x,random1, random2, a, b, kn, kf, alpha_0, alpha_1, alpha_2, mean, mean2, var, var2, sum, tau, fon, dnt
real :: nn, ff, on, phi, y, mean_theory, fano_theory, n2_epsilon, gamma, epsilon, var_epsilon, f2_epsilon, part2, n_epsilon
real :: fano, fano_yone, fano_yinf, f2_yone, f2, q1, q2, r1, r2, n2_yone,variance_yone, n_yone, fano_yone_simplified
real, dimension(1:2, 0:10000) :: dist
real, dimension(1:2, 0:10000) :: dist_norm
real, dimension(0:10000) :: histogram
real :: start, finish


call cpu_time(start)

a = 10
b = 1
kn = .001
kf = .001
e =.01
y = 1+e

gamma = kf/kn
x = a/b
g = kf/kn
epsilon = e

t_max = 200000
dim = 10000

call random_seed

dist_norm = 0
dist = 0


	n = 0
	m = 2

	time = 0.0

	open(2, file = "timeseries.dat", status = "replace")

	do while (time < t_max) 
		
		call random_number (random1)
		call random_number (random2)
		
		if ( n == 0) then
			if ( m == 1) then

				!!! case 1 ---------->

				phi = kn + a/y
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				
				if ( random2 < (a/y)/phi) then
					n = n + 1
				else
					m = m + 1
				end if
			else  

				!!! case 2 ------------>

			    phi = kf + a + a/y
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if (random2 < (a + a/y)/phi) then
					n = n + 1
				else
					m = m - 1
				end if
			end if
		else
			if ( m == 1) then

				!!! case 3 ---------->

				phi = kn + a/(y**(n+1)) + n*b
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				if (random2 < kn/phi) then
					m = m + 1 
				else if (kn/phi < random2 .and. random2 < (kn+a/(y**(n+1)))/phi) then
					n = n + 1
				else 
					n = n - 1 
				end if
			else

				!!! case 4 --------------->

				phi = kf + a + a/(y**(n+1)) + n*b
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				if ( random2 < kf/phi) then
					m = m - 1
				else if (kf/phi < random2 .and. random2 < (kf + a + a/(y**(n+1)))/phi) then 
					n = n + 1
				else
					n = n - 1
				end if 

			end if
		end if

	write(2,*) time,",", m,",", n


	if (n .eq. dim) then
		print*, "n is out of bounds"
	end if


	end do

sum= 0
mean = 0
mean2 = 0
var2 = 0
fon = 0
dnt = 0

! total time
open(3, file = "dist.dat", status = "replace")
do j = 1, 2
 do i = 0, dim
 	sum = sum + dist(j,i)
 	write(3, *) j,",", i, ",", dist(j,i)
 end do
end do
close(3)

! normalize distribution, find mean
open(1, file = "dist_norm.dat", status = "replace")
do j = 1,2
 do i = 0, dim
 	dist_norm(j,i) = real(dist(j,i))/sum
 	write(1, *) i, ",", dist_norm(j,i)
 	mean = mean + i*dist_norm(j,i)
 end do
end do
close(1)

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

open(4, file = "histogram_twornap.dat", status = "replace")
do i = 0,dim
		write(4, *) i, ",", histogram(i)
end do
close(4)


!!! ====== GENERAL RESULTS =======================================
mean_theory = a/b*(1/y*(gamma)/(1+gamma) + (1 + 1/y)/(1+gamma))


!!! ====== Y=1 APPROXIMATION ===================================
f2_yone = a/b*(gamma/(1+ gamma)*(kn + kf)+ 1/(1+gamma)*(b + b*gamma + kn + kf) + b + kn)/(1+ gamma)/(b + kn+ kf)
n_yone = a/b*(2+gamma)/(1+gamma)
n2_yone = a/b*((2+gamma)/(1+gamma)+f2_yone+n_yone)
variance_yone = n2_yone - n_yone**2
fano_yone = (1/(1+gamma)+ f2_yone + n_yone)*((1+gamma)/(2+gamma))-a/b*(2+gamma)/(1+gamma)
fano_yone_simplified = 1 + a/b*(1-(2+gamma)/(1+gamma))+a/b*(2*b+2*b*gamma+2*kn+2*kn*gamma+kf+kf*gamma)/(1+gamma)/(2+gamma)/(b+kn+kf)

!!! ===== EPSILON APPROXIMATION ==================================
n_epsilon = a/b*((2+gamma)/(1+gamma)-epsilon*(1+a/b*(2+gamma)/(1+gamma)))
part2 = x*(b+kn)/(1+g)/(b+kn+kf)-e*x**2*b*(g/(1 + g)*(kn+kf)+1/(1+g)*(kn+kf+b+b*g)+b+kn)/(1 + g)/(b + kn + kf)**2
f2_epsilon =x*(1-e)*(g/(1+g)*(kn+kf)+1/(1+g)*(kn+kf+b+b*g)-e*x*(kn + kf)*(2+g)/(1+g))/(1+g)/(b+kn+kf)+part2
n2_epsilon = x*((2+g)/(1+g)-e*(1+x*(2+g)/(1+g))+f2_epsilon+x*(2+g)/(1+g)*(1-2*e))
var_epsilon = n2_epsilon - n_epsilon**2

!!! ==== Y=INF APPROXIMATION ======================================
fano_yinf = 1 + (b + kn)/(b + kn+ kf) - a/b/(1 + gamma)


call cpu_time(finish)

print*, "     ","simulaton       ", 'analytical      '
print*, "GENERAL RESULTS:----------------------"
print*, "mean =  ", mean, mean_theory, (mean - mean_theory)
print*, "fano =", var2/mean, fano, (var2/mean - fano_theory)
print*, ""
print*, "Y=1 APPROXIMATION---------------------"
print*, "fano_y=1", var2/mean, fano_yone_simplified
print*, ""
print*, "Y=INF APPROXIMATION-------------------"
print*, "fano_y=inf", var2/mean, fano_yinf
print*,""
print*, "EPSILON APPROXIMATION-----------------"
print*, "n_epsilon", mean, n_epsilon
print*, "var_epsilon", var2, var_epsilon
print*, "fano_epsilon", var2/mean, var_epsilon/n_epsilon

print*, '-------------------'

print*, 'total time', sum
print*, 'check norm time', dnt
!print *,"simul time =", finish-start



end program twostate_ssa