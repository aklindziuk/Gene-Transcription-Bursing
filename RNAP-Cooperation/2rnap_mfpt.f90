program twostate_ssa

implicit none

integer :: n, m, i, j, k, dim, counter, c, h, nn, t
real :: time, t_max, g, e, x,random1, random2, a, b, kn, kf, alpha_0, alpha_1, alpha_2, mean, mean2, var, var2, sum, tau, fon, dnt
real :: y, phi,mean_theory, fano_theory, n2_epsilon, gamma, epsilon, var_epsilon, f2_epsilon, part2, n_epsilon
real :: fano, fano_yone, fano_yinf, f2_yone, f2, q1, q2, r1, r2, n2_yone,variance_yone, n_yone, fano_yone_simplified
real :: target, fpmean, O, Q, fact1, fact2
real, dimension(1:2, 0:10000) :: dist
real, dimension(1:2, 0:10000) :: dist_norm
real, dimension(0:10000) :: histogram
real, dimension(0:500) :: fp_time
real, dimension(0:10) :: mfpt
real :: start, finish


call cpu_time(start)

a = 10
b = 1
kn = 1
kf = 1
y = 1.0

x = a/b

gamma = kf/kn

t_max = 300000
dim = 10000

fp_time = 0

mfpt = 0


call random_seed

open(2, file = "timeseries.dat", status = "replace")

open(3, file = "mfpt.dat", status = "replace")

do t = 1, 10 !================iterate through dirrerent targets

	target = real(t)

do c = 1, size(fp_time) !=========== iterate through 100 rounds 
	
	n = 0
	m = 1
	if (c .ge. int(size(fp_time)/2)) then 
		m = 2
	end if


	dist_norm = 0
	dist = 0

	time = 0.0
	counter = 1

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

	if (n == target) then
		if (counter == 1) then
			fp_time(c) = time
		    !print*, c, fp_time(c)
		end if
		counter = counter + 1
	end if

	if (c==5) then
		write(2,*) time, m, n
	end if


	if (n .eq. dim) then
		print*, "n is out of bounds"
	end if


	end do !========= one round over ============

if (counter == 1) then 
	print*, "target not reached in ", t_max, "time units"
	fp_time(c) = t_max
end if 

end do !======== all rounds over ====================

open (5, file = "fptime.dat", status = "replace")
do i = 1, size(fp_time)
	write(5,*), i, fp_time(i)
end do
close(5)

fpmean = 0
do i = 1,size(fp_time)
	fpmean = fpmean + fp_time(i)
end do
fpmean = fpmean/size(fp_time)
!print*, kfs(h)/kns(h), fpmean, target
write(3,*), target,",", fpmean
print*, target, "done"

end do ! =======================finished with all targets-===========

close(3)



call cpu_time(finish)




end program twostate_ssa