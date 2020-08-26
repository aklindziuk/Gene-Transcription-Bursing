program twostate_ssa

implicit none

integer :: n, m, i, j, k, dim, counter, c, h, nn, initail_m, t, rounds
real :: time, t_max, g, e, x,random1, random2, a, b, alpha_0, alpha_1, alpha_2, mean, mean2, var, var2, sum, tau, fon, dnt
real :: y,kf, kn, phi,mean_theory, fano_theory, n2_epsilon, gamma, epsilon, var_epsilon, f2_epsilon, part2, n_epsilon
real :: fano, fano_yone, fano_yinf, f2_yone, f2, q1, q2, r1, r2, n2_yone,variance_yone, n_yone, fano_yone_simplified
real :: target, fpmean, O, Q, fact1, fact2, SD
real, dimension(1:2, 0:10000) :: dist
real, dimension(1:2, 0:10000) :: dist_norm
real, dimension(0:10000) :: histogram
!real, dimension(0:100) :: fp_time
real, dimension (:,:), allocatable :: fp_time 
real, dimension(1:10) :: ys
real :: start, finish
character(len=10) :: file_id
character(len=50) :: file_name
character(20) filename


call cpu_time(start)

rounds = 700

!ys

ys(1) = 1.0
ys(2) = 1.01
ys(3) = 1.03
ys(4) = 1.07
ys(5) = 1.15
ys(6) = 1.31
ys(7) = 1.63
ys(8) = 2.27
ys(9) = 3.55
ys(10) = 6.11


a = 10
b = 1
kn = .1
kf = .1
x = a/b

allocate ( fp_time(1:rounds,1:size(ys)) )    

t_max = 300000
dim = 10000

fp_time = 0

call random_seed


do t = 4, 5 !================iterate through dirrerent targets

	target = real(t)

	write(filename,'(a,i0,a)') 'target', t,'.txt'
	open(3, file = filename , status = "replace")
	write(3,*), "target,","fpmean,", "SD"

	write(filename,'(a,i0,a)') 'fp_time', t,'.txt'
	open(5, file = filename , status = "replace")
	write(5,*), "ys",ys(1),ys(2),ys(3),ys(4),ys(5),ys(6),ys(7),ys(8),ys(9)!,ys(10)

do h = 1, size(ys) !========iterate through different gammas 

	!write(filename,'(a,i0,a,f3.1,a)') 'fptime',t,'-', ys(h),'.txt'
	!open(5, file = filename , status = "replace")

do c = 1, rounds !=========== iterate through 100 rounds 
	
	n = 0
	m = 1
	if (c .ge. int(rounds/2)) then 
		m = 2
	end if

	!Special cases
	if (kf/kn < .001) then 
		m = 2
	end if
	if ( kf/kn > 1000) then
		m = 1
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

				phi = kn + a/ys(h)
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				
				if ( random2 < (a/ys(h))/phi) then
					n = n + 1
				else
					m = m + 1
				end if
			else  

				!!! case 2 ------------>

			    phi = kf + a + a/ys(h)
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau
				if (random2 < (a + a/ys(h))/phi) then
					n = n + 1
				else
					m = m - 1
				end if
			end if
		else
			if ( m == 1) then

				!!! case 3 ---------->

				phi = kn + a/(ys(h)**(n+1)) + n*b
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				if (random2 < kn/phi) then
					m = m + 1 
				else if (kn/phi < random2 .and. random2 < (kn+a/(ys(h)**(n+1)))/phi) then
					n = n + 1
				else 
					n = n - 1 
				end if
			else

				!!! case 4 --------------->

				phi = kf + a + a/(ys(h)**(n+1)) + n*b
				tau = (1/phi)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				if ( random2 < kf/phi) then
					m = m - 1
				else if (kf/phi < random2 .and. random2 < (kf + a + a/(ys(h)**(n+1)))/phi) then 
					n = n + 1
				else
					n = n - 1
				end if 

			end if
		end if

	if (n == target) then
		if (counter == 1) then
			fp_time(c, h) = time
		    !print*, c, fp_time(c)
		end if
		counter = counter + 1
	end if

	if (c==5) then
		!write(2,*) time, m, n
	end if


	if (n .eq. dim) then
		print*, "n is out of bounds"
	end if

	!print*,  y, c


	end do !========= one round over ============

if (counter == 1) then 
	print*, "target not reached in ", t_max, "time units"
	fp_time(c, h) = t_max
end if 

!write(5,*), c, ",", fp_time(c)


end do !======= all 100 rounds over ====================

print*, "done with y=", ys(h)

!close(5)


fpmean = 0
do i = 1,rounds
	fpmean = fpmean + fp_time(i, h)
end do
fpmean = fpmean/real(rounds)

var = 0
SD = 0
do i = 1, rounds
	var = var + (fp_time(i,h)-fpmean)**2
end do
SD = sqrt(var/(real(rounds)-1))

!print*, "var", var
!print*, "size fp_time", rounds
!print*, "SD", SD

!print*, ys(h), fpmean, target
write(3,*), target, ",", ys(h), ",", fpmean, ",", SD


end do ! =======================finished with all ys-===========

close(3)
do i = 1, rounds
	write(5, *), i,",",fp_time(i,1),",",fp_time(i,2),",",fp_time(i,3),",",fp_time(i,4),",", &
	fp_time(i,5),",",fp_time(i,6),",",fp_time(i,7),",",fp_time(i,8),",",fp_time(i,9),",",fp_time(i,10)
end do

end do! ==========finished all targets
close(5)
call cpu_time(finish)




end program twostate_ssa