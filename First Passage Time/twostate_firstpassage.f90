!------------------------------------------------------------------
! Building on the two-state Gillespie simulation for transcription, this 
! code determines the first passage time to reach a given target. The first
! passage time (fp) is the time it takes to reach a target number of mRNAs, 
! designated by n, for the first time. a, b, kn, kf are the kinetic parameters
! of the system. target  is the target number of mRNAs. trails is the number
! of times the system is launeched (aka how many times fp time is collected)
! The output files of interest are fptime.dat, which is a record of the fp
! time recorded during each trial, and the fphist.dat, which gives a histogram 
! of the recorded fp values (there might be a small error (slight shift) in the 
! histogram. The command line output gives mean and variance of the the fp time - 
! these calculations re accurate. 
!-------------------------------------------------------------------


program twostate_ssa

implicit none

integer :: n, m, i, j, k, dim, trials, counter, target
real :: time, t_max, random1, random2, a, b, kn, kf, alpha_0, alpha_1, alpha_2, mean, mean2, var, var2, sum, tau, fon, dnt
real :: nn, ff, on
real, dimension(1:2, 0:250) :: dist
real, dimension(1:2, 0:250) :: dist_norm
real, dimension(0:250) :: histogram
real :: start, finish
real, dimension(1:500) :: fp_time
real, dimension(1:100) :: fp_hist
real :: fpmean, min, max, fpsum, bin, fpvar

call cpu_time(start)

!----------set parameters-----------------
a = 100
b = 1

t_max = 10000
target = 130

kn = 3
kf = 1

dim = 250 ! set to the size of dist (max poss n of system)
trials = 500 ! how many times fp time is sampled
!-----------------------------------------------

open(2, file = "timeseries.dat", status = "replace")

call random_seed

fp_time = 0
fp_hist = 0

do i = 1, size(fp_time) !----------start trials-------
	
	dist_norm = 0
	dist = 0

	n = 0
	m = 1
	time = 0.0
	counter = 1

	do while (time < t_max) !-----------start a random walk-------
		
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

				alpha_0 = kf + a + n*b
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				alpha_1 = kf/alpha_0
				alpha_2 = (kf + a)/alpha_0
				if (random2 < alpha_1) then
					m = m + 1 
				else if (alpha_1 < random2 .and. random2 < alpha_2) then
					n = n + 1
				else 
					n = n - 1 
				end if
			else

				!!! case 4 --------------->

				alpha_0 = kn + n*b
				tau = (1/alpha_0)*log(1/random1)
				if (time > t_max*.20) then
					dist(m,n) = dist(m,n) + tau
				end if
				time = time + tau

				alpha_1 = kn/alpha_0
				if ( random2 < alpha_1) then
					m = m - 1
				else
					n = n - 1
				end if 

			end if
		end if !-------- step taken -----------

		if ( n == target ) then
			if (counter == 1) then
				fp_time(i) = time !if target is reached for the first time, record
			end if
			counter = counter + 1 !counts how many time target is reached
		end if

	if (i == 1) then
		write(2,*) time, m, n !record substate at t_max (during 1st trial)
	end if 

	end do !----------one trial over---------------
	close(2)

	if (counter == 1) then
		print*, "target not reached in ", t_max, "time units"
		fp_time(i) = t_max
	end if

end do !---------all trials over-----------------

open (5, file = "fptime.dat", status = "replace")
do i = 1, size(fp_time)
	write(5,*), i, fp_time(i)
end do
close(5)

!------------------fp mean------------------------

fpmean = 0
do i = 1,size(fp_time)
	fpmean = fpmean + fp_time(i)
end do
fpmean = fpmean/size(fp_time)
print*, 'fp mean', fpmean

!--------------calculate histogram--------------------

min = minval(fp_time)
max = maxval(fp_time)
!print*, "min", min
!print*, "max", max
bin = (max - min)/real(size(fp_hist))
!print*, 'bin size', bin

do i = 1, size(fp_time)
	do j = 1, size(fp_hist)
		if ((fp_time(i).ge.(min+bin*(j-1))).and.(fp_time(i).le.(min+bin*j))) then
			fp_hist(j) = fp_hist(j)+1
		end if
	end do
end do

fpsum = 0
do i = 1,size(fp_hist)
	fpsum = fpsum + fp_hist(i)
end do
!print*, "histogram sum", fpsum

open(7, file = "fphist.dat", status = "replace")
do i = 1,size(fp_hist)
		fp_hist(i) = fp_hist(i)/fpsum
		write(7, *), min + (i + .5)*bin, fp_hist(i)
end do
close(7)

fpsum = 0
do i = 1,size(fp_hist)
	fpsum = fpsum + fp_hist(i)
end do
!print*, "histogram sum 2", fpsum

!----------------calculate variance ----------------------------------

fpmean = 0
do i = 1,size(fp_hist)
		fpmean = fpmean + (min + (i + .5)*bin)*fp_hist(i)
end do
! print*, 'fp mean', fpmean

fpvar = 0
do i = 1, size(fp_hist)
	fpvar = fpvar + fp_hist(i)*(fpmean - (min + (i + .5)*bin))**2
end do
print*, 'fp var', fpvar

print*, ' '


!-----------------------------------------------------

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
 	write(3, *) j, i, dist(j,i)
 end do
end do
close(3)

! normalize distribution, find mean
open(1, file = "dist_norm.dat", status = "replace")
do j = 1,2
 do i = 0, dim
 	dist_norm(j,i) = real(dist(j,i))/sum
 	write(1, *) i, dist_norm(j,i)
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

open(4, file = "histogram.dat", status = "replace")
do i = 0,dim
		write(4, *) i, histogram(i)
end do
close(4)

ff = 1 + ((1 - kn/(kn + kf))*a)/(kn + kf + b)
nn = (kn* a)/((kn + kf)* b)
on = kn/(kn + kf)

call cpu_time(finish)

print*, "results:   ","simulaton       ", 'analytical      ','error'
print*, '-------------------'
print*, "mean =  ", mean, nn, (mean-nn)
print*, "fano2 =", var2/mean, ff, (var2/mean - ff)
print*, "fon =  ", fon, on, (fon- on)
print*, '-------------------'
print*, "variance2", var2
print*, 'total time', sum
print*, 'check norm time', dnt
!print *,"simul time =", finish-start



end program twostate_ssa