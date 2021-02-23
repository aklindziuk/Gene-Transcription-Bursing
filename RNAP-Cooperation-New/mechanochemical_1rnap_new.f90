program twostate_ssa

implicit none

integer :: n, m, i, j, dim, values(1:8)
real :: time, t_max, a, b, kn, kf, y,alpha_0, alpha_1, alpha_2
real :: dist_sum, norm_sum, tau,  mean, m2, var, p_on, p_off, pf_check
real :: t
real, dimension(1:2, 0:10000) :: dist
real, dimension(1:2, 0:10000) :: dist_norm
real, dimension(0:10000) :: histogram
real :: start, finish

character(len=10) :: file_id
character(len=50) :: file_name
character(20) filename
integer, dimension(:), allocatable :: seed
real(8) :: r1, r2


call cpu_time(start)
call date_and_time(values=values) 
call random_seed(size=i)
allocate(seed(1:i))
seed(:) = values(8)
call random_seed(put=seed)

call cpu_time(start)

! ---------------select parameters----------------

a = 15
b = 1
kn = 1
kf = 1
y = 1.0

!--------------------------------------------------
t_max = 100000
dim = 10000
call random_number (r1)
call random_number (r2)
histogram = 0
dist_norm = 0
dist = 0
    n = 0
    m = 1
    j=0
    time = 0.0
    !open(2, file = "timeseries.dat", status = "replace")
    do while (time < t_max) 
        
        if ( m == 1) then
            if ( n == 0) then
                !!! case 1 ----------> ON (n=0)
                alpha_0 = kf + a 
                tau = (1/alpha_0)*log(1/r1)
                if (time > t_max*.20) then
                    dist(m,n) = dist(m,n) + tau
                    time = time + tau
                    !write(2,*) time, 'case 1', m, n, j, tau
                end if
                time = time + tau
                
                if ( r2 < a/alpha_0) then
                    n = n + 1
                else
                    m = m + 1
                end if
            else
                !!! case 2 ------------> ON (n>0)
                alpha_0 = kf + a + n*b
                tau = (1/alpha_0)*log(1/r1)
                if (time > t_max*.20) then
                    dist(m,n) = dist(m,n) + tau
                    time = time + tau
                    !write(2,*) time, 'case 2', m, n, j, tau, r1
                end if
                time = time + tau
                alpha_1 = kf/alpha_0
                alpha_2 = (kf + a)/alpha_0
                if (r2 < alpha_1) then !(if kf)
                    m = m + 1 
                else if (alpha_1 < r2 .and. r2 < alpha_2) then !(if a)
                    n = n + 1
                else !(if nb)
                    n = n - 1 
                end if
            end if
        else 
            if ( n == 0) then
                !!! case 3 ----------> OFF (n=0)
                alpha_0 = kn + a/(y**(j+1))
                tau = (1/alpha_0)*log(1/r1)
                if (time > t_max*.20) then
                    dist(m,n) = dist(m,n) + tau
                    !P(j) = P(j) + tau
                    time = time + tau
                    !write(2,*) time,'case 3', m, n, j, tau, r1
                end if
                time = time + tau
                alpha_1 = kn/alpha_0
                if (r2 < alpha_1) then !(if kn)
                    m = m - 1
                    j = 0 
                else 
                    n = n + 1
                    j = j + 1
                end if
            else
                !!! case 4 ---------------> OFF (n>0)
                alpha_0 = kn + n*b + a/(y**(j+1))
                tau = (1/alpha_0)*log(1/r1)
                if (time > t_max*.20) then
                    dist(m,n) = dist(m,n) + tau
                    !P(j) = P(j) + tau
                    time = time + tau
                    !write(2,*) time,'case 4', m, n, j, tau, r1
                end if
                time = time + tau
                alpha_1 = kn/alpha_0
                alpha_2 = (kn + n*b)/alpha_0
                if ( r2 < alpha_1) then !(if kn)
                    m = m - 1
                    j = 0
                else if (r2>alpha_1 .and. r2<alpha_2) then !(if nb)
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
!print*, "dist sum", dist_sum
!------------normalize distribution, check sum = 1 -------------------
do j = 1, 2
    do i = 0, dim
        dist_norm(j,i) = dist(j,i)/dist_sum
        norm_sum = norm_sum + dist_norm(j,i)
    end do
end do
!print*, "norm sum", norm_sum
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
print*, 'fano', var/mean

end program twostate_ssa