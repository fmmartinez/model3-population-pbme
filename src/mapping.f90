module m_map
implicit none

real(8),parameter :: pi=3.1415926535d0

contains

subroutine get_force_traceless(nmap,lld,kosc,x,c2,rm,pm,f)
implicit none

real(8),dimension(:),allocatable :: c
real(8),dimension(:),intent(in) :: rm,pm,x
real(8),dimension(:),intent(out) :: f

integer :: a,j,n
integer,intent(in) :: nmap

real(8) :: trace,tn
real(8),dimension(:),intent(in) :: kosc,c2
real(8),dimension(:,:),intent(in) :: lld
real(8),dimension(:,:),allocatable :: mdh

allocate(mdh(1:nmap,1:nmap))
allocate(c(1:nmap))

n = size(c2)

f = 0d0
!getting product for faster calculation
do a = 1,nmap
   c(a) = 0.5d0*(rm(a)**2d0 + pm(a)**2d0)
end do

do j = 1, n
   f(j) = -kosc(j)*x(j)
   
   mdh = (lld)*(-c2(j)*2d0)
   
   trace = 0d0
   do a = 1, nmap
      trace = trace + mdh(a,a)
   end do

   tn = trace/nmap
   !for force trace is substracted, in hamiltonian the trace is added (F = -DivV)
   f(j) = f(j) +  tn

   do a = 1, nmap
      f(j) = f(j) + (mdh(a,a) - tn)*c(a)
   end do
end do

end subroutine get_force_traceless

subroutine make_hm_traceless(nmap,hm,tn)
implicit none

real(8),dimension(:,:),intent(inout) :: hm

integer :: i
integer,intent(in) :: nmap

real(8) :: trace
real(8),intent(out) :: tn

trace = 0d0
do i = 1, nmap
   trace = trace + hm(i,i)
end do

trace = trace/nmap

do i = 1, nmap
   hm(i,i) = hm(i,i) - trace
end do

tn = trace

end subroutine make_hm_traceless

subroutine get_totalenergy(nmap,hm,pm,rm,x,p,kosc,etotal)
implicit none

integer :: i,j,n
integer,intent(in) :: nmap

real(8),intent(out) :: etotal
real(8),dimension(:),intent(in) :: x,p,kosc,rm,pm
real(8),dimension(:,:),intent(in) :: hm

n = size(x)

etotal = 0d0
!classical part
do i = 1, n
   etotal = etotal + 0.5d0*(p(i)**2 + kosc(i)*x(i)**2)
end do
!mapping part
do i = 1, nmap
   do j = 1, nmap
      if (i == j) then
         etotal = etotal + 0.5d0*(hm(i,j)*(pm(i)*pm(j) + rm(i)*rm(j) - 1d0))
      else
         etotal = etotal + 0.5d0*(hm(i,j)*(pm(i)*pm(j) + rm(i)*rm(j)))
      end if
   end do
end do

end subroutine get_totalenergy

subroutine get_totalenergy_traceless(nmap,hm,tn,pm,rm,x,p,kosc,etotal,es,eb)
implicit none

integer :: i,j,n
integer,intent(in) :: nmap

real(8),intent(in) :: tn
real(8),intent(out) :: etotal,es,eb
real(8),dimension(:),intent(in) :: x,p,kosc,rm,pm
real(8),dimension(:,:),intent(in) :: hm

n = size(x)

etotal = 0d0
!classical part
do i = 1, n
   etotal = etotal + 0.5d0*(p(i)**2 + kosc(i)*x(i)**2)
end do
eb = etotal
!mapping part
etotal = etotal + tn
do i = 1, nmap
   do j = 1, nmap
      etotal = etotal + 0.5d0*(hm(i,j)*(pm(i)*pm(j) + rm(i)*rm(j)))
   end do
end do
es = etotal - eb - tn
end subroutine get_totalenergy_traceless

subroutine get_coeff(ng,rm,pm,coeff)
implicit none


integer :: i
integer,intent(in) :: ng

real(8) :: z
real(8),intent(out) :: coeff
real(8),dimension(:),intent(in) :: rm,pm
real(8),dimension(:),allocatable :: exp_be,prob

allocate(exp_be(1:ng),prob(1:ng))

exp_be = 0d0
z = 0d0
do i = 1, ng
   !exp_be(i) = exp(-beta*omega*(i - 0.5d0))
   !exp_be(i) = exp(-1.44d0*(i - 0.5d0))
   !the equation beloww corresponds to the one on commented first line on top
   exp_be(i) = exp(-2.29d-1*(i - 0.5d0))
   !z = z + exp_be(i)
end do

prob = exp_be!/z

coeff = 0d0
!only diagonal because in the ngxng matrix by construction their lambdas are
! (1,0,0...), (0,1,0...), (0,0,1...), etc
do i = 1, ng
   coeff = coeff + (rm(i)**2 + pm(i)**2 - 0.5d0)*prob(i)
end do

deallocate(exp_be)
deallocate(prob)
end subroutine get_coeff

subroutine evolve_rm(nmap,dt,hm,pm,rm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt
real(8),dimension(:),intent(in) :: pm
real(8),dimension(:),intent(inout) :: rm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      rm(i) = rm(i) + dt*hm(i,j)*pm(j)
   end do
end do

end subroutine evolve_rm

subroutine evolve_pm(nmap,dt2,hm,rm,pm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt2
real(8),dimension(:),intent(in) :: rm
real(8),dimension(:),intent(inout) :: pm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      pm(i) = pm(i) - dt2*hm(i,j)*rm(j)
   end do
end do

end subroutine evolve_pm

subroutine get_facts_pop_traceless(nmap,coeff,llg,llb,lld,rm,pm,fact1,fact2,fact3)
implicit none

integer :: a,b,i
integer,intent(in) :: nmap

real(8) :: trt
real(8),intent(in) :: coeff
real(8),intent(out) :: fact1,fact2,fact3
real(8),dimension(:),intent(in) :: rm,pm
real(8),dimension(:,:),intent(in) :: llg,llb,lld

trt = 0d0
do i = 1, nmap
   trt = trt + llg(i,i)
end do
fact1 = trt/nmap
do a = 1, nmap
   do b = 1, nmap
      if (a == b) then
         fact1 = fact1 + 0.5d0*(llg(a,b) - trt/nmap)*(rm(a)*rm(b) + pm(a)*pm(b))
      else
         fact1 = fact1 + 0.5d0*(llg(a,b))*(rm(a)*rm(b) + pm(a)*pm(b))
      end if
   end do
end do
fact1 = coeff*fact1

trt = 0d0
do i = 1, nmap
   trt = trt + llb(i,i)
end do
fact2 = trt/nmap
do a = 1, nmap
   do b = 1, nmap
      if (a == b) then
         fact2 = fact2 + 0.5d0*(llb(a,b) - trt/nmap)*(rm(a)*rm(b) + pm(a)*pm(b))
      else
         fact2 = fact2 + 0.5d0*llb(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
      end if
   end do
end do
fact2 = coeff*fact2

trt = 0d0
do i = 1, nmap
   trt = trt + lld(i,i)
end do
fact3 = trt/nmap
do a = 1, nmap
   do b = 1, nmap
      if (a == b) then
         fact3 = fact3 + 0.5d0*(lld(a,b) - trt/nmap)*(rm(a)*rm(b) + pm(a)*pm(b))
      else
         fact3 = fact3 + 0.5d0*lld(a,b)*(rm(a)*rm(b) + pm(a)*pm(b))
      end if
   end do
end do
fact3 = coeff*fact3
end subroutine get_facts_pop_traceless

subroutine get_facts_pop(nmap,ng,nb,coeff,rm,pm,fact1,fact2,fact3)
implicit none

integer :: a
integer,intent(in) :: nmap,ng,nb

real(8),intent(in) :: coeff
real(8),intent(out) :: fact1,fact2,fact3
real(8),dimension(:),intent(in) :: rm,pm

fact1 = 0d0
do a = 1, ng
   fact1 = fact1 + (rm(a)**2 + pm(a)**2 - 1d0)
end do
fact1 = coeff*fact1

fact2 = 0d0
do a = ng+1, ng+nb
   fact2 = fact2 + (rm(a)**2 + pm(a)**2 - 1d0)
end do
fact2 = coeff*fact2

fact3 = 0d0
do a = ng+nb+1, nmap
   fact3 = fact3 + (rm(a)**2 + pm(a)**2 - 1d0)
end do
fact3 = coeff*fact3
end subroutine get_facts_pop

function kronecker_delta(i,j) result (d)
implicit none

integer,intent(in) :: i, j

real(8) :: d

if (i == j) then
   d = 1d0
else
   d = 0d0
end if

end function kronecker_delta

subroutine get_hm2(nmap,ng,nb,mu,et,a1,a2,hs,hm)
implicit none

integer :: i
integer,intent(in) :: nmap,ng,nb

real(8),intent(in) :: mu,et,a1,a2
real(8),dimension(:,:),intent(in) :: hs
real(8),dimension(:,:),intent(out) :: hm

hm = hs
hm(1:ng,ng+1:ng+nb) = hs(1:ng,ng+1:ng+nb)*(-mu*et)
hm(ng+1:ng+nb,1:ng) = hs(ng+1:ng+nb,1:ng)*(-mu*et)
do i = ng+nb+1,nmap
   hm(i,i) = hm(i,i) + a1 + a2
end do
end subroutine get_hm2

subroutine get_preh(ng,nb,nd,eg,eb,ed,delta,omega,hs)
use m_vib
implicit none

integer,parameter :: ip = 10000

character(len=2) :: c_nt
character(len=9) :: fmt1

integer :: i,j,nt
integer,intent(in) :: ng,nb,nd

real(8) :: cg,cb,cd,uint,lint,alpha
real(8),intent(in) :: eg,eb,ed,delta,omega
real(8),dimension(:,:),intent(out) :: hs

nt = ng + nb + nd

cg = 0d0
cb = 2d0*sqrt(10d0)/omega
cd = cb/2d0

uint = cb + 6d0
lint = cb - 6d0

alpha = sqrt(omega)

hs = 0d0

!fill g|g
do i = 1, ng
   hs(i,i) = eg + (i - 0.5d0)*omega
end do
!fill g|b
do i = 1, ng
   do j = ng+1, ng+nb
      hs(i,j) = integrate_t_phiphi(ip,lint,uint,i,cg,j-ng,cb,alpha)
   end do
end do
!fill g|d
!not necessary
!fill b|g
do i = ng+1, ng+nb
   do j = 1, ng
      hs(i,j) = integrate_t_phiphi(ip,lint,uint,i-ng,cb,j,cg,alpha)
   end do
end do
!fill b|b
do i = ng+1, ng+nb
   hs(i,i) = eb + (i -ng - 0.5d0)*omega
end do
!fill b|d
do i = ng+1, ng+nb
   do j = ng+nb+1, nt
      hs(i,j) = delta*integrate_t_phiphi(ip,lint,uint,i-ng,cb,j-ng-nb,cd,alpha)
   end do
end do
!fill d|g
!not necessary
!fill d|b
do i = ng+nb+1, nt
   do j = ng+1, ng+nb
      hs(i,j) = delta*integrate_t_phiphi(ip,lint,uint,i-ng-nb,cd,j-ng,cb,alpha)
   end do
end do
!fill d|d
do i = ng+nb+1, nt
   hs(i,i) = ed + (i -ng -nb - 0.5d0)*omega
end do

if (nt > 9) then
   write(c_nt,'(i2)') nt
else
   write(c_nt,'(i1)') nt
end if

fmt1 = '('//trim(c_nt)//'f10.5)'

!print fmt1, hs

end subroutine get_preh


subroutine iniconq_d(nosc,lumda_d,ome_max,ome,c2,kosc)
implicit none

integer :: i
integer,intent(in) :: nosc

real(8) :: check
real(8),intent(in) :: ome_max,lumda_d
real(8),dimension(:),intent(inout) :: ome,c2,kosc

check=0
do i=1,nosc
!  ome(i)=tan(pi*i/(4*nosc))
   ome(i)=tan(i*datan(ome_max)/nosc)
!  ome(i)=i**2*ome_max/nosc**2
!  rho(i)=sqrt(ome(i)*ome_max)/(1+ome(i)**2)
!  c2(i)=0.5d0*ome(i)*dsqrt(lumda_d/nosc)
   c2(i)=ome(i)*dsqrt(datan(ome_max)*lumda_d/(pi*nosc))
!  c2(i)=ome(i)*sqrt(lumda_d*rho(i)*2/(pi*nosc))
   check=check+c2(i)**2/ome(i)**2
   kosc(i)=ome(i)**2 
end do
write(6,*) check

end subroutine iniconq_d

end module m_map
