program ham
implicit none
integer :: p,q,L,Ny,Nz,Nmax,Nk,i,j,pp,qp,Nkx,i1,j1,offX,offY,m,kpSize
real*8 :: Ly,Lz,pi,Kmax,gam1,gam2,gam3,gamc,meff,Dir,Eph,ScfPerc,Ep0,P0B
real*8 ::  Ep,Eg,angs2bohr,nano2bohr,ev2Hartree,dSO,P0,lb,hb,gam1L,gam2L,gam3L,Eg0
integer,allocatable :: KnY(:),KnZ(:),IWORK(:),Indices(:)
real*8, allocatable :: Ky(:),Kz(:),W(:),RWORK(:),Ener(:),Umat(:,:)
real*8 :: kx,kyy,kyyp,kzz,kzzp,FERMI_LEVEL,CUTOFF,CUTOFF1,k_span,Lx,dx,dE,dy,dz,Efl,Efr,dU
complex*16, allocatable :: Hkp_SO(:,:),WORK(:),SigL(:,:,:),SigR(:,:,:),gL(:,:,:),gR(:,:,:)
integer ::  LRWORK,LWORK,LIWORK,INFO,matSize,cnt,Nx,Nen,Nm,buffReg,flagL,flagR,NphCov
character(len=30) :: fileName
integer :: SpEmTerm,AbsTerm,StEmTerm,NphCovInit
character(len=4) :: Sep
character(len=1) :: Skp
character(len=3) :: SNk

parameter(angs2bohr=1.889)
parameter(nano2bohr=18.89)
parameter(ev2Hartree=0.036749309)
parameter(dSO=0.044*ev2Hartree)
parameter(Nmax=10)
parameter(Ly=2.4*nano2bohr)
parameter(Lz=2.4*nano2bohr)
parameter(gam1L=4.285)
parameter(gam2L=0.339)
parameter(gam3L=1.446)
parameter(meff=0.26)
parameter(Eg0=1.42*ev2Hartree)
parameter(Ep=24.3*ev2Hartree)
parameter(Ep0=24.3*ev2Hartree)
parameter(kpSize=2)
parameter(Nx=140)



pi=4.*atan(1.0)
Ny=1
Nz=1
!gamC=(1./(2.*meff)-Ep/3.*(1./Eg0-1./(Eg0+dSO)))
gamC=0.D0
Eg=1.42*ev2Hartree

gam1=gam1L-Ep/(3.*Eg0)
gam2=gam2L-Ep/(6.*Eg0)
gam3=gam3L-Ep/(6.*Eg0)



!write(*,*) gamC
!pause

P0=dsqrt(Ep/2.)
P0B=dsqrt(Ep0/2.)

StEmTerm=1
SpEmTerm=0
absTerm=1


Eph=1.60*ev2Hartree

dx=0.3*nano2bohr


NphCovInit=100

!NphcovInit=40

ScfPerc=0.1




dy=Ly/Ny
dz=Lz/Nz


!dE=Eph/NphCov




Lx=Nx*dx

!  CUTOFF=1.2*ev2Hartree
!  CUTOFF1=0.7*ev2Hartree


!  CUTOFF=.9*ev2Hartree
!  CUTOFF1=.3*ev2Hartree




dU=1.4*ev2Hartree

if(Eph.ge.dU) then
CUTOFF=.0*ev2Hartree+Eph+0.25*ev2Hartree
CUTOFF1=Eph-(1.42*ev2Hartree-dU)+0.25*ev2Hartree
endif



if(Eph.lt.dU) then
CUTOFF=.0*ev2Hartree+Eph+0.4*ev2Hartree
CUTOFF1=Eph-(1.42*ev2Hartree-dU)+0.4*ev2Hartree
endif




!dU=0.
Dir=1.
BuffReg=4

do i=NphCovInit,NphCovInit+10
NphCov=i
dE=dble(Eph/i)
Nen=int((cutoff1+cutoff)/dE)
if(mod(Nen,2).eq.0) exit
enddo

lb=-Cutoff1
hb=Cutoff

100     format(10000(2x,F18.12))






call getKval(Nmax,Ny,Nz,pi,Nk)




allocate(SigL(Nen,kpSize*Nk,kpSize*Nk),SigR(Nen,kpSize*Nk,kpSize*Nk),gL(Nen,kpSize*Nk,kpSize*Nk),gR(Nen,kpSize*Nk,kpSize*Nk))
allocate(Ener(Nen))
allocate(KnY(Nk),KnZ(Nk),Ky(Nk),Kz(Nk))
allocate(Umat(kpSize*Ny*Nz,kpSize*Nk))

!dE=(cutOff1+cutOFF)/(Nen-1)
do i=1,Nen
Ener(i)=-CutOFF1+dE*(i-1)

write(*,*) Nen,(-CutOFF1+dE*(i-1))/ev2Hartree
enddo
!pause


call fillKval(Nk,Nmax,Ny,Nz,pi,KnY,KnZ,Ky,Kz,Ly,Lz)




call  GetUmat(Nk,KnY,KnZ,Ky,Kz,Ly,Lz,Umat,Ny,Nz,dy,dz,kpSize)



call getFullMatrix(Nk,KnY,KnZ,Ky,Kz,gam1,gam2,gam3,gamC,P0,P0B,Eg,Nx,dx,kpSize,SigL,SigR,Ener,Nen,CUTOFF,CUTOFF1,Umat,Ny,Nz,dU,Dir,BuffReg,dy,dz,NphCov,Eph,lb,hb,absTerm,SpEmTerm,ScfPerc)
deallocate(Umat)
end program ham


subroutine getKval(Nmax,Ny,Nz,pi,Nk)
implicit none
integer p,q,Ny,Nz,Nmax,cnt,modK,Nk
real*8 pi
cnt=0
do p=1,Ny
do q=1,Nz
cnt=cnt+1
enddo
enddo

Nk=cnt
end subroutine getKval




subroutine fillKval(Nk,Nmax,Ny,Nz,pi,KnY,KnZ,Ky,Kz,Ly,Lz)
implicit none
integer p,q,Ny,Nz,Nmax,cnt,modK,Nk,KnY(Nk),KnZ(Nk)
real*8 pi,Ky(Nk),Kz(Nk),Ly,Lz
cnt=0
do p=1,Ny
do q=1,Nz
cnt=cnt+1
KnY(cnt)=(p-1)
KnZ(cnt)=(q-1)
Ky(cnt)=(p-1)*pi/Ly
Kz(cnt)=(q-1)*pi/Lz
enddo
enddo

Nk=cnt
end subroutine fillKval



subroutine getMatrix(gam1,gam2,gam3,gamC,P0,Eg,Hkp4,Wmat,dx,kpSize)

implicit none
integer p,q,pp,qp,kpSize
real*8 Kx,Ky,Kz,KyP,KzP,gam1,gam2,gamc,dl,dqqp,dppp,dppluspp,dqplusqp,P0,gam3,dx
COMPLEX*16 C,Pm,Qm,R,S,Pplus,Pminus,Pz,im,Hkp4(kpSize,kpSize),Wmat(kpSize,kpSize),WmatL(kpSize,kpSize)
real*8 sqr3,Kx2,Ky2,Kz2,prodKz,prodKy,pi,L,M,N,Eg,sign1,sign2


im=(0.d0,1.d0)
pi=4.*Atan(1.d0)
Hkp4=0.
Wmat=0.

!L=0.5*(gam1+4.*gam2)

L=13.04/2.


gamC=8./2.
L=-1./2.



write(*,*) L

!pause
Hkp4(1,1)=gamC*2./(dx*dx)+Eg
Hkp4(2,2)=L*2./(dx*dx)
Hkp4(1,2)=0.
Hkp4(2,1)=0.

Wmat(1,1)=gamC*(-1./(dx*dx))
Wmat(2,2)=(L*-1./(dx*dx))
Wmat(1,2)=P0*1./(2.*dx)
Wmat(2,1)=-P0*1./(2.*dx)



end subroutine getMatrix



subroutine getFullMatrix(Nk,KnY,KnZ,Ky,Kz,gam1,gam2,gam3,gamC,P0,P0B,Eg,Nx,dx,kpSize,SigL,SigR,Ener,Nen,CUTOFF,CUTOFF1,Umat,Ny,Nz,dU,Dir,BuffReg,dy,dz,Nphcov,Eph,lb,hb,absTerm,SpEmTerm,ScfPerc)
implicit none

integer Nk,KnY(Nk),Knz(Nk),i,j,i1,j1,p,q,pp,qp,offX,offY,kpSize,offA,OffSet,Nx,cnt,Nen,NitMax,Nm,Ny,Nz,buff,flag,flag1,flagIt
integer SpEmTerm,AbsTerm,StEmTerm
real*8 ScfPerc,P0B,dEfr
real*8 kyy,kzz,kyyp,kzzp,Ky(Nk),Kz(Nk),gam1,gam2,gam3,gamc,P0,Eg,Kx,dx,Ener(Nen),eps,ev2Hartree,Efl,Efr
real*8 temp1v,temp2v,temp3v,temp4v,temp5v,temp6v,beta,Umat(kpSize*Ny*Nz,kpSize*Nk),V,dU,T(Nen),summa,vol,dy,dz,lb,hb
complex*16 Hkp_SO(kpSize*Nk,kpSize*Nk),Wmat(kpSize,kpSize),Hkp(kpSize,kpSize),Ham(kpSize*Nk,kpSize*Nk)
real*8 dir,EnerProf(kpSize*Nk),NtotNeg,NtotPos,summ2,Enev,currAlt,Nph
complex*16,allocatable:: TEMP1N(:,:),U1(:,:),UNx(:,:)
integer buffreg,cntCB,jj,matSize,TotNegM,TotPosM,Niter,Niterinit,iscf,cntAll
complex*16 im
complex*16,allocatable::gDiag(:,:,:),gUp(:,:,:),gLow(:,:,:),gRTEMP(:,:),gL(:,:,:),gR(:,:,:)
complex*16,allocatable :: gLTEMP(:,:),TEMP(:,:),TEMP1(:,:),TEMP2(:,:),CMAT(:,:),TEMP3(:,:),TEMP4(:,:),TEMP5(:,:,:),TEMP6(:,:,:)
complex*16,allocatable :: Wm(:,:),WmDiager(:,:),SigL(:,:,:),SigR(:,:,:),WmM(:,:,:),WmMD(:,:,:),EigVctr0(:,:),EigVctr1(:,:)
complex*16,allocatable :: Fmat(:,:),Gmat(:,:),H1(:,:),HNx(:,:),SigInL(:,:,:),SigInR(:,:,:)
complex*16,allocatable :: gNdiag(:,:,:),gNlow(:,:,:),gNup(:,:,:),gPdiag(:,:,:),gPlow(:,:,:),gPUp(:,:,:)
real*8, allocatable :: Vmat(:,:,:),VmatR(:,:,:)
complex*16 imPart,summCmplx,summCmplxHole

integer :: LWORKS,INFO,SchD,sdim,Nmix,Nmod,NmodesNeg,NmodesPos,NmodesNeg0,NmodesPos0,diff1,diff2
integer,allocatable:: Indices(:)
parameter(ev2Hartree=0.036749309)
real*8 error,summ,alpha,alpha_init,summ1,EigEn(kpSize*Nk),CUTOFF,CUTOFF1,Ebottom,Eup,kt
real*8 fermiL,fermiR,summEl,summCur,summEl1,currDens(Nx-1),eps_medium,Eph,Avec(3),currDensHole(Nx-1),currElIn(Nx-1),currHIn(Nx-1),ElDensityLayer(Nx-1)
real*8 currDensCB(Nx-1),currDensVB(Nx-1),currElCB(Nx-1),currElVB(Nx-1),CurPrev,CurCur
real*8 currDensBall(Nx-1,Nen),currDensHoleBall(Nx-1,Nen),currDensPhoto(Nx-1),dE,pi
real*8, allocatable :: TempEn(:),RWORKS(:),Theta(:)
complex*16 ModeHam0(kpSize*Nk,kpSize*Nk),EphBase(kpSize,kpSize)
complex*16,allocatable :: Hm(:,:,:),EnerDiag(:,:),Ediag(:,:),EigVctrMode(:,:),gPREVL(:,:),gPREVR(:,:),EigVctrMode1(:,:),HamEph(:,:)
complex*16,allocatable :: HamEphMode(:,:,:),FullSigOut(:,:)
complex*16,allocatable :: rhoN(:,:,:),Fn(:,:,:),diffRhoFn(:,:,:),diffFnFnm1(:,:,:),Unm(:,:),Vn(:),An(:),Hreal(:,:,:),HrealBall(:,:,:)
complex*16,allocatable :: Amat11(:,:),GamL(:,:),GamR(:,:),U1H(:,:),UNxH(:,:),FullSigIn(:,:),Htemp(:,:),Gn(:,:,:),SigInEPh(:,:,:,:),SigREph(:,:,:,:),SigOutEph(:,:,:,:)
complex*16,allocatable :: Gp(:,:,:),FullSigInEph(:,:,:),FullSigOutEph(:,:,:),FullSigREph(:,:,:),HamEphFull(:,:),GnPrev(:,:,:),GpPrev(:,:,:)
integer :: Nscf,NphCov,l,CntEntries
complex*16,allocatable :: gNdiagEn(:,:,:,:),gPDiagEn(:,:,:,:),gDiagEn(:,:,:,:),TEMPOR(:,:),HamEphBaseT(:,:),FullSigInEphPrev(:,:,:)
real*8 pref,absT,emsT,Pref1,Intens,beta3,summEl2,Elow,Ehigh,DivJ(Nx),HamEphNotZero,errorCur

im=(0.d0,1.D0)

write(*,*) Ny,Nz


allocate(Wm(kpSize*Nk,kpSize*Nk),WmDiager(kpSize*Nk,kpSize*Nk))

allocate(EigVctr0(kpSize*Nk,kpSize*Nk),EigVctr1(kpSize*Nk,kpSize*Nk))

allocate(Vmat(Nx,kpSize*Ny*Nz,kpSize*Ny*Nz),VmatR(Nx,kpSize*Nk,kpSize*Nk),TEMP1N(kpSize*Nk,kpSize*Nk))
allocate(H1(kpSize*Nk,kpSize*Nk),HNx(kpSize*Nk,kpSize*Nk))


dE=dabs(Ener(2)-Ener(1))
NitMAx=1255000
Nmix=2
Nmod=4
eps=1.D-8
NtotNeg=1.
NtotPos=1.
alpha_init=0.2
imPArt=4.D-21*im
Niterinit=100
kt=0.0256*ev2Hartree
Nscf=200

vol=Nx*dx*10.*10.*18.89*18.89

eps_medium=12.
Pref=dsqrt(4.*3.14159/(2.*Eph*vol))

Avec(1)=1.*Pref
Avec(2)=0.
Avec(3)=0.


pi=4.*atan(1.0)

!ScfPerc=8.


!################ Intensity in Mw/m^2 #################
Intens=0.5

Pref1=2.1D+7


Nph=Pref1*Intens*vol/(18.89**3.)*dsqrt(eps_medium)/(Eph/0.036749309)


!Nph=.0000794

!ScfPerc=1.
!Nph=0.008

Nph=0.0000001918 !################# Nanowire 20nmX20nm, Iw=10^-5(MW/cm^2)

write(*,*) Nph,Eph,kt,Nk,kpSize
allocate(Unm(Nmix,Nmix),Vn(Nmix),An(Nmix))



open(unit=22, file="dos_Pot_mode_nm.dat", status="unknown")
open(unit=23, file="transmission_mode_nm.dat", status="unknown")
open(unit=24, file="EnProf_mode_nm.dat", status="unknown")
open(unit=25, file="currentEn_mode_nm.dat", status="unknown")
open(unit=26, file="electronDensity_mode_nm.dat", status="unknown")
open(unit=27, file="elDensity_mode2D_nm.dat", status="unknown")
open(unit=28, file="currDensity_nm.dat", status="unknown")
open(unit=29, file="currDensityHole_nm.dat", status="unknown")
open(unit=37, file="PhotoCurrDensity.dat", status="unknown")
open(unit=31, file="curentIntegrated_nm_El.dat", status="unknown")
open(unit=38, file="curentIntegrated_nm_Hole.dat", status="unknown")
Open(unit=52, file="PhotoCurrent2D.dat", status="unknown")





call getMatrix(gam1,gam2,gam3,gamC,P0B,Eg,Hkp,Wmat,dx,kpSize)
Hkp_SO=Hkp
Wm=Wmat
WmDiager=dconjg(transpose(Wm))
TEMP1N=Hkp_SO+Wm+WmDiager

call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)

do j=1,kpSize*Nk
if((EnerProf(j).le.0.).and.(EnerProf(j+1).ge.0.)) then
Efl=EnerProf(j)+0.055*ev2Hartree


endif
enddo

!Efr=Efl
!Efr=Efl


ModeHam0=Hkp_SO+Wm+WmDiager


allocate(gLTEMP(kpSize*Nk,kpSize*Nk),gRTEMP(kpSize*Nk,kpSize*Nk),Ediag(kpSize*Nk,kpSize*Nk))
allocate(gL(Nen,kpSize*Nk,kpSize*Nk),gR(Nen,kpSize*Nk,kpSize*Nk),SigL(Nen,kpSize*Nk,kpSize*Nk),SigR(Nen,kpSize*Nk,kpSize*Nk))
allocate(gPREVL(kpSize*Nk,kpSize*Nk),gPREVR(kpSize*Nk,kpSize*Nk),U1H(kpSize*Nk,kpSize*Nk),UNxH(kpSize*Nk,kpSize*Nk))



call GetVpot(Vmat,Nx,Ny,Nz,kpSize,Efl,Efr,buffReg,Dir,dx,dEfr)

write(*,*) Umat

!pause
do i=1,Nx
VmatR(i,:,:)=matmul(Transpose(Umat),matmul(Vmat(i,:,:),Umat))
enddo

!Efl=Efl+0.01*ev2Hartree
!   Efl=Efl+0.01*ev2Hartree

!Efl=Efl-dU-0.02*ev2Hartree
Efr=Efl+dEfr
!+(1.1-0.01)*ev2Hartree



!Efr=Efl

do i=1,Nx
TEMP1N=ModeHam0+VmatR(i,:,:)
call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)
write(24,'(1000F12.4)') (i-1)*dx,(EnerProf(j)/ev2Hartree,j=1,Nk*kpSize),Efl/ev2Hartree,Efr/ev2Hartree
enddo


Nm=kpSize*Nk



allocate(U1(kpSize*Nk,Nm),UNx(kpSize*Nk,Nm))
allocate(Amat11(Nm,Nm),GamL(Nm,Nm),GamR(Nm,Nm))
allocate(gDiag(Nx,Nm,Nm),gUp(Nx,Nm,Nm),gLow(Nx,Nm,Nm))
allocate(TEMP(Nm,Nm),TEMP1(Nm,Nm),TEMP2(Nm,Nm))
allocate(CMAT(Nm,Nm))
allocate(WmM(Nx,Nm,Nm),WmMD(Nx,Nm,Nm),Hm(Nx,Nm,Nm),EnerDiag(Nx,Nm),EigVctrMode(kpSize*Nk,Nm),EigVctrMode1(kpSize*Nk,Nm),TempEn(Nm))
allocate(TEMP3(Nm,Nm),TEMP4(Nm,Nm))
allocate(SigInL(Nen,Nm,Nm),SigInR(Nen,Nm,Nm),gNdiag(Nx,Nm,Nm),gNlow(Nx,Nm,Nm),gNup(Nx,Nm,Nm),gPdiag(Nx,Nm,Nm),gPlow(Nx,Nm,Nm),gPup(Nx,Nm,Nm))
allocate(SigInEph(Nen,Nx-2*buffReg,Nm,Nm),SigREph(Nen,Nx-2*buffReg,Nm,Nm),SigOutEph(Nen,Nx-2*buffreg,Nm,Nm))
allocate(HamEphBaseT(kpSize*Nk,kpSize*Nk))
allocate(HamEph(kpSize*Nk,kpSize*Nk),HamEphMode(Nx-2*buffReg,Nm,Nm))
allocate(gNdiagEn(Nen,Nx,Nm,Nm),gPDiagEn(Nen,Nx,Nm,Nm),gDiagEn(Nen,Nx,Nm,Nm))



H1=Hkp_SO+VmatR(1,:,:)
HNx=Hkp_SO+VmatR(Nx,:,:)




!##################################### Electron photon matrix elements ###########################3

EphBase=0.

HamEph=0.

do i=1,kpSize
do j=1,kpSize
if(Avec(1).ne.0.) then
EphBase(1,2)=im*P0*Avec(1)
EphBase(2,1)=-im*P0*Avec(1)

!         EphBase(1,2)=P0*Avec(1)
!         EphBase(2,1)=P0*Avec(1)

endif
enddo
enddo



HamEph=EphBase




do i=1,Nx
Ham=Hkp_SO+VmatR(i,:,:)
Hm(i,:,:)=Ham
WmM(i,:,:)=Wm
WmMD(i,:,:)=WmDiager
enddo



SigL=0.
SigR=0.
gL=0.
gR=0.




!############################################## Getting Contact Self-Energies ##########################

SigInEph=0.
SigOutEph=0.
SigREph=0.
flag=0
flag1=0
flagIt=0
cntAll=0


do iSCF=1,Nscf
currElIn=0.
currHIn=0.
currElCB=0.
currElVB=0.
ElDensityLayer=0.


do p=1,Nen
Niter=Niterinit
alpha=alpha_init
Ediag=0.
do i=1,kpSize*Nk
Ediag(i,i)=Ener(p)+imPart
enddo


EnEv=Ener(p)/ev2Hartree

fermiL=1./(1.+dexp((Ener(p)-Efl)/kt))
fermiR=1./(1.+dexp((Ener(p)-Efr)/kt))


!###################### IF SCF #########################
if(p.eq.1) then
gPrevL=0.
gPrevR=0.
else
gPrevL=gL(p-1,:,:)
gPrevR=gR(p-1,:,:)
endif

100 if(flag.eq.0.and.iscf.eq.1) then
gPrevL=0.
gPrevR=0.
endif

write(*,'(A7,2I5,4F20.8)') "P,ISCF",p,iSCF,CurCur,CurPrev,errorCur,ScfPerc

if(iSCF.eq.1) then
call GetConactSelfEnergies(p,H1,kpSize*Nk,gPrevL,Ediag,Wm,WmDiager,gLTEMP,Nmix,Nmod,alpha,NitMax,"gL",eps,flag,Niter)
call GetConactSelfEnergies(p,HNx,kpSize*Nk,gPrevR,Ediag,WmDiager,Wm,gRTEMP,Nmix,Nmod,alpha,NitMax,"gR",eps,flag,Niter)
gL(p,:,:)=gLTEMP
gR(p,:,:)=gRTEMP
SigL(p,:,:)=matmul(WmDiager,matmul(gL(p,:,:),Wm))
SigR(p,:,:)=matmul(Wm,matmul(gR(p,:,:),WmDiager))
write(*,*) "iSCF",isCF
write(*,*) "SigL",CDABS(SUM(TEMP1)),CDABS(SUM(SigL(p,:,:)))
write(*,*) "SigR",CDABS(SUM(TEMP2)),CDABS(SUM(SigR(p,:,:)))

endif


GamL=-2.*Dimag(SigL(p,:,:))
GamR=-2.*Dimag(SigR(p,:,:))


write(*,*) "ISCF",ISCF


!###################### IF SCF #########################


if(iScf.eq.1) then
SigInL(p,:,:)=GamL*im*fermiL
SigInR(p,:,:)=GamR*im*fermiR

endif
!###################### IF SCF #########################


EnerDiag=0.
do i=1,Nx
do j=1,Nm
EnerDiag(i,j)=Ener(p)+imPart
enddo
enddo

write(*,*) "FermiL,FermiR",fermiL,fermiR
write(*,*) "SigInL,SigInR",CDABS(SUM(SigInL(p,:,:))),CDABS(SUM(SigInR(p,:,:)))




write(*,*) "#############################"





! ################################### TEST RECURSIVE ###############################

matSize=Nm*Nx
offSet=Nm


if(p.eq.1.and.iSCF.eq.1.and.cntAll.eq.0) then
allocate(Hreal(Nen,matSize,matSize))
allocate(HrealBall(Nen,matSize,matSize))
allocate(FullSigIn(matSize,matSize))
allocate(FullSigOut(matSize,matSize))
allocate(GnPrev(Nen,matSize,matSize))
allocate(GpPrev(Nen,matSize,matSize))
allocate(Gn(Nen,matSize,matSize))
allocate(Gp(Nen,matSize,matSize))
allocate(FullSigInEph(Nen,matSize,matSize))
allocate(FullSigInEphPrev(Nen,matSize,matSize))
allocate(FullSigOutEph(Nen,matSize,matSize))
allocate(FullSigREph(Nen,matSize,matSize))
allocate(HamEphFull(matSize,matSize))
allocate(TEMPOR(matSize,matSize))
FullSigInEph=0.
FullSigOutEph=0.
HamEphFull=0.
FullSigREph=0.
cntAll=1
endif

HamEphFull=0.
CntEntries=0
!do i=1,Nx
!offA=Nm*(i-1)
!if(i.gt.buffReg.and.i.le.Nx-buffReg-1) then
!TEMP1N=ModeHam0+VmatR(i,:,:)
!call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)
!HamEphBaseT=0.
!Elow=EnerProf(1)
!Ehigh=EnerProf(2)
!if(Ener(p-NphCov).le.Elow.and.Ener(p).ge.Ehigh) then
!HamEphBaseT=HamEph
!CntEntries=CntEntries+1
!endif
!if(Ener(p+NphCov).ge.Ehigh.and.Ener(p).le.Elow) then
!HamEphBaseT=HamEph
!CntEntries=CntEntries+1
!endif

do i=1,Nx
offA=Nm*(i-1)
if(i.gt.buffReg.and.i.le.Nx-buffReg-1) then
TEMP1N=ModeHam0+VmatR(i,:,:)
call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)
HamEphBaseT=0.
Elow=EnerProf(1)
Ehigh=EnerProf(2)
if(Ener(p-NphCov).le.Efl.and.Ener(p).ge.Efr) then
HamEphBaseT=HamEph
CntEntries=CntEntries+1
endif
if(Ener(p+NphCov).gt.Efr.and.Ener(p).le.Efl) then
HamEphBaseT=HamEph
CntEntries=CntEntries+1
endif


!                if(cntEntries.gt.0) then
!                     write(*,'(A24,10F15.4)') "Elow,Ehigh,Ep,E(p+NphCov)",Elow/ev2Hartree,Ehigh/ev2Hartree,Ener(p)/ev2Hartree,Ener(p+NphCov)/ev2Hartree
!                     pause
!               endif

HamEphFull((offA+1):offA+Offset,(offA+1):(offA+Offset))=HamEphBaseT
endif
enddo


if(iSCF.eq.1) then
HrealBall(p,:,:)=0.

HrealBall(p,1:Offset,1:Offset)=-Hm(1,:,:)-SigL(p,:,:)
HrealBall(p,1:offSet,(offset+1):2*OffSet)=-WmM(1,:,:)
HrealBall(p,((Nx-1)*Offset+1):Nx*Offset,((Nx-1)*Offset+1):Nx*Offset)=-Hm(Nx,:,:)-SigR(p,:,:)
HrealBall(p,((Nx-1)*Offset+1):Nx*Offset,((Nx-2)*Offset+1):(Nx-1)*Offset)=-WmMD(Nx-1,:,:)



write(*,*) "MatSize",matSize

do i=2,Nx-1
offA=Nm*(i-1)
HrealBall(p,(offA+1):offA+Offset,(offA+1):(offA+Offset))=-Hm(i,:,:)
HrealBall(p,(offA+1):offA+Offset,(offA+offset+1):(offA+2*offSet))=-WmM(i,:,:)
HrealBall(p,(offA+1):offA+Offset,(offA-offset+1):offA)=-WmMD(i-1,:,:)
enddo
endif


Hreal(p,:,:)=HrealBall(p,:,:)

do i=1,matSize
Hreal(p,i,i)=Hreal(p,i,i)+Ener(p)+imPart
enddo

Hreal(p,:,:)=Hreal(p,:,:)-FullSigREph(p,:,:)

write(*,*) "Before Inversion"

call getInverse(matSize,Hreal(p,:,:))

write(*,*) "After Inversion"

FullSigIn=0.
FullSigOut=0.

FullSigIn(1:Offset,1:Offset)=SigInL(p,:,:)
FullSigIn(((Nx-1)*Offset+1):Nx*Offset,((Nx-1)*Offset+1):Nx*Offset)=SigInR(p,:,:)


FullSigOut(1:Offset,1:Offset)=GamL*(1.-fermiL)*-im
FullSigOut(((Nx-1)*Offset+1):Nx*Offset,((Nx-1)*Offset+1):Nx*Offset)=GamR*(1-fermiR)*-im


FullSigIn=FullSigIn+FullSigInEph(p,:,:)
FullSigOut=FullSigOut+FullSigOutEph(p,:,:)

!  GnPrev(p,:,:)=FullSigInEph(p,:,:)
!  GpPrev(p,:,:)=FullSigOutEph(p,:,:)

Gn(p,:,:)=matmul(Hreal(p,:,:),matmul(FullSigIn,transpose(dconjg(Hreal(p,:,:)))))
Gp(p,:,:)=Hreal(p,:,:)-transpose(dconjg(Hreal(p,:,:)))+Gn(p,:,:)







if((absTerm.eq.1).and.(CntEntries.gt.0)) then
!       write(*,*) "E(p)",EnEv
!        pause
TEMPOR=0.
if(p.gt.NphCov.and.p+NphCov.le.Nen)   TEMPOR=Nph*Gn(p-NphCov,:,:)+Nph*Gn(p+NphCov,:,:)
if(p.le.NphCov) TEMPOR=Nph*Gn(p+NphCov,:,:)
if(p+NphCov.gt.Nen) TEMPOR=Nph*Gn(p-NphCov,:,:)
FullSigInEph(p,:,:)=matmul(HamEphFull,matmul(TEMPOR,HamEphFull))

endif





if((absTerm.eq.1).and.(CntEntries.gt.0)) then
TEMPOR=0.
if(p.gt.NphCov.and.p+NphCov.le.Nen)   TEMPOR=Nph*Gp(p+NphCov,:,:)+Nph*Gp(p-NphCov,:,:)
if(p.le.NphCov) TEMPOR=Nph*Gp(p+NphCov,:,:)
if(p+NphCov.gt.Nen) TEMPOR=Nph*Gp(p-NphCov,:,:)
FullSigOutEph(p,:,:)=matmul(HamEphFull,matmul(TEMPOR,HamEphFull))

endif






if((SpEmTerm.eq.1)) then
TEMPOR=0.
 do i=1,Nen-p
   HamEphFull=0.
     do j=1,Nx
       offA=Nm*(j-1)
          if(j.ge.buffReg.and.j.le.Nx-buffReg) then
           TEMP1N=ModeHam0+VmatR(j,:,:)
           call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)
           HamEphBaseT=0.
           Elow=EnerProf(1)
           Ehigh=EnerProf(2)
            if(Ener(p).le.Elow.and.Ener(p+i).ge.Ehigh)  HamEphBaseT=HamEph
              HamEphFull((offA+1):offA+Offset,(offA+1):(offA+Offset))=HamEphBaseT
           endif
  enddo


    if(SUM(CDABS(HamEphFull)).gt.0)  then
!                                     FullSigInEph(p,:,:)=FullSigInEph(p,:,:)+matmul(HamEphFull,matmul(Gn(p+i,:,:),HamEphFull))*dE
     TEMPOR=TEMPOR+matmul(HamEphFull,matmul(Gn(p+i,:,:),HamEphFull))*dE
!                                       FullSigOutEph(p+i,:,:)=FullSigOutEph(p+i,:,:)+matmul(HamEphFull,matmul(Gp(p,:,:),HamEphFull))*dE
    endif

   enddo

FullSigInEph(p,:,:)=FullSigInEph(p,:,:)+TEMPOR


TEMPOR=0.
   do i=1,p-1
      HamEphFull=0.
         do j=1,Nx
         offA=Nm*(j-1)
          if(j.ge.buffReg.and.j.le.Nx-buffReg) then
            TEMP1N=ModeHam0+VmatR(j,:,:)
            call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)
            HamEphBaseT=0.
            Elow=EnerProf(1)
            Ehigh=EnerProf(2)
              if(Ener(p-i).le.Elow.and.Ener(p).ge.Ehigh)  HamEphBaseT=HamEph
                 HamEphFull((offA+1):offA+Offset,(offA+1):(offA+Offset))=HamEphBaseT
              endif
    enddo


   if(SUM(CDABS(HamEphFull)).gt.0)  then
!                                     FullSigInEph(p,:,:)=FullSigInEph(p,:,:)+matmul(HamEphFull,matmul(Gn(p+i,:,:),HamEphFull))*dE
   TEMPOR=TEMPOR+matmul(HamEphFull,matmul(Gp(p-i,:,:),HamEphFull))*dE
  endif
  enddo


FullSigOutEph(p,:,:)=FullSigOutEph(p,:,:)+TEMPOR


endif




FullSigREph(p,:,:)=0.5D0*(FullSigOutEph(p,:,:)-FullSigInEph(p,:,:))



do i=1,Nx-1
offA=Nm*(i-1)
gUp(i,:,:)=Hreal(p,(offA+1):offA+Offset,(offA+offset+1):(offA+2*offSet))
gLow(i,:,:)=Hreal(p,(offA+offset+1):(offA+2*offSet),(offA+1):offA+Offset)
gNup(i,:,:)=Gn(p,(offA+1):offA+Offset,(offA+offset+1):(offA+2*offSet))
gNlow(i,:,:)=Gn(p,(offA+offset+1):(offA+2*offSet),(offA+1):offA+Offset)
gPup(i,:,:)=Gp(p,(offA+1):offA+Offset,(offA+offset+1):(offA+2*offSet))
gPlow(i,:,:)=Gp(p,(offA+offset+1):(offA+2*offSet),(offA+1):offA+Offset)
enddo

do i=1,Nx
offA=Nm*(i-1)
gDiag(i,:,:)=Hreal(p,(offA+1):offA+Offset,(offA+1):(offA+Offset))
gNdiag(i,:,:)=Gn(p,(offA+1):offA+Offset,(offA+1):(offA+Offset))
gPdiag(i,:,:)=Gp(p,(offA+1):offA+Offset,(offA+1):(offA+Offset))
enddo

write(*,*) "Got gDiag ... etc"


!#################################################################### End Eph coupling ######################3





Amat11=im*(gDiag(1,:,:)-transpose(dconjg(gDiag(1,:,:))))
TEMP1=Amat11-matmul(gDiag(1,:,:),matmul(GamL,transpose(dconjg(gDiag(1,:,:)))))
TEMP2=matmul(GamL,TEMP1)

currDensCB=0.
currDensVB=0.


do j=1,Nx-1
TEMP1=matMul(Wm,DBLE(gNLow(j,:,:)))
TEMP3=matMul(Wm,DBLE(gPLow(j,:,:)))
summCmplx=0.
summCmplxHole=0.

do jj=1,Nm
summCmplx=summCmplx+TEMP1(jj,jj)
summCmplxHole=summCmplxHole+TEMP3(jj,jj)
enddo

TEMP1N=ModeHam0+VmatR(j,:,:)
call   getDiagonalized(kpSize*Nk,TEMP1N,EnerProf)
Elow=EnerProf(1)
Ehigh=EnerProf(2)

if(Ener(p).ge.0.35*ev2Hartree) then
currDensCB(j)=dble(summCmplx)
else
currDensVB(j)=dble(summCmplx)
endif

currDens(j)=dble(summCmplx)
currDensHole(j)=dble(summCmplxHole)

currElIn(j)=currElIn(j)+currDens(j)*dE
currHIn(j)=currHin(j)+currDensHole(j)*dE
currElCB(j)=currElCB(j)+currDensCB(j)*dE
currElVB(j)=currElVB(j)+currDensVB(j)*dE
enddo


TEMP1=im*matmul(gDiag(1,:,:)-transpose(dconjg(gDiag(1,:,:))),SigInL(p,:,:))-matmul(gNDiag(1,:,:),GamL)

summ=0.

do j=1,Nm
summ=TEMP2(j,j)+summ
enddo

T(p)=summ

summ=0.
summ1=0.
summ2=0.
summEl=0.
summCur=0.
summEl2=0.
summEl1=0.

do i=1,Nx
do j=1,Nm
summ=summ+DIMAG(gDiag(i,j,j))
summEl=summEl+dimag(gNdiag(i,j,j))
ElDensityLayer(i)=ElDensityLayer(i)+2.*dimag(gNdiag(i,j,j))/(2.*pi)*dE
enddo
enddo

do i=1,Nm
summ1=summ1+DIMAG(gDiag(2,i,i))
summ2=summ2+DIMAG(gDiag(Nx-1,i,i))
summEl1=summEl1+dimag(gNdiag(Nx/2,i,i))
summEl2=summEl2+dimag(gNdiag(Nx-buffreg+1,i,i))

enddo


do j=1,Nx
summa=0.
do jj=1,Nm
summa=summa+dimag(gNdiag(j,jj,jj))
enddo
write(27,'(10F22.12)') j*dx/18.89, Ener(p)/ev2Hartree, summa

enddo

write(23,*) Ener(p)/ev2Hartree,T(p)
write(22,'(10F22.12)') Ener(p)/ev2Hartree,-summ*2.,-summ1*2.,-summ2*2.
write(26,'(10F22.12)') Ener(p)/ev2Hartree,summEl*2.,summEl1*2.,summEl2*2,summEl/summEl1
write(28,'(10F22.12)') Ener(p)/ev2Hartree,currAlt,currDens(2),currDens(buffReg+(Nx/2-buffreg)/2),currDens(Nx/2),currDens(buffReg+(Nx/2+buffreg)/2),currDens(Nx-buffReg+1)


if(iSCF.gt.1) then
currDensPhoto=currDens-currDensBall(:,p)
write(37,'(10F22.12)') Ener(p)/ev2Hartree,currDensPhoto(2),currDensPhoto(buffReg+(Nx/2-buffreg)/2),currDensPhoto(Nx/2),currDensPhoto(buffReg+(Nx/2+buffreg)/2),currDensPhoto(Nx-buffReg+1)
do j=2,Nx-1
write(52,'(10F22.12)') j*dx/18.89, Ener(p)/ev2Hartree, currDensPhoto(j)
enddo
endif

write(29,'(10F22.12)') Ener(p)/0.036749309,currDensHole(2),currDensHole(buffReg+(Nx/2-buffreg)/2),currDensHole(Nx-buffReg+1)



if(iSCF.eq.1) then
currDensBall(:,p)=currDens
currDensHoleBall(:,p)=currDensHole

endif




enddo


write(27,*)
write(23,*)
write(22,*)
write(26,*)
write(28,*)
write(52,*)
write(37,*)
write(31,'(I5,100F22.12)') iscF,(currElIn(j),j=2,Nx-1)

write(38,'(I5,100F22.12)') iscF,(currHin(j),j=2,Nx-1)

Open(unit=25, file="CurrInLayerElCBVB.dat", status="unknown")
Open(unit=40, file="CurrInLayerEl.dat", status="unknown")
Open(unit=41, file="CurrInLayerHole.dat", status="unknown")
Open(unit=42, file="CurrInLayerSum.dat", status="unknown")
open(unit=46, file="ElDensityLayer.dat", status="unknown")


do j=1,Nx-1
write(25,'(1I4,3F24.12)') j,currElCB(j),currElVB(j),currElCB(j)+currElVB(j)
write(40,*) j,currElIn(j)
write(41,*) j,currHin(j)
write(42,*) j,currElIn(j)+currHin(j)
write(46,*) j,ElDensityLayer(j)

enddo





close(46)
close(40)
close(41)
close(42)
close(25)


if(iSCF.gt.2)  then
if(mod(iSCF,2).ne.0) then
CurPrev=sum(CurrElIn)
else
CurCur=sum(CurrElIn)
endif
endif


if(iSCF.ge.4) then
errorCur=dabs(CurCur-CurPrev)/CurCur*100
write(*,*) "ErrorCurr",errorCur
if(errorCur.le.ScfPerc)  goto 122
endif

enddo

122 continue

deallocate(gDiag,gUp,gLow,gRTEMP)
deallocate(gLTEMP,TEMP,TEMP1,TEMP2)
deallocate(CMAT,Wm,WmMD,WmM,EigVctr0,EigVctr1,Hm,Ediag,EnerDiag,EigVctrMode,EigVctrMode1)
deallocate(gL,gR,SigL,SigR,TempEn,TEMP3,TEMP4,TEMP5)




deallocate(Vmat,VmatR,SigInL,SigInR)


deallocate(gNdiag,gNlow,gNup,gPdiag,gPlow,gPup)


end subroutine getFullMatrix










subroutine  getInverse(N,Mat)
implicit none
integer N,INFO
complex*16 Mat(N,N),WORK(N)
integer IPIV(N)

call zgetrf(N,N,Mat,N,IPIV,INFO)
if(INFO.ne.0) write(*,*) "###############  ERROR IN ZGETRF in getInverse ###############",INFO
call zgetri(N,Mat,N,IPIV,WORK,N,INFO)
if(INFO.ne.0) write(*,*) "###############  ERROR IN ZGETRI in getInverse ###############",INFO

end subroutine getInverse


subroutine  getDiagonalized(N,Mat,EigEn)
implicit none
integer N,INFO
complex*16 Mat(N,N)
real*8 EigEn(N)
integer,allocatable :: IWORK(:)
real*8, allocatable :: RWORK(:)
complex*16, allocatable :: WORK(:)
integer :: LRWORK,LWORK,LIWORK



LRWORK=1+5*N+2*N**2.
LWORK=LRWORK
LIWORK=3+5*N
allocate(RWORK(LRWORK),IWORK(LIWORK),WORK(LWORK))
call zheevd('V','U',N,Mat,N,EigEn,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
if(INFO.ne.0) write(*,*) "###############  ERROR IN ZHEEVD in getDiagonalized  ###############",INFO
deallocate(RWORK,IWORK,WORK)


end subroutine getDiagonalized




subroutine GetConactSelfEnergies(p,Ham,Nm,gPrev,Ediag,WmM,WmMD,gFc,Nmix,Nmod,alpha,NitMax,Fname,eps,flag,Niter)
implicit none
integer ::  p,Nm,i,j,q,Nmix,Nmod,NitMax,ind,ij,iq,flag,Niter
real*8 :: alpha,eps,error
complex*16 :: rhoN(Nmix+4,Nm,Nm),Fn(Nmix+4,Nm,Nm),diffRhoFn(Nmix+4,Nm,Nm),diffFnFnm1(Nmix+4,Nm,Nm)
complex*16 :: TEMP(Nm,Nm),TEMP3(Nm,Nm),TEMP4(Nm,Nm)
complex*16 :: Ham(Nm,Nm),gPrev(Nm,Nm),WmMD(Nm,Nm),WmM(Nm,Nm),Ediag(Nm,Nm),gFc(Nm,Nm),Unm(Nmix,Nmix),An(Nmix),Vn(Nmix)

character(len=4) :: Fname





if(p.eq.1) then
TEMP=Ham
call getInverse(Nm,TEMP)
gFc=TEMP
endif


if(p.gt.1) then
gFc=gPREV
endif

do i=1,Nitmax
TEMP=Ediag-Ham-matmul(WmMD,matmul(gFc,WmM))
call getInverse(Nm,TEMP)
error=maxval(cdabs((TEMP-gFc)))
gFc=TEMP*(1.-alpha)+alpha*gFc
if(error.le.eps) then
write(*,*) "The",Fname, " at error:", error,"converged",p,i
exit
endif
enddo

end subroutine GetConactSelfEnergies




subroutine GetUmat(Nk,KnY,KnZ,Ky,Kz,Ly,Lz,Umat,Ny,Nz,dy,dz,kpSize)
implicit none
integer :: Ny, Nz,indX,indY,pp,qp,j1,i2
integer :: Nk,KnY(Nk),KnZ(Nk),kpSize,p,q,m,n,colN,rowN,offX,offY,i,j,l,i1
real*8 :: Ky(Nk),Kz(Nk),Ly,Lz,Umat(kpSize*Ny*Nz,kpSize*Nk),dy,dz,kp,kq
real*8 :: ym(Ny),zn(Nz)
real*8 kyy,kzz,kyyp,kzzp,summ,factor1,factor2


Umat=0.

do j=1,Nz
zn(j)=dz*(j-1)
enddo


do j=1,Ny
ym(j)=dy*(j-1)
enddo

do i=1,Ny
do j=1,Nz
indY=(i-1)*Nz+j
do l=1,Nk
indX=l
kp=Ky(l)
kq=Kz(l)
do i1=1,kpSize
offX=(indY-1)*kpSize+i1
offY=(indX-1)*kpSize+i1
Umat(offX,offY)=2./dsqrt(dble(Ny*Nz))*dsin(kp*ym(i))*dsin(kq*zn(j))
!                if(kp.eq.0.and.kq.eq.0) Umat(i,i)=1.
enddo
enddo
enddo
enddo


do i=1,kpSize
Umat(i,i)=1.
enddo

end subroutine GetUmat





subroutine GetVpot(Vmat,Nx,Ny,Nz,kpSize,Efl,Efr,buff,Dir,dx,dEfr)
implicit none
integer Nx,Ny,Nz,kpSize,i,j,k,buff,ind,i1,Nread
real*8 dV,V,Efl,Efr,Dir,x0,V0,ev2Hartree,nano2bohr,Value,xx
real*8 Vmat(Nx,kpSize*Ny*Nz,kpSize*Ny*Nz),err,dx,dEfr
character*12 fileNm


real*8,allocatable :: Vread(:),Xread(:)
parameter(ev2Hartree=0.036749309)
parameter(nano2bohr=18.89)

fileNm="pot10.dat"

Open(unit=40, file=fileNm, status="unknown")

read(40,*) Nread,dEfr
read(40,*) x0,V0
close(40)

Open(unit=40, file=fileNm, status="unknown")
read(40,*) Nread,dEfr

allocate(Vread(Nread),Xread(Nread))

do i=1, Nread
read(40,*) Xread(i),Vread(i)
enddo



Xread=(Xread-X0)*nano2bohr
Vread=(Vread-V0)*ev2Hartree


do i=1,Nread
write(*,*) Xread(i)/nano2bohr,Vread(i)/ev2Hartree
enddo
!pause

Vmat=0.
dv=(Efl-Efr)/(Nx-2*buff-2)



do i=1,buff
ind=0
do j=1,Ny
do k=1,Nz
do i1=1,kpSize
ind=ind+1
Vmat(i,ind,ind)=Efl-Efl
enddo
enddo
enddo
enddo

do i=Nx-buff,Nx
ind=0
do j=1,Ny
do k=1,Nz
do i1=1,kpSize
ind=ind+1
Vmat(i,ind,ind)=Efl-Efr
enddo
enddo
enddo
enddo




do i=buff+1,Nx-buff-1
ind=0
do j=1,Ny
do k=1,Nz
do i1=1,kpSize
ind=ind+1
Vmat(i,ind,ind)=dV*(i-buff-1)
enddo
enddo
enddo
enddo


do i=buff+1,Nx-buff-1
ind=0
do j=1,Ny
do k=1,Nz
do i1=1,kpSize
ind=ind+1
!          if(i.ge.buff+5.and.i.le.Nx-buff-6) then
!                  if(i1.eq.1)      Vmat(i,ind,ind)=Vmat(i,ind,ind)+0.1*0.036749309
!                  if(i1.eq.2)      Vmat(i,ind,ind)=Vmat(i,ind,ind)-0.1*0.036749309
!          endif
enddo
enddo
enddo
enddo



do i=1,Nx
ind=0
do i1=1,kpSize
ind=ind+1
xx=(i-1)*dx

call polint(Xread,Vread,Nread,xx,Value,err)
write(*,*) xx/nano2bohr,Value
Vmat(i,ind,ind)=-Value
enddo
enddo


dEfr=dEfr*ev2Hartree
!Vmat=Vmat*-1.


end subroutine GetVpot










SUBROUTINE polint(xa,ya,n,x,y,dy)
implicit none
INTEGER n,NMAX
REAL*8 dy,x,y,xa(n),ya(n),dx
PARAMETER (NMAX=2)
INTEGER i,m,ns
REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
ns=1
dif=dabs(x-xa(1))

do i=1,n
dift=dabs(x-xa(i))
if (dift.lt.dif) then
ns=i
dif=dift
endif
c(i)=ya(i)
d(i)=ya(i)
enddo
y=ya(ns)
ns=ns-1

dx=xa(n)-xa(n-1)
y=(ya(ns+1)-ya(ns))/dx*(x-xa(ns))+ya(ns)

!  do  m=1,n-1
!     do  i=1,n-m
!        ho=xa(i)-x
!        hp=xa(i+m)-x
!        w=c(i+1)-d(i)
!       den=ho-hp
!          if(den.eq.0.) pause 'failure in polint'
!       den=w/den
!       d(i)=hp*den
!       c(i)=ho*den
! enddo
!  if (2*ns.lt.n-m) then
!  dy=c(ns+1)
!   else
!    dy=d(ns)
!    ns=ns-1
!  endif
! y=y+dy
! enddo
!  return

END subroutine polint


















