program galaxy

  implicit none

  double precision,parameter :: tend = 10e3
  double precision,parameter :: yret = 0.1
  double precision,parameter :: yz = 0.099
  double precision,parameter :: ydust = 0.001
  double precision :: tnow,dt
  double precision :: mgas,mstar,mz,mdust,mztot
  double precision :: zdust,zgas,ztot
  double precision :: dmgas,dmstar,dmz,dmdust
  double precision :: getsfr,sfr
  double precision :: dsgas,dsz,dsdust,dism
  double precision :: grow,destroy,gettimestep
  double precision :: grate,drate,safeadd
  double precision :: stardust,totdust
  integer :: i

  tnow = 0.

  ! initialise mgas, mstar, mz and mdust

  mgas = 1.e10
  mstar = 0.
  mz = 0.
  mdust = 0.
  sfr = 0.
  zdust = 0.
  zgas = 0.

  stardust = 0.
  totdust = 0.

  do
     ! get star formation/destruction/growth rates
     sfr = getsfr(tnow,mgas)
     drate = destroy(zgas,zdust,mgas,sfr)
     grate = grow(zgas,zdust,mgas,sfr)
     ! get timestep
     dt = gettimestep(sfr,drate,grate,mgas,mz,mdust)
     ! determine returned gas/metal/dust fractions
     dmstar = (1.-yret)*sfr*dt
     dmgas = -dmstar
     dsgas = yret*sfr*dt
     dsz = (yz*yret - zgas)*sfr*dt
     dsdust = (ydust*yret - zdust)*sfr*dt
     ! determine growth/destruction of preexisting dust
     dmz = dsz + (drate - grate)*dt
     dmdust = dsdust - (drate - grate)*dt
     stardust = stardust + ydust*dsgas
     totdust = totdust + ydust*dsgas + grate*dt
     ! update masses
     mgas = safeadd(mgas,dmgas)
     mstar = safeadd(mstar,dmstar)
     mz = safeadd(mz,dmz)
     mdust = safeadd(mdust,dmdust)
     zgas = mz/mgas
     zdust = mdust/mgas
     mztot = mz + mdust
     ztot = zgas+zdust
     tnow = tnow + dt
     write(*,'(1(ES10.3,2X),4(F9.5,2X),2(ES9.2,2X))') tnow/1e3,ztot/0.0134,zgas/ztot,zdust/ztot,stardust/totdust,mz/grate,mdust/drate
     if (tnow .ge. tend) exit
  end do

end program galaxy

double precision function getsfr(t,mgas)

  implicit none

  double precision,intent(in) :: t,mgas

  getsfr = 2e6*(mgas/1e10)

  if (mgas .le. 0.) getsfr = 0.

end function getsfr

double precision function grow(zgas,zdust,mgas,sfr)

  implicit none

  double precision,parameter :: mp = 1.67e-24
  double precision,parameter :: kb = 1.38e-16
  double precision,parameter :: msun = 1.99e33
  double precision,parameter :: year = 3.15e7
  double precision,intent(in) :: zgas,zdust,mgas,sfr
  double precision :: amass,gtemp,gdens,mu,fdense
  double precision :: vel,adens
  double precision :: dustarea

  ! gas/atom properties
  mu = 2.33*mp
  amass = 12.*mp
  gtemp = 20.
  gdens = 100.
  fdense = 0.015 !* (sfr/2e6)

  ! velocity/density of metals
  vel = sqrt(kb*gtemp/amass)
  adens = mu * gdens * zgas

  ! dust grain area / mass
  dustarea = 4.24e5 ! MRN cm^2 / g
  
  ! growth rate dm/dt ~ area * velocity * density
  ! velocity ~ constant
  ! density ~ zgas
  ! area ~ zdust(?)
  grow = dustarea * zdust*mgas*fdense * vel * adens
  grow = grow * 1e6*year

end function grow

double precision function destroy(zgas,zdust,mgas,sfr)

  implicit none

  double precision,intent(in) :: zgas,zdust,mgas,sfr
  double precision :: fsn
  double precision :: snrate,mshocked,fdest,mdest

  ! number of sn per msun
  fsn = 0.01
  
  snrate = fsn*sfr

  ! dust mass shocked per sn
  mshocked = 1000 * zdust
  mshocked = 8143 * zdust
!  mshocked = 4036 * zdust

  ! destruction efficiency
  fdest = 1.
  fdest = 1. / (1. + zgas/(0.0134*0.14))
!  fdest = 1. / (1. + zgas/(0.0134*0.039))**(0.298)

  ! dust mass destroyed per sn
  mdest = fdest * mshocked

  ! destruction rate dm/dt ~ mdest * snrate
  ! mdest ~ ???? zdust, zgas
  ! mdest ~ mshocked * fdest
  ! snrate ~ sfr
  destroy = mdest * snrate

end function destroy

double precision function safeadd(m1,m2)

  implicit none

  double precision,intent(in) :: m1,m2
  double precision :: add

  add = m1 + m2
  
  if (add .ge. 0.) then
     safeadd = add
  else
     safeadd = 0.
  end if

end function safeadd

double precision function gettimestep(sfr,drate,grate,mgas,mz,mdust)

  implicit none

  double precision,intent(in) :: sfr,drate,grate,mgas,mz,mdust
  double precision,parameter :: tol = 0.01
  double precision :: tstar,tdest,tgrow

  tstar = mgas/sfr
  tdest = mdust/drate
  tgrow = mz/grate

  gettimestep = tol*min(tstar,tdest,tgrow)

  if (isnan(gettimestep)) gettimestep = tol*tstar

end function gettimestep
