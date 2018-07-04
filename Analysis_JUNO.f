
!! 2017/08/23 interpolate the integrated spectra (t,Enu) numerical table
!! 2017/08/24 using PDG formulae (39.78) to calc the expected significance
!! Backgournd of JUNO, spectra maybe??? (mainly reactor neutrios + geo-neutrinos)    
!! considering background uncertainty 
!! 2017/08/28 New reactor & geo-neutrino    
!! JUNO, average efficiency about 0.73;   

!!! Re-scan to find the best time & energy window, when telling mass hierarchies        

!!! 17/09/22 calc p-value & the expected sensitivity 

!!! 17/09/26 using MC simulation      (not adopted)

!!! 17/11/16 focus on advance warning         
      
       module global 
       real*8:: me = 0.511, pi = 3.14159, npyr=7.23d31, Dl=0.2*3.086d21  ! H mass fraction, 12%
       integer,Parameter:: n=701,m=80          !  12,15,20,25 (/1688, 662, 513, 701/)     
       integer:: ierr,il=1, jl=1, ix1(2,m),iy1(2,n),      
     C         ix2(2,m),iy2(2,n)  
       real*8:: a1(4,4,n,m),a2(4,4,n,m),dx(n),dy(m)      
       real*8:: derx,dery,pp = 0.1, vlow = -50.,vint(1:2),  
     C      Event(1:2,1001), Enu(m), Time(n) 
       real*8:: Epk(1:4), tpk(1:4)   ! peak time 
       integer:: StelN    ! stellar model 
       real*8:: enu_ron(153),enu_roff(43),enu_geo1(9),enu_geo2(6),
     C      enu_geo3(12),ron(153),roff(43),geo1(9),geo2(6),geo3(12), 
     C      y2_ron(153),y2_roff(43),y2_geo1(9),y2_geo2(6),y2_geo3(12)   
       real*8:: alpha_sys = 0.1,alpha_sysb=0.1  ! systematic error 
       real*8:: sigalpha = 0.1  ! same as alpha_sys, but somehow used for p_value calculation      
       end module      
   

       program test        
       use global
       implicit none
       real*8, external:: S_b ,Event_react,Event_react_2013, 
     C  Event_geo,CroX, fb,fb2,S_b2,S_b3,signif2, 
     C  Event_geoN, Event_ron, Event_roff, S_b_sys,gammln,pval,likeli2,
     C  xN_3sig_up,xN_3sig_low,pval_low,zbqlnor,p_MC,NH_MC,likeli, 
     C  likeli3,fbb,fb2b,likeli_bk,likeli_bk0,likeli_bk00     

       integer,external:: zbqlpoi
       integer::i,j,k,iT,iE,iiT,ii                 
       real*8:: Spec(1:2,n,200:279)       !  vebar rate: e/mu, T, rhoY  
       real*8:: SpecInt1(n,m),SpecInt2(n,m)      ! for interpolation input  
!       integer:: NT(1:4)  
       real*8::  Tup,Tlow,tmp,Emax,Emin,event_tot(1:7),alpha_sc(10) 
       real*8::xEnu,tmp2,tmp3,tmp4,tmp5,tmp6,xT,xdl,param(0:2),
     C     xevent(1:2),tmp7(100),xT1,xT2   
       character(len=10):: cha(4)
!       real*8:: Epk(1:4), tpk(1:4)   ! peak time 
!       integer:: StelN    ! stellar model  

!       NT2 = (/33887,23620,14130,17005/)     ! Start of Snapshot   
!       NT1 = (/17000,17000,9000,10000/)      ! End of Snapshot    
!       NT(1:4) = (/1688, 662, 513, 701/)      ! stellar models: 12,15,20,25 
!       N = NT(1)    
!c=====================================================================
!c     Read the neutrino spectrum for each thermal process          
!c=====================================================================  
       cha(1:4) = (/'12','15','20','25'/)         
       Spec = 0.                     
       open(31,file='S25_Spectra(tE).dat',status = 'old')       !  12,15,20,25     
       StelN = 4                                                 !  12,15,20,25   
       do i=1,n
        read(31,*) Time(i),Spec(1,i,200:279)       ! eb 
        read(31,*) tmp,Spec(2,i,200:279)           ! xb  
!        write(*,*) 'hah',i,time(i),-log10(Time(i))  
        Time(i) = -log10(Time(i))       
        Spec(1:2,i,200:279) = log10(max(Spec(1:2,i,200:279), 1.d-50 ))  
!        write(*,*) i, Spec(1,i,220)     
!        write(*,*) i,Time(i)    
       end do  
       close(31)
       do i=1,m 
        Enu(i) = 10.**((i-1)*0.02)    
!        write(*,*) i, Enu(i) 
       end do        
       write(*,*) ' finishing reading star_integrated spectra of '
     C     ,cha(stelN)   
 
!!!!===================================================================
!!!      Read reactor & geo-neutrino data for interpolation   
!!!====================================================================

      open(32, file = './Pre-neutrino-figs-paper/Reactor-On.dat', 
     C    status = 'old') 
      open(33, file = './Pre-neutrino-figs-paper/Reactor-Off.dat',
     C    status = 'old') 
      open(34, file = './Pre-neutrino-figs-paper/Geo-1.dat',
     C    status = 'old') 
      open(35, file = './Pre-neutrino-figs-paper/Geo-2.dat',
     C    status = 'old') 
      open(36, file = './Pre-neutrino-figs-paper/Geo-3.dat', 
     C    status = 'old')       

      do i=1,153    ! 9.6 
      read(32,*) enu_ron(i), ron(i)
      end do
      do i=1,43     ! 9. 
      read(33,*) enu_roff(i), roff(i)
      end do
      do i=1,9    ! 2.06146 
      read(34,*) enu_geo1(i), geo1(i) 
      end do
      do i=1,6    ! 2.07348,2.24048 
      read(35,*) enu_geo2(i), geo2(i)
      end do
      do i=1,12  ! 2.27655 - 3.27 
      read(36,*) enu_geo3(i), geo3(i)
      end do 
      call init_back()
 
!c=====================================================================
!c        interpolation the star-integrated spectra 
!c=====================================================================  
! initialization 
      Specint1(1:n,1:m) = Spec(1,1:n,200:279)
      Specint2(1:n,1:m) = Spec(2,1:n,200:279)  
      call rattwo(n,m,Time,Enu,Specint1,vlow,pp,dx,dy,ix1,iy1,a1) 
      call rattwo(n,m,Time,Enu,Specint2,vlow,pp,dx,dy,ix2,iy2,a2)  
!c=====================================================================
!c                    End Interpolation
!c=====================================================================                    
      Tup = 0.1       ! set time window 
      Tlow = 0.01 
      Emax = 4.0      !Enu(m)    ! set energy window  
      Emin = 1.804   
       
!!!! Draw observed spectra, get events rates        
!      open(51,file = 'Observed_Events25.dat')  
      do i = 1, -200    
      tmp = 1.804d0 + 0.02*i  
      call CalcEvent(1.d-6,1.d0,tmp-0.02,tmp,Event_tot)  
      tmp2 = Event_tot(1)*0.681 + Event_tot(2)*0.319 
      tmp2 = tmp2/0.02*0.73         ! per MeV  
      tmp3 = Event_tot(1)*0.022 + Event_tot(2)*0.978   
      tmp3 = tmp3/0.02*0.73
      tmp4 = Event_tot(3)/0.02*0.73   
      tmp5 = Event_tot(4)/0.02*0.73
      tmp6 = Event_tot(5)/0.02*0.73  
!      write(*,*) i,'hah' 
      write(*,'(10G14.5)') tmp,tmp2,tmp3,tmp4,tmp5,tmp6  
      end do 
!      close(51) 
 
      call CalcEvent(1.d-6,1.d0,Emin,4.d0,Event_tot)   
      tmp2 = Event_tot(1)*0.681 + Event_tot(2)*0.319 
      tmp2 = tmp2*0.73         ! per MeV  
      tmp3 = Event_tot(1)*0.022 + Event_tot(2)*0.978 
      tmp3 = tmp3*0.73
      tmp4 = Event_tot(3)*0.73   
!      tmp5 = Event_tot(4)*0.73
      tmp6 = Event_tot(5)*0.73  
!      write(*,'(10G14.5)') tmp2,tmp3,tmp4,tmp6  
      
      do i=1,-200     ! 1.d-3 to 10 days 
      tmp = 1.d-3*10.**(log10(10./1.d-3)/200.*i)    
      call CalcRate(tmp,Emin,4.d0,xEvent)  
      write(*,'(10G14.5)') tmp,xEvent(1)*0.73,xEvent(2)*0.73, 
     C    xEvent(1)/xEvent(2)     
      end do

!c=====================================================================  
!!  for certain event rate, study the significance as a fucntion of 
!!  uncertainty      
!c===================================================================== 
       
      call CalcEvent(1.d-6,1.d0,Emin,4.d0,Event_tot)  
      tmp = (Event_tot(3) + Event_tot(5))*0.73  ! background rate       
      alpha_sc(1:10)=(/0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)  
!      open(77,file='Significance_mu_alpha_2DGrid_IH.dat')  
      do j=1,-101
!      tmp2 = j*10.    ! scan 10 - 1000 events 
      tmp2 = 10.*10.**(0.0230103*(j-1.))     ! 10 -2000   
      do i=1,76       ! 0.05 - 0.8 with step size 0.01  
!      alpha_sys = alpha_sc(i)
      alpha_sys = 0.05 + (i-1.)*0.01    ! 0.01 - 0.8    
      tmp7(i) = likeli2(tmp,tmp2,tmp2/3.4d0)   ! NH
!      tmp7(i) = likeli2(tmp,tmp2/3.4d0,tmp2)    ! IH   
      end do    
      write(*,*) 'finish',j,'of 100'  
      write(*,'(100G14.5)') tmp2,tmp7(1:76)  
!      write(77,'(100G14.5)') tmp2,tmp7(1:76)      
      end do    
!      close(77)        

!!!! compare (NH,20) with (IH 25)
!!!! compare (IH,20) with (NH, 15) 
!!!! 15/20/25: 12.3(3.7),21.0(6.1),25.0(7.2) 
      if(10.>1.) then

      alpha_sys  = 0.1
      tmp2 = 21.0*25.   ! 200 pc
      tmp3 = 7.2*25. 
      write(*,*) likeli2(tmp,tmp2,tmp3)  
      tmp2 = 6.1*25. 
      tmp3 = 12.3*25.   
      write(*,*) likeli2(tmp,tmp2,tmp3)  

      alpha_sys  = 0.18  
      tmp2 = 21.0*25.   ! 200 pc
      tmp3 = 7.2*25. 
      write(*,*) likeli2(tmp,tmp2,tmp3)  
      tmp2 = 6.1*25. 
      tmp3 = 12.3*25.    
      write(*,*) likeli2(tmp,tmp2,tmp3) 

      write(*,*) 'check likeli2 Vs 3 below:'
      do i=1,10 
      tmp = 10.+i*2.
      tmp2 = tmp/3.5 
      write(*,*) tmp,likeli2(17.d0,tmp,tmp2),likeli3(17.d0,tmp,tmp2) 
      end do 

      end if 

!c=====================================================================  
!       scan time-window  (1.,0.2,0.3,0.3) days    
!c=====================================================================  
      Tpk(1:4) = (/1.d0,0.2d0,0.3d0,0.3d0/)   
!      Tpk(1:4) = (/0.3d0,0.3d0,0.3d0,0.3d0/)  
      xdl = 0.5   ! kpc   
!      open(55, file='Significance_12.dat')     ! New E&T-window

!      alpha_sc(1:10)=(/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5/) 
      alpha_sc(1:10)=(/0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8/) 
      do i=2,-2        
      tmp2 = 0.  
      alpha_sys = alpha_sc(i)    ! scan alpha, flux uncertainty   
      do j=1,96   ! scan time window    
      xdl = 0.1 + (j-1)*0.02       ! 0.2 - 2.0         
!      call CalcEvent(1.d-6,0.05d0*j,1.804d0,4.d0,Event_tot)  
      call CalcEvent(1.d-6,Tpk(StelN),1.804d0,4.d0,Event_tot)        
      tmp = (Event_tot(1)*0.022 + Event_tot(2)*0.978)*0.73       ! IH 
      tmp3 =(Event_tot(1)*0.681 + Event_tot(2)*0.319)*0.73       ! NH 
      tmp4 = (Event_tot(3)+Event_tot(5))*0.73
!      write(*,'(10G14.5)') 'T,bcg,IH,NH',0.02*j,tmp4,tmp,tmp3 
      tmp = tmp/xdl**2   
      tmp3 = tmp3/xdl**2 
!      write(*,*) tmp4,tmp,tmp3
      tmp5 = likeli2(tmp4,tmp,tmp3)   ! IH  
      tmp6 = likeli2(tmp4,tmp3,tmp)   ! NH events, how likely fitted by IH
!      write(*,*) 0.05*j,tmp5,tmp6 
!      if(tmp6>tmp2) then 
!          tmp2 = tmp6  
!      else 
!       goto 23    
!      end if    
      write(*,*) xdl,tmp5,tmp6    

!      tmp5 = S_b_sys((Event_tot(3)+Event_tot(5) 
!     C     +tmp)*0.73d0,1) 
!      tmp6 = -S_b_sys((Event_tot(3)+Event_tot(5) 
!     C     +tmp3)*0.73d0,-1)  
!      write(*,*) 0.02*j,tmp6,(tmp3-tmp)/tmp5, (tmp3-tmp)/tmp6    
       end do
!       close(55) 
  23  continue
!      write(*,*) alpha_sys, tmp2   

      end do
       
      tmp3 = 100. 
      do i=1,-10
       tmp = 300.
       tmp2 = tmp - 3.*sqrt(tmp) + 6.*sqrt(tmp)/10.*i   
       write(*,*) tmp2,likeli(tmp3,tmp,tmp2),likeli2(tmp3,tmp,tmp2)      
      end do 

!      open(55,file='Sensi_12_GG.dat') 
      do i=1,-100
      tmp2 = 0.05*10**(0.02*i)  
      tmp = (Event_tot(1)*0.022 + Event_tot(2)*0.978)/tmp2**2     ! IH  
      tmp3 = (Event_tot(1)*0.681 + Event_tot(2)*0.319)/tmp2**2    ! NH   
      
       write(*,'(10G14.5)') tmp2,tmp*0.73,tmp3*0.73,   
     C  tmp*0.73+S_b_sys((Event_tot(3)+Event_tot(5)+tmp)*0.73d0,1),       
     C  tmp3*0.73+S_b_sys((Event_tot(3)+Event_tot(5)+tmp3)*0.73d0,-1),
     C  Event_tot(3)*0.73,S_b_sys(Event_tot(3)*0.73d0,1)     
       end do    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  advance warning (2017/11/16)       !!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      write(*,*) 10**(-time(n)),10.**(-time(1))  

      xdl = 0.5   
      xT1 = 1.5
      xT2 = 3. 
!      call AD_Warning_E(xdl,xT1,xT2)        

!      write(*,*) xdl 
      
      do j=1,-10 
      xdl = 0.5 + 0.05*j      ! 0.55 - 1.0     
      do i=1,60 !60     ! 0.1-20 days      
      tmp = 0.1*10.**(0.04*i)    ! 0.1 - 10.  
      call AD_Warning(xdl,tmp)                     
      end do
      end do  

      do i=1,-100    ! likeli_bk0 is quite good   
      tmp = 10.*i
      write(*,*) tmp,likeli_bk(tmp,300.d0),likeli_bk0(tmp, 300.d0),  
     C    likeli_bk00(tmp,300.d0)                
      end do  

      end  

      
      subroutine AD_Warning(xdl,xT)  ! scan T, fixed 4 MeV as optimal E-window 
      use global
      implicit none
      real*8:: xdl,xT,xTup,xEvent(1:7),tmp1,muNH,muIH,xb,xT0,tmp2,
     C   tmp3,tmp4,tmp5    
      real*8:: Tpkk(4)
      integer:: i,j,k,opt   ! NH or IH  
      real*8,external:: likeli_bk0,likeli_bk,likeli_bk00 
!      Tpkk(1:4) = (/1.d0,1.d0,1.d0,1.d0/)     
!      xT0 = max( Tpkk(StelN), xT)     
      tmp3 = -0.1
      tmp4 = -0.1 
      xT0 = max(xT, 0.2)  
      do i = 1, 30      ! T-Scan, how to    
      xTup = xT0*10.**(0.05*i)      ! 3-6 MeV 
      xTup = min(xTup, 10.**(-time(1)))  
      call CalcEvent(xT,xTup,1.804d0,4.d0,xEvent)  
      muNH = (xEvent(1)*0.681 + xEvent(2)*0.319)/xdl**2*0.73       
      muIH = (xEvent(1)*0.022 + xEvent(2)*0.978)/xdl**2*0.73  
      xb =   (xEvent(3) + xEvent(5))*0.73      ! reactor + geo  
!      xb =  (xEvent(4) + xEvent(5))*0.73      !low reactor mode   
      tmp1 = likeli_bk0(muNH,xb)  
      tmp2 = likeli_bk0(muIH,xb)        
      if(tmp3 < tmp1) then
          tmp3 = tmp1
          tmp4 = tmp2 
          tmp5 = xTup
      end if    
      if( tmp3 > 1.1*tmp1)  goto 13          
!      write(*,'(10G14.5)') i,xT,xTup,muNH,muIH,xb,
!     C    tmp1,tmp2    
!!!!! simply assume background can be well measured           
      end do
13      write(*,'(10G14.5)') StelN,i,xdl,xT,tmp5,tmp3,tmp4              
      end  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! the following 3 functions are for advance warning 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      function likeli_bk(xmu,xb)   !  likelihood over background, systematic error for xmu   
      use global                      
      implicit none 
      real*8:: xmu,xb,xn,likeli_bk,Lmu     
      real*8,external:: fbb,fb2b 
      likeli_bk = 0. 
      xn = xmu + xb 
      Lmu = xn*log(xn) - xn  ! nonsense constant, should somehow add it back
      likeli_bk = fb2b(xb,xmu,xn) + Lmu - xn*log(xb) + xb   ! factorial(n) is cancelled off
      likeli_bk = sqrt(2.*likeli_bk)     ! significance       
      return     
      end       
 
      function likeli_bk0(xmu,xb)  ! no systematic error for mu & xb
      use global
      implicit none 
      real*8:: likeli_bk0,xmu,xb,xn 
      xn = xmu + xb 
      likeli_bk0 = xn*log(xn) - xn - (xn*log(xb)-xb)  
      likeli_bk0 = sqrt(2.*likeli_bk0) 
      return   
      end  

      function likeli_bk00(xmu,xb)  ! systematic error for xb instead 
      use global
      implicit none 
      real*8:: likeli_bk00,xmu,xb,xn 
      real*8,external:: FB_bfluc,FB_bfluc2  ! used for xn is large && xmu/xb >> 1 or <<1     
      likeli_bk00 = 0.  
      xn = xb + xmu  
      likeli_bk00 = FB_bfluc(xb,xmu,xn) - FB_bfluc(xb,0.d0,xn) 
      likeli_bk00 = sqrt(2.*likeli_bk00)  
      return
      end  


!!!!!!!! let's just ignore the systematic error of background 
!!!!!!!! hard to calculate         
       function FB_bfluc(xb,xmu,xn)   ! int db pb b**n e**-b, poisson  
       use global          ! given xn events, xb Vs xn
       implicit none       ! tune xb & xmuIH together as background for NH telling 
       real*8:: xb,xn,fb_bfluc,fb,xmu                
       integer, parameter:: Nnp = 200 
       real*8:: x1,x2,x(Nnp),w1(Nnp),tmp = 0.,sigmab, Lmu  
       integer:: i         
       FB_bfluc = 0.
       fb = 0. 
!       xn = xb + xmu    
       sigmab = alpha_sysb  !*xb    !  let's say 10% uncertainty 
       x1 = max(1. - 5.*sigmab, 0.0001d0)      
       x2 = 1. + 5.*sigmab 
       call gauleg(x1,x2,x,w1,Nnp)    
       do i=1,Nnp
       tmp = 1./(sqrt(2.*pi)*sigmab*xb)
     C       *exp(-(x(i)-1.)**2/(2.*sigmab**2))   ! gaussian distri
       FB = FB + tmp*w1(i)*(x(i)+xmu/xb)**xn*exp(-x(i)*xb-xmu)    ! xb**xn * xb  
       end do 
       fb = log(fb) + (xn+1.)*log(xb)   !  = log(L0)   
       Lmu = xn*log(xn) - xn 
       FB_bfluc = FB - Lmu         ! n! has been factored out 
       return  
       end 
 

       function FB_bfluc2(xb,xmu,xn)   !  gaussian, b is large s     
       use global    
       implicit none 
       real*8:: xb,xn,fb2,fb_bfluc2,xmu   
       integer, parameter:: Nnp = 200, NNp2 = 100.  ! should be large somehow when xb is large 
       real*8:: x1,x2,x(Nnp),w1(Nnp),tmp = 0.,sigmab,Lmu 
       real*8:: xn1,xn2,Wn(NNp2),xnn(NNp2),tmp1,tmp2   
       integer:: i,j 
       fb_bfluc2 = 0. 
       FB2 = 0.
       sigmab = alpha_sysb  !*xb    ! let's say 10% uncertainty 
       x1 = max(1. - 5.*sigmab, 0.0001)  
       x2 = 1. + 5.*sigmab 
       xn1 = xn - 0.5     
       xn2 = xn + 0.5 
       call gauleg(xn1,xn2,xnn,Wn,NNp2)    
       call gauleg(x1,x2,x,w1,Nnp)
       do j=1,Nnp2
       do i=1,Nnp
       tmp2 = x(i)*xmu+xb 
       tmp1 = xnn(j)-tmp2-(xnn(j)+0.5)*log(xnn(j)/tmp2)         
       tmp = 1./(sqrt(2.*pi)*sigmab*xb) 
     C       *exp(-(x(i)-1.)**2/(2.*sigmab**2))    ! gaussian distri
     C       *1./sqrt( 2.*pi*(xb*x(i)+xmu) )   
!     C       *exp( -(xnn(j)-x(i)*xb-xmu)**2/( 2.*x(i)*xb+2.*xmu )  )
     C       *exp(tmp1)           ! modified Gaussian 
       FB2 = FB2 + tmp*w1(i)*wn(j)*xb    
       end do        
       end do 
       Lmu = -0.5*log(2.*pi*xn)  
       FB_bfluc2 = log(FB2) - Lmu  
       return  
       end  


      subroutine AD_Warning_E(xdl,xT1,xT2)   ! scan E & see see 
      use global
      implicit none 
      real*8:: xT1,xT2,xdl,xEvent(1:7),xEmax,muNH,muIH,xb,
     C        tmp1,tmp2,tmp3    
      integer:: j 
      real*8,external:: likeli_bk,likeli_bk0,likeli_bk00   
      do j=1,30
      xEmax = 2.5 + 0.1*j     ! 2.5-5.5 MeV 
      call CalcEvent(xT1,xT2,1.804d0,xEmax,xEvent)  
      muNH = (xEvent(1)*0.681 + xEvent(2)*0.319)/xdl**2*0.73         
      muIH = (xEvent(1)*0.022 + xEvent(2)*0.978)/xdl**2*0.73   
      xb =   ( xEvent(3)+xEvent(5) )*0.73      ! reactor + geo  
!      tmp1 = likeli_bk(muNH, xb)
      tmp2 = likeli_bk0(muNH, xb)    ! no error for back & mu 
!      tmp3 = likeli_bk00(muNH,xb)   
      write(*,'(10G14.5)') xEmax,muNH,xb,tmp2       
      end do 
      end     

      subroutine AD_Warning_old(xdl,xT)   ! advanced warning time, fixed distance(in kpc)   
      use global     
      implicit none 
      real*8:: xT, xdl,xEvent(1:7),xEmax,ratio1,ratio2, 
     C    tmp = 0., tmp2 = 0.,tmp3 = 0., ratio1_IH,ratio2_IH , 
     C    tmp2_IH = 0., tmp3_IH = 0. 
      integer:: i,j,i1,i2  
      real*8,external:: S_b,S_b2        
      
      tmp3 = 0.
      if( xT < 0.5*Tpk(StelN) ) then      ! 0.225     
        i1 = 1     
        i2 = 1
      else 
        i1 = 1  
        i2 = 40
      end if    

      do i=i1, i2  ! time scan   
!      tmp = xT*10.**(0.05*i) 
       tmp = max(Tpk(StelN)/1.122, xT)*10**(0.05*i)      

      tmp2 = 0. 
      do j=1, 30   ! Emax scan  
      xEmax = 2.5 + 0.1*j     ! 3-6 MeV  
      
      call CalcEvent(xT,tmp,1.804d0,xEmax,xEvent)      ! maximize E_Range   
      xEvent(1:2) = xEvent(1:2)/xdl**2                   
      ratio1 = (xEvent(1)*0.681 + xEvent(2)*0.319)/xEvent(6)   ! nu Osc, NH
      ratio1_IH = (xEvent(1)*0.022 + xEvent(2)*0.978)/xEvent(6)  ! IH 
!      write(*,*) j,xEvent(1),xEvent(6),ratio1,tmp2  
      if(ratio1 < tmp2) then       
!         CYCLE 
          goto 33
      else    
         tmp2 = ratio1  
      end if   
!      write(*,*) 'haa',i,xEmax, ratio1
!     C   ,xEvent(1),S_b(xEvent(3)+xEvent(5))    
      end do      
33      ratio2 = tmp2 
!        write(*,*) i,ratio2  
      if(ratio2 < tmp3) then
        ratio2 = tmp3      
        goto 34 
      else 
        tmp3 = ratio2   
      end if
!      write(*,*) i,tmp,xEmax,ratio2   
      end do

34      write(*,*) xdl, xT, tmp, ratio2     

      end      
      

      subroutine CalcEvent(xTlow0,xTup0,xEmin0,xEmax0,xEvent)  ! cacl total events at given T&E window 
      use global     
      integer:: i,j 
      real*8::ttx,tty,xEvent(1:7),xTlow0,xTup0,xTlow,xTup,xEmin,xEmax,
     C   xEmin0,xEmax0,xEnu     
      real*8,external:: Crox,Event_react,Event_react_2013,Event_geo,
     C     S_b,S_b2,Event_ron,Event_roff,Event_GeoN     
!!!!  xEvent (1-5: eb, xb, react, react_2013, geo)
!!!   xEvent(6,7):  3 sigma significance (1:2) React+geo,React_2013+geo       
!!!!  integrated events within a certain T & Enu range       
!!!!  normalizaed to 20 kt (12% mass fraction), 1 kpc       
      xEmax = min(xEmax0, Enu(m))    ! set energy window  
      xEmin = max(1.804, xEmin0)  
      xTup = -log10( min( xTup0, 10**(-time(1)) ) )  
      xTlow = -log10( max( xTlow0, 10**(-time(n)) ) )    ! log scale
      xEvent = 0. 
      Event = 0.
      do i=1,1001   ! Time 
      ttx = (xTlow-xTup)/1000.*(i-1) + xTup 
      do j= 1,1001    ! Enu 
      tty = xEmin*10**( log10(xEmax/xEmin)/1000.*(j-1) )            
       call sp2rat(ttx,tty,vint(1),derx,dery,a1,pp,Time,Enu,dx,dy,  
     &                 ix1,iy1,n,m,il,jl,ierr) 
       call sp2rat(ttx,tty,vint(2),derx,dery,a2,pp,Time,Enu,dx,dy,
     &                 ix2,iy2,n,m,il,jl,ierr)   
      Event(1:2,i) = Event(1:2,i)+croX(tty)
     C     *(10**(log10(xEmax/xEmin)/1000.)-1.) 
     C        *tty*0.0952d-42*Npyr*10**(vint(1:2))     
      end do  
      xEvent(1:2) = xEvent(1:2) + Event(1:2,i) 
     C     *10**(-ttx)*(1.-10**(-(xTlow-xTup)/1000.))*24.*3600.   ! < 0
     C     *20./25.    ! 20 kton & 1 kpc  
!      write(*,*) 'hah',i,xevent(1:2) 
      end do 
!      write(*,*) 'hah',xevent(1:2)
!      write(*,*) 'Time window(days):',xTlow0,xTup0  
!      write(*,*) 'Enu(MeV):',xEmin,xEmax     

      do i=1,1000
      xEnu = (xEmax-xEmin)/1000.*i + xEmin 
      xEvent(3)=xEvent(3)+Event_ron(xEnu)*(xEmax-xEmin)/1000.  
     C   *(npyr/1.d32)*(xTup0-xTlow0)/365.*20. 
      xEvent(4)=xEvent(4)+Event_roff(xEnu)*(xEmax-xEmin)/1000.
     C   *(npyr/1.d32)*(xTup0-xTlow0)/365.*20.  
      xEvent(5)=xEvent(5)+Event_geoN(xEnu)*(xEmax-xEmin)/1000. 
     C   *(npyr/1.d32)*(xTup0-xTlow0)/365.*20.     
      end do 

      xEvent(6) = xEvent(3) + xEvent(5)      ! what's wrong here       
      xEvent(7) = xEvent(4) + xEvent(5)   

!      xEvent(6) = S_b2(xEvent(3) + xEvent(5),1)      what's wrong here       
!      xEvent(7) = S_b2(xEvent(4) + xEvent(5),1)     



!      xEvent(8) = S_b(xEvent(3) + xEvent(5) +        )      ! NH as background 
!      xEvent(9) = S_b(xEvent(3) + xEvent(5) +        )      ! NH as background 


!      write(*,*) 'total event(eb,xb):', xEvent(1:2)  
!      write(*,*) 'background', xEvent(3:5) 
!      write(*,*) '3 sigma', xEvent(6:7)  
!       write(*,*) xEmin0,xEvent(1)/xEvent(6),xEvent(1)/xEvent(7)  
 
      end 

      
      subroutine CalcRate(xT,xEmin0,xEmax0,xEvent) ! cacl Rate as function of time  
      use global     
      integer:: i,j 
      real*8::ttx,tty,xEvent(1:2),xT,xEmin,xEmax,
     C   xEmin0,xEmax0,xEnu     
      real*8,external:: Crox,Event_react,Event_react_2013,Event_geo,
     C     S_b,S_b2,Event_ron,Event_roff,Event_GeoN      
!!!!  xEvent (1-5: eb, xb, react, react_2013, geo)
!!!   xEvent(6,7):  3 sigma significance (1:2) React+geo,React_2013+geo      
!!!!  integrated events within a certain T & Enu range       
!!!!  normalizaed to 20 kt (12% mass fraction), 1 kpc       
      xEmax = min(xEmax0, Enu(m))    ! set energy window  
      xEmin = max(1.804, xEmin0)  
      xEvent = 0.
      ttx = -log10(xT) 
      do j= 1,1001    ! Enu 
      tty = xEmin*10**( log10(xEmax/xEmin)/1000.*(j-1) )            
       call sp2rat(ttx,tty,vint(1),derx,dery,a1,pp,Time,Enu,dx,dy,  
     &                 ix1,iy1,n,m,il,jl,ierr) 
       call sp2rat(ttx,tty,vint(2),derx,dery,a2,pp,Time,Enu,dx,dy,
     &                 ix2,iy2,n,m,il,jl,ierr)   
      xEvent(1:2) = xEvent(1:2)+croX(tty)
     C     *(10**(log10(xEmax/xEmin)/1000.)-1.) 
     C        *tty*0.0952d-42*Npyr*10**(vint(1:2))     
      end do  
      xEvent(1:2) = xEvent(1:2) 
     C     *24.*3600.   ! < 0
     C     *20./25.    ! 20 kton & 1 kpc  
     
      end 


!!!!!!!   likelihood ratio,  L(n, mu1, b), b + mu same systematics        
      function likeli(xb0,xmu1,xmu2)   ! xn = xmu1 + xb as expected   
      use global                        
      implicit none 
      real*8,external:: fb,fb2
      real*8:: xb0, xb,xn, xmu1, xmu2, xlambda1, xlambda2, likeli    
!       xb = xb0 + xmu1      ! as effective background 
       xn = xb0 + xmu1       ! expeted events           
       if(xn>500.) then           
       xlambda1 = fb2(xb0+xmu1,xn)  !xn=xb+(xmu2-xmu1)=(xb0+xmu1)+(xmu2-xmu1)=xb0+xmu2
       xlambda2 = fb2(xb0+xmu2,xn)  
       else
       xlambda1 = fb(xb0+xmu1,xn)   
       xlambda2 = fb(xb0+xmu2,xn)        
       end if
       likeli = sqrt(-2.*(xlambda2-xlambda1))  
      return
      end    

!!!!!!!   likelihood ratio,  L(n, mu1, b), keep b constant     
      function likeli2(xb0,xmu1,xmu2)   ! xn = xmu1 + xb as expected  
!!!!!!!!  good option to proceed            
      use global
      implicit none 
      real*8,external:: fbb,fb2b
      real*8:: xb0, xb,xn, xmu1, xmu2, xlambda1, xlambda2, likeli2    
!       xb = xb0 + xmu1      ! as effective background 
       xn = xb0 + xmu1       ! expeted events           
       if(xn>100.) then            
       xlambda1 = fb2b(xb0,xmu1,xn)  !xn=xb+(xmu2-xmu1)=(xb0+xmu1)+(xmu2-xmu1)=xb0+xmu2
       xlambda2 = fb2b(xb0,xmu2,xn)  
       else 
       xlambda1 = fbb(xb0,xmu1,xn) 
       xlambda2 = fbb(xb0,xmu2,xn)         
       end if
       likeli2 = sqrt(-2.*(xlambda2-xlambda1))   
      return
      end 
       
      function likeli3(xb0,xmu1,xmu2)   ! xn = xmu1 + xb as expected   
      use global
      implicit none 
      real*8,external:: fbb,fb2b
      real*8:: xb0, xb,xn, xmu1, xmu2, xlambda1, xlambda2, likeli3    
!       xb = xb0 + xmu1      ! as effective background 
       xn = xb0 + xmu1       ! expeted events           
!       if(xn>500.) then           
       xlambda1 = fb2b(xb0,xmu1,xn)  !xn=xb+(xmu2-xmu1)=(xb0+xmu1)+(xmu2-xmu1)=xb0+xmu2
       xlambda2 = fb2b(xb0,xmu2,xn)   
!       else
!       xlambda1 = fbb(xb0,xmu1,xn) 
!       xlambda2 = fbb(xb0,xmu2,xn)        
!       end if
       likeli3 = sqrt(-2.*(xlambda2-xlambda1))     
      return
      end 

      
!!!! Telling IH from NH, considering the uncertainty of NH rates 
!!!! seems there is no need to do it separately
!!!! b -> b+S_NH, and sensitivity of  ===> S_IH = S_NH - S_b(2,3)       

!!! Equation solve, from b to s 
       function signif(xs, param)   !equation to be solved for s or n for a given background rate 
       use global 
       implicit none   
       real*8::xs,xb,param(0:1),signif
       xb = param(0) 
       signif = sqrt( 2.*( (xs+xb)*log(1.+xs/xb)-xs ) ) - 3.                        
       return
       end

       function S_b(xb,opt)    ! from b->s, no systematic error  
       implicit none
       real*8::xb,S_b,param(0:2),x1,x2,tol
       real*8,external::signif, zbrent  
       integer:: opt
       param(0) = xb
       if(opt.eq.1) then
       x1 = 0.01
       x2 = 1.d3
       else if(opt.eq.-1) then
       x1 = -0.999*xb   !1.d3
       x2 = 0.
       end if 
       tol = 1.d-10              
       S_b = zbrent(signif,param,x1,x2,tol)
       return 
       end 

!!!!! more complex case, considering the uncertaities of background        
!!!!! involve an integration over background distributions (gauss by default) 
!!!!!  when sb is large, using Gaussian Dis rather than Poisson Dis  
       
       function S_b_sys(xb,opt)    ! systematic error 
           use global     ! opt = 1, up fluctuate; -1 for downward fluctuates
       implicit none 
       real*8:: S_b_sys, xb
       real*8,external:: S_b2,S_b3 
       integer:: opt 
       if(xb<500.) then
        S_b_Sys = S_b2(xb,opt)
       else  
        S_b_Sys = S_b3(xb,opt)   
       end if
       return 
       end 


       function signif2(xs, param)      ! (xb,xn) 
       use global      
       implicit none       
       real*8:: xs,xb,xn,param(0:1),signif2,Lmu,xlambda 
       real*8, external:: FB 
       xb = param(0)  
       xn = xs + xb      
!       Lmu = xn*log(xn) - xn                
!       L0 = log(fb(xb,xn)) + (xn+1.)*log(xb)     
       xlambda = fb(xb,xn)   !  L0-Lmu = fb
       signif2 = sqrt(-2.*xlambda) - 3.  
       return
       end 


       function S_b2(xb,opt)    ! from b->s, valid when b<500   
       implicit none
       real*8::xb,S_b2,param(0:2),x1,x2,tol
       real*8,external::signif2, zbrent  
       integer:: opt 
       param(0) = xb 
       if(opt.eq.1) then
       x1 = 0.001
       x2 = 1.d3 
       else if(opt.eq.-1) then      
       x1 = -0.999*xb
       x2 = 0.
       end if 
       tol = 1.d-1              
       S_b2 = zbrent(signif2,param,x1,x2,tol)  
       return 
       end 


       function signif3(xs, param)      ! (xb,xn) 
       use global      
       implicit none       
       real*8:: xs,xb,xn,param(0:1),signif3,Lmu,xlambda 
       real*8, external:: FB2 
       xb = param(0)  
       xn = xs + xb       ! (xs=NH-IH) & xb = xb0 + IH   
!       Lmu = xn*log(xn) - xn                
!       L0 = log(fb(xb,xn)) + (xn+1.)*log(xb)     
       xlambda = fb2(xb,xn)   ! L0-Lmu = fb
       signif3 = sqrt(-2.*xlambda) - 3.   
       return
       end 
 

       function S_b3(xb,opt)    ! from b->s, when b is large, gaussian 
       implicit none
       real*8::xb,S_b3,param(0:2),x1,x2,tol
       real*8,external::signif3, zbrent 
       integer:: opt 
       param(0) = xb 
       if(opt.eq.1) then
       x1 = 0.1
       x2 = 1.d3     ! should be relative small  
       else if(opt.eq.-1) then     
       x1 = -0.999*xb
       x2 = 0.
       end if
       tol = 1.d-1                
       S_b3 = zbrent(signif3,param,x1,x2,tol)   
       return 
       end 

       function FB(xb,xn)  ! int db pb b**n e**-b, poisson  
       use global          ! given xn events, xb Vs xn
       implicit none       ! tune xb & xmuIH together as background for NH telling 
       real*8:: xb,xn,fb              
       integer, parameter:: Nnp = 200 
       real*8:: x1,x2,x(Nnp),w1(Nnp),tmp = 0.,sigmab, Lmu  
       integer:: i         
       FB = 0.
       sigmab = alpha_sys  !*xb    !  let's say 10% uncertainty 
       x1 = max(1. - 5.*sigmab, 0.0001)   
       x2 = 1. + 5.*sigmab 
       call gauleg(x1,x2,x,w1,Nnp)    
       do i=1,Nnp
       tmp = 1./(sqrt(2.*pi)*sigmab*xb)
     C       *exp(-(x(i)-1.)**2/(2.*sigmab**2))   ! gaussian distri
       FB = FB + tmp*w1(i)*x(i)**xn*exp(-x(i)*xb)    ! xb**xn * xb  
       end do 
       fb = log(fb) + (xn+1.)*log(xb)   !  = log(L0)   
       Lmu = xn*log(xn) - xn 
       FB = FB - Lmu         ! n! has been factored out 

       return 
       end 


!!!!  keep background constant        
       function FBB(xb0,xmu,xn)     ! int db pb b**n e**-b, poisson  
       use global             ! given xn events, xb Vs xn
       implicit none 
       real*8:: xb0,xmu,xn,fbb 
       integer, parameter:: Nnp = 200 
       real*8:: x1,x2,x(Nnp),w1(Nnp),tmp = 0.,sigmab, Lmu 
       integer:: i         
       FBB = 0.
       sigmab = alpha_sys  !*xb    !  let's say 10% uncertainty 
       x1 = max(1. - 5.*sigmab, 0.0001)   
       x2 = 1. + 5.*sigmab 
       call gauleg(x1,x2,x,w1,Nnp)     
       do i=1,Nnp
       tmp = 1./(sqrt(2.*pi)*sigmab)  
     C       *exp(-(x(i)-1.)**2/(2.*sigmab**2))   ! gaussian distri
       FBB=FBB+tmp*w1(i)*((x(i)*xmu+xb0)/xmu)**xn*exp(-x(i)*xmu - xb0)    ! xb**xn * xb  
       end do 
       fbb = log(fbb) + xn*log(xmu)   !  = log(L0)   
       Lmu = xn*log(xn) - xn 
       FBB = FBB - Lmu         ! n! has been factored out   
       return 
       end 

       
       function FB2(xb,xn)   !  gaussian, b is large s     
       use global    
       implicit none 
       real*8:: xb,xn,fb2 
       integer, parameter:: Nnp = 200, NNp2 = 100.  ! should be large somehow when xb is large 
       real*8:: x1,x2,x(Nnp),w1(Nnp),tmp = 0.,sigmab,Lmu 
       real*8:: xn1,xn2,Wn(NNp2),xnn(NNp2) 
       integer:: i,j          
       FB2 = 0.
       sigmab = alpha_sys  !*xb    ! let's say 10% uncertainty 
       x1 = max(1. - 5.*sigmab, 0.0001)  
       x2 = 1. + 5.*sigmab 
       xn1 = xn - 0.5     
       xn2 = xn + 0.5 
       call gauleg(xn1,xn2,xnn,Wn,NNp2)   
       call gauleg(x1,x2,x,w1,Nnp)
       do j=1,Nnp2
       do i=1,Nnp
       tmp = 1./(sqrt(2.*pi)*sigmab*xb) 
     C       *exp(-(x(i)-1.)**2/(2.*sigmab**2))    ! gaussian distri
     C       *1./sqrt(2.*pi*xb*x(i)) 
     C       *exp( -(xnn(j)-x(i)*xb)**2/( 2.*x(i)*xb )  )
       FB2 = FB2 + tmp*w1(i)*wn(j)*xb    
       end do        
       end do 
       Lmu = -0.5*log(2.*pi*xn)  
       FB2 = log(FB2) - Lmu  
       return  
       end  


       function FB2B(xb0,xmu,xn)     ! gaussian, b is large s     
       use global    
       implicit none  
       real*8:: xb0,xmu,xn,fb2b  
       integer, parameter:: Nnp = 200, NNp2 = 100.  ! should be large somehow when xb is large 
       real*8:: x1,x2,x(Nnp),w1(Nnp),tmp = 0.,sigmab,Lmu,tmp1,tmp2 
       real*8:: xn1,xn2,Wn(NNp2),xnn(NNp2) 
       integer:: i,j          
       FB2B = 0.
       sigmab = alpha_sys  !*xb    ! let's say 10% uncertainty 
       x1 = max(1. - 5.*sigmab, 0.0001)   
       x2 = 1. + 5.*sigmab 
       xn1 = xn - 0.5     
       xn2 = xn + 0.5 
       call gauleg(xn1,xn2,xnn,Wn,NNp2)   
       call gauleg(x1,x2,x,w1,Nnp)
       do j=1,Nnp2
       do i=1,Nnp
       tmp2 = x(i)*xmu+xb0 
       tmp1 = xnn(j)-tmp2-(xnn(j)+0.5)*log(xnn(j)/tmp2)  
       tmp = 1./(sqrt(2.*pi)*sigmab)  
     C  *exp(-(x(i)-1.)**2/(2.*sigmab**2))    ! gaussian distri
     C  *1./sqrt(2.*pi*(xmu*x(i)+xb0))      ! modify 
!     C       *exp( -(xnn(j)-x(i)*xmu-xb0)**2/( 2.*x(i)*xmu+2.*xb0) )
     C  *exp(tmp1)                   ! modified Gaussian distri 
       FB2B = FB2B + tmp*w1(i)*wn(j)           ! poi -> gauss 
       end do        
       end do 
       Lmu = -0.5*log(2.*pi*xn)        ! L0 (with n=lambda_0) 
       FB2B = log(FB2B) - Lmu  
       return  
       end            ! valid in large limit of n   

!!!!!!!!!!!!! p-value calculation !!!!!!!!!!!! 
      function pval(xN,xb,xmu)
      use global    
      implicit none     
      real*8:: pval,xN,xb,xmu,tmp1,tmp2,tmp3,mut,tmp4
      integer:: i,j
      integer,parameter :: np=30
      real*8:: xNa(np),xwa(np),xalpha(np),walpha(np)
!      real*8:: sigalpha=0.1
      real*8,external:: gammln 
!! poisson dis for xb, Gauss for xmu with 10% uncertainty 

      tmp2 = 0.
      call gauleg(-6.d0*sigalpha,6.d0*sigalpha,xalpha,walpha,np)  
      do i=1,np
!       mut = xb+(1.+xalpha(i))*xmu
       mut = (1.+xalpha(i))*(xb+xmu)        ! vary back+mu with same pace   
       if( xN < mut+6.d0*sqrt(mut) ) then   !then contributes  
       tmp4=exp(-xalpha(i)**2*0.5/sigalpha**2)
     C     /sqrt(2.*3.14159*sigalpha**2)  !*sigalpha   
       tmp1 = max( mut+6.d0*sqrt(mut), mut+20.) 
       tmp1 = max(tmp1, xN+10.d0) 
       call gauleg(xN,tmp1,xNa,xwa,np)
       do j=1,np
       tmp3 = xna(j)*log(mut)-gammln(xna(j))-mut    ! when n is large 
       tmp3 = exp(tmp3)/xna(j)  
       tmp2=tmp2+walpha(i)*xwa(j) 
     C *tmp4       ! gaussian for alpha   
     C *tmp3        !mut**xna(j)*exp(-mut)/factorial(xna(j))   
       
       end do
       end if 
      end do 
      pval = tmp2  
      return
      end 


!!!! solve xN for pval = 1.35d-3 for 3-sigma C.L. 
      function xN_3sig_up(xb,xmu)
      use global    
      implicit none 
      real*8:: xN_3sig_up,xb,xmu,xNup,xNlow,xN, tmp = 1.35d-3,tmp1 
      integer:: i
      real*8,external:: pval 
!!! Newtonian method 
!!! guess solution: 3.*sqrt( (xb+xmu)+(0.1*xmu)**2 )+xb+xmu   
      xNup = 6.*sqrt(xb+xmu+(sigalpha*xmu)**2)+xb+xmu+20.  
      xNlow = xb+xmu       
      do i=1,100
      xN = 0.5*(xNup+xNlow) 
      tmp1 = pval(xN,xb,xmu)
      if( abs(tmp1-tmp)<1.d-7) goto 12   
      if( tmp1 < tmp) then
          xNup = xN
      else 
          xNlow = xN
      end if    
      end do
12      xN_3sig_up = xN  
      return 
      end 



!!!!!!!!!!!!! p-value calculation !!!!!!!!!!!! 
      function pval_low(xN,xb,xmu)
      use global    
      implicit none     
      real*8:: pval_low,xN,xb,xmu,tmp1,tmp2,tmp3,mut,tmp4 
      integer:: i,j
      integer,parameter :: np=30
      real*8:: xNa(np),xwa(np),xalpha(np),walpha(np)
!      real*8:: sigalpha=0.1
      real*8,external:: gammln 
!! poisson dis for xb, Gauss for xmu with 10% uncertainty 

      tmp2 = 0.
      call gauleg(-6.d0*sigalpha,6.d0*sigalpha,xalpha,walpha,np)  
      do i=1,np
       mut = xb+(1+xalpha(i))*xmu
       
       if( xN > mut-6.d0*sqrt(mut) ) then   ! then contributes  
       tmp4=exp(-xalpha(i)**2*0.5/sigalpha**2)
     C     /sqrt(2.*3.14159*sigalpha**2)  !*sigalpha   
       tmp1 = max( mut-6.d0*sqrt(mut), 0.)  
       call gauleg(tmp1,xN,xNa,xwa,np)
       do j=1,np
       tmp3 = xna(j)*log(mut)-gammln(xna(j))-mut    ! when n is large 
       tmp3 = exp(tmp3)/xna(j)  
       tmp2=tmp2+walpha(i)*xwa(j) 
     C *tmp4      ! gaussian for alpha   
     C *tmp3      !mut**xna(j)*exp(-mut)/factorial(xna(j))   
       
       end do
       end if 
      end do 
      pval_low = tmp2  
      return
      end 


      function xN_3sig_low(xb,xmu)
      use global    
      implicit none 
      real*8:: xN_3sig_low,xb,xmu,xNup,xNlow,xN,tmp = 1.35d-3,tmp1 
      integer:: i
      real*8,external:: pval_low 
!!! Newtonian method 
!!! guess solution: 3.*sqrt( (xb+xmu)+(0.1*xmu)**2 )+xb+xmu   
      xNup = xb+xmu  
      xNlow = max(0.d0,-6.*sqrt(xb+xmu+(sigalpha*xmu)**2)+xb+xmu) 
      do i=1,100
      xN = 0.5*(xNup+xNlow) 
      tmp1 = pval_low(xN,xb,xmu) 
      if( abs(tmp1-tmp)<1.d-7) goto 12  
      if( tmp1 < tmp) then
          xNlow = xN 
      else 
          xNup = xN 
      end if    
      end do
12      xN_3sig_low = xN  
      return 
      end 



      FUNCTION gammln(xx)   ! gamma function 
      implicit none     
      REAL*8 gammln,xx
!      Returns the value ln[Γ(xx)] for xx > 0. 
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
! Internal arithmetic will be done in double precision, a nicety 
! that you can omit if five-figure accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
       x=xx
       y=x
       tmp=x+5.5d0
       tmp=(x+0.5d0)*log(tmp)-tmp
       ser=1.000000000190015d0
       do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
       end do 
       gammln=tmp+log(stp*ser/x) 
       return
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  MC simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      function p_MC(xmuNH,xmuIH,xb,np)          
      use global    
      implicit none 
      integer,external:: zbqlpoi 
      real*8:: xmuNH,xmuIH,xb,p_MC
      integer:: nIH,nNH,nb1,nb2,i,ii,np
!      call zbqlini(0)
      ii = 0
!      np = 1000000 
      do i=1,np
      nb1 = zbqlpoi(xb)
      nIH = zbqlpoi(xmuIH) 
      nb2 = zbqlpoi(xb)
      nNH = zbqlpoi(xmuNH)   
      if(nIH+nb1-nNH-nb2.ge.0) then
          ii = ii+1
      end if    
      end do
      p_MC = 1.*ii/np 
      return 
      end 

!!!! find the upbound of IH       
      function NH_MC(xmuIH,xb,np) 
!      subroutine xNH_MC(xmuIH,xb,np)         
      use global
      implicit none 
      real*8:: xmuIH,xb,xmuNH,NH_MC,xmuNHup,xmuNHlow
      integer:: i,np 
      real*8,external:: p_MC

       xmuNHup = max(30.d0,xmuIH+10.d0*sqrt(xmuIH))   
       xmuNHlow = xmuIH
       xmuNH = xmuNHup 
      do i=1,500
        if( abs(p_MC(xmuNH,xmuIH,xb,np)/1.32d-3 - 1.)<1.d-2 ) goto 12
        if( p_MC(xmuNH,xmuIH,xb,np) > 1.32d-3) then 
          xmuNHlow = xmuNH
        else
          xmuNHup = xmuNH
        end if
       xmuNH = 0.5*(xmuNHup + xmuNHlow)
!       write(*,*) xmuNH, p_MC(xmuNH,xmuIH,xb,np)
      end do
12    NH_MC = xmuNH
!       write(*,*) xmuNH, p_MC(xmuNH,xmuIH,xb,np)
      return
      end 



      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END


      FUNCTION zbrent(func,param,x1,x2,tol)
      INTEGER ITMAX
      REAL*8 zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,param(50) 
      a=x1
      b=x2
      fa=func(a,param)
      fb=func(b,param)
!      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause 
!     c 'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,param)
11    continue
!      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END



      function croX(xE)   ! IBD cross section
      use global     
      implicit none           ! 1.804 = 0.511 + 1.293 
      real*8:: croX, xE    
      if(xE<1.804) then 
          croX = 0.
      else 
          croX = sqrt( (xE-1.293)**2 - me**2 )*(xE-1.293)  
      end if 

      return 
      end 
       

!!!! reactor neutrino and geo-neutrinos for JUNO 

      subroutine Init_back()    
      use global     
      real*8:: yp1,yp2 
      integer:: np
!      153,43,9,6,12 

      np = 153  
      yp1 = (ron(2)-ron(1))/(Enu_ron(2)-Enu_ron(1))
      ypn = (ron(NP)-ron(NP-1))/(Enu_ron(NP)-Enu_ron(NP-1))  
      call spline(enu_ron,ron,153,yp1,ypn,y2_ron)  

      np = 43  
      yp1 = (roff(2)-roff(1))/(Enu_roff(2)-Enu_roff(1))
      ypn = (roff(NP)-roff(NP-1))/(Enu_roff(NP)-Enu_roff(NP-1))  
      call spline(enu_roff,roff,43,yp1,ypn,y2_roff) 

      np = 9  
      yp1 = (geo1(2)-geo1(1))/(Enu_geo1(2)-Enu_geo1(1))
      ypn = (geo1(NP)-geo1(NP-1))/(Enu_geo1(NP)-Enu_geo1(NP-1))  
      call spline(enu_geo1,geo1,9,yp1,ypn,y2_geo1)  

      np = 6  
      yp1 = (geo2(2)-geo2(1))/(Enu_geo2(2)-Enu_geo2(1))
      ypn = (geo2(NP)-geo2(NP-1))/(Enu_geo2(NP)-Enu_geo2(NP-1))  
      call spline(enu_geo2,geo2,6,yp1,ypn,y2_geo2)   

      np = 12   
      yp1 = (geo3(2)-geo3(1))/(Enu_geo3(2)-Enu_geo3(1))
      ypn = (geo3(NP)-geo3(NP-1))/(Enu_geo3(NP)-Enu_geo3(NP-1))  
      call spline(enu_geo3,geo3,9,yp1,ypn,y2_geo3) 

      end 

      function Event_geoN(xxEnu) 
      use global
      real*8:: Event_geoN, xxEnu, y 
!      153,43,9,6,12 
      if(xxEnu < 2.07) then
      call splint(Enu_geo1,geo1,y2_geo1,9,xxEnu,y)  
      else if(xxEnu < 2.26) then
      call splint(Enu_geo2,geo2,y2_geo2,6,xxEnu,y)   
      else if(xxEnu < 3.27) then
      call splint(Enu_geo3,geo3,y2_geo3,12,xxEnu,y) 
      else 
      y = 0.
      end if 
      Event_geoN = y*100.    ! TNU/MeV
      return
      end 

      function Event_Ron(xxEnu)
      use global
      real*8:: xxEnu, Event_Ron, y   
      call splint(Enu_ron,ron,y2_ron,153,xxEnu,y) 
      Event_Ron = y*100.  
      return
      end 

      function Event_Roff(xxEnu)
      use global
      real*8:: xxEnu, Event_Roff, y   
      call splint(Enu_roff,roff,y2_roff,43,xxEnu,y) 
      Event_Roff = y*100.  
      return
      end 


       function Event_geo(xxEnu)
       integer,parameter:: NP = 23 
       real*8:: xxEnu, xEnu(NP), xF(NP), Event_geo     
       real*8:: yp1,ypn,y2(NP),y   
       if(xxEnu > 3.286 .or. xxEnu < 1.804) then
           Event_geo = 0.
       else      
        xEnu(1:23)=(/1.81324,1.82403,1.85615,1.89364,1.94712,2.00596,
     C        2.05945,2.11822,2.17702,2.25716,2.26764,2.27278,2.27776,
     C        2.2779,2.28298,2.40576,2.52321,2.68333,2.88082,3.07829,
     C        3.27038,3.28097,3.28621/)                
        xF(1:23) = (/0.02416,0.11206,0.18796,0.27983,0.36769,0.46753,
     C          0.55539,0.60727,0.67914,0.743,0.70311,0.13855,0.11044,
     C         0.1099,0.09562,0.11947,0.14333,0.15916,0.17095,0.17475,
     C         0.15058,0.14479,0.13993/)
!!! interpolation then 
      yp1 = (xF(2)-xF(1))/(xEnu(2)-xEnu(1))
      ypn = (xF(NP)-xF(NP-1))/(xEnu(NP)-xEnu(NP-1))  
      call spline(xEnu,xF,NP,yp1,ypn,y2) 
      call splint(xEnu,xF,y2,NP,xxEnu,y)  
      Event_geo = y*100.      ! TNU/MeV
      end if  
      return 
      end 


       function Event_react(xxEnu)  ! old 
       implicit none 
       integer,parameter:: NP = 46
       real*8:: xxEnu, xEnu(NP), xF(NP), Event_react    
       real*8:: yp1,ypn,y2(NP),y
       if(xxEnu > 9.0 .or. xxEnu < 1.804) then
           Event_react = 0.
       else   
         xEnu(1:NP)= (/1.80791,1.81348,1.81378,1.82476,
     C           1.84107,1.86811,1.8898,1.91145,1.94369,1.9599,
     C           1.9708,2.00296,2.04589,2.1102,2.19053,2.25482,
     C           2.41458,2.54767,2.69148,2.87271,3.14516,3.43912,
     C           3.70111,4.00582,4.35857,4.77515,5.02586,5.28718,
     C           5.60689,5.89993,6.07595,6.16649,6.27286,6.38987,
     C           6.53898,6.7202,6.89076,7.01327,7.1304,7.30613,
     C           7.47129,7.73782,8.05773,8.31373,8.72987,8.99668/) 
         xF(1:NP) = (/0.02816,0.414,0.43175,0.65552,
     C           0.88328,1.15098,1.4067,1.63046,1.79826,1.9421,
     C           2.11791,2.21778,2.39356,2.59729,2.805,2.98476,
     C           2.72886,2.47698,2.26105,2.09303,2.29656,2.62394,
     C          2.97932,3.36263,3.7459,3.97325,3.88109,3.74497,
     C          3.36503,2.9971,2.925,2.78505,2.50922,2.2094,
     C          1.9655,1.79349,1.62149,1.44554,1.23763,0.94175, 
     C           0.7298,0.49377,0.26568,0.14154,0.03723,0.01299/)    
      yp1 = (xF(2)-xF(1))/(xEnu(2)-xEnu(1))
      ypn = (xF(NP)-xF(NP-1))/(xEnu(NP)-xEnu(NP-1))  
      call spline(xEnu,xF,NP,yp1,ypn,y2) 
      call splint(xEnu,xF,y2,NP,xxEnu,y)
      Event_react = y*100.      ! TNU/MeV   
      end if 
      return
      end 


       function Event_react_2013(xxEnu)    ! React/2013, old  
       implicit none 
       integer,parameter:: NP = 35
       real*8:: xxEnu, xEnu(NP), xF(NP), Event_react_2013    
       real*8:: yp1,ypn,y2(NP),y
       if(xxEnu > 8.0 .or. xxEnu < 1.804) then
           Event_react_2013 = 0.
       else   
        xEnu(1:NP) = (/1.81855,1.89329,1.96268,2.05349,2.15499,
     C      2.26176,2.39513,2.51246,2.6245,2.7847,2.88089,
     C      3.00911,3.14266,3.18539,3.32952,3.53223,3.63354,
     C      3.73489,3.8682,3.99619,4.23097,4.46041,4.73263,
     C      5.02088,5.22903,5.5119,5.79475,6.10427,6.40844,
     C      6.72861,7.20354,7.50238,7.7692,7.99335,8.29754/) 
        xF(1:NP) = (/0.00418,0.02408,0.03201,0.09185,0.16368,
     C      0.19155,0.14746,0.0874,0.06331,0.13108,0.22689,
     C      0.33066,0.42643,0.44637,0.4742,0.40606,0.33803,
     C      0.30196,0.21791,0.14586,0.10966,0.07746,0.10916,
     C      0.14883,0.16061,0.17232,0.16404,0.15174,0.12746,
     C      0.09118,0.05075,0.03047,0.01821,0.01399,0.0057/)  
      yp1 = (xF(2)-xF(1))/(xEnu(2)-xEnu(1))
      ypn = (xF(NP)-xF(NP-1))/(xEnu(NP)-xEnu(NP-1))  
      call spline(xEnu,xF,NP,yp1,ypn,y2) 
      call splint(xEnu,xF,y2,NP,xxEnu,y)
      Event_react_2013 = y*100.      ! TNU/MeV    
      end if 
      return
      end 


!!!! interpolation for 1D 
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit none     
      INTEGER n
      REAL*8 x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h 
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))
     C   *(h**2)/6.
      return
      END

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit none    
      INTEGER:: n 
      REAL*8:: yp1,ypn,x(n),y(n),y2(n)
      integer, PARAMETER:: NMAX=500 
      INTEGER:: i,k
      REAL*8:: p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1) = -0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-
     C    (y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then 
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k) 
12    continue     
      return
      END

!!!! interpolation for 2D 
      include 'sp2rat.f'    
      include 'sp2ratco.f'      
      include 'randgen.f'     ! random number generator:ZBQLINI(SEED),ZBQLNOR(MU,SIGMA),ZBQLPOI(MU)

