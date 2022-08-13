
module sim_parameters
  implicit none
  double precision,parameter :: pai = 3.141592653589793D+00
  double precision,parameter :: dt=0.001    !A step width of dynamics
  
  integer,parameter :: N=1000               !Number of KaiC hexamers
  integer,parameter :: ifft=0               !ifft=1 for the FFT calculation and =0 for others
  
  integer,parameter :: no=ifft*(100000+100*2**15)+(1-ifft)*300000   !Max number of simulation steps

  
  integer,parameter :: itempshift=100                 !Step no. at which temperature is shifted
  integer,parameter :: itempshift2=90000

  integer,parameter :: iinput=0
  integer,parameter :: ioutput=0
  integer,parameter :: ioutstepa=61400  
  integer,parameter :: ioutstepb=71400
 
  double precision,parameter :: cT=30.0 +273.15       !Temperature at the beginning
  double precision,parameter :: cT2=25.0 +273.15      !Temperature later than the itempshift step
  double precision,parameter :: cT3=25.0 +273.15      !Temperature later than the itempshift2 step
  character(len=64),parameter :: fname='u30.dat'
  
  
  double precision,parameter :: sn1=1.0               !scaling factor
  double precision,parameter :: sn2=1.0               !scaling factor
  double precision,parameter :: sn3=1.0               !scaling factor
  
!Structure modulation
  double precision,parameter :: tc0=2.0               !Shifting the average W
  double precision,parameter :: tc1=3.0*sn1           !Negative feedback strength to determine W
  double precision,parameter :: tc2=5.0*sn2           !Positive feedback strength to determine W
  double precision,parameter :: tc2kaiA=5.0*sn2       !Positive feedback strength to determine W  

!KaiAB abiudance  
  double precision,parameter :: arA=1.0               !Modulation coefficient of rates on W  
  double precision,parameter :: arB=1.0               !Modulation coefficient of rates on W
  
  double precision,parameter :: Atot0=1.0*N           !Reference total number of KaiA dimers
  double precision,parameter :: Btot0=6.0*N           !Reference total number of KaiB monomers  
  double precision,parameter :: Atot=1.0*N            !Total number of KaiA dimers
  double precision,parameter :: Btot=6.0*N            !Total number of KaiB monomers 

!****** From here, rates should be scaled in a parallel way when unit of time is changed ********

!Soft spin amplitude   
  double precision,parameter :: tam=0.5                !Amplitude of soft-spin potential
  
!Phoshorylation/dephosphorylation 
  double precision,parameter :: trk1=0.18              !Rate of phosphorylation
  double precision,parameter :: tdk1=0.18              !Rate of dephosphorylation


!Binding/unbinding rates of KaiABC       
  double precision,parameter :: thA0=0.5/(Atot0)       !Rate coefficient of binding KaiA to KaiC
  double precision,parameter :: tfA0=1.0               !Rate of unbinding KaiA from KaiC
  double precision,parameter :: thB0=0.30/(Btot0)      !Rate coefficient of binding KaiB to KaiC
  double precision,parameter :: tfB0=2.0               !Rate of unbinding KaiB from KaiC
  double precision,parameter :: thBA=6/(Atot0)         !Rate coefficient of binding KaiA to KaiCB
  double precision,parameter :: tfBA=1.0               !Rate of unbinding KaiA from KaiCB

  
!********** Uo to here, rates and a duration  ******************************************************
  
!Parameters for ATP hydrolysis 
  double precision,parameter :: taw=0.5               !Modulation coefficient of duration on W    
  real(kind=8):: aw    
  double precision,parameter :: taq=(2.0/6)*sn1       !Strength of perturbation due to the Pi release

  
!ATP hydrolysis rate(frequency) and duration (Notice that duration should be changed inversely!!!)
    double precision,parameter :: tomega=1.0*sn3      !Basic frequency of Pi release
    double precision,parameter :: tdelta0=1.0*sn3     !Duration steps of perturbation due to the Pi release

    integer,parameter :: iatpchange=0
    integer,parameter :: iatpdrop=75000
    integer,parameter :: iatprise=iatpdrop+6000    
    double precision,parameter :: atpratio=4.0




  integer,parameter :: mesh=100             !# of mesh for P1 and P2 to calculate mutual information

  integer,parameter :: imode=2              !imode=1 for soft-spin dynamics
                                            !imode=2 for chemical kinetics dynamics
    
  integer,parameter :: isynch=3  !isynch=1 to calculate the dispersion among N KaiC molecules
                                 !isynch=2 to calculate the mutual information
                                 !isynch=3 to average individual mutual information



   
end module sim_parameters




!****************************************************
program main
  
  use sim_parameters
  implicit none
  
  integer i,j,k,ki,kj,istep,iseed,ireact,ireact2,irdiff
  double precision s,PN,time,Amol,Bmol,sfA,sPha2,Utot,Wtot,areact
  double precision Uave,disp,D2,DD2,umut
  
  double precision r8_uniform_01,r8_normal_01

  integer,dimension(6,N) :: jstep         !Memory of Pi release at ith subunit of kth KaiC hexamer
  double precision,dimension(N) :: U      !Phosphorylation level of kth KaiC hexamer
  double precision,dimension(N) :: W      !Structure state of kth KaiC hexamer
  double precision,dimension(N) :: fw     !Frequency of Pi release modulated by W
  double precision,dimension(N) :: Pha2   !Probability of kth KaiC hexamer to bind KaiA dimer
  double precision,dimension(0:6,N) :: PB
  double precision,dimension(N) :: PBt  
  double precision,dimension(N) :: q      !Structure moduation in kth KaiC due to Pi release
  double precision,dimension(6,N) :: qq   !Structure modulation in ith subunit of kth KaiC
  double precision,dimension(6,N) :: qq2

  double precision c0,c1,c2,c2kaiA,am,aq,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T
  double precision qdiff,qt,ha1,gkBA


 
  double precision,dimension(0:mesh,0:mesh,10,10) :: P2
  double precision,dimension(0:mesh,10) ::P1a
  double precision,dimension(0:mesh,10) ::P1b  
  integer,dimension(10) :: ipair1,ipair2
  
  double precision,dimension(N) :: Ua
  double precision,dimension(N) :: Wa
  double precision,dimension(N) :: Pha2a
  double precision,dimension(N) :: Ub
  double precision,dimension(N) :: Wb
  double precision,dimension(N) :: Pha2b
  double precision,dimension(6,N) :: qqa
  double precision,dimension(6,N) :: qq2a
  double precision,dimension(6,N) :: qqb
  double precision,dimension(6,N) :: qq2b
  double precision,dimension(0:6,N) :: PBa
  double precision,dimension(0:6,N) :: PBb
  integer,dimension(6,N) :: jstepa
  integer,dimension(6,N) :: jstepb
  double precision Amola,Amolb,Bmola,Bmolb,timea,timeb
  integer iseedb,iouta,ioutb
  
  !============= Initial values  =======================
      iseed=2341571
      ireact=0
      ireact2=0

      timea=ioutstepa*dt
      timeb=ioutstepb*dt

      D2=0.0
      DD2=0.0

      if(iinput==0) then
         Amol=0.2*Atot
         Bmol=0.6*Btot
         ha1 = thBA*Amol
         gkBA= ha1/tfBA
      
         do k=1,N
            s=r8_normal_01(iseed)
            U(k)=0.0+s*0.2
            s=r8_normal_01(iseed)
            W(k)=0.0+s*0.2
            do i=1,6
              jstep(i,k)=0
              s=r8_uniform_01(iseed)*0
              if(s.le.0.5) then
                qq(i,k)=0.0
              else
                qq(i,k)=-aq
              endif
              qq2(i,k)=qq(i,k)
            enddo
        
            Pha2(k)=abs(sin(0.01*k*pai+30))
            PBt(k)=0.0
            Do i=0,6
               PB(i,k)=abs(sin(0.02*(k+i)*pai+10))
               PBt(k)=PBt(k)+PB(i,k)
            Enddo
            PN=Pha2(k)+PBt(k)
            Pha2(k)=Pha2(k)/PN
            Do i=0,6
              PB(i,k)=PB(i,k)/PN
            Enddo
         Enddo
      endif

      
      if(iinput==1) then
            open(unit=81, file='r0a.dat')
            read(81,*) iseed,iouta,Amol,Bmol
            do k=1,N
               read(81,*) U(k),W(k),Pha2(k)
               read(81,*) PB(0,k)
               do i=1,6
                  read(81,*) jstepa(i,k),PB(i,k),qq(i,k),qq2(i,k)
                  jstep(i,k)=jstepa(i,k)-iouta
               enddo
            enddo
            ha1 = thBA*Amol
            gkBA= ha1/tfBA
            close(unit=81)
      endif     

      if(iinput==2) then        
            open(unit=81, file='r0b.dat')
            read(81,*) iseed,iouta,Amola,Bmola
            do k=1,N
               read(81,*) Ua(k),Wa(k),Pha2a(k)
               read(81,*) PBa(0,k)
               do i=1,6
                  read(81,*) jstepa(i,k),PBa(i,k),qqa(i,k),qq2a(i,k)
               enddo
            enddo
            close(unit=81)
            
            open(unit=82, file='r0a.dat')
            read(82,*) iseedb,ioutb,Amolb,Bmolb
            do k=1,N
               read(82,*) Ub(k),Wb(k),Pha2b(k)
               read(82,*) PBb(0,k)
               do i=1,6
                  read(82,*) jstepb(i,k),PBb(i,k),qqb(i,k),qq2b(i,k)
               enddo
            enddo
            close(unit=82)
            do k=1,N/2
               U(k)=Ua(k)
               W(k)=Wa(k)
               Pha2(k)=Pha2a(k)
               PB(0,k)=PBa(0,k)
               do i=1,6
                  jstep(i,k)=jstepa(i,k)-iouta
                  PB(i,k)=PBa(i,k)
                  qq(i,k)=qqa(i,k)
                  qq2(i,k)=qq2a(i,k)
               enddo
            enddo
            do k=N/2+1,N
               U(k)=Ub(k)
               W(k)=Wb(k)
               Pha2(k)=Pha2b(k)
               PB(0,k)=PBb(0,k)
               do i=1,6
                  jstep(i,k)=jstepb(i,k)-ioutb
                  PB(i,k)=PBb(i,k)
                  qq(i,k)=qqb(i,k)
                  qq2(i,k)=qq2b(i,k)
               enddo                                   
            enddo
            Amol=(Amola+Amolb)*0.5
            Bmol=(Bmola+Bmolb)*0.5
            ha1 = thBA*Amol
            gkBA= ha1/tfBA
      endif     
         
      do i=1,10
         ipair1(i)=i*60
         ipair2(i)=i*60+25
         do j=0,mesh
            P1a(j,i)=0.0
            P1b(i,i)=0.0
         enddo
      enddo
      do i=1,10
         do j=1,10
            do ki=0,mesh
               do kj=0,mesh
                  P2(ki,kj,i,j)=0.0
               enddo
            enddo
         enddo
      enddo
      

!============= File opening  =======================
      
      open(unit=30, file=fname, status='replace')
      open(unit=31, file='w.dat', status='replace')
      open(unit=33, file='a.dat', status='replace')
      open(unit=34, file='ba.dat', status='replace')
      open(unit=35, file='amol.dat', status='replace')

      open(unit=40, file='u0.dat')
      open(unit=41, file='u1.dat')
      open(unit=42, file='u2.dat')
      open(unit=43, file='u3.dat')
      open(unit=44, file='u4.dat')

      open(unit=49, file='qt.dat')
      open(unit=50, file='qdiff.dat')
      open(unit=51, file='q1.dat')

      open(unit=61, file='w1.dat')
      
      open(unit=71, file='r0a.dat') 
      open(unit=72, file='r0b.dat')
      open(unit=73, file='u0a.dat')
      open(unit=74, file='u0b.dat')

      
      open(unit=39, file='itemps')
      write(39,*) itempshift     
!============= Loop of dynamics =======================

  T=cT/cT
  call Tmodulation(c0,c1,c2,c2kaiA,am,aq,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
  call ATPhyd1(U,W,fw,aq,q,qq,qq2,qt,Pha2,jstep,istep,iseed,ireact,c0,am,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
  call binding(sfA,sPha2,W,Amol,Bmol,Pha2,PB,q,istep,c0,c1,c2,c2kaiA,am,rk1,dk1, &
                      omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
  call diffmodel(Pha2,PB,U,W,q,Utot,Wtot,c0,c1,c2,c2kaiA,am,rk1,dk1,omega,delta0, &
                       hA0,fA0,hB0,fB0,hBA,fBA,T,iseed,aq)  
  time=0
            Uave=Utot/dble(N)
            write(30,*) time,(Uave+1)*0.5
            flush 30
            if(ifft.eq.0) then 
               write(31,*) time,(Wtot/dble(N)+1)*0.5
               write(33,*) time,sPha2/dble(N)
               write(34,*) time,sfA/dble(N)
               write(35,*) time,Amol/Atot           
               flush 31
               flush 33
               flush 34
               flush 35
            
               write(40,*) time,(U(1)+1)*0.5
               write(41,*) time,(U(100)+1)*0.5
               write(42,*) time,(U(200)+1)*0.5
               write(43,*) time,(U(300)+1)*0.5
               write(44,*) time,(U(400)+1)*0.5
               flush 40
               flush 41
               flush 42
               flush 43
               flush 44

            irdiff=ireact-ireact2
            ireact2=ireact
            qdiff=real(irdiff)/(6.0*N)/2*10
               write(49,*) time,-qt
               write(50,*) time,qdiff
               write(51,*) time,-q(100)/(6.0*aq)
               flush 49
               flush 50
               flush 51
             endif      

  
  do istep=1,no
    
         if(istep.eq.itempshift) then
            T=cT2/cT
            call Tmodulation(c0,c1,c2,c2kaiA,am,aq,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
         endif
         if(istep.eq.itempshift2) then
            T=cT3/cT
            call Tmodulation(c0,c1,c2,c2kaiA,am,aq,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
         endif         
         if(iatpchange==1) then
            if(istep==iatpdrop) delta0=atpratio*delta0
            if(istep==iatprise) delta0=delta0/atpratio
         endif      


         
         call ATPhyd1(U,W,fw,aq,q,qq,qq2,qt,Pha2,jstep,istep,iseed,ireact,c0,am,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)  
         call binding(sfA,sPha2,W,Amol,Bmol,Pha2,PB,q,istep,c0,c1,c2,c2kaiA,am,rk1,dk1, &
                      omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
         call diffmodel(Pha2,PB,U,W,q,Utot,Wtot,c0,c1,c2,c2kaiA,am,rk1,dk1,omega,delta0, &
                       hA0,fA0,hB0,fB0,hBA,fBA,T,iseed,aq)        


         areact=dble(ireact)/(6.0*N)/(istep*dt)



         Uave=Utot/dble(N)
         if(istep.gt.100000) then
            
            if(isynch.eq.1) then
               call dispersion(U,Uave,disp)
               DD2=DD2+disp
               D2=DD2/dble(istep-100000)
            endif
            
            if(isynch.ge.2) then
               call contourplot(U,P2,P1a,P1b,ipair1,ipair2)
            endif           
            
         endif 
         
         

         
            
  if(ifft.eq.0) then
     ! if(mod(istep,200).eq.0) then !.and.istep.lt.100000) then
     !        irdiff=ireact-ireact2
    !        ireact2=ireact
    !        qdiff=real(irdiff)/(6.0*N)/2*10
    !   endif 
         
      if(mod(istep,200).eq.0) then !.and.istep.ge.100000) then
            
            time=istep*dt
            write(30,*) time,(Uave+1)*0.5
            flush 30
            if(ioutput==2) then
               write(73,*) time-timea,(Uave+1)*0.5
               write(74,*) time-timeb,(Uave+1)*0.5
            endif
            
               write(31,*) time,(Wtot/dble(N)+1)*0.5
               write(33,*) time,sPha2/dble(N)
               write(34,*) time,sfA/dble(N)
               write(35,*) time,Amol/Atot           
               flush 31
               flush 33
               flush 34
               flush 35
            
               write(40,*) time,(U(1)+1)*0.5
               write(41,*) time,(U(100)+1)*0.5
               write(42,*) time,(U(200)+1)*0.5
               write(43,*) time,(U(300)+1)*0.5
               write(44,*) time,(U(400)+1)*0.5
               flush 40
               flush 41
               flush 42
               flush 43
               flush 44             

            irdiff=ireact-ireact2
            ireact2=ireact
            qdiff=real(irdiff)/(6.0*N)/2*10
            
               write(49,*) time,-qt
               write(50,*) time,qdiff
               write(51,*) time,-q(100)/(6.0*aq)
               flush 49
               flush 50
               flush 51

               write(61,*) time,(W(100)+1)*0.5
               flush 61


            endif
         else
            if(mod(istep,100).eq.0) then
               time=istep*dt
               write(30,*) time,(Uave+1)*0.5
               flush 30           
            endif
         endif

         if(mod(istep,1000).eq.0.and.ifft.eq.0) then
             if(isynch.eq.3) then
                call mutinf2(P2,P1a,P1b,umut,istep)
!               write(6,*) istep,Amol,Bmol,areact,umut
                write(6,*) istep,T,areact,umut
            endif           
            if(isynch.eq.2) then
               call mutinf(P2,P1a,P1b,umut,istep)
               write(6,*) istep,Amol,Bmol,areact,umut
            endif            
            if(isynch.eq.1) write(6,*) istep,Amol,Bmol,areact,D2
            if(isynch.eq.0) write(6,*) istep,Amol,Bmol,areact
         endif

         
         if(istep==ioutstepa.and.ioutput.ge.1) then
            write(71,*) iseed,ioutstepa,Amol,Bmol
            do k=1,N
               write(71,*) U(k),W(k),Pha2(k)
               write(71,*) PB(0,k)
               do i=1,6
                  write(71,*) jstep(i,k),PB(i,k),qq(i,k),qq2(i,k)
               enddo
            enddo
            flush 71
         endif
         if(istep==ioutstepb.and.ioutput==2) then
            write(72,*) iseed,ioutstepb,Amol,Bmol
            do k=1,N
               write(72,*) U(k),W(k),Pha2(k)
               write(72,*) PB(0,k)
               do i=1,6
                  write(72,*) jstep(i,k),PB(i,k),qq(i,k),qq2(i,k)
               enddo
            enddo
            flush 72
         endif        
         
  enddo


stop
end program main


!*****************************************************
subroutine Tmodulation(c0,c1,c2,c2kaiA,am,aq,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
    use sim_parameters
    implicit none
    
    double precision ap,ap2,ap3,ap4,ap5,ap6,ap7,ap8,ap10,ap12,ap15,ap20,gp
    double precision c0,c1,c2,c2kaiA,am,aq,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T
 

  
    ap=exp(1-1/T)
       ap2=ap*ap
       ap3=ap**3
       ap4=ap3*ap
       ap5=ap4*ap
       ap6=ap3*ap3
       ap7=ap3*ap4
       ap8=ap4*ap4
       ap10=ap5*ap5
       ap12=ap7*ap5
       ap15=ap10*ap5
       ap20=ap10*ap10
     
  


!     ap increases as T increases
!     Enhancement of positive feedback or Reduction of negative feedback
!     increases the oscillation period

       

!       gp=12.0            !shifting the average W
!       c0=tc0*(1.0+gp*(1.0-T))       !Shifting the average W       
!       c1=tc1*(1.0+gp*(1.0-T))       !Effects of P/dP on W (Negative f.b.)
!       aq=taq*(1.0+gp*(1.0-T))       !Effects of ATPase on W  

       c0=tc0           !Shifting the average W   
       c1=tc1/ap7       !Effects of P/dP on W (Negative f.b.)
       aq=taq/ap7       !Effects of ATPase on W
       
       c2=tc2          !Effects of KaiB binding on W (Positive f.b.)
       c2kaiA=tc2kaiA    !Effects of KaiA binding on W (Positive f.b.)

         
       am=tam*ap10             !Amplitude of soft-spin potential
      
       rk1=trk1*ap10           !Rate of phosphorylation
       dk1=tdk1*ap10           !Rate of dephosphorylation
  
       omega=tomega*ap10           !Basic frequency of Pi release
       delta0=tdelta0!/ap10         !Duration steps of perturbation due to the Pi release  
       aw=taw!/ap20
       
       hA0=thA0*ap10          !Rate coefficient of binding KaiA to KaiC
       fA0=tfA0*ap10          !Rate of unbinding KaiA from KaiC
       hB0=thB0*ap10*ap10       !Rate coefficient of binding KaiB to KaiC
       fB0=tfB0*ap10*ap10*ap2      !Rate of unbinding KaiB from KaiC
       hBA=thBA*ap10          !Rate coefficient of binding KaiA to KaiCB 
       fBA=tfBA*ap10          !Rate of unbinding KaiA from KaiCB

       
  return
end subroutine Tmodulation








!************************************************************************
subroutine ATPhyd1(U,W,fw,aq,q,qq,qq2,qt,Pha2,jstep,istep,iseed,ireact,c0,am,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
  
  use sim_parameters
  implicit none
  
  integer i,k,istep,iseed,iwait,ireact
  double precision s,r8_uniform_01,wq,r8_normal_01,twq,aq
  
  integer,dimension(6,N) :: jstep
  double precision,dimension(N) :: U
  double precision,dimension(N) :: W  
  double precision,dimension(N) :: fw
  double precision,dimension(N) :: Pha2  
  double precision,dimension(N) :: q  
  double precision,dimension(6,N) :: qq
  double precision,dimension(6,N) :: qq2
  
  double precision c0,am,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T
  double precision qt,qrt,srt,qrtf
  
  twq=delta0/dt
  qt=0.0
  
  do k=1,N

!              qrtf=1/cosh(aw*0.5*(W(k)-1.0))**2
              qrt=1-tanh(aw*W(k))
            
 !             fw(k)=omega*dt*qrtf
              fw(k)=omega*dt             
              srt=(twq*qrt)**0.5
         
              q(k)=0.0
              do i=1,6
                s=r8_uniform_01(iseed)
                if(s.le.fw(k).and.qq(i,k).ge.0.0) then
                   qq(i,k)= -aq
                   jstep(i,k)=istep
!                   ireact=ireact+1
                endif
                s=r8_normal_01(iseed)   
                iwait=int( twq*qrt+srt*s  )

                
                if(istep.gt.jstep(i,k)+iwait) qq(i,k)=0.0
                if(qq2(i,k).lt.0.0.and.qq(i,k).ge.0.0) then
                   ireact=ireact+1
                endif
                qq2(i,k)=qq(i,k)
                q(k)=q(k)+qq(i,k)
             enddo
             qt=qt+q(k)
  enddo
  qt=qt/(6.0*N*aq)  
return
end subroutine ATPhyd1


!************************************************************************
subroutine binding(sfA,sPha2,W,Amol,Bmol,Pha2,PB,q,istep,c0,c1,c2,c2kaiA,am,rk1,dk1, &
  omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T)
  
  use sim_parameters
  implicit none
  
  integer i,j,k,istep,iloop
  double precision Amol,Bmol,ha2,hb1      
  double precision pt,gkBA,gkA,A0,B0,Abind,Bbind
  double precision sfA,sPha2,shB,sfB,art,brt
  double precision g1,amm,bmm

  double precision,dimension(N) :: W
  double precision,dimension(N) :: q 
  double precision,dimension(N) :: Pha2
  double precision,dimension(N) :: PBt
  double precision,dimension(0:6,N) :: PB
  double precision,dimension(0:6,N) :: PBold    
  double precision,dimension(N) :: Ph0
  double precision,dimension(N) :: hA
  double precision,dimension(N) :: fA
  double precision,dimension(N) :: hB
  double precision,dimension(N) :: fB
  double precision,dimension(6) :: ff
  
  double precision c0,c1,c2,c2kaiA,am,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T
 

  Do k=1,N
     Do i=0,6
        PBold(i,k)=PB(i,k)
     Enddo
  Enddo   
  do k=1,N
       art = tanh(W(k)*arA)
       brt = tanh(W(k)*arB)
       hA(k) = hA0*( 1.0+art )
       fA(k) = fA0*( 1.0-art ) 
       hB(k) = hB0*( 1.0-brt )
       fB(k) = fB0*( 1.0+brt )
  enddo  

  shB=0.0
  sfB=0.0     
  sfA=0.0
  sPha2=0.0
  
  gkBA= hBA/fBA
  g1=gkBA*Amol/(1+gkBA*Amol)

  
  
    do k=1,N

       hb1 = hB(k)*Bmol
       gkA=hA(k)/fA(k)
       
       Do i=1,5
!          PB(i,k)=PBold(i,k)+dt*( -fB(k)*PBold(i,k)-hb1*PBold(i,k) +fB(k)*PBold(i+1,k) + hb1*PBold(i-1,k) )
          PB(i,k)=PBold(i,k)+dt*( -i*fB(k)*PBold(i,k)-(6-i)*hb1*PBold(i,k) +(i+1)*fB(k)*PBold(i+1,k) +(7-i)* hb1*PBold(i-1,k) )
       Enddo
!       PB(6,k)=PBold(6,k)+dt*(-fB(k)*PBold(6,k)+hb1*PBold(5,k))
!       PB(0,k)=PBold(0,k)+dt*(-hb1*PBold(0,k)+fB(k)*PBold(1,k))
       PB(6,k)=PBold(6,k)+dt*(-6*fB(k)*PBold(6,k)+hb1*PBold(5,k))
       PB(0,k)=PBold(0,k)+dt*(-6*hb1*PBold(0,k)+fB(k)*PBold(1,k))      
       PBt(k)=0.0
       Do i=1,6
          PBt(k)=PBt(k)+PB(i,k)
       Enddo      
       Ph0(k)=PB(0,k)

       
       Pha2(k)=Ph0(k)*gkA*Amol/(1+gkA*Amol)       
       pt=PBt(k)+Ph0(k)

       Do i=0,6
          PB(i,k)=PB(i,k)/pt
       Enddo
       PBt(k)=PBt(k)/pt
       Ph0(k)=Ph0(k)/pt
       Pha2(k)=Pha2(k)/pt
    enddo

    A0=0
    B0=0
    Do k=1,N
       gkA=hA(k)/fA(k)
       A0= A0 + Ph0(k)*gkA/(1+gkA*Amol)
       Do i=1,6
          B0 = B0 + real(i)* PB(i,k)
       Enddo
       shB=shB+PBt(k)
    Enddo
    Bbind=B0*g1
    Abind=A0*Amol
    
    if(Atot.ne.0) then
       amm=gkBA*(1+A0)
       bmm=1+A0+gkBA*(B0 -Atot)
       Amol=(-bmm +(bmm**2+4*Atot*amm)**0.5)/(2*amm)
    else
       Amol=0.0
    endif
    Bmol=Btot-B0   

    if(Amol.lt.0.0) Amol=0.0
    if(Amol.gt.Atot) Amol=Atot
    If(Bmol.lt.0.0) Bmol=0.0
    if(Bmol.gt.Btot) Bmol=Btot
    

    sfA=Bbind
    sfB=Bbind
    sPha2=Abind

return
end subroutine binding
    
!************************************************************************

subroutine diffmodel(Pha2,PB,U,W,q,Utot,Wtot,c0,c1,c2,c2kaiA,am,rk1,dk1,omega,delta0, &
  hA0,fA0,hB0,fB0,hBA,fBA,T,iseed,aq)
 
  use sim_parameters
  implicit none
  
  integer i,j,k,iseed
  double precision sb,Utot,Wtot

  double precision,dimension(N) :: U
  double precision,dimension(N) :: W
  double precision,dimension(N) :: q
  double precision,dimension(N) :: Pha2
  double precision,dimension(0:6,N) :: PB  

  double precision c0,c1,c2,c2kaiA,am,rk1,dk1,omega,delta0,hA0,fA0,hB0,fB0,hBA,fBA,T,P6,aq,aq6
  double precision rkk,dkk,skk,r8_normal_01



      Utot=0.0
      Wtot=0.0
      aq6=6*aq
      
      do k=1,N
         if(imode.eq.1) then
            U(k)= U(k)+dt*( rk1*Pha2(k) -dk1*(1.0-Pha2(k))+ am*(0.5*U(k)-U(k)**3))
         endif
         if(imode.eq.2) then
            P6=Pha2(k)
            rkk=dt*rk1*P6/(0.1+P6)*(1-U(k))*0.5      
            dkk=dt*dk1*0.1/(0.1+P6)*(1+U(k))*0.5          
            skk=sqrt(abs(rkk+dkk))*r8_normal_01(iseed)*0.0
            U(k)= U(k)+rkk-dkk+skk
         endif
         if(imode.eq.3) then
            rkk=(1+W(k))*0.5
!            dkk=(1-w(k))*0.5
            U(k)= U(k)+dt*( rk1*rkk*Pha2(k) -dk1*(1.0-Pha2(k))+ am*(0.5*U(k)-U(k)**3))
         endif       


         
         sb=0.0
         Do i=1,6
            sb=sb+PB(i,k)
         Enddo
!         W(k) = tanh(( c0 -c1*U(k) +c2kaiA*Pha2(k) -c2*sb +q(k) -(6*(-aq)-q(k)) )/T)
         W(k) = tanh(( c0 -c1*U(k) +c2kaiA*(Pha2(k)-0.0*Pha2(k)**2) -c2*sb  &
             +q(k)*(1+W(k))/2 -(-6*aq-q(k))*(1-W(k))/2 )/T)        
         Utot=Utot+U(k)
         Wtot=Wtot+W(k)
     enddo

return
end subroutine diffmodel

!***********************************************************************
subroutine dispersion(U,Uave,disp)

  use sim_parameters
  implicit none

  integer i
  double precision Uave,disp
  double precision,dimension(N) :: U

  disp=0.0
  do i=1,N
     disp=disp+(U(i)-Uave)**2
  enddo
  disp=disp/dble(N)

  

return
end subroutine dispersion


!***********************************************************************
subroutine contourplot(U,P2,P1a,P1b,ipair1,ipair2)

  use sim_parameters
  implicit none

  integer i,j,ipa,ipb,ima,imb
  double precision,dimension(N) :: U
  double precision,dimension(0:mesh,0:mesh,10,10) :: P2
  double precision,dimension(0:mesh,10) ::P1a
  double precision,dimension(0:mesh,10) ::P1b  
  integer,dimension(10) :: ipair1,ipair2

  do i=1,10
     ipa=ipair1(i)
     ima=int((U(ipa)+1)*0.5*mesh)
     if(ima.le.0) ima=0
     if(ima.ge.mesh) ima=mesh
     P1a(ima,i)=P1a(ima,i)+1.0

     ipb=ipair2(i)
     imb=int((U(ipb)+1)*0.5*mesh)
     if(imb.le.0) imb=0
     if(imb.ge.mesh) imb=mesh
     P1b(imb,i)=P1b(imb,i)+1.0
  enddo

  do i=1,10
     ipa=ipair1(i)
     ima=int((U(ipa)+1)*0.5*mesh)
     if(ima.le.0) ima=0
     if(ima.ge.mesh) ima=mesh
     do j=1,10
        ipb=ipair2(j)
        imb=int((U(ipb)+1)*0.5*mesh)
        if(imb.le.0) imb=0
        if(imb.ge.mesh) imb=mesh
        P2(ima,imb,i,j)=P2(ima,imb,i,j)+1.0
     enddo
  enddo 

return
end subroutine contourplot


!***********************************************************************

subroutine mutinf(P2,P1a,P1b,umut,istep)

  use sim_parameters
  implicit none

  integer i,j,istep,imesh,jmesh
  double precision umut
  double precision,dimension(0:mesh,0:mesh,10,10) :: P2
  double precision,dimension(0:mesh,10) :: P1a
  double precision,dimension(0:mesh,10) :: P1b  
  double precision,dimension(0:mesh,0:mesh) :: PT2
  double precision,dimension(0:mesh) :: PT1a
  double precision,dimension(0:mesh) :: PT1b

  umut=0.0
  do imesh=0,mesh
     PT1a(imesh)=0.0
     PT1b(imesh)=0.0
     do jmesh=0,mesh
        PT2(imesh,jmesh)=0.0
     enddo
  enddo


  do imesh=0,mesh
     do i=1,10
        PT1a(imesh)=PT1a(imesh)+P1a(imesh,i)*0.1/dble(istep-100000)
        PT1b(imesh)=PT1b(imesh)+P1b(imesh,i)*0.1/dble(istep-100000)
     enddo
     do jmesh=0,mesh
        do i=1,10
           do j=1,10
              PT2(imesh,jmesh)=PT2(imesh,jmesh)+P2(imesh,jmesh,i,j)*0.01/dble(istep-100000)
           enddo
        enddo
     enddo    
  enddo

  do imesh=0,mesh
     do jmesh=0,mesh
        if(PT1a(imesh).gt.0.0.and.PT1b(jmesh).gt.0.0.and.PT2(imesh,jmesh).gt.0.0) then
           umut=umut+PT2(imesh,jmesh)*log(PT2(imesh,jmesh)/PT1a(imesh)/PT1b(jmesh))
        endif
     enddo
  enddo

  open(unit=71, file='pt2.dat')
  do imesh=0,mesh
     do jmesh=0,mesh
        write(71,*) imesh,jmesh,PT2(imesh,jmesh)
     enddo
     write(71,*)
!                write(71,*) imesh,PT1b(imesh)
  enddo
  close(71)

  
  return
end subroutine mutinf

!***********************************************************************

subroutine mutinf2(P2,P1a,P1b,umut,istep)

  use sim_parameters
  implicit none

  integer i,j,istep,imesh,jmesh
  double precision umut
  double precision,dimension(0:mesh,0:mesh,10,10) :: P2
  double precision,dimension(0:mesh,10) :: P1a
  double precision,dimension(0:mesh,10) :: P1b  
  double precision PT2,PT1a,PT1b
  umut=0.0

  do i=1,10
     do imesh=0,mesh
        PT1a=P1a(imesh,i)/dble(istep-100000)
        
        do j=1,10
           do jmesh=0,mesh
              PT1b=P1b(jmesh,j)/dble(istep-100000)
              PT2=P2(imesh,jmesh,i,j)/dble(istep-100000)
              if(PT1a.gt.0.0.and.PT1b.gt.0.0.and.PT2.gt.0.0) then             
                 umut=umut + PT2*log(PT2/PT1a/PT1b)
              endif
           enddo
        enddo
        
     enddo    
  enddo


  umut=umut*0.01

  
  return
end subroutine mutinf2


!************************************************************************
!
      function r8_normal_01 (seed)
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.

!    Input/output, integer SEED, a seed for the random number generator.
!    Output, double precision R8_NORMAL_AB, a sample of the normal PDF.

      implicit none

      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      double precision x

      r1 = r8_uniform_01 ( seed )
      r2 = r8_uniform_01 ( seed )
      r8_normal_01 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )


      return
      end

!*********************************************************************
!
      function r8_uniform_01 ( seed )
!
!  R8_UNIFORM_01 returns a unit pseudorandom R8.

!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702

!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
        

        
      implicit none

      integer i4_huge 
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10


      return
      end
    
