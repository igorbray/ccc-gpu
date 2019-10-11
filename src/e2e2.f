      subroutine E2E2
     : (lg,gk,npk,wk,vt,nchtop,ovlpn,csfile,slowery,phasel,phaseq)

C Integration of the Born ME with half-on-shell T-matrix
      
      include 'par.f'
      parameter (knm3 = knm * 3, ncs = (lamax+1)*3)
                                                             
      COMMON/BORN/ br(0:lamax,kmax,nchan),bv(0:lamax,kmax,nchan)
      COMMON/BM2/ br2(0:lamax,kmax,nchan),bv2(0:lamax,kmax,nchan)
      COMMON/MJ/   MJ
      COMMON/BIM2/ Hr2(0:lamax,nchan),  Hv2(0:lamax,nchan)
      COMMON/channel/ iopen(nchan)    !Indicates open channels (=1) 
      COMMON /faza/ ph(0:lamax,nchan)
      complex ph,phasel(kmax,nchan),phase, phaseq(knm)
      COMMON /moment/ qq,omega
                                                             
      character*3 chan(knm), chs(0:lamax)
      character*3 stat
      
      common /charchan/ chan
      common /chanen/ enchan(knm)
      common /pspace/ nabot(0:lamax),labot,natop(0:lamax),latop,
     >   ntype,ipar,nze,ninc,linc,lactop,nznuc,zasym
      character lh(0:10), nh(0:9),csfile*(*), file*50
      data lh / 's','p','d','f','g','h','n','o','x','x','x'/
      data nh / '0','1','2','3','4','5','6','7','8','9'/

      integer npk(nchtop+1)
      dimension gk(kmax,nchan), temp(maxr), ovlpn(knm), cs(knm,3),
     >   tcs(0:lamax,3),tics(0:lamax,3),tnbcs(0:lamax,3)

      complex wk(kmax*nchan)    !Greens function
      complex vt(npk(nchtop+1)-1,nchtop) !T-matrix
      complex GF,T,ci,TR,TV,dTR,dTV
      complex Hr2, Hv2
      complex*16 coulphase
      real*8 deta

      data ltop,ry,cs,tnbcs,tics,tcs/0,13.6058,knm3*0.0,
     >   ncs*0.0, ncs*0.0, ncs*0.0/

      if(nznuc.eq.2)GS  =  5.807 !He  atom GS energy in Ry
      if(nznuc.eq.3)GS  = 14.559 !Li+  ion GS energy in Ry
      if(nznuc.eq.4) GS = 2.0236 !BeIII-BeI   energy in Ry Radzig, Smirnov
      if(nznuc.eq.8)GS  =118.312 !O^6+ ion GS energy in Ry
      if(nznuc.eq.1)GS  =1.05544 !H-   ion GS energy in Ry
              
      pi = acos(-1.0)                                       
      ci = cmplx(0.0, 1.0)

      
      i = 1
      do while (csfile(i:i).ne.'_'.and.csfile(i:i).ne.' ')
         i = i + 1
      enddo 

      if (i.eq.8) then
        open(42,file='L2borncs.J='//nh(lg)//csfile(i:60))
        open(43,file='V2borncs.J='//nh(lg)//csfile(i:60))
!        open(42,file='L2borncs.J='//nh(lg)//nh(MJ))
      end if
      
      write(6,1000)
C Final channel loop
      
    
      ntmax = 0
      do nchf = 1, nchtop
         call getchinfo (nchf,nt,lg,temp,maxpsif,ef,lf,nf,llf)
         EI=GS+ef                                        
         if (nchf.eq.1) then
            etot = ef + gk(1,nchf)**2
            om = EI + gk(1,nchf)**2
               write(42,51)qq,om*ry, etot*ry
               write(43,51)qq,om*ry, etot*ry
               write (42,52) slowery*ry
               write (43,52) slowery*ry
         endif 

 51    format('Transferred momentum:  ' ,f7.4,' au'/,
     >        'Transferred energy:    ' ,f7.2,' eV'/,
     >        'Electron pair energy:  ' ,f7.2,' eV')

 52    format('Inner electron energy: ' ,f7.2,' eV'//,
     >        'Transition:    J   Lin Lout   Re(T)       Im(T) ')


         ph(lg,nchf) = phasel(1,nchf)
         phi = atan(aimag(ph(lg,nchf))/real(ph(lg,nchf)))
         

         nqm = npk(nchf+1) - npk(nchf)
         write(6,1954) nchf,EI,nf,lh(lf),lh(LLf),nqm,phi
         write(6,'(11x,A,F9.4)') 'Bound state energy',ef
         ipseudo=0              !Eigenstates
         if(ef.gt.0.) ipseudo=1 !Pseudostates
         if(IOPEN(nchf).eq.0)then
            write(6,'(11x,A,//)') 'Channel closed'
            goto 113
         end if
                                                             
         TR =(0.0,0.0)                                          
         TV =(0.0,0.0)                                          
         
c$$$         write(6,35)                                            

C Initial channel loop
      
         nv = 0
         do nchi = 1, nchtop
            nqm = npk(nchi+1) - npk(nchi)
            dTR =(0.,0.)
            dTV =(0.,0.)
            do k = 1, nqm
               nv = nv + 1
               q0 = gk(1,nchf)
               q  = gk(k,nchi)
               if (q.lt.0.0) then
                  tmp = 1.0
               else
                  tmp = q
               end if 
               if(k.eq.1.and.q.lt.0.0)goto 112
               if(q.gt.0.0)then                              
                  rn=sqrt(2./q)                              
               else                                          
                  rn=1.0                                     
               end if                                        
               r =br2 (lg,k,NCHI)*rn                             
               v =bv2 (lg,k,NCHI)*rn
               GF=  wk(nv)* tmp**2
               T =vt(nv,nchf)/q0/tmp
c$$$*               write(6,40) q,gf,T,r,v,x
               if(k.eq.1.and.NCHI.eq.NCHF)then !Diagonal term f-f
                  p =q
                  ro=r 
                  vo=v
               end if           
               dTR=dTR+r*GF*T                                
               dTV=dTV+v*GF*T
 112           continue                                         
            end do
            
!            write(6,41) NCHI,dTR
!            if(nchi.eq.nchf) then !Diagonal terms only
               TR = (TR + dTR)
               TV = (TV + dTV)
!            end if
         end do

!         write(6,42) 0, ro,0d0
!         write(6,43) TR+ro
         
         TR = (TR+ro) 
         TV = (TV+vo)

C Store integrated Born ME
C Restore normalization
C Multiply by phase factor and overlap
         
         phase = phasel(1,nchf)
         if (ef.gt.0.0) then
            deta = - (zasym + 1.0) / sqrt(ef)
c$$$            phase = phase * coulphase(deta,lf) * ovlpn(nt)
            phase = phase * phaseq(nt) * ovlpn(nt)
c$$$            print*,'phases:',phaseq(nt),coulphase(deta,lf)
         endif

         Hr2(lg,nchf) = TR*phase*sqrt(p/2)
         Hv2(lg,nchf) = TV*phase*sqrt(p/2)

         if (nt.gt.ntmax) ntmax = nt
         if (lf.gt.ltop) ltop = lf
         
 113     continue
      end do

      print'(A)', 'STEST'
      
      do nchf = 1, nchtop
         call getchinfo (nchf,nt,lg,temp,maxpsif,ef,lf,nf,llf)
         if(IOPEN(nchf).eq.1)then
      write(42,'(I2,A1,A,A1,5x,3i4,4x,1p,2e12.4,0p,2F9.4)')
     >      nf,lh(lf),'->E',lh(llf),  lg,lf,llf,Hr2(lg,nchf),ef*Ry,
     >      ovlpn(nt)   
      write(43,'(I2,A1,A,A1,5x,3i4,4x,1p,2e12.4,0p,2F9.4)')
     >      nf,lh(lf),'->E',lh(llf),  lg,lf,llf,Hv2(lg,nchf),ef*Ry,
     >      ovlpn(nt)   
      end if
      end do



 35   format(/'Initial',8x,'Length',19x,'Velocity'/,
     >   'channel',3x,
     >   2('Re            Im          ')/)     
 40   format(10(E11.4,2x))                                             
 41   format(I3,2x,6(E11.4,2x))                                        
 42   format(/I3,2x,6(E11.4,2x))                                       
 43   format(' Sum ',6(E11.4,2x))                                      
 100  format(10(3x,E11.4) )                                            
 1000 FORMAT(//15X,'INTEGRATION OF ',                                  
     >   'THE 2nd BORN ME WITH HALF ON-SHELL T-MATRIX'/,15x,59('-')/)    
 1954 format(/11x,'FINAL CHANNEL ',I2/,                                
     >        11x,'Ion energy',F7.3/,                                  
     >        11x,'Bound electron',I2,A1/,                            
     >        11x,'Free  electron k', A1/,                             
     >        11x,'Total k     ', I2/,
     >        11x,'Coulomb phase   ', F9.4/)

 1955 format(/13x,'1st Cross Section',             
     >       11x,'TM Cross Section',/             
     >       5x,2(6x,'Length',7x,'Velocty',2x)/,61('-'))

      close(42)
      close(43)
 
      return
      end

