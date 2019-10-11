      subroutine E2E
     : (lg,gk,npk,wk,vt,nchtop,ovlpn,csfile,slowery,phasel,phaseq)

C Integration of the Born ME with half-on-shell T-matrix
      
      include 'par.f'
      parameter (knm3 = knm * 3, ncs = (lamax+1)*3)
                                                             
      COMMON/BORN/ br(0:lamax,kmax,nchan),bv(0:lamax,kmax,nchan)
      COMMON/BORN2r/ Hr(0:lamax,nchan)
      COMMON/BORN2v/ Hv(0:lamax,nchan)
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
      complex Hr,Hv
      complex*16 coulphase
      real*8 deta

      data ltop,ry,cs,tnbcs,tics,tcs/0,13.6058,knm3*0.0,
     >   ncs*0.0, ncs*0.0, ncs*0.0/

      if(nznuc.eq.2)GS  =  5.807 !He  atom GS energy in Ry
      if(nznuc.eq.3)GS  = 14.559 !Li+  ion GS energy in Ry
      if(nznuc.eq.8)GS  =118.312 !O^6+ ion GS energy in Ry
      if(nznuc.eq.1)GS  =1.05544 !H-   ion GS energy in Ry
      if(nznuc.eq.12)GS = 1.6670 !MgIII-MgI   energy in Ry Radzig, Smirnov
      pi = acos(-1.0)                                       
      ci = cmplx(0.0, 1.0)

      
      open(42,file='Lborncs.J='//nh(lg)//csfile(8:60))
      open(43,file='Vborncs.J='//nh(lg)//csfile(8:60))

      write(6,1000)
C Final channel loop
      
      sigma_l =0.0
      sigma_v =0.0
      sigma1_l =0.0
      sigma1_v =0.0
    
      ntmax = 0
      do nchf = 1, nchtop
         call getchinfo (nchf,nt,lg,temp,maxpsif,ef,lf,nf,llf)
         EI=GS+ef                                        
         if (nchf.eq.1) then
            etot = ef + gk(1,nchf)**2
            om = EI + gk(1,nchf)**2
            do ng = 1, 2
               write(41+ng,51)qq,om*ry, etot*ry
               write (41+ng,52) slowery*ry
*               if (slowery.gt.0.0) write (41+ng,52) slowery*ry
            enddo
         endif 

 51    format('Transferred momentum:  ' ,f7.4,' au'/,
     >        'Transferred energy:    ' ,f7.2,' eV'/,
     >        'Electron pair energy:  ' ,f7.2,' eV')

 52    format('Inner electron energy: ' ,f7.2,' eV'//,
     >        'Transition:    J   Lin Lout   Re(T)       Im(T) ')


         ph(lg,nchf) = phasel(1,nchf)
         phi = atan(aimag(ph(lg,nchf))/real(ph(lg,nchf)))
         

         nqm = npk(nchf+1) - npk(nchf)
*         write(6,1954) nchf,EI,nf,lh(min(lf,10)),lh(min(LLf,10)),
*     >      nqm,phi
*         write(6,'(11x,A,F9.4)') 'Bound state energy',ef
         ipseudo=0              !Eigenstates
         if(ef.gt.0.) ipseudo=1 !Pseudostates
         if(IOPEN(nchf).eq.0)then
*            write(6,'(11x,A,//)') 'Channel closed'
            goto 113
         end if
                                                             
         TR =(0.0,0.0)                                          
         TV =(0.0,0.0)                                          
         
c$$$         write(6,35)                                            

C Initial channel loop
      
         nv = 0
         const=8.067 !  4 pi^2 * alpha * a0^2 * 10^{18}, alpha = 1/137 

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
               r =br(lg,k,NCHI)*rn                             
               v =bv(lg,k,NCHI)*rn                             
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
            
c$$$            write(6,41) NCHI,dTR,dTV
            TR = (TR + dTR)
            TV = (TV + dTV)
         end do

c$$$         write(6,42) 0, ro,0d0,vo,0d0
c$$$         write(6,43) TR+ro,TV+vo
         
C Restore normalization

         const=const*p/2                                        
         om=EI+p**2

         TR = (TR+ro) 
         TV = (TV+vo) 
         
c$$$         write(6,'(/A,F9.4/)') 'OMEGA=', om*ry

         CLo=  const*2/3*om*abs(ro)**2                          
         CVo=4*const*2/3/om*abs(vo)**2                          
         CL =  const*2/3*om*abs(TR)**2                          
         CV =4*const*2/3/om*abs(TV)**2                          

C Store integrated Born ME
C Restore normalization
C Multiply by phase factor and overlap
         
         phase = phasel(1,nchf)
         if (ef.gt.0.0) then
            deta = - 2.0 / sqrt(ef)
c$$$            phase = phase * coulphase(deta,lf) * ovlpn(nt)
            phase = phase * phaseq(nt) * ovlpn(nt)
         endif

         Hr(lg,nchf) = TR*phase*sqrt(p/2)
         Hv(lg,nchf) = TV*phase*sqrt(p/2)

C  Here to write out the double ionization matrix elements.
C  The asymptotic charge seen by the inner target-space electron is +2. OVLP is
C  the overlap between the EF-energy state with the true continuum function
C  of same energy. This contains the continuum normalisation sqrt(2/pi), but
C  NOT division by sqrt(ef).
            
*         if (ef.gt.0.0.and.abs(ef-slowery)/ef.lt.1e-1) then
*            write(42,70) lg,lf,llf,Hr(lg,nchf)
*            write(43,70) lg,lf,llf,Hv(lg,nchf)
*         endif
* 70      format(12x,3i4,4x,1p,2e12.4)
                  

c$$$         write(6,1955)                                          
c$$$         write(6,'(I2,A1,2x,10(3x,E11.4) )')nf,lh(lf),
c$$$     >                CLo,CVo,CL,CV, ovlpn(nt)

         sigma_l = sigma_l+CL
         sigma_v = sigma_v+Cv

         if (nt.gt.ntmax) ntmax = nt
         if (lf.gt.ltop) ltop = lf
         
         if(ipseudo.eq.0)then       !Negative energy states only
            proj = ovlpn(nt)
            if (proj.lt.0.0.or.proj.gt.1.0001)
     >         print*, 'CAUTION: OVLPN for CHANNEL:',proj,nt
            sigma1_l = sigma1_l+CL * proj
            sigma1_v = sigma1_v+Cv * proj
         end if
 113     continue
      end do

*      do ngauge = 1, 2
*         write(41+ngauge,53)
* 53      format(52('=')/,'Transition: ' ,
*     >        '   J   Lin Lout   Re(T)       Im(T) ')
*      end do
      do nchf = 1, nchtop
         call getchinfo (nchf,nt,lg,temp,maxpsif,ef,lf,nf,llf)
         if(IOPEN(nchf).eq.1)then
      write(42,'(I2,A1,A,A1,5x,3i4,4x,1p,2e12.4,0p,F9.4)')
     >      nf,lh(min(lf,10)),'->E',lh(min(10,llf)),  lg,lf,llf,
     >         Hr(lg,nchf),ef*Ry
         write(43,'(I2,A1,A,A1,5x,3i4,4x,1p,2e12.4,0p,F9.4)')
     >      nf,lh(min(lf,10)),'->E',lh(min(10,llf)),  lg,lf,llf,
     >   Hv(lg,nchf),ef*Ry
      end if
      end do



      asymcs = 0.0
      do l = 0, ltop
         chs(l) = ' + '
         if (l.eq.0) chs(l) = ' = '
      enddo 
      
      ein = gk(1,1)**2*ry
!      ddl = (sigma_l-sigma1_l)/sigma1_l*100.  
!      ddv = (sigma_v-sigma1_v)/sigma1_v*100.  
!      ddx = (sigma_x-sigma1_x)/sigma1_x*100.  
      
!      write(20,'(F7.2,3F11.4)') om*ry,sigma_l, sigma_v, sigma_x   
!      write(21,'(F7.2,3F11.4)') om*ry,sigma1_l,sigma1_v, sigma1_x
!      write(22,'(''#ph energy  double L single     double V single'',
!     >   ''     double A single   el energy'')') 
!      write(22,'(F8.2,1x,1p,6E10.3,0p,f10.3)') om*ry, sigma_l-sigma1_l,
!     >   sigma1_l,sigma_v-sigma1_v,sigma1_v,sigma_x-sigma1_x,sigma1_x,
!     >   ein
!      write(22,'(F8.2,3F11.4,2x,a20)') om*ry, ddl, ddv,ddx,file
 35   format(/'Initial',8x,'Length',19x,'Velocity'/,
     >   'channel',3x,
     >   2('Re            Im          ')/)     
 40   format(10(E11.4,2x))                                             
 41   format(I3,2x,6(E11.4,2x))                                        
 42   format(/I3,2x,6(E11.4,2x))                                       
 43   format(' Sum ',6(E11.4,2x))                                      
 100  format(10(3x,E11.4) )                                            
 1000 FORMAT(//15X,'INTEGRATION OF ',                                  
     >   'THE BORN ME WITH HALF ON-SHELL T-MATRIX'/,15x,56('-')/)    
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

