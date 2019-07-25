c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program testrot  --  derivative test; torsional version  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "testrot" computes and compares the analytical and numerical
c     gradient vectors of the potential energy function with respect
c     to rotatable torsional angles
c
c
      program testrot
      use sizes
      use domega
      use energi
      use inform
      use iounit
      use math
      use omega
      use zcoord
      implicit none
      integer i
      real*8 e,e0,etot,energy
      real*8 delta,delta0,eps
      real*8 eb0,ea0,eba0,eub0
      real*8 eaa0,eopb0,eopd0
      real*8 eid0,eit0,et0
      real*8 ept0,ebt0,eat0
      real*8 ett0,ev0,ec0,ecd0
      real*8 ed0,em0,ep0,er0
      real*8 es0,elf0,eg0,ex0
      real*8, allocatable :: derivs(:)
      real*8, allocatable :: nderiv(:)
      logical exist,query
      character*240 string
c
c
c     set up the molecular mechanics calculation
c
      call initial
      call getint
      call mechanic
      call initrot
c
c     get the stepsize for numerical gradient calculation
c
      delta = -1.0d0
      delta0 = 1.0d-3
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  delta
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)  delta0
   20    format (/,' Enter Finite Difference Stepsize [',d8.1,
     &              ' Deg] :  ',$)
         read (input,30,err=10)  delta
   30    format (f20.0)
      end if
      if (delta .le. 0.0d0)  delta = delta0
      eps = -delta / radian
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(nomega))
c
c     make the call to get analytical torsional derivatives
c
      call gradrot (etot,derivs)
c
c     print the total potential energy of the system
c
      if (digits .ge. 8) then
         write (iout,40)  etot
   40    format (/,' Total Potential Energy :',8x,f20.8,' Kcal/mole')
      else if (digits .ge. 6) then
         write (iout,50)  etot
   50    format (/,' Total Potential Energy :',8x,f18.6,' Kcal/mole')
      else
         write (iout,60)  etot
   60    format (/,' Total Potential Energy :',8x,f16.4,' Kcal/mole')
      end if
c
c     print the energy breakdown over individual components
c
      write (iout,70)
   70 format (/,' Potential Energy Breakdown by Individual',
     &           ' Components :')
      if (digits .ge. 8) then
         write (iout,80)
   80    format (/,'  Energy',7x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &           /,'  Terms',8x,'EAA',13x,'EOPB',12x,'EOPD',12x,'EID',
     &           /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
     &           /,15x,'EAT',13x,'ETT',13x,'EV',14x,'EC',
     &           /,15x,'ECD',13x,'ED',14x,'EM',14x,'EP',
     &           /,15x,'ER',14x,'ES',14x,'ELF',13x,'EG',
     &           /,15x,'EX')
         write (iout,90)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
     &                    eit,et,ept,ebt,eat,ett,ev,ec,ecd,
     &                    ed,em,ep,er,es,elf,eg,ex
   90    format (/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,
     &              /,6x,4f16.8,/,6x,4f16.8,/,6x,f16.8)
      else if (digits .ge. 6) then
         write (iout,100)
  100    format (/,'  Energy',6x,'EB',12x,'EA',12x,'EBA',11x,'EUB',
     &              11x,'EAA',
     &           /,'  Terms',7x,'EOPB',10x,'EOPD',10x,'EID',
     &              11x,'EIT',11x,'ET',
     &           /,14x,'EPT',11x,'EBT',11x,'EAT',11x,'ETT',11x,'EV',
     &           /,14x,'EC',12x,'ECD',11x,'ED',12x,'EM',12x,'EP',
     &           /,14x,'ER',12x,'ES',12x,'ELF',11x,'EG',12x,'EX')
         write (iout,110)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
     &                     eit,et,ept,ebt,eat,ett,ev,ec,ecd,
     &                     ed,em,ep,er,es,elf,eg,ex
  110    format (/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,/,6x,5f14.6,
     &              /,6x,5f14.6)
      else
         write (iout,120)
  120    format (/,'  Energy',6x,'EB',10x,'EA',10x,'EBA',9x,'EUB',
     &              9x,'EAA',9x,'EOPB',
     &           /,'  Terms',7x,'EOPD',8x,'EID',9x,'EIT',9x,'ET',
     &              10x,'EPT',9x,'EBT',
     &           /,14x,'EAT',9x,'ETT',9x,'EV',10x,'EC',10x,'ECD',
     &              9x,'ED',
     &           /,14x,'EM',10x,'EP',10x,'ER',10x,'ES',10x,'ELF',
     &              9x,'EG',
     &           /,14x,'EX')
         write (iout,130)  eb,ea,eba,eub,eaa,eopb,eopd,eid,
     &                     eit,et,ept,ebt,eat,ett,ev,ec,ecd,
     &                     ed,em,ep,er,es,elf,eg,ex
  130    format (/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4,
     &               /,6x,f12.4)
      end if
c
c     print a header for the gradients of individual potentials
c
      write (iout,140)
  140 format (/,' Torsional Gradient Breakdown by Individual',
     &              ' Components :')
      if (digits .ge. 8) then
         write (iout,150)
  150    format (/,2x,'Atom',9x,'d EB',12x,'d EA',12x,'d EBA',
     &              11x,'d EUB',
     &           /,2x,'Axis',9x,'d EAA',11x,'d EOPB',10x,'d EOPD',
     &              10x,'d EID',
     &           /,2x,'Type',9x,'d EIT',11x,'d ET',12x,'d EPT',
     &              11x,'d EBT',
     &           /,15x,'d EAT',11x,'d ETT',11x,'d EV',12x,'d EC',
     &           /,15x,'d ECD',11x,'d ED',12x,'d EM',12x,'d EP',
     &           /,15x,'d ER',12x,'d ES',12x,'d ELF',11x,'d EG',
     &           /,15x,'d EX')
      else if (digits .ge. 6) then
         write (iout,160)
  160    format (/,2x,'Atom',8x,'d EB',10x,'d EA',10x,'d EBA',
     &              9x,'d EUB',9x,'d EAA',
     &           /,2x,'Axis',8x,'d EOPB',8x,'d EOPD',8x,'d EID',
     &              9x,'d EIT',9x,'d ET',
     &           /,2x,'Type',8x,'d EPT',9x,'d EBT',9x,'d EAT',
     &              9x,'d ETT',9x,'d EV',
     &           /,14x,'d EC',10x,'d ECD',9x,'d ED',10x,'d EM',
     &              10x,'d EP',
     &           /,14x,'d ER',10x,'d ES',10x,'d ELF',9x,'d EG',
     &              10x,'d EX')
      else
         write (iout,170)
  170    format (/,2x,'Atom',6x,'d EB',8x,'d EA',8x,'d EBA',
     &              7x,'d EUB',7x,'d EAA',7x,'d EOPB',
     &           /,2x,'Axis',6x,'d EOPD',6x,'d EID',7x,'d EIT',
     &              7x,'d ET',8x,'d EPT',7x,'d EBT',
     &           /,2x,'Type',6x,'d EAT',7x,'d ETT',7x,'d EV',
     &              8x,'d EC',8x,'d ECD',7x,'d ED',
     &           /,12x,'d EM',8x,'d EP',8x,'d ER',8x,'d ES',
     &              8x,'d ELF',7x,'d EG',
     &           /,12x,'d EX')
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (nderiv(nomega))
c
c     get numerical derivatives for each of the rotatable torsions
c
      do i = 1, nomega
         ztors(zline(i)) = ztors(zline(i)) + delta/2.0d0
         call makexyz
         e0 = energy ()
         eb0 = eb
         ea0 = ea
         eba0 = eba
         eub0 = eub
         eaa0 = eaa
         eopb0 = eopb
         eopd0 = eopd
         eid0 = eid
         eit0 = eit
         et0 = et
         ept0 = ept
         ebt0 = ebt
         eat0 = eat
         ett0 = ett
         ev0 = ev
         ec0 = ec
         ecd0 = ecd
         ed0 = ed
         em0 = em
         ep0 = ep
         er0 = er
         es0 = es
         elf0 = elf
         eg0 = eg
         ex0 = ex
         ztors(zline(i)) = ztors(zline(i)) - delta
         call makexyz
         e = energy ()
         ztors(zline(i)) = ztors(zline(i)) + delta/2.0d0
         nderiv(i) = (e-e0) / eps
c
c     print analytical gradients of each energy term for each atom
c
         if (digits .ge. 8) then
            write (iout,180)  iomega(2,i),teb(i),tea(i),teba(i),teub(i),
     &                        iomega(1,i),teaa(i),teopb(i),teopd(i),
     &                        teid(i),teit(i),tet(i),tept(i),tebt(i),
     &                        teat(i),tett(i),tev(i),tec(i),tecd(i),
     &                        ted(i),tem(i),tep(i),ter(i),tes(i),
     &                        telf(i),teg(i),tex(i)
  180       format (/,i6,4f16.8,/,i6,4f16.8,/,' Anlyt',4f16.8,
     &                 /,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,f16.8)
         else if (digits .ge. 6) then
            write (iout,190)  iomega(2,i),teb(i),tea(i),teba(i),teub(i),
     &                        teaa(i),iomega(1,i),teopb(i),teopd(i),
     &                        teid(i),teit(i),tet(i),tept(i),tebt(i),
     &                        teat(i),tett(i),tev(i),tec(i),tecd(i),
     &                        ted(i),tem(i),tep(i),ter(i),tes(i),
     &                        telf(i),teg(i),tex(i)
  190       format (/,i6,5f14.6,/,i6,5f14.6,/,' Anlyt',5f14.6,
     &                 /,6x,5f14.6,/,6x,5f14.6)
         else
            write (iout,200)  iomega(2,i),teb(i),tea(i),teba(i),teub(i),
     &                        teaa(i),teopb(i),iomega(1,i),teopd(i),
     &                        teid(i),teit(i),tet(i),tept(i),tebt(i),
     &                        teat(i),tett(i),tev(i),tec(i),tecd(i),
     &                        ted(i),tem(i),tep(i),ter(i),tes(i),
     &                        telf(i),teg(i),tex(i)
  200       format (/,i6,6f12.4,/,i6,6f12.4,/,' Anlyt',6f12.4,
     &                 /,6x,6f12.4,/,6x,f12.4)
         end if
c
c     print numerical gradients of each energy term for each atom
c
         if (digits .ge. 8) then
            write (iout,210)  iomega(2,i),(eb-eb0)/eps,(ea-ea0)/eps,
     &                        (eba-eba0)/eps,(eub-eub0)/eps,
     &                        iomega(1,i),(eaa-eaa0)/eps,
     &                        (eopb-eopb0)/eps,(eopd-eopd0)/eps,
     &                        (eid-eid0)/eps,(eit-eit0)/eps,
     &                        (et-et0)/eps,(ept-ept0)/eps,
     &                        (ebt-ebt0)/eps,(eat-eat0)/eps,
     &                        (ett-ett0)/eps,(ev-ev0)/eps,(ec-ec0)/eps,
     &                        (ecd-ecd0)/eps,(ed-ed0)/eps,(em-em0)/eps,
     &                        (ep-ep0)/eps,(er-er0)/eps,(es-es0)/eps,
     &                        (elf-elf0)/eps,(eg-eg0)/eps,(ex-ex0)/eps
  210       format (/,i6,4f16.8,/,i6,4f16.8,/,' Numer',4f16.8,
     &                 /,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8,/,6x,f14.8)
         else if (digits .ge. 6) then
            write (iout,220)  iomega(2,i),(eb-eb0)/eps,(ea-ea0)/eps,
     &                        (eba-eba0)/eps,(eub-eub0)/eps,
     &                        (eaa-eaa0)/eps,iomega(1,i),
     &                        (eopb-eopb0)/eps,(eopd-eopd0)/eps,
     &                        (eid-eid0)/eps,(eit-eit0)/eps,
     &                        (et-et0)/eps,(ept-ept0)/eps,
     &                        (ebt-ebt0)/eps,(eat-eat0)/eps,
     &                        (ett-ett0)/eps,(ev-ev0)/eps,(ec-ec0)/eps,
     &                        (ecd-ecd0)/eps,(ed-ed0)/eps,(em-em0)/eps,
     &                        (ep-ep0)/eps,(er-er0)/eps,(es-es0)/eps,
     &                        (elf-elf0)/eps,(eg-eg0)/eps,(ex-ex0)/eps
  220       format (/,i6,5f14.6,/,i6,5f14.6,/,' Numer',5f14.6,
     &                 /,6x,5f14.6,/,6x,5f14.6)
         else
            write (iout,230)  iomega(2,i),(eb-eb0)/eps,(ea-ea0)/eps,
     &                        (eba-eba0)/eps,(eub-eub0)/eps,
     &                        (eaa-eaa0)/eps,(eopb-eopb0)/eps,
     &                        iomega(1,i),(eopd-eopd0)/eps,
     &                        (eid-eid0)/eps,(eit-eit0)/eps,
     &                        (et-et0)/eps,(ept-ept0)/eps,
     &                        (ebt-ebt0)/eps,(eat-eat0)/eps,
     &                        (ett-ett0)/eps,(ev-ev0)/eps,(ec-ec0)/eps,
     &                        (ecd-ecd0)/eps,(ed-ed0)/eps,(em-em0)/eps,
     &                        (ep-ep0)/eps,(er-er0)/eps,(es-es0)/eps,
     &                        (elf-elf0)/eps,(eg-eg0)/eps,(ex-ex0)/eps
  230       format (/,i6,6f12.4,/,i6,6f12.4,/,' Numer',6f12.4,
     &                 /,6x,6f12.4,/,6x,f12.4)
         end if
      end do
c
c     print a header for the analytical vs. numerical comparison
c
      if (digits .ge. 8) then
         write (iout,240)
  240    format (//,5x,'Torsion',19x,'Anlyt Deriv',9x,'Numer Deriv',/)
      else if (digits .ge. 6) then
         write (iout,250)
  250    format (//,5x,'Torsion',18x,'Anlyt Deriv',7x,'Numer Deriv',/)
      else
         write (iout,260)
  260    format (//,5x,'Torsion',17x,'Anlyt Deriv',5x,'Numer Deriv',/)
      end if
c
c     print comparison of analytical and numerical derivatives
c
      if (digits .ge. 8) then
         do i = 1, nomega
            write (iout,270)  iomega(2,i),iomega(1,i),derivs(i),
     &                        nderiv(i)
  270       format (1x,i5,'-',i5,10x,2f20.8)
         end do
      else if (digits .ge. 6) then
         do i = 1, nomega
            write (iout,280)  iomega(2,i),iomega(1,i),derivs(i),
     &                        nderiv(i)
  280       format (1x,i5,'-',i5,10x,2f18.6)
         end do
      else
         do i = 1, nomega
            write (iout,290)  iomega(2,i),iomega(1,i),derivs(i),
     &                        nderiv(i)
  290       format (1x,i5,'-',i5,10x,2f16.4)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      deallocate (nderiv)
c
c     perform any final tasks before program exit
c
      call final
      end
