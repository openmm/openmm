C  ******************************************************************************
C                                     DISCLAIMER
C                    Everything that is free comes with NO warranty!!!
C                                   Pedro E M Lopes
C  ******************************************************************************
      program patchfind
      implicit none

      integer mxatm
      parameter (mxatm=1000)
      integer mxcon
      parameter (mxcon=8)
      integer i,j,k,l
      integer n,nunit
      integer ia,ib,ic,id
      integer ncryst
      integer valence(mxatm)
      integer ncon(mxatm),nval(mxatm)
      integer icon(mxatm,mxcon),cel(mxatm,mxcon)
      integer nbnd,nang,ntor
      integer ibond(2000,2),bcel(2000,2)
      integer iangl(2000,3),acel(2000,3)
      integer itor(2000,4),tcel(2000,4)
      integer indcel,jtrans
      integer indx
      integer celind(2000,8)
      integer nprim,nxx,nxy,nyy,nxxyy,nxxxy,nyyxy,nxxyyxy
      double precision radii(mxatm)
      double precision dist, rcut
      double precision x(mxatm),y(mxatm),z(mxatm)
      character*4 name(mxatm)
      logical primitive,similar
      logical xx,yy,xy
      logical skip,skip1,skip2,skip3,skip4
c
c     read the unit cell coordinates                                 
c
      open(3,file='input.xyz',status='old')

      read(3,*) n
      do i = 1, n
         read(3,*) name(i),x(i),y(i),z(i)
      end do
c
c     create a supercell by replicating the unit cell along X, Y and XY   PEML
c
      nunit = n      ! NUNIT is the number of atoms in the unit cell   PEML

      call crystal(n,name,valence,x,y,z)   ! CRYSTAL replicates the unit cell along X, Y and XY   PEML 

      ncryst = n

      n = nunit

c      do i = 1, ncryst
c         write(6,*) name(i),x(i),y(i),z(i)
c      end do

c
c     determine valences of atoms 
c
      do i = 1, n     
         if ((name(i)(1:1) .eq. 'S') .or. (name(i)(1:1) .eq. 's')) 
     $             valence(i) = 4
         if ((name(i)(1:1) .eq. 'O') .or. (name(i)(1:1) .eq. 'o')) 
     $             valence(i) = 2
         if ((name(i)(1:1) .eq. 'H') .or. (name(i)(1:1) .eq. 'h'))
     $             valence(i) = 1
      end do
c
      do i = 1, ncryst
         if ((name(i)(1:1) .eq. 'S') .or. (name(i)(1:1) .eq. 's'))
     $             radii(i) = 1.09d0
         if ((name(i)(1:1) .eq. 'O') .or. (name(i)(1:1) .eq. 'o'))
     $             radii(i) = 0.66d0
         if ((name(i)(1:1) .eq. 'H') .or. (name(i)(1:1) .eq. 'h'))
     $             radii(i) = 0.30d0
      end do
c
c
      do i = 1, mxatm 
         ncon(i) = 0
         nval(i) = 0
         do j = 1, mxcon
            icon(i,j) = 0
         end do
      end do
c         
      primitive = .false.
      similar = .false.
      do i = 1, n-1
         do j = i+1, n
            rcut = 1.35d0*(radii(i)+radii(j))
            dist =sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2) 
            if (dist .le. rcut) then
               ncon(i) = ncon(i) + 1
               nval(i) = nval(i) + 1
               ncon(j) = ncon(j) + 1
               nval(j) = nval(j) + 1
               icon(i,ncon(i)) = j
               cel(i,ncon(i)) = 1
               icon(j,ncon(j)) = i
               cel(j,ncon(j)) = 1
            end if
         end do
      end do
c      
      do i = 1, n   
         do j = n+1, ncryst
            rcut = 1.35d0*(radii(i)+radii(j))
            dist =sqrt((x(i)-x(j))**2 + 
     $                 (y(i)-y(j))**2 + (z(i)-z(j))**2)
            if (dist .le. rcut) then
               indcel = int(j/n) + 1
               jtrans = j - (indcel-1)*n                     
               ncon(i) = ncon(i) + 1
               nval(i) = nval(i) + 1
c               ncon(jtrans) = ncon(jtrans) + 1
c               nval(jtrans) = nval(jtrans) + 1
               icon(i,ncon(i)) = jtrans
               cel(i,ncon(i)) = indcel
c               icon(jtrans,ncon(jtrans)) = i
c               cel(jtrans,ncon(jtrans)) = 1
            end if
         end do
      end do
c
      do i = 1, 2000
         ibond(i,1) = 0
         ibond(i,2) = 0
         iangl(i,1) = 0
         iangl(i,2) = 0
         iangl(i,3) = 0
         itor(i,1) = 0
         itor(i,2) = 0
         itor(i,3) = 0
         itor(i,4) = 0
      end do
c         
      nbnd = 0
      do i = 1, n
         ia = i
         if (ncon(ia) .gt. 0) then
            do j = 1, ncon(ia)
               ib = icon(ia,j)
               nbnd = nbnd + 1
               ibond(nbnd,1) = ia
               ibond(nbnd,2) = ib
               bcel(nbnd,1) = 1
               bcel(nbnd,2) = cel(ia,j)
               
               if (nbnd .ge. 2) then
                  skip = .false.
                  skip1 = .false.
                  skip2 = .false.
                  do k = 1, nbnd-1
                     if ((ibond(k,1) .eq. ibond(nbnd,2)) .and.
     $                  (bcel(k,1) .eq. bcel(nbnd,2))) skip1 = .true.
                     if ((ibond(k,2) .eq. ibond(nbnd,1)) .and.
     $                  (bcel(k,2) .eq. bcel(nbnd,1))) skip2 = .true.
                     if (skip1 .and. skip2) skip = .true.
                     if (skip) then
                        nbnd = nbnd - 1
                        skip = .false.
                        skip1 = .false.
                        skip2 = .false.
                        goto 550
                     end if
                  end do
                  skip = .false.
                  skip1 = .false.
                  skip2 = .false.
                  do k = 1, nbnd-1
                     if ((ibond(k,1) .eq. ibond(nbnd,1)) .and.
     $                  (bcel(k,1) .eq. bcel(nbnd,1))) skip1 = .true.
                     if ((ibond(k,2) .eq. ibond(nbnd,2)) .and.
     $                  (bcel(k,2) .eq. bcel(nbnd,2))) skip2 = .true.
                     if (skip1 .and. skip2) skip = .true.
                     if (skip) then
                        nbnd = nbnd - 1
                        skip = .false.
                        skip1 = .false.
                        skip2 = .false.
                        goto 550
                     end if
                  end do
  550             continue
               end if
            end do
         end if
      end do
c
      nang = 0
      do i = 1, nbnd
         ia = ibond(i,1)
         ib = ibond(i,2)
         if (ncon(ib) .gt. 0) then
            do j = 1, ncon(ib)
               ic = icon(ib,j)
               if (ic .ne. ia) then
                  nang = nang + 1
                  iangl(nang,1) = ia
                  iangl(nang,2) = ib
                  iangl(nang,3) = ic
                  acel(nang,1) = bcel(i,1)
                  acel(nang,2) = bcel(i,2)
                  acel(nang,3) = bcel(i,2) + cel(ib,j) - 1
               end if
            end do
         end if
c
         if (ncon(ia) .gt. 0) then
            do j = 1, ncon(ia)
               ic = icon(ia,j)
               if (ic .ne. ib) then
                  nang = nang + 1
                  iangl(nang,1) = ic
                  iangl(nang,2) = ia
                  iangl(nang,3) = ib
                  acel(nang,1) = bcel(i,1) + cel(ia,j) - 1
                  acel(nang,2) = bcel(i,1)                
                  acel(nang,3) = bcel(i,2)
               end if
            end do
         end if
      end do
c
      do i = 1, nang-1 
         do k = i+1, nang
            skip = .false.
            skip1 = .false.
            skip2 = .false.
            skip3 = .false.
            if ((iangl(i,1) .eq. iangl(k,1)) .and.
     $         (acel(i,1) .eq. acel(k,1))) skip1 = .true.
            if ((iangl(i,2) .eq. iangl(k,2)) .and.
     $         (acel(i,2) .eq. acel(k,2))) skip2 = .true.
            if ((iangl(i,3) .eq. iangl(k,3)) .and.
     $         (acel(i,3) .eq. acel(k,3))) skip3 = .true.
            if (skip1 .and. skip2 .and. skip3) skip = .true.
            if (skip) then
               do l = k, nang-1
                  iangl(l,1) = iangl(l+1,1)
                  acel(l,1) = acel(l+1,1)
                  iangl(l,2) = iangl(l+1,2)
                  acel(l,2) = acel(l+1,2)
                  iangl(l,3) = iangl(l+1,3)
                  acel(l,3) = acel(l+1,3)
               end do
               nang = nang - 1
            end if
         end do
      end do
      do i = 1, nang-1
         do k = i+1, nang
            skip = .false.
            skip1 = .false.
            skip2 = .false.
            skip3 = .false.
            if ((iangl(i,1) .eq. iangl(k,3)) .and.
     $         (acel(i,1) .eq. acel(k,3))) skip1 = .true.
            if ((iangl(i,2) .eq. iangl(k,2)) .and.
     $         (acel(i,2) .eq. acel(k,2))) skip2 = .true.
            if ((iangl(i,3) .eq. iangl(k,1)) .and.
     $         (acel(i,3) .eq. acel(k,1))) skip3 = .true.
            if (skip1 .and. skip2 .and. skip3) skip = .true.
            if (skip) then
               do l = k, nang-1
                  iangl(l,1) = iangl(l+1,1)
                  acel(l,1) = acel(l+1,1)
                  iangl(l,2) = iangl(l+1,2)
                  acel(l,2) = acel(l+1,2)
                  iangl(l,3) = iangl(l+1,3)
                  acel(l,3) = acel(l+1,3)
               end do
               nang = nang - 1
            end if
         end do
      end do
c
      ntor = 0
      do i = 1, nang
         ia = iangl(i,1)
         ib = iangl(i,2)
         ic = iangl(i,3)
         if (ncon(ic) .gt. 0) then
            do j = 1, ncon(ic)
               id = icon(ic,j)
               if (id .ne. ib) then
                  ntor = ntor + 1
                  itor(ntor,1) = ia
                  itor(ntor,2) = ib
                  itor(ntor,3) = ic
                  itor(ntor,4) = id
                  tcel(ntor,1) = acel(i,1)
                  tcel(ntor,2) = acel(i,2)
                  tcel(ntor,3) = acel(i,3)
                  tcel(ntor,4) = acel(i,3) + cel(ic,j) - 1
               end if
            end do
         end if
c
         if (ncon(ia) .gt. 0) then
            do j = 1, ncon(ia)
               id = icon(ia,j)
               if (id .ne. ib) then
                  ntor = ntor + 1
                  itor(ntor,1) = id
                  itor(ntor,2) = ia
                  itor(ntor,3) = ib
                  itor(ntor,4) = ic
                  tcel(ntor,1) = acel(i,1) + cel(ia,j) - 1
                  tcel(ntor,2) = acel(i,1)
                  tcel(ntor,3) = acel(i,2)
                  tcel(ntor,4) = acel(i,3) 
               end if
            end do
         end if   
      end do
c
      do i = 1, ntor-1
         do k = i+1, ntor
            skip = .false.
            skip1 = .false.
            skip2 = .false.
            skip3 = .false.
            skip4 = .false.
            if ((itor(i,1) .eq. itor(k,1)) .and.
     $         (tcel(i,1) .eq. tcel(k,1))) skip1 = .true.
            if ((itor(i,2) .eq. itor(k,2)) .and.
     $         (tcel(i,2) .eq. tcel(k,2))) skip2 = .true.
            if ((itor(i,3) .eq. itor(k,3)) .and.
     $         (tcel(i,3) .eq. tcel(k,3))) skip3 = .true.
            if ((itor(i,4) .eq. itor(k,4)) .and.
     $         (tcel(i,4) .eq. tcel(k,4))) skip4 = .true.
            if (skip1 .and. skip2 .and. skip3 .and. skip4) skip = .true.
            if (skip) then
               do l = k, ntor-1
                  itor(l,1) = itor(l+1,1)
                  tcel(l,1) = tcel(l+1,1)
                  itor(l,2) = itor(l+1,2)
                  tcel(l,2) = tcel(l+1,2)
                  itor(l,3) = itor(l+1,3)
                  tcel(l,3) = tcel(l+1,3)
                  itor(l,4) = itor(l+1,4)
                  tcel(l,4) = tcel(l+1,4)
               end do
               ntor = ntor - 1
            end if
         end do
      end do
      do i = 1, ntor-1
         do k = i+1, ntor
            skip = .false.
            skip1 = .false.
            skip2 = .false.
            skip3 = .false.
            skip4 = .false.
            if ((itor(i,1) .eq. itor(k,4)) .and.
     $         (tcel(i,1) .eq. tcel(k,4))) skip1 = .true.
            if ((itor(i,2) .eq. itor(k,3)) .and.
     $         (tcel(i,2) .eq. tcel(k,3))) skip2 = .true.
            if ((itor(i,3) .eq. itor(k,2)) .and.
     $         (tcel(i,3) .eq. tcel(k,2))) skip3 = .true.
            if ((itor(i,4) .eq. itor(k,1)) .and.
     $         (tcel(i,4) .eq. tcel(k,1))) skip4 = .true.
            if (skip1 .and. skip2 .and. skip3 .and. skip4) skip = .true.
            if (skip) then
               do l = k, ntor-1
                  itor(l,1) = itor(l+1,1)
                  tcel(l,1) = tcel(l+1,1)
                  itor(l,2) = itor(l+1,2)
                  tcel(l,2) = tcel(l+1,2)
                  itor(l,3) = itor(l+1,3)
                  tcel(l,3) = tcel(l+1,3)
                  itor(l,4) = itor(l+1,4)
                  tcel(l,4) = tcel(l+1,4)
               end do
               ntor = ntor - 1
            end if
         end do
      end do

c
c     PRINT SECTION
c
c
c     print bonds 
c

      write(6,'(///,a)')'Bonds inside unit cell'
      do i = 1,nbnd
         if (bcel(i,2) .eq. 1) then
            write(6,200) ibond(i,1),ibond(i,2),bcel(i,1),bcel(i,2)
         end if
      end do

      write(6,'(//,a)')'Bonds along XX'
      do i = 1,nbnd
         if (bcel(i,2) .eq. 2) then
            write(6,200)ibond(i,1),ibond(i,2),bcel(i,1),bcel(i,2)
         end if
      end do
      write(6,'(/,a)')'Bonds along YY'
      do i = 1,nbnd
         if (bcel(i,2) .eq. 3) then
            write(6,200)ibond(i,1),ibond(i,2),bcel(i,1),bcel(i,2)
         end if
      end do
      write(6,'(/,a)')'Bonds along XY'
      do i = 1,nbnd
         if (bcel(i,2) .eq. 4) then
            write(6,200)ibond(i,1),ibond(i,2),bcel(i,1),bcel(i,2)
         end if
      end do
 200  format (2i5,4x,2i3)

c      write(6,*)'Nang=',nang

      do i = 1, 2000
         celind(i,1) = 0
         celind(i,2) = 0
         celind(i,3) = 0
         celind(i,4) = 0
         celind(i,5) = 0
         celind(i,6) = 0
         celind(i,7) = 0
         celind(i,8) = 0
      end do
                             
      nprim = 0
      nxx = 0
      nxy = 0
      nyy = 0
      nxxyy = 0
      nxxxy = 0
      nyyxy = 0
      do i = 1, nang
         xx = .false.
         xy = .false.
         yy = .false.
         do j = 1, 3
            if (acel(i,j) .eq. 2) xx = .true.
            if (acel(i,j) .eq. 3) yy = .true.
            if (acel(i,j) .eq. 4) xy = .true.
         end do

         if (.not. xx .and. .not. yy .and. .not. xy) then
            nprim = nprim + 1
            celind(nprim,1) = i
         end if

         if (xx .and. .not. yy .and. .not. xy) then
            nxx = nxx + 1
            celind(nxx,2) = i
         end if

         if (yy .and. .not. xx .and. .not. xy) then
            nyy = nyy + 1
            celind(nyy,3) = i
         end if
c
         if (xy .and. .not. xx .and. .not. yy) then
            nxy = nxy + 1
            celind(nxy,4) = i
         end if
c
         if (xx .and. yy .and. .not. xy) then
            nxxyy = nxxyy + 1
            celind(nxxyy,5) = i
         end if
c
         if (xx .and. .not. yy .and. xy) then
            nxxxy = nxxxy + 1
            celind(nxxxy,6) = i
         end if
c
         if (.not. xx .and.  yy .and. xy) then
            nyyxy = nyyxy + 1
            celind(nyyxy,7) = i
         end if
      end do
c
c     print angles
c

c      write(6,'(///,a)')'Angles inside unit cell'
c      do i = 1, nprim
c         indx = celind(i,1)
c         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3)
c      end do

      write(6,'(//,a)')'Angles along XX'
      do i = 1, nxx
         indx = celind(i,2)
         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3)
      end do
      write(6,'(/,a)')'Angles along YY'
      do i = 1, nyy
         indx = celind(i,3)
         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3)
      end do
      write(6,'(/,a)')'Angles along XY'
      do i = 1, nxy
         indx = celind(i,4)
         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3)
      end do
      write(6,'(/,a)')'Angles along XX and YY'
      do i = 1, nxxyy
         indx = celind(i,5)
         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3)
      end do
      write(6,'(/,a)')'Angles along XX and XY'
      do i = 1, nxxxy
         indx = celind(i,6)
         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3) 
      end do
      write(6,'(/,a)')'Angles along YY and XY'
      do i = 1, nyyxy
         indx = celind(i,7)
         write(6,201)(iangl(indx,k),k=1,3),(acel(indx,k),k=1,3)
      end do

 201  format (3i5,4x,3i3)

      nang = nang - (nprim+nxx+nyy+nxy+nxxyy+nxxxy+nyyxy)

      if (nang .eq. 0) then
         write(6,'(/,a)') 'All angles are accounted for'
      else 
         write(6,202)'Angles missing: ',nang
         write(6,'(/,a)')'You have to change the code, print all angles' 
         write(6,'(a)')'and see which ones are missing'
      end if
 202  format (/,a16,i6)

c      write(6,*)'ntor=',ntor

      do i = 1, 2000
         celind(i,1) = 0
         celind(i,2) = 0
         celind(i,3) = 0
         celind(i,4) = 0
         celind(i,5) = 0
         celind(i,6) = 0
         celind(i,7) = 0
         celind(i,8) = 0
      end do

      nprim = 0
      nxx = 0
      nxy = 0
      nyy = 0
      nxxyy = 0
      nxxxy = 0
      nyyxy = 0
      nxxyyxy = 0
      do i = 1, ntor
         xx = .false.
         yy = .false.
         xy = .false.
         do j = 1, 4
            if (tcel(i,j) .eq. 2) xx = .true.
            if (tcel(i,j) .eq. 3) yy = .true.
            if (tcel(i,j) .eq. 4) xy = .true.
         end do
c
         if (.not. xx .and. .not. yy .and. .not. xy) then
            nprim = nprim + 1
            celind(nprim,1) = i
         end if
c
         if (xx .and. .not. yy .and. .not. xy) then 
            nxx = nxx + 1
            celind(nxx,2) = i
         end if

         if (yy .and. .not. xx .and. .not. xy) then
            nyy = nyy + 1
            celind(nyy,3) = i
         end if
c
         if (xy .and. .not. xx .and. .not. yy) then
            nxy = nxy + 1
            celind(nxy,4) = i
         end if
c
         if (xx .and. yy .and. .not. xy) then
            nxxyy = nxxyy + 1
            celind(nxxyy,5) = i
         end if
c
         if (xx .and. .not. yy .and. xy) then
            nxxxy = nxxxy + 1
            celind(nxxxy,6) = i
         end if
c
         if (.not. xx .and.  yy .and. xy) then
            nyyxy = nyyxy + 1
            celind(nyyxy,7) = i
         end if
c
         if (xx .and.  yy .and. xy) then
            nxxyyxy = nxxyyxy + 1
            celind(nxxyyxy,8) = i
         end if
      end do

c
c     print torsions
c

c      write(6,'(///,a)')'Torsions inside unit cell'
c      do i = 1, nprim
c         indx = celind(i,1)
c         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
c      end do

      write(6,'(//,a)')'Torsions along XX'
      do i = 1, nxx
         indx = celind(i,2)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do
      write(6,'(/,a)')'Torsions along YY'
      do i = 1, nyy
         indx = celind(i,3)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do
      write(6,'(/,a)')'Torsions along XY'
      do i = 1, nxy
         indx = celind(i,4)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do
      write(6,'(/,a)')'Torsions along XX and YY'
      do i = 1, nxxyy
         indx = celind(i,5)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do
      write(6,'(/,a)')'Torsions along XX and XY'
      do i = 1, nxxxy
         indx = celind(i,6)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do
      write(6,'(/,a)')'Torsions along YY and XY'
      do i = 1, nyyxy
         indx = celind(i,7)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do
      write(6,'(/,a)')'Torsions along XX and XY and XY'
      do i = 1, nxxyyxy
         indx = celind(i,8)
         write(6,203)(itor(indx,k),k=1,4),(tcel(indx,k),k=1,4)
      end do

 203  format (4i5,4x,4i3)

      ntor = ntor - (nprim+nxx+nyy+nxy+nxxyy+nxxxy+nyyxy+nxxyyxy)

      if (ntor .eq. 0) then
         write(6,'(//,a)') 'All torsions are accounted for'
      else
         write(6,204)'Torsions missing: ',ntor
         write(6,'(/,a)')'You have to change the code, print all angles'
         write(6,'(a)')'and see which ones are missing'
      end if
 204  format (//,a18,i6)
 
      end 
