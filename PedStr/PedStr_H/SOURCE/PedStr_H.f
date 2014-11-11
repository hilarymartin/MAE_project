      program pedstr
      parameter (np=32768, np1=512)
      include "blk.incl"
      integer*1 sex(np), mmbit
      integer*2 sexp(nt), nkr(np), nnbt(32)
      integer*4 ix(np), id(np), pid(np,2), iz(np), ia(np), ir(np), ns
      integer*4 ipd(np), nn, na, nr(np), nbt(np), iz0(np)
      integer*4 i2(np), ip2(np,2), n2, nbits2, n, pere(nt), mere(nt)
      integer*4 n20(np), i20(np,np1), ip20(np,np1,2), nrr(np,nt), it
      integer*4 nz(np), nmx, nx, inmx
      real*8 rr0(np,np1), rr(nt), ri(np), rf(np)
      logical izp(nt), nd(np)
      character*25 flnm(3)
      common /com12/ id
      common /com123/ pid, iz, nn, ipd
      common /csex/ sex
      common /com123/ ip2, i2, n2, nbits2, mmbit
      common /cprnt/ rr, pere, mere, sexp, izp, n

      open (unit=9, file='pedstr.ini')
      do 30 i = 1, 2
  30  read(9,*) flnm(i)
      read(9,*) mmbit
      close (9)
      open (unit=1, file=flnm(1))
      i = 1
      mix = 0
      na = 0
c      read (1, *) flnm(3)
c  10  read(1,*,end=11) ix(i), id(i), pid(i,1),pid(i,2),sex(i),nn,iz(i)
  10  read(1,*,end=11) ix(i), id(i), pid(i,1), pid(i,2), sex(i), iz(i)
      if (iz(i).ne.0) na = na + 1
      i = i + 1
      if (i.gt.np) stop 'Too much data!'
      goto 10
  11  close (1)
c      print *, 'izmerennyh ', na
c      pause
      nn = i - 1
      na = 0

      open (unit=2, file=flnm(2))

      do 64 i = 1, nn
      iz0(i) = 0
      rf(i) = .0
      nbits2 = -1
      if (iz(i).eq.0) goto 64
      n2 = 0
      print *, 'Proceed ', i, ' of ', nn
      call bcount(i)

      do 12 j = 1, n2
      i20(i,j) = i2(j)
      ip20(i,j,1) = ip2(j,1)
      ip20(i,j,2) = ip2(j,2)
      n20(i) = n2
  12  continue

      do 17 j = 1, n2
      rr(j) = .0
      izp(j) = .false.
      if (iz(i2(j)).ne.0) izp(j) = .true.
      if (ip20(i,j,1).eq.0) then
      pere(j) = 0
      mere(j) = 0
      goto 17
      end if
      do 16 k = 1, n2
      if (ip20(i,j,1).eq.i2(k)) pere(j) = k
      if (ip20(i,j,2).eq.i2(k)) mere(j) = k
  16  continue
  17  continue
      n = n2
      call parente

      do 32 j = 1, n2
      rr0(i,j) = rr(j)
      rf(i) = rf(i) + rr(j)
      if (rr(j).gt.ri(i2(j))) then
      ri(i2(j)) = rr(j)
      nr(i2(j)) = i
      end if
  32  continue
      if (rf(i).eq.0) n20(i) = 0
      if (rf(i).ne.0) nbt(i) = nbits2
  64  continue

ccccc search for absolute equal shapes
      do 13 ii1 = 1, nn
      if (n20(ii1).eq.0) goto 13
      do 14 ii2 = 1, nn
      if (ii1.eq.ii2) goto 14
      if (n20(ii2).eq.0) goto 14
      if ((nbt(ii1).ne.nbt(ii2)).or.(rf(ii1).ne.rf(ii2))) goto 14
      do 15 jj1 = 1, n20(ii1)
      do 48 jj2 = 1, n20(ii2)
      if (i20(ii1,jj1).eq.i20(ii2,jj2)) goto 15
  48  continue
      goto 14
  15  continue
      do 38 jj1 = 1, n20(ii1)
      if (nr(i20(ii1,jj1)).eq.ii1) then
      nr(i20(ii1,jj1)) = ii2
      end if
  38  continue
      n20(ii1) = 0
  14  continue
  13  continue

c      do 34 ii1 = 1, nn
c      if (n20(ii1).eq.0) goto 34
c      do 35 ii2 = 1, nn
c      if (n20(ii2).eq.0) goto 35
c      if (ii1.eq.ii2) goto 35
c      do 36 jj1 = 1, n20(ii1)
c      if (iz(i20(ii1,jj1)).eq.0) goto 36
c      do 37 jj2 = 1, n20(ii2)
c      if (i20(ii1,jj1).eq.i20(ii2,jj2)) goto 36
c  37  continue
c      goto 35
c  36  continue
c      do 38 jj1 = 1, n20(ii1)
c      if (nr(i20(ii1,jj1)).eq.ii1) then
c      print *, i20(ii1,jj1), rr0(ii1,jj1), rr0(ii2,jj2)
c      nr(i20(ii1,jj1)) = ii2
c      end if
c  38  continue
c      n20(ii1) = 0
c      goto 34
c  35  continue
c  34  continue

      nmx = 0
      do 19 i = 1, nn
      nz(i) = 0
      nkr(i) = 0
      nd(i) = .false.
      do 25 j = 1, n20(i)
      it = i20(i,j)
      if (iz(it).eq.0) then
      nz(i) = nz(i) + 1
      goto 25
      end if
      if (rr0(i,j).eq.ri(it)) then
      if (nrr(i,j).eq.2) goto 27
      nrr(i,j) = 1
      do 26 i1 = i+1, nn
      do 50 j1 = 1, n20(i1)
      if ((it.eq.i20(i1,j1)).and.(rr0(i1,j1).eq.ri(it))) then
      nrr(i1,j1) = 2
      nrr(i,j) = 2
      end if
  50  continue
  26  continue
  27  continue
      nkr(i) = nkr(i) + 1
      end if
  25  continue
      if (nkr(i).gt.nmx) nmx = nkr(i)
  19  continue

      do 21 i = 1, nmx
      do 22 j1 = 1, nn
      if (nkr(j1).ne.i) goto 22
      do 29 k1 = 1, n20(j1)
      if (nrr(j1,k1).ne.2) goto 29
      do 31 j2 = 1, nn
      if (j1.eq.j2) goto 31
      do 33 k2 = 1, n20(j2)
      if ((nrr(j2,k2).eq.2).and.(i20(j1,k1).eq.i20(j2,k2))) then
      nrr(j1,k1) = 0
      nkr(j1) = nkr(j1) - 1
      nr(i20(j1,k1)) = j2
      do 39 j3 = 1, nn
      if ((j3.eq.j1).or.(j3.eq.j2)) goto 39
      do 45 k3 = 1, n20(j3)
      if ((nrr(j3,k3).eq.2).and.(i20(j1,k1).eq.i20(j3,k3))) goto 29
  45  continue
  39  continue
      nrr(j2,k2) = 1
      goto 29
      end if
  33  continue
  31  continue
  29  continue
  22  continue
  21  continue

      do 40 i = 0, nmx-1
      nx = nmx - i
      do 41 j1 = 1, nn
      if (nkr(j1).ne.nx) goto 41
      do 42 k1 = 1, n20(j1)
      if (nrr(j1,k1).eq.2) then
      do 43 j2 = 1, nn
      if (j1.eq.j2) goto 43
      do 44 k2 = 1, n20(j2)
      if ((nrr(j2,k2).eq.2).and.(i20(j1,k1).eq.i20(j2,k2))) then

      if (nkr(j1).gt.nkr(j2)) then
      nrr(j2,k2) = 0
      nkr(j2) = nkr(j2) - 1
      nr(i20(j2,k2)) = j1
      na = 1
      if (nkr(j2).eq.0) n20(j2) = 0
      goto 43
      end if

      if (nkr(j1).eq.nkr(j2)) then
      if (nz(j1).lt.nz(j2)) then
      nrr(j2,k2) = 0
      nkr(j2) = nkr(j2) - 1
      nr(i20(j2,k2)) = j1
      na = 1
      if (nkr(j2).eq.0) n20(j2) = 0
      goto 43
      else
      nrr(j1,k1) = 0
      nkr(j1) = nkr(j1) - 1
      nr(i20(j1,k1)) = j2
      na = 1
      if (nkr(j1).eq.0) then
      n20(j1) = 0
      goto 41
      end if
      goto 42
      end if
      end if

      end if
  44  continue
  43  continue
      end if
  42  continue
  41  continue
  40  continue

      nmx = 0
      inmx = 0
      do 28 i = 1, nn
      if (nkr(i).gt.nmx) then
      nmx = nkr(i)
      inmx = i
      end if
  28  continue
c      if (inmx.eq.0) goto 666
      i = inmx
c      do 18 j = 1, n20(i)
c      write(2,*) i20(i,j), ip20(i,j,1), ip20(i,j,2),
c     * sex(i20(i,j)), iz(i20(i,j)), rr0(i,j), nrr(i,j)
c      iz(i20(i,j)) = 0
c      pid(i20(i,j),1) = 0
c      pid(i20(i,j),2) = 0
c  18  continue
c      write(2,*) 'next'
c      goto 001

      ns = 0
      do 47 i = 1, mmbit
      nnbt(i) = 0
  47  continue
c      open (unit=1, file='r.out')
      do 49 i = 1, nn
      if (iz(i).eq.0) goto 49
c      write (1,*) i, nkr(i), nr(i), ri(i), nd(i)
c      write (1,*) i, nr(i)
cc      write (1,*) i,nkr(i),nr(i),ri(i),rf(i),nbt(i),ndb(i,1)
cc      if(ndb(i,1).gt.1)write (1,*)i,ndb(i,1)-1,(ndb(i,j),j=2,ndb(i,1))
cc      if (nkr(i).ne.0) nnbt(nbt(i)) = nnbt(nbt(i)) + 1
cc      print *, i, nkr(i), nr(i), ri(i), rf(i), nbt(i), nnbt(nbt(i))
      if (nkr(i).gt.0) then
      ns = ns + 1
      do 20 j = 1, n20(i)
      iz0(i20(i,j)) = 1
c      write(2,*) ns, i20(i,j), ip20(i,j,1), ip20(i,j,2),
c     * sex(i20(i,j)), iz(i20(i,j)), rr0(i,j), nrr(i,j)
c      write(2,*) i20(i,1), i20(i,j), ip20(i,j,1), ip20(i,j,2),
c     * sex(i20(i,j)), iz(i20(i,j)), rr0(i,j), nrr(i,j)
      write(2,*) ns, i20(i,j), ip20(i,j,1), ip20(i,j,2),
     * sex(i20(i,j)), iz(i20(i,j))
  20  continue
      end if
  49  continue
      close (2)
c      close (1)
      
      do 70 i = 1, nn
      if ((iz(i).ne.0).and.(iz0(i).eq.0)) print *,i,' is not included!'
  70  continue
      
      stop
      end


      subroutine bcount(i)
      parameter (np=32768, tnp=512)
      integer*1 sex(np), tsex, gnr(tnp), ign(2), mmbit
      integer*4 id(np), pid(np,2), iz(np), ipd(np)
      integer*4 tid(tnp), ti, tpid(tnp,2)
      integer*4 i2(np), ip2(np,2), n2, nbits2, nbits1
      common /com12/ id
      common /com123/ pid, iz, nn, ipd
      common /csex/ sex
      common /bcnt/ tid, tpid, gnr, ti, nbits1
      common /com123/ ip2, i2, n2, nbits2, mmbit

c      print *, i
      tid(1) = i
      ti = 1
      gnr(1) = 0
      if (sex(i).eq.1) tsex = 2
      if (sex(i).eq.2) tsex = 1
      tpid(1,1) = pid(i,1)
      tpid(1,2) = pid(i,2)

c dobavljaem roditelej
      if (pid(i,1).ne.0) then
      do 65 j = 1, 2
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = pid(i,j)
      gnr(ti) = 1
      tpid(ti,1) = 0
      tpid(ti,2) = 0
  65  continue
      end if
      if (ti.gt.1) call zachist(i)

c nabiraem potomkov
      do 24 j = 1, nn
      if (pid(j,sex(i)).eq.i) then
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = j
      gnr(ti) = -1
      tpid(ti,1) = pid(j,1)
      tpid(ti,2) = pid(j,2)
      do 28 k = 1, ti
      if (pid(j,tsex).eq.tid(k)) goto 29
  28  continue
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = pid(j,tsex)
      gnr(ti) = 0
      tpid(ti,1) = 0
      tpid(ti,2) = 0
  29  continue
      call zachist(i)
      if (nbits1.ge.mmbit) return
      end if
  24  continue
c      print *, ti
c      do 99 j = 1, ti
c      print *, 'nachalo', gnr(j), j, tid(j), tpid(j,1), tpid(j,2)
c  99  continue
c      call tcheck

c razshirenie
      do 54 j = 1, nn
c      print *, j
c      pause
      if (j.gt.ti) return
      if (gnr(j).eq.0) goto 54
c dlja predkov
      if (gnr(j).gt.0) then
c      print *, 'dlja predkov'

c dobavljaem potomkov
      do 67 k = 1, nn
      if (pid(k,sex(tid(j))).eq.tid(j)) then
      do 68 l = 1, ti
      if (id(k).eq.tid(l)) goto 67
  68  continue
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = id(k)
      gnr(ti) = gnr(j) + 1
      tpid(ti,1) = pid(k,1)
      tpid(ti,2) = pid(k,2)
c      print *,'potomka predkam',gnr(ti),ti,tid(ti),tpid(ti,1),tpid(ti,2)

      if (sex(tid(j)).eq.1) tsex = 2
      if (sex(tid(j)).eq.2) tsex = 1
      do 69 l = 1, ti
      if (pid(k,tsex).eq.tid(l)) goto 70
  69  continue
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = pid(k,tsex)
      gnr(ti) = 0
      tpid(ti,1) = 0
      tpid(ti,2) = 0
  70  continue
c      call tcheck
      call zachist(i)
      if (nbits1.ge.mmbit) return
      end if
  67  continue

c dobavljaem roditelej
      if (pid(tid(j),1).ne.0) then
      do 66 k = 1, 2
      do 71 l = 1, ti
c      if (tpid(j,k).eq.tid(l)) goto 66
      if (pid(tid(j),k).eq.tid(l)) then
      tpid(j,k) = pid(tid(j),k)
      goto 66
      end if
  71  continue
      ti = ti + 1
      if (ti.gt.tnp) return
      tpid(j,k) = pid(tid(j),k)
      tid(ti) = tpid(j,k)
      gnr(ti) = gnr(j) + 1
      tpid(ti,1) = 0
      tpid(ti,2) = 0
  66  continue
c      print *,'roditel predkam',gnr(ti),ti,tid(ti),tpid(ti,1),tpid(ti,2)
c      call tcheck
      call zachist(i)
      if (nbits1.ge.mmbit) return
      end if
      end if

c dlja potomkov
      if (gnr(j).lt.0) then
c      print *, 'dlja potomkov'
c dobavljaem potomkov
      do 72 k = 1, nn
      if (pid(k,sex(tid(j))).eq.tid(j)) then
      do 73 l = 1, ti
      if (id(k).eq.tid(l)) goto 72
  73  continue
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = id(k)
      gnr(ti) = gnr(j) - 1
      tpid(ti,1) = pid(k,1)
      tpid(ti,2) = pid(k,2)
c      print *,'potomka potomku',gnr(ti),ti,tid(ti),tpid(ti,1),tpid(ti,2)

      if (sex(tid(j)).eq.1) tsex = 2
      if (sex(tid(j)).eq.2) tsex = 1
      do 74 l = 1, ti
      if (pid(k,tsex).eq.tid(l)) goto 75
  74  continue
      ti = ti + 1
      if (ti.gt.tnp) return
      tid(ti) = pid(k,tsex)
      gnr(ti) = 0
      tpid(ti,1) = 0
      tpid(ti,2) = 0
  75  continue
c      call tcheck
      call zachist(i)
      if (nbits1.ge.mmbit) return

      end if
  72  continue
      end if
  54  continue

c      print *, ti
  
      return
      end


      subroutine zachist(ii)
      parameter (np=32768, tnp=512)
      integer*1 sex(np), tsex, gnr(tnp), ign(2), mmbit
      integer*4 id(np), pid(np,2), iz(np), ipd(np)
      integer*4 tid(tnp), ti, tpid(tnp,2), tp(2)
      integer*4 i1(np), ip1(np,2), i2(np), ip2(np,2), n1, n2
      integer*4 nbits1, nbits2
      common /com12/ id
      common /com123/ pid, iz, nn, ipd
      common /csex/ sex
      common /bcnt/ tid, tpid, gnr, ti, nbits1
      common /com123/ ip2, i2, n2, nbits2, mmbit
      common /pchk/ ip1, i1, n1

      nbits1 = 0
      n1 = ti
      do 14 i = 1, ti
c      if (ii.eq.1482) print *, i, ti
      i1(i) = tid(i)
      ip1(i,1) = tpid(i,1)
      ip1(i,2) = tpid(i,2)
c      print *, ii, i, ti, i1(i), ip1(i,1), ip1(i,2)
  14  continue

cc Zachistka ot neobsledovannyh
  48  continue
      if (n1.eq.0) return
      isk = 0
      do 27 i = 1, n1
      if (i1(i).eq.0) then
      n1 = n1 - 1
      do 28 j = 1, n1
      i1(j) = i1(j+1)
      ip1(j,1) = ip1(j+1,1)
      ip1(j,2) = ip1(j+1,2)
  28  continue
      goto 48
      end if
  27  continue
c      print *, 'Next!'
c      call pcheck

      tp(1) = 0
      tp(2) = 0
cc ot potomkov
      do 41 j = 1, n1
c      if ((ip1(j,1).eq.0).or.(ip1(j,2).eq.0)) then
c      if (ip1(j,1).ne.ip1(j,2)) then
c      print *, 'oshibka', i1(j), ip1(j,1), ip1(j,2)
c      pause
c      end if
c      end if
      if (iz(i1(j)).ne.0) goto 41
      do 42 k = 1, n1
      if ((ip1(k,1).eq.i1(j)).or.(ip1(k,2).eq.i1(j))) goto 41
  42  continue
      isk = j
      n1 = n1 - 1
      tp(1) = ip1(j,1)
      tp(2) = ip1(j,2)
      do 43 k = j, n1
      i1(k) = i1(k+1)
      ip1(k,1) = ip1(k+1,1)
      ip1(k,2) = ip1(k+1,2)
  43  continue
      i1(n1+1) = 0
      ip1(n1+1,1) = 0
      ip1(n1+1,2) = 0
      goto 33
  41  continue
      goto 34
c otvalivshiesja roditeli
  33  continue
      do 15 i = 1, 2
      do 16 j = 1, n1
      if (ip1(j,i).eq.tp(i)) goto 15
  16  continue
      do 17 j = 1, n1
      if (i1(j).eq.tp(i)) then
      if (ip1(j,1).eq.0) then
      n1 = n1 - 1
      do 18 k = j, n1
      i1(k) = i1(k+1)
      ip1(k,1) = ip1(k+1,1)
      ip1(k,2) = ip1(k+1,2)
  18  continue
      i1(n1+1) = 0
      ip1(n1+1,1) = 0
      ip1(n1+1,2) = 0
      end if
      goto 15
      end if
  17  continue
  15  continue
      if (isk.ne.0) then
c      print *, 'ot potomkov'
c      call pcheck
      end if
      if (isk.ne.0) goto 48
  34  continue
c      print *, 'pered roditeljami'
c      call pcheck
cc ot roditelej
      do 19 i = 1, n1
      if (ip1(i,1).eq.0) goto 19
      tp(1) = ip1(i,1)
      tp(2) = ip1(i,2)
      if ((iz(tp(1)).ne.0).or.(iz(tp(2)).ne.0)) goto 19
      do 20 j = 1, n1
      if (i.eq.j) goto 20
      if ((ip1(j,1).eq.tp(1)).or.(ip1(j,2).eq.tp(2))) goto 19
  20  continue
      do 21 j = 1, 2
      do 22 k = 1, n1
      if (i1(k).eq.tp(j)) then
      if (ip1(k,1).ne.0) then
      tp(j) = 0
      else
      tp(j) = k
      end if
      goto 21
      end if
  22  continue
  21  continue
      if ((tp(1).eq.0).or.(tp(2).eq.0)) goto 19
      if (tp(2).gt.tp(1)) tp(2) = tp(2) - 1
      ip1(i,1) = 0
      ip1(i,2) = 0
      do 23 j = 1, 2
      do 25 k = tp(j), n1 
      i1(k) = i1(k+1)
      ip1(k,1) = ip1(k+1,1)
      ip1(k,2) = ip1(k+1,2)
  25  continue
      i1(n1+1) = 0
      ip1(n1+1,1) = 0
      ip1(n1+1,2) = 0
  23  continue
      isk = 1
      n1 = n1 - 2
c      print *, 'ot roditelej'
c      call pcheck
      goto 48
  19  continue
c      if (isk.ne.0) goto 48

c bits counter

c      call pcheck

      do 26 i = 1, n1
c      print *, i1(i), ip1(i,1), ip1(i,2)
      if (ip1(i,1).eq.0) nbits1 = nbits1 - 1
      if (ip1(i,1).ne.0) nbits1 = nbits1 + 2
  26  continue

      if (nbits1.le.mmbit) then
      if (nbits1.gt.nbits2)then
      nbits2 = nbits1
      n2 = n1
      do 31 i = 1, n2
      i2(i) = i1(i)
      ip2(i,1) = ip1(i,1)
      ip2(i,2) = ip1(i,2)
  31  continue
      end if
      end if
  
c      print *, 'Bity = ', nbits1, ' na ', n1

      return
      end


c      subroutine pcheck
c      parameter (np=32768, tnp=1024)
c      integer*4 i1(np), ip1(np,2), n1
c      common /pchk/ ip1, i1, n1

c      do 36 i = 1, n1
c      print *, 'pcheck:', i, i1(i), ip1(i,1), ip1(i,2)
c  36  continue

c      do 10 i = 1, n1
c      do 11 j = 1, n1
c      if (i.eq.j) goto 11
c      if (i1(i).eq.i1(j)) then
c      print *, 'dubl!'
c      pause
c      end if
c  11  continue
c  10  continue

c      do 33 i = 1, n1
c      if ((ip1(i,1).eq.0).or.(ip1(i,2).eq.0)) then
c      if (ip1(i,1).ne.ip1(i,2)) then
c      print *, 'Parents error!', i1(i), ip1(i,1), ip1(i,2)
c      pause
c      end if
c      end if
c      do 34 j = 1, 2
c      if (ip1(i,j).ne.0) then
c      do 35 k = 1, n1
c      if (i1(k).eq.ip1(i,j)) goto 34
c  35  continue
c      print *, 'Error!:', i1(i), ip1(i,j)
c      pause
c      end if
c  34  continue
c  33  continue
c      return
c      end

c      subroutine tcheck
c      parameter (np=32768, tnp=1024)
c      integer*1 gnr(tnp)
c      integer*4 tid(tnp), ti, tpid(tnp,2), nbits1
c      common /bcnt/ tid, tpid, gnr, ti, nbits1

c      do 36 i = 1, ti
c      print *, 'tcheck:', i, tid(i), tpid(i,1), tpid(i,2)
c  36  continue

c      do 10 i = 1, ti
c      do 11 j = 1, ti
c      if (i.eq.j) goto 11
c      if (tid(i).eq.tid(j)) then
c      print *, 'dubl!', i, ti, tid(i)
c      pause
c      end if
c   11 continue
c   10 continue

c      do 33 i = 1, ti
c      if ((tpid(i,1).eq.0).or.(tpid(i,2).eq.0)) then
c      if (tpid(i,1).ne.tpid(i,2)) then
c      print *, 'Parents error!', tid(i), tpid(i,1), tpid(i,2)
c      pause
c      end if
c      end if
c      do 34 j = 1, 2
c      if (tpid(i,j).ne.0) then
c      do 35 k = 1, ti
c      if (tid(k).eq.tpid(i,j)) goto 34
c  35  continue
c      print *, 'Error!:', tid(i), tpid(i,j)
c      pause
c      end if
c  34  continue
c  33  continue
c      return
c      end


      subroutine parente
      implicit none
      include "blk.incl"
      integer*2 sex(nt)
      integer*4 i, j, k, n, ii, jj, kk, is, il, ne, na, ie, ip, im
      integer*4 je, nd, iad, iaf, ifail, ig, ig1, ig2, eff
      integer*4 lliste(nlistx,2), nlist(2), iopt(4)
      integer*4 pere(nt), mere(nt), ped(2,nt), point(nt)
      integer*4 ord(nt), rord(nt), ic(nt), ndes(nt)
      real*8 statis(106), x, f(0:nt), l(nt), d(nt), rr(nt)
      logical ts, ts2(2,2), chy, izp(nt)
      character*128 sti, stc, jour
      common /cprnt/ rr, pere, mere, sex, izp, n
      include "format.incl"

      ts=.true.
c      if (stc.eq.'no') ts=.false.
      do i=1, 2
        do j=1, 2
          ts2(i,j)=ts
        end do
      end do
      nlist(1)=0
      nlist(2)=0
	
      do 1 i = 1, n
      if (izp(i)) then
        nlist(1) = nlist(1) + 1
c	print *, n, nlist(1)
        lliste(nlist(1),1) = i
c	print *, n, lliste(nlist(1),1)
      end if
  1   continue

      do ig=1, 2
       do i=1, nlist(ig)
        if (lliste(i,ig).gt.n) then
          print 8512,lliste(i,ig)
          stop
          end if
        end do
       end do

c *** renumerotation du plus vieux au plus jeune
      call comp_d (n, pere, mere, ped, ord, rord)

      do ig=1, 2
       do i=1, nlist(ig)
        ndes(ord(lliste(i,ig)))=1
        end do
       end do

      do i=1,n
        point(i)=0
        l(i)=0.
        d(i)=0.
        end do
      call meuw(n, ped, f, d, l, point,ndes)

      if (ts) then
cTolya        open (3,file=stc,form='formatted')
        do ig=1, 2
          do i=1, nlist(ig)
            j=lliste(i,ig)
      if (ts2(ig,ig)) then
c      print 8888, j, j, f(ord(j)), ig, ig
      rr(j) = f(ord(j))
      end if
            end do
          end do
        end if
 8888 format (2i8,e14.6,2i2)

      n=n + 1
      do i=1, n
        point(i)=0
        l(i)=0.
        end do

      do ig1=1, 2
       do ig2=ig1, 2
        call statt(x,statis,0)
c        print *,'Groups ',ig1,' x ',ig2
        do i=1,nlist(ig1)
          ii=1
          if (ig1.eq.ig2) ii=i+1
          do j=ii, nlist(ig2)
           ip=ord(lliste(i,ig1))
           im=ord(lliste(j,ig2))
           ped(1,n)=ip
           ped(2,n)=im
           call consang(n,ped,f,d,l,point)
      if (ts2(ig1,ig2)) then
c      print 8888, lliste(i,ig1), lliste(j,ig2), f(n), ig1, ig2
      rr(lliste(i,ig1)) = rr(lliste(i,ig1)) + f(n)
      rr(lliste(j,ig2)) = rr(lliste(j,ig2)) + f(n)
      end if
           call statt(f(n),statis,1)
           end do
         end do
        call statt(x,statis,2)
        end do
       end do
cTolya       if (ts) close(3)
      end

      subroutine statt(x,statis,io)
      implicit none
      integer*4 io, nm,i,ii,j,jj,k,l
      real*8  statis(*),x
      character*1 etoile

      include "format.incl"

      if (io.eq.0) then
        do i=1, 106
          statis(i)=0.d0
          end do
        return

      else if (io.eq.1) then
        statis(3)=statis(3) + 1.d0
        statis(4)=statis(4) + x
        statis(5)=statis(5) +x*x
        k=6+int(100.d0 * x)
        statis(k)=statis(k)+1.d0
        return

      else if (io.eq.2) then
        etoile='*'
        if (statis(3).gt.0.d0) then
          statis(4)=statis(4)/statis(3)
          statis(5)=statis(5) - statis(3)*statis(4)*statis(4)
          statis(5)=dsqrt(statis(5)/statis(3))
          end if
c        print 8517,statis(3),statis(4),statis(5)
        do i=1, 100
          k=int(.5 + 100.*statis(5+i)/statis(3))
c          if (statis(5+i).gt.0)
c     *    print 8518, i-1, i, statis(5+i), k, (etoile,l=1,k)
          end do
        end if
      return
      end

      subroutine comp_d (n,sire, dam,ped,ord,rord)
      implicit none
      integer*4 n,sire(*),dam(*),ped(2,*),nbit,k,i,j,ks,kd
      integer*4 ord(*), rord(*)

      include "format.incl"

	nbit=0
	do i=1, n
	  ord(i)=0
	  end do
	k=0
	do while (k.lt.n .and. nbit.le.20)
	  nbit=nbit + 1
	  do i=1, n
	    if (ord(i).eq.0) then
		if (sire(i).le.0 .or. ord(sire(i)).ne.0) then
	   	 if (dam(i) .le.0 .or. ord( dam(i)).ne.0) then
		  k=k+1
		  ord(i)=k
                rord(k)=i
	         end if
	        end if
	    end if
	  end do
	end do

c  Test qu'il n'y a pas de boucle dans le pedigree
      j=0
      if (k.ne.n) then
        do i=1,n
          if (ord(i).eq.0) then
            j=j+1
            if (j.lt.100) print '(3i10)',i,sire(i),dam(i)
            sire(i)=0
            dam(i)=0
            k=k+1
            ord(i)=k
            rord(k)=i
          end if
       end do
      end if
      if (j.gt.0) print 900,j

      do i=1, n
	  j=rord(i)
	  ped(1,i)=sire(j)
	  ped(2,i)=dam(j)
	  if (ped(1,i).gt.0) ped(1,i)=ord(sire(j))
	  if (ped(2,i).gt.0) ped(2,i)=ord(dam(j))
        end do

      DO i=1,n
        if (i.le.ped(1,i) .or. i.le.ped(2,i)) then
	    print *,'Problem in coding pedigree'
          print *,i, ped(1,i),ped(2,i)
          stop
          end if
        ks=ped(1,i)
        kd=ped(2,i)
        ped(1,i)=max(ks,kd)
        ped(2,i)=min(ks,kd)
        end do
	return
	end


C*** Methode de Meuwissen
      subroutine meuw(n, ped, f, d, l ,point,ndes)
      implicit none
      integer*4 n, ped(2,*), point(*),ndes(*),np,npar
      integer*4 ninbr, i, j,k, ik, is, id, ks, kd
      real*8 f(0:n), d(*), l(*),r, fi
      ninbr=0
      f(0)=-1.d0
      DO i=1,n
        point(i)=0
        end do
      DO i=1,n
        if (ped(1,i).gt.0) ndes(ped(1,i))=ndes(ped(1,i))+1
        if (ped(2,i).gt.0) ndes(ped(2,i))=ndes(ped(2,i))+1
        end do
      npar=0
      do i=1, n
        if (ndes(i).gt.0) npar=npar+1
        end do
      if (npar.eq.0) return
      DO i=1,n
       if (ndes(i).gt.0) then
        is=ped(1,i)
        id=ped(2,i)
        ped(1,i)=max(is,id)
        ped(2,i)=min(is,id)
        d(i)=.5d0 - .25d0*(f(is)+f(id))
        if (is.eq.0.or.id.eq.0) then
           f(i)=0.d0
        else
          np=0
          fi=-1.d0
          l(i)=1.d0
          j=i
          do while(j.ne.0)
            k=j
            r=.5d0 * l(k)
            ks=ped(1,k)
            kd=ped(2,k)
            if (ks.gt.0) then
              l(ks)=l(ks) + r
              do while(point(k).gt.ks)
                k=point(k)
                end do
              if (ks.ne.point(k)) then
                point(ks)=point(k)
                point(k)=ks
                end if
              if (kd.gt.0) then
                l(kd)=l(kd) + r
                do while(point(k).gt.kd)
                  k=point(k)
                  end do
                if (kd.ne.point(k)) then
                  point(kd)=point(k)
                  point(k)=kd
                  end if
                end if
              end if
              fi=fi + l(j)*l(j)*d(j)
              l(j)=0.d0
              k=j
              j=point(j)
              point(k)=0
              np=np+1
              end do
            f(i)=fi
            if (fi.gt.0.000001d0) ninbr=ninbr + 1
            end if
           end if
          end do
       RETURN
       END

C*** Methode de Meuwissen
      subroutine consang(i, ped, f, d, l ,point)
      implicit none
      integer*4 ped(2,*), point(*)
      integer*4 i, j,k, ik, is, id, ks, kd
      real*8 f(0:i), d(*), l(*),r, fi
        is=ped(1,i)
        id=ped(2,i)
        ped(1,i)=max(is,id)
        ped(2,i)=min(is,id)
        d(i)=.5d0 - .25d0*(f(is)+f(id))
        if (is.eq.0.or.id.eq.0) then
           f(i)=0.d0
           return
           end if
        fi=-1.d0
        l(i)=1.d0
        j=i
        do while(j.ne.0)
          k=j
          r=.5d0 * l(k)
          ks=ped(1,k)
          kd=ped(2,k)
          if (ks.gt.0) then
            l(ks)=l(ks) + r
            do while(point(k).gt.ks)
              k=point(k)
              end do
            if (ks.ne.point(k)) then
              point(ks)=point(k)
              point(k)=ks
              end if
            if (kd.gt.0) then
              l(kd)=l(kd) + r
              do while(point(k).gt.kd)
                k=point(k)
                end do
              if (kd.ne.point(k)) then
                point(kd)=point(k)
                point(k)=kd
                end if
              end if
            end if
          fi=fi + l(j)*l(j)*d(j)
          l(j)=0.d0
          k=j
          j=point(j)
          point(k)=0
          end do
        f(i)=fi
       RETURN
       END
