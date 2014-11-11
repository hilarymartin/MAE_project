      program pedstr
      parameter (np=32768, np1=64)
      include "blk.incl"
      integer*1 sex(np), mmbit, mrad
      integer*4 nbst, ibst(np), pid(np,2), iz(np), ns, nn, nr(np)
      integer*4 i2(nt), ip2(nt,2), n2, nbits2, n, pere(ntp), mere(ntp)
      integer*4 n20(np), i20(np,nt), ip20(np,nt,2), nrr(np,nt), it
      integer*4 nz(np), nx, inmx, nbt(np), nnn, m, mt
      real*8 rr0(np,nt), rr(ntp), ri(np), rf(np), rm, rmm, rmn, rnn
      real*8 rt, rtm, rbst
      logical izp(ntp)
      character*25 flnm(3)
      common /com123a/ pid, iz, nn
      common /csex/ sex
      common /cmrad/ mrad
      common /com123b/ ip2, i2, n2, mt, nbits2, mmbit
      common /cprnt/ rr, pere, mere, izp, n

      open (unit=9, file='pedstr.ini')
      do 30 i = 1, 2
  30  read(9,*) flnm(i)
      read(9,*) mmbit
      read(9,*) mrad
      if (mmbit.gt.30) stop 'Too much bits!'
      close (9)
      open (unit=1, file=flnm(1))
      i = 1
      mix = 0
  10  read(1,*,end=11) nn, nn, pid(i,1), pid(i,2), sex(i), iz(i)
      if (nn.ne.i) stop 'ID is not consecutive!'
      i = i + 1
      if (i.gt.np) stop 'Too much data!'
      goto 10
  11  close (1)
      nn = i - 1
      nnn = 0
      open (unit=2, file=flnm(2))
 001  continue
      call zzach()
      nnn = nnn + 1
      m = 0
      do 51 i = 1, nn
      if (iz(i).ne.0) m = m + 1
  51  continue
      print *,'It is the ',nnn,'th shape; stil ',m,' observed persons'
      do 005 i = 1, nt
      pere(i) = 0
      mere(i) = 0
      rr(i) = .0
      i2(i) = 0
      ip2(i,1) = 0
      ip2(i,2) = 0
 005  continue
      do 002 i = 1, np
      nr(i) = 0
      ri(i) = .0
      rf(i) = .0

      do 003 j = 1, nt
      i20(i,j) = 0
      ip20(i,j,1) = 0
      ip20(i,j,2) = 0
      rr0(i,j) = .0
 003  continue
      do 004 j = 1, nt
      nrr(i,j) = 0
 004  continue
 002  continue

      do 64 i = 1, nn
      n20(i) = 0
      nbt(i) = 0
      rf(i) = .0
      nbits2 = -1
      if (iz(i).eq.0) goto 64
      n2 = 0
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

      do 19 i = 1, nn
      nz(i) = 0
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
      end if
  25  continue
  19  continue

      inmx = 0
      rmm = .0
      nbst = 0
      do 28 i = 1, nn
      if (n20(i).eq.0) goto 28
      rm = .0
      do 52 j = 1, n20(i)
      if (iz(i20(i,j)).ne.0) then
      rm = rm + rr0(i,j)
      end if
  52  continue
      if (rm.eq.rmm) then
      nbst = nbst + 1
      ibst(nbst) = i
      end if
      if (rm.gt.rmm) then
      nbst = 1
      ibst(1) = i
      rmm = rm
      end if
  28  continue

      if (nbst.eq.0) goto 666
      if (nbst.eq.1) then
      inmx = ibst(1)
      goto 125
      end if

      rbst = .0
      do 124 i = 1, nbst
      do 120 ii = 1, nn
      rr(ii) = .0
      pere(ii) = pid(ii,1)
      mere(ii) = pid(ii,2)
      izp(ii) = .false.
      if (iz(ii).ne.0) then
      izp(ii) = .true.
      do 122 j = 1, n20(ibst(i))
      if (i20(ibst(i),j).eq.ii) then
      izp(ii) = .false.
      goto 123
      end if
  122 continue
      end if
  123 continue
  120 continue
      n = nn
      call parente
      rmn = .0
      do 121 ii = 1, nn
      rmn = rmn + rr(ii)
  121 continue
      if (rmn.gt.rbst) then
      rbst = rmn
      inmx = ibst(i)
      end if
  124 continue
  125 continue

      i = inmx
      k = 0
      kj = 0
      do 23 j = 1, n20(i)
      if (iz(i20(i,j)).ne.0) then
      k = k + 1
      kj = j
      end if
  23  continue
      if (k.eq.1) then
      iz(i20(i,kj)) = 0
      nnn = nnn - 1
      goto 001
      end if

      do 18 j = 1, n20(i)
c      write(2,*) nnn, i20(i,j), ip20(i,j,1), ip20(i,j,2),
c     * sex(i20(i,j)), iz(i20(i,j)), rr0(i,j)
      write(2,*) nnn, i20(i,j), ip20(i,j,1), ip20(i,j,2),
     * sex(i20(i,j)), iz(i20(i,j))
      if (iz(i20(i,j)).ne.0) iz(i20(i,j)) = 0
  18  continue
      goto 001

      close (2)
 666  continue
      stop
      end


      subroutine zzach()
      parameter (np=32768, tnp=4096)
      integer*4 pid(np,2), iz(np), i, j, k, nn
      logical ft
      common /com123a/ pid, iz, nn
  12  continue
      ft = .false.
      do 10 i = 1, nn
      if (pid(i,1).eq.0) goto 10
      if ((iz(pid(i,1)).ne.0).or.(iz(pid(i,2)).ne.0)) goto 13
      if ((pid(pid(i,1),1).ne.0).or.(pid(pid(i,2),1).ne.0)) goto 13
      do 14 j = 1, nn
      if (i.eq.j) goto 14
      if ((pid(i,1).eq.pid(j,1)).or.(pid(i,2).eq.pid(j,2))) goto 13
  14  continue
      pid(i,1) = 0
      pid(i,2) = 0
      ft = .true.
      goto 10
  13  continue
      if (iz(i).ne.0) goto 10
      do 11 j = 1, nn
      if ((pid(j,1).eq.i).or.(pid(j,2).eq.i)) goto 10
  11  continue
      pid(i,1) = 0
      pid(i,2) = 0
      ft = .true.
  10  continue
      if (ft) goto 12
      return
      end


      subroutine bcount(i)
      parameter (np=32768, tnp=4096)
      include "blk.incl"
      integer*1 sex(np), tsex, gnr(tnp), gn(tnp), mmbit, mrad, izc(tnp)
      integer*1 tgn(tnp), tgnr(tnp)
      integer*4 pid(np,2), iz(np), i, j, k, ii, jj, kk, mt, kt
      integer*4 tid(tnp), ti, tpid(tnp,2), tpd(tnp,2), nof
      integer*4 i2(nt), ip2(nt,2), n2, nbits2, nbits1, n3, mr
      real*8 rt, rf
      logical trfl(tnp,2), ngt(tnp)
      common /com123a/ pid, iz, nn
      common /csex/ sex
      common /cmrad/ mrad
      common /bcnt/ tid, tpid, ti, nbits1
      common /com123b/ ip2, i2, n2, mt, nbits2, mmbit
      common /comadd/ gn, gnr
      common /comrf/ rf
      common /cizc/ izc

      do 19 ki = 1, tnp
      trfl(ki,1) = .false.
      trfl(ki,2) = .false.
      ngt(ki) = .false.
      izc(ki) = 0
  19  continue
      tid(1) = i
      ti = 1
      gnr(ti) = 1
      gn(ti) = 1
      tpid(ti,1) = 0
      tpid(ti,2) = 0
      if (iz(i).eq.0) izc(ti) = 1
      n3 = 0
      mr = 0

      do 10 j = 1, mmbit+2

      mt = 0
      call zachist()
      rt = rf
      mt = ti
      kt = 0
      do 15 k = 1, mt
      tpd(k,1) = tpid(k,1)
      tpd(k,2) = tpid(k,2)
      tgnr(k) = gnr(k)
      tgn(k) = gn(k)
  15  continue

      do 16 k = 1, ti
      if ((gn(k).ne.j).or.(izc(k).gt.mrad)) goto 16
      call addoff(j,k)
      if ((gnr(k).gt.0).and.(pid(tid(k),1).ne.0)
     *.and.(izc(k).le.mrad)) then
      call addpar(j,k)
      end if
  16  continue
      call zachist()
      if (nbits1.eq.mmbit) return
      if (n3.eq.n2) then
      mr = mr + 1
      if (mr.eq.mrad) return
      end if
      if (n2.gt.n3) then
      n3 = n2
      mr = 0
      end if
      if (nbits1.lt.mmbit) goto 10
      ti = mt
      do 17 k = 1, mt
      tpid(k,1) = tpd(k,1)
      tpid(k,2) = tpd(k,2)
      gnr(k) = tgnr(k)
      gn(k) = tgn(k)
  17  continue

  01  continue
      mt = 0
      call zachist()
      rt = rf
      mt = ti
      kt = 0
      do 12 k = 1, mt
      tpd(k,1) = tpid(k,1)
      tpd(k,2) = tpid(k,2)
      tgnr(k) = gnr(k)
      tgn(k) = gn(k)
  12  continue

      do 11 k = 1, ti
      if (gn(k).ne.j) goto 11
      if ((trfl(k,1)).or.(izc(k).gt.mrad)) goto 18
      call addoff(j,k)
      call zachist()
      if (nbits1.eq.mmbit) return
      if (nbits1.gt.mmbit) then
      ngt(k) = .true.
      trfl(k,1) = .true.
      ti = mt

      do 25 ki = 1, mt
      tpid(ki,1) = tpd(ki,1)
      tpid(ki,2) = tpd(ki,2)
      gnr(k) = tgnr(k)
      gn(k) = tgn(k)
  25  continue

      goto 18
      end if
      if (rf.ge.rt) then
      rt = rf
      kt = k
      end if
      ti = mt

      do 13 ki = 1, mt
      tpid(ki,1) = tpd(ki,1)
      tpid(ki,2) = tpd(ki,2)
      gnr(k) = tgnr(k)
      gn(k) = tgn(k)
  13  continue
  18  continue
      if (trfl(k,2)) goto 11

      if ((gnr(k).gt.0).and.(pid(tid(k),1).ne.0)
     *.and.(izc(k).le.mrad)) then
      call addpar(j,k)
      call zachist()
      if (nbits1.eq.mmbit) return
      if (nbits1.gt.mmbit) then
      trfl(k,2) = .true.
      ti = mt
      do 26 ki = 1, mt
      tpid(ki,1) = tpd(ki,1)
      tpid(ki,2) = tpd(ki,2)
      gnr(k) = tgnr(k)
      gn(k) = tgn(k)
  26  continue
      goto 11
      end if
      if (rf.ge.rt) then
      rt = rf
      kt = -k
      end if
      ti = mt
      do 14 ki = 1, mt
      tpid(ki,1) = tpd(ki,1)
      tpid(ki,2) = tpd(ki,2)
      gnr(k) = tgnr(k)
      gn(k) = tgn(k)
  14  continue
      end if
  11  continue
      if (kt.eq.0) return
      if ((kt.gt.0).and.(izc(kt).le.mrad)) then
      call addoff(j,kt)
      call zachist()
      trfl(kt,1) = .true.
      goto 01
      end if
      if ((kt.lt.0).and.(izc(-kt).le.mrad)) then
      call addpar(j,-kt)
      call zachist()
      trfl(-kt,2) = .true.
      goto 01
      end if

      return
  10  continue

      return
      end


      subroutine addpar(j,k)
      parameter (np=32768, tnp=4096)
      include "blk.incl"
      integer*1 sex(np), tsex, gnr(tnp), gn(tnp), mmbit, mrad, izc(tnp)
      integer*4 pid(np,2), iz(np), i, j, k, ii, jj, kk, mt
      integer*4 tid(tnp), ti, tpid(tnp,2)
      integer*4 i2(nt), ip2(nt,2), n2, nbits2, nbits1
      common /com123a/ pid, iz, nn
      common /csex/ sex
      common /bcnt/ tid, tpid, ti, nbits1
      common /com123b/ ip2, i2, n2, mt, nbits2, mmbit
      common /comadd/ gn, gnr
      common /cmrad/ mrad
      common /cizc/ izc
      do 66 jj = 1, 2
      do 71 l = 1, ti
      if (pid(tid(k),jj).eq.tid(l)) then
      tpid(k,jj) = tid(l)
      if (gnr(l).le.0) then
      gn(ti) = j + 1
      gnr(ti) = j + 1
      end if
      goto 66
      end if
  71  continue
      ti = ti + 1
      tpid(k,jj) = pid(tid(k),jj)
      tid(ti) = tpid(k,jj)
      gn(ti) = j + 1
      gnr(ti) = j + 1
      tpid(ti,1) = 0
      tpid(ti,2) = 0
      if (iz(tid(ti)).eq.0) izc(ti) = izc(k) + 1
  66  continue
      return
      end


      subroutine addoff(j,k)
      parameter (np=32768, tnp=4096)
      include "blk.incl"
      integer*1 sex(np), tsex, gnr(tnp), gn(tnp), mmbit, mrad, izc(tnp)
      integer*4 pid(np,2), iz(np), i, j, k, ii, jj, kk, mt
      integer*4 tid(tnp), ti, tpid(tnp,2)
      integer*4 i2(nt), ip2(nt,2), n2, nbits2, nbits1
      common /com123a/ pid, iz, nn
      common /csex/ sex
      common /bcnt/ tid, tpid, ti, nbits1
      common /com123b/ ip2, i2, n2, mt, nbits2, mmbit
      common /comadd/ gn, gnr
      common /cmrad/ mrad
      common /cizc/ izc
      do 12 ii = 1, nn
      if ((pid(ii,1).ne.tid(k)).and.(pid(ii,2).ne.tid(k))) goto 12
      do 13 jj = 1, ti
      if (tid(jj).eq.ii) then
      if (gn(jj).eq.0) then
      gn(jj) = j + 1
      gnr(jj) = -gn(jj)
      end if
      if ((tpid(jj,1).eq.0).or.(tpid(jj,2).eq.0)) then
      tpid(jj,1) = pid(ii,1)
      tpid(jj,2) = pid(ii,2)
      end if
      goto 14
      end if
  13  continue
      ti = ti + 1
      tid(ti) = ii
      gn(ti) = j + 1
      gnr(ti) = -gn(ti)
      tpid(ti,1) = pid(ii,1)
      tpid(ti,2) = pid(ii,2)
      if (iz(tid(ti)).eq.0) izc(ti) = izc(k) + 1
  14  continue
      if (sex(tid(k)).eq.1) tsex = 2
      if (sex(tid(k)).eq.2) tsex = 1
      do 15 jj = 1, ti
      if (tid(jj).eq.pid(ii,tsex)) goto 12
  15  continue
      ti = ti + 1
      tid(ti) = pid(ii,tsex)
      gn(ti) = 0
      gnr(ti) = 0
      tpid(ti,1) = 0
      tpid(ti,2) = 0
      if (iz(tid(ti)).eq.0) izc(ti) = izc(k) + 1
  12  continue
      return
      end


      subroutine zachist()
      parameter (np=32768, tnp=4096)
      include "blk.incl"
      integer*1 mmbit
      integer*4 pid(np,2), iz(np), pere(ntp), mere(ntp)
      integer*4 tid(tnp), ti, tpid(tnp,2), tp(2), n
      integer*4 i1(np), ip1(np,2), i2(nt), ip2(nt,2), n1, n2
      integer*4 nbits1, nbits2, mt
      real*8 rr(ntp), rf
      logical izp(ntp)
      common /com123a/ pid, iz, nn
      common /bcnt/ tid, tpid, ti, nbits1
      common /com123b/ ip2, i2, n2, mt, nbits2, mmbit
      common /cprnt/ rr, pere, mere, izp, n
      common /comrf/ rf

      nbits1 = 0
      n1 = ti
      do 14 i = 1, ti
      i1(i) = tid(i)
      ip1(i,1) = tpid(i,1)
      ip1(i,2) = tpid(i,2)
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

      tp(1) = 0
      tp(2) = 0
cc ot potomkov
      do 41 j = 1, n1
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
      end if
      if (isk.ne.0) goto 48
  34  continue
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
      goto 48
  19  continue

      do 26 i = 1, n1
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

      rf = .0
      if (mt.eq.n2) return

      do 11 j = 1, n2
      rr(j) = .0
      izp(j) = .false.
      if (iz(i2(j)).ne.0) izp(j) = .true.
      if (ip2(j,1).eq.0) then
      pere(j) = 0
      mere(j) = 0
      goto 11
      end if
      do 10 k = 1, n2
      if (ip2(j,1).eq.i2(k)) pere(j) = k
      if (ip2(j,2).eq.i2(k)) mere(j) = k
  10  continue
  11  continue
      n = n2
      call parente
      do 12 j = 1, n2
      rf = rf + rr(j)
  12  continue
      
  
      return
      end


      subroutine parente
      implicit none
      include "blk.incl"
      integer*4 i, j, k, n, ii, jj, kk, is, il, ne, ie, ip, im
      integer*4 je, nd, iad, iaf, ifail, ig, ig1, ig2, eff
      integer*4 lliste(nlistx), nlist, iopt(4)
      integer*4 pere(ntp), mere(ntp), ped(2,ntp), point(ntp)
      integer*4 ord(ntp), rord(ntp), ndes(ntp)
      real*8 statis(ntp), x, f(0:ntp), l(ntp), d(ntp), rr(ntp)
      logical ts, ts2, chy, izp(ntp)
      character*128 sti, stc, jour
      common /cprnt/ rr, pere, mere, izp, n
      include "format.incl"

      ts=.true.
      ts2 = ts
      nlist = 0
	
      do 1 i = 1, n
      if (izp(i)) then
        nlist = nlist + 1
        lliste(nlist) = i
      end if
  1   continue

      do i = 1, nlist
      if (lliste(i).gt.n) then
          print 8512, lliste(i)
          stop
        end if
      end do

c *** renumerotation du plus vieux au plus jeune
      call comp_d (n, pere, mere, ped, ord, rord)

      do i = 1, nlist
        ndes(ord(lliste(i))) = 1
      end do

      do i = 1, n
        point(i)=0
        l(i)=0.
        d(i)=0.
      end do
      call meuw(n, ped, f, d, l, point,ndes)

      if (ts) then
        do i = 1, nlist
          j = lliste(i)
          if (ts2) then
            rr(j) = f(ord(j))
          end if
        end do
      end if
 8888 format (2i8,e14.6,2i2)

      n=n + 1
      do i=1, n
        point(i)=0
        l(i)=0.
      end do

      call statt(x,statis,0)
      do i = 1, nlist
        ii = 1
        ii = i + 1
        do j = ii, nlist
          ip = ord(lliste(i))
          im = ord(lliste(j))
          ped(1,n) = ip
          ped(2,n) = im
          call consang(n,ped,f,d,l,point)
          if (ts2) then
            rr(lliste(i)) = rr(lliste(i)) + f(n)
            rr(lliste(j)) = rr(lliste(j)) + f(n)
          end if
          call statt(f(n),statis,1)
        end do
      end do
      call statt(x,statis,2)
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
        do i=1, 100
          k=int(.5 + 100.*statis(5+i)/statis(3))
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
