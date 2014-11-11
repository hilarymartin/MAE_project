      program parente
      implicit none

      include "blk.incl"
      integer*4 i, j, k, n,ii,jj,kk,is,il,ne,eff,iopt(4),na,ie,ip,im
      integer*4 je,nd,nlist(2),iad,iaf,ifail,ig,ig1,ig2
      character*128 str,sti,stc,jour

      integer*4 lliste(nlistx,2)
      integer*4 pere(nt),mere(nt),ped(2,nt), point(nt)
      integer*2 sex(nt),annai(nt)
      integer*4 ord(nt),rord(nt),id(nt),ic(nt),ndes(nt)
      real*8 statis(106),x
      real*8 f(0:nt),l(nt),d(nt)
      logical ts,ts2(2,2),chy

      include "format.incl"

      call fdate(jour)
      PRINT 990,jour
      print 8001
      print 1001
      read(5,*) str
      print 1002,str

      print 8501
      read (5,*) sti
      print 8502,sti

      print 1003
      print *,'no if no output file'
      read (5,*) stc
      print 1004,stc
      ts=.true.
      if (stc.eq.'no') ts=.false.
      do i=1, 2
        do j=1, 2
          ts2(i,j)=ts
	  end do
	end do

      nlist(1)=0
      nlist(2)=0
      open (9,file=sti,form='formatted')
   18 read (9,*,end=19) i,ii
        if (ii.lt.1.or.ii.gt.2) then
          print 8004,i,ii
          stop
          end if
        nlist(ii)=nlist(ii)+1
        if (i.lt.1) stop 'number < 1 in the list'
        if (nlist(ii).gt.nlistx)
     *    stop 'Increase parameter nlistx in blk.incl and recompile !'
        lliste(nlist(ii),ii)=i
        goto 18
   19   close (9)

      if (nlist(2).eq.0) then
        print 8005,nlist(1)
      else
        print 8006,nlist(1)
        print 8007,nlist(2)
        end if
	
c COMMENTED TO OUTPUT LARGE KINSHIPS BY FAN, YURII
c limitation of output	
c      do i=1, 2
c        do j=1, 2
c	   if (nlist(i)*nlist(j).gt.1000000) ts2(i,j)=.false.
c	   end do
c	end do


c lecture du pedigree
      n=0
      chy=.false.
      open (1,file=str,form='FORMATTED')

    1 read (1,*,end=2) i,j,k,jj,is
      n=n+1
      if (jj.le.11) then
             jj=jj+2000
             chy=.true.
      else if (jj.le.100) then
             jj=jj+1900
             chy=.true.
             end if
      if (i.ne.n) then
         print 102,n,i
         stop
         end if
      if (i.gt.nt.or.j.gt.nt.or.k.gt.nt) then
         print 101
         stop
         end if
      pere(n)=j
      mere(n)=k
      if (j.lt.0) pere(n)=0
      if (k.lt.0) mere(n)=0
      annai(n)=jj
      sex(n)=is
      goto 1

    2 if (chy) print 103
      print 900, n
      close(1)
      

      do ig=1, 2
       do i=1, nlist(ig)
        if (lliste(i,ig).gt.n) then
          print 8512,lliste(i,ig)
          stop
          end if
        end do
       end do
      do i=1, n
        if (pere(i).gt.n) then ; print 8513; stop; end if
        if (mere(i).gt.n) then ; print 8514; stop; end if
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
        open (3,file=stc,form='formatted')
        do ig=1, 2
          do i=1, nlist(ig)
            j=lliste(i,ig)
            if (ts2(ig,ig)) write (3,8888) j,j,f(ord(j)),ig,ig
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
        print *
        print *
        print *,'Groups ',ig1,' x ',ig2
        print *,'****************'
        do i=1,nlist(ig1)
          ii=1
          if (ig1.eq.ig2) ii=i+1
          do j=ii, nlist(ig2)
           ip=ord(lliste(i,ig1))
           im=ord(lliste(j,ig2))
           ped(1,n)=ip
           ped(2,n)=im
           call consang(n,ped,f,d,l,point)
           if (ts2(ig1,ig2)) write (3,8888) lliste(i,ig1),
     *      lliste(j,ig2),f(n),ig1,ig2
           call statt(f(n),statis,1)
           end do
         end do
        call statt(x,statis,2)
        end do
       end do
       if (ts) close(3)

      call fdate(jour)
      print *
      PRINT 991,jour
	stop
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
        print 8517,statis(3),statis(4),statis(5)
        do i=1, 100
          k=int(.5 + 100.*statis(5+i)/statis(3))
          if (statis(5+i).gt.0) print 8518,i-1,i,statis(5+i),k,
     * 	    (etoile,l=1,k)
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
c        print *,nbit,k,n
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
c      print *, 'Coefficients initiaux calcules (parents + candidats) : ',npar
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


c      PRINT 3000,ninbr
c 3000 FORMAT (' Nb de parents consanguins : ',I8)
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
