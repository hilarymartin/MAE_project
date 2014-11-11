c selection recursive des genealogies utiles

      program selgen
      implicit none
      integer nt,ncarx
      parameter(nt=2000000,ncarx=5)
      character*10 anim(nt),pere(nt),mere(nt),base(nt),inc,trav
      character*1 elim
      integer prc(nt), mrc(nt), ndes(nt),ord(nt),rord(nt)
      integer*4 car(nt,ncarx), gen(nt)
      integer maxg,mgm(2),nc,ifail,necr,ncv,ncar
      integer i,j,k,l,n, nl, nbit,ip,im,ii,jj,npop,ng1, ng2,nn ,ngg
      logical t, pop(nt)
      character*128 str,sts,jour,stc
      external hpsort,hpperm,ipperm
      character*80 fmt
      character*1 ftx


      call fdate(jour)
      PRINT *, 'date and time  ',jour

      print *
      print *,'Name of pedigree file'
      read(5,*) str
      print *,'Pedigree file :',str
      inc='0000000000'
      print *,'Unknown parents are coded ',inc

      print *
      print *,'Name of file of reference animals'
      read (5,*) stc
      print *,'File of reference animals : ',stc

      print *
      print *,'Name of output file'
      read(5,*) sts
      print *,'Output file :',sts

      print *
      print *,'Maximum Nb of generations traced ? '
      read (5,*) maxg
      print *,'Maximum Nb of generations traced : ',maxg

      print *
      print *,'Number of other parameters to be read in the'
      print *,' pedigree file and written into output file ? '
      read (5,*) ncar
      if (ncar.gt.ncarx) then
         print *,'Error : increase parameter ncarx and recompile'
         stop
         end if
      print *,'Number of parameters to be read and written : ',ncar

      print *
      print *,'Elimination of useless pedigree ?'
      print *,' ie ancestors not in the ref list, without parents',
     *   ' and with<2 progeny (y / n) '
      read (5,*) elim
      if (elim.eq.'Y') elim='y'
      if (elim.ne.'y') elim='n'
      print *,'Elimination of useless pedigree  : ',elim

c lecture du fichier choix

        open (1,file=stc,form='FORMATTED')
        nc=1
   10   read (1,'(a10)',end=20) base(nc)
        nc=nc+1
        if (nc.gt.nt) stop 'Increase nt'
        goto 10
   20   close (1)
        nc=nc-1
        call hpsort(base,nc,1,10,ord,2,trav,ifail)

c lecture du fichier pedigree
      nl=1
      necr=0
      open (1,file=str,form='FORMATTED')
   11 read (1,*,end=21) anim(nl),pere(nl),mere(nl),
     *    (car(nl,i),i=1,ncar)

c   11 read (1,99,end=21) anim(nl),pere(nl),mere(nl)
c   99 format(a10,1x,a10,1x,a10,i5,3i2)
      if (necr.lt.10) then
        print '(i8,2x,3(a10,1x),5i8)',
c     *    nl,anim(nl),pere(nl),mere(nl)
     *    nl,anim(nl),pere(nl),mere(nl),(car(nl,i),i=1,ncar)
        necr=necr+1
        end if
      NL=NL+1
      if (nl.gt.nt) stop 'Increase parameter nt'
      GOTO 11
   21 close (1)
      nl=nl-1

c tri du fichier genealogie
      call hpsort(anim,nl,1,10,ord,2,trav,ifail)
      call hpperm(pere,nl,ord,trav,ifail)
      call hpperm(mere,nl,ord,trav,ifail)
      do ii=1, ncar
        call ipperm(car(1,ii),nl,ord,ifail)
        end do

c selection de la base
      j=1
      ncv=0
      do i=1, nl
        pop(i)=.false.
        gen(i)=-1
        end do
      do i=1, nc
        do while(j.lt.nl.and.anim(j).lt.base(i))
           j=j+1
           end do
        if (anim(j).eq.base(i)) then
           pop(j)=.true.
           gen(j)=1
	else
	   ncv=ncv+1
           end if
        end do
      if (ncv.gt.0) then
        print *,ncv,' Individuals in the list',
     * 	' but not in the pedigree file'
        print *,'They are not in the output file'
        end if

c recodif des peres et meres
      call hpsort(pere,nl,1,10,ord,2,trav,ifail)

      j=1
      do i=1, nl
        if (pere(i).eq.inc) then
           prc(ord(i))=0
        else
           do while(pere(i).gt.anim(j).and.j.lt.nl)
             j=j+1
             end do
           if (pere(i).eq.anim(j)) then
              prc(ord(i))=j
           else
              prc(ord(i))=0
              end if
           end if
        end do
      do i=1, nl
         rord(ord(i))=i
         end do
      call hpperm(pere,nl,rord,trav,ifail)

      call hpsort(mere,nl,1,10,ord,2,trav,ifail)
      j=1
      do i=1, nl
        if (mere(i).eq.inc) then
           mrc(ord(i))=0
        else
           do while(mere(i).gt.anim(j).and.j.lt.nl)
             j=j+1
             end do
           if (mere(i).eq.anim(j)) then
               mrc(ord(i))=j
           else
              mrc(ord(i))=0
              end if
           end if
        end do

      do i=1, nl
         rord(ord(i))=i
         end do
      call hpperm(mere,nl,rord,trav,ifail)

      print *,'Nb of individuals in the pedigree file :' , nl
      print *,'Nb of individuals in the list          :' , nc
c a ce stade les recodif sont terminees

c selection des ancetres des individus de base
      print *
      print *,'Selection of ancestors'
      print *,'----------------------'
      j=1
      nbit=0
      do while(nbit.lt.20.and.(j.gt.0.or.k.gt.0))
         nbit=nbit+1
         j=0
         k=0
         do i=1, nl
           if (pop(i)) then
              ip=prc(i)
              if (ip.gt.0) then
                if (.not.pop(ip)) then
                  pop(ip)=.true.
                  gen(ip)=gen(i)+1
                  j=j+1
                else if (gen(ip).gt.gen(i)+1) then
                  gen(ip)=gen(i)+1
                  k=k+1
                  end if
                end if
              im=mrc(i)
              if (im.gt.0) then
                if (.not.pop(im)) then
                  pop(im)=.true.
                  gen(im)=gen(i)+1
                  j=j+1
                else if (gen(im).gt.gen(i)+1) then
                  gen(im)=gen(i)+1
                  k=k+1
                  end if
                end if
              end if
           end do
         print *,'At iteration ',nbit,', ',j,' ancestors selected'
	 print *,'                   and ',k,' generation Nr updated'
         end do

c restriction a mgm generations
      j=0
      do i=1, nl
        if (gen(i).gt.maxg) then
	  if (prc(i).ne.0) then
            pop(prc(i))=.false.
            prc(i)=0
	    j=j+1
	    end if
	  if (mrc(i).ne.0) then
            pop(mrc(i))=.false.
            mrc(i)=0
	    j=j+1
	    end if	
          end if
        end do
      if (j.gt.0) then
        print *
	print *,'Elimination due to high number of generations : ',j
	end if	


c restriction aux individus utiles
      if (elim.eq.'y') then
      nbit=0
      j=1
      k=nl
      print *
      print *,'Elimination of useless ancestors'
      print *,'--------------------------------'

      do while(j.gt.0 .and. nbit.le.20)
        nbit=nbit+1
        do i=1, nl
           ndes(i)=0
           end do
        do i=1, nl
          if (pop(i).and.prc(i).ne.0) ndes(prc(i))=ndes(prc(i)) + 1
          if (pop(i).and.mrc(i).ne.0) ndes(mrc(i))=ndes(mrc(i)) + 1
          end do
        j=0
        do i=1, nl
          ip=prc(i)
          im=mrc(i)
          if (ip.eq.0 .and. im.eq.0 .and.ndes(i).lt.2
     *       .and. gen(i).ne.1 .and. pop(i)) then
            pop(i)=.false.
            j=j+1
            end if
          end do
        do i=1, nl
          ip=prc(i)
          im=mrc(i)
          if (ip.gt.0 .and. .not.pop(ip)) prc(i)=0
          if (im.gt.0 .and. .not.pop(im)) mrc(i)=0
          end do
        k=k-j
        print *,'At iteration ',nbit,', elimination of ',j,' ancestors'
        print *,'          There remains ',k,' individuals in the list'
      end do
      end if

c recodification des individus choisis
      j=0
      do i=1, nl
        ord(i)=0
        if (pop(i)) then
          j=j+1
          ord(i)=j
          end if
        end do

c recodification des peres et meres et svg
      necr=0
      open (2,file=sts,form='formatted')
      if (ncar.eq.0) then
        fmt='(3i8,3(1x,a10),i8,3i3)'
      else
        write(ftx,'(i1)') ncar
        fmt='(3i8,'//ftx//'i8,3(1x,a10),i8,3i3)'
	end if
      print *,'Output format : ',fmt
      do i=1, nl
        if (pop(i)) then
          if (prc(i).gt.0) prc(i)=ord(prc(i))
          if (mrc(i).gt.0) mrc(i)=ord(mrc(i))

          write (2,fmt)
     *      ord(i),prc(i),mrc(i),(car(i,j),j=1,ncar),
     *      anim(i),pere(i),mere(i),ndes(i),gen(i)
          if (necr.lt.5) print fmt,
     *      ord(i),prc(i),mrc(i),(car(i,j),j=1,ncar),
     *      anim(i),pere(i),mere(i)
          necr=necr+1
          end if
        end do
      print *	
      print *,necr,' records written in the file ',sts
      print *,'END of SELECTION'
      stop
      END


*DECK HPPERM
      SUBROUTINE HPPERM (HX, N, IPERM, WORK, IER)
C***BEGIN PROLOGUE  HPPERM
C***PURPOSE  Rearrange a given array according to a prescribed
C            permutation vector.
C***LIBRARY   SLATEC
C***CATEGORY  N8
C***TYPE      CHARACTER (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
C***KEYWORDS  APPLICATION OF PERMUTATION TO DATA VECTOR
C***AUTHOR  McClain, M. A., (NIST)
C           Rhoads, G. S., (NBS)
C***DESCRIPTION
C
C         HPPERM rearranges the data vector HX according to the
C         permutation IPERM: HX(I) <--- HX(IPERM(I)).  IPERM could come
C         from one of the sorting routines IPSORT, SPSORT, DPSORT or
C         HPSORT.
C
C     Description of Parameters
C         HX - input/output -- character array of values to be
C                 rearranged.
C         N - input -- number of values in character array HX.
C         IPERM - input -- permutation vector.
C         WORK - character variable which must have a length
C                   specification at least as great as that of HX.
C         IER - output -- error indicator:
C             =  0  if no error,
C             =  1  if N is zero or negative,
C             =  2  if work array is not long enough,
C             =  3  if IPERM is not a valid permutation.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   901004  DATE WRITTEN
C   920507  Modified by M. McClain to revise prologue text and to add
C           check for length of work array.
C***END PROLOGUE  HPPERM
      INTEGER N, IPERM(*), I, IER, INDX, INDX0, ISTRT
      CHARACTER*(*) HX(*), WORK
C***FIRST EXECUTABLE STATEMENT  HPPERM
      IER=0
      IF(N.LT.1)THEN
         IER=1
         CALL XERMSG ('SLATEC', 'HPPERM',
     +    'The number of values to be rearranged, N, is not positive.',
     +    IER, 1)
         RETURN
      ENDIF
      IF(LEN(WORK).LT.LEN(HX(1)))THEN
         IER=2
         CALL XERMSG ('SLATEC', 'HPPERM',
     +    'The length of the work variable, WORK, is too short.',IER,1)
         RETURN
      ENDIF
C
C     CHECK WHETHER IPERM IS A VALID PERMUTATION
C
      DO 100 I=1,N
         INDX=ABS(IPERM(I))
         IF((INDX.GE.1).AND.(INDX.LE.N))THEN
            IF(IPERM(INDX).GT.0)THEN
               IPERM(INDX)=-IPERM(INDX)
               GOTO 100
            ENDIF
         ENDIF
         IER=3
         CALL XERMSG ('SLATEC', 'HPPERM',
     +    'The permutation vector, IPERM, is not valid.', IER, 1)
         RETURN
  100 CONTINUE
C
C     REARRANGE THE VALUES OF HX
C
C     USE THE IPERM VECTOR AS A FLAG.
C     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
C
      DO 330 ISTRT = 1 , N
         IF (IPERM(ISTRT) .GT. 0) GOTO 330
         INDX = ISTRT
         INDX0 = INDX
         WORK = HX(ISTRT)
  320    CONTINUE
         IF (IPERM(INDX) .GE. 0) GOTO 325
            HX(INDX) = HX(-IPERM(INDX))
            INDX0 = INDX
            IPERM(INDX) = -IPERM(INDX)
            INDX = IPERM(INDX)
            GOTO 320
  325    CONTINUE
         HX(INDX0) = WORK
  330 CONTINUE
C
      RETURN
      END


*DECK HPSORT
      SUBROUTINE HPSORT (HX, N, STRBEG, STREND, IPERM, KFLAG, WORK, IER)
C***BEGIN PROLOGUE  HPSORT
C***PURPOSE  Return the permutation vector generated by sorting a
C            substring within a character array and, optionally,
C            rearrange the elements of the array.  The array may be
C            sorted in forward or reverse lexicographical order.  A
C            slightly modified quicksort algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A1C, N6A2C
C***TYPE      CHARACTER (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
C***KEYWORDS  PASSIVE SORTING, SINGLETON QUICKSORT, SORT, STRING SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Rhoads, G. S., (NBS)
C           Sullivan, F. E., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   HPSORT returns the permutation vector IPERM generated by sorting
C   the substrings beginning with the character STRBEG and ending with
C   the character STREND within the strings in array HX and, optionally,
C   rearranges the strings in HX.   HX may be sorted in increasing or
C   decreasing lexicographical order.  A slightly modified quicksort
C   algorithm is used.
C
C   IPERM is such that HX(IPERM(I)) is the Ith value in the
C   rearrangement of HX.  IPERM may be applied to another array by
C   calling IPPERM, SPPERM, DPPERM or HPPERM.
C
C   An active sort of numerical data is expected to execute somewhat
C   more quickly than a passive sort because there is no need to use
C   indirect references. But for the character data in HPSORT, integers
C   in the IPERM vector are manipulated rather than the strings in HX.
C   Moving integers may be enough faster than moving character strings
C   to more than offset the penalty of indirect referencing.
C
C   Description of Parameters
C      HX - input/output -- array of type character to be sorted.
C           For example, to sort a 80 element array of names,
C           each of length 6, declare HX as character HX(100)*6.
C           If ABS(KFLAG) = 2, then the values in HX will be
C           rearranged on output; otherwise, they are unchanged.
C      N  - input -- number of values in array HX to be sorted.
C      STRBEG - input -- the index of the initial character in
C               the string HX that is to be sorted.
C      STREND - input -- the index of the final character in
C               the string HX that is to be sorted.
C      IPERM - output -- permutation array such that IPERM(I) is the
C              index of the string in the original order of the
C              HX array that is in the Ith location in the sorted
C              order.
C      KFLAG - input -- control parameter:
C            =  2  means return the permutation vector resulting from
C                  sorting HX in lexicographical order and sort HX also.
C            =  1  means return the permutation vector resulting from
C                  sorting HX in lexicographical order and do not sort
C                  HX.
C            = -1  means return the permutation vector resulting from
C                  sorting HX in reverse lexicographical order and do
C                  not sort HX.
C            = -2  means return the permutation vector resulting from
C                  sorting HX in reverse lexicographical order and sort
C                  HX also.
C      WORK - character variable which must have a length specification
C             at least as great as that of HX.
C      IER - output -- error indicator:
C          =  0  if no error,
C          =  1  if N is zero or negative,
C          =  2  if KFLAG is not 2, 1, -1, or -2,
C          =  3  if work array is not long enough,
C          =  4  if string beginning is beyond its end,
C          =  5  if string beginning is out-of-range,
C          =  6  if string end is out-of-range.
C
C     E X A M P L E  O F  U S E
C
C      CHARACTER*2 HX, W
C      INTEGER STRBEG, STREND
C      DIMENSION HX(10), IPERM(10)
C      DATA (HX(I),I=1,10)/ '05','I ',' I','  ','Rs','9R','R9','89',
C     1     ',*','N"'/
C      DATA STRBEG, STREND / 1, 2 /
C      CALL HPSORT (HX,10,STRBEG,STREND,IPERM,1,W)
C      PRINT 100, (HX(IPERM(I)),I=1,10)
C 100 FORMAT (2X, A2)
C      STOP
C      END
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified by John A. Wisniewski to use the Singleton
C           quicksort algorithm.
C   811001  Modified by Francis Sullivan for string data.
C   850326  Documentation slightly modified by D. Kahaner.
C   870423  Modified by Gregory S. Rhoads for passive sorting with the
C           option for the rearrangement of the original data.
C   890620  Algorithm for rearranging the data vector corrected by R.
C           Boisvert.
C   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
C   920507  Modified by M. McClain to revise prologue text.
C   920818  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
C***END PROLOGUE  HPSORT
C     .. Scalar Arguments ..
      INTEGER IER, KFLAG, N, STRBEG, STREND
      CHARACTER * (*) WORK
C     .. Array Arguments ..
      INTEGER IPERM(*)
      CHARACTER * (*) HX(*)
C     .. Local Scalars ..
      REAL R
      INTEGER I, IJ, INDX, INDX0, IR, ISTRT, J, K, KK, L, LM, LMT, M,
     +        NN, NN2
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT, LEN
C***FIRST EXECUTABLE STATEMENT  HPSORT
      IER = 0
      NN = N
      IF (NN .LT. 1) THEN
         IER = 1
         CALL XERMSG ('SLATEC', 'HPSORT',
     +    'The number of values to be sorted, N, is not positive.',
     +    IER, 1)
         RETURN
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         IER = 2
         CALL XERMSG ('SLATEC', 'HPSORT',
     +    'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.',
     +    IER, 1)
         RETURN
      ENDIF
C
      IF(LEN(WORK) .LT. LEN(HX(1))) THEN
         IER = 3
         CALL XERMSG ('SLATEC',' HPSORT',
     +    'The length of the work variable, WORK, is too short.',
     +    IER, 1)
         RETURN
      ENDIF
      IF (STRBEG .GT. STREND) THEN
         IER = 4
         CALL XERMSG ('SLATEC', 'HPSORT',
     +    'The string beginning, STRBEG, is beyond its end, STREND.',
     +    IER, 1)
         RETURN
      ENDIF
      IF (STRBEG .LT. 1 .OR. STRBEG .GT. LEN(HX(1))) THEN
         IER = 5
         CALL XERMSG ('SLATEC', 'HPSORT',
     +    'The string beginning, STRBEG, is out-of-range.',
     +    IER, 1)
         RETURN
      ENDIF
      IF (STREND .LT. 1 .OR. STREND .GT. LEN(HX(1))) THEN
         IER = 6
         CALL XERMSG ('SLATEC', 'HPSORT',
     +    'The string end, STREND, is out-of-range.',
     +    IER, 1)
         RETURN
      ENDIF
C
C     Initialize permutation vector
C
      DO 10 I=1,NN
         IPERM(I) = I
   10 CONTINUE
C
C     Return if only one value is to be sorted
C
      IF (NN .EQ. 1) RETURN
C
C     Sort HX only
C
      M = 1
      I = 1
      J = NN
      R = .375E0
C
   20 IF (I .EQ. J) GO TO 70
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location L
C
      IJ = I + INT((J-I)*R)
      LM = IPERM(IJ)
C
C     If first element of array is greater than LM, interchange with LM
C
      IF (HX(IPERM(I))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND)) THEN
         IPERM(IJ) = IPERM(I)
         IPERM(I) = LM
         LM = IPERM(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than LM, interchange with LM
C
      IF (HX(IPERM(J))(STRBEG:STREND) .LT. HX(LM)(STRBEG:STREND)) THEN
         IPERM(IJ) = IPERM(J)
         IPERM(J) = LM
         LM = IPERM(IJ)
C
C        If first element of array is greater than LM, interchange
C        with LM
C
         IF (HX(IPERM(I))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND))
     +      THEN
               IPERM(IJ) = IPERM(I)
               IPERM(I) = LM
               LM = IPERM(IJ)
         ENDIF
      ENDIF
      GO TO 50
   40 LMT = IPERM(L)
      IPERM(L) = IPERM(K)
      IPERM(K) = LMT
C
C     Find an element in the second half of the array which is smaller
C     than LM
C
   50 L = L-1
      IF (HX(IPERM(L))(STRBEG:STREND) .GT. HX(LM)(STRBEG:STREND))
     +    GO TO 50
C
C     Find an element in the first half of the array which is greater
C     than LM
C
   60 K = K+1
      IF (HX(IPERM(K))(STRBEG:STREND) .LT. HX(LM)(STRBEG:STREND))
     +   GO TO 60
C
C     Interchange these elements
C
      IF (K .LE. L) GO TO 40
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 80
C
C     Begin again on another portion of the unsorted array
C
   70 M = M-1
      IF (M .EQ. 0) GO TO 110
      I = IL(M)
      J = IU(M)
C
   80 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   90 I = I+1
      IF (I .EQ. J) GO TO 70
      LM = IPERM(I+1)
      IF (HX(IPERM(I))(STRBEG:STREND) .LE. HX(LM)(STRBEG:STREND))
     +   GO TO 90
      K = I
C
  100 IPERM(K+1) = IPERM(K)
      K = K-1
C
      IF (HX(LM)(STRBEG:STREND) .LT. HX(IPERM(K))(STRBEG:STREND))
     +    GO TO 100
      IPERM(K+1) = LM
      GO TO 90
C
C     Clean up
C
  110 IF (KFLAG .LE. -1) THEN
C
C        Alter array to get reverse order, if necessary
C
         NN2 = NN/2
         DO 120 I=1,NN2
           IR = NN-I+1
           LM = IPERM(I)
           IPERM(I) = IPERM(IR)
           IPERM(IR) = LM
  120    CONTINUE
      ENDIF
C
C     Rearrange the values of HX if desired
C
      IF (KK .EQ. 2) THEN
C
C        Use the IPERM vector as a flag.
C        If IPERM(I) < 0, then the I-th value is in correct location
C
         DO 140 ISTRT=1,NN
            IF (IPERM(ISTRT) .GE. 0) THEN
               INDX = ISTRT
               INDX0 = INDX
               WORK = HX(ISTRT)
  130          IF (IPERM(INDX) .GT. 0) THEN
                  HX(INDX) = HX(IPERM(INDX))
                  INDX0 = INDX
                  IPERM(INDX) = -IPERM(INDX)
                  INDX = ABS(IPERM(INDX))
                  GO TO 130
               ENDIF
               HX(INDX0) = WORK
            ENDIF
  140    CONTINUE
C
C        Revert the signs of the IPERM values
C
         DO 150 I=1,NN
            IPERM(I) = -IPERM(I)
  150    CONTINUE
C
      ENDIF
C
      RETURN
      END

*DECK IPPERM
      SUBROUTINE IPPERM (IX, N, IPERM, IER)
C***BEGIN PROLOGUE  IPPERM
C***PURPOSE  Rearrange a given array according to a prescribed
C            permutation vector.
C***LIBRARY   SLATEC
C***CATEGORY  N8
C***TYPE      INTEGER (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
C***KEYWORDS  APPLICATION OF PERMUTATION TO DATA VECTOR
C***AUTHOR  McClain, M. A., (NIST)
C           Rhoads, G. S., (NBS)
C***DESCRIPTION
C
C         IPPERM rearranges the data vector IX according to the
C         permutation IPERM: IX(I) <--- IX(IPERM(I)).  IPERM could come
C         from one of the sorting routines IPSORT, SPSORT, DPSORT or
C         HPSORT.
C
C     Description of Parameters
C         IX - input/output -- integer array of values to be rearranged.
C         N - input -- number of values in integer array IX.
C         IPERM - input -- permutation vector.
C         IER - output -- error indicator:
C             =  0  if no error,
C             =  1  if N is zero or negative,
C             =  2  if IPERM is not a valid permutation.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   900618  DATE WRITTEN
C   920507  Modified by M. McClain to revise prologue text.
C***END PROLOGUE  IPPERM
      INTEGER IX(*), N, IPERM(*), I, IER, INDX, INDX0, ITEMP, ISTRT
C***FIRST EXECUTABLE STATEMENT  IPPERM
      IER=0
      IF(N.LT.1)THEN
         IER=1
         CALL XERMSG ('SLATEC', 'IPPERM',
     +    'The number of values to be rearranged, N, is not positive.',
     +    IER, 1)
         RETURN
      ENDIF
C
C     CHECK WHETHER IPERM IS A VALID PERMUTATION
C
      DO 100 I=1,N
         INDX=ABS(IPERM(I))
         IF((INDX.GE.1).AND.(INDX.LE.N))THEN
            IF(IPERM(INDX).GT.0)THEN
               IPERM(INDX)=-IPERM(INDX)
               GOTO 100
            ENDIF
         ENDIF
         IER=2
         CALL XERMSG ('SLATEC', 'IPPERM',
     +    'The permutation vector, IPERM, is not valid.', IER, 1)
         RETURN
  100 CONTINUE
C
C     REARRANGE THE VALUES OF IX
C
C     USE THE IPERM VECTOR AS A FLAG.
C     IF IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
C
      DO 330 ISTRT = 1 , N
         IF (IPERM(ISTRT) .GT. 0) GOTO 330
         INDX = ISTRT
         INDX0 = INDX
         ITEMP = IX(ISTRT)
  320    CONTINUE
         IF (IPERM(INDX) .GE. 0) GOTO 325
            IX(INDX) = IX(-IPERM(INDX))
            INDX0 = INDX
            IPERM(INDX) = -IPERM(INDX)
            INDX = IPERM(INDX)
            GOTO 320
  325    CONTINUE
         IX(INDX0) = ITEMP
  330 CONTINUE
C
      RETURN
      END

      subroutine xermsg (c1,c2,c3,j,i)
      character*10 c1,c2
      character*80 c3
      integer i,j

      print *,'Message d erreur'
      print *,'****************'
      print *,'Subroutine : ',c2
      print *,c3
      if (i.eq.1) stop
      return
      end


