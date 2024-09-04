!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains miscellaneous utility routines.  They need to be revised
!  and reworked, and none of them should be regarded as publicly callable at this
!  time.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module utils

interface prin2
   module procedure prin2_0
   module procedure prin2_1
   module procedure prin2_2
   module procedure prin2_3
end interface prin2

interface prind
   module procedure prind_0
   module procedure prind_1
   module procedure prind_2
end interface prind

interface princ
   module procedure princ_0
   module procedure princ_1
   module procedure princ_2
end interface princ

interface prinz
   module procedure prinz_0
   module procedure prinz_1
   module procedure prinz_2
end interface prinz

interface prini
   module procedure prini_0
   module procedure prini_1
   module procedure prini_2
end interface prini

interface prinl
   module procedure prinl_0
   module procedure prinl_1
   module procedure prinl_2
end interface prinl

contains

subroutine prin2_0(str,a)
implicit double precision (a-h,o-z)

double precision a
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(6(2x,e15.7))") a

end subroutine

subroutine prin2_1(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(6(2x,e15.7))") a

end subroutine

subroutine prin2_2(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(8(2x,e15.7))") a

end subroutine


subroutine prin2_3(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,e15.7))",a

write (13,*) str
write (13,"(6(2x,e15.7))") a

end subroutine

subroutine prind_0(str,a)
implicit double precision (a-h,o-z)

double precision :: a
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine

subroutine prind_1(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine


subroutine prind_2(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine



subroutine princ_0(str,a)
implicit double precision (a-h,o-z)
double complex               :: a
character(len=*), intent(in) :: str

print *,str
print "(3(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(3(d15.7,',',d15.7,2X))") a

end subroutine


subroutine princ_1(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(3(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(3(d15.7,',',d15.7,2X))") a

end subroutine


subroutine princ_2(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(3(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(3(d15.7,',',d15.7,2X))") a

end subroutine


subroutine prinz_0(str,a)
implicit double precision (a-h,o-z)
double complex               :: a
character(len=*), intent(in) :: str

print *,str
print "(1(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(1(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prinz_1(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(1(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(1(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prinz_2(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(1(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(1(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prini_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I9))",a

write (13,*) str
write (13,"(8(2x,I9))") a

end subroutine

subroutine prini_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine



subroutine prini_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine

subroutine prina(str)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str

print *,str
write (13,*) str

end subroutine


subroutine prinl_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine



subroutine elapsed(t)
implicit double precision (a-h,o-z)
integer*8 i,irate
call system_clock(i,irate)

dd = i
dd = dd/irate
t = dd
return
end subroutine




subroutine insort(k,a)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
double precision, intent (inout) :: a(k)

if (k .le. 1) return

do i=2,k
val=a(i)
j=i-1
do while (j .ge. 1 .AND. a(j) .gt. val) 
a(j+1)=a(j)
j=j-1
end do
a(j+1)=val
end do
end subroutine


subroutine insorti(k,ia)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)

if (k .le. 1) return

do i=2,k
ival=ia(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)=ia(j)
j=j-1
end do
ia(j+1)=ival
end do
end subroutine


subroutine insorti2(k,ia,idxs)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),idxs(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
idxval = idxs(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
idxs(j+1) = idxs(j)
j=j-1
end do
ia(j+1)   = ival
idxs(j+1) = idxval
end do

end subroutine

subroutine insorti3(k,ia,ia2,ia3)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),ia2(k),ia3(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
ival2  = ia2(i)
ival3  = ia3(i)

j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
ia2(j+1)  = ia2(j)
ia3(j+1)  = ia3(j)
j=j-1
end do
ia(j+1)   = ival
ia2(j+1)  = ival2
ia3(j+1)  = ival3
end do

end subroutine


subroutine insort3(k,ia,vals)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)
double complex                   :: vals(k)
double complex                   :: val

if (k .le. 1) return

do i=2,k
ival   = ia(i)
val    = vals(i)

j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
vals(j+1) = vals(j)
j=j-1
end do
ia(j+1)   = ival
vals(j+1) = val
end do
end subroutine


subroutine quicksort(n,vals)
implicit double precision (a-h,o-z)
dimension istack(2,20000)
dimension vals(1),idxs(1)
!
!       Sort a list of double precision numbers.
!
k        = 60

if (n .lt. k) then
call insort(n,vals)
return
endif

maxstack = 10000

m = 1
istack(1,1) = 1
istack(2,1) = n
!
 1000 continue
if (m .eq. 0) goto 1100
i1 = istack(1,m)
i2 = istack(2,m)
m=m-1
!
l = i2-i1+1
if (l .le. k) then
call insort(l,vals(i1))
goto 1000
endif
!
! Otherwise perform quicksort.
!
call quicksort01(vals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
! if (m+2 .ge. maxstack) then
! print *,"quicksort out of memory"
! stop
! endif
!
!  Make sure the smaller half is processed first to reduce storage
!  to O(logn).
!             
n1 = i3-i1+1
n2 = i2-i3
!
if (n2 .lt. n1) then
!
m = m+1
istack(1,m) = i1
istack(2,m) = i3
!
m = m+1
istack(1,m) = i3+1
istack(2,m) = i2
!
else
!
m = m+1
istack(1,m) = i3+1
istack(2,m) = i2
!
m = m+1
istack(1,m) = i1
istack(2,m) = i3
!
endif
!
goto 1000
 1100 continue 
end subroutine


        subroutine quicksort01(vals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        val  = vals(ipiv)
!
!       Swap the pivot element and the last element.
!
        vals(ipiv) = vals(i2)
        vals(i2)   = val
!
       i3 = i1
!
        do 1000 i=i1,i2-1
!
        if( vals(i) .lt. val) then
        d  = vals(i)
!
        vals(i)  = vals(i3)
        vals(i3) = d       
!
        i3=i3+1
        endif
!
 1000 continue
!
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
!
        end subroutine



        subroutine quicksorti(n,ivals)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti(n,ivals)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti(l,ivals(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti0(ivals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti0(ivals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv) = ivals(i2)
        ivals(i2)   = ival
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
!
        ivals(i)  = ivals(i3)
        ivals(i3) = id
!
        i3=i3+1
        endif
 1000 continue
!
        id        = ivals(i3)
        ivals(i3) = ivals(i2)
        ivals(i2) = id
!
        end subroutine



        subroutine quicksorti2(n,ivals,ivals2)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1),ivals2(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti2(n,ivals,ivals2)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti2(l,ivals(i1),ivals2(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti20(ivals,ivals2,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti20(ivals,ivals2,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1),ivals2(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
        ival2 = ivals2(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv)  = ivals(i2)
        ivals2(ipiv) = ivals2(i2)
        ivals(i2)    = ival
        ivals2(i2)   = ival2
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
        id2 = ivals2(i)
!
        ivals(i)   = ivals(i3)
        ivals2(i)  = ivals2(i3)
        ivals(i3)  = id
        ivals2(i3) = id2
!
        i3=i3+1
        endif
 1000 continue
!
        id         = ivals(i3)
        id2        = ivals2(i3)
        ivals(i3)  = ivals(i2)
        ivals2(i3) = ivals2(i2)
        ivals(i2)  = id
        ivals2(i2) = id2
        
!
        end subroutine



        subroutine quicksorti3(n,ivals,ivals2,ivals3)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000),ivals(1),ivals2(1),ivals3(1)
!
!       Sort a list of integers.
!
        if (n .lt. 60) then
        call insorti3(n,ivals,ivals2,ivals3)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
!
!
        l = i2-i1+1
        if (l .le. k) then
        call insorti3(l,ivals(i1),ivals2(i1),ivals3(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksorti30(ivals,ivals2,ivals3,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
!
        if (n2 .lt. n1) then
!
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        else
!
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
!
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
!
        endif
!
        goto 1000
 1100 continue
        end subroutine


        subroutine quicksorti30(ivals,ivals2,ivals3,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1),ivals2(1),ivals3(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        ival  = ivals(ipiv)
        ival2 = ivals2(ipiv)
        ival3 = ivals3(ipiv)
!
!       Swap the pivot element and the last element.
!
        ivals(ipiv)  = ivals(i2)
        ivals2(ipiv) = ivals2(i2)
        ivals3(ipiv) = ivals3(i2)
        ivals(i2)    = ival
        ivals2(i2)   = ival2
        ivals3(i2)   = ival3
!
        i3 = i1
!
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
        id2 = ivals2(i)
        id3 = ivals3(i)
!
        ivals(i)   = ivals(i3)
        ivals2(i)  = ivals2(i3)
        ivals3(i)  = ivals3(i3)
        ivals(i3)  = id
        ivals2(i3) = id2
        ivals3(i3) = id3
!
        i3=i3+1
        endif
 1000 continue
!
        id         = ivals(i3)
        id2        = ivals2(i3)
        id3        = ivals3(i3)
        ivals(i3)  = ivals(i2)
        ivals2(i3) = ivals2(i2)
        ivals3(i3) = ivals3(i2)
        ivals(i2)  = id
        ivals2(i2) = id2
        ivals3(i2) = i3
        
!
        end subroutine



        subroutine iremove(n,ia,m,ib)
        implicit double precision (a-h,o-z)
        dimension ia(n),ib(m)
!
!       Remove from the list ia of length n all integers appearing in 
!       the list ib of length m.  Both the list ia and the list ib
!       must be sorted before this call is made.  The results will
!       also be sorted.
!

!        call quicksorti(n,ia)
!        call quicksorti(m,ib)

        isrc = 1
        itar = 1
        ii   = 1
 1000 continue

        if (ii .gt. m)   goto 2000
        if (isrc .gt. n) goto 3000
        if (ia(isrc) .gt. ib(ii)) then
        ii=ii+1
        goto 1000
        endif

        if (ia(isrc) .lt. ib(ii)) then         
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 1000
        endif
        isrc=isrc+1
        goto 1000

 2000 continue
        if (isrc .gt. n) goto 3000
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 2000

 3000 continue
        n = itar-1
        return        
        
        end



function eye(n) result(a)
implicit double precision (a-h,o-z)

double precision, allocatable :: a(:,:)
integer n

!
!  Return an identity matrix which is dimensioned (n,n).
!

allocate(a(n,n))
a = 0
do i=1,n
a(i,i)=1.0d0
end do

end function



subroutine plot_functions(scriptname,filename,nfuns,xlabel,ylabel,iflogx,iflogy,x1_,x2_,y1_,y2_, &
   n,xs,ys,legend_loc,legends,styles)
implicit double precision (a-h,o-z)

character(len=*)               :: scriptname,filename,legends,legend_loc,styles
character(len=*)               :: xlabel,ylabel
double precision               :: xs(n),ys(nfuns,n)
character(len=:), allocatable  :: command, legend, style

!
!  Produce a python script which generates a PDF file containing a plot of one or more
!  functions of one variable suitable for inclusion in a paper.
!
!  Input parameters:
!    scriptname - the name of the python script
!    filename - the name of the PDF file to produce
!
!    xlabel - a label for the x-axis (no label will be present if the string is empty)
!    ylabel - a label for the y-axis (no label will be present if the string is empty)
!
!    iflogx - an integer parameter specifying whether the x-axis should be logarithm
!      scale or not
!    iflogy - an integer parameter specifying whether the y-axis should be logarithm
!      scale or not
!
!    (x1,x2,y1,y2) - extents of the axes ---- if x2 <= x1 then these will be set
!      automatically; likewise if y2 <= y1
!
!    legend_loc - a string specifying the location for the legend --- this can
!       be "upper right", "upper left", "lower right", etc  OR a blank string if 
!       no legend is to be included OR "best" for automatic placement
!
!    legends - a string specifying the legend labels for each function ...
!       each should be terminated by an asterik "*" so "1*2*3*" specifies
!       the label 1 for the first function, 2 for the second, and 3 for the
!       third
!
!       this is ignored if legend_loc is blank
!
!    styles - a string specifying styles for each function ... separated as in
!       the legends string ... ignored if empty
!     
!    n - the number of point in the graph of the function to be specified
!    xs - the x-coordinates of the points to plot
!    ys - an (nfuns,n) array whose jth column gives the y-coordinates on
!        the graph of the jth function
!
!  Output parameters:
!    N/A
!


x1 = x1_
x2 = x2_
y1 = y1_
y2 = y2_

if (x2 .le. x1) then
x1 =  1d300
x2 = -1d300
do i=1,n
x1 = min(x1,xs(i))
x2 = max(x2,xs(i))
end do

endif

if (y2 .le. y1) then
y1 =  1d300
y2 = -1d300

do i=1,n
do j=1,nfuns
y1 = min(y1,ys(j,i))
y2 = max(y2,ys(j,i))
end do
end do


if (iflogy .eq. 0) then

if (y1 .gt. 0) then
y1 = y1 * 0.98d0
else
y1 = y1*1.02d0
endif

if (y2 .gt. 0) then
y2 = y2 * 1.02d0
else
y2 = y2 * 0.98d0
endif


endif


endif


iw = 1001
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "import warnings"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A)") 'warnings.filterwarnings("ignore")'
write(iw,"(A,I5,A)") "xs  = np.zeros(",n,")"
write(iw,"(A,I5,A,I5,A)") "ys = np.zeros((",nfuns,",",n,"))"

write(iw,"(A)") "plt.rcParams.update({'font.size': 17, 'figure.autolayout': True})"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"]  = ",xs(i)
end do

do i=1,n
do j=1,nfuns
write(iw,"(A,I5,A,I5,A,ES30.18E3)") "ys[",j-1,",",i-1,"] = ",ys(j,i)
end do
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif



!allocate(legend(0),style(0))

idx1 = 1
idx2 = 1

idx3 = 1 
idx4 = 1

do j=1,nfuns


! find the legend string
if (len(legend_loc) .gt. 0) then
do while (legends(idx2:idx2) .ne. '*') 
idx2 = idx2+1
end do
ll = idx2-idx1
allocate(character(ll) :: legend )
legend(1:ll) = legends(idx1:idx2-1)
idx1 = idx2+1
idx2 = idx1
else
allocate(character(0) :: legend )
endif




! find the style string
if (len(styles) .gt. 0) then
do while (styles(idx4:idx4) .ne. '*') 
idx4 = idx4+1
end do
ll = idx4-idx3
allocate(character(ll) :: style )
style(1:ll) = styles(idx3:idx4-1)
idx3 = idx4+1
idx4 = idx3
else
allocate(character(0) :: style  )
endif



write(iw,"(A,I5,A,A,A,A,A)") 'ax.plot(xs,ys[',j-1,',:],"',style,'",label="',legend,'")'

deallocate(style)
deallocate(legend)

end do


if (len(legend_loc) .gt. 0) then
write(iw,"(A,A,A)") 'plt.legend(loc="',legend_loc,'")'
endif

if (iflogx .eq. 1) then
write(iw,"(A)") 'plt.xscale("log")'
elseif (iflogx .eq. 2) then
write(iw,"(A)") 'ax.set_xscale("log", base=2)'
endif

if (iflogy .eq. 1) then
write(iw,"(A)") 'plt.yscale("log")'
elseif (iflogy .eq. 2) then
write(iw,"(A)") 'ax.set_yscale("log", base=2)'
endif

write(iw,"(A,E24.15,A,E25.15,A)") "plt.xlim([",x1,",",x2,"])"
write(iw,"(A,E24.15,A,E25.15,A)") "plt.ylim([",y1,",",y2,"])"
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

close(iw)

allocate(character(7+len(scriptname)) :: command )
write(command,"(A,A)") "python ",scriptname
call system(command)

end subroutine



subroutine plot_function(filename,xlabel,ylabel,n,xs,ys)
implicit double precision (a-h,o-z)

character(len=*)             :: filename,xlabel,ylabel
double precision             :: xs(n),ys(n)
!
!  Use Python and the matplotlib package to produce a plot a function
!  of one variable.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!  product of this code.
!  
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!


iw = 20
open(iw,FILE="plotscript.py")
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"


write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"] = ",xs(i)
write(iw,"(A,I5,A,ES30.18E3)") "ys[",i-1,"] = ",ys(i)
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif

write(iw,"(A)") "ax.plot(xs,ys)"
!write (iw,"(A)") "plt.ylim([-36,1])"
!write (iw,"(A)") "plt.xlim([0,1z])"

write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

! write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
! write(iw,"(A)") 'plt.show()'

close(iw)

call system("python plotscript.py")
! call system("rm -f plotscript.py")

end subroutine


subroutine plot_points(ifshow,filename,n,xs,ys)
implicit double precision (a-h,o-z)

character(len=*)             :: filename
double precision             :: xs(n),ys(n)

!
!  Use Python and the pyplot package to produce a scatter plot.
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!

character(len=12)          :: scriptname
character(len=21)          :: command

call random_number(dd)
ii = dd * (10000)

write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
write(command,"(A,A,I5.5,A)")  "python ","plot",ii,".py &"


iw = 20
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib.pyplot as plt"


write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i)
write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i)
end do

write(iw,"(A)") "fig, ax = plt.subplots()"
write(iw,"(A)") "ax.scatter(xs,ys,s=10)"
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'


if (ifshow .eq. 1) then
write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
write(iw,"(A)") 'plt.show()'
endif

close(iw)

 call system(command)

!call system("python plotscript.py")
!call system("rm -f plotscript.py")

end subroutine


        subroutine quicksort2(n,vals,idxs)
        implicit double precision (a-h,o-z)
        dimension istack(2,20000)
        dimension vals(1),idxs(1)
!
!       Sorts the list of double precision numbers in vals and keep track of
!       indices.
!
        if (n .lt. 100) then
        call insort2(n,vals,idxs)
        return
        endif
!
        maxstack = 10000
        k        = 60
!
        m = 1
        istack(1,1) = 1
        istack(2,1) = n
!
 1000 continue
        if (m .eq. 0) goto 1100
        i1 = istack(1,m)
        i2 = istack(2,m)
        m=m-1
!
        l = i2-i1+1
        if (l .le. k) then
        call insort2(l,vals(i1),idxs(i1))
        goto 1000
        endif
!
!       Otherwise perform quicksort.
!
        call quicksort21(vals,idxs,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
        if (m+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
!
!       Make sure the smaller half is processed first to ensure storage
!       requirements are O(log(n))
!             
        n1 = i3-i1+1
        n2 = i2-i3
!
        if (n2 .lt. n1) then
!
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
!
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
!
        else
!
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
!
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
!
        endif
!
        goto 1000
 1100 continue
        end


        subroutine quicksort21(vals,idxs,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1),idxs(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
        ipiv = i1+(i2-i1)/2
!
        val  = vals(ipiv)
        ival = idxs(ipiv)
!
!       Swap the pivot element and the last element.
!
        vals(ipiv) = vals(i2)
        vals(i2)   = val
!
        idxs(ipiv) = idxs(i2)
        idxs(i2)   = ival
!
        i3 = i1
!
        do 1000 i=i1,i2-1
!
        if( vals(i) .lt. val) then
        d  = vals(i)
        id = idxs(i)
!
        vals(i)  = vals(i3)
        vals(i3) = d       
!
        idxs(i)  = idxs(i3)
        idxs(i3) = id
        i3=i3+1
        endif
 1000 continue
!
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
!
        idd = idxs(i3)
        idxs(i3) = idxs(i2)
        idxs(i2) = idd
!
        end


        subroutine insort2(k,a,ib)
        implicit double precision (a-h,o-z)
        dimension a(1),ib(1)
        if (k .le. 1) return
        do 1000 i=2,k
        val=a(i)
        ival=ib(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val) 
        a(j+1)=a(j)
        ib(j+1)=ib(j)
        j=j-1
 1100 continue
        a(j+1)=val
        ib(j+1)=ival
 1000 continue
        end






! subroutine equispaced_intervals(nints,a,b,ab)
! implicit double precision (a-h,o-z)

! integer, intent(in)            :: nints
! double precision, intent(in)   :: a,b
! double precision, allocatable, intent(out) :: ab(:,:)

! allocate(ab(2,nints))
! do int=1,nints
! ab(1,int) = a + (b-a) * (int-1.0d0)/(nints)
! ab(2,int) = a + (b-a) * (int+0.0d0)/(nints)
! end do

! end subroutine


! subroutine bisected_intervals(nints,a,b,ab)
! implicit double precision (a-h,o-z)

! integer, intent(in)            :: nints
! double precision, intent(in)   :: a,b
! double precision, allocatable, intent(out) :: ab(:,:)

! allocate(ab(2,nints))

! do int=1,nints
! ab(1,int) = 2.0d0**(-nints+int-1) 
! ab(2,int) = 2.0d0**(-nints+int)
! end do

! ab = a + (b-a)*ab
! end subroutine




end module
