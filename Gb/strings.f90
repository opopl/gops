
MODULE STRINGS

USE PRECISION

private :: value_dr,value_sr,value_di,value_si
private :: write_dr,write_sr,write_di,write_si
private :: writeq_dr,writeq_sr,writeq_di,writeq_si

interface value  ! Generic operator for converting a number string to a 
                 ! number. Calling syntax is 'call value(numstring,number,ios)' 
                 ! where 'numstring' is a number string and 'number' is a 
                 ! real number or an integer (single or double precision).         
   module procedure value_dr
   module procedure value_sr
   module procedure value_di
   module procedure value_si
end interface

interface writenum  ! Generic  interface for writing a number to a string. The 
                    ! number is left justified in the string. The calling syntax
                    ! is 'call writenum(number,string,format)' where 'number' is
                    ! a real number or an integer, 'string' is a character string
                    ! containing the result, and 'format' is the format desired, 
                    ! e.g., 'e15.6' or 'i5'.
   module procedure write_dr
   module procedure write_sr
   module procedure write_di
   module procedure write_si
end interface

interface writeq  ! Generic interface equating a name to a numerical value. The
                  ! calling syntax is 'call writeq(unit,name,value,format)' where
                  ! unit is the integer output unit number, 'name' is the variable
                  ! name, 'value' is the real or integer value of the variable, 
                  ! and 'format' is the format of the value. The result written to
                  ! the output unit has the form <name> = <value>.
   module procedure writeq_dr
   module procedure writeq_sr
   module procedure writeq_di
   module procedure writeq_si
end interface


!**********************************************************************

contains

!**********************************************************************

SUBROUTINE PARSE(STR,DELIMS,ARGS,NARGS)
! {{{
! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

! sub
CHARACTER(LEN=*) :: STR,DELIMS
CHARACTER(LEN=*),DIMENSION(:) :: ARGS
INTEGER NARGS

! loc
CHARACTER(LEN=LEN_TRIM(STR)) :: STRSAV
INTEGER i,NA,K,LENSTR

STRSAV=STR
CALL COMPACT(STR)
NA=SIZE(ARGS)
DO I=1,NA
  ARGS(I)=' '
END DO  
NARGS=0
LENSTR=LEN_TRIM(STR)
IF(LENSTR==0) RETURN
K=0

DO
   IF(LEN_TRIM(STR) == 0) EXIT
   NARGS=NARGS+1
   CALL SPLIT(STR,DELIMS,ARGS(NARGS))
   CALL REMOVEBKSL(ARGS(NARGS))
END DO   
STR=STRSAV
! }}}
END SUBROUTINE PARSE

!**********************************************************************

SUBROUTINE COMPACT(STR)
! {{{
! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

! sub
CHARACTER(LEN=*):: STR

! loc
CHARACTER(LEN=1):: CH
CHARACTER(LEN=LEN_TRIM(STR)):: OUTSTR

INTEGER LENSTR,ISP,K,ICH,i

STR=ADJUSTL(STR)
LENSTR=LEN_TRIM(STR)
OUTSTR=' '
ISP=0
K=0

DO I=1,LENSTR
  CH=STR(I:I)
  ICH=IACHAR(CH)
  
  SELECT CASE(ICH)
  
    CASE(9,32)     ! SPACE OR TAB CHARACTER
      IF(ISP==0) THEN
        K=K+1
        OUTSTR(K:K)=' '
      END IF
      ISP=1
      
    CASE(33:)      ! NOT A SPACE, QUOTE, OR CONTROL CHARACTER
      K=K+1
      OUTSTR(K:K)=CH
      ISP=0
      
  END SELECT
  
END DO

STR=ADJUSTL(OUTSTR)
! }}}
END SUBROUTINE COMPACT

!**********************************************************************

subroutine removesp(str)

! Removes spaces, tabs, and control characters in string str

! sub
character(len=*):: str

! loc
character(len=1):: ch
character(len=len_trim(str))::outstr
integer lenstr,k,i,ich

str=adjustl(str)
lenstr=len_trim(str)
outstr=' '
k=0

do i=1,lenstr
  ch=str(i:i)
  ich=iachar(ch)
  select case(ich)    
    case(0:32)  ! space, tab, or control character
         cycle       
    case(33:)  
      k=k+1
      outstr(k:k)=ch
  end select
end do

str=adjustl(outstr)

end subroutine removesp

!**********************************************************************

subroutine value_dr(str,rnum,ios)

! Converts number string to a double precision real number

! sub
character(len=*)::str
real(kr8)::rnum
integer :: ios
! loc
INTEGER ::  ilen,ipos

ilen=len_trim(str)
ipos=scan(str,'Ee')
if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
   ios=3
   return
end if
read(str,*,iostat=ios) rnum

end subroutine value_dr

!**********************************************************************

subroutine value_sr(str,rnum,ios)

! Converts number string to a single precision real number

! sub
character(len=*)::str
real(kr4) :: rnum
integer :: ios

! local
real(kr8) :: rnumd 

call value_dr(str,rnumd,ios)
if( abs(rnumd) > huge(rnum) ) then
  ios=15
  return
end if
if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
rnum=rnumd

end subroutine value_sr

!**********************************************************************

subroutine value_di(str,inum,ios)

! Converts number string to a double precision integer value

character(len=*)::str
integer(ki8) :: inum
integer ios


real(kr8) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki8)

end subroutine value_di

!**********************************************************************

subroutine value_si(str,inum,ios)

! Converts number string to a single precision integer value

character(len=*)::str
integer(ki4) :: inum
integer :: ios
real(kr8) :: rnum

call value_dr(str,rnum,ios)
if(abs(rnum)>huge(inum)) then
  ios=15
  return
end if
inum=nint(rnum,ki4)

end subroutine value_si

!**********************************************************************

subroutine shiftstr(str,n)
 
! Shifts characters in in the string 'str' n positions (positive values
! denote a right shift and negative values denote a left shift). Characters
! that are shifted off the end are lost. Positions opened up by the shift 
! are replaced by spaces.

character(len=*):: str
integer, intent(in) :: n
integer nabs,lenstr

lenstr=len(str)
nabs=iabs(n)
if(nabs>=lenstr) then
  str=repeat(' ',lenstr)
  return
end if
if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
return

end subroutine shiftstr

!**********************************************************************

subroutine insertstr(str,strins,loc)

! Inserts the string 'strins' into the string 'str' at position 'loc'. 
! Characters in 'str' starting at position 'loc' are shifted right to
! make room for the inserted string. Trailing spaces of 'strins' are 
! removed prior to insertion

character(len=*):: str,strins
character(len=len(str))::tempstr
integer loc,lenstrins

lenstrins=len_trim(strins)
tempstr=str(loc:)
call shiftstr(tempstr,lenstrins)
tempstr(1:lenstrins)=strins(1:lenstrins)
str(loc:)=tempstr
return

end subroutine insertstr

!**********************************************************************

subroutine delsubstr(str,substr)

! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.

! sub
character(len=*):: str,substr
! loc
integer lensubstr,ipos

lensubstr=len_trim(substr)
ipos=index(str,substr)
if(ipos==0) return
if(ipos == 1) then
   str=str(lensubstr+1:)
else
   str=str(:ipos-1)//str(ipos+lensubstr:)
end if   
return

end subroutine delsubstr

!**********************************************************************

subroutine delall(str,substr)

! Deletes all occurrences of substring 'substr' from string 'str' and
! shifts characters left to fill holes.

! sub
character(len=*):: str,substr

! loc
integer lensubstr,ipos

lensubstr=len_trim(substr)
do
   ipos=index(str,substr)
   if(ipos == 0) exit
   if(ipos == 1) then
      str=str(lensubstr+1:)
   else
      str=str(:ipos-1)//str(ipos+lensubstr:)
   end if
end do   
return

end subroutine delall

!**********************************************************************

function uppercase(str) result(ucstr)

! convert string to upper case

! sub
character (len=*):: str
character (len=len_trim(str)):: ucstr

! loc
integer ilen,ioffset,iquote,i,iav,iqc

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')     
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do
return

end function uppercase

!**********************************************************************

function lowercase(str) result(lcstr)

! convert string to lower case

! sub
character (len=*):: str
character (len=len_trim(str)):: lcstr

! local
integer ilen,ioffset,iquote,i,iav,iqc

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do
return

end function lowercase

!**********************************************************************

subroutine readline(nunitr,line,ios)

! Reads line from unit=nunitr, ignoring blank lines
! and deleting comments beginning with an exclamation point(!)

! sub
character (len=*):: line
integer :: ios,nunitr

! local 
integer ipos

do  
  read(nunitr,'(a)', iostat=ios) line      ! read input line
  if(ios /= 0) return
  line=adjustl(line)
  ipos=index(line,'!')
  if(ipos == 1) cycle
  if(ipos /= 0) line=line(:ipos-1)
  if(len_trim(line) /= 0) exit
end do
return

end subroutine readline

!**********************************************************************

subroutine match(str,ipos,imatch)

! Sets imatch to the position in string of the delimiter matching the delimiter
! in position ipos. Allowable delimiters are (), [], {}, <>.

! subroutine
character(len=*) :: str
integer, intent(in) :: ipos
integer, intent(out) :: imatch

! local 
character :: delim1,delim2,ch
integer :: lenstr,i,idelim2,iend,istart,inc,isum

lenstr=len_trim(str)
delim1=str(ipos:ipos)
select case(delim1)
!{{{
   case('(')
      idelim2=iachar(delim1)+1
      istart=ipos+1
      iend=lenstr
      inc=1
   case(')')
      idelim2=iachar(delim1)-1
      istart=ipos-1
      iend=1
      inc=-1
   case('[','{','<')
      idelim2=iachar(delim1)+2
      istart=ipos+1
      iend=lenstr
      inc=1
   case(']','}','>')
      idelim2=iachar(delim1)-2
      istart=ipos-1
      iend=1
      inc=-1
   case default
      write(*,*) delim1,' is not a valid delimiter'
      return
      ! }}}
end select
if(istart < 1 .or. istart > lenstr) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if
delim2=achar(idelim2) ! matching delimiter

isum=1
do i=istart,iend,inc
   ch=str(i:i)
   if(ch /= delim1 .and. ch /= delim2) cycle
   if(ch == delim1) isum=isum+1
   if(ch == delim2) isum=isum-1
   if(isum == 0) exit
end do
if(isum /= 0) then
   write(*,*) delim1,' has no matching delimiter'
   return
end if   
imatch=i

return

end subroutine match

!**********************************************************************

subroutine write_dr(rnum,str,fmt)

! Writes double precision real number rnum to string str using format fmt

real(kr8) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_dr

!***********************************************************************

subroutine write_sr(rnum,str,fmt)

! Writes single precision real number rnum to string str using format fmt

real(kr4) :: rnum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) rnum
str=adjustl(str)

end subroutine write_sr

!***********************************************************************

subroutine write_di(inum,str,fmt)

! Writes double precision integer inum to string str using format fmt

integer(ki8) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_di

!***********************************************************************

subroutine write_si(inum,str,fmt)

! Writes single precision integer inum to string str using format fmt

integer(ki4) :: inum
character(len=*) :: str,fmt
character(len=80) :: formt

formt='('//trim(fmt)//')'
write(str,formt) inum
str=adjustl(str)

end subroutine write_si

!***********************************************************************

subroutine trimzero(str)

! Deletes nonsignificant trailing zeroes from number string str. If number
! string ends in a decimal point, one trailing zero is added.

character(len=*) :: str

! local
character :: ch
character(len=10) :: exp
integer i,ipos,lstr

ipos=scan(str,'eE')
if(ipos>0) then
   exp=str(ipos:)
   str=str(1:ipos-1)
endif
lstr=len_trim(str)
do i=lstr,1,-1
   ch=str(i:i)
   if(ch=='0') cycle          
   if(ch=='.') then
      str=str(1:i)//'0'
      if(ipos>0) str=trim(str)//trim(exp)
      exit
   endif
   str=str(1:i)
   exit
end do
if(ipos>0) str=trim(str)//trim(exp)

end subroutine trimzero

!**********************************************************************

subroutine writeq_dr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr8) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_dr

!**********************************************************************

subroutine writeq_sr(unit,namestr,value,fmt)

! Writes a string of the form <name> = value to unit

real(kr4) :: value
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr

call writenum(value,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_sr

!**********************************************************************

subroutine writeq_di(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki8) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_di

!**********************************************************************

subroutine writeq_si(unit,namestr,ivalue,fmt)

! Writes a string of the form <name> = ivalue to unit

integer(ki4) :: ivalue
integer :: unit
character(len=*) :: namestr,fmt
character(len=32) :: tempstr
call writenum(ivalue,tempstr,fmt)
call trimzero(tempstr)
write(unit,*) trim(namestr)//' = '//trim(tempstr)

end subroutine writeq_si

!**********************************************************************

FUNCTION IS_LETTER(CH) RESULT(RES)
! {{{
! Returns .true. if ch is a letter and .false. otherwise

CHARACTER :: CH
LOGICAL :: RES

SELECT CASE(CH)
CASE('a':'z','A':'Z')
  RES=.TRUE.
CASE DEFAULT
  RES=.FALSE.
END SELECT
RETURN
! }}}
END FUNCTION IS_LETTER

!**********************************************************************

FUNCTION IS_DIGIT(CH) RESULT(RES)
! {{{
! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

CHARACTER :: CH
LOGICAL :: RES

SELECT CASE(CH)
CASE('0':'9')
  RES=.TRUE.
CASE DEFAULT
  RES=.FALSE.
END SELECT
RETURN
! }}}
END FUNCTION IS_DIGIT

!**********************************************************************

SUBROUTINE SPLIT(STR,DELIMS,BEFORE,SEP)
! {{{
! declarations  {{{

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

! subroutine 
CHARACTER(LEN=*) :: STR,DELIMS,BEFORE
CHARACTER,OPTIONAL :: SEP

! local 
LOGICAL :: PRES
CHARACTER :: CH,CHA
INTEGER LENSTR,K,IBSL,I,IPOS,IPOSA
! }}}


PRES=PRESENT(SEP)
STR=ADJUSTL(STR)
CALL COMPACT(STR)
LENSTR=LEN_TRIM(STR)
IF(LENSTR == 0) RETURN        ! STRING STR IS EMPTY
K=0
IBSL=0                        ! BACKSLASH INITIALLY INACTIVE
BEFORE=' '
DO I=1,LENSTR
   CH=STR(I:I)
   IF(IBSL == 1) THEN          ! BACKSLASH ACTIVE
      K=K+1
      BEFORE(K:K)=CH
      IBSL=0
      CYCLE
   END IF
   IF(CH == '\\') THEN          ! BACKSLASH WITH BACKSLASH INACTIVE
      K=K+1
      BEFORE(K:K)=CH
      IBSL=1
      CYCLE
   END IF
   IPOS=INDEX(DELIMS,CH)         
   IF(IPOS == 0) THEN          ! CHARACTER IS NOT A DELIMITER
      K=K+1
      BEFORE(K:K)=CH
      CYCLE
   END IF
   IF(CH /= ' ') THEN          ! CHARACTER IS A DELIMITER THAT IS NOT A SPACE
      STR=STR(I+1:)
      IF(PRES) SEP=CH
      EXIT
   END IF
   CHA=STR(I+1:I+1)            ! CHARACTER IS A SPACE DELIMITER
   IPOSA=INDEX(DELIMS,CHA)
   IF(IPOSA > 0) THEN          ! NEXT CHARACTER IS A DELIMITER
      STR=STR(I+2:)
      IF(PRES) SEP=CHA
      EXIT
   ELSE
      STR=STR(I+1:)
      IF(PRES) SEP=CH
      EXIT
   END IF
END DO
IF(I >= LENSTR) STR=''
STR=ADJUSTL(STR)              ! REMOVE INITIAL SPACES
RETURN
! }}}
END SUBROUTINE SPLIT

!**********************************************************************

SUBROUTINE REMOVEBKSL(STR)
! {{{
! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

CHARACTER(LEN=*):: STR
CHARACTER(LEN=1):: CH
CHARACTER(LEN=LEN_TRIM(STR))::OUTSTR

INTEGER  K,IBSL,I,LENSTR

STR=ADJUSTL(STR)
LENSTR=LEN_TRIM(STR)
OUTSTR=' '
K=0
IBSL=0                        ! BACKSLASH INITIALLY INACTIVE

DO I=1,LENSTR
  CH=STR(I:I)
  IF(IBSL == 1) THEN          ! BACKSLASH ACTIVE
   K=K+1
   OUTSTR(K:K)=CH
   IBSL=0
   CYCLE
  END IF
  IF(CH == '\\') THEN          ! BACKSLASH WITH BACKSLASH INACTIVE
   IBSL=1
   CYCLE
  END IF
  K=K+1
  OUTSTR(K:K)=CH              ! NON-BACKSLASH WITH BACKSLASH INACTIVE
END DO

STR=ADJUSTL(OUTSTR)
! }}}
END SUBROUTINE REMOVEBKSL

!**********************************************************************

end module strings  


