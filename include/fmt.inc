
! formats and strings                                                            {{{

! strings {{{
character(40) s_stars
s_stars="*****************************************"
! }}}
! 1 - string
1 format(a40,2a20,a10)
10 format(a)
11 format(a40)
12 format(a,a20)
120 format(a,a20,a100)
121 format(a,a20)
13 format(a20,a40)

14 format(3a) 
15 format(a20,a50) 
                ! 2 - real
2 format(a40,a20,f20.5)
21 format(a40,a20,e20.1)
                ! 3 - integer
3 format(a40,a20,i20)
                ! 4 - logical
4 format(a40,a20,l20)
40 format(a40,a20,l20,i10)
41 format(a20,l20)
! }}}

