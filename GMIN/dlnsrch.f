C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,cockup)
      INTEGER n
      LOGICAL check,cockup
      DOUBLE PRECISION f,fold,stpmax,g(n),p(n),x(n),xold(n),ALF
     *,TOLX,grad(n)
      PARAMETER (ALF=1.d-4,TOLX=1.d-7)
      INTEGER i
      DOUBLE PRECISION a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2
     *,slope,sum,temp,
     *test,tmplam

      PRINT '(A)',' Please supply the lnsrch routine!'
      STOP

      END
