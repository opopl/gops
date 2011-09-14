#   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
#   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
#   This file is part of CONNECT module. CONNECT module is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

PPATH=$(PWD)/../
INCPATH = $(PPATH)/../include/
include $(INCPATH)/inc.mk

VPATH = ../NEB:..:../CHARMM
LIBNAME = ../libnc.a

.SUFFIXES:
.SUFFIXES: .o .f90 

.f90.o :
	${FC} -c ${FFLAGS} $<

new_connect: ${LIBNAME}

${LIBNAME}: ${OBJS} 
	ar sr ${LIBNAME} ${OBJS}

clean:
	rm -f ${LIBNAME} ${OBJS} *.mod ifc?????? *.il *.stb *.o

keyconnect.o: keyutils.o keydecide.o charutils.o keyneb.o gsdata.o  key.o
decide.o:          porfuncs.o ncdata.o ncutils.o keydecide.o
dodumpsp.o:      ncdata.o
ido.o:      ncdata.o ncutils.o keyconnect.o key.o porfuncs.o keyconnect.o ncutils.o 
keydecide.o: charutils.o
keyutils.o: key.o
ncutils.o:  ncdata.o ntc.o keyutils.o key.o commons.o intcommons.o keydecide.o ntc.o keyconnect.o amhglobals.o syminf.o modhess.o modcharmm.o porfuncs.o modunres.o 
newconnect.o: ncdata.o keyconnect.o ncutils.o modcharmm.o decide.o tryconnect.o ido.o porfuncs.o amhglobals.o key.o modneb.o commons.o
tryconnect.o: newneb.o ncdata.o keyconnect.o ncutils.o ntc.o keyneb.o key.o modguess.o modunres.o modcharmm.o modamber9.o modmec.o keyutils.o commons.o porfuncs.o gsdata.o amhglobals.o
