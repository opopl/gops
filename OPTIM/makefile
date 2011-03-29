
# GPL License Info {{{
# 
#   OPTIM: A program for optimizing geometries and calculating reaction pathways
#   Copyright (C) 1999-2011 David J. Wales
#   This file is part of OPTIM.
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
# }}}

#

#ROOTPATH=$(PWD)/../

PPATH=$(PWD)
include $(PPATH)/make.inc

# object files from source files generated during compilation

OBJSrt := porfuncs.o display_version.o
OBJS+= $(OBJSrt)

VPATH = .:AMH:NEB:CONNECT:.

LDFLAGS = -L.

LIBS = -lnn  -lnc $(LLIBS) -lnc

include $(INCPATH)/targets
include $(INCPATH)/mkf

libnn.a: SAT-Ghost
	cd NEB; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}"
libnc.a: SAT-Ghost
	cd CONNECT; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}"
libamh.a: SAT-Ghost
	cd AMH; make FC="${FC}" FFLAGS="${FFLAGS} ${SEARCH_PATH}" 

$(PROG): libnn.a libnc.a $(LLIBS)

libnn.a: modcharmm.o key.o commons.o modunres.o porfuncs.o modmec.o modguess.o efol.o growstring.o gsdata.o intcommons.o intcoords.o amhglobals.o 
libnc.a: libnn.a key.o syminf.o modhess.o amhglobals.o
libamh.a: amhglobals.o

include deps
