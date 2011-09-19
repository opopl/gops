
BASE_SOURCE:= $(shell cat base_source.in) 
BLN_SOURCE:= $(shell cat bln_source.in) 

BASE_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(BASE_SOURCE))))
BLN_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(BLN_SOURCE))))
