
BASE_SOURCE:= $(shell cat base_source.in) 
BLN_SOURCE:= p46merdiff.f g46merdiff.f BLN.f

BASE_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(BASE_SOURCE))))
BLN_OBS := $(patsubst %.F,%.o,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(BLN_SOURCE))))
