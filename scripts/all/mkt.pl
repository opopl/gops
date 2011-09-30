#!/usr/bin/perl 
# header {{{
#===============================================================================
#
#         FILE:  mkt.pl
#
#        USAGE:  ./mkt.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  YOUR NAME (), 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  30/09/11 21:41:28
#     REVISION:  ---
#===============================================================================
#}}}
# use ... dhelp() {{{
use strict;
use warnings;

if ( !@ARGV ){ 
	dhelp();
	exit 1;
}

dhelp(){
#{{{
print << "HELP";
=========================================================
PURPOSE: Generate targets for a Fortran project
USAGE: 
	$this_script FILE
		FILE is output file with dependencies 
SCRIPT LOCATION:
	$0
=========================================================

HELP
#}}}
}
#}}}
# vars {{{
my %libs=(
  	"CONNECT" 	=> "libnc.a",
	"NEB"		=> "libnn.a",
	"AMH"		=> "libamh.a",
	"Bowman"	=> "libbowman.a"
);

my @libdirs=keys %libs;
#}}}

foreach my $ldir (@libdirs){

print << "eof";
cd $ldir
make FC=\${FC} INCPATH=$INCPATH FFLAGS=$FFLAGS
eof

}
