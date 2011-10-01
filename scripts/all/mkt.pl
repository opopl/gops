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
use File::Basename;

my $this_script=&basename($0);
my $shd=&dirname($0);

if ( !@ARGV ){ 
	&dhelp();
	exit 1;
}
#}}}
# vars {{{

# fortran compiler
my $FC;
# include path
my $INCPATH;

my %libs=(
  	"CONNECT" 	=> "libnc.a",
	"NEB"		=> "libnn.a",
	"AMH"		=> "libamh.a",
	"Bowman"	=> "libbowman.a"
);

my @libdirs=keys %libs;
#}}}
# subs {{{
# gettime () {{{

sub gettime(){
my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $time = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
return $time;
}

#}}}
# dhelp() {{{
sub dhelp(){
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
}
#}}}
# }}}

my $time=&gettime();
print << "head";
# Make include file with targets
# Time: $time
head

foreach my $ldir (@libdirs){

print << "eof";
$libs{$ldir}:
	cd $ldir
	# Firstly, specify the compiler
	make $FC
	# Then, the actual making
	make INCPATH=$INCPATH
eof

}
