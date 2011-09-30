#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  cpge.pl
#
#        USAGE:  ./cpge.pl  
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
#      CREATED:  09/29/11 21:39:32
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use File::Copy;

my ($jn,$gf,$d,@f);
$gf="g.sh";
$d="/home/op226/data/28-09-11/";

open(F,"<$gf");

while(<F>){
	@f=split ;
	$jn=$f[5];

}
close F;
