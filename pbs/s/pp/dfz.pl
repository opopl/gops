#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  dfz.pl
#
#        USAGE:  ./dfz.pl  
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
#      CREATED:  09/28/11 16:32:36
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my($dfz,$f,$ext,$force);

$f=$ARGV[0];
open(F,"<$f") or die $!; 

while(<F>){
	split; 
	$dfz=$_[6]; 
	$force=$_[0];
	if ( $force != 0  ){
		$ext=abs($dfz/$force);
		$ext= sprintf("%10.5e", $ext);
		$force= sprintf("%10.5e", $force);
		print "$force $ext\n" ;
	}
}
close F;
