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

my $f=$ARGV[0];
open(F,"<$f") or die $!; 

my($dfz,$f,$ext);

while(<F>){
	split; 
	$dfz=$_[6]; 
	$f=$_[0];
	if ( $f != 0 ){
		$ext=abs($dfz/$f);
		$ext= sprintf("%10.5e", $ext);
		$f= sprintf("%10.5e", $f);
		print "$f $ext\n" ;
	}
}
close F;
