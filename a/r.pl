#!/usr/bin/perl 

use strict;
use warnings;
use boolean;

my %rex=(
# {{{
	"dr\(i,i+1,([:0-9^,]+)\)"		=>	"bvr(i,\1)",
	"dr\(i+1,i+2,([:0-9^,]+)\)"		=>	"bvr(i+1,\1)",
	"dr\(i-1,i,([:0-9^,]+)\)"		=>	"bvr(i-1,\1)",
	"dr\(i-2,i-1,([:0-9^,]+)\)"		=>	"bvr(i-2,\1)",
	"dot_prod"						=>	"dpd",
	"x_prod"						=>	"xpd_2",
	"bond_angle\(([^()]*)\)"		=>	"ang(\1,1)",
	"tor_angle\(([^()]*)\)"			=>	"ang(\1,2)",
	"fba_x\([^()]*\)"				=>	"fba(\1,1)",
	"fba_y\([^()]*\)"				=>	"fba(\1,2)",
	"fba_z\([^()]*\)"				=>	"fba(\1,3)",
	"fta_x\([^()]*\)"				=>	"fta(\1,1)",
	"fta_y\([^()]*\)"				=>	"fta(\1,2)",
	"fta_z\([^()]*\)"				=>	"fta(\1,3)",
	"<++>"							=>	"<++>"
	"<++>"							=>	"<++>"
# }}}
);

# Isolate and extract IF ELSEIF statements
~= /(if|elseif)\([^\(\)]*\)//i;
