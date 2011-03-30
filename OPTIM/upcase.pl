#/usr/bin/perl -w

eval 'exec /usr/bin/perl -S $0 ${1+"$@"}'
if 0; #$running_under_some_shell

use strict;
use File::Find ();

# Set the variable $File::Find::dont_use_nlink if you're using AFS,
# since AFS cheats.

# for the convenience of &wanted calls, including -eval statements:
use vars qw/*name *dir *prune/;
*name   = *File::Find::name;
*dir    = *File::Find::dir;
*prune  = *File::Find::prune;

sub wanted;

# Traverse desired filesystems
File::Find::find({wanted => \&wanted}, '.');
#exit;

my @ffiles;
my $f;

sub wanted {
    /^.*\.(f90|f|F)\z/s && push(@ffiles,"$name");
}

foreach $f (@ffiles) {
system("cat $f | tr '[a-z]' '[A-Z]' > n; mv n $f")  
}

