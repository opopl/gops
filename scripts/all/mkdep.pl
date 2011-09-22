#!/usr/bin/perl -w

# start {{{
# Changelog:
#
# Original makemake utility - Written by Michael Wester <wester@math.unm.edu> December 27, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# 14:19:31 (Sat, 26-Mar-2011):
#
# mkdep - put under git control by op
#

eval 'exec /usr/bin/perl -S $0 ${1+"$@"}'
if 0; #$running_under_some_shell

#use strict;
use File::Find ();
use Getopt::Std;
use Cwd;
use File::Basename;

# Set the variable $File::Find::dont_use_nlink if you're using AFS, since AFS cheats.
# for the convenience of &wanted calls, including -eval statements:
use vars qw/*name *dir *prune/;
*name   = *File::Find::name;
*dir    = *File::Find::dir;
*prune  = *File::Find::prune;

#}}}
# vars {{{

# %excluded %libs {{{
my %excluded=(
	"mpif.h" =>	"",
	"MPIF.H" =>	""
);

my %libs=(
  	"CONNECT" 	=> "libnc.a",
	"NEB"		=> "libnn.a",
	"AMH"		=> "libamh.a",
	"libbowman.a"	=> "libbowman.a"
      );
  #}}}

my @libdirs=keys %libs;

# array for not-used fortran files which are taken from file nu.mk
my @nused;
# does the nu.mk file exist? 0 for no, 1 for yes
my $nu_exist;

my $nu_dir;

# current project full path, e.g., /home/op226/gops/G
my $PPATH=&cwd();
# current program name, e.g., G
my $PROGNAME=&basename($PPATH);
# scripts/all/
my $SAPATH=&dirname($0);
# $HOME/gops/
my $ROOTPATH="$SAPATH/../../";
my $INCPATH="$ROOTPATH/include/";
my $F_NU="$INCPATH/nu_$PROGNAME.mk";
my $F_DP="$PPATH/$ARGV[0]";
open(DP, ">$F_DP") or die $!; 
print DP "# Project dir: $PPATH\n";
print DP "# Program name: $PROGNAME\n";

#}}}
# here-doc{{{
my %opt;
@ARGV > 0 and getopts('n:s:m:', \%opt) and not (keys %opt > 1) or die 
+<< "USAGE";
=========================================================
PURPOSE: Generate dependencies for a Fortran project
USAGE: $0 FILE
=========================================================
USAGE
#}}}
# Time stuff  {{{

my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
my $year = 1900 + $yearOffset;
my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

#}}}
# Subroutines  {{{

# wanted() get_unused() {{{

sub wanted;
sub get_unused;

sub get_unused{
	#{{{
	$nu_exist=1;
  open(NUF, $F_NU) or $nu_exist=0; 
  if ($nu_exist eq 1){
	  while(<NUF>){
		chomp;
		if( ! /^(#|d:)/ ){
			push(@nused,$_);
		}
		elsif( /^d:(\w+)/g ){
			push(@nu_dirs,$1);
			print DP "# Not used dir: $1\n";
		}
}
	foreach(@nused) { s/^\s+//; s/\s+$//; }
	close(NUF);
	}
# }}}
}

sub wanted {
	#{{{
    my ($dev,$ino,$mode,$nlink,$uid,$gid);
	( my $dirname = $File::Find::dir ) =~ s/$PPATH//g;
	$dirname =~ s/^(\.|\/)+//g;

    if ( ( /^.*\.(f90|f|F)\z/s) &&
    ( $nlink || (($dev,$ino,$mode,$nlink,$uid,$gid) = lstat($_)) ) &&
    ( ! /^.*\.(save|o|old|ref)\..*\z/s ) &&
	( ! grep { $_ eq $dirname } @nu_dirs )
	) {
		$name =~ s/^\.\///; 
		if ( ! grep { $_ eq $name } @nused ) {
				push(@fortranfiles,"$name"); 
			}
	}
#}}}
}

#}}}
#Subs: MakeDepends PrintWords LanguageCompiler toLower uniq {{{

sub MakeDepends{
# {{{
# &MakeDepends(language pattern, include file sed pattern); --- dependency maker
  my $subname = (caller(0))[3];
   local(@incs);
   local($lang) = $_[0];
   local($pattern) = $_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
	 /$pattern/i && push(@incs, $1);
	 }
      if (defined @incs) {
	 $file =~ s/\.[^.]+$/.o/;
	 print DP "$file: ";
	 &PrintWords(length($file) + 2, 0, @incs);
	 print DP "\n";
	 undef @incs;
	 }
      }
# }}}
}

sub PrintWords {
  # {{{
# &PrintWords(current output column, extra tab?, word list); --- print words nicely
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print DP $_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
	 print DP " $word";
	 $columns -= $wordlength + 1;
	 }
      else {
	 #
	 # Continue onto a new line
	 #
	 if ($extratab) {
	    print DP " \\\n\t\t$word";
	    $columns = 62 - $wordlength;
	    }
	 else {
	    print DP " \\\n\t$word";
	    $columns = 70 - $wordlength;
	    }
	 }
      }
# }}}
}

# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
sub LanguageCompiler {
  # {{{
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
	 grep(/^$compiler$/, ("fc", "f77")) &&
	    do { $compiler = "FC"; last CASE; };
	 grep(/^$compiler$/, ("cc", "c"))   &&
	    do { $compiler = "CC"; last CASE; };
	 $compiler = "F90";
	 }
      }
   else {
      CASE: {
	 grep(/\.f90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
	 grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
	 grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
	 $compiler = "???";
	 }
      }
   $compiler;
 #}}}
}

# &toLower(string); --- convert string into lower case
sub toLower {
   local($string) = $_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
}

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
# {{{
   local(@words);
   foreach $word (@_) {
      #if ( ( (defined($word)) && ($word ne $words[$#words])) || ( $#words==0 )) {
      if ( $word ne $words[$#words]) {
	 	push(@words, $word);
	 	}
      }
   @words;
#}}}
}

#}}}

sub MakeDependsf90 {
#{{{
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
  my $subname = (caller(0))[3];
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it {{{
   
   foreach $file (@fortranfiles) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      # get extension from the $file
       if ( $file =~ /.*\.(f|F|f90)$/ ){
	 $fext=$1;
       }
       # get the object name for the module
      while (<FILE>) {
	 /^\s*module\s+([^\s!]+)/i &&
	    ($filename{&toLower($1)} = $file) =~ s/\.$fext$/.o/;
	 }
      }
   # }}}
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   # {{{
   foreach $file (@fortranfiles) {
      open(FILE, $file);
	  @flines=<FILE>;
	  #get_use();
      foreach (@flines) {
		#{{{
	  		chomp;	  
	  		#include files {{{	  
	 		if ( $_ =~ /^\s*include\s+["\']([^"\']+)["\']/i ){ 
	   			if ( !exists $excluded{$1} ){
					# included files themselves not printed 
					# in the DP file
                    # push(@incs, $1); 
			 		# now use those statements in that file
			 		# specified after the include ... statement
			 		my $include_file=$1;
			 		open(IFILE,$include_file);
			 		while(<IFILE>){
	 					/^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
			 		}
					close IFILE;
	   			}
	 		}
     		#}}}
	 	/^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
 	 	#}}}
	 }
      if (defined @incs || defined @modules) {
	 ($objfile = $file) =~ s/\.[^.]+$/.o/;
	 # for each object file, print relevant dependencies 
	 # ONLY those object files which are in the root directory
	 if ( ( $objfile !~ /^(.*)\// ) 
		 && ( $objfile !~ /.*\.(inc|i)\..*/ ) # don't print include files
		 && ( $objfile !~ /.*\.(old|o)\..*/ ) # don't print old files
		 && ( $objfile !~ /.*\.(save)\..*/ ) # don't use save files
		 && ( $objfile !~ /.*\.(ref)\..*/ ) # don't print "reference" files, i.e. files included from other code
	 ) {
	 print DP "$objfile: ";
	 undef @dependencies;
	 foreach $module (@modules) {
	   if (defined $filename{$module})
	   {
	    $mo=$filename{$module};
	  } else
	  {
		 ( $mo=$module ) =~ s/$/\.o/;
	  } 
	    foreach $libdir (@libdirs){
	        $libname=$libs{$libdir};
		$mo =~ s/^$libdir\/.*/$libname/;
	    }
		    if ( ! exists $excluded{$mo} ){
		      #print "$mo \n";
			    push(@dependencies, $mo);
		    }
	  }
	 @dependencies = &uniq(sort(@dependencies));
	 &PrintWords(length($objfile) + 2, 0,
		     @dependencies, &uniq(sort(@incs)));
	 print DP "\n";
	 }
	 #
	 undef @incs;
	 undef @modules;
	       }
      }
    #}}}
# }}}
}

# }}}
# main {{{

&get_unused;
File::Find::find({wanted => \&wanted}, '.')  ;
@f90 = uniq("F","f","f90");
foreach (@f90) { s/^/*./ };

print DP << "head";
# $F_DP
# Fortran dependency file
# Created: $theTime
head

&MakeDependsf90();

close DP;
#}}}
