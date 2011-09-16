#!/usr/bin/perl
#
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
#}}}

# from find2perl {{{
# Set the variable $File::Find::dont_use_nlink if you're using AFS,
# since AFS cheats.

# for the convenience of &wanted calls, including -eval statements:
use vars qw/*name *dir *prune/;
*name   = *File::Find::name;
*dir    = *File::Find::dir;
*prune  = *File::Find::prune;

sub wanted;
sub get_unused;

sub get_unused{
	#{{{
	my @nused;
	$nused_file="nu.mk";
 if ( open(NUF, $nused_file) ) {
	@nused=<NUF>;
	foreach(@nused) { s/^\s+//; s/\s+$//; }
	chomp(@nused);
	close(NUF);
	}
#}}}
}

sub wanted {
	#{{{
    my ($dev,$ino,$mode,$nlink,$uid,$gid);

    if ( ( /^.*\.(f90|f|F)\z/s) &&
    ( $nlink || (($dev,$ino,$mode,$nlink,$uid,$gid) = lstat($_)) ) &&
    ( ! /^.*\.(save|o|old|ref)\..*\z/s )
	) {
		$name =~ s/^\.\///; 
		if ( ! grep { $_ eq $name } @nused ) {
				push(@fortranfiles,"$name"); 
		}	
	}
#}}}
}

#}}}

#here-doc{{{
my %opt;
@ARGV > 0 and getopts('n:s:m:', \%opt) and not (keys %opt > 1) or die 
+<< "USAGE";
=========================================================
PURPOSE: Generate dependencies for a Fortran project
USAGE: $0 FILE
=========================================================
USAGE
#}}}

my $depsfile="$ARGV[0]";
open(MAKEFILE, ">$depsfile") or die $!; 

get_unused;
File::Find::find({wanted => \&wanted}, '.')  ;
#
# Allow Fortran 90 module source files to have extensions other than .f90
#
@f90 = uniq(F,f,f90);

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
my @libdirs=keys %libs;

foreach (@f90) { s/^/*./ };
# Dependency listings

# Time stuff  {{{

@months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
@weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
$year = 1900 + $yearOffset;
$theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

print MAKEFILE "# $depsfile \n";
print MAKEFILE "# Fortran dependency file \n";
print MAKEFILE "# Created: $theTime \n";

#}}}

&MakeDependsf90();
#&MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

#system("sort $depsfile > n; mv n $depsfile");

#Subroutines  {{{

# &PrintWords(current output column, extra tab?, word list); --- print words nicely
sub PrintWords {
  # {{{
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
	 print MAKEFILE " $word";
	 $columns -= $wordlength + 1;
	 }
      else {
	 #
	 # Continue onto a new line
	 #
	 if ($extratab) {
	    print MAKEFILE " \\\n\t\t$word";
	    $columns = 62 - $wordlength;
	    }
	 else {
	    print MAKEFILE " \\\n\t$word";
	    $columns = 70 - $wordlength;
	    }
	 }
      }
# }}}
}

#Subs: LanguageCompiler toLower uniq {{{

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
   local($string) = @_[0];
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
      if ($word ne $words[$#words]) {
	 push(@words, $word);
	 }
      }
   @words;
#}}}
}

#}}}

# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
# {{{
  my $subname = (caller(0))[3];
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
	 /$pattern/i && push(@incs, $1);
	 }
      if (defined @incs) {
	 $file =~ s/\.[^.]+$/.o/;
	 print MAKEFILE "$file: ";
	 &PrintWords(length($file) + 2, 0, @incs);
	 print MAKEFILE "\n";
	 undef @incs;
	 }
      }
# }}}
}

# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
sub MakeDependsf90 {
#{{{
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
      while (<FILE>) {
	 if ( $_ =~ /^\s*include\s+["\']([^"\']+)["\']/i ){ 
	   if ( !exists $excluded{$1} ){
     	     push(@incs, $1); 
	   }
	 }
	 /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
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
	 print MAKEFILE "$objfile: ";
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
	 print MAKEFILE "\n";
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
