#!/usr/bin/perl -w

use strict;
use Pod::Usage;
use Getopt::Long;

our ($opt_c, $opt_l);
our ($opt_help, $opt_man);
GetOptions("help", "man", "l=i", "c")
  or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

my $line_len = 72;
if ($opt_l) {
    $line_len = $opt_l;
}
my $convert;
if ($opt_c) {
    $convert = 1;
}
#
# Loop over all statements
#
my $lookahead = <>;
my($left, $parexp, $after);
statement: while ($_ = get_line()) {
#
# Delete trailing spaces
#
    s/ *$//;
#
# Skip comments and blank lines
#
    if (/^[*c#!]|^$/i) {
        print;
        next;
    }
#
# Replace tabs with spaces
#
    s/\t/        /g;
#
# Put space around =
#
    s/\s*=\s*/ = /;
#
# Look for if-then-else and do loops 
#
    my($before, $middle, $begin, $end, $statement);
    if (/^(.{6} *(else)? *if) *\(/i) {
        $_ = $';
        $before = $1;
        $left = 1;
        find_match();
        if ($left != 0) { die "Illegal if statement"; }
    # add spaces inside parentheses
        $statement = '';
        $_ = $parexp;
        while (/( *)\.(gt|lt|eq|ne|le|ge|and|or|not|eqv|neqv)\.( *)/im) {
            $begin = $`;
            $middle = '.'.$2.'.';
            $end = $';
            if ($convert) {
MATCH:      {
                if ($middle =~ /\.gt\./m) { $middle = '>'; last MATCH; }
                if ($middle =~ /\.ge\./m) { $middle = '>='; last MATCH; }
                if ($middle =~ /\.lt\./m) { $middle = '<'; last MATCH; }
                if ($middle =~ /\.le\./m) { $middle = '<='; last MATCH; }
                if ($middle =~ /\.eq\./m) { $middle = '=='; last MATCH; }
                if ($middle =~ /\.ne\./m) { $middle = '/='; last MATCH; }
            } }
            $statement .= $begin.' '.$middle.' ';
            $_ = $end;
        }
    $statement .= $_;
    } else {
        print;
        next;
    }
#
# Put string back together
#
    $after =~ s/^ *//m;
    $_ = $before.' ('.$statement.' '.$after;
#
# Check for falling off the edge of the world
#
    my @lines = split(/^/);
    my $line;
    foreach $line (@lines) {
    if (length($line) > $line_len+1 &&    # +1 for '\n'
            !(substr($_,0,73) =~ /!/)) {      # don't complain about long comments
            print STDERR "Line ", $.-1, " too long. Fix it\n";
        }
#
# Delete trailing spaces, then print
#
        $line =~ s/ *$//;
        print $line;
    }
}
#
# end main, begin subs
# get a line, combining continuation lines
#
sub get_line {
    my $thisline = $lookahead;
    if ($lookahead) {
        line: while ($lookahead = <>) {
            if ($lookahead =~ /^     \S/) {
#            if ($lookahead =~ /^     \S|^$|^[*c#!]/i) {
                $thisline .= $lookahead;
            }
            else {
                last line;
            }
        }
    }
    $thisline;
}
#
# find matching parentheses
#
sub find_match {
    $parexp = '';
    while (/[()]/m) {
        $parexp .= $`;
        $parexp .= $&;
        $_ = $';
        if ($& eq "(") { $left++; }
        else           { $left--; }
        if ($left == 0) { last; }
    }
    $after = $_;
}
__END__

ifspace - Add space around Fortran if statements

=head1 SYNOPSIS

ifspace [--help] [--man] [-l length] [-c] [file ...] > outfile

=head1 DESCRIPTION

Add space around 'if' statements in a Fortran program, for instance
changing:

      if(i.lt.20.or.j.gt.30)then

to

      if (i .lt. 20 .or. j .gt. 30) then

It also puts one space on either side of an equals sign.

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<-help>

Print more details about the arguments.

=item I<-man>

Print a full man page.

=item I<-l length>

Specifies the maximum length of your Fortran lines (default is 72).

=item I<-c>

Convert to f90 comparisons, e.g., ".gt." becomes ">".

=back

=head1 EXAMPLE


=head1 AUTHOR

   Kate Hedstrom
   kate@arsc.edu

=head1 BUGS

       It checks to see if lines are too long, but complains
    to STDERR instead of fixing them.
       It doesn't always do exactly what you want for multi-line if
    statements.  For instance, 

         if ((plot_s .or. plot_t .or. plot_rho .or. plot_curl)
        &           .and. .not. bigpix) call frame

    becomes

         if ((plot_s .or. plot_t .or. plot_rho .or. plot_curl)
        & .and.  .not. bigpix) call frame

=cut
