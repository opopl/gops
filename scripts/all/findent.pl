#!/usr/bin/perl -w

use Pod::Usage;
use Getopt::Long;

our ($opt_help, $opt_man);
our($opt_l, $opt_n, $opt_c);
GetOptions("help", "man", "c", "l=i", "n=i")
  or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

my $indent = 0;
my $section_no = 0;
my $s_pref = $section_no."_";
my $line_len = 72;
my $delta1 = 2;

if ($opt_l) {
    $line_len = $opt_l;
}
if ($opt_n) {
    $delta1 = $opt_n;
}
#
# first, remove all indentation
#
my $lookahead = <>;
statement: while ($_ = &get_line()) {
#
# Delete trailing spaces
#
    s/ *$//;
#
# Skip comments and blank lines
#
    if (/^[*c#!]|^$/i) {
        if ($opt_c) { s/^[*c!] */!  /; }
        print;
        next;
    }
#
# Replace tabs with spaces
#
    s/\t/        /g;
#
# Check for new section (function or subroutine)
#
    if (/function|subroutine/is && $` !~ /'/) {
        $section_no++;
        $s_pref = $section_no."_";
    }
#
# Remove current indentation
#
    $match = /^(.{6}) */;
    if ($match) {
        $before = $1;
        $after = $';
    } else {
        die "Illegal Fortran statement ";
    }
#
# Look for if-then-else and do loops 
#
    $outy = 0;
    $_ = $after;
    if (/^ *if *\(/i) {
        $_ = $';
        $parexp = &find_match(1);
        if (/^ *then|^ *\n     \S *then/is) {
            $indent++;
            $outy = 1;
        }
    }
    elsif (/^ *else/i) {
        $outy = 1;
    }
    elsif (/^ *end *if|^ *end *do/i) {
        $indent--;
        if ($indent < 0) { die "Too many endif/enddo's"; }
    }
    elsif (/^ *do *\w+ *= *[-+*\/\w() ]+ *, *[-+*\/\w() ]+/i) {
        $indent++;
        $outy = 1;
    }
    elsif (/^ *do *while *\(/i) {
        $indent++;
        $outy = 1;
    }
    elsif (/^ *do *([0-9]+) *\w+ *= *[-+*\/\w() ]+ *, *[-+*\/\w() ]+/i) {
        $type{$s_pref.$1} = "enddo";
        $count{$s_pref.$1}++;
        $outy = 1;
        if ($count{$s_pref.$1} == 1) {
            $indent++;
        }
    }
    elsif (/^ *continue/i) {
        if ($before =~ /^ *([0-9]+)/) {
            if ($type{$s_pref.$1} eq "enddo") {
                $indent--;
                if ($indent < 0) {
                    die "Too many endif/enddo's";
                }
            }
        }
    }
#
# Put string back together
#
    $_ = $before.' ' x ($delta1*($indent-$outy)).$after;
#
# Check for falling off the edge of the world
#
    @lines = split(/^/);
    foreach $line (@lines) {
        if (length($line) > $line_len+1 &&    # +1 for '\n'
            !(substr($_,0,$line_len+1) =~ /!/)) {      # don't complain about long comments
            print STDERR "Line ", $.-1, " too long. Fix it\n";
        }
        print $line;
    }
}
#
# end of main program
#
# get a line, combining continuation lines
#
sub get_line {
    $thisline = $lookahead;
    if ($lookahead) {
        line: while ($lookahead = <>) {
#               if ($lookahead =~ /^     \S|^$|^[\*c#!]/i) {
            if ($lookahead =~ /^     \S/) {
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
    my $left = shift;
    my $parexp = '';

    while (/[()]/) {
        $parexp .= $`;
        $parexp .= $&;
        $_ = $';
        if ($& eq "(") { $left++; }
        else           { $left--; }
        if ($left == 0) { last; }
    }
    if ($left != 0) { die "Illegal if statement"; }
    $parexp;
}
__END__

findent - Fortran 77 indenter

=head1 SYNOPSIS

findent [--help] [--man] [-c] [-l length] [-n num] [file ...] > outfile

=head1 DESCRIPTION

Indent a Fortran 77 program. It also strips off trailing blanks.

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<-help>

Print more details about the arguments.

=item I<-man>

Print a full man page.

=item I<-c>

Formats comments with 2 spaces after the !.

=item I<-n num>

Sets the indent increment to num (default is 2).

=item I<-l length>

Specifies the maximum length of your Fortran lines (default is 72).

Version 2 - It can now deal with (most) multi-line 'if' statements.

=back

=head1 EXAMPLE

    % findent -n 3 file.f > file.indent

=head1 AUTHOR

    Kate Hedstrom
    kate@arsc.edu

With some tips from Sverre Froyen's relabel program.

do while patch by himanshu (hoberoi@limerick.cbs.umn.edu)

=head1 BUGS

    It assumes you have already run your code through the
    redo program or at least don't have any 'do' loops like:

       do 210 i=1,10
  210      mm(i) = i

    It checks to see if lines are too long, but complains
    to STDERR instead of fixing them.

    It doesn't like multiline 'if' tests with comments or blank
    lines interspersed.

=cut
