#!/usr/bin/perl

use strict;
use Pod::Usage;
use Getopt::Long;

our ($opt_help, $opt_man);
GetOptions("help", "man")
  or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

my $cont_char = '&';
my $do_count = 0;
#
# set these up to not conflict with existing labels
#
my $label_no_start = 2000;
my $label_no_incr = 10;
my $label_no = $label_no_start;
my $lookahead;
my @label;
#
# get a line, combining continuation lines
#
sub get_line {
    my $thisline = $lookahead;
    line: while ($lookahead = <>) {
        if ($lookahead =~ s/^     \S/     $cont_char/) {
            $thisline .= $lookahead;
        }
        else {
            last line;
        }
    }
    $thisline;
}
#
# main loop - look for do-enddo loops
#
while (<>) {
#
# Skip comments and blank lines
#
    if (/^[c#]|^$/i) {
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
    if (/function|subroutine/i && $` !~ /'/) {
        if ($do_count != 0) {
            die "Unmatched do-enddo";
        }
        $label_no = $label_no_start;
    }
#
# Do some simple tests to see if line needs further processing
# (to speed things up)
#
    if ($_ !~ /^ *do|^ *end *do/i) {
        print;
        next;
    }
#
# found a do loop
#
    if (/^( *do *)(\w+ *= *[-+*\/\w() ]+ *, *[-+*\/\w() ]+)/i)
    {
        $do_count++;
        $label[$do_count] = $label_no;
        $label_no += $label_no_incr;
        print $1;
        printf (" %s ", $label[$do_count]);
        print $2;
        print $';
        next;
    }
#
# look for enddo's
#
        if (/^ *end *do */i) {
        printf ("%5d continue\n", $label[$do_count]);
        $do_count--;
        next;
        }
    print $_;
}
__END__

unenddo - Fortran 77 reformatting tool

=head1 SYNOPSIS

unenddo [--help] [--man] [file ...] > outfile

=head1 DESCRIPTION

Change Fortran do-enddo loops to do-continue loops.

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<-help>

Print more details about the arguments.

=item I<-man>

Print a full man page.

=back

=head1 BUGS

It doesn't check to see if new labels are already being used.

It doesn't indent the way I would (but see findent).

=head1 AUTHOR

Kate Hedstrom
kate@arsc.edu

=cut
