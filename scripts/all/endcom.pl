#!/usr/bin/perl -w

use strict;
use Pod::Usage;
use Getopt::Long;

# From the Pod::Usage man page
our ($opt_help, $opt_man);
GetOptions("help", "man")
  or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;
my $comment_char = '!';
#
# main loop
#
while (<>) {
#
# Skip comments and blank lines
#
    if (/^[c#!]|^$/i) {
        print;
        next;
    }
#
# Replace tabs with spaces
#
    s/\t/        /g;
#
# Skip short lines
#
    if (length() <= 73) {     # 72 + '\n'
        print;
        next;
    }
#
# Found a long line - but skip if comment
#
    if (substr($_,0,72) =~ /$comment_char/) {
	print;
	next;
    }
#
# Stick in the comment character plus a space
#
    my $rest = length() - 72;
    my($before, $after) = unpack("a72 a$rest", $_);
    print $before . $comment_char . " " . $after;
}
__END__

endcom - Make comments at end of line in Fortran

=head1 SYNOPSIS

endcom [--help] [--man] [file ...] > outfile

=head1 DESCRIPTION

Add a ! before the comments at the end of a line (past column 72).

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<-help>

Print more details about the arguments.

=item I<-man>

Print a full man page.

=back

=head1 AUTHOR

Kate Hedstrom
kate@arsc.edu

=cut
