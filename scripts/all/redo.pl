#!/usr/bin/perl
#
# With fixes by Robert Apthorpe to make it more portable.
# February, 2010.
#

#use strict;
use File::Spec;
use Getopt::Long;
use Pod::Usage;

my $enddo;
my ($opt_help, $opt_man);
GetOptions('enddo' => \$enddo,
           'help' => \$opt_help,
           'man' => \$opt_man)
  or pod2usage("Try '$0 --help' for more information");
pod2usage(-verbose => 1) if $opt_help;
pod2usage(-verbose => 2) if $opt_man;

# Temporary file name
my $tmpdir = File::Spec->tmpdir();
my $tmp = File::Spec->catfile($tmpdir, qq{redo.$$});

my $section_no = 0;
my $s_pref = $section_no."_";
my $cont_char = '&';
my %type;
my %count;
my %goto;
my ($left);
#
# first pass - find all "enddo" labels
#
open(TMP,">$tmp") || die "Can't open tmp file";
my $lookahead = <>;                # Get first line
while ($_ = get_line('ARGV')) {
#
# Skip comments and blank lines
#
    if (/^[c#!*]|^$/i) {
        print TMP;
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
# Do some simple tests to see if line needs further processing
# (to speed things up)
#
    if ($_ !~ /^ *do|^ *if|^ *go *to/is) {
        print TMP;
        next;
    }
#
# found a do loop
#
    if (/^( *do) *(\d+) *(\w+ *= *[-+*\/\w() ]+ *, *[-+*\/\w() ]+)/is)
    {
        $type{$s_pref.$2} = "enddo";
        $count{$s_pref.$2}++;
        if ($enddo) {
            print TMP $1;
            print TMP ' ';
            print TMP $3;
            $_ = $';
        }
    }
    print TMP;
#
# look for goto's
# Must first skip past `if (...)' constructs
#
        $_ = substr($_,6,999);
        if (/^ *if *\(/is) {
                $_ = $';
                $left = 1;
                find_match();
                if ($left != 0) { die "Illegal if statement"; }
        }
#
# goto
#
        if (/^ *go *to */is) {
                $_ = $';
                if (/^([0-9]+) *$/i) {        # only simplest type
            $goto{$s_pref.$1} = "true";
                        $_ = $';
                }
                else { die "Illegal goto"; }
        }
}
close(TMP);
#
# Second pass - check for continues on do loops
#
open(TMP,"$tmp") || die "Can't open tmp file - second pass";
$lookahead = <TMP>;                # Get first line
$section_no = 0;
$s_pref = $section_no."_";
while ($_ = &get_line('TMP')) {
#
# Skip comments and blank lines
#
    if (/^[c#*]|^$/i) {
        print;
        next;
    }
#
# Check for new section (function or subroutine)
#
    if (/function|subroutine/is && $` !~ /'/) {
        $section_no++;
        $s_pref = $section_no."_";
    }
#
# Check for numeric label field
#
    my $label = substr($_,0,5);
    my $label_field;
    if ($label =~ s/^[ 0]*([1-9][0-9]*) */$1/) {
        $label_field = $s_pref.$label;
        if ($type{$label_field} eq "enddo") {
            my $rest = substr($_,6,999);
            if ($rest =~ /^ *format/i) {    # Label type
                $type{$label_field} = "format";
                close(TMP);
                unlink $tmp;
                die "wrong label type $label_field";
            }
            elsif ($rest =~ /^ *continue/i) {
                $type{$label_field} = "continue";
                if ($enddo) {
                    $_ = "";
                }
            }
            else {
                substr($_,0,5) = "     ";
            }
        }
    }
    print;
    if (($type{$label_field} eq "enddo") && !$enddo) {
        printf ("%5d continue\n", $label);    # Print label
    }
    if ($enddo && ($type{$label_field} eq "enddo" ||
                    $type{$label_field} eq "continue")) {
        for(my $i=1;$i<=$count{$label_field};$i++){
            print "      enddo\n";
        }
        if ($goto{$label_field} eq "true") {
            printf ("%5d continue\n", $label);
        }
    }
    $type{$label_field} = "";        # until next label
}
#
# cleanup
#
close(TMP);
unlink $tmp;
#
# end main, start subroutines
# get a line, combining continuation lines
#
sub get_line {
    my $file = shift;
    local($_);

    my $thisline = $lookahead;
    if ($lookahead) {
        line: while ($lookahead = <$file>) {
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
# Find matching parenthesis
#
sub find_match {
    my $parexp = '';
    while (/[()]/) {
        $parexp .= $`;
        $parexp .= $&;
        $_ = $';
        if ($& eq "(") { $left++; }
        else           { $left--; }
        if ($left == 0) { last; }
    }
}
__END__

redo - Reformat Fortran do loops

=head1 SYNOPSIS

redo [--help] [--man] [--enddo] [file ...] > outfile

=head1 DESCRIPTION

Update do loops that don't end in a continue or enddo (see example).

=head1 OPTIONS AND ARGUMENTS

=over 4

=item I<--help>

Print more details about the arguments.

=item I<--man>

Print a full man page.

=item I<--enddo>

Put enddo at the end of the loop instead of continue.

=back

=head1 EXAMPLE

Add continue statements to "do" loops in a Fortran program:

      do 10 i=1,20
  10     sum = sum+i

goes to:

      do 10 i=1,20
         sum = sum+i
  10  continue

The --enddo option causes it to produce the following instead:

      do i=1,20
         sum = sum+i
      enddo

=head1 BUGS

   The --enddo option will break your code if you jump to do loop
      labels in any way other than a simple goto.
   It doesn't indent the way I would.

=head1 AUTHOR

Kate Hedstrom
kate@arsc.edu

With many pieces from Sverre Froyen's relabel program.

=cut
