#!/usr/bin/perl 
# header {{{
#===============================================================================
#
#         FILE:  mkt.pl
#
#        USAGE:  ./mkt.pl  
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
#      CREATED:  30/09/11 21:41:28
#     REVISION:  ---
#===============================================================================
#}}}
# use ... dhelp(); this_script, shd, bdir {{{
use strict;
use warnings;
use File::Basename;
use File::Path qw(mkpath);

my $this_script=&basename($0);
my $shd=&dirname($0);
my $bdir="$shd/../../";
my $incdir="$bdir/inc";

if ( !@ARGV ){ 
	&dhelp();
	exit 1;
}
#}}}
# vars {{{ 
my(%fold);

# folds {{{
%fold=( 
	o => "{{{",
	c => "}}}",
);

my $fo=$fold{o};
my $fc=$fold{c};
#}}}

my($dxfile,$dxdir,$dxlog,$pdir);
my(@nused,$F_NU);
$dxlog="$shd/dx.log";

my @prj=qw( $ARGV[0] );

#}}}
# subs {{{
# gettime () {{{

sub gettime(){
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $time = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	return $time;
}

#}}}
# dhelp() {{{
sub dhelp(){
print << "HELP";
=========================================================
PURPOSE: 
	Generate targets for a Fortran project
USAGE: 
	$this_script FILE
		FILE is output file with dependencies 
SCRIPT LOCATION:
	$0
=========================================================

HELP
}
#}}}
# }}}

sub main(){
#{{{
foreach my $p (@prj){
# p-dependent definitions: pdir dxdir dxfile ... {{{
$pdir="$bdir/$p";
$dxdir="$bdir/dx/$p";
$dxfile="$bdir/$p/Doxyfile";
$F_NU="$incdir/nu_$p.mk";
#}}}
# open files: NU DX {{{
open(DX,">$dxfile");
open(NU,"<$F_NU");
#}}}
# process NU {{{

while(<NU>){
	chomp;
	my $line =~ s/^\s*#.*$//g;
	push(@nused,split($line));
}

# }}}
# INPUT EXCLUDE {{{
print << "eof";
Project: 
	$p
Input project directory:
	$pdir
Output doxygen directory:
	$dxdir
eof
print DX "INPUT EXCLUDE $fo\n";
print DX "INPUT $pdir\n";
foreach my $file (@nused){
	print DX "EXCLUDE $file\n";
}
print DX << "dx";
EXCLUDE_PATTERNS       = 
EXCLUDE_SYMBOLS        =
EXCLUDE_SYMLINKS       = 
dx
print DX "$fc\n";
# }}}
# Print the rest of the Doxyfile {{{
print DX << "dx";
# Doxyfile template {{{
# A-E{{{
ABBREVIATE_BRIEF       =
ALIASES                = "input=\\par Input parameters:\n" 
ALLEXTERNALS           = NO
ALPHABETICAL_INDEX     = YES
ALWAYS_DETAILED_SEC    = NO
BINARY_TOC             = NO
BRIEF_MEMBER_DESC      = YES
BUILTIN_STL_SUPPORT    = NO
CALLER_GRAPH           = YES
CALL_GRAPH             = YES
CASE_SENSE_NAMES       = YES
CHM_FILE               =
CHM_INDEX_ENCODING     =
CLASS_DIAGRAMS         = YES
CLASS_GRAPH            = YES
COLLABORATION_GRAPH    = YES
COLS_IN_ALPHA_INDEX    = 5
COMPACT_LATEX          = NO
COMPACT_RTF            = NO
CPP_CLI_SUPPORT        = NO
CREATE_SUBDIRS         = YES
DIRECTORY_GRAPH        = YES
DISABLE_INDEX          = NO
DISTRIBUTE_GROUP_DOC   = NO
DOCSET_BUNDLE_ID       = org.doxygen.Project
DOCSET_FEEDNAME        = "Doxygen generated docs"
DOCSET_PUBLISHER_ID    = org.doxygen.Publisher
DOCSET_PUBLISHER_NAME  = Publisher
DOTFILE_DIRS           =
DOT_CLEANUP            = YES
DOT_FONTNAME           = FreeSans.ttf
DOT_FONTPATH           =
DOT_FONTSIZE           = 10
DOT_GRAPH_MAX_NODES    = 50
DOT_IMAGE_FORMAT       = png
DOT_MULTI_TARGETS      = NO
DOT_NUM_THREADS        = 0
DOT_PATH               =
DOT_TRANSPARENT        = NO
DOXYFILE_ENCODING      = UTF-8
ECLIPSE_DOC_ID         = org.doxygen.Project
# }}}
# E {{{
ENABLED_SECTIONS       = YES
ENABLE_PREPROCESSING   = YES
ENUM_VALUES_PER_LINE   = 4
EXAMPLE_PATH           =
EXAMPLE_PATTERNS       =
EXAMPLE_RECURSIVE      = NO
EXPAND_AS_DEFINED      =
EXPAND_ONLY_PREDEF     = NO
EXTENSION_MAPPING      =
EXTERNAL_GROUPS        = YES
EXTRACT_ALL            = YES
EXTRACT_ANON_NSPACES   = NO
EXTRACT_LOCAL_CLASSES  = YES
EXTRACT_LOCAL_METHODS  = NO
EXTRACT_PRIVATE        = NO
EXTRACT_STATIC         = NO
EXTRA_PACKAGES         =
EXT_LINKS_IN_WINDOW    = NO
#}}}
# F {{{
FILE_PATTERNS          = *.f *.f90 *.F
FILE_VERSION_FILTER    =
FILTER_PATTERNS        =
FILTER_SOURCE_FILES    = NO
FORCE_LOCAL_INCLUDES   = NO
FORMULA_FONTSIZE       = 10
FORMULA_TRANSPARENT    = YES
FULL_PATH_NAMES        = YES
# }}}
# G {{{
GENERATE_AUTOGEN_DEF   = NO
GENERATE_BUGLIST       = YES
GENERATE_CHI           = NO
GENERATE_DEPRECATEDLIST= YES
GENERATE_DOCSET        = NO
GENERATE_ECLIPSEHELP   = NO
GENERATE_HTML          = YES
GENERATE_HTMLHELP      = NO
GENERATE_LATEX         = YES
GENERATE_LEGEND        = YES
GENERATE_MAN           = YES
GENERATE_PERLMOD       = NO
GENERATE_QHP           = NO
GENERATE_RTF           = NO
GENERATE_TAGFILE       =
GENERATE_TESTLIST      = YES
GENERATE_TODOLIST      = YES
GENERATE_TREEVIEW      = YES
GENERATE_XML           = NO
GRAPHICAL_HIERARCHY    = YES
GROUP_GRAPHS           = YES
#}}}
# H {{{
HAVE_DOT               = NO
HHC_LOCATION           =
HIDE_FRIEND_COMPOUNDS  = NO
HIDE_IN_BODY_DOCS      = NO
HIDE_SCOPE_NAMES       = NO
HIDE_UNDOC_CLASSES     = NO
HIDE_UNDOC_MEMBERS     = NO
HIDE_UNDOC_RELATIONS   = YES
HTML_ALIGN_MEMBERS     = YES
HTML_COLORSTYLE_GAMMA  = 80
HTML_COLORSTYLE_HUE    = 220
HTML_COLORSTYLE_SAT    = 100
HTML_DYNAMIC_SECTIONS  = NO
HTML_FILE_EXTENSION    = .html
HTML_FOOTER            =
HTML_HEADER            =
HTML_OUTPUT            = html
HTML_STYLESHEET        =
HTML_TIMESTAMP         = YES
# }}}
# I-L  {{{
IDL_PROPERTY_SUPPORT   = YES
IGNORE_PREFIX          =
IMAGE_PATH             =
INCLUDED_BY_GRAPH      = YES
INCLUDE_FILE_PATTERNS  =
INCLUDE_GRAPH          = YES
INCLUDE_PATH           =
INHERIT_DOCS           = YES
INLINE_INFO            = YES
INLINE_INHERITED_MEMB  = NO
INLINE_SOURCES         = YES
INPUT_ENCODING         = UTF-8
INPUT_FILTER           =
INTERNAL_DOCS          = NO
JAVADOC_AUTOBRIEF      = NO
LATEX_BATCHMODE        = NO
LATEX_CMD_NAME         = latex
LATEX_HEADER           =
LATEX_HIDE_INDICES     = NO
LATEX_OUTPUT           = latex
LATEX_SOURCE_CODE      = NO
LAYOUT_FILE            =
#}}}
# M - Q {{{
MACRO_EXPANSION        = NO
MAKEINDEX_CMD_NAME     = makeindex
MAN_EXTENSION          = .3
MAN_LINKS              = YES
MAN_OUTPUT             = man
MAX_DOT_GRAPH_DEPTH    = 0
MAX_INITIALIZER_LINES  = 30
MSCFILE_DIRS           =
MSCGEN_PATH            =
MULTILINE_CPP_IS_BRIEF = NO
OPTIMIZE_FOR_FORTRAN   = YES
OPTIMIZE_OUTPUT_FOR_C  = NO
OPTIMIZE_OUTPUT_JAVA   = NO
OPTIMIZE_OUTPUT_VHDL   = NO
OUTPUT_DIRECTORY       = $dxdir
OUTPUT_LANGUAGE        = English
PAPER_TYPE             = a4
PDF_HYPERLINKS         = YES
PERLMOD_LATEX          = NO
PERLMOD_MAKEVAR_PREFIX =
PERLMOD_PRETTY         = YES
PERL_PATH              = /usr/bin/perl
PREDEFINED             =
PROJECT_NAME           = $p
QCH_FILE               =
QHG_LOCATION           =
QHP_CUST_FILTER_ATTRS  =
QHP_CUST_FILTER_NAME   =
QHP_NAMESPACE          = org.doxygen.Project
QHP_SECT_FILTER_ATTRS  =
QHP_VIRTUAL_FOLDER     = doc
QT_AUTOBRIEF           = NO
QUIET                  = NO
 #}}}
# R - S {{{
RECURSIVE              = YES
REFERENCED_BY_RELATION = NO
REFERENCES_LINK_SOURCE = YES
REFERENCES_RELATION    = NO
REPEAT_BRIEF           = YES
RTF_EXTENSIONS_FILE    =
RTF_HYPERLINKS         = NO
RTF_OUTPUT             = rtf
RTF_STYLESHEET_FILE    =
SEARCHENGINE           = YES
SEARCH_INCLUDES        = YES
SEPARATE_MEMBER_PAGES  = NO
SERVER_BASED_SEARCH    = NO
SHORT_NAMES            = NO
SHOW_DIRECTORIES       = NO
SHOW_FILES             = YES
SHOW_INCLUDE_FILES     = YES
SHOW_NAMESPACES        = YES
SHOW_USED_FILES        = YES
SIP_SUPPORT            = NO
SKIP_FUNCTION_MACROS   = YES
SORT_BRIEF_DOCS        = NO
SORT_BY_SCOPE_NAME     = NO
SORT_GROUP_NAMES       = NO
SORT_MEMBERS_CTORS_1ST = NO
SORT_MEMBER_DOCS       = YES
SOURCE_BROWSER         = YES
STRIP_CODE_COMMENTS    = NO
STRIP_FROM_INC_PATH    =
STRIP_FROM_PATH        =
SUBGROUPING            = YES
SYMBOL_CACHE_SIZE      = 0
#}}}
# T-X{{{
TAB_SIZE               = 8
TAGFILES               =
TEMPLATE_RELATIONS     = NO
TOC_EXPAND             = NO
TREEVIEW_WIDTH         = 250
TYPEDEF_HIDES_STRUCT   = NO
UML_LOOK               = NO
USE_HTAGS              = NO
USE_INLINE_TREES       = NO
USE_PDFLATEX           = YES
VERBATIM_HEADERS       = YES
WARNINGS               = YES
WARN_FORMAT            = "\$file:\$line: \$text"
WARN_IF_DOC_ERROR      = YES
WARN_IF_UNDOCUMENTED   = YES
WARN_LOGFILE           =
WARN_NO_PARAMDOC       = NO
XML_DTD                =
XML_OUTPUT             = xml
XML_PROGRAMLISTING     = YES
XML_SCHEMA             =
# }}}
# }}}
dx
#}}}
# close files {{{
close DX;
close NU;
# }}}

exec "doxygen 2>> $dxlog";
}
#}}}
}

&main();
