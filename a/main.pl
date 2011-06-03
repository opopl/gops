#!/usr/bin/perl 

use strict;
use warnings;
use boolean;

my @progs=( GMIN OPTIM PATHSAMPLE );
my $prog=GMIN;

# Initial part {{{

my $natoms;
my @screenc;

my %fname=(
	# {{{
	logfile			=>	$prog.log,
	coords			=>	coords,
	best			=>	best.dat,
	markov			=>	markov.dat,
	rmsd			=>	rmsd.dat,
	pairdist		=>	pairdists,
	data			=>	data
	# }}}
)

#kwd {{{

#my %kwdh=(
## {{{		=>,
#2D		=><++>,
#ACCEPTRATIO		=><++>,
#ACE		=><++>,
#ACKLAND		=><++>,
#ALGLUE		=><++>,
#CISTRANS		=><++>,
#AMBERENERGIES		=><++>,
#AMCHPMAX		=><++>,
#AMCHPMIN		=><++>,
#AMCHNMAX		=><++>,
#AMCHNMIN		=><++>,
#AMH		=><++>,
#NINT_AMH		=><++>,
#HARM_AMH		=><++>,
#ANGSTROM		=><++>,
#ARGON		=><++>,
#ARM		=><++>,
#ARNO		=><++>,
#AVOID		=><++>,
#AXTELL		=><++>,
#BFGS		=><++>,
#BHPT		=><++>,
#BINARY		=><++>,
#BINSTRUCTURES		=><++>,
#BLJCLUSTER		=><++>,
#BLN		=><++>,
#BSMIN		=><++>,
#BSPT		=><++>,
#BSPTDUMPFRQ		=><++>,
#BSPTQRANGE		=><++>,
#BSPTRESTART		=><++>,
#BSWL		=><++>,
#CAPSID		=><++>,
#CENTRE		=><++>,
#CENTREXY		=><++>,
#SETCENTRE		=><++>,
#QUCENTRE		=><++>,
#CG		=><++>,
#AMBERMDSTEPS		=><++>,
#AMBER9		=><++>,
#MODEL1		=><++>,
#MOVABLEATOMS		=><++>,
#LIGMOVE		=><++>,
#CHANGEACCEPT		=><++>,
#FREEZEGROUP		=><++>,
#ROTAMER		=><++>,
#GROUPROTATION		=><++>,
#CHARMM		=><++>,
#CHARMMENERGIES		=><++>,
#CHARMMTYPE		=><++>,
#CHECKMARKOV		=><++>,
#CHMD		=><++>,
#CHPMAX		=><++>,
#CHPMIN		=><++>,
#CHNMAX		=><++>,
#CHNMIN		=><++>,
#NOPHIPSI		=><++>,
#TOMEGA		=><++>,
#CHFREQ		=><++>,
#CHRIGIDTRANS		=><++>,
#CHRIGIDROT		=><++>,
#FIXEDEND		=><++>,
#OSASA		=><++>,
#ODIHE		=><++>,
#OEINT		=><++>,
#ORGYR		=><++>,
#CNORANDOM		=><++>,
#CPERMDIHE		=><++>,
#COLDFUSION		=><++>,
#COMPRESS		=><++>,
#COOPMOVE		=><++>,
#CPMD		=><++>,
#CSM		=><++>,
#CUTOFF		=><++>,
#D5H		=><++>,
#DBRENT		=><++>,
#DEBUG		=><++>,
#DECAY		=><++>,
#DFTB		=><++>,
#DGUESS		=><++>,
#DIELEC		=><++>,
#DIFFRACT		=><++>,
#DIPOLES		=><++>,
#DL_POLY		=><++>,
#DONTMOVE		=><++>,
#DONTMOVEGROUP		=><++>,
#DONTMOVEALL		=><++>,
#DONTMOVERES		=><++>,
#DUMP		=><++>,
#DUMPINT		=><++>,
#DZUGUTOV		=><++>,
#EAMAL		=><++>,
#EAMLJ		=><++>,
#EDIFF		=><++>,
#EQUILIBRATION		=><++>,
#EVSTEP		=><++>,
#EXPFAC		=><++>,
#FAKEWATER		=><++>,
#FAL		=><++>,
#FIXBOTH		=><++>,
#FIXCOM		=><++>,
#FIXD		=><++>,
#FIXSTEP		=><++>,
#FIXTEMP		=><++>,
#FNI		=><++>,
#FRAUSI		=><++>,
#FREEZE		=><++>,
#FREEZERES		=><++>,
#FREEZEALL		=><++>,
#FS		=><++>,
#G46		=><++>,
#GAMMA		=><++>,
#GAUSS		=><++>,
#GROUND		=><++>,
#GUIDE		=><++>,
#GUPTA		=><++>,
#HISTSMOOTH		=><++>,
#HISTTEMP		=><++>,
#IH		=><++>,
#INTMIN		=><++>,
#JM		=><++>,
#JUMPMOVE		=><++>,
#LB2		=><++>,
#LOCALSAMPLE		=><++>,
#LJMF		=><++>,
#MAKEOLIGO		=><++>,
#MAXBFGS		=><++>,
#MAXERISE		=><++>,
#MAXIT		=><++>,
#MGGLUE		=><++>,
#MORSE		=><++>,
#MPI		=><++>,
#MSORIG		=><++>,
#MSTRANS		=><++>,
#MULLERBROWN		=><++>,
#MULTIPLICITY		=><++>,
#MYSD		=><++>,
#NATB		=><++>,
#NEON		=><++>,
#NEWCONF		=><++>,
#NEWJUMP		=><++>,
#NEWRESTART		=><++>,
#NMAX		=><++>,
#NMIN		=><++>,
#NOCHIRALCHECKS		=><++>,
#UACHIRAL		=><++>,
#NOCISTRANS		=><++>,
#NOCISTRANSDNA		=><++>,
#NOCISTRANSRNA		=><++>,
#NOFREEZE		=><++>,
#NORESET		=><++>,
#OH		=><++>,
#OTP		=><++>,
#P46		=><++>,
#PACHECO		=><++>,
#PAH		=><++>,
#PAIRDIST		=><++>,
#PARALLEL		=><++>,
#PBGLUE		=><++>,
#PERIODIC		=><++>,
#PERMDIST		=><++>,
#PERMOPT		=><++>,
#PLUS		=><++>,
#PMAX		=><++>,
#PMIN		=><++>,
#POWER		=><++>,
#PROJI		=><++>,
#PROJIH		=><++>,
#PRTFRQ		=><++>,
#PTMC		=><++>,
#PULL		=><++>,
#Q4		=><++>,
#QCUTOFF		=><++>,
#QDTEST2		=><++>,
#QDTEST		=><++>,
#QUAD		=><++>,
#QUENCHDOS		=><++>,
#RADIUS		=><++>,
#RANSEED		=><++>,
#RCUTOFF		=><++>,
#READMSB		=><++>,
#RENORM		=><++>,
#RESIZE		=><++>,
#RESTART		=><++>,
#RESTORE		=><++>,
#RGCL2		=><++>,
#RINGROTSCALE		=><++>,
#RKMIN		=><++>,
#RMS		=><++>,
#SAVE		=><++>,
#SAVEINTE		=><++>,
#SC		=><++>,
#SEED		=><++>,
#SHELLMOVE		=><++>,
#SHIFTCUT		=><++>,
#SIDESTEP		=><++>,
#SORT		=><++>,
#CSQUEEZE		=><++>,
#STAR		=><++>,
#STEP		=><++>,
#STEEREDMIN		=><++>,
#TRACKDATA		=><++>,
#STEPOUT		=><++>,
#STEPS		=><++>,
#STICKY		=><++>,
#LJCOUL		=><++>,
#STOCK		=><++>,
#CHECKD		=><++>,
#DBP		=><++>,
#DBPTD		=><++>,
#DBLPY		=><++>,
#DMBLM		=><++>,
#LINROD		=><++>,
#LWOTP		=><++>,
#NCAP		=><++>,
#NPAH		=><++>,
#NTIP		=><++>,
#PATCHY		=><++>,
#ASAOOS		=><++>,
#STOCKAA		=><++>,
#MSSTK		=><++>,
#MSBIN		=><++>,
#MULTPAHA		=><++>,
#TDHD		=><++>,
#WATERDC		=><++>,
#WATERKZ		=><++>,
#GB		=><++>,
#GBD		=><++>,
#GBDP		=><++>,
#GEM		=><++>,
#PAHA		=><++>,
#PAHW99		=><++>,
#PYG		=><++>,
#PYGDP		=><++>,
#MSGB		=><++>,
#MSPYG		=><++>,
#MULTISITEPY		=><++>,
#GAYBERNE		=><++>,
#PARAMONOV		=><++>,
#PYOVERLAPTHRESH		=><++>,
#PYGPERIODIC		=><++>,
#LJCAPSID		=><++>,
#EXTRALJSITE		=><++>,
#EXTRALJSITEATTR		=><++>,
#LJSITECOORDS		=><++>,
#SWAPMOVES		=><++>,
#PYBINARY		=><++>,
#GAYBERNEDC		=><++>,
#STRAND		=><++>,
#SUPERSTEP		=><++>,
#SW		=><++>,
#SETCHIRAL		=><++>,
#SYMMETRISE		=><++>,
#SYMMETRISECSM		=><++>,
#TABOO		=><++>,
#TARGET		=><++>,
#TD		=><++>,
#TEMPERATURE		=><++>,
#TETHER		=><++>,
#THOMSON		=><++>,
#THRESHOLD		=><++>,
#TIP		=><++>,
#CTN		=><++>,
#TOLBRENT		=><++>,
#TOLD		=><++>,
#TOLE		=><++>,
#TOSI		=><++>,
#TSALLIS		=><++>,
#TWOPLUS		=><++>,
#UPDATES		=><++>,
#VGW		=><++>,
#VGWCPS		=><++>,
#VGWCPF		=><++>,
#VGWTOL		=><++>,
#VISITPROP		=><++>,
#TSTAR		=><++>,
#WELCH		=><++>,
#WENZEL		=><++>,
#ZETT1		=><++>,
#ZETT2		=><++>,
##}}}		=><++>,
#);

# default values 
my %kwddef=(
# {{{
	acceptratio 	=>	"0.5",
	changeaccept	=>	"50",
	dguess			=>	"0.1",
	ediff			=>	"0.02",
	maxbfgs			=>	"0.4",
	maxit			=>	"500 500",
	orgyr			=>	"",
	parallel		=>	"1",
	save			=>	"10",
	pull			=>	"1 46 0.01",
	sloppyconv 		=>	"10E-3",
	step			=>	"step astep ostep block",
	steps			=>	"1000 1.0",
	temperature		=>	"0.35",
	tightconv 		=>	"10E-10",
	permdist		=>	"",
# }}}
);

my %kwdh=(
#{{{
	"#"				=>	"",
	acceptratio 	=>	"accrat",
	changeaccept	=>	"naccept",
	dguess			=>	"",
	dump			=>	"",
	ediff			=>	"econv",
	maxbfgs			=>	"maxbfgs",
	maxit			=>	"maxit",
	orgyr			=>	"",
	parallel		=>	"npar",
	pull			=>	"atom1 atom2 force",
	sloppyconv 		=>	"bqmax",
	step			=>	"step astep ostep block",
	steps			=>	"mcsteps tfac",
	temperature		=>	"temp",
	tightconv 		=>	"cqmax",
	permdist		=>	"",
	target			=>	"target",
	<++>			=>	"<++>",
	<++>			=>	"<++>",
#}}}
);	

my @kwds=keys %kwdh;

# }}}

# }}}

# main part  {{{
open LFILE, ">$logfile" or die "$!";

# read in the "data" file
&keyword;
&mc;
&finalq;
&finalio;

close LFILE or die "can't close $: $!";;

# }}}

# subroutines {{{

sub finalio {
# {{{

open LOW,">lowest" or die "$!";

for ( $iconf = 0; $iconf < $nsave; $iconf++){
	print LOW, "$natoms\n";
	print LOW, "Energy of minimum $iconf = $qmin[$iconf]  first found at step $ff[$iconf] \n";
         WRITE(MYUNIT2,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30       FORMAT('LA ',3F20.10)
	print LOW, "LA $qminp \n";
}

close LOW or die "can't close $: $!";;

# }}}
}

sub finalq {
# {{{
# }}}
}

sub keyword {
# {{{
open FILE,"<data" or die "$!";

while(<FILE>){
	my( $k, @vars )= split $_;
	foreach $kwd (@kwds){
		( lc($k) eq $kwd ) && &handle_kwd($kwd);
		}
}

close FILE or die "can't close $: $!";;
# }}}
}

sub mc {
# {{{
my ($nsteps,$scalefac,@screenc)=@_;

print LFILE, "Calculating initial energy\n";
print LFILE, "Qu $nq E= $potel steps=$iterations RMS=$rms Markov E=$potel t=($time-$tstart) \n"; 

$nqtot=$nqtot+1;

print LFILE, "Starting MC run of",$nsteps," steps \n";
print LFILE, "Temperature will be multiplied by ",$scalefac," at every step\n";
&bhloop;
#}}}
}

#}}}
