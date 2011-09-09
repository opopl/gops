
export pref_default=""
export logf=out/r.log
#export pref="$pref_default"

export ps2pdf_opts="-dOptimize=true -dUseFlateCompression=true \
	       -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 \
		-dSubsetFonts=true -dEmbedAllFonts=true \
		-dAutoFilterColorImages=false \
		-dAutoFilterGrayImages=false \
		-dColorImageFilter=/FlateEncode \
		-dGrayImageFilter=/FlateEncode \
		-dModoImageFilter=/FlateEncode "	

export dvips_opts="-Pamz -Pcmz"

export pdfv="gv"
export psv="gv"

dvitopdf(){

#dvips $dvips_opts -o $1.ps $1.dvi
#ps2pdf $ps2pdf_opts $1.ps $1.pdf
dvips  -o $1.ps $1.dvi
ps2pdf $1.ps $1.pdf

}


date_dmy_hms(){

date +"%D   %H:%M:%S" 

}

date_dm_hm(){

date +"%d%m-%H%M" 

}

date_hms(){

date +\%H:\%M:\%S 

}

# decompose time in days, hours, minutes, and seconds, and display this

time_dhms(){

te=$1
if [ ! -z "$2" ]; then 
	tb=$2
else
	tb=0
fi

let t=$(($te-$tb))
# number of days
let days=$(($t/86400))

# number of seconds in the last (incomplete) day
let t=$(($t%86400))
let hrs=$(($t/3600))

let t=$(($t%3600))
let mins=$(($t/60))

let t=$(($t%60))
let secs=$t

days_s=""
if [ $days -gt 0 ]; then
 	days_s="$days (days)"
fi	

hrs_s=""
if [ $hrs -gt 0 ]; then
 	hrs_s="$hrs (hrs)"
fi	

mins_s=""
if [ $mins -gt 0 ]; then
 	mins_s="$mins (mins)"
fi	

secs_s=""
if [ $secs -gt 0 ]; then
 	secs_s="$secs (secs)"
fi	

echo "$days_s $hrs_s $mins_s $secs_s"

}

# return date in seconds

date_in_secs(){

secs=` date +"%S" `
mins=` date +"%M" `
hrs=` date +"%H" `
days=` date +"%d" `

echo " $( echo " $secs+60*$mins+3600*$hrs+86400*$days " | bc ) "

}

eo(){

  case "$1" in
   	"-n") echo "" >& $logf ;;
        *) echo "$pref> $1" >> $logf    ;;
  esac
}

esl(){

echo "**************************************************"

}

eo_s(){

echo "@@ $1" >> $outf

}

ps_merge(){

gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$1 $2 $3 -c quit

}

pdf_merge(){

pdftk $2 $3 cat output $1.new
mv $1.new $1

}

pdf_add(){

if [ -f $1 ]; then
pdf_merge $1 $1 $2
else
cp $2 $1
fi

}

ps_add(){

  if [ -f $1 ]; then
	ps_merge $1 $1 $2
      else
	cp $2 $1
fi 

}

cpl(){

while [ ! -z "$1" ]; do
  	dr=` date_dm_hm `
	#scp -r $1 op226@leonov:~/gmGof/$dr/
	#scp -r $1 op226@leonov:~/gmGof/
	scp -r $1 0502790@venus.phys.uu.nl:~/gmGof/
	shift
done

}

