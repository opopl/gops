
cat >& $r/t.jmol << EOF

load trajectory "$xyzfile"
#frame $i_frame

frame 1

width = 640;
height = 480;

set axes on; color axes white
axes scale 2.0

axes unitcell

rotate x $x_angle
rotate y $y_angle
rotate z $z_angle

zoom on
zoom $zoom

background green;
select all; spacefill 200

#set perspectiveModel 10;

labels on
color labels black
color bonds blue

set fontSize 14


# BLN sequence is B9-N3-(LB)4-N3-B9-N3-(LB)5-L  {{{

var Bc="red"
var Nc="blue"
var Lc="yellow"

select { atomno>0 and atomno <10 } ; color @Bc; label "%i"
select { atomno>9 and atomno <13 } ; color @Nc; label "%i"

select atomno=13; color @Lc; label "%i"
select atomno=14; color @Bc; label "%i"
select atomno=15; color @Lc; label "%i"
select atomno=16; color @Bc; label "%i"
select atomno=17; color @Lc; label "%i"
select atomno=18; color @Bc; label "%i"
select atomno=19; color @Lc; label "%i"
select atomno=20; color @Bc; label "%i"

select { atomno>20 and atomno < 24 } ; color @Nc; label "%i"
select { atomno>23 and atomno < 33 } ; color @Bc; label "%i"
select { atomno>32 and atomno < 36 } ; color @Nc; label "%i"

select atomno=36; color @Lc; label "%i"
select atomno=37; color @Bc; label "%i"
select atomno=38; color @Lc; label "%i"
select atomno=39; color @Bc; label "%i"
select atomno=40; color @Lc; label "%i"
select atomno=41; color @Bc; label "%i"
select atomno=42; color @Lc; label "%i"
select atomno=43; color @Bc; label "%i"
select atomno=44; color @Lc; label "%i"
select atomno=45; color @Bc; label "%i"

select atomno=46; color @Lc; label "%i"

for (var an=1; an <=45; an=an+1)
  		ann=an+1
	  	select atomno=an,atomno=ann; connect (selected)
end for

# }}}

set echo top left; font echo $font_size serif bolditalic ; color echo white
echo "$jtext"

#echo this is|myecho; set echo myecho center

var filename = "$pic" 
write IMAGE 800 600 $pic_ext @filename

EOF
#vim:set fdm=marker
