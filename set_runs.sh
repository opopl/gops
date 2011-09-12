
temp=0.03
te=-14.0
force=0.1
#te=-0.2936531215

opts="-q s1 -nap -te $te -radius 2.0 -nrg 1 -nsteps 100000 -temp $temp"

runs[0]="-q s1 -radius 2.0 -nrg 1 -nsteps 100000 -temp 0.03 --no_force_loop -force 0.0 -od r0 "
runs[1]="-q s1 -nap -sys P46 -te $te -radius 100.0 -nrg 1 -nsteps 1000000 -temp $temp --no_force_loop -force $force -od r1"

runs[2]="-q s1 -nap -sys G46 -te $te -radius 2.0 -nrg 1 -nsteps 10000 -temp $temp -if 0.0 0.01 0.001 -od 0"
runs[3]="-q s1 -nap -sys P46 -te $te -radius 2.0 -nrg 1 -nsteps 100000 -temp $temp -if 0.0 0.01 0.001 -od chfwt"

runs[4]="$opts -sys G46 -if 0.01 0.1 0.01 -od chfgo2"
runs[5]="$opts -sys G46 -if 0.1 1.0 0.1 -od chfgo3"





