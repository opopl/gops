
#runs[0]="-q s1 -radius 10.0 -nrg 10 -nsteps 10000 -temp 0.1" 

runs[0]="-q s1 -radius 1.8 -nrg 1 -nsteps 1000 -temp 0.03 -ap --no_force_loop -force 0.0" 
runs[1]="-radius 100.0 -nrg 100" 
runs[3]="-temp 2.0" 

#

runs[2]="--no_force_loop -force 0.0 -nrg 1 -radius 2.0" 
runs[4]="-q s1 -radius 2.0 -nrg 1 -nsteps 10000 -temp 0.03 -ap -if 0.0 0.02 0.001" 

#
runs[5]="-q s1 -radius 2.0 -nrg 1 -nsteps 100000 -temp 0.03 -ap -if 0.0 0.02 0.001" 

