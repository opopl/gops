
force=0.0

force_min=0.000
force_i=0.0001
force_max=0.01

# number of starting random geometries
nrg=100
nrg=10

# random number seed
seed=1

# container radius
radius=100.0
radius=10.0

# exe
exe=stats.csh

pull_start=1
pull_end=46

nsteps=30
nsteps=1000
nsteps=10000
#nsteps=9000000

gmin_conv_always=1
adj_target=0
#
# print_ef (string) 	determines what to write into ef.dat
#		Possible values are:
#			a: 	force (E-E0)/E0 [Several lowest energies]
#
print_ef="a"
print_gmin_warnings=1
print_acceptance_ratios=1 
remove_dn=1
write_t_dat=1

n_en=5
temp=0.03

target_energy=-0.2936531215
gs_energy=-0.2936531215

sys="G46"

do_force_loop=0
