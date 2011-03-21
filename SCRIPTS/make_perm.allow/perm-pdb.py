#! /usr/bin/env python
import sys

###################################################################################
## version of perm.py programm written by May Khalili modified by Edyta Malolepsza
##
## 1. Use on PDB file
## 2. PDB should be in all-atom representation
## 3. PDB file which should contain END on the end of sequence
## 4. Permutation is carried out only for atoms belonging to proteins and nucleic
##
## use as ./perm-pdb-allatom.py NAME.PDB
##
###################################################################################

inp = open(sys.argv[1])

print '\n---------------------------------------------------------------------------'
print 'Please check if names of terminal residues for proteins contain \'N\' and \'C\''
print 'for N and C terminus, respectively. If not, please change them, e.g. ALA -> '
print 'NALA or CALA. Otherwise they will be treated as not terminal residues'
print '---------------------------------------------------------------------------\n'

out=open('perm.allow','w')

class Atom:
  name=""
  index =0
  acidname=""
  def __str__(this):
    return this.name
  def __init__(this):
    this.name = ''
  def __str__(this):
    return this.name+this.acidname+" %d "%(this.index)
  def __init__(this):
    this.name = ''
    this.index = 0
    this.acidname=''
 
  def copy(this):
    n = Atom()
    n.name=this.name
    n.acidname=this.acidname
    n.atom=this.index
    return n

# This function reads and stores atom name, its corresponding residue name and its index
def readatom(s):
  a = Atom()
  l = s.split()
  a.acidname=l[3]
  a.name = l[2]
  a.index = int(l[1])
  return a

prev=0
finals=[]
for each in inp:
  els = each.split()
  if len(els)==1 and prev==0: continue
  if (els[0]=='TER'): continue
  if prev==0: 
    prev=els[4]
    ATMlist=[]
  if ((els[0] != 'END') and (els[4]==prev) and (prev != 1)):
    AM=readatom(each)
    ATMlist.append(AM)
     
  else:
    atnum=[]
    atnum2=[]
    atnum3=[]
    atnum4=[]
    atnum5=[]
    atnum6=[]
    group=[]     # for group permutation
    count=2      # number of permutable atoms
    groupcount=0 # number of permutable atoms in group
    swap=0       # number of other pairs of atoms that must swap if the first pair is permuted
    group2=[]
    swap2=0
    groupcount2=0

    ###############
    ## amino acids
    ###############
    if(ATMlist[0].acidname=='GLN'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
      count3=2
      atnum3.append(ATMlist[13].index)
      atnum3.append(ATMlist[14].index)
    elif(ATMlist[0].acidname=='NGLN'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
      count3=2
      atnum3.append(ATMlist[15].index)
      atnum3.append(ATMlist[16].index)
      count4=3
      atnum4.append(ATMlist[1].index)
      atnum4.append(ATMlist[2].index)
      atnum4.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CGLN'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
      count3=2
      atnum3.append(ATMlist[13].index)
      atnum3.append(ATMlist[14].index)
      count4=2
      atnum4.append(ATMlist[16].index)
      atnum4.append(ATMlist[17].index)
    elif(ATMlist[0].acidname=='GLH'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
    elif(ATMlist[0].acidname=='ARG'): 
      swap=2
      groupcount=2
      group.append(ATMlist[16].index)
      group.append(ATMlist[19].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[20].index)
      group.append(ATMlist[18].index)
      group.append(ATMlist[21].index)
      atnum.append(ATMlist[17].index)
      atnum.append(ATMlist[18].index)
      count2=2
      atnum2.append(ATMlist[20].index)
      atnum2.append(ATMlist[21].index)
      count3=2
      atnum3.append(ATMlist[5].index)
      atnum3.append(ATMlist[6].index)
      count4=2
      atnum4.append(ATMlist[8].index)
      atnum4.append(ATMlist[9].index)
      count5=2
      atnum5.append(ATMlist[11].index)
      atnum5.append(ATMlist[12].index)
    elif(ATMlist[0].acidname=='NARG'): 
      swap=2
      groupcount=2
      group.append(ATMlist[18].index)
      group.append(ATMlist[21].index)
      group.append(ATMlist[19].index)
      group.append(ATMlist[22].index)
      group.append(ATMlist[20].index)
      group.append(ATMlist[23].index)
      atnum.append(ATMlist[19].index)
      atnum.append(ATMlist[20].index)
      count2=2
      atnum2.append(ATMlist[22].index)
      atnum2.append(ATMlist[23].index)
      count3=2
      atnum3.append(ATMlist[7].index)
      atnum3.append(ATMlist[8].index)
      count4=2
      atnum4.append(ATMlist[10].index)
      atnum4.append(ATMlist[11].index)
      count5=2
      atnum5.append(ATMlist[13].index)
      atnum5.append(ATMlist[14].index)
      count6=3
      atnum6.append(ATMlist[1].index)
      atnum6.append(ATMlist[2].index)
      atnum6.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CARG'): 
      swap=2
      groupcount=2
      group.append(ATMlist[16].index)
      group.append(ATMlist[19].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[20].index)
      group.append(ATMlist[18].index)
      group.append(ATMlist[21].index)
      atnum.append(ATMlist[17].index)
      atnum.append(ATMlist[18].index)
      count2=2
      atnum2.append(ATMlist[20].index)
      atnum2.append(ATMlist[21].index)
      count3=2
      atnum3.append(ATMlist[5].index)
      atnum3.append(ATMlist[6].index)
      count4=2
      atnum4.append(ATMlist[8].index)
      atnum4.append(ATMlist[9].index)
      count5=2
      atnum5.append(ATMlist[11].index)
      atnum5.append(ATMlist[12].index)
      count6=2
      atnum6.append(ATMlist[23].index)
      atnum6.append(ATMlist[24].index)
    elif(ATMlist[0].acidname=='VAL'):
      swap=3
      groupcount=2
      group.append(ATMlist[6].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[7].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[8].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[11].index)
      count=3
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      atnum.append(ATMlist[9].index)
      count2=3
      atnum2.append(ATMlist[11].index)
      atnum2.append(ATMlist[12].index)
      atnum2.append(ATMlist[13].index)
    elif(ATMlist[0].acidname=='NVAL'):
      swap=3
      groupcount=2
      group.append(ATMlist[8].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[15].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[14].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[13].index)
      count=3
      atnum.append(ATMlist[9].index)
      atnum.append(ATMlist[10].index)
      atnum.append(ATMlist[11].index)
      count2=3
      atnum2.append(ATMlist[13].index)
      atnum2.append(ATMlist[14].index)
      atnum2.append(ATMlist[15].index)
      count3=3
      atnum3.append(ATMlist[1].index)
      atnum3.append(ATMlist[2].index)
      atnum3.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CVAL'):
      swap=3
      groupcount=2
      group.append(ATMlist[6].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[7].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[8].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[11].index)
      count=3
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      atnum.append(ATMlist[9].index)
      count2=3
      atnum2.append(ATMlist[11].index)
      atnum2.append(ATMlist[12].index)
      atnum2.append(ATMlist[13].index)
      count3=2
      atnum3.append(ATMlist[15].index)
      atnum3.append(ATMlist[16].index)
    elif(ATMlist[0].acidname=='ASP'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
    elif(ATMlist[0].acidname=='NASP'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
      count3=3
      atnum3.append(ATMlist[1].index)
      atnum3.append(ATMlist[2].index)
      atnum3.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CASP'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
      count3=2
      atnum3.append(ATMlist[11].index)
      atnum3.append(ATMlist[12].index)
    elif(ATMlist[0].acidname=='ASH'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='LEU'):
      swap=3
      groupcount=2
      group.append(ATMlist[9].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[16].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[15].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[14].index)
      count=3
      atnum.append(ATMlist[10].index)
      atnum.append(ATMlist[11].index)
      atnum.append(ATMlist[12].index)
      count2=3
      atnum2.append(ATMlist[14].index)
      atnum2.append(ATMlist[15].index)
      atnum2.append(ATMlist[16].index)
      count3=2
      atnum3.append(ATMlist[5].index)
      atnum3.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='NLEU'):
      swap=3
      groupcount=2
      group.append(ATMlist[11].index)
      group.append(ATMlist[15].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[18].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[14].index)
      group.append(ATMlist[16].index)
      count=3
      atnum.append(ATMlist[12].index)
      atnum.append(ATMlist[13].index)
      atnum.append(ATMlist[14].index)
      count2=3
      atnum2.append(ATMlist[16].index)
      atnum2.append(ATMlist[17].index)
      atnum2.append(ATMlist[18].index)
      count3=2
      atnum3.append(ATMlist[7].index)
      atnum3.append(ATMlist[8].index)
      count4=3
      atnum4.append(ATMlist[1].index)
      atnum4.append(ATMlist[2].index)
      atnum4.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CLEU'):
      swap=3
      groupcount=2
      group.append(ATMlist[9].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[16].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[15].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[14].index)
      count=3
      atnum.append(ATMlist[10].index)
      atnum.append(ATMlist[11].index)
      atnum.append(ATMlist[12].index)
      count2=3
      atnum2.append(ATMlist[14].index)
      atnum2.append(ATMlist[15].index)
      atnum2.append(ATMlist[16].index)
      count3=2
      atnum3.append(ATMlist[5].index)
      atnum3.append(ATMlist[6].index)
      count4=2
      atnum4.append(ATMlist[18].index)
      atnum4.append(ATMlist[19].index)
    elif(ATMlist[0].acidname=='GLU'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
      count3=2
      atnum3.append(ATMlist[11].index)
      atnum3.append(ATMlist[12].index)
    elif(ATMlist[0].acidname=='NGLU'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
      count3=2
      atnum3.append(ATMlist[13].index)
      atnum3.append(ATMlist[14].index)
      count4=3
      atnum4.append(ATMlist[1].index)
      atnum4.append(ATMlist[2].index)
      atnum4.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CGLU'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
      count3=2
      atnum3.append(ATMlist[11].index)
      atnum3.append(ATMlist[12].index)
      count4=2
      atnum4.append(ATMlist[14].index)
      atnum4.append(ATMlist[15].index)
    elif(ATMlist[0].acidname=='ASN'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
    elif(ATMlist[0].acidname=='NASN'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=2
      atnum2.append(ATMlist[12].index)
      atnum2.append(ATMlist[13].index)
      count3=3
      atnum3.append(ATMlist[1].index)
      atnum3.append(ATMlist[2].index)
      atnum3.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CASN'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
      count3=2
      atnum3.append(ATMlist[13].index)
      atnum3.append(ATMlist[14].index)
    elif(ATMlist[0].acidname=='LYS'):
      count=3
      atnum.append(ATMlist[17].index)
      atnum.append(ATMlist[18].index)
      atnum.append(ATMlist[19].index)			
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[8].index)
      atnum3.append(ATMlist[9].index)
      count4=2
      atnum4.append(ATMlist[11].index)
      atnum4.append(ATMlist[12].index)
      count5=2
      atnum5.append(ATMlist[14].index)
      atnum5.append(ATMlist[15].index)
    elif(ATMlist[0].acidname=='NLYS'):
      count=3
      atnum.append(ATMlist[19].index)
      atnum.append(ATMlist[20].index)
      atnum.append(ATMlist[21].index)			
      count2=2
      atnum2.append(ATMlist[7].index)
      atnum2.append(ATMlist[8].index)
      count3=2
      atnum3.append(ATMlist[10].index)
      atnum3.append(ATMlist[11].index)
      count4=2
      atnum4.append(ATMlist[13].index)
      atnum4.append(ATMlist[14].index)
      count5=2
      atnum5.append(ATMlist[16].index)
      atnum5.append(ATMlist[17].index)
      count6=3
      atnum6.append(ATMlist[1].index)
      atnum6.append(ATMlist[2].index)
      atnum6.append(ATMlist[3].index)			
    elif(ATMlist[0].acidname=='CLYS'):
      count=3
      atnum.append(ATMlist[17].index)
      atnum.append(ATMlist[18].index)
      atnum.append(ATMlist[19].index)			
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[8].index)
      atnum3.append(ATMlist[9].index)
      count4=2
      atnum4.append(ATMlist[11].index)
      atnum4.append(ATMlist[12].index)
      count5=2
      atnum5.append(ATMlist[14].index)
      atnum5.append(ATMlist[15].index)
      count6=2
      atnum6.append(ATMlist[21].index)
      atnum6.append(ATMlist[22].index)
    elif(ATMlist[0].acidname=='LYN'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)
      count3=2
      atnum3.append(ATMlist[11].index)
      atnum3.append(ATMlist[12].index)
      count4=2
      atnum4.append(ATMlist[14].index)
      atnum4.append(ATMlist[15].index)
    elif(ATMlist[0].acidname=='TYR'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      swap=3
      groupcount=2
      group.append(ATMlist[8].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[15].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[18].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[16].index)
#      count2=2                          #can i permute these two hydrogen????
#      atnum2.append(ATMlist[9].index)
#      atnum2.append(ATMlist[18].index)
#      count3=2
#      atnum3.append(ATMlist[11].index)
#      atnum3.append(ATMlist[16].index)
    elif(ATMlist[0].acidname=='NTYR'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      swap=3
      groupcount=2
      group.append(ATMlist[10].index)
      group.append(ATMlist[19].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[20].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[18].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CTYR'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      swap=3
      groupcount=2
      group.append(ATMlist[8].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[15].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[18].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[16].index)
      count2=2
      atnum2.append(ATMlist[20].index)
      atnum2.append(ATMlist[21].index)
    elif(ATMlist[0].acidname=='PHE'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      swap=3
      groupcount=2
      group.append(ATMlist[8].index)
      group.append(ATMlist[16].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[14].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[15].index)
#      count2=2
#      atnum2.append(ATMlist[9].index)
#      atnum2.append(ATMlist[17].index)
#      count3=2
#      atnum3.append(ATMlist[11].index)
#      atnum3.append(ATMlist[15].index)
    elif(ATMlist[0].acidname=='NPHE'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      swap=3
      groupcount=2
      group.append(ATMlist[10].index)
      group.append(ATMlist[18].index)
      group.append(ATMlist[12].index)
      group.append(ATMlist[16].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[19].index)
      group.append(ATMlist[13].index)
      group.append(ATMlist[17].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)			
    elif(ATMlist[0].acidname=='CPHE'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      swap=3
      groupcount=2
      group.append(ATMlist[8].index)
      group.append(ATMlist[16].index)
      group.append(ATMlist[10].index)
      group.append(ATMlist[14].index)
      group.append(ATMlist[9].index)
      group.append(ATMlist[17].index)
      group.append(ATMlist[11].index)
      group.append(ATMlist[15].index)
      count2=2
      atnum2.append(ATMlist[19].index)
      atnum2.append(ATMlist[20].index)
    elif(ATMlist[0].acidname=='MET'):
      count=3
      atnum.append(ATMlist[12].index)
      atnum.append(ATMlist[13].index)
      atnum.append(ATMlist[14].index)			
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[8].index)
      atnum3.append(ATMlist[9].index)
    elif(ATMlist[0].acidname=='NMET'):
      count=3
      atnum.append(ATMlist[14].index)
      atnum.append(ATMlist[15].index)
      atnum.append(ATMlist[16].index)			
      count2=2
      atnum2.append(ATMlist[7].index)
      atnum2.append(ATMlist[8].index)
      count3=2
      atnum3.append(ATMlist[10].index)
      atnum3.append(ATMlist[11].index)
      count4=3
      atnum4.append(ATMlist[1].index)
      atnum4.append(ATMlist[2].index)
      atnum4.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CMET'):
      count=3
      atnum.append(ATMlist[12].index)
      atnum.append(ATMlist[13].index)
      atnum.append(ATMlist[14].index)			
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[8].index)
      atnum3.append(ATMlist[9].index)
      count4=2
      atnum4.append(ATMlist[16].index)
      atnum4.append(ATMlist[17].index)
    elif(ATMlist[0].acidname=='ILE'):
      count=3
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      atnum.append(ATMlist[9].index)			
      count2=2
      atnum2.append(ATMlist[11].index)
      atnum2.append(ATMlist[12].index)
      count3=3
      atnum3.append(ATMlist[14].index)
      atnum3.append(ATMlist[15].index)
      atnum3.append(ATMlist[16].index)
    elif(ATMlist[0].acidname=='NILE'):
      count=3
      atnum.append(ATMlist[9].index)
      atnum.append(ATMlist[10].index)
      atnum.append(ATMlist[11].index)			
      count2=2
      atnum2.append(ATMlist[13].index)
      atnum2.append(ATMlist[14].index)
      count3=3
      atnum3.append(ATMlist[16].index)
      atnum3.append(ATMlist[17].index)
      atnum3.append(ATMlist[18].index)
      count4=3
      atnum4.append(ATMlist[1].index)
      atnum4.append(ATMlist[2].index)
      atnum4.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CILE'):
      count=3
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      atnum.append(ATMlist[9].index)			
      count2=2
      atnum2.append(ATMlist[11].index)
      atnum2.append(ATMlist[12].index)
      count3=3
      atnum3.append(ATMlist[14].index)
      atnum3.append(ATMlist[15].index)
      atnum3.append(ATMlist[16].index)
      count4=2
      atnum4.append(ATMlist[18].index)
      atnum4.append(ATMlist[19].index)
    elif(ATMlist[0].acidname=='SER'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='NSER'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CSER'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
    elif(ATMlist[0].acidname=='GLY'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
    elif(ATMlist[0].acidname=='NGLY'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CGLY'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=2
      atnum2.append(ATMlist[6].index)
      atnum2.append(ATMlist[7].index)
    elif(ATMlist[0].acidname=='TRP'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='NTRP'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CTRP'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[23].index)
      atnum2.append(ATMlist[24].index)
    elif(ATMlist[0].acidname=='HIS' or ATMlist[0].acidname=='HIE' or ATMlist[0].acidname=='HID' or ATMlist[0].acidname=='HIP'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='NHIS' or ATMlist[0].acidname=='NHIE' or ATMlist[0].acidname=='NHID' or ATMlist[0].acidname=='NHIP'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)			
    elif(ATMlist[0].acidname=='CHIS' or ATMlist[0].acidname=='CHIE' or ATMlist[0].acidname=='CHID' or ATMlist[0].acidname=='CHIP'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[16].index)
      atnum2.append(ATMlist[17].index)
    elif(ATMlist[0].acidname=='ALA'):
      count=3
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      atnum.append(ATMlist[7].index)			
    elif(ATMlist[0].acidname=='NALA'):
      count=3
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      atnum.append(ATMlist[3].index)			
      count2=3
      atnum2.append(ATMlist[7].index)
      atnum2.append(ATMlist[8].index)
      atnum2.append(ATMlist[9].index)			
    elif(ATMlist[0].acidname=='CALA'):
      atnum.append(ATMlist[9].index)
      atnum.append(ATMlist[10].index)			
      count2=3
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      atnum2.append(ATMlist[7].index)			
    elif(ATMlist[0].acidname=='THR'):
      count=3
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      atnum.append(ATMlist[9].index)
    elif(ATMlist[0].acidname=='NTHR'):
      count=3
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      atnum.append(ATMlist[3].index)
      count2=3
      atnum2.append(ATMlist[9].index)
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
    elif(ATMlist[0].acidname=='CTHR'):
      count=3
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      atnum.append(ATMlist[9].index)
      count2=2
      atnum2.append(ATMlist[13].index)
      atnum2.append(ATMlist[14].index)
    elif(ATMlist[0].acidname=='CYS'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='NCYS'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CCYS'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[10].index)
      atnum2.append(ATMlist[11].index)
    elif(ATMlist[0].acidname=='CYM'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='CYX'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
    elif(ATMlist[0].acidname=='NCYX'):
      atnum.append(ATMlist[7].index)
      atnum.append(ATMlist[8].index)
      count2=3
      atnum2.append(ATMlist[1].index)
      atnum2.append(ATMlist[2].index)
      atnum2.append(ATMlist[3].index)
    elif(ATMlist[0].acidname=='CCYX'):
      atnum.append(ATMlist[5].index)
      atnum.append(ATMlist[6].index)
      count2=2
      atnum2.append(ATMlist[9].index)
      atnum2.append(ATMlist[10].index)
    elif(ATMlist[0].acidname=='PRO'):
      atnum.append(ATMlist[2].index)
      atnum.append(ATMlist[3].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[8].index)
      atnum3.append(ATMlist[9].index)
    elif(ATMlist[0].acidname=='NPRO'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[4].index)
      atnum2.append(ATMlist[5].index)
      count3=2
      atnum3.append(ATMlist[7].index)
      atnum3.append(ATMlist[8].index)
      count4=2
      atnum4.append(ATMlist[10].index)
      atnum4.append(ATMlist[11].index)
    elif(ATMlist[0].acidname=='CPRO'):
      atnum.append(ATMlist[2].index)
      atnum.append(ATMlist[3].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[8].index)
      atnum3.append(ATMlist[9].index)
      count4=2
      atnum4.append(ATMlist[13].index)
      atnum4.append(ATMlist[14].index)
    elif(ATMlist[0].acidname=='ACE'):
      count=3
      atnum.append(ATMlist[0].index)
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
    elif(ATMlist[0].acidname=='NME'):
      count=3
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      atnum.append(ATMlist[5].index)

    #####################
    ## nucleic residues
    #####################
    elif(ATMlist[0].acidname=='DA' or ATMlist[0].acidname=='DA3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[19].index)
      atnum3.append(ATMlist[20].index)
      count4=2
      atnum4.append(ATMlist[29].index)
      atnum4.append(ATMlist[30].index)
    elif(ATMlist[0].acidname=='DA5' or ATMlist[0].acidname=='DAN'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=2
      atnum2.append(ATMlist[17].index)
      atnum2.append(ATMlist[18].index)
      count3=2
      atnum3.append(ATMlist[27].index)
      atnum3.append(ATMlist[28].index)
    elif(ATMlist[0].acidname=='DC' or ATMlist[0].acidname=='DC3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[19].index)
      atnum3.append(ATMlist[20].index)
      count4=2
      atnum4.append(ATMlist[27].index)
      atnum4.append(ATMlist[28].index)
    elif(ATMlist[0].acidname=='DC5' or ATMlist[0].acidname=='DCN'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=2
      atnum2.append(ATMlist[17].index)
      atnum2.append(ATMlist[18].index)
      count3=2
      atnum3.append(ATMlist[25].index)
      atnum3.append(ATMlist[26].index)
    elif(ATMlist[0].acidname=='DG' or ATMlist[0].acidname=='DG3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[23].index)
      atnum3.append(ATMlist[24].index)
      count4=2
      atnum4.append(ATMlist[30].index)
      atnum4.append(ATMlist[31].index)
    elif(ATMlist[0].acidname=='DG5' or ATMlist[0].acidname=='DGN'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=2
      atnum2.append(ATMlist[21].index)
      atnum2.append(ATMlist[22].index)
      count3=2
      atnum3.append(ATMlist[28].index)
      atnum3.append(ATMlist[29].index)
    elif(ATMlist[0].acidname=='DT' or ATMlist[0].acidname=='DT3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=3
      atnum3.append(ATMlist[17].index)
      atnum3.append(ATMlist[18].index)
      atnum3.append(ATMlist[19].index)
      count4=2
      atnum4.append(ATMlist[29].index)
      atnum4.append(ATMlist[30].index)
    elif(ATMlist[0].acidname=='DT5' or ATMlist[0].acidname=='DTN'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=3
      atnum2.append(ATMlist[15].index)
      atnum2.append(ATMlist[16].index)
      atnum2.append(ATMlist[17].index)
      count4=2
      atnum4.append(ATMlist[27].index)
      atnum4.append(ATMlist[28].index)
    elif(ATMlist[0].acidname=='RA' or ATMlist[0].acidname=='RA3' or ATMlist[0].acidname=='RC' or ATMlist[0].acidname=='RC3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[19].index)
      atnum3.append(ATMlist[20].index)
    elif(ATMlist[0].acidname=='RA5' or ATMlist[0].acidname=='RAN' or ATMlist[0].acidname=='RC5' or ATMlist[0].acidname=='RCN'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=2
      atnum2.append(ATMlist[17].index)
      atnum2.append(ATMlist[18].index)
    elif(ATMlist[0].acidname=='RG' or ATMlist[0].acidname=='RG3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
      count3=2
      atnum3.append(ATMlist[23].index)
      atnum3.append(ATMlist[24].index)
    elif(ATMlist[0].acidname=='RG5' or ATMlist[0].acidname=='RGN'):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
      count2=2
      atnum2.append(ATMlist[21].index)
      atnum2.append(ATMlist[22].index)
    elif(ATMlist[0].acidname=='RU' or ATMlist[0].acidname=='RU3'):
      atnum.append(ATMlist[1].index)
      atnum.append(ATMlist[2].index)
      count2=2
      atnum2.append(ATMlist[5].index)
      atnum2.append(ATMlist[6].index)
    elif((ATMlist[0].acidname=='RU5') or (ATMlist[0].acidname=='RUN')):
      atnum.append(ATMlist[3].index)
      atnum.append(ATMlist[4].index)
    else:
      print 'Neither amino acid nor nucleic residue - please check residue %s' % ATMlist[0].acidname


    if(els[0]!='END'):
      ATMlist=[]
      prev=els[4]
      AM=readatom(each)
      ATMlist.append(AM)

    if(len(group)!=0):
      s=str(groupcount)+' '+str(swap)
      finals.append(s)
      s=''
      for i in range(0,len(group)):
        s=s+' '+str(group[i])
      finals.append(s)

    if(len(group2)!=0):
      s=str(groupcount2)+' '+str(swap2)
      finals.append(s)
      s=''
      for i in range(0,len(group2)):
        s=s+' '+str(group2[i])
      finals.append(s)

    if(len(atnum)!=0):
      s=str(count)+' 0'
      finals.append(s)
      s=''
      for i in range(0,len(atnum)):
        s=s+' '+str(atnum[i])
      finals.append(s)

    if(len(atnum2)!=0):
      s=str(count2)+' 0'
      finals.append(s)
      s=''
      for i in range(0,len(atnum2)):
        s=s+' '+str(atnum2[i])
      finals.append(s)	

    if(len(atnum3)!=0):
      s=str(count3)+' 0'
      finals.append(s)
      s=''
      for i in range(0,len(atnum3)):
        s=s+' '+str(atnum3[i])
      finals.append(s)	

    if(len(atnum4)!=0):
      s=str(count4)+' 0'
      finals.append(s)
      s=''
      for i in range(0,len(atnum4)):
        s=s+' '+str(atnum4[i])
      finals.append(s)	

    if(len(atnum5)!=0):
      s=str(count5)+' 0'
      finals.append(s)
      s=''
      for i in range(0,len(atnum5)):
        s=s+' '+str(atnum5[i])
      finals.append(s)	

    if(len(atnum6)!=0):
      s=str(count6)+' 0'
      finals.append(s)
      s=''
      for i in range(0,len(atnum6)):
        s=s+' '+str(atnum6[i])
      finals.append(s)	

totalperm=len(finals)/2
print >> out, totalperm
for i in finals:
  print >> out, i
    
out.close()
