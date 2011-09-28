
SYS=P46 G46
sh=$(wildcard *.sh )

default: a
a: a.pdf

clean: c
c:
	rm -rf *.pdf *.ps

a.pdf: $(sh)
	ef.sh $(SYS) 
	scp a.pdf $(ser):~/data/plots/
