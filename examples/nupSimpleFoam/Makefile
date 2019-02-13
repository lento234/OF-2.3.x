# File			: Makefile
# Created		: Oct 23 2018
# Description	: Makefile for run the case
# Copyright 2018 Lento Manickathan.

# cluster name
# (hypatia | daint | euler)
#CLUSTER = hypatia
CLUSTER = hypatia

# submit flags
SJOBFLAGS	= 
ENVTARGETS 	=

# Leaf area density mapping
SOURCE = ladfieldmesh
FIELD  = a

# submit file
FILE = submit.$(CLUSTER)

# Folder name
FOLDERNAME = $(shell basename $(CURDIR))

#HOSTNAME = $(shell hostname)

# setup cluster specific variables
ifeq ($(CLUSTER),hypatia)
	# job submit script
	SJOB = qsub
	SCRATCHPATH = /mnt/scratch/male
	JOBNAMEPATTERN =PBS -N
else
	SJOB = sbatch
	ENVTARGETS += daintenv
	SCRATCHPATH = /scratch/snx3000/mlento
	JOBNAMEPATTERN =SBATCH --job-name=
endif



# Commands
.PHONY: all purge clean reconstruct mesh run daintenv push pull initialize map

#hostname:
#	echo $(HOSTNAME)

all: purge run

purge: clean
	rm -f slurm-*.out *.e[0-9]* log.*
	rm -rf constant/polyMesh/{boundary,faces,neighbour,owner,points}
	rm -rf [1-9]* postProcessing* logs VTK
	#rm -f 0/{y,yPlus,a,phi.*}
	rm -rf 0

map:
	#cp -v 0.org/$(FIELD) 0/$(FIELD)
	#mapFields $(SOURCE) -case . -fields '($(FIELD))'
	initLAD

decompose:
	decomposePar

vtk:
	foamToVTK

clean:
	rm -rf processor*

reconstruct: $(ENVTARGETS)
	reconstructPar

mesh: $(ENVTARGETS)
	blockMesh
	cp -rv 0.org 0

initialize:
	cp -v 0.org/U 0/U
	potentialFoam > log.potentialFoam

daintenv:
	source ~/loadopenfoam

kill:
	qstat -u male | grep wtcomp | grep 'R' | cut -d " " -f1 | xargs qdel

pull:
	rsync -avhu $(CLUSTER):$(SCRATCHPATH)/$(FOLDERNAME)/ . --exclude='processor*'

deletepull:
	rsync -avhu --delete $(CLUSTER):$(SCRATCHPATH)/$(FOLDERNAME)/ . --exclude='processor*'

sample:
	sample -latestTime

push:
	rsync -avhu . $(CLUSTER):$(SCRATCHPATH)/$(FOLDERNAME)

run: $(ENVTARGETS) $(FILE)
	sed -i "/#$(JOBNAMEPATTERN)/c #$(JOBNAMEPATTERN) $(FOLDERNAME)" $(FILE)
	$(SJOB) $(SJOBFLAGS) $(FILE)
