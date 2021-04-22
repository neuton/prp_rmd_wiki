#!/bin/bash
#SBATCH --job-name=equilibration
#SBATCH --gres=gpu:1
cd $SLURM_SUBMIT_DIR

# README
# Usage: ./equilibrate.sh n
# Input: unfolded/unfolded_n.pdb em.mdp nvt2.mdp npt.mdp charmm36-mar2019.ff unfolded/unfolded_*.pdb

# load GROMACS environment:
module load gromacs/2020.4

# setup the box:
gmx editconf -f unfolded/unfolded_$1.pdb -o newbox.gro -c -d 1.4 -bt dodecahedron

# solvate:
cp '#topol.top.1#' newtopol.top
gmx solvate -cp newbox.gro -cs spc216.gro -o newsolvated.gro -p newtopol.top

# add ions:
gmx grompp -f em.mdp -c newsolvated.gro -p newtopol.top -o newions.tpr -maxwarn 1
echo SOL | gmx genion -s newions.tpr -o newsolv_ions.gro -p newtopol.top -pname NA -nname CL -neutral -conc 0.15

# minimize energy:
gmx grompp -f em.mdp -c newsolv_ions.gro -p newtopol.top -o newem.tpr
gmx mdrun -v -deffnm newem

# NVT equilibrate:
gmx grompp -f nvt2.mdp -c newem.gro -r newem.gro -p newtopol.top -o newnvt.tpr
gmx mdrun -v -deffnm newnvt -update gpu -bonded gpu

# NPT equilibrate:
gmx grompp -f npt.mdp -c newnvt.gro -r newnvt.gro -t newnvt.cpt -p newtopol.top -o newnpt.tpr
gmx mdrun -v -deffnm newnpt -update gpu -bonded gpu

# save output for rMD:
mkdir -p equilibrated/c$1
cp newnpt.gro equilibrated/c$1/npt.gro
cp newtopol.top equilibrated/c$1/topol.top