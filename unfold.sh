#!/bin/bash
#SBATCH --job-name=General
#SBATCH --gres=gpu:1

module load gromacs/2020.4

cd $SLURM_SUBMIT_DIR

# README !!!
# Before executing, prepare the following dirs/files inside current directory:
#  native.pdb em.mdp nvt.mdp md.mdp nvt2.mdp npt.mdp charmm36-mar2019.ff

# generate topology:
gmx pdb2gmx -f native.pdb -o start.gro -ignh -ff CHARMM36 -water tip3p

# define simulation box:
gmx editconf -f start.gro -o box.gro -c -d 4.0 -bt dodecahedron

# add solvent:
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

# add ions:
gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
# !!!!!!!!!!!!! solvent group 13 ?

# minimize energy:
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# equlibrate with NVT:
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -update gpu -bonded gpu

# unfold:
gmx grompp -f md.mdp -c nvt.gro -p topol.top -t nvt.cpt -o md.tpr
gmx mdrun -v -deffnm md -update gpu -bonded gpu

