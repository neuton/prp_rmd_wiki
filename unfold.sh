#!/bin/bash
#SBATCH --job-name=unfolding
#SBATCH --gres=gpu:1
cd $SLURM_SUBMIT_DIR

# README
# Input: native.pdb em.mdp nvt.mdp md.mdp charmm36-mar2019.ff

# load GROMACS environment:
module load gromacs/2020.4

# generate topology:
gmx pdb2gmx -f native.pdb -o start.gro -ignh -ff charmm36-mar2019 -water tip3p

# define simulation box:
gmx editconf -f start.gro -o box.gro -c -d 4.0 -bt dodecahedron

# add solvent:
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

# add ions:
gmx grompp -f em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# minimize energy:
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# equlibrate with NVT:
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt -update gpu -bonded gpu

# unfold:
gmx grompp -f md.mdp -c nvt.gro -p topol.top -t nvt.cpt -o md.tpr
gmx mdrun -v -deffnm md -update gpu -bonded gpu

# calculate RMSD:
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg

# save frames with high RMSD:
mkdir -p unfolded
echo 1 1 | gmx trjconv -f md.xtc -s md.tpr -pbc mol -center -sep -o unfolded/unfolded_.pdb -skip 3 -drop rmsd.xvg -dropunder 2.2