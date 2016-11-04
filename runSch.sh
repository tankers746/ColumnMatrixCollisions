#PBS -l nodes=17:ppn=3
source /etc/bash.bashrc
cd ~/ColumnMatrixCollisions
mpirun main
