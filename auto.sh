#PBS -l nodes=29:ppn=12
 source /etc/bash.bashrc
 cd ~/ColumnMatrixCollisions
 mpirun main 29 12