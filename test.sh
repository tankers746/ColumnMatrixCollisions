for i in 2 5 8 11 14 17 20 23 26 29      ### Nodes ###
do
    for j in 1 3 6 9 12 ### PPN ###
    do
	printf "#PBS -l nodes=$i:ppn=$j\n source /etc/bash.bashrc\n cd ~/ColumnMatrixCollisions\n mpirun main $i $j" > auto.sh
	sleep 1
	qsub auto.sh
    done
done
