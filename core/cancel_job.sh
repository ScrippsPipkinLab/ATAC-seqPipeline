for job in "$(squeue -u snagaraja | awk 'NR > 1{ print $1 }')" 
do scancel $job 
done