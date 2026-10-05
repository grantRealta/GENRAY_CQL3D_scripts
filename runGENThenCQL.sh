  # expects the target directory to be passed in in which to find the .pbs file 
  genray_id=$(sbatch --parsable $1/genr.slurm)
  cql_id=$(sbatch --parsable --dependency=afterok:$genray_id $1/cql.slurm)