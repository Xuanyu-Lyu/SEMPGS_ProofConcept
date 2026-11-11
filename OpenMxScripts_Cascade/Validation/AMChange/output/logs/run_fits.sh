#!/bin/zsh
# sequential fit chain for the AMChange validation (3 workers x 2 threads at a time)
cd "$(dirname "$0")/../.."
echo "START mvn continuous $(date)";  Rscript 02_fit_models.R --mode mvn --binary 0 --workers 3 > output/logs/02_mvn_continuous_b.log 2>&1
echo "START reduction $(date)";       Rscript 02_fit_models.R --mode reduction --exhaustive 0 --workers 3 > output/logs/02_reduction.log 2>&1
echo "START sim continuous $(date)";  Rscript 02_fit_models.R --mode sim --binary 0 --branch all --workers 3 > output/logs/02_sim_continuous.log 2>&1
echo "START mvn binary $(date)";      Rscript 02_fit_models.R --mode mvn --binary 1 --workers 3 > output/logs/02_mvn_binary.log 2>&1
echo "START sim binary $(date)";      Rscript 02_fit_models.R --mode sim --binary 1 --workers 3 > output/logs/02_sim_binary.log 2>&1
echo "DONE $(date)"
