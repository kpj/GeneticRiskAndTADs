bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -q normal.24h \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  -pr \
  --use-conda \
  --restart-times 3 \
  --cores 100 \
  --local-cores 1 \
  --latency-wait 30 \
  "$@"
