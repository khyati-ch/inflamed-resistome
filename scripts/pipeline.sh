#!/usr/bin/env bash
# usage: bash pipeline.sh manifest.tsv
# manifest (tab-separated, with header): run_accession  fastq_ftp  fastq_md5  diagnosis

BASE=/mnt/d/ibd_project
DATA=$BASE/data
REF=$HOME/ref/GRCh38_noalt_as
CARD=$BASE/card
SCRATCH=/mnt/d/scratch
CONDA_SH=$HOME/miniconda3/etc/profile.d/conda.sh
THREADS=2
SKIP_HUMAN=1   # 1 = skip bowtie2 human-read removal, 0 = run it
CAP=5000000        # max read pairs per sample (first N pairs, via fastp)
MIN_FREE_GB=20     # stop if D: has less free space than this
MIN_FREE_HOME=10   # stop if Ubuntu home disk has less free space than this

MANIFEST=$1
mkdir -p $DATA/raw $DATA/qc $DATA/logs $DATA/amr $SCRATCH
SUMMARY=$DATA/pipeline_summary.tsv
[ -f $SUMMARY ] || printf "run\tdiagnosis\traw_pairs\ttrimmed_pairs\tnonhuman_pairs\thuman_pct\trgi_mapped_reads\tminutes\n" > $SUMMARY
source $CONDA_SH

pairs() { seqkit stats -T "$1" | awk 'NR==2{print $4}'; }
fail()  { echo "FAILED [$RUN] at: $1" | tee -a $DATA/logs/failed.txt; rm -f $SCRATCH/$RUN*; return 1; }

process() {
  RUN=$1; FTP=$2; MD5=$3; DIAG=$4
  T0=$(date +%s)
  [ -s $DATA/amr/$RUN.gene_mapping_data.txt ] && { echo "[$RUN] already done, skipping"; return 0; }

  freeD=$(df --output=avail -BG /mnt/d | tail -1 | tr -dc '0-9')
  freeH=$(df --output=avail -BG $HOME | tail -1 | tr -dc '0-9')
  [ "$freeD" -ge "$MIN_FREE_GB" ] && [ "$freeH" -ge "$MIN_FREE_HOME" ] || { echo "Low disk (D: ${freeD} GB, home ${freeH} GB). Stopping."; exit 2; }

  conda activate ibd

  # 1. download and verify md5
  IFS=';' read -ra URLS <<< "$FTP"
  IFS=';' read -ra SUMS <<< "$MD5"
  R1=""; R2=""
  for i in "${!URLS[@]}"; do
    f=$(basename "${URLS[$i]}")
    case "$f" in
      *_1.fastq.gz) R1=$DATA/raw/$f ;;
      *_2.fastq.gz) R2=$DATA/raw/$f ;;
      *) continue ;;
    esac
    ok=0
    for try in 1 2 3 4 5 6; do
      aria2c -c -x8 -s8 --file-allocation=none --console-log-level=error --summary-interval=0 --max-tries=10 --retry-wait=10 -d $DATA/raw "https://${URLS[$i]}"
      if [ -f $DATA/raw/$f ] && [ "$(md5sum $DATA/raw/$f | cut -d' ' -f1)" = "${SUMS[$i]}" ]; then ok=1; break; fi
      echo "[$RUN] md5 mismatch for $f (try $try), redownloading"; [ -f $DATA/raw/$f.aria2 ] || rm -f $DATA/raw/$f; sleep 30
    done
    [ $ok = 1 ] || { fail "download $f"; return 1; }
  done
  [ -s "$R1" ] && [ -s "$R2" ] || { fail "paired files missing"; return 1; }

  N_RAW=$(pairs $R1)
  [ "$N_RAW" = "$(pairs $R2)" ] || { fail "raw R1/R2 counts differ"; return 1; }

  # 2. trim (first CAP pairs only)
  fastp -i $R1 -I $R2 -o $SCRATCH/${RUN}_1.trim.fastq.gz -O $SCRATCH/${RUN}_2.trim.fastq.gz \
    --detect_adapter_for_pe --cut_tail --cut_mean_quality 20 --length_required 50 \
    --reads_to_process $CAP --thread $THREADS \
    -h $DATA/qc/${RUN}_fastp.html -j $DATA/qc/${RUN}_fastp.json > $DATA/logs/${RUN}_fastp.log 2>&1 || { fail fastp; return 1; }
  N_TRIM=$(pairs $SCRATCH/${RUN}_1.trim.fastq.gz)
  [ "$N_TRIM" = "$(pairs $SCRATCH/${RUN}_2.trim.fastq.gz)" ] || { fail "trimmed counts differ"; return 1; }
  rm -f $R1 $R2

  # 3. human-read removal (SKIP_HUMAN=1 skips it)
  if [ "$SKIP_HUMAN" = "1" ]; then
    NH1=$SCRATCH/${RUN}_1.trim.fastq.gz; NH2=$SCRATCH/${RUN}_2.trim.fastq.gz
    N_NH=$N_TRIM; HUMAN_PCT=NA
  else
  # 3. remove human reads
  bowtie2 -p $THREADS --sensitive -x $REF \
    -1 $SCRATCH/${RUN}_1.trim.fastq.gz -2 $SCRATCH/${RUN}_2.trim.fastq.gz \
    --un-conc-gz $SCRATCH/${RUN}_nh_%.fastq.gz -S /dev/null 2> $DATA/logs/${RUN}_bowtie2.log || { fail bowtie2; return 1; }
  rm -f $SCRATCH/${RUN}_1.trim.fastq.gz $SCRATCH/${RUN}_2.trim.fastq.gz
  NH1=$SCRATCH/${RUN}_nh_1.fastq.gz; NH2=$SCRATCH/${RUN}_nh_2.fastq.gz
  N_NH=$(pairs $NH1)
  [ "$N_NH" = "$(pairs $NH2)" ] || { fail "non-human counts differ"; return 1; }
  HUMAN_PCT=$(grep "overall alignment rate" $DATA/logs/${RUN}_bowtie2.log | cut -d'%' -f1)

  fi

  # 4. resistance genes (must run from the CARD folder because of localDB)
  conda activate rgi
  ( cd $CARD && rgi bwt -1 $NH1 -2 $NH2 -a kma -n $THREADS -o $SCRATCH/$RUN --local --clean ) > $DATA/logs/${RUN}_rgi.log 2>&1
  [ -s $SCRATCH/$RUN.gene_mapping_data.txt ] || { fail rgi; return 1; }
  for s in allele_mapping_data.txt overall_mapping_stats.txt reference_mapping_stats.txt artifacts_mapping_stats.txt gene_mapping_data.txt; do
    cp $SCRATCH/$RUN.$s $DATA/amr/
  done
  MAPPED=$(awk '/^Mapped reads:/{print $3}' $DATA/amr/$RUN.overall_mapping_stats.txt)

  mkdir -p $DATA/trimmed_keep; [ "$SKIP_HUMAN" = "1" ] && cp $NH1 $NH2 $DATA/trimmed_keep/
  rm -f $SCRATCH/$RUN*
  MIN=$(( ($(date +%s) - T0) / 60 ))
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $RUN $DIAG $N_RAW $N_TRIM $N_NH "$HUMAN_PCT" "$MAPPED" $MIN >> $SUMMARY
  return 0
}

FAILS=0
while IFS=$'\t' read -r run ftp md5 diag; do
  [ -n "$run" ] || continue
  if process "$run" "$ftp" "$md5" "$diag" < /dev/null >> $DATA/logs/pipeline_$run.log 2>&1; then
    echo "$(date +%H:%M)  $run  OK"; FAILS=0
  else
    echo "$(date +%H:%M)  $run  FAILED (see data/logs/failed.txt)"; FAILS=$((FAILS+1)); [ $FAILS -ge 3 ] && { echo "3 failures in a row, stopping"; exit 3; }
  fi
done < <(tail -n +2 "$MANIFEST" | tr -d '\r')
echo "All samples in $MANIFEST processed."
