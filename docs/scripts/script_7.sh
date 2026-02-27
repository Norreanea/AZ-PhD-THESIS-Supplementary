#!/usr/bin/env bash
set -euo pipefail
#set -x
# ---------- SM paths ----------
BASEDIR="/mnt/d/scRNA-seq/small_RNA"
DATADIR="${BASEDIR}/Data"
REFDIR="${BASEDIR}/reference"
REF_FASTA="${REFDIR}/SM_ncRNA_filtered.fa"
BOWTIE_PREFIX="${REFDIR}/SM_ncRNA_filtered_bowtie2_index"  
#GTF="${REFDIR}/final_anno_only_ncRNA.gtf"                       
THREADS=8
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
MINLEN=14
# ----------------------------------------

# ---------- Conda env  ----------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate srna_env

# ---------- Directory layout ----------
OUT="${BASEDIR}/GenXPro_our"
RAWQC="${OUT}/qc/raw"
TRIMQC="${OUT}/qc/trimmed"
TRIM1="${OUT}/trimmed_adapter_q"       # pass1
CLEAN="${OUT}/cleaned_poly"            # pass2
EXTRACT="${OUT}/umi_extracted"         # UMI after trimming
DEDUP="${OUT}/dedup_pre_align"         # pre-align dedup
MAP="${OUT}/map_bowtie2"
COUNT="${OUT}/counts"
TPM="${OUT}/tpm"
LOGS="${OUT}/logs"
mkdir -p "$OUT" "$RAWQC" "$TRIMQC" "$TRIM1" "$CLEAN" "$EXTRACT" "$DEDUP" "$MAP" "$COUNT" "$TPM" "$LOGS"

# ---------- Helper: relaxed UMI extractor (8nt 5′ UMI, optional 4nt 3′ UMI) ----------
# Keeps reads lacking a recognizable 3′ UMI by assigning UMI3=NNNN
set -x
echo "UMI extraction (fast; parallel)"
n_in=$(ls -1 "${CLEAN}"/*.clean.fastq.gz 2>/dev/null | wc -l); [ "$n_in" -gt 0 ] || { echo "No inputs in ${CLEAN}"; exit 1; }

EXTRACT_PY="${OUT}/_umi_extract_relaxed.py"
cat > "$EXTRACT_PY" << 'PY'
import sys, gzip, os
from Bio.SeqIO.QualityIO import FastqGeneralIterator

inp, outp = sys.argv[1], sys.argv[2]
os.makedirs(os.path.dirname(outp), exist_ok=True)

def open_in(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

with open_in(inp) as fh, gzip.open(outp, "wb") as out:
    for h, s, q in FastqGeneralIterator(fh):
        if len(s) < 9:
            continue
        umi5 = s[:8]; body = s[8:]
        if len(body) >= 5:
            umi3 = body[-4:]; ins = body[:-4]; qins = q[8:-4]
        else:
            umi3 = "NNNN"; ins = body; qins = q[8:8+len(ins)]
        if not ins:
            continue
        out.write(f"@{h.rstrip()} UMI:{umi5}-{umi3}\n{ins}\n+\n{qins}\n".encode())
PY




#echo "UMI extraction (Biopython; 5' UMI required, 3' UMI optional)"
#for fq in "${CLEAN}"/*.clean.fastq.gz; do
#  base=$(basename "$fq" .clean.fastq.gz)
#  python -u "${EXTRACT_PY}" "$fq" "${EXTRACT}/${base}.umi.fastq.gz" \
#    2> "${LOGS}/${base}.umi_extract.stderr" || { echo "UMI extract failed for $base"; exit 1; }
#  ls -lh "${EXTRACT}/${base}.umi.fastq.gz" | cat
#  break
#done

for f in "${EXTRACT}"/*.umi.fastq.gz; do
  tmp="${f%.umi.fastq.gz}.umi.fixed.fastq.gz"
  gzip -cd "$f" | awk 'NR%4==1{$0="@"$0}1' | gzip > "$tmp" && mv "$tmp" "$f"
done



# ---------- Helper: FASTQ dedup by exact (sequence + UMI pair) ----------
DEDUP_PY="${OUT}/_dedup_by_seq_umi.py"
cat > "$DEDUP_PY" << 'PY'
import sys, gzip, hashlib, os

inp, outp, logp = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(os.path.dirname(outp), exist_ok=True)

def opn_r(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")
def opn_w(p): return gzip.open(p, "wb") if p.endswith(".gz") else open(p, "wb")

def umi_from_header(h):
    # header starts with '@'; UMI stored like "... UMI:NNNNNNNN-NNNN"
    i = h.rfind("UMI:")
    return "" if i == -1 else h[i+4:].split()[0]

seen=set(); kept=dup=0
with opn_r(inp) as r, opn_w(outp) as w, open(logp,"w") as lg:
    while True:
        h = r.readline()
        if not h: break
        s = r.readline().rstrip("\n")
        plus = r.readline()
        q = r.readline().rstrip("\n")

        if not h.startswith("@"):    # guard against malformed records
            continue
        umi = umi_from_header(h.rstrip())
        key = hashlib.md5((s + "|" + umi).encode()).hexdigest()

        if key in seen:
            dup += 1
            continue
        seen.add(key); kept += 1

        w.write(h.encode())
        w.write((s + "\n").encode())
        w.write(plus.encode())
        w.write((q + "\n").encode())

    lg.write(f"kept\t{kept}\nremoved_duplicates\t{dup}\n")
PY


# ---------- FastQC raw ----------
# fastqc -t $THREADS -o "$RAWQC" "${DATADIR}"/*.fastq.gz || true

# ---------- Cutadapt pass 1: adapter + quality (as in vendor report) ----------
# echo "Cutadapt pass 1 (adapter+qtrim)"
# for fq in "${DATADIR}"/*.fastq.gz; do
  # base=$(basename "$fq" .fastq.gz)
  # cutadapt -e 0.1 -O 3 -q 20 -m ${MINLEN} -n 8 \ #check!!!!!!!
           # -a "${ADAPTER}" \
           # -o "${TRIM1}/${base}.trim1.fastq.gz" "$fq" \
           # > "${TRIM1}/${base}.trim1.log"
# done

# # ---------- Cutadapt pass 2: homopolymer cleaning ----------
# echo "Cutadapt pass 2 (homopolymer cleaning)"
# for fq in "${TRIM1}"/*.trim1.fastq.gz; do
  # base=$(basename "$fq" .trim1.fastq.gz)
  # cutadapt -a 'A{10};o=10' -a 'T{10};o=10' -a 'C{10};o=10' -a 'G{10};o=10' \ #check!!!!!!!
           # -n 3 -m ${MINLEN} \
           # -o "${CLEAN}/${base}.clean.fastq.gz" "$fq" \
           # > "${CLEAN}/${base}.clean.log"
# done

# ---------- UMI extraction AFTER trimming (relaxed) ----------
echo "UMI extraction (fast; 5' UMI required, 3' UMI optional)"
# sanity: inputs present?
# n_in=$(ls -1 "${CLEAN}"/*.clean.fastq.gz 2>/dev/null | wc -l)
# if [ "$n_in" -eq 0 ]; then
  # echo "ERROR: No inputs in ${CLEAN}/*.clean.fastq.gz" >&2; exit 1
# fi

# # parallel if available, else serial
# if command -v parallel >/dev/null 2>&1; then
  # ls "${CLEAN}"/*.clean.fastq.gz \
  # | sed 's#.*/##; s/.clean.fastq.gz$//' \
  # | parallel -j ${THREADS} '
      # python -u "'"${EXTRACT_PY}"'" \
        # "'"${CLEAN}"'"/{}.clean.fastq.gz \
        # "'"${EXTRACT}"'"/{}.umi.fastq.gz \
      # 2> "'"${LOGS}"'"/{}.umi_extract.stderr
    # '
# else
  # for fq in "${CLEAN}"/*.clean.fastq.gz; do
    # base=$(basename "$fq" .clean.fastq.gz)
    # python -u "${EXTRACT_PY}" "$fq" "${EXTRACT}/${base}.umi.fastq.gz" \
      # 2> "${LOGS}/${base}.umi_extract.stderr"
  # done
# fi


# ---------- Deduplicate BEFORE mapping: exact (UMI pair + insert) ----------
echo "Deduplicate (pre-align) by (UMI pair + insert sequence)"
for fq in "${EXTRACT}"/*.umi.fastq.gz; do
  base=$(basename "$fq" .umi.fastq.gz)
  python "${DEDUP_PY}" "$fq" "${DEDUP}/${base}.dedup.fastq.gz" "${DEDUP}/${base}.dedup.stats.txt"
done

# ---------- Build Bowtie2 index if missing ----------
if [ ! -e "${BOWTIE_PREFIX}.1.bt2" ] && [ ! -e "${BOWTIE_PREFIX}.1.bt2l" ]; then
  echo "Building Bowtie2 index for ${REF_FASTA}"
  bowtie2-build "${REF_FASTA}" "${BOWTIE_PREFIX}"
fi

# ---------- Map to ncRNA with Bowtie2 --sensitive --local ----------
echo "Bowtie2 mapping (--sensitive --local) to ncRNA"
mkdir -p "${MAP}"
for fq in "${DEDUP}"/*.dedup.fastq.gz; do
  base=$(basename "$fq" .dedup.fastq.gz)
  echo $base
  bowtie2 --threads ${THREADS} --sensitive --local \
          -x "${BOWTIE_PREFIX}" -U "$fq" \
    2> "${MAP}/${base}.bowtie2.log" \
  | samtools view -b -F 4 - \
  | samtools sort -@4 -o "${MAP}/${base}.sorted.bam"
  samtools index "${MAP}/${base}.sorted.bam"
done

# ---------- GTF (autogenerate if absent) ----------
if [ ! -s "${GTF}" ]; then
  echo "No GTF found at ${GTF}. Autogenerating from FASTA headers (single exon per record)."
  GTF="${REFDIR}/autogen_ncRNA_from_fasta.gtf"
  python - "$REF_FASTA" "$GTF" << 'PYCODE'
import sys, gzip
fa, gtf = sys.argv[1], sys.argv[2]
def op(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
with op(fa) as fh, open(gtf,'w') as out:
    name=None; seq=[]
    def flush(nm, seqlen):
        if not nm: return
        out.write(f"{nm}\tgenxpro\ttranscript\t1\t{seqlen}\t.\t+\t.\ttranscript_id \"{nm}\"; gene_id \"{nm}\";\n")
    for line in fh:
        if line.startswith('>'):
            if name is not None:
                flush(name, len(''.join(seq)))
            name=line[1:].strip().split()[0]
            seq=[]
        else:
            seq.append(line.strip())
    if name is not None:
        flush(name, len(''.join(seq)))
PYCODE
fi

# ---------- htseq-count ----------
echo "htseq-count on BAMs"
mkdir -p "${COUNT}"
GTF="${REFDIR}/autogen_ncRNA_from_fasta.gtf"
for bam in "${MAP}"/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)
  echo $base
  htseq-count -f bam -r pos -s no -a 0 -t transcript -i transcript_id \
    --nonunique=fraction \
    "$bam" "$GTF" > "${COUNT}/${base}.counts.txt"
  #htseq-count \
  #  -f bam \
  #  -r pos \                # coordinate-sorted BAM
  #  -s no \                 # unstranded
  #  -a 0 \                  # no min AQual filter
  #  -t transcript \         # your GTF uses 'transcript'
  #  -i transcript_id \      # attribute you wrote
  #  "$bam" "$GTF" > "${COUNT}/${base}.counts.txt"
done



featureCounts -T ${THREADS} -s 0 -t transcript -g transcript_id -M --fraction \
  -a "${REFDIR}/autogen_ncRNA_from_fasta.gtf" \
  -o "${COUNT}/featureCounts.transcript.txt" \
  ${MAP}/*.sorted.bam

# ---------- Merge counts and compute TPM ----------
echo "Merge counts and compute TPM"
samtools faidx "${REF_FASTA}"
mkdir -p "${TPM}"
cut -f1,2 "${REF_FASTA}.fai" > "${TPM}/lengths.tsv"

python - "${COUNT}" "${TPM}/lengths.tsv" "${TPM}" << 'PYCODE'
import sys, os, glob, csv
from collections import defaultdict
count_dir, len_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
lengths = {}
with open(len_path) as f:
    for line in f:
        tid, ln = line.rstrip().split('\t')[:2]
        lengths[tid] = float(ln)
samples=[]; counts=defaultdict(dict)
for fn in sorted(glob.glob(os.path.join(count_dir, "*.counts.txt"))):
    s=os.path.basename(fn).replace(".counts.txt",""); samples.append(s)
    with open(fn) as fh:
        for row in fh:
            if row.startswith("__"): continue
            tid,c=row.rstrip().split('\t'); counts[tid][s]=float(c)
all_ids=list(counts.keys())
with open(os.path.join(out_dir,"counts_matrix.tsv"),"w",newline="") as out:
    w=csv.writer(out,delimiter='\t'); w.writerow(["transcript_id"]+samples)
    for tid in all_ids: w.writerow([tid]+[int(counts[tid].get(s,0)) for s in samples])
def tpm_for_sample(s):
    rpk={}; 
    for tid in all_ids:
        ln=lengths.get(tid,0.0); c=counts[tid].get(s,0.0)
        rpk[tid]=0.0 if ln<=0 else c/(ln/1000.0)
    denom=sum(rpk.values()) or 1.0
    return {tid:(v/denom)*1e6 for tid,v in rpk.items()}
per={s:tpm_for_sample(s) for s in samples}
with open(os.path.join(out_dir,"tpm_matrix.tsv"),"w",newline="") as out:
    w=csv.writer(out,delimiter='\t'); w.writerow(["transcript_id"]+samples)
    for tid in all_ids: w.writerow([tid]+[f"{per[s][tid]:.6f}" for s in samples])
PYCODE

# ---------- MultiQC ----------
multiqc "${OUT}" -o "${OUT}/multiqc"

echo "Done."
echo "Outputs:"
echo "  - Trimmed:       ${TRIM1}"
echo "  - Cleaned:       ${CLEAN}"
echo "  - UMI-extracted: ${EXTRACT}"
echo "  - Deduped:       ${DEDUP}"
echo "  - BAMs:          ${MAP}"
echo "  - Counts:        ${COUNT}"
echo "  - TPM:           ${TPM}"
