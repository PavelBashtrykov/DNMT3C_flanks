#!/usr/bin/bash -i

################################################################################
# This script is written by Pavel Bashtrykov
# pavel.bashtrykov@ibtb.uni-stuttgart.de
# pavel.bashtrykov@gmail.com
################################################################################

# PROVIDE ABSOLUTE PATH TO THE GENOME AND WGBS2FLANKS SCRIPTS
PATH_GENOME_FASTA=~/HDD/genomes/mm9/mm9.fa
PATH_WGBS2FLANKS=~/tools/wgbs_draft/

# GET DATA FROM GEO
SAMPLES=(GSM3239884_SKO_mES_BS_seq_mCG.txt GSM4809269_SKO_Dnmt1_KO_BS_seq_mCG.txt)

wget -c https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3239nnn/GSM3239884/suppl/GSM3239884_SKO_mES_BS_seq_mCG.txt.gz
gzip -dk GSM3239884_SKO_mES_BS_seq_mCG.txt.gz

wget -c https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4809nnn/GSM4809269/suppl/GSM4809269_SKO_Dnmt1_KO_BS_seq_mCG.txt.gz
gzip -dk GSM4809269_SKO_Dnmt1_KO_BS_seq_mCG.txt.gz


# FLANKING SEQUENCE ANALYISIS

analyse_flanks(){
    python ${2}bsmap2bismark.py ${1}

    python ${2}wgbs2bed.py\
        --infile ${1}.bis.txt\
        --outfile ${1}_coordinates.bed

    bedtools getfasta -tab -s -name\
        -fi ${3}\
        -bed ${1}_coordinates.bed\
        -fo ${1}_coordinates_sequences.txt

    python ${2}compute_flanks.py --depth 10\
        --infile ${1}_coordinates_sequences.txt\
        --outfile ${1}_flanks_methylation_d10.csv

    python ${2}compute_flanks.py --depth 5\
        --infile ${1}_coordinates_sequences.txt\
        --outfile ${1}_flanks_methylation_d5.csv
}

for i in ${SAMPLES[@]}; do analyse_flanks "$i" "$PATH_WGBS2FLANKS" "$PATH_GENOME_FASTA"; done

# GENERATE FILES CONTAINING CG COORDINATES AND FLANK SEQUENCES FOR LOCAL CORRELATION ANALYSIS

WRITE_COORDINATES_FLANKS=$(cat <<EOF

infiles = [
    "GSM3239884_SKO_mES_BS_seq_mCG.txt_coordinates_sequences.txt",
    "GSM4809269_SKO_Dnmt1_KO_BS_seq_mCG.txt_coordinates_sequences.txt"
    ]

for infile in infiles:
    outfile = infile.split(".")[0] + "_flanks.txt"
    nametag = infile.split(".")[0]
    outfile = nametag + "_flanks.txt"

    with open(infile, "r") as fh, open(outfile, "w") as rh:
        rh.write("\t".join(["chr", "start", "end", "5mC", "totalC", "sequence\n"]))
        line = fh.readline()
        while line:
            newline = []
            info, sequence = line.split("\t")
            elements = info.split("::")
            coordinates = elements[1].split(":")
            newline.append(coordinates[0])
            coor = coordinates[1].split("(")[0].split("-")
            newline.append(coor[0])
            newline.append(coor[1])
            newline.append(info.split(";")[0])
            newline.append(info.split(";")[1])
            newline.append(sequence.upper())
            rh.write("\t".join(newline))
            line = fh.readline()

    outfile2 = nametag + "_flanks_d5.txt"
    with open(outfile, "r") as fh, open(outfile2, "w") as rh:
        rh.write("\t".join(["chr", "start", "end", "5mC", "totalC", "sequence\n"]))
        line = fh.readline()
        line = fh.readline()
        while line:
            if int(line.split("\t")[4])>=5:
                rh.write(line)
            line = fh.readline()

EOF
)

python3 -c "$WRITE_COORDINATES_FLANKS"
