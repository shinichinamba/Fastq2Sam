#! /bin/bash
#$ -S /bin/bash
#$ -l s_vmem=5G,m_mem_free=5G

# @ST-E00104:1016:HVJLMCCXY:1:1101:11302:1309 1:N:0:NTATGGAT+NAATACAG
# /work9/SG_PVWGSfastq/201812_SG_HC_OS_00116_00158_SG_PV_TM_00059_00072/1901AHX-0002_hdd3/SG_HC_OS_00158/190129_HVJLM_SG_HC_OS_00158_L001_R1.fastq.gz
# /work9/SG_PVWGSfastq/201812_SG_HC_OS_00116_00158_SG_PV_TM_00059_00072/1901AHX-0002_hdd3/SG_HC_OS_00158/190129_HVJLM_SG_HC_OS_00158_L001_R2.fastq.gz

export LD_LIBRARY_PATH="/work1/home/snamba/conda_pack/gxx13/lib:$LD_LIBRARY_PATH"
SCAN() {
  local FASTQ_PFX="$1"
  local IDX="${2:-0}"
  \time \
    -f '==========\ntime: %esec (system:%S user:%U %P)\nmaxmem: %Mkb\ninput:%I (%r sockets) output:%O (%s sockets)\ncommand: %C\nexit: %x' \
    $NH/Projects/Fastq2Sam/build/fastq2sam -1 "${FASTQ_PFX}1.fastq.gz" -2 "${FASTQ_PFX}2.fastq.gz" -n sample1 -i "$IDX"
}

SCAN "/work9/SG_PVWGSfastq/201812_SG_HC_OS_00116_00158_SG_PV_TM_00059_00072/1901AHX-0002_hdd3/SG_HC_OS_00158/190129_HVJLM_SG_HC_OS_00158_L001_R"
SCAN "/work28/public_WGS/ENA/PRJEB11455/ERR1813570_" 1
SCAN "/work28/public_WGS/SRA_fastq_1/ERR2868229/ERR2868229_" 1
