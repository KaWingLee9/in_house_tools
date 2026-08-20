

## Download all sra files under the Bioproject and convert to fastq file  
``` bash
# Before using, install sra-tools in home directory: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit#the-installation-processes-for-mac-os-x-and-the-two-linux-distributions-are-roughly-identical

export PATH=~/sratoolkit.3.4.1-alma_linux64/bin:$PATH
Accession=PRJNA1147209

# download all files under Bioproject
python DownloadBioproject.py ${Accession} --workers 6 --output-root ./

# convert the sra files to fastq files
python ConvertSraToFastq.py ${Accession} \
  --sra-dir ${Accession}/sra \
  --output-dir ${Accession}/fastq \
  --temp-dir ${Accession}/tmp \
  --metadata-dir ${Accession}/metadata \
  --jobs 2 \
  --threads-per-job 16
```

## generate md5 files in batch
ls ./ | while read file; do md5sum ${file} > ${file}.md5; done
