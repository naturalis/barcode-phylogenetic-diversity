 ## Metabarcoding data

### Sequences: the fastq.gz files

The workflow is designed to visualize the alpha and beta diversity of the sequences in this directory. They are the 
fungal ITS sequences acquired from soil sampled in seven locations in the Netherlands. Their format is demultiplexed 
FASTQ Phred33 paired-end, indicating they are not multiplexed anymore, contain quality information, were the product of 
Illumina sequencing, and require a 'Manifest' file (see below). Another feature of note is that every sample has a 
forward (R1) and reverse (R2) pair that will be merged eventually, producing fourteen initial files. They were zipped 
to reduce storage, which resulted in the .gz extension.

### Manifest: _4854_manifest.csv_

During the pipeline, the sequences will have to be imported into QIIME2 inside the Galaxy environment. To do so, 
the 'qiime2 tools import' tool needs a Manifest file in CSV format, containing the sample id, the absolute filepath 
(based on the UUID inside Galaxy), and the direction (forward/reverse). One of the bin scripts, _paste_and_manifest.py_,
pastes the data inside Galaxy, creates the Manifest file, and also pastes this Manifest file into Galaxy. It is meant 
to be a temporary file, but it is included here for clarity. 
