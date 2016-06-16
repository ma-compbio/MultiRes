# MultiRes

MultiRes is a pipeline for finding ancestral gene orders with duplications when
given a phylogenetic tree, a set of synteny blocks on the extant species, and a
set of contiguous ancestral regions on these synteny blocks.

MultiRes is distributed as a bash script which runs several Python routines in
sequence. MultiRes is free and open-source software. Code included in it may be
used in part or full by others.

## Usage
<code> . MultiRes.sh `<gene families>` `<synteny blocks>` `<tree>` 
              `<CAR file>` `<copy number threshold>` 
              `<segment length>` `<window length>` `<output dir>`
</code>

## Requirements
MultiRes uses the NetworkX module to work with graphs. It also includes the
weight computation module from ANGES (Jones et al. 2012), which requires numpy.


## Input
MultiRes takes as input the following data:
<ol>
<li> A phylogenetic tree in NHX format.</li>
<li> A set of gene families in the inferCARs format.</li>
<li> A set of synteny blocks in the inferCARs format.</li>
<li> Contiguous ancestral regions at the ancestor defined in terms of the
     synteny blocks.</li>
<li> An integer copy number threshold.</li>
<li> An integer segment length parameter L.</li>
<li> An integer window length paramter l, with l <= L.</li>
<li> An output directory.</li>
</ol>

Detailed explanation:
1. The tree must be in NHX format, with clearly labelled extant species. Branch
   lengths may be specified.  The ancestor must be given the name '@'.

2. The gene families must be specified in the inferCARs format, which is as
   follows.

   `><gene_id>`

   `<species_name>.<chromosome_name>:<start_position>-<stop_position> <orientation>`
   
    .
   
    .
     
    Here <gene_id>,<start_position> and <stop_position>  must be integers, and
    <orientation> must be one of '+' or '-'. The rest are strings.

3. The synteny blocks must also be specified in the inferCARs format (see gene
   family format for specification.  The lengths of the blocks (i.e.
   stop_position-start_position+1) are expected to be at least as long as those
   of the genes. 

4. Contiguous ancestral regions (CARs) may be specified in the PQ-tree format
   output by ANGES (Jones et al. 2012), or in the format output by MGRA
   (Aleksyev and Pevzner 2009). The CARs must be halved, i.e. each synteny block
   must be directed.

5. The copy number threshold is used to filter the gene families by excluding 
   those families which have high ancestral copy numbers, as computed by 
   Sankoff-Roussea parsimony (Csuros 2008).

6. Segment length
7. Window length

   Segment length and window length are parameters used by the main method to
   partition the CARs. Larger parameters cause a dramatic increase in time
   complexity for mammalian genomes, without significantly affecting the
   accuracy. The user is advised to try the method with small parameters
   relative to the length of the longest CAR. 

8. Self-explanatory


## Output
The final output of MultiRes, written in the output directory, is a file named
FINAL_GENE_ORDERS. The file contains a set of CARs, with names identical to the
CAR names in the original CAR file. The CARs consist of sequences of gene family
extremities: a gene family i has tail and head extremities identified by 2*i-1
and 2*i respectively. Furthermore, each extremity is labelled by a
`localization' label. For eg. 

\#CAR2

1212.0 1211.0 1208.0 1207.0 1206.0 1205.0 1204.0 1203.0 705.0 706.0 | 1202.0
1201.0 1200.0 1199.0 1198.0 1197.0 1034.0 1033.0 1196.0 1195.0 1194.0 1193.0
1192.0 1191.0 1190.0 1189.0 1188.0 1187.0 1186.0 1185.0 | .. $.

Here, 1212.0 1211.0 are extremities of the same gene family, identified by the
localization label '0'. A subsequence which is bordered by '|' on both sides is
understood to be in the correct position, but with missing adjacencies
connecting it to the two subsequences flanking it. The '$' marks the end of the
CAR.

## References

<ol>
<li> Alekseyev, Max, and Pavel A. Pevzner. "Breakpoint graphs and ancestral
genome reconstructions." Genome research (2009): gr-082784. </li>
<li> Csűrös, Miklós. "Ancestral reconstruction by asymmetric Wagner parsimony
over continuous characters and squared parsimony over distributions."
Comparative Genomics. Springer Berlin Heidelberg, 2008. 72-86. </li>
<li> Jones, Bradley R., et al. "ANGES: reconstructing ANcestral GEnomeS maps."
Bioinformatics 28.18 (2012): 2388-2390. </li>
<li> Maňuch, Ján, et al. "Linearization of ancestral multichromosomal genomes."
BMC bioinformatics 13.Suppl 19 (2012): S11. </li>
</ol>
