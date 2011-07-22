#!/usr/bin/perl 
#===============================================================================
# /home/dsth/.vim/perl-support/templates/file-description.template
# dsth@cantab.net
use strict;
use warnings;

One of generated cumpolsory external data file:

1) NucleotideAceessions2NucleotideSeqIDs - some WGS_SCFLDs gain annotation and so their accessions and seqids
are not available in previous p2g files - generate from exisiting release - i.e. to WGS_SCFLD level - i.e. via
entrez web-access not ftp

this is a one-off file that should never change between submissions - i.e. they should only change if 
submitting a 'new' assembly (e.g. even just modifying any AGP data and clearly modifying contigs themselves)

    generate with?

Every submission generated cumpolsory external data files:

1) short_introns_file : Small Intron file - this need is required for genes with introns < 5/10bp? - gives 
the sequence that needs to be put into a .pep file and removes the intron. uses exception 
with supporting evidence if an INSD entry name is provided else uses poor quality sequence region qualifier

THIS FILE WILL IN PROBABILITY REQUIRE REGENERATION EACH SUBMISSION

    generate list with?

    get the evidence with...

2) p2g : TranslationStables2PeptideAccessionsNSeqIDs - allows mapping from 'vb' translation stable ids to 
Gb peptide Accessions and SeqIds - this is generated from p2g file - in past they have been 'inconsistent'
as to what to use so may be 'interesting'

=================================================================================

* With any re-submitted peptide you must submit the Gb Accession.

* With newer submissions the submitted SeqID (i.e. when an old one is not found its allocated) is 
  the translation stable_id or nuc_accession-stable_id when the gene is split across submission level 
  contigs e.g. for pest

    grep -P '^AGAP' AAAB01.May2011.p2g | head -2
    AGAP012829-PA   EDO64865        157021149       1       343     null-gene-symb \
        AgaP_AGAP012829 null-gene-syn   null-gene-id       CRA_x9P1GAV4PJ4  \
        AAAB01000388    19594990        1       AGAP012829-PA [Anopheles gambiae str. PEST]
    AGAP012928-PA   EDO64864        157021148       1       160     null-gene-symb  \
        AgaP_AGAP012928 null-gene-syn   null-gene-id       CRA_x9P1GAV4PJH \ 
        AAAB01000394    19595002        1       AGAP012928-PA [Anopheles gambiae str. PEST]

    thus in these cases its simple col1->col2 mapping

* However, in past they have generally NOT been translation stable_ids (in aegypti they are transcript 
  stable_ids?!? and in pest some examples are:

    grep -P '^e' AAAB01.May2011.p2g | head -2   
    ebiP2536        EAA01830        116133526       4       125     null-gene-symb  \
        AgaP_AGAP012847 null-gene-syn   null-gene-id       CRA_x9P1GAV4NYY AAAB01000105  \
        19594432        1       AGAP012847-PA [Anopheles gambiae str. PEST]
    ebiP2973        EAA01834        55247439        2       484     null-gene-symb  \
        AgaP_AGAP012800 null-gene-syn   null-gene-id       CRA_x9P1GAV4P6Q AAAB01000212  \
        19594640        1       AGAP012800-PA [Anopheles gambiae str. PEST]

    grep -P '^E' AAAB01.May2011.p2g | head -2   
    ENSANGP00000027317      EAL42448        116133531       2       87      null-gene-symb  \
        AgaP_AGAP012915 null-gene-syn      null-gene-id    CRA_x9P1GAV4NTC AAAB01000034  \
        19594291        1       AGAP012915-PA [Anopheles gambiae str. PEST]
    ENSANGP00000029855      EAU78074        116133530       1       55      null-gene-symb \
        AgaP_AGAP012805 null-gene-syn      null-gene-id    CRA_x9P1GAV4NUW AAAB01000057  \
        19594337        1       AGAP012805-PA [Anopheles gambiae str. PEST]
        
    grep -P '^[^eEA]' AAAB01.May2011.p2g | head -2
    agCP9167        EAA01835        55247437        2       115     null-gene-symb \
        AgaP_AGAP012596 null-gene-syn   null-gene-id       CRA_x9P1GAV4P7E AAAB01000223 \
        19594661        1       AGAP012596-PA [Anopheles gambiae str. PEST]
    agCP9255        EAA01837        157021153       4       176     null-gene-symb \
        AgaP_AGAP012597 null-gene-syn   null-gene-id       CRA_x9P1GAV4PC9 \ 
        AAAB01000300    19594815        1       AGAP012597-PA [Anopheles gambiae str. PEST]

=================================================================================

    Thus can either search for AGAP-[RP] each time or just have a here-doc in the p2g parser that
    prints the old versions and then appends all new ones - i.e. anything where the FIRST column
    matches a -P[A-Z]...

    in theory we now use the translation stable_id AS seqid - or when the gene is split across contigs 
    contig_name-stable_id (thus w/o multiple-translations) so its simple to pull the accessions from the 
    last p2g file - i.e. ALWAYS need accessions so the file must be generated each time - BUT the old 
    annoyingly formated ones should be one-off - thus perhaps separate them into two files

        1x one of file of poorly choosen seqids allowing translation stable_id mapping to accession/seqid

        1x every submission generated file from last p2g that just checks for appropriate seqids and pulls
        their accessions
    
pest required files - one-off files 'should' be same accross submissions

1) p2o : TranslationStables2OldLocusIds - the 'gene' locus IDs underwent changes in the history of pest so this is required for the appropriate old locus tags to be added

2) n2h : NucleotideAccessions2Haplotype - pest assembly was generated from heterozygous 'pest' hybrid?!? - there are alternative haplotypes that require adding 

by having blacklist on genes that're f'up can simply append old annotation into appropriate file instead of 
substituting it...

=head1 template file for tbl2asn

Submit-block ::= {
  contact {
    contact {
      name
        name {
          last "Hughes" ,
          first "Daniel" ,
          initials "S. T." ,
          suffix "" } ,
      affil
        std {
          affil "European Bioinformatics Institute" ,
          div "VectorBase" ,
          city "Cambridge" ,
          country "United Kingdom" ,
          street "Hinxton" ,
          email "dsth@cantab.net" ,
          phone "+44 01223-494678" ,
          postal-code "CB10 1SA" } } } ,
  cit {
    authors {
      names
        std {
          {
            name
              name {
                last "Hughes" ,
                first "Daniel" ,
                initials "S. T." } } } ,
      affil
        std {
          affil "European Bioinformatics Institute" ,
          div "VectorBase" ,
          city "Cambridge" ,
          country "United Kingdom" ,
          street "Hinxton" ,
          postal-code "CB10 1SA" } } ,
    date
      std {
        year 2011 ,
        month 01 ,
        day 28 } } ,
  subtype new  }

=cut

=head1 generate asn/gb format files

generate both asn (.sqn) files and genbank format files to use for qc...?!?

    #!/bin/bash
    /home/dsth/project/tools/parser/tbl2asn/linux.tbl2asn -t GenBank-Update-Template.sbt -p $1 -J -F p -U -V vb -Z discrep-report.txt

=cut

=head1 check discrepency report

check the Summary section in 'Discrepancy Report Results' file

=cut

=head1 get summaries of problems in .val files

get quick summaries...

=head2 general numbers of notes, errors, warnings...

    cat anopheles_gambiae_11May31/*.val | perl -MData::Dumper -0 -ne 'my $set; s/^(\w+)/$set{$1}++/egms;print Dumper \%set'; 

    $VAR1 = {
          'NOTE' => 761,
          'ERROR' => 574,
          'WARNING' => 13380
    };

=head2 categories of error/note/warnings etc., and number of occurences

    cat anopheles_gambiae_11May31/*.val | perl -MData::Dumper -0 -ne 'my $set; s/\w+?\s+\[(S[\w\.]*?)\]/$set{$1}++/egms;print Dumper \%set'

    $VAR1 = {
          'SEQ_FEAT.CDSmRNAXrefLocationProblem' => 12,
          'SEQ_FEAT.NotSpliceConsensusDonor' => 288,
          'SEQ_FEAT.InternalStop' => 227,
          'SEQ_FEAT.SeqLocOrder' => 4,
          'SEQ_FEAT.RareSpliceConsensusDonor' => 290,
          'SEQ_FEAT.NoStop' => 116,
          'SEQ_INST.InternalNsInSeqRaw' => 7374,
          'SEQ_INST.StopInProtein' => 227,
          'SEQ_FEAT.NotSpliceConsensusAcceptor' => 145,
          'SEQ_FEAT.ReplicatedGeneSequence' => 2,
          'SEQ_FEAT.CollidingGeneNames' => 13,
          'SEQ_FEAT.PartialProblem' => 5655,
          'SEQ_FEAT.MultipleGeneOverlap' => 112,
          'SEQ_FEAT.TransLen' => 2,
          'SEQ_INST.HighNContentPercent' => 28,
          'SEQ_FEAT.CDSmRNArange' => 220
    };

=head2 categories by type i.e. note/error/warning...

    cat anopheles_gambiae_11May31/*.val | perl -MData::Dumper -0 -ne 'my $set; s/^(\w+.+?\[S[\w\.]*?)\]/$set{$1}++/egms;print Dumper \%set'

    $VAR1 = {
          'NOTE: valid [SEQ_FEAT.RareSpliceConsensusDonor' => 290,
          'NOTE: valid [SEQ_FEAT.PartialProblem' => 469,
          'WARNING: valid [SEQ_FEAT.TransLen' => 2,
          'ERROR: valid [SEQ_FEAT.SeqLocOrder' => 4,
          'ERROR: valid [SEQ_INST.StopInProtein' => 227,
          'WARNING: valid [SEQ_FEAT.CDSmRNArange' => 220,
          'WARNING: valid [SEQ_FEAT.CDSmRNAXrefLocationProblem' => 12,
          'NOTE: valid [SEQ_FEAT.ReplicatedGeneSequence' => 2,
          'WARNING: valid [SEQ_INST.HighNContentPercent' => 28,
          'WARNING: valid [SEQ_FEAT.CollidingGeneNames' => 13,
          'WARNING: valid [SEQ_INST.InternalNsInSeqRaw' => 7374,
          'ERROR: valid [SEQ_FEAT.InternalStop' => 227,
          'WARNING: valid [SEQ_FEAT.NotSpliceConsensusDonor' => 288,
          'ERROR: valid [SEQ_FEAT.NoStop' => 116,
          'WARNING: valid [SEQ_FEAT.PartialProblem' => 5186,
          'WARNING: valid [SEQ_FEAT.MultipleGeneOverlap' => 112,
          'WARNING: valid [SEQ_FEAT.NotSpliceConsensusAcceptor' => 145
    };

=head2 ignore the pointless stuff 

    cat anopheles_gambiae_11May31/*.val | grep -v -P 'SEQ_FEAT.NotSpliceConsensusDonor|SEQ_INST.InternalNsInSeqRaw|SEQ_FEAT.RareSpliceConsensusDonor|SEQ_FEAT.PartialProblem|SEQ_FEAT.NotSpliceConsensusAcceptor' | perl -MData::Dumper -0 -ne 'my $set; s/\w+?\s+\[(S[\w\.]*?)\]/$set{$1}++/egms;print Dumper \%set'

    $VAR1 = {
          'SEQ_FEAT.CDSmRNAXrefLocationProblem' => 12,
          'SEQ_FEAT.InternalStop' => 227,
          'SEQ_FEAT.SeqLocOrder' => 4,
          'SEQ_FEAT.CollidingGeneNames' => 13,
          'SEQ_FEAT.NoStop' => 116,
          'SEQ_INST.StopInProtein' => 227,
          'SEQ_FEAT.MultipleGeneOverlap' => 112,
          'SEQ_FEAT.ReplicatedGeneSequence' => 2,
          'SEQ_FEAT.TransLen' => 2,
          'SEQ_INST.HighNContentPercent' => 28,
          'SEQ_FEAT.CDSmRNArange' => 220
    };



=cut

=head1 split exon genes

for pest there are split exon genes due to the fact that the agp is nonsense - just replace the appropriate
entries. in both cases they seem to just get one entry - i.e. not sure if CDS spans boundary?!? - just
used old entry from AgamP3.4.

the full files with the replaced sections for AgamP3.5 submission:

    /home/dsth/projects/GenBank/dev/GenBank_SplitExonGenes/AAAB01008980.tbl
    /home/dsth/projects/GenBank/dev/GenBank_SplitExonGenes/AAAB01008982.tbl

=head2 AGAP003293 from AAAB01008982 of AgamP3.4:

    <170479 >168776 gene
                            locus_tag       AgaP_AGAP003293
    <170479 >168776 mRNA
                            locus_tag       AgaP_AGAP003293
    170479  >168776 CDS
                            locus_tag       AgaP_AGAP003293

=head2 AGAP009563 from AAAB01008980 of AgamP3.4:

    116994  >114704 gene
                            locus_tag       AgaP_AGAP009563
    116994  116717  mRNA
    116643  115618
    115484  115314
    115229  115022
    114949  >114704
                            locus_tag       AgaP_AGAP009563
    116815  116717  CDS
    116643  115618
    115484  115314
    115229  115022
    114949  >114704
                            locus_tag       AgaP_AGAP009563

=cut

=head1 QC

=head2 comparison with previous submission

in order to compare with the previous submission download the current Gb entry to the annotation (WGS_SCFLD)
level - i.e. via Entrez web interface and not ftp (has raw data - with or without AGP - i.e. full version gives
you sequence and 'none' full version give coordinates for scaffold in contig cs).

generate file of sequences for peptides:

    perl PreviousSbmPeptideSeqExtraction.pl

run comparison script:

    perl Compare_current-ensembl_old-Gb_new-Asn_translations.pl

=cut

for particular errors just grep them from the submisssion dir and fix

