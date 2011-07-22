my $nt = qq{\n}.' 'x9;

=head1 Classes

=head2 Processing

=cut

=pod

B<Interfaces>

=cut

=head3 GffDoc::Activity::Role 

    Basic interface for GffDoc manipulating modules.

=cut

package GffDoc::Activity::Role;
use Mouse::Role;
has 'gffdoc' => (is => 'ro', isa => 'GffDoc', required => 1);
has 'log4' => (is => 'ro', required => 1);
has 'prob_genes' => (is => 'ro', isa => 'ArrayRef');
requires qw/run/; 

1;

=head3 GffDoc::Activity::Minimal::Role 

    Minimal interface for GffDoc manipulating modules (i.e. those ones that don't quite fit).

=cut

package GffDoc::Activity::Minimal::Role;
use Mouse::Role;
requires qw/run/; 

1;

=head2 Instances

    These are all defined in the separate modules where they are used e.g. GffDoc::Parser, GffDoc::Build...

=cut

#b/##### containers #####

=head2 Containers

=cut 

=pod

B<Main Container Convinience Class>

=cut

=head3 GffDoc

    The class whose instance forms the structure for holding gff feature objects and eventually gffGene objects. 
    To simply things and ensure some degree of typesafety they rather inelegantly have individual 
    getters/setters (where are Cpp Templates when you want them?!?) for the different types of array objects 
    it may hold (see next).

=cut

######

package GffDoc;
use Mouse;
has 'geneArray' => (is => 'ro', isa => 'GffDoc::Feature::Array::gene', writer => '_geneArray');
has 'mRNAArray' => (is => 'ro', isa => 'GffDoc::Feature::Array::mRNA', writer => '_mRNAArray');
has 'exonArray' => (is => 'ro', isa => 'GffDoc::Feature::Array::exon', writer => '_exonArray');
has 'CDSArray' => (is => 'ro', isa => 'GffDoc::Feature::Array::CDS', writer => '_CDSArray');
has 'GffGenes' => (is => 'ro', isa => 'GffDoc::Array::Gff::Gene', writer => '_GffGenes');
#y make this meta?!?
has 'polypeptideArray' => (is => 'ro', isa => 'GffDoc::Feature::Array::polypeptide', writer => '_polypeptideArray');
has 'GbAccession' => (is => 'ro', isa => 'Str', writer => '_GbAccession');
has 'GbSeqId' => (is => 'ro', isa => 'Str', writer => '_GbSeqId');
has 'SbmObj' => ( is => 'ro', isa => 'GenBank::Submission' );
has 'chromosome_location_str' => ( is => 'ro', default => q{} );
has 'sequence' => ( is => 'rw' );

=head1 dumpGff3 method?!? - perhaps multiple seq_regions at a time - or an append mode for all in one?!?

=cut 

#=head1 dumpTbl method - only ever for one seq region at a time

sub from_slice {

    my ($class, $slice, $SbmObj) = @_;

    Exception::GffDoc->throw( 
        stage => 'Type', 
        type => 'Genbank::Submission', 
        error => 'from_slice method requires a Genbank::Submission instance',
    ) if (!$SbmObj);

    #print qq{\nthe e! seq name is: }, $SbmObj->e_seq_region_name();

    #$hash{_NucleotideIds} = {
    #    toSeqId => \%NucleotideAccessions2NucleotideSeqIDs,
    #    toAccession => \%NucleotideSeqIDs2NucleotideAccessions,

    #y upwards projection

    #y seems to be problem here!?!?

    my $chr = q{};
    if ($SbmObj->chromosomes()) {
        $chr = $slice->project('chromosome')->[0]->[2]->seq_region_name ne 'UNKN' 
            ? '[chromosome='.$slice->project('chromosome')->[0]->[2]->seq_region_name.']' 
            : '[chromosome=Unknown]';
    }

    #/ argh!?!
    my $seq = $slice->seq;
    chomp($seq);

    my $gffdoc = GffDoc->new(SbmObj => $SbmObj, chromosome_location_str => $chr, sequence => $seq);

    #my $method = $SbmObj->e_seq_region_name();

    if ($SbmObj->e_seq_region_name() eq 'GbAccession') {
        $gffdoc->_GbAccession($slice->seq_region_name);
        $gffdoc->_GbSeqId($SbmObj->_NucleotideIds->{toSeqId}->{$slice->seq_region_name});
    } else {
        $gffdoc->_GbAccession($SbmObj->_NucleotideIds->{toAccession}->{$slice->seq_region_name});
        $gffdoc->_GbSeqId($slice->seq_region_name);
    }

#    print qq{\ne_seq_...: }.$SbmObj->e_seq_region_name();
#    print qq{\naccession: }.$SbmObj->_NucleotideIds->{toAccession}->{$slice->seq_region_name};
#    print qq{\nseqid: }.$SbmObj->_NucleotideIds->{toSeqId}->{$slice->seq_region_name};
#        print qq{\nseq_region_name: }.$slice->seq_region_name();

    $gffdoc->_GffGenes(GffDoc::Array::Gff::Gene->new());

    $gffdoc->GbAccession() || die qq{\nwe must have the genbank accession availble in the appropriate file};
    $gffdoc->GbSeqId() || die qq{\nwe must have the genbank SeqId availble in the appropriate file};

    my $ga = $slice->adaptor()->db()->get_GeneAdaptor();

    #r perhaps should re-sort ALL genes - i.e. including split ones with standard...

    #/ they are handled separately as the split gene processing has too much overhead to use generally
    if ($SbmObj 
          && $SbmObj->split_genes() 
          #&& $SbmObj->split_genes()->{toGeneId}->{$gffdoc->SeqId()}
          && $SbmObj->split_genes()->{toGeneId}->{$slice->seq_region_name()}
        ) {

        for my $split_gene_id (

          #@{$SbmObj->split_genes()->{toGeneId}->{$gffdoc->SeqId()}}
          @{$SbmObj->split_genes()->{toGeneId}->{$slice->seq_region_name()}}
            ) {

            next if (exists $SbmObj->blist()->{$split_gene_id});

            my $g = $ga->fetch_by_stable_id($split_gene_id) or die; # 4 small - strand

#        $hash{split_genes} = {
#            toSeqId => \%MultiSegmentGenes_map,
#            toGeneId => \%MultiSegmentGenes_rev_map,
#        };

#            my $time = time-$start;
#            print qq{\n\n[}.int(100*($split_gene_count+$gene_count)/$total).'% / '.$time.qq{ sec]\n} 
#            if (($split_gene_count+$gene_count)%100==0);
#            &_split_gene_processing($slice_name, $tbl_fh, $split_gene_stableid);
#            $split_gene_count++;
            #my ($class, $split_gene, $submission_cs, $seqid_with_part_of_gene) = @_;
            
            my $gffgene = GffDoc::Gff3::Gene->_partial_from_split_eGene(
              $g,$SbmObj->coordsystem,$slice->seq_region_name()
            );

#r/            print qq{\nshoving in split gene }, $g->stable_id();
#print qq{\nshoving in split gene [}.$g->strand.'] ', $g->stable_id();

            #/ ignore invalid genes
            $gffdoc->GffGenes()->add($gffgene) if $gffgene;

        }

    }

    for my $gene (sort {$a->start <=> $b->start} @{$slice->get_all_Genes}) {

            next if (exists $SbmObj->blist()->{$gene->stable_id()});
#        my $time = time-$start;
#        print qq{\n\n[}.int(100*($split_gene_count+$gene_count)/$total).'% / '.$time.qq{ sec]\n} 
#        if (($split_gene_count+$gene_count)%100==0);
            my $gffgene = GffDoc::Gff3::Gene->from_eGene($gene);
#        $gene_count++;

#r/            print qq{\nshoving in gene }, $gene->stable_id();
            $gffdoc->GffGenes()->add($gffgene);

    }

    return $gffdoc;
}

#mtbl

sub print_ContigTbl {

    my $self = shift;
    my $n = $Tools::n;
    my $t = $Tools::t;
    my $qual_kv_pair = $Tools::qual_kv_pair;

    my $SbmObj = $self->SbmObj();
    $SbmObj || die;

    #print qq{\n> processing slice corresponding to genbank accession: }, $self->GbAccession();

=head1

    #w pest this is accession for WGS level (submission level)
    #w aegypti this is Gb SeqID for WGS_SCFLD level (submission level)

    #b/ WGS e.g. aedes: WGS: AAGE02000001:AAGE02036206 #Scaffolds: CH477186:CH479178, CH899794:CH902558     
    # grabAegyptiAccessionsFromSeqIDnames.pl

    #y seemingly Gb Wgs take two levels: WGS and WGS_SCFLD. 
    #y in standard Wgs projects these correspond to WGS=contigs / WGS_SCFLD=supercontigs
    #y presumably send Gb the raw contigs with names (turned into Gb SeqIds) and the AGP (again with names...):
    #y both sets get assigned accession ids - but - in the case of the contigs they are given the Wgs 4 digit prefix
    #y of the project while the WGS_SCFL gets the standard 'arbitrary' Gb style accessions.
    #y annotation is provided at the WGS_SCFLD level

    #b/ PEST WGS:AAAB01000001:AAAB01069724 #Scaffold:CM000356:CM000360 

    #r PEST is old-style project w/ chromosomes-scaffold-chunks
    #r but it's a mess!?! but some scaffolds should really be single scaffolds
    #r further, the two levels seemingly submitted where chromosome and scaffold!?!
    #r thus scaffold have the standard 4 digit key accessions which the chromosomes have arbitrary ones
    #r the contigs don't figure into it at all - thus (1) you submit the raw scaffolds - i.e. with long runs of Ns
    #r (2) unlike with WGS as this is already a superstructure you submit annotation at the WGS and not the WGS_SCFLD
    #r level!?! - but still cos of the mess that is the scaffolds there are split genes!?! chunks have nothing to do
    #r with Gb submissions

    #b/ WGS: WGS=contigs                PEST: WGS=scaffold
    #b/ WGS: WGS_SCFLD=supercontigs     PEST: WGS_SCFLD=chromosomes

    #b/ annotation_submission: WGS: WGS_SCFLD       PEST: WGS

    #b/ e! pest scaffold seqregion_names: WGA_Accessions     aegypti: supercontig seqregion_names: supercont... (Gb SeqIDs)

    #b/ thus pest: need to look up Gb SeqIds - i.e. Celera names and for aegypti have Gb SeqIds (vb names) need to pull WGA_SCFLD accessions

    #w pest the scaffold names are the accession and we look up vb names (Gb seqid: celera...)
    #w aegypti supercontig names are vb names (supercont... - i.e. Gb seqid) so need to look up accession
    #w hence they are reversed!?! pest: accession->seqid / aegypti seqid->accession

=cut

#    ($slice_alternate_name,$slice_name) = ($slice_name,$slice_alternate_name) if ($proj_prefix ne 'AAAB'); # if (!$wgs); ?!?

    #y annotation submissions: don't submit anything without annotation
#    next if ((scalar @{$slice->get_all_Genes} == 0) && (!exists $MultiSegmentGenes_rev_map{$slice_name}));

#    print qq{\nprocessing contig $slice_name :}, $slice->strand;

    open my $fa_fh, '>', $SbmObj->_outdir().'/'.$self->GbAccession().'.fsa' or die;
    open my $tbl_fh, '>', $SbmObj->_outdir().'/'.$self->GbAccession().'.tbl' or die;
#    $tbl_fh = *STDOUT;
#    my $fa_fh = *STDOUT;

    my $spc_str = ' [organism='.$self->SbmObj()->species_strain().'] [tech=wgs] '.$self->chromosome_location_str();

    my $string = 'gnl|WGS:'.$self->SbmObj()->wgs_proj_prefix().'|'.$self->GbSeqId().'|gb|'.$self->GbAccession();

#   'AAAB01008839' => [
#                              {
#                                'comment' => 'a probable alternative assembly of this region is represented by part of scaffold AAAB01008980',
#                                'end' => '109822',
#                                'start' => '1'
#                              },
#                              {
#                                'comment' => 'a probable alternative assembly of this region is represented by scaffold AAAB01008837',
#                                'end' => '367413',
#                                'start' => '250843'
#                              },
#                              {
#                                'comment' => 'a probable alternative assembly of this region is represented by scaffold AAAB01003176',
#                                'end' => '2027030',
#                                'start' => '1994213'
#                              }
#                            ]
#        };

    my $hap_str = q{};
    #/ forgot to do this!?!
    if ($self->SbmObj()->n2h() && exists $self->SbmObj()->n2h->{$self->GbAccession()}) {
        for my $i (@{$self->SbmObj()->n2h->{$self->GbAccession()}}) {
            $hap_str .= $i->{start}.$t.$i->{end}.$t.'misc_feature'.$n;
            $hap_str .= $qual_kv_pair.'note'.$t.$i->{comment}.$n;
        }
    }

    #y print the main GB record entry for the contig/scaffold
    print $tbl_fh '>Feature '.$string.$n.$hap_str;
    print $fa_fh '>'.$string.$spc_str.$n.$self->sequence();

    close $fa_fh;

    for my $gffgene (@{$self->GffGenes()->features()}) {


#/debug!?! 
#next if ($gffgene->id ne 'AGAP001127');

        print $tbl_fh $gffgene->toTblString($self->SbmObj(), $self->GbAccession());
        $SbmObj->_increment_processed_genes();

    }

    close $tbl_fh;

    #print qq{\nfinished processing slice\n};

    return;
 
}

__PACKAGE__->meta->make_immutable();

1;

=pod

B<Interfaces>

=cut

=head3 GffDoc::Feature::Array::Role;

    Minimal interface for the different 'type' arrays that hold the various feature types. As it's a mouse role
    it also provides the meat of some methods e.g. add, distinctAsHashRef, findParentsDestructive, 
    findParentsNonDestructive, length. 

    Because the 'add' method circumvents the individual auto-generated constructors for the different class instances
    and therefore slips past the normal mouse typesafety checks, before is used to perform the appropriate checks in 
    each implementing class.

=cut

package GffDoc::Feature::Array::Role;
use Mouse::Role; 
requires qw/add/; 
has 'features' => (is => 'ro', isa => qq{ArrayRef}, init_arg => undef, writer => '_features'); # auto_deref => 1,);

#y appalingly inefficient...
#sub shift {
#    #my $self = shift;
#    my ($self) = @_;
#    return shift @{$self->features()};
#}

sub length {
    return scalar @{shift->features()};
}

sub distinctAsHashRef {

    my ($self,$what) = @_;

    my $type = ref $self;
    $type =~ s/Array/Gff3/;
    my $meta = $type->meta();
    my $ref = ref $self;

    Exception::GffDoc->throw( 
        stage => 'Type', 
        type => 'NoSuchAttrib', 
        error => 'what on earth is a '.$what.'?',
    ) if (!$meta->has_attribute($what));

    Exception::GffDoc->throw ( 
        stage => 'Type', 
        type  => 'Build::EmptyFeatureArray', 
        error => 'Empty feature array: '.$ref
    ) if (!$self->features());

    Exception::GffDoc->throw ( 
        stage => 'Type', 
        type  => 'Build::EmptyFeatureArray', 
        error => 'Empty feature array: '.$ref
    ) if (@{$self->features()} == 0);
    
    my %distinct;
    @distinct{map { $_->$what } @{$self->features()} } = ();
    return \%distinct;
}

sub add { 
    push @{$_[0]->{features}}, $_[1];
}

sub findParentsDestructive {

    my ($self,$parentid) = @_;
    Exception::GffDoc->throw ( 
        stage => 'Type', 
        type  => 'Build::EmptyFeatureArray', 
        error => 'feature array '.ref $self.' is empty!' 
    ) if (!$self->features());
    Exception::GffDoc->throw ( 
        stage => 'Type', 
        type  => 'Build::EmptyFeatureArray', 
        error => 'feature array '.ref $self.' is empty!',
    ) if (@{$self->features()} == 0);

    my @features = grep { $_->parent ne $parentid } @{$self->features()};
    my @return = grep { $_->parent eq $parentid } @{$self->features()};

    #y need to allow for pseudogenes
    #Exception::Build::FeatureLine::Dangling::NoChildren->throw ( error => qq{Feature $parentid exists but has no }
    #  .'lower level features referencing it. Are you sure that ID= and Parent= identifiers are consistent?!' )
    #  if (scalar @return == 0);
    
    $self->_features(\@features);
    return \@return;
}

sub findParentsNonDestructive {

    my ($self,$parentid) = @_;
    Exception::Build::EmptyFeatureArray->throw ( error => 'feature array '.ref $self.' is empty!' ) 
      if (!$self->features());
    Exception::Build::EmptyFeatureArray->throw ( error => 'feature array '.ref $self.' is empty!' ) 
      if (@{$self->features()} == 0);

    return [grep { $_->parent eq $parentid } @{$self->features()}];
}

1;

=pod

B<Instances>

=cut

=head3 GffDoc::Feature::Array::gene

    Array container class for holding GffDoc::Feature::Gff3::gene instances 'with 'GffDoc::Feature::Array::Role'. Provides a shift method.

=cut

package GffDoc::Feature::Array::gene;
use Mouse;
with 'GffDoc::Feature::Array::Role';
{ my $type = 'GffDoc::Feature::Gff3::gene';
before 'add' => sub {
    Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::WrongType', error => 'i want a '.$type.' object' ) 
      if (!$_[1]->isa($type));
}; }

sub shift {
    my $self = shift;
    return shift @{$self->features()};
}

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Array::mRNA

    Array container class for holding GffDoc::Feature::Gff3::mRNA instances 'with 'GffDoc::Feature::Array::Role'.

=cut

package GffDoc::Feature::Array::mRNA;
use Mouse;
with 'GffDoc::Feature::Array::Role';
{ my $type = 'GffDoc::Feature::Gff3::mRNA';
before 'add' => sub {
    Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::WrongType', error => 'i want a '.$type.' object' ) 
      if (!$_[1]->isa($type));
}; }

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Array::exon

    Array container class for holding GffDoc::Feature::Gff3::exon instances 'with 'GffDoc::Feature::Array::Role'.

=cut

package GffDoc::Feature::Array::exon;
use Mouse;
with 'GffDoc::Feature::Array::Role';
{ my $type = 'GffDoc::Feature::Gff3::exon';
before 'add' => sub {
    Exception::GffDoc->throw( 
        stage => 'Type', 
        type  => 'Type',
        error => 'play fair. '.$_[1].' is not a '.$type,
    ) if (!$_[1]->isa($type));
}; }

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Array::CDS

    Array container class for holding GffDoc::Feature::Gff3::CDS instances 'with 'GffDoc::Feature::Array::Role'.

=cut

package GffDoc::Feature::Array::CDS;
use Mouse;
with 'GffDoc::Feature::Array::Role';
{ my $type = 'GffDoc::Feature::Gff3::CDS';
before 'add' => sub {
    Exception::GffDoc->throw( 
        stage => 'Type', 
        type  => 'Type',
        error => 'play fair. '.$_[1].' is not a '.$type,
        ) if (!$_[1]->isa($type));
}; }

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Array::polypeptide

    Array container for holding GffDoc::Feature::Gff3::polypeptide instances. 

=cut

package GffDoc::Feature::Array::polypeptide;
use Mouse;
with 'GffDoc::Feature::Array::Role';
{ my $type = 'GffDoc::Feature::Gff3::polypeptide';
before 'add' => sub {
    Exception::Type->throw( error => 'play fair. '.$_[1].' is not a '.$type ) 
      if (!$_[1]->isa($type));
}; }

__PACKAGE__->meta->make_immutable();

1;

#b/##### features #######

=head2 Feature Classes

=cut

=pod

B<Interfaces>

=cut

=head3 GffDoc::Feature::Role

    Generic interface for feature class core functionality. Generally stores the main feature 
    data as read-only accessors: seqid, source, strand, start, end, phase, score, biotype, id.

    The constructor performs a few basic checks like that start is <= end, and strand is '+' or '-' etc..

    It also provides a copy2Type method. This simply takes a Feature object and deep copies its contents
    (the copy process isn't hard coded due to lazyness - i.e. i don't want to wonder why the code breaks
    when new feature attribs are added to a class...).

=cut

package GffDoc::Feature::Role;
use Mouse::Role; 
#use Mouse::Util::TypeConstraints;
has [qw/seqid/] => (is => 'ro', required => 1, isa => 'Str', writer => '_seqid_for_partial');
has 'source' => (is => 'ro', isa => 'Str');
has 'strand' => (is => 'ro', required => 1);
has 'start' => (is => 'ro', required => 1, isa => 'Int', writer => '_start' ); #y completely ignore phase anyway
has 'end' => (is => 'ro', required => 1, isa => 'Int', writer => '_end' ); #y completely ignore phase anyway
#has 'phase' => (is => 'rw', isa => enum([qw/ 1 2 3 . /]) );
has 'phase' => (is => 'rw'); #y completely ignore phase anyway
has 'score' => (is => 'ro'); #y completely ignore...
has 'biotype' => (is => 'ro', writer => '_biotype');
has 'id' => (is => 'ro', isa => 'Str', writer => '_id');

our $additional_types_string = q{};

around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    my %hash = @_;
    
    Exception::GffDoc::FeatureConstructor->throw( 
        error => 'hmmm '.$hash{id},
        hash => \%hash,
        stage => 'Parser', 
        type => 'FeatureLine::WTF', 
    ) if (!exists $hash{id} && !exists $hash{ID});
    #y unnecessary, but just making sure things throw early...
    #Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Coords', error => 'start must be less than end! feature: '.Dumper \%hash) 
    # Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Coords', error => 'start must be less than end! feature: '.$hash{id}) 
    Exception::GffDoc::FeatureConstructor->throw( 
        error => 'hmmm '.$hash{id},
        hash => \%hash,
        stage => 'Parser', 
        type => 'FeatureLine::WTF', 
    ) if ($hash{seqid} eq q{});

#    Exception::GffDoc::FeatureConstructor->throw( 
#        error => 'hmmm '.$hash{id},
#        hash => \%hash,
#        stage => 'Parser', 
#        type => 'FeatureLine::WTF', 
#    ) if ($hash{id} eq q{});

#y mechanistically all these are pseudogenes... - i.e. only safe thing is force explicit activation of protein_coding procedures...
# transposable_element_gene=gene:transposable_element -type pseudogenic_exon=exon -type protein=ignore -type pseudogenic_transcript=mRNA:pseudogene -type ncRNA=mRNA:ncRNA -type tRNA=mRNA:tRNA -type miRNA=mRNA:miRNA

    Exception::GffDoc::FeatureConstructor->throw( 
        error => 'have not yet implemented biotype '.$hash{biotype}.' for '.$hash{id}.$nt
        .'If you want to allow for additional non_protein_coding biotypes use -non_protein_coding_types with a list of comma-separated types',
        hash => \%hash,
        stage => 'Parser', 
        type => 'FeatureLine::Biotype', 
    ) if ($hash{biotype} !~ /^(protein_coding|pseudogene|rRNA|ncRNA|tRNA|miRNA${additional_types_string})$/);

#    Exception::GffDoc::FeatureConstructor->throw( 
#        error => 'start must be less than end!', 
#        hash => \%hash,
#        stage => 'Parser', 
#        type => 'FeatureLine::Coords', 
#    ) if ($hash{start} >= $hash{end});
    
    Exception::GffDoc::FeatureConstructor->throw( 
        error => 'strand must be + or - feature: '.$hash{id},
        hash => \%hash,
        stage => 'Parser', 
        type => 'FeatureLine::Strand', 
    ) if ($hash{strand} ne '+' && $hash{strand} ne '-');

    return $class->$orig(@_);
};

sub copy2Type {
    my ($self,$type) = @_;
    Exception::Type->throw( error => 'what is a '.$type )
      if ($type ne 'gene' && $type ne 'mRNA' && $type ne 'exon' && $type ne 'CDS');
    my $this = do { +{%{$self}} }; # ~= my $this = {%{$self}};
    my $copy = &_deep_copy($this);
    my $class = 'GffDoc::Feature::Gff3::'.$type;
    return bless $copy, $class;
}

sub _deep_copy {
    my $this = shift;

    if (not ref $this) { 
        $this;
    } elsif (ref $this eq "ARRAY") {
        [map &_deep_copy($_), @$this];
    } elsif (ref $this eq "HASH") {
        +{map { $_ => &_deep_copy($this->{$_}) } keys %$this};
    } else { die "copy: what type is $_?" }
}

#no Mouse::Util::TypeConstraints;
1;

=head3 GffDoc::Feature::Gff::Leaf::Role

    Poorly named interface for terminal gff 'or gtf feature types (CDS/exon). 

    Provides parent and parentsparent attributes.

    The constructor handles id, parent's id and parentsparent id extraction - you hand the 
    constructor the attribs column string as 'id' and the 'gtf' format parameter (true/false) 
    and the constructor does the rest - i don't want code handling gtf vs. gff tokens messing up
    bit that that are just trying to extract features from a doc hence this is done here (yeah, 
    yeah, prolly should be 'string2names' method, DIY!).

=cut

package GffDoc::Feature::Gff::Leaf::Role;
use Mouse::Role;
has 'parent' => (is => 'ro', isa => 'Str', writer => '_parent'); 
has 'parentsparent' => (is => 'ro', isa => 'Str', writer => '_parentsparent');

around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    my %hash = @_;
    #y if you modify id in situ first you fuck the rest... duh.

    #g getting rid of separate class, using meta to add attrib and buildargs to modify construction

    if (!$hash{gtf}) {

        #y this gets manually reset with exploded exons
        if ($hash{id} =~ /Parent=\s?"?(\S+?)"?\s?;/) { $hash{parent} = $1;
        } else { 
            Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::ParentId', 
              error => 'gff CDS/exon parent identifiers must be of gff3 format (Parent=blah)):'
              .$nt.$hash{id}.' - if they have gtf identifiers you should use the -gtf option'
            ) 
        }

        #y don't really care about id for exon/CDS
        if ($hash{id} =~ /ID=\s?"?(\S+?)"?\s?;/) { $hash{id} = $1; }
        return $class->$orig(%hash);

    } else {

        if ($hash{id} =~ /transcript_id\s?"?(\S+?)"?\s?;/) { 
            $hash{parent} = $1; 
        } else { 

            Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::ParentId', 
              error => 'gtf CDS/exon parent identifiers must be of gtf format (transcript_id "blah"): '
              .$nt.$hash{id}
        ) 

        } if ($hash{id} =~ /gene_id\s?"?(\S+?)"?\s?;/) { 
            $hash{parentsparent} = $1; 
        } else { 

            Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::ParentsParentId', 
              error => 'gtf CDS/exon parentsparent identifiers must be of gtf format (gene_id "blah"): '
              .$nt.$hash{id}) 
        }
        #y don't give a...
        if ($hash{id} =~ /ID=\s?"?(\S+?)"?\s?;/) { $hash{id} = $1; 
        } else { $hash{id} = 0; }
        #} else { $hash{id} = undef; }
        #} else { $hash{id} = 'NULL'; }

        return $class->$orig(%hash);

    }
};

1;

=pod

B<Instances>

=cut

=head3 GffDoc::Feature::Gff3::gene

    Class representing gff gene features 'with GffDoc::Feature::Role'

    As with GffDoc::Feature::Gff::Leaf::Role constructor is passed the attrib column as 'id' parameter
    and extracts the actual feature id.

=cut

package GffDoc::Feature::Gff3::gene;
use Mouse;
with 'GffDoc::Feature::Role';

around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    my %hash = @_;
    if ($hash{ID}) { 
        $hash{id} = $hash{ID};
    } elsif ($hash{id} =~ /ID=\s?(\S+?)\s?;/) { 
        $hash{id} = $1; 
    } else { 
        Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Identifier', error => 'gene identifiers must be of gff3 format (ID=blah): '.$nt.$hash{id}) 
    }
    return $class->$orig(%hash);
};

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Gff3::mRNA

    Class representing gff mRNA features 'with GffDoc::Feature::Role'

    As with GffDoc::Feature::Gff::Leaf::Role constructor is passed the attrib column as 'id' parameter
    and extracts the actual feature id 'and' parent id.

=cut

package GffDoc::Feature::Gff3::mRNA;
use Mouse;
with 'GffDoc::Feature::Role';
has 'parent' => (is => 'ro', isa => 'Str', writer => '_parent'); 
#has 'polypeptide_seq' => (is => 'rw', isa => 'Str');
#has 'gene' => (is => 'rw', isa => 'GffDoc::Feature::Gff3::Gene', required => 1);
#has 'exons' => (is => 'rw', isa => 'ArrayRef[GffDoc::Feature::Gff3::exon]', required => 1);
#has 'CDS' => (is => 'rw', isa => 'ArrayRef[GffDoc::Feature::Gff3::CDS]', required => 1);

around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    my %hash = @_;
    #y if you modify id in situ first you fuck the rest... duh.
    if ($hash{id} =~ /(Parent=|gene_id)\s?"?(\S+?)"?\s?;/) { $hash{parent} = $2; 
    } else { Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::ParentId', error => 'mRNA parent identifiers  must be of gff3 format (Parent=blah) or gtf format (gene_id "blah"): '.$nt.$hash{id}) }
    if ($hash{id} =~ /(ID=|transcript_id)\s?"?(\S+?)"?\s?;/) { $hash{id} = $2; 
    } else { Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Identifier', error => 'mRNA identifiers must be of gff3 format(ID=blah) or gtf format (transcript_id "blah"): '.$nt.$hash{id}) }
    return $class->$orig(%hash);
};

__PACKAGE__->meta->make_immutable();

1;
 
=head3 GffDoc::Feature::Gff3::exon

    Class representing gff exon features 'with GffDoc::Feature::Role' and 'with GffDoc::Feature::Gff::Leaf::Role'

    Also holds ensembl equivalent of exon phase (for what it's worth...).

=cut

package GffDoc::Feature::Gff3::exon;
use Mouse;
with 'GffDoc::Feature::Role';
with 'GffDoc::Feature::Gff::Leaf::Role';
#y these are daft, but that's the system...
has 'ephase_start' => (is => 'ro', isa => 'Int', writer => '_ephase_start');
has 'ephase_end' => (is => 'ro', isa => 'Int', writer => '_ephase_end');
has 'eExon' => (is => 'ro', isa => 'Bio::EnsEMBL::Exon', writer => '_eExon');
#has 'mRNA' => (is => 'rw', isa => 'GffDoc::Feature::Gff3::mRNA', required => 1);

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Gff3::CDS

    Class representing gff CDS features 'with GffDoc::Feature::Role' and 'with GffDoc::Feature::Gff::Leaf::Role'

=cut

package GffDoc::Feature::Gff3::CDS;
use Mouse;
with 'GffDoc::Feature::Role';
with 'GffDoc::Feature::Gff::Leaf::Role';
#has 'mRNAs' => (is => 'rw', isa => 'ArrayRef[GffDoc::Feature::Gff3::mRNA]', required => 1);

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Feature::Gff3::CDS

    Class representing a minimal set of info for integrating polypeptide data when required.
    At the moment it is only used for when CDS data is implicitly defined via polypeptides, but it 
    also has id and seq attributes for if-and-when external modules are written to allow for 
    translation names to be pulled this way and/or seqedits to be directly integrated at import time.

=cut

package GffDoc::Feature::Gff3::polypeptide;
use Mouse;
has 'parent' => (is => 'ro', isa => 'Str', writer => '_parent'); 
has 'id' => (is => 'ro', isa => 'Str', writer => '_id'); 
has 'start' => (is => 'ro', isa => 'Str', writer => '_start'); 
has 'end' => (is => 'ro', isa => 'Str', writer => '_end'); 
has 'seq' => (is => 'ro', isa => 'Str', writer => '_seq'); 

around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    my %hash = @_;
    if ($hash{id} =~ /ID=\s?(\S+?)\s?;/) { $hash{id} = $1; 
    } else { Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Identifier', error => 'gene identifiers must be of gff3 format (ID=blah): '.$nt.$hash{id}) }
    return $class->$orig(%hash);
};

__PACKAGE__->meta->make_immutable();

1;

#b/##### composed classes for final genes

=head3 GffDoc::Arrayinternacional/::Gff::Gene

    Class' sole existence is to give a different structure to put 'finished' gff genes into - so that mouse bits me when 
    i tiredly hand it the wrong thing. with GffDoc::Feature::Array::Role.

=cut

package GffDoc::Array::Gff::Gene;
use Mouse;
with 'GffDoc::Feature::Array::Role';
{ my $type = 'GffDoc::Gff3::Gene';
before 'add' => sub {
    Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::WrongType', error => 'i want a '.$type.' object' ) 
      if (!$_[1]->isa($type));
}; }

__PACKAGE__->meta->make_immutable();

1;

package FakeExon;

sub new {
    my ($caller) = @_;
    #my ($caller,$start,$end) = @_;
    my $class = ref($caller) || $caller;
    my $self = bless [], $class;
    #my $self = bless [$start,$end], $class;
    return $self;
}

#y accessors/setters to emulate the exon objects
sub start {
    my ($self,$start) = @_;
    if ($start) {
        $self->[0] = $start;
        return;
    }
    else {
        return $self->[0];
    }
}

sub end {
    # my $self = shift;
    # return $self->[1];
    my ($self,$end) = @_;
    if  ($end) {
        $self->[1] = $end;
        return;
    }
    else {
        return $self->[1];
    }
}

1;

=head3 GffDoc::Gff3::Gene

    Composed class whose sole reason for existence is to have gene feature objects as different type to built gene objects.
    with GffDoc::Feature::Role.

=cut

#######

package GffDoc::Gff3::Gene;
use Mouse;

######################## hashofshorts should be map of id => seq - so never have to store sequence itself
######################## hashofshorts should be map of id => seq - so never have to store sequence itself
######################## hashofshorts should be map of id => seq - so never have to store sequence itself
######################## hashofshorts should be map of id => seq - so never have to store sequence itself
######################## hashofshorts should be map of id => seq - so never have to store sequence itself
######################## hashofshorts should be map of id => seq - so never have to store sequence itself
######################## hashofshorts should be map of id => seq - so never have to store sequence itself

with 'GffDoc::Feature::Role';
#has 'mRNAs' => (is => 'rw', isa => 'ArrayRef[GffDoc::Feature::Gff3::mRNA]', required => 1, auto_deref => 1);
has 'mRNAs' => (is => 'rw', isa => 'ArrayRef[GffDoc::Gff3::mRNA]', auto_deref => 1);
has 'five_utr' => (is => 'ro', default => q{}, writer => '_five_utr');
has 'three_utr' => (is => 'ro', default => q{}, writer => '_three_utr');
has 'external_name' => (is => 'ro', default => q{}, writer => '_external_name');
has 'description' => (is => 'ro', default => q{}, writer => '_description');
has 'split' => (is => 'ro', default => 0, writer => '_split');

#mgene

sub toTblString {

    my ($self, $SbmObj, $GbAccession) = @_;

    #r either put reference to parent gffdoc in gene - ewww or just pass value if split
    #r or just pass it...

    #my ($self, $SbmObj,$split_flag) = @_;
    $SbmObj||=0;
    my $split_flag = $self->split() ? 1 : 0;

    die if ($SbmObj && !$SbmObj->isa('GenBank::Submission'));

    my $full_string = q{};

    my $nucAccession = $self->seqid();
    my $stable = $self->id;
    my $strand = $self->strand;
    my $start = $self->start;
    my $end = $self->end;
    my $biotype = $self->biotype;
    my $prefix = $SbmObj ? $SbmObj->locus_tag_prefix() : 'ensemblgenomes_';
    
    #my $ncRNA_types = $SbmObj->ncRNAs();
    my $ncRNA_types = qr{antisense_RNA|autocatalytically_spliced_intron  
        |ribozyme|hammerhead_ribozyme|RNase_P_RNA     
        |RNase_MRP_RNA|telomerase_RNA|guide_RNA       
        |rasiRNA|scRNA|siRNA|miRNA|piRNA|snoRNA|snRNA
        |SRP_RNA|vault_RNA|Y_RNA}x; # must have this format with x modifier...

    my $seqid = $self->seqid();
    my $id = $self->id();

    #r values only ever set once...
    my $n = $Tools::n;
    my $t = $Tools::t;
    my $qual_kv_pair = $Tools::qual_kv_pair;

    if ($biotype =~ 'protein_coding|pseudogene|tRNA|rRNA|misc_RNA|tmRNA' 
      || $biotype =~ /^$ncRNA_types$/x) {

    #mysql> select distinct biotype from gene;
    #| miRNA          |
###    #| snRNA          |
###    #| snoRNA         |

        my $gstring = $self->strand() eq '+' 
        # my $gstring = $self->strand() == 1 
          ? $self->five_utr().$start.$t.$self->three_utr().$end.$t.'gene'.$n 
          : $self->five_utr().$end.$t.$self->three_utr().$start.$t.'gene'.$n;

        $gstring .= $qual_kv_pair.'pseudo'.$n if ($self->biotype eq 'pseudogene');

        my $locus_tag = $qual_kv_pair.'locus_tag'.$t.$prefix.$stable.$n;
        $gstring .= $locus_tag;

        my $rRNA_product = q{};

        #print qq{\nbiotype: }.$biotype;

        if ($self->external_name()) {
            #g make sure we don't give same name more than once in a project...
            #r this is real submission dependent!?! - should we pass it to package or just inject it into the namespace?!?
            my $qual = $biotype eq 'rRNA' ? 'note' : 'gene';
            my $name = $SbmObj ? $SbmObj->AvoidRepeatNames($self->external_name()) : $self->external_name();
            $gstring .= $qual_kv_pair.$qual.$t.$name.$n; 

            # rna/AAAB01000012.tbl:                   note    5_8S_rRNA_a
            # rna/AAAB01000148.tbl:                   note    SSU_rRNA_5_a
            # rna/AAAB01008967.tbl:                   product 5.8S ribosomal RNA
            # rna/AAAB01008984.tbl:                   product 5S ribosomal RNA
            if ($biotype eq 'rRNA') {

                my $tprod = $self->external_name() eq '5_8S_rRNA'   ? '5.8S ribosomal RNA'
                          : $self->external_name() eq 'SSU_rRNA_5'  ? '5S ribosomal RNA'
                          :                                           'ribosomal RNA';

                $rRNA_product = $qual_kv_pair.'product'.$t.$tprod.$n 
            }
        }

        #g to avoid issues of multiple SeqIDs for the same translation but on different contigs we use compound names...
        my $compound_name = q{};
        my $old_protein_ref;

        if ($SbmObj) {
            #print qq{\nobject: }.$SbmObj; die;
            #use Data::Dumper; print Dumper $SbmObj->OldLocusTag($self->seqid().$SbmObj->_field_sep().$self->id()); die;
            # $compound_name = $seqid.$SbmObj->_field_sep().$id.'-PA';
            
            #$gstring .= $SbmObj->OldLocusTag($compound_name);
            $gstring .= $SbmObj->OldLocusTag($self->seqid(),$self->id());



            #$SbmObj->vectorbase_dbxref() && ($gstring .= $qual_kv_pair.'db_xref'.$t.'VectorBase:'.$id.$n);
            $gstring .= $qual_kv_pair.'db_xref'.$t.'VectorBase:'.$id.$n if ($SbmObj->vectorbase_dbxref());
        }

        $full_string .= $gstring;

        for my $mRNA (sort {$a->id cmp $b->id} @{$self->mRNAs}) {
            
            #r/ this needs revamping - i.e. short intron hash needs to have ->{seq} and ->{evidence}

            #r does it have a short intron?
            my $short_intron_support = ($SbmObj 
              && $mRNA->biotype eq 'protein_coding' 
              && $SbmObj->CheckHashOfShortIntrons($mRNA->trnsl_name())) ? $SbmObj->CheckHashOfShortIntrons($mRNA->trnsl_name()) : 0;
              #&& $SbmObj->CheckHashOfShortIntrons($mRNA->trnsl_name())) ? 1 : 0;

#            next if (!$short_intron_support); # TEMP DEBUG

            my $ncRNA_qualifier = q{};
            my $feature_type = q{};
            my $transcript_peptide_qualifier = q{};
            my $seqid_accession_valuestring;

            #y get the qualifiers for protein_coding and pseudogenes
            if ($biotype eq 'protein_coding' || $biotype eq 'pseudogene') {
                ($transcript_peptide_qualifier,$seqid_accession_valuestring) = &_protein_coding_string($mRNA,$SbmObj,$split_flag); 
                $feature_type = 'mRNA';
            }
            #r/ rRNA isn't correct!?!
            elsif ($biotype eq 'tRNA' || $biotype eq 'rRNA' || $biotype eq 'misc_RNA' || $biotype eq 'tmRNA') {
                $feature_type = $biotype;
            }
            #r/ ncRNAs are a big grouping of RNA types
            elsif ($biotype =~ /^$ncRNA_types$/x) {
                $feature_type = 'ncRNA';
                $ncRNA_qualifier = $qual_kv_pair.'ncRNA_class'.$t.$biotype.$n;
            } else {

                die qq{\nWTF: unrecognised type: $biotype};
            }

            # need locally scoped copy to stop it just getting constantly appended to...
            my $full_qualifiers = $locus_tag;
            $full_qualifiers .= $transcript_peptide_qualifier;

            my $mRNA_id = $mRNA->id();

#            my $exons_to_use_mrna = $short_intron_support ? &_fake_big_exon($mRNA) : $mRNA->exons();

            #y print the mRNA section
            my $coords_full_qualifiers = &_grab_coords(

                $strand,

                $short_intron_support ? &_fake_big_exon($mRNA):$mRNA->exons(),
#               $exons_to_use_mrna,

                #!$trans->five_prime_utr ? '<' : '',
                #!$trans->three_prime_utr ? '>' : '',
                $mRNA->five_utr,
                $mRNA->three_utr,
                $feature_type,
            );

            $full_string .= $coords_full_qualifiers;

            $full_string .= $ncRNA_qualifier if $ncRNA_qualifier;

            #y print mRNA qualifiers
            $full_string .= $full_qualifiers;

            #r product qualifier
            if ($biotype eq 'tRNA') { 

                #y a load of aegypti transcripts have nothing as opposed to null desc?!?
                my $desc = $mRNA->description() ? $mRNA->description() : q{};
                if ($desc =~ /tRNA-\w{3}$/) { 
                    $full_string .= $qual_kv_pair.'product'.$t.$desc.$n; 
                } else { 
                    $full_string .= $qual_kv_pair.'product'.$t.'tRNA-Xxx'.$n; 
                }
            } else { 
                $full_string .= $rRNA_product;
                $full_string .= $qual_kv_pair.'product'.$t.$mRNA->id().$n; 
            }

            $full_string .= $qual_kv_pair.'db_xref'.$t.'VectorBase:'.$mRNA_id.$n if ($SbmObj && $SbmObj->vectorbase_dbxref());

            #r CDS info if appropriate...
            if ($biotype eq 'pseudogene') { 
                $full_string .= $qual_kv_pair.'pseudo'.$n; 
            } elsif ($biotype eq 'protein_coding') { #b just anal
            #y print CDS section and qualifiers

                my $TRNSL_id = $mRNA->trnsl_name();

                my $splstr = $split_flag ? 'part of ' : q{};
                $full_string .= $qual_kv_pair.'note'.$t.$splstr.$mRNA_id.' coding for '.$TRNSL_id.$n;

                if ($mRNA->external_name()) { 
                    $full_string .= $qual_kv_pair.'note'.$t.$mRNA->external_name.$n; 
                }

                #my $full_cds_seq = $mRNA->seq();

                #y print the actual CDS details...


#                my $exons_to_use_cds = $short_intron_support ? &_fake_big_exon($mRNA) : $mRNA->exons();

                my $coords_full_qualifiers = &_grab_coords(

                  $strand,

                  $short_intron_support ? &_fake_big_exon($mRNA):$mRNA->CDSs(),
#                 $exons_to_use_cds,

                  #$full_cds_seq =~ /^M/ ? '' : '<', 
                  #substr($trans->spliced_seq,$trans->cdna_coding_end-3,3) 
                  #  =~ /(TAG|TAA|TGA)/ || $full_cds_seq =~ /\*$/ ? '' : '>', 
                  $mRNA->trnsl_N_terminal(),
                  $mRNA->trnsl_C_terminal(),
                  'CDS'
                );

                $full_string .= $coords_full_qualifiers;

                $full_string .= $qual_kv_pair.'codon_start'.$t.$mRNA->genbank_codon_offset().$n
                  if $mRNA->genbank_codon_offset(); 

                $full_string .= $full_qualifiers;

                #r/ this needs populating upstream... - i.e. either external id OR translation_stable
                # my $product = $mRNA->product() ? $mRNA->product : 
                # $full_string .= $mRNA->product : $qual_kv_pair.'product'.$t.$mRNA->product.$n;
                $full_string .= $mRNA->product() ? $qual_kv_pair.'product'.$t.$mRNA->product().$n : q{};

                $full_string .= $SbmObj && $SbmObj->vectorbase_dbxref() ? $qual_kv_pair.'db_xref'.$t.'VectorBase:'.$TRNSL_id.$n : q{};

                if ($short_intron_support) {

                    #delete $short_intron_support->{insd};

                    my $pepfile = $SbmObj->_outdir().'/'.$GbAccession.'.pep';
                    #unlink $pepfile if (-e $pepfile);
                    open my $pepfh, '>>', $pepfile or die;
                    print $pepfh '>'.$seqid_accession_valuestring.$short_intron_support->{seq}.qq{\n};
                    close $pepfh or die;

                    #r hmmmmm?!?
                    if ($short_intron_support->{insd}) {
                    # if (exists $short_intron_support->{insd}) {
                        $full_string .= $qual_kv_pair.'exception'.$t
                        .'annotated by transcript or proteomic data'.$n
                        .$qual_kv_pair.'inference'.$t
                        .'similar to AA sequence:'.$short_intron_support->{insd}.$n;
                    } elsif ($SbmObj->heteregeneous_population()) {
                        $full_string .= $qual_kv_pair.'exception'.$t.'heterogenous population sequenced'.$n;
                    } else {
                        $full_string .= $qual_kv_pair.'exception'.$t.'low-quality sequence region'.$n;
                    }
                }

                $full_string .= $qual_kv_pair.'note'.$t.$splstr.$TRNSL_id.' encoded by '.$mRNA_id.$n;
            }
        }

        #$coding_pseudo_genes++;
#        print qq{\n\n};
#        print $tbl_fh $full_string;
        return $full_string;
        #die;
    } else {

        die qq{\nencountered previously unknown e! biotype - need to tell me how to handle it?!?};

    }


    #y things should no be ncRNA in the db - if we have some we need to get a class for them
#    elsif ($biotype eq 'ncRNA') { 
#        push @skipped_ncRNA_list, $gene->stable_id; 
#    } else { 
#        push @skipped_list, [$gene->stable_id, $gene->biotype];
#    } 
}

sub _print_coords {

    #y this is not elegant, but...
    my ($strand,$list_ref,$fiveprimeNterminal,$threeprimeCterminal,$type,$tbl_fh) = @_;

    my $n = $Tools::n;
    my $t = $Tools::t;
    my $qual_kv_pair = $Tools::qual_kv_pair;

    #y screw it, no need for this...
    #my $fiveprimeNterminal = &start_coderef->()...

    my @list = @{$list_ref};

    #y single exons genes 
    if (scalar @list == 1) {
        print $tbl_fh $strand == 1 ? 
          $fiveprimeNterminal.$list[0]->start.$t.$threeprimeCterminal.$list[0]->end.$t.$type.$n : 
          $fiveprimeNterminal.$list[0]->end.$t.$threeprimeCterminal.$list[0]->start.$t.$type.$n;
          #print $tbl_fh qq{\n----\n}; # freakish logic test is after this?!?
    }
    else {
        for my $i (0..$#list) {
        # for my $i (0..$len-1) {
        #use Data::TreeDraw; draw($list[$i]); 
            if ($i == 0) {
                  #print $fiveprimeNterminal.$list[$i]->start.$t.$list[0]->end.$t.$type.$n;
                print $tbl_fh $strand == 1 ? 
                  $fiveprimeNterminal.$list[$i]->start.$t.$list[0]->end.$t.$type.$n : 
                  $fiveprimeNterminal.$list[$i]->end.$t.$list[0]->start.$t.$type.$n;
            }
            elsif ($i == $#list) {
                #  print $list[$i]->start.$t.$threeprimeCterminal.$list[$#list]->end.$t.$n;
                print $tbl_fh $strand == 1 ? 
                  $list[$i]->start.$t.$threeprimeCterminal.$list[$#list]->end.$t.$n : 
                  $list[$i]->end.$t.$threeprimeCterminal.$list[$#list]->start.$t.$n;
            }
            else {
                #  print $list[$i]->start.$t.$list[$i]->end.$t.$n;
                print $tbl_fh $strand == 1 ? 
                  $list[$i]->start.$t.$list[$i]->end.$t.$n : 
                  $list[$i]->end.$t.$list[$i]->start.$t.$n;
            }
        }
    }

    #y don't forget to return!?!
    return;
}

sub _grab_coords {

    #y this is not elegant, but...
    my ($strand,$list_ref,$fiveprimeNterminal,$threeprimeCterminal,$type) = @_;

    my $n = $Tools::n;
    my $t = $Tools::t;
    my $qual_kv_pair = $Tools::qual_kv_pair;

    #y screw it, no need for this...
    #my $fiveprimeNterminal = &start_coderef->()...

    my @list = @{$list_ref};

    my $coords_string = q{};

    #y single exons genes 
    if (scalar @list == 1) {
        $coords_string .= $strand eq '+' ? 
          $fiveprimeNterminal.$list[0]->start.$t.$threeprimeCterminal.$list[0]->end.$t.$type.$n : 
          $fiveprimeNterminal.$list[0]->end.$t.$threeprimeCterminal.$list[0]->start.$t.$type.$n;
    }
    else {
        for my $i (0..$#list) {
        # for my $i (0..$len-1) {
        #use Data::TreeDraw; draw($list[$i]); 
            if ($i == 0) {
                  #print $fiveprimeNterminal.$list[$i]->start.$t.$list[0]->end.$t.$type.$n;
                $coords_string .= $strand eq '+' ? 
                # $coords_string .= $strand == 1 ? 
                  $fiveprimeNterminal.$list[$i]->start.$t.$list[0]->end.$t.$type.$n : 
                  $fiveprimeNterminal.$list[$i]->end.$t.$list[0]->start.$t.$type.$n;
            }
            elsif ($i == $#list) {
                #  print $list[$i]->start.$t.$threeprimeCterminal.$list[$#list]->end.$t.$n;
                $coords_string .= $strand eq '+' ? 
                #$coords_string .= $strand == 1 ? 
                  $list[$i]->start.$t.$threeprimeCterminal.$list[$#list]->end.$t.$n : 
                  $list[$i]->end.$t.$threeprimeCterminal.$list[$#list]->start.$t.$n;
            }
            else {
                #  print $list[$i]->start.$t.$list[$i]->end.$t.$n;
                $coords_string .= $strand eq '+' ? 
                #$coords_string .= $strand == 1 ? 
                  $list[$i]->start.$t.$list[$i]->end.$t.$n : 
                  $list[$i]->end.$t.$list[$i]->start.$t.$n;
            }
        }
    }

    #y don't forget to return!?!
    return $coords_string;
}

sub _protein_coding_string {
#y create the string to print for each translation - poorly named - invoked for all CDS

    #/ split not yet integrated...
    my ($mRNA,$SbmObj,$split_flag) = @_;
    $split_flag||=0;

    my $TRNSL_id = $mRNA->trnsl_name();
    my $nucAccession = $mRNA->seqid();
    my $n = $Tools::n;
    my $t = $Tools::t;
    my $qual_kv_pair = $Tools::qual_kv_pair;

    my $peptideSeqID;
    my $old_peptide_suffix = q{};
        
    my $proj_prefix = $SbmObj ? $SbmObj->wgs_proj_prefix() : 'ENSEMBLGENOMES';

    my $TranslationStables2PeptideAccessionsNseqIDs = $SbmObj ? $SbmObj->_TranslationStables2PeptideAccessionsNseqIDs() : +{};

    #y if its an old peptide grab old_peptide suffix, accession, SeqID

    print qq{\nDEBUG: undefined value for }.$mRNA->id() if (!defined $TRNSL_id || !$TRNSL_id);
    
    if (exists $TranslationStables2PeptideAccessionsNseqIDs->{$nucAccession.'_'.$TRNSL_id}) {

        my $peptideAccession = $TranslationStables2PeptideAccessionsNseqIDs->{$nucAccession.'_'.$TRNSL_id}->{protAccession}; 
        $peptideSeqID = $TranslationStables2PeptideAccessionsNseqIDs->{$nucAccession.'_'.$TRNSL_id}->{protSeqId}; 
        $old_peptide_suffix = '|gb|'.$peptideAccession;
        #$old_proteins++;
    }
    #y/ for new peptides only use compound names as seqids if they are split between contigs
    else {
        if ($split_flag) {
            $peptideSeqID  = $nucAccession.'_'.$TRNSL_id;
        } else {
            $peptideSeqID  = $TRNSL_id;
        }
        #$new_proteins++;
    }

    print qq{\nDEBUG: undefined value for }.$mRNA->id() if (!defined $peptideSeqID || !$peptideSeqID);

    #r transcript_id - just usual format that is cross-referenced to transl
    my $string = $qual_kv_pair.'transcript_id'.$t.'gnl|WGS:'.$proj_prefix.'|'.$mRNA->id().$n;
    #r if a resubmitted protein must have same details - in suffix
    $string .= $qual_kv_pair.'protein_id'.$t.'gnl|WGS:'.$proj_prefix.'|'.$peptideSeqID.$old_peptide_suffix.$n;

    my $seqid_accession_valuestring = 'gnl|WGS:'.$proj_prefix.'|'.$peptideSeqID.$old_peptide_suffix.$n;

    print qq{\nDEBUG: undefined value for }.$mRNA->id() if (!defined $TRNSL_id || !defined $peptideSeqID);

    return ($string,$seqid_accession_valuestring);
}

sub _silly_vb_description_mess {
    my $translation = shift;
    my $product = undef;
    #my $product = q{};
    for my $x (@{$translation->get_all_DBEntries}) {
        $product = $x->primary_id if ($x->dbname eq 'VB_Community_Symbol');
    }
    #/ duh - it's set to empty string so it ALWAYS defined!?! - either use simple if or set to undef!?!
    return defined $product ? $product : $translation->stable_id;
}

#/ atm. strictly for coding genes - easily extended...
sub _partial_from_split_eGene {

    my ($class, $split_gene, $submission_cs, $seqid_with_part_of_gene) = @_;

    #print qq{\nfuck!?! here: $submission_cs / $seqid_with_part_of_gene}; 
    #r Gb submission so has nothing to do with the convinience of canonical. we need to build mini transcripts 
    #r for each transcript corresponding to parts on each seqid_with_part_of_gene - ought to do all in one go but too lazy...
    #r so we go through each transcript...
    #b instead of printing each transcript at a time we just build a GffDoc gene per transcript...

    #y project from default cs (i.e. should be toplevel) to submission level
    my @list_of_projsegs = @{$split_gene->project($submission_cs)};

    #print qq{\ngggg: }.$_->to_Slice()->seq_region_name() for (@list_of_projsegs);


    #y grab terminal projection segments
    my $five_prime_end_seqid_with_part_of_gene = $list_of_projsegs[0]->to_Slice->seq_region_name;
    my $three_prime_end_seqid_with_part_of_gene = $list_of_projsegs[$#list_of_projsegs]->to_Slice->seq_region_name;
    #my $three_prime_end_seqid_with_part_of_gene = $list_of_projsegs[$#list_of_projsegs]->[2]->seq_region_name;

    my $strand = $split_gene->strand;

    my $print_gene_flag = 1;
    my $g_string = q{};

    #/ we allow the info to be incorrect wrt gene/mrna start/end etc as this is solely for Gb which uses
    #/ exons/CDS only

    #r/ need to re-write gene start and end!?!

    my $gffGene = bless $split_gene->toGeneDocString, 'GffDoc::Gff3::Gene';
    $gffGene->_split(1);

    #//// as this is a split gene we've retrived from a higher cs - so next to translate start/end
    #//// but we don't yet know what that is - i.e. need to use maximal extent of mRNAs on this seqid
    #//// exons are always returns 5'-3' hence strand dependency of code - thus just ignore and check at end?!?

    $gffGene->_seqid_for_partial($seqid_with_part_of_gene);
    #my $gffGene = bless {}, 'GffDoc::Gff3::Gene';
    #$gffGene->_start(undef);

    $gffGene->_external_name($split_gene->external_name()) if ($split_gene->external_name());
    $gffGene->_description($split_gene->description()) if ($split_gene->description());

    my $v_first_exon = 1;
    my $no_exons = 1;

    ######################################
    #b get the local - i.e. current transcript extent
    my $local_five_gene;
    my $local_three_gene;

    #///// move these out of transcript loop so they get set just the once!?!
    #b/ why are there two flags for same thing - i.e. just opposites!?!
    my $first_exon_on_this_contig_flag = 1;
    my $past_first_exon_on_this_contig_flag = 0;
    ######################################

    TR:
    #/ this is not designed for gff/gtf printing - we aren't bothering with phase info... - as need offset...
    for my $tr (@{$split_gene->get_all_Transcripts}) {

        my $gffmRNA = bless $tr->toGeneDocString, 'GffDoc::Gff3::mRNA';
        $gffmRNA->_seqid_for_partial($seqid_with_part_of_gene);

        #my $gffmRNA = bless {}, 'GffDoc::Gff3::mRNA';
        $gffmRNA->_external_name($tr->external_name()) if ($tr->external_name());
        $gffmRNA->_description($tr->description()) if ($tr->description());

        my @gffCDSs;
        my @gffexons;
        
        my $five_prime_exon_of_transcript = 0;
        my $three_prime_exon_of_transcript = 0;

        #r exons ergo the slices will be returned 5'->3' however, the slices will return coords independent of strand...
        my @exons = @{$tr->get_all_Exons};

        #r/ why fake exons when you could just use the originals?!?
        ############## go through exons - essential build exons for this transcript ###############
        for my $p (0..$#exons) {

            #y project exon down to submission level
            my $pp = $exons[$p]->project($submission_cs);


#            if (scalar @{$pp} >> 1) {
#                print qq{\nTHIS GENE HAS A SPLIT EXON!?! IT IS A MESS!?! - skipping\n};
#                print ref $pp->[0]; 
#                print @{$pp}; 
#                exit;
#            }
            
            
            #/ exons should only have 1 projection segment hence [0]
#            print qq{\nprojction part: }.$pp->[0]->to_Slice()->seq_region_name;
            #print qq{\nprojction part: }.$pp->[0][2]->seq_region_name;
#            print qq{\nwhat we want: }.$seqid_with_part_of_gene;
            #y ignore exon if not on the particular contig we're processing
            next if ($pp->[0]->to_Slice()->seq_region_name ne $seqid_with_part_of_gene);
            #next if ($pp->[0][2]->seq_region_name ne $seqid_with_part_of_gene);


$no_exons = 0;

            #y why is this down here?!?
            $five_prime_exon_of_transcript = 1 if ($p == 0);
            $three_prime_exon_of_transcript = 1 if ($p == $#exons);

    #b/ this is important - i.e. if this transcript does not include the 5' terminal the CDSs may be out of phase
    #b/ in such cases without stating the out of phaseness of the transcript Gb will translate it out of phase!?!
    #r/ we add this phase info to stop it breaking such transcripts - phase 0 (=codon_start 1) is irrelevant as its in phase
    #r/ this must be added as a flag?!? e.g. split_cds_transcript_start_phase method!?! - i.e. print if set...

            #r translatable returns truncated exons - i.e. cds so they will have e! == 0
            #y default codon_start is 1 so if e! phase is 0 or -1 just leave
            if ($first_exon_on_this_contig_flag && $exons[$p]->phase == 2) { #r not sure why 
                $gffmRNA->_genbank_codon_offset(2); 
            } elsif ($first_exon_on_this_contig_flag && $exons[$p]->phase == 1) { 
                $gffmRNA->_genbank_codon_offset(3); 
            }

            #r/ the two first exon flags are redundant!?!
            undef $first_exon_on_this_contig_flag;
            
            my $gffexon = bless {}, 'GffDoc::Feature::Gff3::exon';

    #b/ since we have projected the exon down to submission level to find out where it is we not
    #b/ get the values from the projected slice - we ASSUME that the exon is projected in its 
    #b/ entirity to new contig - i.e. it is NOT split - this breaks a gene in pest where the exon 
    #b/ actually spans a contig gap - as there is no gap between them!?!
            
            if ($strand == 1) { 
                $gffexon->_start($pp->[0][2]->start);
                $gffexon->_end($pp->[0][2]->end);
                $local_five_gene = $pp->[0][2]->end if (!$past_first_exon_on_this_contig_flag); 
                $past_first_exon_on_this_contig_flag = 1;
                $local_three_gene = $pp->[0][2]->start;
            }
            else {
                $gffexon->_start($pp->[0][2]->end);
                $gffexon->_end($pp->[0][2]->start);
                $local_three_gene = $pp->[0][2]->start if (!$past_first_exon_on_this_contig_flag); 
                $past_first_exon_on_this_contig_flag = 1;
                $local_five_gene = $pp->[0][2]->end;
            }

            my $mrna_end;
            my $cds_end;

            #y for now only giving a shit about 3' end as it complains...
            if ($three_prime_exon_of_transcript) {

                #y use standard rules...
                $gffGene->_three_utr(!$tr->three_prime_utr ? '' : '>');
                $gffmRNA->_three_utr(!$tr->three_prime_utr ? '' : '>');

                $gffmRNA->_trnsl_C_terminal(substr($tr->spliced_seq,$tr->cdna_coding_end-3,3) 
                  =~ /(TAG|TAA|TGA)/ ? '' : '>');
                #$gffmRNA->_trnsl_C_terminal($tr->translation()->seq() =~ /\*$/ ? '' : '>');

            }
            else {
                #y just not really very bothered...
                $gffGene->_three_utr('>');

                $gffmRNA->_three_utr('>');
                $gffmRNA->_trnsl_C_terminal('>');
            }

            my $mrna_starting;
            my $cds_starting;

            if ($five_prime_exon_of_transcript) {

                #y just not really very bothered...
                $gffGene->_five_utr(!$tr->five_prime_utr ? '<' : '');

                #y use standard rules...
                $gffmRNA->_five_utr(!$tr->five_prime_utr ? '<' : '');
                # !$tr->five_prime_utr ? '' : '<';
                $gffmRNA->_trnsl_N_terminal($tr->spliced_seq =~ /^M/ ? '' : '<');
                # $tr->spliced_seq =~ /^M/ ? '' : '<';
            }
            else {

                #y just not really very bothered...
                $gffGene->_five_utr(!$tr->five_prime_utr ? '<' : '');

                $gffmRNA->_five_utr('<'); 
                # $mrna_starting = '<';
                $gffmRNA->_trnsl_N_terminal('<');
                # $cds_starting = '<';
            }

            push @gffexons, $gffexon;
                    
        }

        ############## go through CDS - build CDS for this transcript #################

        for my $p (@{$tr->get_all_translateable_Exons}) {

            my $pp = $p->project($submission_cs);

            my $gffCDS = bless {}, 'GffDoc::Feature::Gff3::CDS';

            if ($strand == 1) { 
                $gffCDS->_start($pp->[0][2]->start);
                $gffCDS->_end($pp->[0][2]->end);
            }
            else {
                $gffCDS->_start($pp->[0][2]->end); 
                $gffCDS->_end($pp->[0][2]->start);
            }

            push @gffCDSs, $gffCDS if ($pp->[0][2]->seq_region_name eq $seqid_with_part_of_gene);

        }

        $gffmRNA->_exons(\@gffexons);
        $gffmRNA->_CDSs(\@gffCDSs);
        $gffmRNA->_gffGene($gffGene);
        $gffGene->addmRNA($gffmRNA);

        #use Data::Dumper; print Dumper $gffGene; die;


        #r/ really nasty but i'm tired and bored so find gene extent here!?!
#        my $gene_start;
#        my $gene_end;
#        my ($l,$h);
#        my $first_exon = 1;
#        print qq{\nentering loop for gene!};

#       my $no_exons = 1;

#        for my $t (@{$gffGene->mRNAs()}) {
#            print qq{\ntranscript loop};
#            for my $e (@{$t->exons()}) {
#                print qq{\nexon loop};
#                print qq{\nexon: }.$e->start.'-'.$e->end.' - '.$gffGene->strand;
#                $no_exons = 0;
#            }
#        }
        return 0 if $no_exons;

=head1

#                if ($first_exon) {
#                    if ($gffGene->strand() eq '+') { 
#                        $gene_start = $e->start();
#                        $gene_end = $e->end();
#                    } else {
#                        $gene_start = $e->start();
#                        $gene_end = $e->end();
#                        $l = $e->start() >= $e->end() ? $e->end() : $e->start();
#                        $h = $e->start() <= $e->end() ? $e->end() : $e->start();
#                    }
#                        $first_exon = 0;
#                } else {
#                    if ($gffGene->strand() eq '+') { 
#                        $gene_start = $e->start() if ($gene_start > $e->start());
#                        $gene_end = $e->end() if ($gene_end < $e->end());
#                    } else {
#                        $gene_start = $e->end() if ($gene_start < $e->end());  
                        #$gene_start = $e->end() if ($gene_start < $e->end());   xxx - 4008
                        #$gene_start = $e->end() if ($gene_start > $e->end()); xxx - 721
#                        $gene_end = $e->start() if ($gene_end > $e->start());  
                        #$gene_end = $e->start() if ($gene_end < $e->start());  4099  -
                        #$gene_end = $e->start() if ($gene_end > $e->start());  453  -
#                        $l = $e->end() if ($l > $e->end());  
#                        $l = $e->start() if ($l > $e->start());  
#                        $h = $e->end() if ($l < $e->end());  
#                        $h = $e->start() if ($l < $e->start());  
#                    }
#                }
#            }
#        }

        #/ ignore any gene that has no exons on the contig!?! 
        return 0 if $no_exons;



        #print qq{\nwe have high and low: }.$l.' / '.$h;



        #$gffGene->_start($gene_start);
        # $gffGene->_end($gene_end);
        # $gffGene->_start($l);
        #$gffGene->_end($h);
        #print qq{\nhigh: $h\nlow: $l};
#        $gffGene->_start($h);
#        $gffGene->_end($l);

=cut

        #my $transl = 
        #  ($trans->translation && ($trans->biotype() eq 'protein_coding')) ? $trans->translation : 0; 
        #if ($transl) {
        $gffmRNA->_trnsl_name($tr->translation()->stable_id());
        #$gffmRNA->_trnsl_id($trans->translation()->stable_id());
        $gffmRNA->_product(&_silly_vb_description_mess($tr->translation()));
        #}

    }

    $gffGene->_start($local_five_gene);
    $gffGene->_end($local_three_gene);

    return $gffGene;
}

sub from_eGene {

    my ($class, $eGene) = @_;

    #y are we really going to be accepting many args - i.e. can remove those two and force rest 
    #y to hash for named args. or just have a boolean on third position?!?

#    my ($class, $eGene, $remove_short_introns?!?) = splice (@_, 0, 3);

    #while (my ($k,$v) = each %hash) { print qq{\nk: $k => v: $v }; }

    die if (!$eGene->isa('Bio::EnsEMBL::Gene'));

    #y circumvent the constructor - i.e. really no need...
    #my $gffGene = GffDoc::Gff3::Gene->new(%{$eGene->toGeneDocString()});

    #r/ check for royally fucked genes!?!
    my $gffGene;
    eval {

    $gffGene = bless $eGene->toGeneDocString, 'GffDoc::Gff3::Gene';
    # my $gffGene = bless $eGene->toGeneDocString, 'GffDoc::Gff3::Gene';

#/debug!?! 
#next if ($gffGene->id ne 'AGAP001127');

                #!$trans->five_prime_utr ? '<' : '',
                #!$trans->three_prime_utr ? '>' : '',
    $gffGene->_five_utr(!$eGene->canonical_transcript->five_prime_utr ? '<' : '');
    $gffGene->_three_utr(!$eGene->canonical_transcript->three_prime_utr ? '>' : '');
    $gffGene->_external_name($eGene->external_name()) if ($eGene->external_name());
    $gffGene->_description($eGene->description()) if ($eGene->description());
    };
    if ($@) {
        print qq{THERE'S AN ISSUE WITH GENERATING gffGene from }.$eGene->stable_id();
        next;
    }

    #/ clearly instead we could use if(biotype eq 'protein_coding' and then use get_all_translateable_Exons separately...
    #/ but then we're also stuck with incorrect phase info if its in the db etc. - i.e. that method uses fact that some 
    #/ values are stored in the db...
    
    my $strand = $eGene->strand == 1 ? '+' : '-'; 

    for my $trans (@{$eGene->get_all_Transcripts}) {

        #r don't bother with constructor - we already know these bits fit together so can be much sloppier
        my $gffmRNA = bless $trans->toGeneDocString, 'GffDoc::Gff3::mRNA';

        $gffmRNA->_external_name($trans->external_name()) if ($trans->external_name());
        $gffmRNA->_description($trans->description()) if ($trans->description());

        #r damned lazy...
        $gffmRNA->_five_utr(!$trans->five_prime_utr ? '<' : '');
        $gffmRNA->_three_utr(!$trans->three_prime_utr ? '>' : '');

        my @exons;
        my @CDS;

        my $cds_count = 0;
        my $phase_minusOne;   
        my $length_minusOne;
        my $ending = 0;

        #y define scalars to store important anchors
        my $first_transl_exon;
        my $last_transl_exon; 
        my $transl_start; 
        my $transl_end;
        my $physical_start;
        my $physical_end;
                
        #y not actually checking biotype checking for translation presence - prolly ought to make it both!?!
        #y set translation flag - i.e. it has CDS processing and get translation object
        my $transl = ($trans->translation && ($trans->biotype() eq 'protein_coding')) ? $trans->translation : 0; 
        #my $transl = $trans->translation ? $trans->translation : 0; 

        #y/ prolly ought to refuse biotype protein_coding without translation

        if ($transl) {
            $physical_start = $strand eq '+' ? $trans->start : $trans->end;
            $physical_end = $strand eq '+' ? $trans->end : $trans->start;
            $first_transl_exon  = $transl->start_Exon;
            $last_transl_exon    = $transl->end_Exon;
            $transl_start     = $transl->start;
            $transl_end       = $transl->end;
                  
            #/ if check?!? 
            $gffmRNA->_trnsl_id($trans->translation()->stable_id());
            
            #print qq{\nhey: }.$trans->seq()->seq().qq{\n};
            #print qq{\nhey: }.$trans->translation()->seq().qq{\n};

            #/ can use ^M with peptide of ATG with cDNA...
            $gffmRNA->_trnsl_N_terminal($trans->translation()->seq =~ /^M/ ? '' : '<');

            #/ api cleaves of terminal stops so do at cDNA level?!?
            $gffmRNA->_trnsl_C_terminal(substr($trans->spliced_seq,$trans->cdna_coding_end-3,3) 
             =~ /(TAG|TAA|TGA)/ ? '' : '>');
            #$gffmRNA->_trnsl_C_terminal($trans->translation()->seq() =~ /\*$/ ? '' : '>');

            #y sheer laziness...
            $gffmRNA->_trnsl_name($trans->translation()->stable_id());
    
            $gffmRNA->_product(&_silly_vb_description_mess($trans->translation()));
        }

        # my ($exon_list = $fake_short_introns ? &_fake_big_exon($trans) : $trans->get_all_Exons();

        for my $exon (@{$trans->get_all_Exons})  {

            my @exon_position = ($exon->start, $exon->end);

            push @exons, bless $exon->toGeneDocString, 'GffDoc::Feature::Gff3::exon';

######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!
            if ($transl && !$ending) {
#            if ($transl) {
######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!

                my  $alt_strand = $eGene->strand == 1 ? '-' : '+';
                my  $switch = $eGene->strand == 1 ? 0 : 1 ;

                next if ($exon ne $first_transl_exon && !@CDS);

######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!
                last if ($ending);
######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!
######## fuckwit you are putting CDS in that our non-coding!?!


                #y should put into single conditional
                #my $five_utr = q{};
                #$five_utr = $exon_position[$switch] if ($exon eq $first_transl_exon);

                my $l = $exon->length;

                $exon_position[$switch] = eval(qq{$exon_position[$switch]$strand$transl_start${alt_strand}1}) 
                  if ($exon eq $first_transl_exon);

                #if ($exon eq $first_transl_exon) {
                #    $five_utr =  $five_utr == $exon_position[$switch] ? '<' : q{};
                #}

                $exon_position[1-$switch] = eval(qq{$exon_position[1-$switch]$strand($transl_end-$l)}) 
                  if ($exon eq $last_transl_exon);

                $l = $exon_position[1] - $exon_position[0] + 1 
                  if ($exon eq $first_transl_exon || $exon eq $last_transl_exon);

                #f these are ORDERED EXONS!?! so its always the same irrespective of strand
                my $p = $cds_count == 0 ? 0 : ((3 -(($length_minusOne%3)-$phase_minusOne) ) % 3);

                push @CDS, bless +{ %{$exon->toGeneDocString},
                    #y these values over-write
                    start => $exon_position[0], end => $exon_position[1], phase => $p }, 'GffDoc::Feature::Gff3::CDS';

                #y stash values for phase calculations for next iteration
                $length_minusOne = $l;
                $phase_minusOne = $p;
                $cds_count++;

                $ending = 1 if ($exon eq $last_transl_exon);
            } 
        }

        #use Data::Dumper; print Dumper \@CDS;  die;
        $gffmRNA->_exons(\@exons);
        $gffmRNA->_CDSs(\@CDS);
        $gffmRNA->_gffGene($gffGene);
        $gffGene->addmRNA($gffmRNA);

        #y using constructor so have to give gene at same time...
        #my $gffmRNA = GffDoc::Gff3::mRNA->new(%{$trans->toGeneDocString()}, gffGene => $gffGene);
        #my $gffGene = bless $eGene->toGeneDocString, 'GffDoc::Gff3::Gene';

    }

    return $gffGene;

}

sub _fake_big_exon {

    my ($mRNA,$cds) = @_;
    #print qq{\n>DEBUG FAKING A BIG EXON for }.$mRNA->id();
    $cds||=0;

    my $exon_list = $cds ? $mRNA->CDSs() : $mRNA->exons();
    my $strand = $mRNA->strand();
    my @exons;
    my $c = 0;
    my $start_exon_in_fake;
    my $end_exon_in_fake;
    my $state = 0;

    my $grr = 0;

    my $last_exon;
    my $last_end = undef;

    for my $current_exon (@{$exon_list}) {
    # for my $current_exon (@{$trans->exons()}) {
    
    # for my $current_exon (@{$t->get_all_Exons()}) {

        if (!defined $last_end) {
            $last_end = $current_exon->end();
        } else {

            my $intron_length = $strand eq '+' ? $current_exon->start()-$last_exon->end() : $last_exon->start()-$current_exon->end();
            #my $intron_length = $strand == 1 ? $current_exon->start()-$last_exon->end() : $last_exon->start()-$current_exon->end();
            $intron_length--;

#            print qq{\nnew intron length: $intron_length\n};
            if ($intron_length < 10) {
#                print qq{\nSMALL [}.++$grr.']';
            }
            if ($intron_length < 10 && !$state) { 
                $state = 1;
                $start_exon_in_fake = $last_exon;
#                print qq{\nstart: $start_exon_in_fake};
                $end_exon_in_fake = $current_exon;
#                print qq{\nend: $end_exon_in_fake};
            #y over-write end any time we see another small...
            } elsif ($intron_length < 10) {
                $end_exon_in_fake = $current_exon;
#                print qq{\nend: $end_exon_in_fake};
            }
        }
        $last_exon = $current_exon;
    }

    #r so not elegant... but now have the start and end of fake exon so need strand dependency and then replace it in list

    #y/ fake it
    
    my $fakeexon = $cds ? bless { id => 'FakeExon', phase => '.' }, 'GffDoc::Feature::Gff3::CDS' : bless { id => 'FakeExon', phase => '.' }, 'GffDoc::Feature::Gff3::exon';
    
    # my $fakeexon = bless {}, 'GffDoc::Gff3::Feature::exon';
    # my $fakeexon = FakeExon->new();
    if ($strand eq '+') {
    # if ($strand == 1) {
        $fakeexon->_start($start_exon_in_fake->start());
        $fakeexon->_end($end_exon_in_fake->end());
    } else {
        $fakeexon->_end($start_exon_in_fake->end());
        $fakeexon->_start($end_exon_in_fake->start());
    }
    #print qq{\n\n >>> fake exon: }, $fakeexon->start(), '---', $fakeexon->end();
    my $before = 1;
    my $after = 0;
    my $replaced = 0;

    #y/ replace them
    my $e = 0;
    for my $exon (@{$exon_list}) {
        $e++;
#        print qq{\n > [$e] EXON: }.$exon->start().'-'.$exon->end();
        if ($exon == $start_exon_in_fake) {
            $before = 0;
#            print qq{\n# switching before off: before[$before] after[$after]};
#            print qq{\n  > PUTTING IN FAKE: }.$fakeexon->start().'-'.$fakeexon->end();
            push @exons, $fakeexon;
        }
        elsif ($exon == $end_exon_in_fake) {
#            print qq{\n# switching after on: before[$before] after[$after]};
            $after = 1;
        }
        #/ duh we want before OR after!?!
        #elsif ($before || !$after) {
        elsif ($before || $after) {
#            print qq{\n# before is true or after is false: before[$before] after[$after]};
#            print qq{\n  > shoving standard exon in}; 
            push @exons, $exon;
        }
        #y we replace the problem region...
    }
    #use Data::Dumper; print Dumper \@exons; 
    return \@exons;
}

sub toGff3String {

    my ($self, $short_intron) = @_;

#b/ small intron handling.    
#    $short_intron = 1;

    my $t = qq{\t};
    my $n = qq{\n};
    my $seqid = $self->seqid();
    my $source = $self->source() ? $self->source() : 'ensemblgenomes';
    my $strand = $self->strand();
    my $biotype = $self->biotype() eq 'protein_coding' ? 'mRNA' : $self->biotype();

    my $str = $seqid.$t.$source.$t.'gene'.$t.$self->start().$t.$self->end()
      .$t.q{.}.$t.$strand.$t.q{.}.$t.'ID='.$self->id().';'.$n;

    for my $mRNA (@{$self->mRNAs}) {

        $str .= $seqid.$t.$source.$t.$biotype.$t.$mRNA->start().$t.$mRNA->end()
          .$t.q{.}.$t.$strand.$t.q{.}.$t.'ID='.$mRNA->id().';Parent='.$self->id().';'.$n;

#b/ fake merge exons with small introns...
#        my $exons = &_fake_big_exon($mRNA) if $short_intron;
#        for my $e (@{$exons}) {
        for my $e (@{$mRNA->exons()}) {

            $str .= $seqid.$t.$source.$t.'exon'.$t.$e->start().$t.$e->end()
              .$t.q{.}.$t.$strand.$t.q{.}.$t.'ID='.$e->id().';Parent='.$mRNA->id().';'.$n;
        }    

#b/ fake merge exons with small introns...
#        my $CDSs = &_fake_big_exon($mRNA,1) if $short_intron;
#        for my $c (@{$CDSs}) {
        for my $c (@{$mRNA->CDSs()}) {

            $str .= $seqid.$t.$source.$t.'CDS'.$t.$c->start().$t.$c->end()
              .$t.q{.}.$t.$strand.$t.$c->phase().$t.'ID='.$c->id().';Parent='.$mRNA->id().';'.$n;
        }    
    }

    return $str;
}

sub toGtfString {

    my ($self,$exons) = splice(@_,0,2);
    $exons||= 0;
    my $t = qq{\t};
    my $n = qq{\n};
    my $seqid = $self->seqid();
    my $source = $self->source() ? $self->source() : 'ensemblgenomes';
    my $strand = $self->strand();

    my $str = q{};

    for my $mRNA (@{$self->mRNAs}) {

        if ($exons) {

            for my $e (@{$mRNA->exons()}) {

                $str .= $seqid.$t.$source.$t.'exon'.$t.$e->start().$t.$e->end()
                  .$t.q{.}.$t.$strand.$t.q{.}.$t.'gene_id "'.$self->id().'"; transcript_id "'.$mRNA->id().'";'.$n;
            }     

        }

        for my $c (@{$mRNA->CDSs()}) {

            $str .= $seqid.$t.$source.$t.'CDS'.$t.$c->start().$t.$c->end()
              .$t.q{.}.$t.$strand.$t.$c->phase().$t.'gene_id "'.$self->id().'"; transcript_id "'.$mRNA->id().'";'.$n;
        }    
    }

    return $str;
}

around BUILDARGS => 
sub {
    my $orig  = shift;
    my $class = shift;
    my %hash = @_;

    if (exists $hash{eGene}) {
        print qq{\n > Not as yet implemented!\n};
        die;
    }

    return $class->$orig(%hash);
};

sub addmRNA { 
    my $self = shift;
    Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::WrongType', error => 'i want a GffDoc::Gff3::mRNA object' ) 
      if (!$_[0]->isa('GffDoc::Gff3::mRNA'));
    push @{$self->{mRNAs}}, @_;
}

__PACKAGE__->meta->make_immutable();

1;

=head3 GffDoc::Gff3::mRNA

    Composed class whose reason for existence is to have mRNA feature objects as different type to those that are part of 
    built gene objects. with GffDoc::Feature::Role.

    Instances of this class also hold the info for building ensembl-style translation objects (i.e. relative to transcript
    rather than as a separate fully-fledged feature as in gff). Consequently, it also provides the methods for doing this
    conversion: orderNPhaseCalc, eTrnslConvert.

=cut

package GffDoc::Gff3::mRNA;
use Mouse;
use Scalar::Util 'refaddr';
with 'GffDoc::Feature::Role';
#has 'parent' => (is => 'ro', isa => 'Str', writer => '_parent');  #y not relevant anymore
#y to map up to the gene info - thus need consistency check as use gene's slice, strand, biotype...
has 'gffGene' => (is => 'ro', isa => 'GffDoc::Gff3::Gene', required => 1, writer => '_gffGene');
has 'exons' => (is => 'ro', isa => 'ArrayRef[GffDoc::Feature::Gff3::exon]', required => 1, writer => '_exons');
has 'CDSs' => (is => 'ro', isa => 'ArrayRef[GffDoc::Feature::Gff3::CDS]', required => 1, writer => '_CDSs');
#y either record the appropriate exon info here - i.e. by number or reference or store it in the exon itself 
#y makes no diff - i.e. either boolean: exon->is_start() vs. address exon eq mRNA->start_exon - use Scalar::Util 'refaddr';
has 'trnsl_start_exon' => (is => 'ro', isa => 'GffDoc::Feature::Gff3::exon', writer => '_trnsl_start_exon');
has 'trnsl_start_base' => (is => 'ro', isa => 'Int', writer => '_trnsl_start_base');
has 'trnsl_end_exon' => (is => 'ro', isa => 'GffDoc::Feature::Gff3::exon', writer => '_trnsl_end_exon');
has 'trnsl_end_base' => (is => 'ro', isa => 'Int', writer => '_trnsl_end_base');
has 'trans_start_exon' => (is => 'ro', isa => 'GffDoc::Feature::Gff3::exon', writer => '_trans_start_exon');
has 'trans_end_exon' => (is => 'ro', isa => 'GffDoc::Feature::Gff3::exon', writer => '_trans_end_exon');
has 'trnsl_N_terminal' => (is => 'ro', default => q{}, writer => '_trnsl_N_terminal');
has 'trnsl_C_terminal' => (is => 'ro', default => q{}, writer => '_trnsl_C_terminal');
has 'trnsl_id' => (is => 'ro', default => q{}, writer => '_trnsl_id');
has 'five_utr' => (is => 'ro', default => q{}, writer => '_five_utr');
has 'three_utr' => (is => 'ro', default => q{}, writer => '_three_utr');
has 'genbank_codon_offset' => (is => 'ro', default => q{}, writer => '_genbank_codon_offset');

#b/ only for seqedit addition on print out and .pep for fake_exons
has 'trnsl_name' => (is => 'ro', default => q{}, writer => '_trnsl_name');
has 'external_name' => (is => 'ro', default => q{}, writer => '_external_name');
has 'description' => (is => 'ro', default => q{}, writer => '_description');
has 'product' => (is => 'ro', default => q{}, writer => '_product');

sub orderNPhaseCalc {

    #/ since phase is so unreliable it just ignores them anyway

    my $self = shift;

    my @exons  = sort { $a->start() <=> $b->start() } @{$self->exons};
    my @CDSs  = sort { $a->start() <=> $b->start() } @{$self->CDSs};
           
    if ($self->gffGene->strand()  eq '-') { 
        @exons = reverse @exons;
        @CDSs  = reverse @CDSs;
    }

    my $phase_minusOne;
    my $length_minusOne;
    my @temp;

    CDSPHASE:
    for my $i (0..$#CDSs) {

        #o ensembl phase = (length+ensphase)%3 = (length-gffphase)%3
        #o i.e. ens-phase = (3-gffphase)%3

        my $phase = $i == 0 ? 0 : ((3 -(($length_minusOne%3)-$phase_minusOne) ) % 3);
        my $length = $CDSs[$i]->end()-$CDSs[$i]->start()+1; #o already checked for acceptable values of strand so no probs
        my $cds_numer = $i+1;

        if ($CDSs[$i]->phase() eq q{.}) {
              #r complain or not?!?
              $CDSs[$i]->phase($phase); #y put in value
        }
        elsif ($CDSs[$i]->phase() != $phase) {
              #r complain or not?!?
              $CDSs[$i]->phase($phase); #y put in correct value
        }

        #y store info for next iteration
        $length_minusOne = $length;
        $phase_minusOne = $phase;
    }

    #$self->_CDSs([@CDSs]);
    #$self->_exons([@exons]);
    $self->_CDSs(\@CDSs);
    $self->_exons(\@exons);
}

#b unnecessary overkill but for clarity using refaddr 

sub eTrnslConvert {

    my $self = shift;
    my @exons  = @{$self->exons};
    my @CDSs  = @{$self->CDSs};
    my $strand = $self->gffGene()->strand() eq '+' ? 1 : -1;

    my $processed_exons = 0;
    my $processed_cds = 0;
    my $sumed_cds_length = 0; #y keep track of the length of the coding region - i.e. should be len%3=0
    my $cds_num = 0; #y create counter for CDS looping within exon loop
    my $eExon_transcription_start_Exon; #y create vars to store important exons for later API assignment
    my $eExon_transcription_end_Exon;
    my $eExon_translation_start_Exon;
    my $eExon_translation_end_Exon;

    my $biotype = $self->gffGene()->biotype();
    #my $biotype = $self->biotype();

    EXONCDSLOOP:
    for my $exon_num (0..$#exons) {

        $processed_exons++;

        # $exon_stable_id||=q{};

        #g/ this is pretty pointless its by definition $exons[0]... - i.e. take this outside of loop
        #f/ putting things into the gff objects
        $self->_trans_start_exon($exons[$exon_num]) if ($exon_num == 0);
        $self->_trans_end_exon($exons[$exon_num]) if ($exon_num == $#exons); 

        #/ use gene's seqid
        my $landmark = $exons[$exon_num]->seqid();

#        if (!exists $slices->{$landmark}) 
        Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Build', type => 'Undefined::exon', 
          error => 'there is a strange referential problem for mRNA '.$self->id(),
          list => \@exons,
        ) if (!$exons[$exon_num]);
          #if (!$exons[$exon_num]->isa('GffDoc::Feature::Gff3::exon'));
        Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Build', type => 'Undefined::CDS', 
          error => 'there is a strange referential problem for mRNA '.$self->id(),
          list => \@CDSs,
        ) if (!$CDSs[$cds_num] && $self->biotype() eq 'protein_coding');
          #y silly, silly
          #if (!$CDSs[$cds_num]);
          #if (!$CDSs[$cds_num]->isa('GffDoc::Feature::Gff3::CDS'));

        #/ i.e. if @CDS == 0 / not protein_coding            
        if ($biotype ne 'protein_coding') { #y give -1 phase to exons if there are no CDS features 
            $exons[$exon_num]->_ephase_start(-1);
            $exons[$exon_num]->_ephase_end(-1);
        }

        #y match each CDS to an exon - i.e. if this doesn't match the exon on non-coding - 
        #y technically we're matching exons to CDS but there may be CDS-less exons
        if ( 
            ($biotype eq 'protein_coding') && ($exons[$exon_num]->{start} <=  $CDSs[$cds_num]->{end}) &&
            ($CDSs[$cds_num]->{start}   <=  $exons[$exon_num]->{end}) 
        ) {

            $processed_cds++;
            $sumed_cds_length += ($CDSs[$cds_num]->{end} - $CDSs[$cds_num]->{start} +1 ); #y keep check of the sumed CDS length for later

            #y check appropriate settings for phase calculations
            my $x = $strand == 1 ? 'start' : 'end';
            my $y = $strand == 1 ? 'end'   : 'start';

#            $self->_cdsExon_overlap_checks($cds_num, $exon_num, \@CDSs, \@exons, $x, $y, $eTrans->stable_id)

            #y assign ensembl start of exon phase (just 3 - Gffphase)%3
            my $ensembl_phase = 
            #w if starts don't overlap means its first CDS so starting phase of exon is -1
            $exons[$exon_num]->{$x} == $CDSs[$cds_num]->{$x} 
            ? (3-$CDSs[$cds_num]->{phase})%3 : -1;

            $exons[$exon_num]->_ephase_start($ensembl_phase);

            $exons[$exon_num]->_ephase_end(
                $exons[$exon_num]->{$y} == $CDSs[$cds_num]->{$y} 
                ? (
                #w if they have the same end and ensembl start phase is -1 then it must be the first CDS ->0
                (($ensembl_phase == -1 ? 0 
                #w the ensembl start is not -1 so we calculate the phase properly (ens_start_phase + length)%3
                : $ensembl_phase)
                + ($CDSs[$cds_num]->{end}-$CDSs[$cds_num]->{start}+1) ) % 3
                    )
                #w they don't end so it so this exon must be partially UTR at this end -> -1
                : -1
            );
                    
            $self->_trnsl_start_exon($exons[$exon_num]) if ($cds_num == 0);
            $self->_trnsl_end_exon($exons[$exon_num]) if ($cds_num == $#CDSs);

            #y absolute last thing we do CDS-exon matcher is increment CDS for next exon match 
            $cds_num++ if ($cds_num != $#CDSs); #w try/do next cds

                #w end of if clause on overlapping cds/exons 
        } else {

            #y doesn't have a matching CDS so its non-coding
            $exons[$exon_num]->_ephase_start(-1);
            $exons[$exon_num]->_ephase_end(-1);
        }

    } #w end exon loop

    #y cds QC and saving EnsTranslation details
    if ($biotype eq 'protein_coding') { #if ($CDS_loop) 

        Exception::GffDoc->throw( stage => 'Build', type => 'Undefined::TRSNLstartExon', error => 'Unable to identify translation start exon for mRNA: '
          .$self->id().qq{\n\t}.'Are the exons/CDS properly paired for this mRNA?'
        ) if (!$self->trnsl_start_exon());
        Exception::GffDoc->throw( stage => 'Build', type => 'Undefined::TRSNLendExon', error => 'Unable to identify translation end exon for mRNA: '
          .$self->id().qq{\n\t}.'Are the exons/CDS properly paired for this mRNA?'
        ) if (!$self->trnsl_end_exon());

        #y check for CDS that haven't been matched to a gene
        if ($cds_num != $#CDSs) {
            my $mapped_cds_num = $cds_num + 1; # i.e. adjust for 0..-1index
            my $n = scalar @CDSs;
    #        . qq{* Unable to allocate all CDS regions to exons for transcript $stable_id }
    #        . qq{ (mapped $mapped_cds_num CDS out a total of $n. Skipping transcript).});
        } else {

            #y set first translated exon - i.e. first exon to map to CDS (arrays were ordered by strand orientation already)

            $self->_trnsl_start_base(
                _exon_coord(
                    $self->trnsl_start_exon(),
                    $strand == 1 ? $CDSs[0]->{start} : $CDSs[0]->{end}, 
                    $strand,
                )
            );

            $self->_trnsl_end_base(
                _exon_coord(
                    $self->trnsl_end_exon(),
                    $strand == 1 ? $CDSs[$#CDSs]->{end} : $CDSs[$#CDSs]->{start}, 
                    $strand,
                )
            );

        }
    }

}

sub _exon_coord {
    my ($exon,$coord,$strand) = @_;

    #if ($exon->strand == 1) {
    if ($strand == 1) {
        my $start = $exon->start;
        return $coord - $start + 1;
    } 
    else {
        my $end = $exon->end;
        return $end - $coord + 1;
    }    
}

__PACKAGE__->meta->make_immutable();

1;

#g could always instantiate classes otf Class::MOP::Class->create($package_name, %options),
#g $metaclass->get_method($method_name), $metaclass->get_method_list, $metaclass->add_method($method_name, $method)
#g $metaclass->get_all_methods$metaclass->get_all_methods, $metaclass->get_attribute($attribute_name)

package Bio::EnsEMBL::Exon;

sub toGeneDocString {
    
    my $thing = shift;

    return +{
        seqid => $thing->slice()->seq_region_name(),
        strand => $thing->strand() == 1 ? '+' : '-',
        start => $thing->start(),
        end => $thing->end(),
        id => $thing->stable_id,
    };
}

1;

package Bio::EnsEMBL::Transcript;

sub toGeneDocString {
    return &Tools::_stringifyGene(shift);
}

1;

package Bio::EnsEMBL::Gene;
#y careful, circumventing constructor and therefore type controls etc.                    
#my $gffGene = bless $gene, 'GffDoc::Gff3::Gene';                                          
#                                                                                            
#$log4->info('Reconstructing gene: '.$gene->id());                                         
#                                                                                            
#for my $kiddymRNA (@{$kiddymRNAs}) {                                                      
#                                                                                            
#    my $gffmRNA = bless $kiddymRNA, 'GffDoc::Gff3::mRNA';                                 
#                                                                                            
#    $gffmRNA->_gffGene($gffGene);  two-way mapping                                       
#                                                                                            
#    my $kiddyname = $kiddymRNA->id();                                                     
#    my $kiddyexons = $gffdoc->exonArray->findParentsDestructive($kiddyname);              
#    my $kiddyCDSs = $gffdoc->CDSArray->findParentsDestructive($kiddyname);                
#                                                                                            
#     we allow for no CDS for pseudogenes - till the genes are constructed we don't know which is which
#    / have to make this biotype dependent?!?                                             
#    Exception::GffDoc::Gene->throw( stage => 'Build', type => 'FeatureLine::Dangling::NoChildren',
#        error => 'we do not tolerate mRNAs without exons perhaps the mRNA ID= and '         
#        .'exon Parent= identifiers are not consistent?!? mRNA id='.$kiddyname               
#        . '. See Build log.', $gene->id(),                                                  
#    ) if (scalar @{$kiddyexons} == 0);                                                    
#                                                                                            
#    / build the mRNA                                                                     
#    print qq{\nexons: }, scalar @{$kiddyexons};                                          
#    print qq{\nCDSs: }, scalar @{$kiddyCDSs};                                            
#    $gffmRNA->_exons($kiddyexons);                                                        
#    $gffmRNA->_CDSs($kiddyCDSs);                                                          
#    $gffGene->addmRNA($gffmRNA);                                                          
#}                                                                                         
#$gffdoc2->GffGenes()->add($gffGene);    

sub toGeneDocString {
    return &Tools::_stringifyGene(shift);
}
    
#y prolly out to be ->new_from_eGene?!?

sub toGeneDocGene {

    print qq{\n > experimental!};
    print qq{\n here: }, $_ for (@_);
    print qq{\n};

    return;
}



1;

package Tools;

#/ can only access global vars with scope-resolution - so to access lexicals use sub returns...

#my $n = qq{\n};
#my $t = qq{\t};
#my $qual_kv_pair = qq{\t\t\t};
our $n = qq{\n};
our $t = qq{\t};
our $qual_kv_pair = qq{\t\t\t};


sub _stringifyGene {
    
    my $thing = shift;
    return +{
        seqid => $thing->slice()->seq_region_name(),
        strand => $thing->strand() == 1 ? '+' : '-',
        start => $thing->start(),
        end => $thing->end(),
        #has 'phase' => (is => 'rw', isa => enum([qw/ 1 2 3 . /]) );
        biotype => $thing->biotype,
        #ID => $thing->stable_id,
        id => $thing->stable_id,
    };
}



1;







######## toTblString ######### best to clean up that code before trying to really integrate it!?!

