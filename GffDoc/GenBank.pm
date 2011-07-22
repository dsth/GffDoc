package GenBank::Submission;
use GffDoc::Exceptions;
use GffDoc::Types;
use Mouse;

#/ type checking only enforced AFTER buildargs so no problem in generating appropriate structures from wrong types...
has 'block_repeat_names' => (is => 'ro' );
has 'blist' => (is => 'ro' );
has 'gene_count' => (is => 'ro' );
has 'genes_processed' => (is => 'ro', default => 0, writer => '_genes_processed');
has '_repeat_names' => (is => 'ro', isa => 'HashRef', auto_deref => 1);
has '_TranslationStables2OldLocusTags' => (is => 'ro', isa => 'HashRef');
has '_TranslationStables2PeptideAccessionsNseqIDs' => (is => 'ro', isa => 'HashRef');
has '_NucleotideIds' => (is => 'ro', isa => 'HashRef');
has 'locus_tag_prefix' => (is => 'ro', isa => 'Str', required => 1);
has 'wgs_proj_prefix' => (is => 'ro', isa => 'Str', required => 1);
has 'coordsystem' => (is => 'ro', isa => 'Str', required => 1);
has '_outdir' => (is => 'ro', isa => 'Str');
has '_short_introns' => (is => 'ro');
has 'split_genes' => (is => 'ro');
has 'n2h' => (is => 'ro');
#has 'alternative_haplotype' => ( is => 'ro' );
has 'vectorbase_dbxref' => (is => 'ro');
has 'heteregeneous_population' => (is => 'ro');
has '_ncRNAs' => (is => 'ro');
has '_field_sep' => (is => 'ro', isa => 'Str');
has 'e_seq_region_name' => (is => 'ro', isa => 'Str'); # states which of Gb accession/SeqId is used as internal seqregion id - use as symbolic for accession vs. SeqId 
has 'chromosomes' => (is => 'ro');
has 'species_strain' => ( is => 'ro', required => 1 );
#has 'species_file' => (is => 'ro', isa => 'Str', required => 1);
#has 'species' => (is => 'ro', isa => 'Str');
    
my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @a = localtime; 
my $year = $a[5]-100;
my $date = $year.$abbr[$a[4]].$a[3];
my $field_sep = '_';

#r values only ever set once...
my $n = $Tools::n;
my $t = $Tools::t;
my $qual_kv_pair = $Tools::qual_kv_pair;

#/ this routine assumes genes are anchored to toplevel and thus are already there when called?!?
sub _grab_split_genes {

    my ($gene_adapt,$cs) = @_;

    my @genes = @{$gene_adapt->fetch_all};

    my %MultiSegmentGenes;
    my %MultiSegmentGenes_revmap;

#    if ($dev) {

    for my $gene (@genes) {

        if (scalar @{$gene->project($cs)} > 1) {

            ## project to the clone coordinate system
            $MultiSegmentGenes{$gene->stable_id} 
              = [map { $_->[2]->seq_region_name } @{$gene->project($cs)}];

            push @{$MultiSegmentGenes_revmap{$_}}, $gene->stable_id 
              for (map { $_->[2]->seq_region_name } @{$gene->project($cs)});
        }
    }

    #print qq{\n\nmulti-segment genes: }, scalar (keys %MultiSegmentGenes);
#    if ($dev == 2) {
#        store(\%MultiSegmentGenes, 'MultiSegmentGenes.frz') 
#            or die "Can't store MultiSegmentGenes!\n";
#        store(\%MultiSegmentGenes_revmap, 'MultiSegmentGenes_revmap.frz') 
#            or die "Can't store MultiSegmentGenes_revmap!\n";
#    }
#    }
#    else {
#        my $MultiSegmentGenes_ref = retrieve('MultiSegmentGenes.frz') 
#          or die "Can't rertieve \$MultiSegmentGenes\n";
#        %MultiSegmentGenes = %{$MultiSegmentGenes_ref};
#        my $MultiSegmentGenes_revmap_ref = retrieve('MultiSegmentGenes_revmap.frz') 
#          or die "Can't rertieve \$MultiSegmentGenes_revmap\n";
#        %MultiSegmentGenes_revmap = %{$MultiSegmentGenes_revmap_ref};
#    }
    return (\%MultiSegmentGenes, \%MultiSegmentGenes_revmap, scalar @genes);
}

sub _increment_processed_genes {
    my ($self) = @_;
    $self->_genes_processed($self->genes_processed+1);
    return;
}
#sub new {
#    my $class = shift;
#    return bless {
#        }, $class;
#}

around BUILDARGS => 
sub {

    my $orig  = shift;
    my $class = shift;
    my %hash = @_;

    #g never use these again after construction and dbad will cache loads of crap stuff and bloat everything...
    #has 'p2g' => (is => 'ro', isa => 'Str', required => 1);
    #has 'n2g' => (is => 'ro', isa => 'Str', required => 1);
    #has 'dbname' => (is => 'ro', isa => 'Str', required => 1);
    #has 'dbad' => (is => 'ro', required => 1);
    #has 'short_introns_file' => (is => 'ro');
    my $check   = !exists $hash{p2g}    ?  'p2g'
                : !exists $hash{n2g}    ?  'n2g'
                : !exists $hash{dbname} ?  'dbname'
                : !exists $hash{dbad}   ?  'dbad'
                :                          q{};

    #r/ only relevant when locus tag has changed!?!
    #: !exists $hash{p2o}    ?  'p2o'

    Exception::GffDoc->throw ( 
        stage => 'GenBank', 
        type  => 'Option::'.$check, 
        error => 'you must supply the option '.$check,
    ) if ($check);

    $hash{_field_sep} = $field_sep;

    my $ncRNA_types = qr{antisense_RNA|autocatalytically_spliced_intron  
        |ribozyme|hammerhead_ribozyme|RNase_P_RNA     
        |RNase_MRP_RNA|telomerase_RNA|guide_RNA       
        |rasiRNA|scRNA|siRNA|miRNA|piRNA|snoRNA|snRNA
        |SRP_RNA|vault_RNA|Y_RNA}x; # must have this format with x modifier...

    $hash{_ncRNAs} = $ncRNA_types;

    if ($hash{dbname} =~ /([a-zA-Z]+_[a-zA-Z]+)_core/) { 
        $hash{_outdir} = $1;
    } else { 
        Exception::GffDoc->throw ( 
            stage => 'GenBank', 
            type  => 'DbName', 
            error => 'database name should be of the form genus_species%' 
        );
    }

######################################
#/ mappings to Gb Accession and SeqIDs
######################################

    my %TranslationStables2PeptideAccessionsNseqIDs;

############ Those generated with translation stable_id = SeqID

    #r for some odd reason this only works with 'or' not '||'
    open my $pfh, '<', $hash{p2g} or Exception::GffDoc->throw ( 
            stage => 'GenBank', 
            type  => 'File::p2g', 
            error => 'cannot read '.$hash{p2g},
    );

    while (<$pfh>) {
    #trnslStableId=protSeqId  nucAccession
        next if (/^\s*#/);
        chomp;
        my @fields = split qq{\t};
        $TranslationStables2PeptideAccessionsNseqIDs{$fields[1].$field_sep.$fields[0]} 
          = { protAccession => $fields[1], protSeqId => $fields[0] };
    }

    close $pfh || die;

############ For earlier genomes a variety of things were used as SeqIDs bring those in too

    #r can't be bothered to remove those entries that already have decent SeqIds as they will only get over-written by newer entries...
    if ($hash{p2g_relic}) {

        open $pfh, '<', $hash{p2g_relic} or Exception::GffDoc->throw ( 
                stage => 'GenBank', 
                type  => 'File::p2g', 
                error => 'cannot read '.$hash{p2g_relic},
        );

        while (<$pfh>) {
        #trnslStableId  nucAccession    protAccession   protSeqId       nucSeqId
            next if (/^\s*#/);
            chomp;
            my @fields = split qq{\t};
            $TranslationStables2PeptideAccessionsNseqIDs{$fields[1].$field_sep.$fields[0]} 
            = { protAccession => $fields[2], protSeqId => $fields[3] };
        }
        close $pfh;

    }

    $hash{_TranslationStables2PeptideAccessionsNseqIDs} = \%TranslationStables2PeptideAccessionsNseqIDs;

####################
    open my $ofh, '<', $hash{p2o} or Exception::GffDoc->throw ( 
            stage => 'GenBank', 
            type  => 'File::p2o', 
            error => 'cannot read '.$hash{p2o},
    );

    my %TranslationStables2OldLocusTags;
    while (<$ofh>) {
    #trnslStableId  nucAccession    protAccession   protSeqId       nucSeqId
        next if (/^\s*#/);
        chomp;
        my @fields = split qq{\t};
        my @oldids = split q{,}, $fields[2];
        $TranslationStables2OldLocusTags{$fields[1].$field_sep.$fields[0]} 
          = [@oldids];
    }
    close $pfh;

    $hash{_TranslationStables2OldLocusTags} = \%TranslationStables2OldLocusTags;
####################
    #b nucleotide accessions to nucleotide seqid mapping is done separately from p2g as this is peptide-centric
    #b and thus excludes those nucleotide sequences that have gained annotation between submissions

    open my $nfh, '<', $hash{n2g} or Exception::GffDoc->throw ( 
            stage => 'GenBank', 
            type  => 'File::n2g', 
            error => 'cannot read '.$hash{n2g},
    );
    my %NucleotideAccessions2NucleotideSeqIDs;
    my %NucleotideSeqIDs2NucleotideAccessions;
    while (<$nfh>) {
    #nucSeqId       nucAccession
        next if (/^\s*#/);
        chomp;
        my @fields = split qq{\t};
        $NucleotideAccessions2NucleotideSeqIDs{$fields[1]} = $fields[0];
        $NucleotideSeqIDs2NucleotideAccessions{$fields[0]} = $fields[1];
    }
    close $nfh;
    
    $hash{_NucleotideIds} = {
        toSeqId => \%NucleotideAccessions2NucleotideSeqIDs,
        toAccession => \%NucleotideSeqIDs2NucleotideAccessions,
    };

    #r we use UniParcRetrivalForShirtIntrons.pl internally...
    if ($hash{short_introns_file}) {

        open my $sifh, '<', $hash{short_introns_file} or Exception::GffDoc->throw ( 
                stage => 'GenBank', 
                type  => 'File::n2g', 
                error => 'cannot read '.$hash{short_introns_file},
        );
        my %hashOfShort;
        while (<$sifh>) {
        #trnslStableId  PepSeq  INSDentry       Method
            next if (/^\s*#/);
            chomp;
            my @fields = split qq{\t};
            #my $insd = $fields[2] ? qq{insd => $fields[2],} : q{};
            $hashOfShort{$fields[0]} = { seq => $fields[1], insd => $fields[2] };
        }
        close $sifh;
    
        $hash{_short_introns} = \%hashOfShort;
    }

    if ($hash{n2h}) {

        open my $hfh, '<', $hash{n2h} or Exception::GffDoc->throw ( 
                stage => 'GenBank', 
                type  => 'File::n2h', 
                error => 'cannot read '.$hash{n2h},
        );

        my %NucleotideAccession2AlternativeHaplotype;
        while (<$hfh>) {
            next if (/^\s*#/);
            chomp;
            my @fields = split qq{\t};
            #push @{$NucleotideAccession2AlternativeHaplotype{$fields[0]}, [@fields[1.3]];
            push @{$NucleotideAccession2AlternativeHaplotype{$fields[0]}}, 
              {start => $fields[1], end => $fields[2], comment => $fields[3]};
        }
        close $hfh;
        $hash{n2h} = \%NucleotideAccession2AlternativeHaplotype;
    }

    if ($hash{blist}) {
        open my $bfh, '<', $hash{blist} or Exception::GffDoc->throw ( 
                stage => 'GenBank', 
                type  => 'File::blist', 
                error => 'cannot read '.$hash{blist},
        );
        my %blist;
        while (<$bfh>) {
            next if (/^\s*#/);
            chomp;
            $blist{$_} = 1;
        }
        close $bfh;
        $hash{blist} = \%blist;
    }

    if ($hash{block_repeat_names}) {

        my $gene_adapt = $hash{dbad}->get_GeneAdaptor;

        my %repeat_names;
        for my $g (@{$gene_adapt->fetch_all()}) {
            my $en = $g->external_name();
            next if (!$en);
            $repeat_names{$en}++;
        }
        my @list = grep { $repeat_names{$_} > 1 } (keys %repeat_names);

        #y being silly...
        undef %repeat_names;
        @repeat_names{@list} = ();

        $hash{_repeat_names} = \%repeat_names;
    }

    my @MultiSegmentGenes;
    my %MultiSegmentGenes_map;
    my %MultiSegmentGenes_rev_map;

    if ($hash{split_genes}) {
    # identify split genes - as yet this is pest specific as only pest is submitted at a funny level - however, projection builds?!?

        my $gene_adapt = $hash{dbad}->get_GeneAdaptor;

        Exception::GffDoc->throw ( 
                stage => 'GenBank', 
                type  => 'Args::coordsystem', 
                error => 'you must supply a coordsystem value',
        ) if (!exists $hash{coordsystem});

        print qq{\n> Identifying split genes\n};
        my @MultiSegmentGenes = &_grab_split_genes($gene_adapt,$hash{coordsystem}); #r dev setting
        my %MultiSegmentGenes_map = %{$MultiSegmentGenes[0]};
        my %MultiSegmentGenes_rev_map = %{$MultiSegmentGenes[1]};

        $hash{split_genes} = {
            toSeqId => \%MultiSegmentGenes_map,
            toGeneId => \%MultiSegmentGenes_rev_map,
        };
    }

    $hash{_outdir} = $hash{_outdir}.'_'.$date;

    #/// have this as an option?!?

    if(-d $hash{_outdir}) { 
        #print qq{\nwiping previous dir }.$hash{_outdir}.' contents';
        # unlink $outdir.'/*' or die; 
        system 'rm -f '.$hash{_outdir}.'/*' || die;
    } else { 
        print qq{\ncreating output directory }.$hash{_outdir};
        mkdir $hash{_outdir} or die; 
        $hash{_outdir} = $hash{_outdir}.'/';
    };

    Exception::GffDoc->throw( 
        stage => 'Type', 
        type => 'Option::e_seq_region_name', 
        error => 'You must supply e_seq_region_name and it must be one of GbAccession or GbSeqId',
    ) if (!exists $hash{e_seq_region_name} 
        || (exists $hash{e_seq_region_name} 
        && ($hash{e_seq_region_name} ne 'GbAccession' 
        && $hash{e_seq_region_name} ne 'GbSeqId'))
    );

    my $gene_adapt = $hash{dbad}->get_GeneAdaptor;
    my $gene_count = scalar @{$gene_adapt->fetch_all()};
    $hash{gene_count} = $gene_count;

    return $class->$orig(%hash);
};

sub CheckHashOfShortIntrons {
    my ($self,$id) = @_;
    # return exists $hashOfShort{shift->} ? 1 : 0;
    # return exists $self->_short_introns->{$id} ? 1 : 0;
    return exists $self->_short_introns()->{$id} ? $self->_short_introns()->{$id} : 0;
}

#/put in subs to get all appropriate info to keep actual GffDoc classes clean

sub AvoidRepeatNames {

    my ($self,$id) = @_;

    my %repeat_names = $self->_repeat_names;

    if (exists $repeat_names{$id}) {
        $repeat_names{$id}++;
        $id .= sprintf(q{_%c},$repeat_names{$id}+96);
    } else {}

    # return $qual_kv_pair.'gene'.$t.$id.$n; 
    return $id; 

}

sub OldLocusTag {

    my ($self,$seqid,$id) = @_;

    my $prefix = $self->locus_tag_prefix();

    my $str = q{};

    #/ this is prolly an error - hard coding PA!?! - wipe it from file?!? or ignore it in GbSbm construction
    my $compound_name = $seqid.$self->_field_sep().$id.'-PA';
    #my $compound_name = $seqid.$self->_field_sep().$id;

    my $TranslationStables2OldLocusTags = $self->_TranslationStables2OldLocusTags();

#    use Data::Dumper; print Dumper $TranslationStables2OldLocusTags;
#print qq{\nsearching for }.$seqid.' '.$id.qq{\nthat is: }.$compound_name; die
#;
    #y the for old genes we put in previous locus names
    if (exists $TranslationStables2OldLocusTags->{$compound_name}) { 
        for my $old_tag (@{$TranslationStables2OldLocusTags->{$compound_name}}) { 

            next if ($prefix.$id eq $old_tag);
            $str .= $qual_kv_pair.'old_locus_tag'.$t.$old_tag.$n;
        }
    }
    else {}
    return $str;
}

1;
