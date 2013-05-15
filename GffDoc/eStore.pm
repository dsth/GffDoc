package GffDoc::eSlices;
use Mouse;
use GffDoc::Exceptions qw/moan/;
with 'GffDoc::Activity::Role';
has 'slices' => (is => 'ro', writer => '_slices');
has 'dbad' => (is => 'ro', required => 1);
has 'cs_name' => (is => 'ro', required => 1);
has 'cs_version' => (is => 'ro', required => 1);

#g need proper exception handling in slice routines...

sub run {

    my $self = shift;
    my %seqids;

    my $cs_version  = $self->cs_version() ? $self->cs_version() : undef;

    #y/ just need ALL seq ids - not all element types exist at this point
    my @features;
    push @features, @{$self->gffdoc->geneArray()->features()} if $self->gffdoc->geneArray()->features();
    push @features, @{$self->gffdoc->mRNAArray()->features()} if $self->gffdoc->mRNAArray()->features();
    push @features, @{$self->gffdoc->CDSArray()->features()} if $self->gffdoc->CDSArray()->features();
    push @features, @{$self->gffdoc->exonArray()->features()} if $self->gffdoc->exonArray()->features();

    #y only to grab slices early and avoid waiting ages to have an inane issue like that kill the thing
    for my $feat (@features) {

        #print qq{\n}, ref $feat;
        $seqids{$feat->seqid()} = 1;
    }

    my @seqids = (keys %seqids);

    my $db_slice_adaptor = $self->dbad()->get_adaptor("Slice") || Exception::GffDoc->throw( stage => 'ensembl', type => 'Adaptor::Slice',
        error => 'Could not generate slice adaptor for e! db'
    );

    my %slices;
    for my $landmark (@seqids) {

        if (my $slice = $db_slice_adaptor->fetch_by_region($self->cs_name, $landmark, undef, undef, undef, $cs_version)) {
            $slices{$landmark}  = $slice;
        } else {
            Exception::GffDoc->throw( stage => 'ensembl', type => 'Slice::Lookup',
                error => 'Could not fetch slice for seqid '.$landmark.' in cs '.$self->cs_name()
            );
        }
    }

    #$self->_slices(\%slices);
    return \%slices;
}

__PACKAGE__->meta->make_immutable();

1;

package GffDoc::eStore;
use Mouse;
use GffDoc::Exceptions qw/moan moan_e $nt/;
with 'GffDoc::Activity::Role';

#g put in to_eObj() method into gene and just pass type and the appropriate stable for it to return e! object!?!
#g or have it in role and each type sets its e! type and so the appropriate one gets called

has 'logicname' => (is => 'ro', isa => 'Str', required => 1);
has 'version' => (is => 'ro', isa => 'Str', required => 1);
has 'analysis' => (is => 'rw');
has 'slices' => (is => 'ro', writer => '_slices', required => 1);
has 'dbad' => (is => 'ro', required => 1);
has 'cs_name' => (is => 'ro', required => 1);
has 'non_coding_cds' => (is => 'rw');

#g not bothering...
#sub eObject {
#    my ($self, $type, $gene, @rest) = @_;
#    my $e = 'Bio::EnsEMBL::'.$type;
#    return $e->new(
#        -ANALYSIS   =>  $self->analysis(),
#        -VERSION    =>  $self->version(),
#        -STRAND     =>  $gene->strand,
#        -BIOTYPE    =>  $gene->biotype,
#        @rest
#    );
#}

sub run {

    my $exception_missingslice = 0;
    my $gene_nosave = 0;
    my $gene_save = 0;
    my $non_coding = 0;

    my %exon_names;

    my $self = shift;

    my $db_ad       = $self->dbad(); # handle is [2][1] # slices are [2][2][0]
    my $gffdoc      = $self->gffdoc();
    my $log4        = $self->log4();
    my $prob_genes  = $self->prob_genes();
    my $slices      = $self->slices();

    my %eException_msgs;

    $log4->info('Starting commitment.');

    my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $self->logicname());

    $self->analysis($analysis);

    my $genes = scalar @{$gffdoc->GffGenes()->features()};

    my $c = 0;

    GFFGENE:
    for my $gffGene (@{$gffdoc->GffGenes()->features()}) {

        $log4->info('Storing gene '.$gffGene->id());
        my $gid = $gffGene->id();

#o/
eval {
#o/

        #/################### wrap this in eval and catch at the end of the loop with a next GFFGENE...

        if ($c == 0 || $c % 10 == 0) {
            my $p = sprintf(q{%.2f},100*$c/$genes);
            #GffDoc::print::print_clean(9);
            GffDoc::print::print_number($p.'%');
            GffDoc::print::print_clean(6);
        }
        $c++;

        #y should really already have been picked up on...
        Exception::GffDoc::Gene::MissingSlice->throw( stage => 'ensembl', type => 'MissingSlice', error => 'cannot find
            slice '.$gffGene->seqid().' for gene '.$gffGene->id(), id => $gffGene->id(),
        )

        if (!exists $slices->{$gffGene->seqid()});

        my $slice = $slices->{$gffGene->seqid()},
        my $strand = $gffGene->strand() eq '+' ? 1 : -1,

        my $biotype = $gffGene->biotype();

        my $eGene = Bio::EnsEMBL::Gene->new(
                -ANALYSIS       =>  $self->analysis(),
                -VERSION        =>  $self->version(),
                -SLICE          =>  $slice,
                -STRAND         =>  $strand,
                -START          =>  $gffGene->start(), # really no need the API bins it and calculates it from transcripts
                -END            =>  $gffGene->end(),
                -STABLE_ID      =>  $gffGene->id(),
                -BIOTYPE        =>  $biotype,
                -CREATED_DATE   =>  time,
                -MODIFIED_DATE  =>  time,
            );

        for my $mRNA (@{$gffGene->mRNAs}) {

             my $eTrans = Bio::EnsEMBL::Transcript->new(
                    -ANALYSIS       =>  $self->analysis(),
                    -VERSION        =>  $self->version(),
                    -SLICE          =>  $slice,
                    -STRAND         =>  $strand,
                    -START          =>  $mRNA->start(),
                    -END            =>  $mRNA->end(),
                    -BIOTYPE        =>  $mRNA->biotype(),
                    -STABLE_ID      =>  $mRNA->id(),
                    -CREATED_DATE   =>  time,
                    -MODIFIED_DATE  =>  time,
                );

			# Edit by JA to use polypeptide ids, if we have 'em
			my $canonical_TRNSL;
			my $alt_TRNSLs = 0;
            if ($biotype eq 'protein_coding') {
				my @translations;
				if (defined $mRNA->polypeptides() && scalar(@{$mRNA->polypeptides()}) > 0) {
            		for my $polypeptide (@{$mRNA->polypeptides()}) {
						push @translations, $polypeptide->id();
					}
				} else {
                	my $TRNSLid = $mRNA->id();
                	$TRNSLid =~ s/-R(\w)$/-P$1/;
					push @translations, $TRNSLid;
					# Show warning if $TRNSLid eq $mRNA->id() ?
                }

				# Health warning: if polypeptides > 1, we assume that
				# the first is the canonical translation, but have
				# done nothing to ensure that is the case.
				my $canonical = 1;
				foreach my $TRNSLid (@translations) {
            		my $eTRNSL; #y peptide so caplitalise...
                    $eTRNSL = Bio::EnsEMBL::Translation->new(
                        -ANALYSIS       =>  $self->analysis(),
                        -VERSION        =>  $self->version(),
                        -SLICE          =>  $slice,
                        -STABLE_ID      =>  $TRNSLid,
                        -STRAND         =>  $strand,
						-CREATED_DATE   =>  time,
						-MODIFIED_DATE  =>  time,
                    );

                	#y bind translation-transcript
					if ($canonical) {
                		$eTrans->translation($eTRNSL);
						$canonical_TRNSL = $eTRNSL;
						$canonical = 0;
					} else {
						$eTrans->add_alternative_translation($eTRNSL);
						$alt_TRNSLs = 1;
					}
				}
            }

            ENSEXON:
            my @exons = @{$mRNA->exons()};

            for my $n (0..$#exons) {

                #y/ exon unique names
                # JA: The previous method didn't work if more than two exons had the same
                # location/ID but different phases.
                my $exonid = $exons[$n]->id() ? $exons[$n]->id() : $gffGene->id().'-E'.($n+1);
                if (exists $exon_names{$exonid}) {
                    $exonid .= '-'.$exon_names{$exonid};
                }
                $exon_names{$exonid}++;

                 my $eExon = Bio::EnsEMBL::Exon->new(
                        -ANALYSIS       =>  $self->analysis(),
                        -VERSION        =>  $self->version(),
                        -STABLE_ID      =>  $exonid,
                        -SLICE          =>  $slice,
                        -STRAND         =>  $strand,
                        -START          =>  $exons[$n]->start(),
                        -END            =>  $exons[$n]->end(),
                        -PHASE          =>  $exons[$n]->ephase_start(),
                        -END_PHASE      =>  $exons[$n]->ephase_end(),
                        -CREATED_DATE   =>  time,
                        -MODIFIED_DATE  =>  time,
                    );

                    $exons[$n]->_eExon($eExon); # lazy but still let system be extendible to multiple translations

                    $eTrans->add_Exon($eExon);

            }

            $eTrans->start_Exon($mRNA->trans_start_exon()->eExon());
            $eTrans->end_Exon($mRNA->trans_end_exon()->eExon());

            if ($biotype eq 'protein_coding') {

#                if ($sumed_cds_length % 3 != 0) {
#                if ($cds_num != $#{$CDS_refs}) {

                $canonical_TRNSL->start_Exon($mRNA->trnsl_start_exon()->eExon());
                $canonical_TRNSL->start($mRNA->trnsl_start_base());
                $canonical_TRNSL->end_Exon($mRNA->trnsl_end_exon()->eExon());
                $canonical_TRNSL->end($mRNA->trnsl_end_base());

				if ($alt_TRNSLs) {
					foreach my $TRNSL (@{$eTrans->get_all_alternative_translations()}) {
						$TRNSL->start_Exon($mRNA->trnsl_start_exon()->eExon());
                		$TRNSL->start($mRNA->trnsl_start_base());
                		$TRNSL->end_Exon($mRNA->trnsl_end_exon()->eExon());
                		$TRNSL->end($mRNA->trnsl_end_base());
					}
				}

                my $seq = $canonical_TRNSL->seq();

                if(!$seq || $seq eq '' || $seq =~ /\*/) {

                    if ($self->non_coding_cds()) {
                        $eTrans->biotype('nontranslating_CDS');
                        $eTrans->translation(undef);
                        $eGene->biotype('nontranslating_CDS');
                        $non_coding++;
                        $log4->warn('Translation for transcript '.$mRNA->id().' contains stop codons. Re-classifying biotype as nontranslating_CDS');
                    } else {
                        Exception::GffDoc::StopCodons->throw( stage => 'ensembl', type => 'TRNSL::StopCodons',
                            error => 'Stop codons are not tolerated within Translations unless you run with -non_coding_cds option!'
                        );
                    }

                    #print qq{\n\nseq }, $seq;
                }

            }
            $eGene->add_Transcript($eTrans);

        }


        $db_ad->get_GeneAdaptor->store($eGene) ;

#o/
};
#o/

my $e;
#y/ print for every new msg
#y have to use e! exception system here. just pull out bits and pieces - i.e. message as hash - so that we don't give
#y the same one twice - and put the rest in the log...
if ( $e = Exception::Class->caught('Exception::GffDoc::StopCodons') ) {
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught('Exception::GffDoc::Gene::MissingSlice') ) {
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string ) if (!$exception_missingslice);
    $exception_missingslice = 1;
    $log4->error('Slice missing for gene '.$gffGene->id());
    push @{$prob_genes}, $e->id;
    $gene_nosave++;
    next GFFGENE;
} elsif ( $e = Exception::Class->caught() ) {
    my $msg = $@;
    my $trace;
    if ($msg =~ /MSG:\s+(.*?)(STACK.*?)Ensembl API version/s) { $msg = $1; $trace = $2 }
    chomp $msg;
    if (!exists $eException_msgs{$msg}) {
                moan ( 'EnsEMBL API', 'Gene '.$gid.': '.$msg, $trace );
        $eException_msgs{$msg} = 1;
    }
    $log4->error('Gene '.$gid.': '.$msg);
    push @{$prob_genes}, $gid;
    $gene_nosave++;
    next GFFGENE;
} else {
    $gene_save++;
}

    }

    return [$gene_save, $gene_nosave, $non_coding];
}

sub _sum {
    my $s = 0;
    $s += $_ for (@_);
    return $s;
}

__PACKAGE__->meta->make_immutable();

1;
