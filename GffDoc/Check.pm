#/ this module is responsible for doing pre-commitment checks, but 'also' generate info for e! db storage
#/ of genes - this can be turned off?!? - or modified?!?
package GffDoc::Check;
use GffDoc::Exceptions qw/moan moan_e $nt moan_ref/;
use GffDoc::print;
use Mouse;
with 'GffDoc::Activity::Role';
has 'biotype_resolutions' => (is => 'ro', isa => 'HashRef');
has 'mRNA_biotype_dominant' => (is => 'ro', isa => 'Int');
has 'gene_biotype_dominant' => (is => 'ro', isa => 'Int');

sub run {

    my $self = shift;

    my $gffdoc      = $self->gffdoc();
    my $log4        = $self->log4();
    my $prob_genes  = $self->prob_genes();

    my $genes = scalar @{$gffdoc->GffGenes()->features()};

    my $c = 0;
    my $seqid_exception = 0;
    my $strand_exception = 0;

    GFFGENE:
    for my $gffGene (@{$gffdoc->GffGenes()->features()}) {

        $log4->info('Checking gene '.$gffGene->id());

        #/ for seqedits, translation name and start/end coords for dumbass gff formats... - use -polypeptide flag now?!?

        #/ for rebuilding we'd use -M and bring in polypeptides in the proceeding step...

#o/
eval {
#o/
        if ($c == 0 || $c % 10 == 0) {  
            my $p = sprintf(q{%.2f},100*$c/$genes);
            #GffDoc::print::print_clean(9);
            GffDoc::print::print_number($p.'%');
            GffDoc::print::print_clean(6);
        }
        $c++;

        #my $strand = $gffGene->strand() eq '+' ? 1 : -1,
        my $seqid = $gffGene->seqid(),
        my $strand = $gffGene->strand();

        #y we must have decided upon a biotype before we go down this route...
        #$biotypes{$gffGene->biotype()} = 1;

        my %biotypes;
        my %mRNA_biotypes;
        my $final_biotype;
        
        my $biotype_g = $gffGene->biotype();
        @biotypes{(map { $_->biotype() } @{$gffGene->mRNAs}), $biotype_g} = ();
        #@biotypes{(map { $_->biotype() } @{$gffGene->mRNAs}), $gffGene->biotype()} = ();
        @mRNA_biotypes{map { $_->biotype() } @{$gffGene->mRNAs}} = ();


        #r THIS WAS WRITTEN IN A HURRY COS OF TAIR... SHOULD BE A TERNARY, BUT CAN'T DECIDE ON SEQUENCE OF EVENTS and want an exception thrown...
        #y first check is unnecessary but... pseudogene always over-rides...
        #if (exists $mRNA_biotypes{'pseudogene'} || $gffGene->biotype() eq 'pseudogene') {
        if (scalar (keys %biotypes) == 1) {
            ($final_biotype) = (keys %biotypes);
            #($final_biotype) = (keys %mRNA_biotypes);
        #y Tair break puts biotype at gene sometimes so we can't just have mRNA biotype dominate...
        #} elsif (scalar (keys %mRNA_biotypes) == 1) {
        } elsif (exists $biotypes{'pseudogene'}) {
            $final_biotype = 'pseudogene';
        } elsif (scalar (keys %mRNA_biotypes) == 1 && $self->mRNA_biotype_dominant()) {
            ($final_biotype) = (keys %mRNA_biotypes);
        } elsif ($self->gene_biotype_dominant()) {
            $final_biotype = $biotype_g;
            die qq{\nhey: $final_biotype} if $final_biotype eq 'miRNA';
        } else {

            if ($self->biotype_resolutions()) {

                my @bts = sort { $a cmp $b } (keys %biotypes);
                my $str = join ':', @bts;

                if (exists $self->biotype_resolutions()->{$str}) {

                        $final_biotype = $self->biotype_resolutions()->{$str};

                } else {
            
                    Exception::GffDoc->throw( 
                    #Exception::GffDoc::Gene::InconsistentBiotype->throw( 
                      type => 'InconsistentBiotype',
                      error => 'Inconsistent biotypes across gene '.$gffGene->id(), 
                    );
                }
            }
        }    

        #y mRNAs get their biotype from their parent gene so we set that
        $gffGene->_biotype($final_biotype);

        for my $mRNA (@{$gffGene->mRNAs}) {

            #y override biotype
            #$mRNA->_biotype($final_biotype);

            #y these routines generate info for e! genes. 
            $mRNA->orderNPhaseCalc();
            $mRNA->eTrnslConvert();

            Exception::GffDoc::Gene::SeqRegion->throw( 
                error   => 'mRNA must be located on same SeqRegion as gene: '.$gffGene->id().'. Skipping gene',
                id      => $gffGene->id(),
                stage   => 'Check', 
                type    => 'SeqRegion::Mismatch', 
            ) if ($mRNA->seqid() ne $seqid);

            Exception::GffDoc::Gene::Strand->throw( 
                error => 'mRNA must have same sense as gene: '.$gffGene->id().'. Skipping gene',
                id      => $gffGene->id(),
                stage => 'Check', 
                type => 'Strand::Mismatch', 
            ) if ($mRNA->strand() ne $strand);

            #$biotypes{$mRNA->biotype()} = 1; #y should have sorted biotype from comments at object instantiation...
            
        }

#o/
};
#o/

#//// all next GFFGENE here need to add gene id to blacklist and print to other file if option set!?!
my $e;
#/ catch at low level once, then generic level later...
if ( $e = Exception::Class->caught('Exception::GffDoc::Gene::SeqRegion') ) { 
    $log4->error($e->error);
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string ) if !$seqid_exception;
    $seqid_exception = 1;
    push @{$prob_genes}, $e->id();
    next GFFGENE;
} elsif ( $e = Exception::Class->caught('Exception::GffDoc::Gene::Strand') ) { 
    $log4->error($e->error);
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string ) if !$strand_exception;
    $strand_exception = 1;
    push @{$prob_genes}, $e->id();
    next GFFGENE;
#y catch any others?!?
} elsif ( $e = Exception::Class->caught('Exception::GffDoc::ReferentialIntegrity') ) { 
    moan_ref ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string,
        scalar @{$e->list}, @{$e->list}[0..3] 
    );
    $log4->error($e->error);
    exit;
} elsif ( $e = Exception::Class->caught('Exception::GffDoc') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    $log4->error($e->error);
    exit;
} elsif ( $e = Exception::Class->caught('Exception') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing shoulbe be caught at this level.', $@ );
    exit;
}

    }    

    return $gffdoc;

} 

__PACKAGE__->meta->make_immutable();

1;
