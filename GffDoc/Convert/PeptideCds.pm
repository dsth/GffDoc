package GffDoc::Convert::PeptideCds;
use GffDoc::Exceptions qw/moan moan_e $nt/;
use Data::Dumper;
use Mouse;
with 'GffDoc::Activity::Role';

sub run {

    my ($self, $convert) = @_;

    my ($before,$after);
    if ($convert) {
        my @temp = split ':', $convert;

        Exception::GffDoc->throw( stage => 'Convert', type => 'IncorrectUsage', 
          error => 'Name conversion requires two parameters'
        ) if (scalar @temp != 2);

        ($before,$after) = @temp;
    }

    my $gffdoc = $self->gffdoc();
    my $log4 = $self->log4();

    my $meta = $gffdoc->meta();

    my $excess = 0;
    my $missing_e = 0;
    my $missing_c = 0;

    $missing_e = 1 if (!$gffdoc->exonArray());
    $missing_c = 1 if (!$gffdoc->CDSArray());
    $missing_e = 1 if (!$gffdoc->exonArray()->features());
    $missing_c = 1 if (!$gffdoc->CDSArray()->features());

    for my $attr ($meta->get_attribute_list()) {

        next if (!$gffdoc->$attr());
        next if (!$gffdoc->$attr()->features());
        
        $log4->info('There are ', $gffdoc->$attr()->length, ' features within the ', $attr);

        next if ($attr eq 'CDSArray' || $attr eq 'exonArray');

        $excess = 1;
    }




#    if ($mode eq 'e') {
#        Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
#          error => 'There should be both exon and CDS features AND no more to run this conversion module'
#            .qq{\n\t}.'ANY EXONS LACKING CDS WILL BE MADE PSEUDOGENES'
#        ) if ( $excess || $missing_c || $missing_e );
#    } elsif ($mode eq 's') {
#        Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
#          error => 'There should be only CDS features AND no more to run this conversion module'
#        ) if ( $excess || $missing_c || !$missing_e );
#    } else {
#        Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
#          error => 'Need to run this module in either standard -M GffDoc::Convert::Gtf:s or extended '
#          .'mode -M GffDoc::Convert::Gtf:e'
#        );
#    }


    my %mRNA_coords;
    my %stuff;

    my %genes;
    my %mRNAs;

    #y this shouldn't be destructive - duh!?: while (my $exon = $gffdoc->exonArray()->shift()) or regular findParents
    for my $peptide (@{$gffdoc->polypeptideArray()->features()}) {

        my $mRNAid = $peptide->id();
        if ($convert) { $mRNAid =~ s/$before/$after/; }
 #       print qq{\n},$mRNAid;

        #/ really unnecessary, just can't be arsed to make polypeptides fully-features features...
        my ($mRNA) = grep { $_->id() eq $mRNAid } @{$gffdoc->mRNAArray()->features()};

        if ($mRNA) {

            $genes{$mRNA->parent()} = 1 if (!exists $genes{$mRNA->parent()}); #y no actual reason for condition
            $mRNAs{$mRNA->id()} = 1 if (!exists $genes{$mRNA->id()}); 

            my @exon_objs = grep { $_->parent() eq $mRNAid } @{$gffdoc->exonArray()->features()};


            @exon_objs = sort { $a->start() <=> $b->start() } @exon_objs;

            #print qq{\npep start: }, $peptide->start();
            #print qq{\npep end: }, $peptide->end();

            for my $exon (@exon_objs) {

                next if (($exon->end() < $peptide->start()) || ($exon->start() > $peptide->end()));

                my $CDS_start = ($exon->start <= $peptide->start && $exon->end >= $peptide->start) ? $peptide->start : $exon->start;
                my $CDS_end = ($exon->start <= $peptide->end && $exon->end >= $peptide->end) ? $peptide->end : $exon->end;

                #print qq{\n\texon start: }, $exon->start(), '    -    cds: ', $CDS_start;
                #print qq{\n\texon end: }, $exon->end(), '    -     cds: ',$CDS_end;

                my $Feature = GffDoc::Feature::Gff3::CDS->new(
                    biotype => 'protein_coding',
                    end     => $CDS_end,
                    id      => 'ID='.$mRNAid.'-E;Parent='.$mRNAid.';',
                    phase   => '.', 
                    seqid   => $mRNA->seqid(), 
                    source  => 'GffDoc::Convert::PeptideCds',
                    start   => $CDS_start, 
                    strand => $mRNA->strand(), 
                    #parent => $from_thing->parent(),
                );

                $gffdoc->CDSArray->add($Feature);

                #use Data::Dumper;
                #print Dumper $Feature;
                #print Dumper $exon;
            }

        } else { 
            Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', error => 'cannot find mRNA' );
        }

    }

    #y now that we've made cds we re-assign anything without polypeptides to psuedogenes

    my $g = 0;
    for my $gene (@{$gffdoc->geneArray()->features()}) {

        next if (exists $genes{$gene->id()});

        $gene->_biotype('pseudogene');
        $g++;

    }

    my $m = 0;
    for my $mRNA (@{$gffdoc->mRNAArray()->features()}) {

        next if (exists $mRNAs{$mRNA->id()});

        $mRNA->_biotype('pseudogene');
        $m++;

    }
    $log4->info('Generated ', $gffdoc->CDSArray->length(), ' features within the CDSArray');
    $log4->info('Changed biotype of ', $g, ' gene features within geneArray to pseudogene');
    $log4->info('Changed biotype of ', $m, ' mRNA features within mRNAArray to pseudogene');

}

__PACKAGE__->meta->make_immutable();

1;
