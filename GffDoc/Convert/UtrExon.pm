package GffDoc::Convert::UtrExon;
use GffDoc::Exceptions qw/moan moan_e $nt/;
use Data::Dumper;
use Mouse;
with 'GffDoc::Activity::Role';

sub run {

    my ($self) = @_;

    my $gffdoc = $self->gffdoc();
    my $log4 = $self->log4();

    my $missing = 0;

    Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', error => 'buttmunch' )
      if (!$gffdoc->mRNAArray());
    Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', error => 'buttmunch' )
      if (!$gffdoc->mRNAArray()->features());

    $missing = 1 if (!$gffdoc->geneArray());
    $missing = 1 if (!$gffdoc->mRNAArray());
    $missing = 1 if (!$gffdoc->CDSArray());
    $missing = 1 if (!$gffdoc->geneArray()->features());
    $missing = 1 if (!$gffdoc->mRNAArray()->features());
    $missing = 1 if (!$gffdoc->CDSArray()->features());

    Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
      error => 'This module requires gene, mRNA and CDS feature types!'
    ) if ($missing);

    $log4->info('There are ', $gffdoc->geneArray()->length, ' features within the geneArray');
    $log4->info('There are ', $gffdoc->mRNAArray()->length, ' features within the mRNAArray');
    $log4->info('There are ', $gffdoc->CDSArray()->length, ' features within the CDSArray');

    my %cds_coords_by_mrna;

    for my $CDS (@{$gffdoc->CDSArray()->features()}) {

        my $parent = $CDS->parent();

        if (!exists $cds_coords_by_mrna{$parent}) {
            $cds_coords_by_mrna{$parent} = [];
        }
        push @{$cds_coords_by_mrna{$parent}}, [ $CDS->start(), $CDS->end() ];

    }

    #my @shorter = (keys %mrna_coords);
    #my ($prob1) = grep { !defined $mrna_coords{$_} } (keys %cds_coords_by_mrna);
    #my ($prob2) = grep { !defined $cds_coords_by_mrna{$_} } (keys %mrna_coords);

    while (my ($k,$v) = each %cds_coords_by_mrna) {

        #y not sure why bother re-ordering yet, but makes little difference...
        my @ordered_list = sort { $a->[0] <=> $b->[0] } @{$v};
        #my @list = sort { $b->[0] <=> $a->[0] } @{$v};

        if (scalar @ordered_list == 0) { 

            Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', error => 'this is a mess' );

        #r grab the coords of the mRNA and print as exon
        } elsif (scalar @ordered_list == 1) {

            my ($mRNA) = grep { $_->id() eq $k } @{$gffdoc->mRNAArray()->features()};

            if ($mRNA) {

                my $Feature = GffDoc::Feature::Gff3::exon->new(
                    biotype => 'protein_coding',
                    end     => $mRNA->end(),
                    id      => 'ID='.$mRNA->id().'-E;Parent='.$mRNA->id().';',
                    phase   => '.', 
                    seqid   => $mRNA->seqid(), 
                    source  => 'GffDoc::Convert::UtrExon',
                    start   => $mRNA->start(), 
                    strand => $mRNA->strand(), 
                    #parent => $from_thing->parent(),
                );

                $Feature->_parentsparent($mRNA->parent());

                $gffdoc->exonArray->add($Feature);

            } else { 
                Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', error => 'cannot find mRNA' );
            }
        #r need to iterate through the ordered list
        } else {

            my ($mRNA) = grep { $_->id() eq $k } @{$gffdoc->mRNAArray()->features()};

            if ($mRNA) {

                
                for my $i (0..$#ordered_list) {

                    my $cds = $ordered_list[$i];

#                    my $exon_start;
#                    my $exon_end;
#                    if ($i == 0) {
#                        $exon_start = $mRNA->start();
#                        $exon_end = $cds->[1];
#                    }
#                    elsif ($i == $ordered_list) {
#                        $exon_start = $cds->[0];
#                        $exon_end = $mRNA->end;
#                    }
#                    else {
#                        $exon_start = $cds->[0];
#                        $exon_end = $cds->[1];
#                    }

                    my $exon_start = $i == 0 ? $mRNA->start() : $cds->[0];
                    my $exon_end = $i == $#ordered_list ? $mRNA->end() : $cds->[1];

                    my $Feature = GffDoc::Feature::Gff3::exon->new(
                        biotype => 'protein_coding',
                        end     => $exon_end,
                        id      => 'ID='.$mRNA->id().'-E;Parent='.$mRNA->id().';',
                        phase   => '.', 
                        seqid   => $mRNA->seqid(), 
                        source  => 'GffDoc::Convert::UtrExon',
                        start   => $exon_start, 
                        strand => $mRNA->strand(), 
                        #parent => $from_thing->parent(),
                    );

                    $Feature->_parentsparent($mRNA->parent());

                    $gffdoc->exonArray->add($Feature);

                } 
                
            } else { 
                Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', error => 'cannot find mRNA' );
            }

        }

    }

    $log4->info('Generated ', $gffdoc->exonArray()->length(), ' features within the exonArray');

}

__PACKAGE__->meta->make_immutable();

1;


