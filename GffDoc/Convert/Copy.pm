package GffDoc::Convert::Copy;
use GffDoc::Exceptions qw/moan moan_e $nt/;
use Data::Dumper;
use Mouse;
with 'GffDoc::Activity::Role';

#/ make it also copy mRNA to genes... should be pretty much the same except need to modify id for gene and mRNA...

sub run {

    my ($self, $copy) = @_;
    # $conversion||=q{};
    # my ($from,$to) = split ':', $conversion;

    my $gffdoc = $self->gffdoc();
    my $log4 = $self->log4();


    my $from = $copy eq 'exon' ? 'exon' : 'CDS';
    my $to = $copy eq 'exon' ? 'CDS' : 'exon';
    my $method_from = $from.'Array';
    my $method_to = $to.'Array';
    my $class = 'GffDoc::Feature::Gff3::'.$to;
    my $missing = 0;

    $missing = 1 if (!$gffdoc->$method_from);
    $missing = 1 if (!$gffdoc->$method_from->features());
#    $missing_c = 1 if (!$gffdoc->CDSArray());
#    $missing_c = 1 if (!$gffdoc->CDSArray()->features());
#    my $excess = 0;
#    my $missing_c = 0;
#    my $meta = $gffdoc->meta();
#    for my $attr ($meta->get_attribute_list()) {
#        next if (!$gffdoc->$attr());
#        next if (!$gffdoc->$attr()->features());
#        next if ($attr eq 'CDSArray' || $attr eq 'exonArray');
#        $excess = 1;
#    }

    if ($copy ne 'exon' && $copy ne 'CDS') {
        Exception::GffDoc->throw( stage => 'Convert', type => 'IncorrectUsage', 
          error => 'This module simple copies exons->CDS or CDS->exons and requires exon or CDS as an argument'
        );
    } elsif ($missing) {
        Exception::GffDoc->throw( stage => 'Convert', type => 'IncorrectUsage', 
          error => 'There must be '.$from.' in order to copy them'
        );
    }

    $log4->info('There are ', $gffdoc->$method_from()->length, ' features within the ', $method_from);

    for my $from_thing (@{$gffdoc->$method_from->features()}) {
        

        my $Feature = $class->new(
            biotype => 'protein_coding',
            end     => $from_thing->end(),
            id      => 'ID='.$from_thing->id().';Parent='.$from_thing->parent().';',
            phase   => '.', 
            seqid   => $from_thing->seqid(), 
            source  => 'GffDoc::Convert::Copy',
            start   => $from_thing->start(), 
            strand  => $from_thing->strand(), 
            #parent  => $from_thing->parent(),
        );

        $Feature->_parentsparent($from_thing->parentsparent()) if $from_thing->parentsparent();

        $gffdoc->$method_to->add($Feature);

    }

    $log4->info('Generated ', $gffdoc->$method_to->length(), ' features within the '.$method_to);

}

__PACKAGE__->meta->make_immutable();

1;
