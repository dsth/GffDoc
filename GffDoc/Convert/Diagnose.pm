package GffDoc::Convert::Diagnose;
use Data::Dumper;
use Mouse;
with 'GffDoc::Activity::Role';

sub run {

    my $self = shift;

    my $gffdoc = $self->gffdoc();

    my $meta = $gffdoc->meta();

    for my $attr ($meta->get_attribute_list()) {

        #next if ($attr eq 'GffGenes' || $attr eq 'polypeptideArray');
        if (!$gffdoc->$attr()) {
            print qq{\nThere is no array of type $attr};
            next;
        }
        if (!$gffdoc->$attr()->features()) {
            print qq{\nThere is an array of type: $attr but it has no features};
            next;
        }

        print qq{\nThere are }, scalar @{$gffdoc->$attr()->features()}, ' features within the ', $attr;
    }

    print qq{\n};

    #/ this is where we put all the referential integrity checks if desired?!?

    Exception::Diagnostic->throw( error => 'done checking' );
    
}

__PACKAGE__->meta->make_immutable();

1;



