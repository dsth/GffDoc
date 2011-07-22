package GffDoc::Convert::Dummy;
use Data::Dumper;
use Mouse;
with 'GffDoc::Activity::Role';

sub run {

    my $self = shift;
    #print Dumper $self->gffdoc();
    #die; #/ not catching this higher up
    return;
    
}

__PACKAGE__->meta->make_immutable();

1;



