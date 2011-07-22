package GffDoc::print;
use strict;
use warnings;

my $l;

sub print_number {
    my $str = shift;
    #print $str;
    system qq{echo -n $str};
    $l = length($str);
    return;
}

sub print_clean {
    my $r = shift;
    print qq{\x08} x $l;
    return;
}


1;

