package GffDoc::Exceptions;
@ISA = qw(Exporter);
@EXPORT_OK = qw(moan moan_e moan_c moan_ref $nt);
our $nt = qq{\n}.' 'x9;
our $ntt = qq{\n}.' 'x11;

my $line = qq{\n==================================================\n};

#y this is carnage. loads of sub-types with no differential handling - its a horribly inefficient form of message logging...
#y hence rarely is anything actually caught at low-level namespaces...
use Exception::Class (
    'Exception',

    'Exception::GffDoc'                     => { isa => 'Exception',
        fields => ['type','stage'],
    },

    'Exception::GffDoc::NonUniq'                     => { isa => 'Exception', },

    'Exception::GffDoc::FeatureConstructor' => { isa => 'Exception::GffDoc',
        fields => ['type','stage','hash'], # don't need to re-def type/stage...
    },

    'Exception::GffDoc::ReferentialIntegrity' => { isa => 'Exception::GffDoc',
        fields => ['type','stage','list'], # don't need to re-def type/stage...
    },

    'Exception::GffDoc::StopCodons'         => { isa => 'Exception::GffDoc' },

    'Exception::GffDoc::mRNA'               => { isa => 'Exception::GffDoc::ReferentialIntegrity' },
    'Exception::GffDoc::exon'               => { isa => 'Exception::GffDoc::ReferentialIntegrity' },
    'Exception::GffDoc::CDS'                => { isa => 'Exception::GffDoc::ReferentialIntegrity' },

    'Exception::GffDoc::Gene'               => { isa => 'Exception::GffDoc',
        fields => ['type','stage', 'id'],
    },


    'Exception::GffDoc::Gene::Strand'       => { isa => 'Exception::GffDoc::Gene' },
    'Exception::GffDoc::Gene::SeqRegion'    => { isa => 'Exception::GffDoc::Gene' },

    'Exception::GffDoc::Gene::MissingSlice' => { isa => 'Exception::GffDoc::Gene' },
);

sub moan { 
    #/ make traces an option
    #print $line; 
    #print qq{\n$_} for (@_);
    #print $line; }
    #my $thing = $_[2];
    $_[2] =~ s/\n(\S)/\n\n* $1/g;
    $_[2] =~ s/called/\n> called/g;
    #$_[2] =~ s/^(.*line.*line\s\d+)\n/$1/s;
    #$thing =~ s/^(.*?line.*?line)/$1/s;
    #$thing =~ s/^(.*\n.*\n).*/$1/m;
    #$_[2] =~ s/(.*\n.*)\n/$1/m;
    printf( "\n"
      . "-------------------- EXCEPTION --------------------\n"
      . "Class:  %s\n"
      . "Msg:    %s\n"
      . "Trace:  %s"
      . "---------------------------------------------------\n",
    @_);
}


sub moan_e { 
    my ($c,$s,$ty,$e,$p,$f,$l,$tr,$n,$xxxxxx) = @_;
    chomp $tr;
    $tr =~ s/\n(\S)/\n\n* $1/g;
    $tr =~ s/called/\n> called/g;
#moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
#    printf( "\n"
#      . "-------------------- EXCEPTION --------------------\n"
#      . "Class:  %s\n"
#      . "Stage:  %s\n"
#      . "Type:   %s\n"
#      . "yyyy:   %s\n"
#      . "Msg:    %s\n"
#      . "Trace:  %s"
#      . "---------------------------------------------------\n",
#    $c, $s, $ty, $e, $n, $tr);
    my $msg =<<HMM;
-------------------- EXCEPTION --------------------
Class:   $c
Stage:   $s
Package: $p
Type:    $ty
Msg:     $e
File:    $f
Line:    $l
Trace:   $tr
---------------------------------------------------
HMM

    print $msg;

}

sub moan_c { 
    my ($c,$s,$ty,$e,$p,$f,$l,$tr,$feature) = @_;
    chomp $tr;
    $tr =~ s/\n(\S)/\n\n* $1/g;
    $tr =~ s/called/\n> called/g;
    my $feat_str = q{};
    while (my ($k, $v) = each %{$feature}) { $feat_str .= qq{${ntt}$k => $v}; }
    #$feat_str .= qq{\n};

    my $msg =<<HMM;
-------------------- EXCEPTION --------------------
Class:   $c
Stage:   $s
Package: $p
Type:    $ty
Msg:     $e
File:    $f
Line:    $l
Feature: $feat_str
---------------------------------------------------
HMM

    print $msg;

}

sub moan_ref { 
    my ($c,$s,$ty,$e,$p,$f,$l,$tr,$num,@examples) = @_;
    chomp $tr;
    $tr =~ s/\n(\S)/\n\n* $1/g;
    $tr =~ s/called/\n> called/g;
    @examples = grep { defined } @examples;
    @examples[$_] .= ', ' for (0..$#examples-1);
    my $msg =<<HMM;
-------------------- EXCEPTION --------------------
Class:   $c
Stage:   $s
Package: $p
Type:    $ty
Msg:     $e
Exmpl(s) There is/are $num instances incl.: @examples - See log file.
File:    $f
Line:    $l
---------------------------------------------------
HMM

    print $msg;

}

1;
