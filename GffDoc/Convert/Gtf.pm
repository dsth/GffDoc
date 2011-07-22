package GffDoc::Convert::Gtf;
use GffDoc::Exceptions qw/moan moan_e $nt/;
use Data::Dumper;
use Mouse;
with 'GffDoc::Activity::Role';

sub run {

    my ($self, $mode, @args) = @_;
    #print qq{\n$_} for (@args); die;
    $mode||=q{};

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

    if ($mode eq 'e') {
        Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
          error => 'There should be both exon and CDS features AND no more to run this conversion module'
            .qq{\n\t}.'ANY EXONS LACKING CDS WILL BE MADE PSEUDOGENES'
        ) if ( $excess || $missing_c || $missing_e );
    } elsif ($mode eq 's') {
        Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
          error => 'There should be only CDS features AND no more to run this conversion module'
        ) if ( $excess || $missing_c || !$missing_e );
    } else {
        Exception::GffDoc->throw( stage => 'Convert', type => 'WrongModule', 
          error => 'Need to run this module in either standard -M GffDoc::Convert::Gtf:s or extended '
          .'mode -M GffDoc::Convert::Gtf:e'
        );
    }

    my %mRNA_coords;
    my %stuff;

    #/ quick hack cos i'm not well to allow it to handle extended gtf with utrs
    my $method = $mode eq 'e' ? 'exonArray' : 'CDSArray';

    my $CDS_num = $gffdoc->$method->length; 

    $log4->info('Processing '.$method.' features (I of II)');
    my $c = 0;

    #y this shouldn't be destructive - duh!?: while (my $exon = $gffdoc->exonArray()->shift()) or regular findParents
    for my $CDS (@{$gffdoc->$method->features()}) {
    # for my $CDS (@{$gffdoc->CDSArray()->features()}) {

        if ($c == 0 || $c % 50 == 0) {  
        my $p = sprintf(q{%.2f},100*$c/$CDS_num);
            GffDoc::print::print_number($p.'%');
            GffDoc::print::print_clean(6);
        }
        GffDoc::print::print_clean();
        $c++;

        #/ are there CDS with the same parent? - i.e. is it protein_coding
        my $biotype = @{$gffdoc->CDSArray->findParentsNonDestructive($CDS->parent())} == 0 ? 'pseudogene' : 'protein_coding';

        #y check once...
        $stuff{$CDS->parent()} = [$CDS->parentsparent(),$biotype] if (!exists $stuff{$CDS->parent()});

        if (!exists $mRNA_coords{$CDS->parent()}) {
            $mRNA_coords{$CDS->parent()} = [$CDS->start(),$CDS->end(),$CDS->seqid(),$CDS->strand()];
            #$mRNA_coords{$CDS->parent()} = [$CDS->start(),$CDS->end(),$CDS->seqid(),$CDS->strand(),$biotype];
        }
        elsif ($CDS->start() < $mRNA_coords{$CDS->parent()}->[0]) {
            $mRNA_coords{$CDS->parent()}->[0] = $CDS->start();
        }
        elsif ($CDS->end() > $mRNA_coords{$CDS->parent()}->[1]) {
            $mRNA_coords{$CDS->parent()}->[1] = $CDS->end();
        }

        if ($mode eq 's') {
            my $exon = $CDS->copy2Type('exon');
            $gffdoc->exonArray->add($exon);
        }

    }

#    while (my ($k, $v) = each %biotypes) { print qq{\n$k => $v}; }

    #/ this uses a fudge. genes are just abstract units and the api actually ignores the coords you give it in
    #/ recalculates start/end based on most extreme positions of constituent mRNAs so don't bother doing it here...

    my $mRNA_num = scalar (keys %mRNA_coords); 
    $log4->info('Generating gene/mRNA features (II of II)');
    $c = 0;

    my %set;

    while (my ($mRNAid,$v) = each %mRNA_coords) {

        if ($c == 0 || $c % 50 == 0) {  
        my $p = sprintf(q{%.2f},100*$c/$mRNA_num);
            GffDoc::print::print_number($p.'%');
            GffDoc::print::print_clean(6);
        }
        $c++;

        #y we don't care about phase as Check.pm generates them

        if (!exists $set{$stuff{$mRNAid}->[0]}) { #b make a gene first time

            my $geneFeature = GffDoc::Feature::Gff3::gene->new(
                biotype => $stuff{$mRNAid}->[1],
                end     => $v->[1],
                id      => 'ID='.$stuff{$mRNAid}->[0].';',
                phase   => '.', 
                seqid   => $v->[2], 
                source  => 'GffDoc::Convert::Broad',
                start   => $v->[0], 
                strand  => $v->[3], 
            );

            $gffdoc->geneArray->add($geneFeature);
        }

        #y using exists not defined so how gives a proverbial...
        $set{$stuff{$mRNAid}->[0]} = 0; #b duh, this was not the gene name (was biotype)...
        
        my $mRNAFeature = GffDoc::Feature::Gff3::mRNA->new(
            biotype => $stuff{$mRNAid}->[1],
            end     => $v->[1],
            id      => 'ID='.$mRNAid.';Parent='.$stuff{$mRNAid}->[0].';',
            phase   => '.', 
            seqid   => $v->[2], 
            source  => 'GffDoc::Convert::Broad',
            start   => $v->[0], 
            strand  => $v->[3], 
        );

        $gffdoc->mRNAArray->add($mRNAFeature);

    }

    $log4->info('Generated ', $gffdoc->mRNAArray->length(), ' features within the mRNAArray');
    $log4->info('Generated ', $gffdoc->geneArray()->length(), ' features within the geneArray');
    $log4->info('Generated ', $gffdoc->exonArray()->length(), ' features within the exonArray') if ($mode eq 's');

    #use Data::Dumper; print Dumper @{$gffdoc->mRNAArray->features}[0];
    #use Data::Dumper; print Dumper @{$gffdoc->geneArray->features}[0]; 

}

__PACKAGE__->meta->make_immutable();

1;
