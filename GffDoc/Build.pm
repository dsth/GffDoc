package GffDoc::Build;
use Mouse;
use GffDoc::Exceptions qw/moan moan_e $nt moan_ref/;
use GffDoc::print;
with 'GffDoc::Activity::Role';
has 'polypeptide_build' => (is => 'ro');
has 'ignore_ref_check'  => (is => 'rw', isa => 'Int');
my $nt = qq{\n\t};

#y/ as yet this does not grab names or seq edits from polypetides with polypeptide option - the option is really only
#y/ there for the polypeptide_build option atm.

sub run {

    my $self = shift;

    my $gffdoc      = $self->gffdoc();
    my $gffdoc2     = GffDoc->new(GffGenes => GffDoc::Array::Gff::Gene->new());
    my $log4        = $self->log4();
    my $prob_genes  = $self->prob_genes();

    my $c = 0;
    my $nochildren_exception = 0;

#o/
eval {
#o/
    Exception::GffDoc->throw( stage => 'Build', type => 'NoGeneArray',error => 'There appear to be no genes to build! - disabled preparser incorrectly?')
      if (!$gffdoc->geneArray()->features());
    Exception::GffDoc->throw( stage => 'Build', type => 'NoGenes',error => 'There appear to be no genes to build! - disabled preparser incorrectly?')
      if (scalar @{$gffdoc->geneArray()->features()} == 0);
#o/
};
#o/
my $e;
#y cannot put other conditions in catch without breaking it?!? thus embed them...
if ( $e = Exception::Class->caught('Exception') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string ) if !$nochildren_exception;
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing shoulbe be caught at this level.', $@ );
    exit;
}
    my $genes = scalar @{$gffdoc->geneArray()->features()};

    GENE: 
    while (my $gene = $gffdoc->geneArray()->shift()) {

#o/
eval {
#o/

        if ($c == 0 || $c % 5 == 0) {  
            my $p = sprintf(q{%.2f},100*$c/$genes);
            #GffDoc::print::print_clean(9);
            GffDoc::print::print_number($p.'%');
            GffDoc::print::print_clean(6);
        }
        $c++;

        my $kiddymRNAs = $gffdoc->mRNAArray->findParentsDestructive($gene->id());

        #y this is not done in findParents so as to allow for pseudogenes
        Exception::GffDoc::Gene->throw( stage => 'Build', type => 'FeatureLine::Dangling::NoChildren',
          error => 'we do not tolerate orphan genes - i.e. genes without associated mRNAs '.$nt 
          .'perhaps the gene ID= and mRNA Parent= identifiers are not consistent?!? gene id='.$gene->id()
          . '. See Build log.', id => $gene->id(),
        ) if (scalar @{$kiddymRNAs} == 0);

        #y careful, circumventing constructor and therefore type controls etc.
        my $gffGene = bless $gene, 'GffDoc::Gff3::Gene';

        $log4->info('Reconstructing gene: '.$gene->id());

        for my $kiddymRNA (@{$kiddymRNAs}) {

            my $gffmRNA = bless $kiddymRNA, 'GffDoc::Gff3::mRNA';

            $gffmRNA->_gffGene($gffGene); # two-way mapping

            my $kiddyname = $kiddymRNA->id();
            my $kiddyexons = $gffdoc->exonArray->findParentsDestructive($kiddyname);
            my $kiddyCDSs = $gffdoc->CDSArray->findParentsDestructive($kiddyname);

            # we allow for no CDS for pseudogenes - till the genes are constructed we don't know which is which
            #/ have to make this biotype dependent?!?
            Exception::GffDoc::Gene->throw( stage => 'Build', type => 'FeatureLine::Dangling::NoChildren',
              error => 'we do not tolerate mRNAs without exons perhaps the mRNA ID= and '
              .'exon Parent= identifiers are not consistent?!? mRNA id='.$kiddyname
              . '. See Build log.', $gene->id(),
            ) if (scalar @{$kiddyexons} == 0);
            
            #/ build the mRNA
            #print qq{\nexons: }, scalar @{$kiddyexons};
            #print qq{\nCDSs: }, scalar @{$kiddyCDSs};
            $gffmRNA->_exons($kiddyexons);
            $gffmRNA->_CDSs($kiddyCDSs);
            $gffGene->addmRNA($gffmRNA);
        }

        $gffdoc2->GffGenes()->add($gffGene);
#o/
};
#o/
my $e;
#y cannot put other conditions in catch without breaking it?!? thus embed them...
if ( $e = Exception::Class->caught('Exception::GffDoc::Gene') ) { 
    $log4->error($e->error);
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string ) if !$nochildren_exception;
    $nochildren_exception = 1;
    push @{$prob_genes}, $e->id;
    next GENE;
} elsif ( $e = Exception::Class->caught('Exception::GffDoc') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing shoulbe be caught at this level.', $@ );
    exit;
}

    }

    # print qq{\nleft over genes: }, scalar @{$gffdoc->GeneArray()->features()}; #y this is clearly not an option!?!

#/ horrible but can't really be bothered to clean it up - as most have been checked for in RefCheck these aren't recoverable - i.e. have you been screwing around?!?

#o/ 
eval {
    Exception::GffDoc->throw( 
        error => 'Feature mRNAArray is empty!',
        stage => 'Build',
        type  => 'EmptyFeatureArray',
    ) if (!$gffdoc->mRNAArray()->features());
    
    Exception::GffDoc::mRNA->throw( 
        error => 'mRNA feature without '
          .'corresponding gene features - perhaps ID= and Parent= identifiers aren not consistent?! See Build log.',
        stage => 'Build',
        type  => 'FeatureLine::Dangling::mRNA::NoParent',
        list  => [map {$_->parent() } @{$gffdoc->mRNAArray()->features()}],
    ) if (scalar @{$gffdoc->mRNAArray()->features()} > 0);

    Exception::GffDoc->throw( 
        error => 'Feature exonArray is empty!',
        stage => 'Build',
        type  => 'EmptyFeatureArray',
    ) if (!$gffdoc->exonArray()->features());
    
    Exception::GffDoc::exon->throw( 
        error => 'exon feature without '
          .'corresponding mRNA features - perhaps ID= and Parent= identifiers aren not consistent?! See Build log.',
        stage => 'Build',
        type  => 'FeatureLine::Dangling::exon::NoParent',
        list  => [map {$_->parent() } @{$gffdoc->exonArray()->features()}],
    ) if (scalar @{$gffdoc->exonArray()->features()} > 0);

    Exception::GffDoc->throw( 
        error => 'Feature array CDSArray is empty!',
        stage => 'Build',
        type  => 'EmptyFeatureArray',
    ) if (!$gffdoc->CDSArray()->features());
    
    Exception::GffDoc::CDS->throw( 
        error => 'CDS feature without '
          .'corresponding mRNA features - perhaps ID= and Parent= identifiers aren not consistent?! See Build log.',
        stage => 'Build',
        type  => 'FeatureLine::Dangling::CDS::NoParent',
        list  => [map {$_->parent() } @{$gffdoc->CDSArray()->features()}],
    ) if (scalar @{$gffdoc->CDSArray()->features()} > 0);
};
if ( $e = Exception::Class->caught('Exception::GffDoc::mRNA') ) { 
    moan_ref ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string,
      scalar @{$e->list}, @{$e->list}[0..3]
    );
    $log4->error('mRNA '.$_->id().' without corresponding gene: '.$_->parent()) for (@{$gffdoc->mRNAArray()->features()});
    #$log4->error('mRNA without corresponding gene: '.$_) for (map { $_->id() } @{$gffdoc->mRNAArray()->features()});
    #/ make it possible to disable this at their own risk - why are undefined values getting past?!?
    #/ i.e. during standard gff3 - when using external module don't disable
    exit if (!$self->ignore_ref_check());
} elsif ( $e = Exception::Class->caught('Exception::GffDoc::exon') ) { 
    moan_ref ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string,
      scalar @{$e->list}, @{$e->list}[0..3]
    );
    $log4->error('exon '.$_->id().' without corresponding mRNA: '.$_->parent()) for (@{$gffdoc->exonArray()->features()});
    exit if (!$self->ignore_ref_check());
} elsif ( $e = Exception::Class->caught('Exception::GffDoc::CDS') ) { 
    moan_ref ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string,
      scalar @{$e->list}, @{$e->list}[0..3]
    );
    $log4->error('CDS '.$_->id().' without corresponding mRNA: '.$_->parent()) for (@{$gffdoc->CDSArray()->features()});
    exit if (!$self->ignore_ref_check());
} elsif ( $e = Exception::Class->caught('Exception::GffDoc') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing shoulbe be caught at this level.', $@ );
    exit;
}

    #y carry over the polypeptides if using them for construction of genes with dodgy gff format
    $gffdoc2->_polypeptideArray($gffdoc->polypeptideArray) if ($self->polypeptide_build);

    return $gffdoc2;

}

__PACKAGE__->meta->make_immutable();

1;

