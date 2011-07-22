package GffDoc::RefCheck;
use Data::Dumper;
use GffDoc::Exceptions qw/$nt/;
use Mouse;
with 'GffDoc::Activity::Role';
has 'non_protein_coding' => (is => 'rw');

#y as modules may be mixed-n-matched many things are re-checked here and re-checked later (i.e. parts may have bypassed some module checks..
sub run {


    my $self = shift;
    my $gffdoc = $self->gffdoc();
    my $log4 = $self->log4();

    my $meta = $gffdoc->meta();
   
    #y complete overkill but will shut up the hmmm...

    my $gene_num = $gffdoc->geneArray()->length();
    my %gene_ids = %{$gffdoc->geneArray()->distinctAsHashRef('id')}; #/ undef values so use exists...
    my $mRNA_num = $gffdoc->mRNAArray()->length();
    my %mRNA_ids = %{$gffdoc->mRNAArray()->distinctAsHashRef('id')};
    my $exon_num = $gffdoc->exonArray()->length();
    my $CDS_num = $gffdoc->CDSArray()->length();

    $log4->info('There are '.$gene_num.' genes in the geneArray');
    $log4->info('There are '.scalar (keys %gene_ids).' distinct gene ids in the geneArray');
    $log4->info('There are '.$mRNA_num.' mRNAs in the mRNAArray');
    $log4->info('There are '.scalar (keys %mRNA_ids).' distinct mRNA ids in the mRNAArray');

    my %mRNA_parents = %{$gffdoc->mRNAArray()->distinctAsHashRef('parent')};
    my %exon_parents = %{$gffdoc->exonArray()->distinctAsHashRef('parent')};
    my %CDS_parents = %{$gffdoc->CDSArray()->distinctAsHashRef('parent')};

    $log4->info('There are '.$exon_num.' exons in the exonArray');
    $log4->info('There are '.$CDS_num.' CDSs in the CDSArray');
    $log4->info('There are '.scalar (keys %mRNA_parents).' distinct parent genes refered to in the mRNAArray');
    $log4->info('There are '.scalar (keys %exon_parents).' distinct parent mRNAs refered to in the exonArray');
    $log4->info('There are '.scalar (keys %CDS_parents).' distinct parent mRNAs refered to in the CDSArray');

    #g turn off with -pseudogene option

    Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Check', type => 'Referential::Missing::CDS', 
      error => 'There are mRNA(s) referenced by exon that are not referenced by CDS. This is not legal for protein_coding genes'
      .$nt.'You must explicitly permit non-protein coding genes with the -non_protein_coding option',
      list => [(grep { !exists $CDS_parents{$_} } (keys %exon_parents))],
    ) if (!$self->non_protein_coding && (scalar (grep { !exists $CDS_parents{$_} } (keys %exon_parents)) > 0 )); #y scalar > 0 just for clarity

    Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Check', type => 'Referential::Missing::exon', 
      error => 'There are mRNA(s) referenced by CDS that are not referenced by exons. This is not legal Gff',
      list => [(grep { !exists $exon_parents{$_} } (keys %CDS_parents))],
    ) if (scalar (grep { !exists $exon_parents{$_} } (keys %CDS_parents)) > 0 ); #y scalar > 0 just for clarity

    my %piece_parents = (%exon_parents, %CDS_parents);
    $log4->info('There are '.scalar (keys %piece_parents).' distinct parent mRNA refered to by exonArray and CDSArray');

    Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Check', type => 'Referential::Missing::mRNA', 
      error => 'There are mRNA(s) refered to by exon/CDS that do not exist as explicit feature entities. This is not legal Gff',
      list => [(grep { !exists $mRNA_ids{$_} } (keys %piece_parents))],
    ) if (scalar (grep { !exists $mRNA_ids{$_} } (keys %piece_parents)) > 0 ); 

    Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Check', type => 'Referential::Missing::mRNA::Components', 
      error => 'There are orphan mRNA(s) features that are not refered to by any exon/CDS. This is not legal Gff',
      list => [(grep { !exists $piece_parents{$_} } (keys %mRNA_ids))],
    ) if (scalar (grep { !exists $piece_parents{$_} } (keys %mRNA_ids)) > 0 ); 

    Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Check', type => 'Referential::Missing::gene', 
      error => 'There are gene(s) refered to by mRNA that do not exist as explicit feature entities. This is not legal Gff',
      list => [(grep { !exists $gene_ids{$_} } (keys %mRNA_parents))],
    ) if (scalar (grep { !exists $gene_ids{$_} } (keys %mRNA_parents)) > 0 );

    #/ clearly already thrown if the exon/CDS are declared and just have mRNA missing
    Exception::GffDoc::ReferentialIntegrity->throw( stage => 'Check', type => 'Referential::Missing::gene::Components', 
      error => 'There are orphan gene(s) features that are not refered to by any mRNA. This is not legal Gff',
      list => [(grep { !exists $mRNA_parents{$_} } (keys %gene_ids))],
    ) if (scalar (grep { !exists $mRNA_parents{$_} } (keys %gene_ids)) > 0 ); 

    Exception::GffDoc::NonUniq->throw( stage => 'Check', type => 'Referential::gene::NonUniqID',
      error => 'There are genes with non-unique ids'
    ) if ($gene_num > (keys %gene_ids));

    Exception::GffDoc::NonUniq->throw( stage => 'Check', type => 'Referential::mRNA::NonUniqID',
      error => 'There are mRNAs with non-unique ids'
    ) if ($mRNA_num > (keys %mRNA_ids));

    return;
    
}

sub repeats {

    my $list_ref = shift;

    my %uniq;
    my %non_uniq;

    for my $id (@{$list_ref}) {

        if (exists $uniq{$id}) { 
            #push @non_unique, $id; 
            $non_uniq{$id} = 1;
        } else { 
            $uniq{$id} = 1;
        }

    }

    return [(keys %non_uniq)];

}


__PACKAGE__->meta->make_immutable();

1;

=head2

genes need a method to return ALL ids
mRNA needs a method to return all ids and all parents
exons/cds need a method to return ALL parents (perhaps also parentsparents?!?)

- get extras etc., BEFORE - running reconstruction - this throws an uncaught exception
unless you run with a specific option to force it to continue...

my @simpsons=("homer","bart","marge","maggie","lisa");
my @females=("lisa","marge","maggie","maude");
my %simpsons=map{$_ =>1} @simpsons;
my %females=map{$_=>1} @females;
# the intersection of @females and @simpsons: my @female_simpsons = grep( $simpsons{$_}, @females );
# proof it works print "Female Simpson:\t$_\n" foreach (@female_simpsons);
# the difference of @females and @simpsons my @male_simpsons=grep(!defined $females{$_}, @simpsons);
# proof it works print "Male Simpson:\t$_\n" foreach (@male_simpsons);
my %union = ();
# the union of @females and @simpsons foreach(@females,@simpsons){ $union{$_}=1; } my @union2 = keys %union; i
# or just do this # my @union = (@females, @simpsons);

=cut


