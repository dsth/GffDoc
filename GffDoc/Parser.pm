package GffDoc::Preparser;
use Mouse;
use GffDoc::Types;
use GffDoc::Exceptions qw/moan moan_e $nt moan_c/;
with 'GffDoc::Activity::Minimal::Role';
has [qw/ignore_id_duplicates leafnonunique gff_feature_num polypeptide/] => (is => 'rw');
has [qw/gff3file/]  => (is => 'ro', required => 1,);
has 'types_regexp' => (is => 'ro');
has 'modules' => (is => 'rw', isa => 'Int');

sub run {

    my $self = shift;
    
    my $infile = $self->gff3file();

    open (my $GFF_hndl, '<', $infile) || Exception::GffDoc->throw( stage => 'PreReqs', type => 'GffFile', error => 'Cannot open '.$infile.' file for reading: '.$!);

    my $c = 0;
    my $gff = 0;
    my $gtf = 0;
    my $nine_col_format = 0;
    my $no_final_semicolon = 0;
    my $unrecognised_types = 0;
    my %ids;
    my %type_set;

    LINE:
    while (my $line = <$GFF_hndl>) {

        chomp $line;

        next LINE if ($line =~ /^[\x20\x09]*$/x); # ignore empty lines
        next LINE if ($line =~ /^[\x20\x09]*\x23/x); # ignore comment lines
        
        $no_final_semicolon = 1 if ($line !~ /\x3B\x20*$/); # check lines end with ';' 

        #y pointless considering the following...
        $nine_col_format = 1 if ($line !~/^([^\t]+\t){8}[^\t]+$/);

        if ($line =~/^([^\t]+\t){2}([^\t]+)\t([^\t]+\t){5}([^\t]+)$/) {
            my $lt = $2;
            $type_set{$lt} = 1;
            my $comment = $4; #y each new match resets $2/$4...
            $gff = 1 if ($comment =~ /(ID|Parent)=\S+?\s?\;/);
            $gtf = 1 if ($comment =~ /(gene|transcript)_id\s+?"\S+?"\s?\;/);

            if ($comment =~ /ID=(\S+?)\s?\;/) {

                #g never disable for higher level structures
                unless(($self->leafnonunique() && ($lt eq 'exon' || $lt eq 'CDS')) || $self->ignore_id_duplicates()) {

                    if (exists $ids{$1}) {
                        Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::NonUniqID',
                            error => 'Feature IDs must be unique: '.$1.' at line '.$..'.'
                            .$nt.'This requirement can be disabled for exons and CDS with -leafnonunique option'
                            .$nt.'Alternatively, if this concerns types irrelevant to actual parsing ifor preparser to ignore with -ignore_id_duplicates'
                        );
                    } else {
                        $ids{$1} = 1;
                    }    
                }
            }
        }    

        $c++;
        if ($c % 1000 == 0) {
            #GffDoc::print::print_number(sprintf("%d",++$c));
            GffDoc::print::print_number(sprintf("%d",$c));
            GffDoc::print::print_clean(length($c));
        }

    }

    $self->gff_feature_num($c);
    
    return 1 if ($self->modules()); #y only v. partial preparse with modules option

    my $permitted_types = $self->types_regexp();

    my @not_permitted = grep { $_ !~ /^$permitted_types$/ } (keys %type_set);

    #/b separate these blocks to make it easier to read...

    #b block anything that doesn't adhere to the accepted types or is overtly ignored
    if (scalar @not_permitted == 1 && $not_permitted[0] eq 'polypeptide') {
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::Polypeptides', 
          error => 'To parse polypeptide features you must use the -polypeptide option'
        ) unless $self->polypeptide();
    }
    elsif (scalar @not_permitted > 0) {
        $unrecognised_types = 1;
    #b so we now have just gene/pseudogene, mRNA, exon, CDS - and any combination of them - don't care about pseudogenes for these purposes
    }
    
    #y we want these warnings whether or not we throw below so deal with this crap here in horrible fashion
    unless (($gtf || $gff) && !$no_final_semicolon && !$nine_col_format && !$unrecognised_types) { 
        my $n = qq{\n};
        print $n
          .'=============== PREPARSER PROBLEMS ================'.$n.$n;
        print &_thereDocNoSemiColon($infile) if $no_final_semicolon;
        print &_thereDocNameErrors if (!$gtf && !$gff);
        print &_thereDoc9Col if $nine_col_format;
        print &_thereDocUnrecognisedTypes(@not_permitted) if $unrecognised_types;
        print '==================================================='.$n;
    }
    
    if ( 
      exists $type_set{'gene'} 
      && exists $type_set{'mRNA'} 
      && exists $type_set{'exon'} 
      && exists $type_set{'CDS'} 
    ) { 
        #f gene, mRNA, exon, CDS
        #y completely normal gff3... 
    } elsif (
      !exists $type_set{'gene'} 
      && !exists $type_set{'mRNA'}  ###### force use of options?!? ######
      && !exists $type_set{'exon'} 
      && exists $type_set{'CDS'} 
    ) { 
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::GtfStandard', 
          error => 'This looks like basic Gtf (CDS only) run GTF module in standard mode with "-M GffDoc::Convert::Gtf:s"'
          .$nt.'If the gene/transcript identifiers adhere to gtf conventions you will also need the -gtf option'
        );
    } elsif (
      !exists $type_set{'gene'} 
      && !exists $type_set{'mRNA'} 
      && exists $type_set{'exon'} 
      && exists $type_set{'CDS'} 
    ) {  
        #f exon, CDS
        #f so with extended option -M GffDoc::Convert::GTF:e - should really use exons to build genes?!?
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::GtfExtended', 
          error => 'This looks like extended Gtf (CDS & exons) run GTF module in extended mode with "-M GffDoc::Convert::Gtf:e"'
          .$nt.'If the gene/transcript identifiers adhere to gtf conventions you will also need the -gtf option'
        );
    } elsif (
      !exists $type_set{'gene'} 
      && !exists $type_set{'mRNA'} 
      && exists $type_set{'exon'} 
      && !exists $type_set{'CDS'} 
    ) { 
        #f exon - 
        #y extended sort of Gtf - i.e. has non-coding regions too
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::GtfTypeError', 
          error => 'Looks like these exons are misnamed CDS and this is basic GTF format'.$nt
            .'Either overtly change the type of the exons to CDS or tell parser to process them as CDS with -exon=CDS'.$nt
            .'You will also need to use the GTF module in standard mode with "-M GffDoc::Convert::Gtf:s"'
        );
    } elsif (
      !exists $type_set{'gene'} 
      && exists $type_set{'mRNA'} 
      && exists $type_set{'exon'} 
      && !exists $type_set{'CDS'} 
    ) { 
        #f mRNA, exon
        #y extended sort of Gtf - i.e. has non-coding regions too
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::PlaceFiller', 
          error => 'Looks like these exons are misnamed CDS and this is a mess - rename exons to CDS...'
        );
    } elsif (
      !exists $type_set{'gene'} 
      && exists $type_set{'mRNA'} 
      && !exists $type_set{'exon'} 
      && exists $type_set{'CDS'} 
    ) { 
        #f mRNA, CDS
        #y extended sort of Gtf - i.e. has non-coding regions too
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::PlaceFiller', 
          error => 'Silly hybrid - easily dealt with using as yet unimplemented external module with -M option'
        );
    } elsif (
      !exists $type_set{'gene'} 
      && exists $type_set{'mRNA'} 
      && exists $type_set{'exon'} 
      && exists $type_set{'CDS'} 
    ) { 
        #f mRNA, exon, CDS
        #y can't recall where, but have seen this. either bin the mRNAs and process as an extended Gtf
        #y or use the appropriate module - See POD
        print &_thereDocTempPlaceFiller;
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::PlaceFiller', 
          error => 'Silly hybrid - easily dealt with using as yet unimplemented external module with -M option'
        );
    } elsif (
      exists $type_set{'gene'} 
      && exists $type_set{'mRNA'} 
      && !exists $type_set{'exon'} 
      && exists $type_set{'CDS'} 
    ) { 
        #f gene, mRNA, CDS
        print &_thereDocMissingExons($infile);
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::MissingExons', 
          error => 'Silly hybrid - easily dealt with using as yet unimplemented external module with -M option'
        );

    } elsif (
      exists $type_set{'gene'} 
      && exists $type_set{'mRNA'} 
      && exists $type_set{'exon'} 
      && !exists $type_set{'CDS'} 
    ) { 
        #f gene, mRNA, exon
        print &_thereDocNotValid($infile);
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::IllegalCombination', 
          error => 'this combination of features (gene, mRNA, exon) without CDS is not a valid gene model format - Read the above message!'
        );
    #y not really necessary, but...
    } else { # what genes & mRNA?!?...
        my $types_str = q{};
        $types_str .= qq{${nt} > $_} for (keys %type_set);
        Exception::GffDoc->throw( stage => 'PreParser', type => 'Types::UnrecognisedCombination', 
          error => 'This file has a combination of feature types seems to be meaningless.'.$nt
          .'Gff should have at least types gene, mRNA, CDS, exon and Gtf should have CDS (and possibly exon).'.$nt
          .'This file has:'.$types_str.$nt.' if this is a simple naming error use -type optin to re-classify e.g. -type bad_name=proper_name'
        );
    }

    if (($gtf || $gff) && !$no_final_semicolon && !$nine_col_format && !$unrecognised_types) { 
        return 1;
    } else {
        return 0;
    }
}

sub _thereDocNoSemiColon {

    my $file = shift;
my $message = <<"PREFIX";
[*] Gtf/Gff features lines must terminate with ';'. For a quick fix try:

  > perl -i -pe 's/(.*[^;]\\s*)\\n/\${1};\\n/' $file

PREFIX

return $message;

}

sub _thereDocNameErrors {

my $message = <<"PREFIX";
[*] Gff/Gtf features must adhere to proper naming conventions e.g. 

    GFF3: contigX   VB    mRNA  ...   ID=mRNA00001;Parent=gene00001;

    GTF:  ContigY   SNAP  CDS   ...   gene_id "001"; transcript_id "001.1";

see: http://www.sequenceontology.org/gff3.shtml and http://mblab.wustl.edu/GTF22.html 
respectively.

PREFIX

return $message;

}

sub _thereDoc9Col {

my $message = <<"PREFIX";
[*] Gff/Gtf features must adhere 9 column format e.g. 

    GFF3: <seqname> <source> <feature> <start> <end> <score> <strand> <phase> [attributes]

    GTF:  <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

In which fields are separated by the ASCI horizontal tab character (0x09). Make sure all 9 columns are
present and that there are no extra tab characters concealed. See: http://www.sequenceontology.org/gff3.shtml 
and http://mblab.wustl.edu/GTF22.html respectively.

PREFIX

return $message;

}

sub _thereDocUnrecognisedTypes {

    my @types = @_;
    my $n = qq{\n};

my $message = '[*] Unrecognised types! Do not recognise the following type(s):'.$n.$n;

    @types = grep { $_ ne 'polypeptide' } @types;
    my $sep = ' 'x8;
    $message .= $sep.q{* '}.$_.q{'}.$n for (@types);

$message .= <<"PREFIX";

    Only gene|pseudogene, mRNA, exon and CDS types are actually used by the 
    parser. Beyond these 5 types: introns, polyA_sequence, polyA_site, 
    five_prime_UTR, three_prime_UTR are tolerated, but completely ignored. 
    No 'other' types are permitted within a gff3 file unless you explicitly 
    hand these types as options to the parser i.e. 

        -type <daftly_named_feature>=<what_to_do_with_it>

    <what_to_do_with_it> may take the values of: ignore, gene, mRNA,exon, CDS.

    thus if there are imaginitively named utr sequences that should be ignored 
    by the parser use:

        -type annoying_alternative_name_for_utr=ignore

    or if someone has decided to call CDSs something else, tell the parser 
    about this with:

        -funky_name_for_CDS=CDS

PREFIX

return $message;

}

sub _thereDocNotValid {

    my $file = shift;
my $message = <<"PREFIX";
===================================================
INVALID GFF FORMAT: Read the message.
===================================================

This file seems to have genes, mRNA and exons features without the presence of
CDS features. This is clearly not valid Gff. There are a few likely explanations:

[1] There is no protein coding information in this file - i.e. all models are 
    pseudogenes or the models are incomplete. Either of these options seems
    unlikely. However, if its the former you can turn off preparsing by 
    running in external module mode with a dummy module (add: 
    -M GffDoc::Convert::Dummy) to commandline. Alternatively, if it is the 
    latter, then the gff file models are junk gene model descriptions.

[2] The coding sequence boundaries are encoded elsewhere e.g. polypeptide features
    give the offset of the coding start relative to the mRNAs - this is the case 
    with sources such as AphidBase - you'll need to run the parser with the 
    appropriate external module to fix this using the -M option.

[3] The exons are 'misnamed' CDS and there are no UTRs - that is to say the exons 
    are identical to the CDS. This is fixed easily by cloning exon features of a
    file (be careful not to clone exons that are part of pseudogenes as this 
    will either promote them to protein_coding if the pseudogene status is not
    marked-up correctly, or will cause exceptions to be thrown due to the presence
    of CDS for pseudogenes (not as yet permitted). To clone all exons as CDS to 
    make valid gff you can try something simple like:

        > perl -i -pe 'my \$s = \$_; if(s/(([^\\t]+\\t){2})exon/\${1}CDS/) { print \$s }' $file

    Alternatively, you can inform the parser that these are incorrectly name CDS with
    
        -type exon=CDS
       
    and run using the appropriate external module with -M option.

[4] The exons are 'misnamed' CDS and there are UTRs. Clearly, all intron/exon 
    organisation of non-coding regions is lost (is it really worth keeping the gene/mRNA 
    info and not processing it as if it were Gtf?). If you wish to process them with the 
    non-coding regions, a 'fudge' is to assume that all non-coding regions are contained 
    within the terminal 5' and 3' exons. This is clearly similar to situation [2]. Again,
    you can inform the parser that these are incorrectly named CDS with
    
        -type exon=CDS

    and use the appropriate external module to fix this using the -M option.

PREFIX

return $message;

}

sub _thereDocMissingExons {

    my $file = shift;
my $message = <<"PREFIX";
===================================================
Missing exons
===================================================

This file seems to be lacking exon features. This means that either:

[1] There are no UTR sequences and CDS and exons are identical in which case you 
    can either copy the the CDS to make exons with something like:

        > perl -i -pe 'my \$s = \$_; if(s/(([^\\t]+\\t){2})CDS/\${1}exon/) { print \$s }' $file

    To make valid gff3 or you can load the file using the appropriate external
    module with the -M option - if you're not certain this is the case use the
    same module as in [2].

[2] There are UTR sequences and  all intron/exon organisation of non-coding regions is 
    lost (is it really worth keeping the gene/mRNA info and not just processing it as 
    if it were Gtf?). If you wish to process them with the non-coding regions, a 'fudge' 
    is to assume that all non-coding regions are contained within the terminal 5' and 3' 
    exons and use the appropriate external module to fix this with the -M option.

PREFIX

return $message;

}

sub _thereDocTempPlaceFiller {

    my $file = shift;
my $message = <<"PREFIX";

This combination requires the use of additional external modules. See the POD.

PREFIX

return $message;

}

__PACKAGE__->meta->make_immutable();

1;

package GffDoc::Parser;
use GffDoc::Exceptions qw/moan moan_e $nt moan_c/;
use GffDoc::print;
use GffDoc::Types;
use Mouse;
use Carp;
use List::MoreUtils qw(any);
with 'GffDoc::Activity::Minimal::Role';

#y gtf has start and stop codons - prolly ought to ignore them and modify start/stop appropriately...
has [qw/gtf pseudo_string leafnonunique gff_feature_num polypeptide/] => (is => 'rw');
#has 'gff_gene_num' => (is => 'rw'); 
has [qw/gff3file/]  => (is => 'ro', required => 1,);
has 'types' => (is => 'rw', isa => 'HashRef', auto_deref => 1);
has 'types_regexp' => (is => 'ro');
has 'modules' => (is => 'rw', isa => 'Int');
has 'biotype_key' => (is => 'rw', isa => 'Str');
has 'biotype_conversions' => (is => 'rw', isa => 'HashRef');

sub run {

    my $self = shift;

    my $infile = $self->gff3file();
    my $biotype_key = $self->biotype_key() ? $self->biotype_key() : q{};
    my $biotype_conversions = $self->biotype_conversions() if $self->biotype_conversions();

    open (my $GFF_hndl, '<', $infile) || Exception::GffFile->throw( 
      error => 'Cannot open '.$infile.' file for reading: '.$!
    );

    my $gffdoc = GffDoc->new();
    
    #/ pseudo genes get classified as genes with biotype type pseudogene here.

    #b/ for gtf have option that generates those features and then also reconstructs - i.e. brings in -M automatically

    my @types = qw/gene mRNA exon CDS/;
    push @types, 'polypeptide' if $self->polypeptide();
    #my $polypeptide = $self->polypeptide() ? 'polypeptide' : q{};

    for my $type (@types) {
    
        #y instantiate each class array
        my $classFeat = qq{GffDoc::Feature::Gff3::${type}};
        my $classArray = qq{GffDoc::Feature::Array::${type}};

        my $FeatArray = $classArray->new();
        my $method = qq{_${type}Array};
        $gffdoc->$method($FeatArray);
        
    }

#    #y/ could just hard code it... modify the actual feature class if necessary
#    GffDoc::Feature::Gff3::CDS->meta->add_attribute(
#        parentsparent => (is => 'ro', isa => 'Str', writer => '_parentsparent')
#    ) if $self->gtf();


#o/
eval {
#o/

    my $c = 0;
    my $total = $self->gff_feature_num();

    #/ just explode exons with multiple parents - i.e. create a separate entry for each one 
    #/ reference counting would require more effort to do safely than i can be bothered with
    LINE:
    while (my $line = <$GFF_hndl>) {

        $c++;
        if ($c % 5000 == 0) {
            my $p = sprintf(q{%.1f},100*$c/$total);
            GffDoc::print::print_number($p.'%');
            GffDoc::print::print_clean(5);
        }

        chomp $line;

        next LINE if ($line =~ /^[\x20\x09]*$/x); # ignore empty lines
        next LINE if ($line =~ /^[\x20\x09]*\x23/x); # ignore comment lines

        Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Format', 
            error => qq{Malformed feature line $. (lines must terminate with ';') [run preparser]}
        ) if ($line !~ /\x3B\x20*$/);

        if ( 

          my ($landmark, $source, $type, $start, $end, $score, $strand, $phase, $attribs) = 
          ($line =~/^\x20*([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t(.+)\x20*$/x)) {

            #y/ change types to the acceptable ones and/or biotypes...
            $type = exists $self->types()->{$type} ? $self->types()->{$type} : $type;

            next LINE if ($type =~ /^(ignore|intron|polyA_sequence|polyA_site|five_prime_UTR|three_prime_UTR|start_codon|stop_codon)$/);

            my $permitted_types = $self->types_regexp();

            Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::UnrecognisedType: '.$type, error => 'Run the preparser (turn off the -M option)' )
              if ($type !~ /^$permitted_types$/);

            #r default biotype - can allow for any biotype as only protein coding get CDS etc.
            my $biotype = 'protein_coding';

            #/ v. basic atm., with type pseudogene made a gene with biotype pseudogene
            #/ alternatively you can supply a regexp to apply to the attribs column to do this
            my $pseudo = $self->pseudo_string();
            if ($type eq 'pseudogene' || ($pseudo && $attribs =~ /$pseudo/)) {
                $type = 'gene';
                $biotype = 'pseudogene';
            #y if it's gene/mRNA we check for biotype and apply
            } elsif ($type =~ /^(gene|mRNA)\x3a(\w+)$/) { 
                $type = $1;
                $biotype = $2;
            } elsif ($biotype_key && $attribs =~ /${biotype_key}\s?=\s?(\S+?);/) {
                $biotype = $1;
            }

            $biotype = exists $biotype_conversions->{$biotype} ? $biotype_conversions->{$biotype} : $biotype;

            #y/ at this point they must have been forced into one of the 4/5 streams properly so just in case
            Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Type', error => 'Doing something funny?!? type='.$type )
              if ($type !~ /^(gene|mRNA|exon|CDS|polypeptide)$/);

            #r/ explode exons/CDS - for simplicity make independent objects and not to leave it to reference counting?!?
            my $parent_str = q{};

#y/ combine these and put in an option for polypeptides.

#/ make it handle two ways multiple parents may be declared

            my @parents = ($attribs =~ /Parent=\s?(\S+?)\s?;/g);
 
            if ((scalar @parents == 1) && ($parents[0] =~ /,/)) {
                @parents = split(q{,}, $parents[0]);              
            } elsif (scalar @parents == 1){
                $parent_str = $parents[0]; # too lazy to change the vars below
            } else {
                $parent_str = q{};
            }

            #if ($attribs =~ /Parent=\s?(\S+?)\s?;/) { $parent_str = $1; }

            if ((scalar @parents > 1) && ($type eq 'exon' || $type eq 'CDS')) {
            #if (($parent_str =~ /,/) && ($type eq 'exon' || $type eq 'CDS')) {

                #my @parents = split(q{,}, $parent_str);

                for my $parent (@parents) {

                    my $Feature = qq{GffDoc::Feature::Gff3::${type}}->new(
                        biotype => $biotype,
                        end     => $end,
                        gtf     => $self->gtf(),
                        id      => $attribs,
                        phase   => $phase, 
                        seqid   => $landmark, 
                        source  => $source,
                        start   => $start, 
                        strand  => $strand, 
                        type    => $type,
                    );

                    #y there are multiple parents so we override the constructor parent recognition from id/attribs
                    $Feature->_parent($parent);

                    #y use var as method name
                    my $method = $type.'Array';
                    $gffdoc->$method->add($Feature);
                
                }

            } else {

                #y use var as class name
                my $Feature = qq{GffDoc::Feature::Gff3::${type}}->new(
                    biotype => $biotype,
                    end     => $end,
                    gtf     => $self->gtf(),
                    id      => $attribs,
                    phase   => $phase, 
                    seqid   => $landmark, 
                    source  => $source,
                    start   => $start, 
                    strand  => $strand, 
                    type    => $type,
                );
                
                #$Feature->_id($id) if $id;
                #$Feature->_parent($parent_str) if $parent_str;
                #$Feature->_attribs($attribs) if $attribs;

                #y use var as method name
                my $method = $type.'Array';
                $gffdoc->$method->add($Feature);

                my $t = qq{GffDoc::Feature::Gff3::${type}};

                #$t->meta->add_attribute(name =>(is => 'ro', isa => 'Str'));
                #$Feature->meta->add_attribute(name =>(is => 'ro', isa => 'Str'));

            }

        } else { 

            Exception::GffDoc->throw( stage => 'Parser', type => 'FeatureLine::Format', error => qq{Malformed feature line $. (lines must adhere to nine-column, tab-delimited format)});
            #f throw incorrectly formatted line exception
        }    

    } # end of line processing using diamond operator

#o/
};
my $e;
# my ($c,$s,$ty,$e,$p,$f,$l,$tr,$feature) = @_;
if ( $e = Exception::Class->caught('Exception::GffDoc::FeatureConstructor') ) {
    moan_c ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string, $e->hash );
    
    #while (my ($k, $v) = each %{$e->hash()}) { print qq{\n$k => $v}; }

    exit;
} elsif ( $e = Exception::Class->caught('Exception') ) {
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
#y some exceptions in types.pm are in the wrong namespace atm....
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing should be caught at this level.', $@ );
    exit;
}
#o/

    close $GFF_hndl;

    return $gffdoc;
}

__PACKAGE__->meta->make_immutable();

1;

