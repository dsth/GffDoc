#/ this has been done with no regard to efficiency at all. there are some fairly simple changes that can be made to improve speed 
#r/ either way the process should be faster as it gets furhter along - i.e. shorter lists left...

#={{{1 use statements
use strict;
use warnings;
use GffDoc::Exceptions qw/moan moan_e moan_ref $nt/;
use GffDoc::Types; # not needed 
use GffDoc::Parser;
use GffDoc::Build;
use GffDoc::Check;
use GffDoc::RefCheck;
use GffDoc::eStore;
use DBI;
use Log::Log4perl qw(get_logger :levels);
use Getopt::Long;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use version; our $VERSION = qv('0.0.1');
#  class_type 'DateTime';
#=}}}


#/ exceptions are generally severly abused. they are either uncaught or caught at high-level
#/ - i.e. the sub-types are generally there to make messages clearer not for differential handling

#b/ options/log4/dbad setup... ########################################={{{1

my @time = localtime;
my @mon = qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/;
$time[5] -= 100;
my $tstamp = $time[3].$mon[$time[4]].$time[5];

my $cs_name = q{};
my $cs_version = q{};
my $db = 'blah_db'; 
my $file = 'blah.gff';
my $gtf = 0;
my $host = 'localhost';
my $leafnonunique = 0;
my $ignore_id_duplicates = 0;
my $logic_name = 'ensemblgenomes';
my $non_coding_cds = 0;
my $non_protein_coding;
my $non_protein_coding_types = q{};
my $polypeptide = 0;
my $polypeptide_build = 0;
my $port = '3027';
my $pw = q{};
my $user = 'dsth';
my $validate = 0;
my $version = '1';
my @biotype_resolutions;
my @biotype_conversions;
my $mRNA_biotype_dominant = 0;
my $gene_biotype_dominant = 0;
my @modules = ();
my @types;
my $pseudo_string = q{};
my $no_grep_file = 0;
my $ignore_ref_check1 = 0;
my $ignore_ref_check2 = 0;
my $biotype_key = q{};

&GetOptions( 
    #'M=s@'              => \@modules,
    #'types=s@'          => \@types,
    'M=s'                           => \@modules,
    'biotype_convert=s'                        => \@biotype_conversions,
    'biotype_key=s'                 => \$biotype_key,
    'biotype_resolution=s'          => \@biotype_resolutions,
    'coordsystem:s'                 => \$cs_name,
    'coordsystem_version:s'         => \$cs_name,
    'dbname:s'                      => \$db,
    'file:s'                        => \$file,
    'gene_biotype_dominant'         => \$gene_biotype_dominant,
    'gtf_names'                     => \$gtf,
    'host:s'                        => \$host,
    'ignore_ref_check1'              => \$ignore_ref_check1,
    'ignore_ref_check2'              => \$ignore_ref_check2,
    'ignore_id_duplicates'              => \$ignore_id_duplicates,
    'leafnonunique'                 => \$leafnonunique,
    'logicname=s'                   => \$logic_name,
    'mRNA_biotype_dominant'         => \$mRNA_biotype_dominant,
    'nogrep'                        => \$no_grep_file,
    'non_coding_cds'                => \$non_coding_cds,
    'non_protein_coding'            => \$non_protein_coding,
    'non_protein_coding_types:s'    => \$non_protein_coding_types,
    'password:s'                    => \$pw,
    'polypeptide'                   => \$polypeptide,
    'polypeptide_build'             => \$polypeptide_build,
    'port:i'                        => \$port,
    'pseudo_string:s'                      => \$pseudo_string,
    'type=s'                        => \@types,
    'user:s'                        => \$user,
    'validate'                      => \$validate,
    'version:s'                     => \$version,
);

#y put this into lower-level namespace to stop messages propagating up log chain... can't be bothered to read documnetation atm
my $log4 = get_logger('GffDoc::Run');
$log4->level($DEBUG);
my $appender1 = Log::Log4perl::Appender->new('Log::Dispatch::Screen');
my $layout = Log::Log4perl::Layout::PatternLayout->new("%p %d %M > %m%n");
#my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %p> %F{1}:%L %M - %m%n");
$appender1->layout($layout);
$log4->add_appender($appender1);

#g do this now so that we don't do all the earlier stuff to have it die on submission...
my %types;
my %biotype_resolutions;
my %biotype_conversions;
my ($permitted_types, $db_ad);

eval {

    $db_ad = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -dbname     => $db, 
        -host       => $host, 
        -pass       => $pw,
        -port       => $port, 
        -user       => $user,
    ) || Exception::GffDoc->throw( stage => 'PreReqs', type => 'eDbAdaptor', error => 'Could not connect to ensembl database '.$db );

    for my $type (@types) {
        my ($t, $a) = split '=', $type;
        $a =~ /^(ignore|gene(\x3a\w+)?|pseudogene|mRNA(\x3a\w+)?|exon|CDS)$/ || Exception::GffDoc->throw( stage => 'Options', type => 'IlegalTypeReclassification',
            error => 'The type reclassifier must be one of ignore, gene, pseudogene, mRNA, exon or CDS: '.$t.'='.$a 
        );
        $types{$t} = $a;
    }

    #y add start_codon etc.?!?
    $permitted_types = 'gene(\x3a\w+)?|mRNA(\x3a\w+)?|CDS|exon|pseudogene|intron'
    .'|polyA_sequence|polyA_site|five_prime_UTR|three_prime_UTR'
    .'|start_codon|stop_codon';

    #y using peptides for build process clearly needs them to be parsed to turn on this option
    $polypeptide = 1 if $polypeptide_build;

    my $other_types = q{};
    if (scalar (keys %types) > 0) {
        my %keys;
        @keys{(keys %types)} = ();
        $other_types = join '|', (keys %keys);
        $permitted_types .= qq{$permitted_types|${other_types}};
        #$permitted_types .= qq{$permitted_types|${other_types}${polypeptide}};
    }

    #$polypeptide = $polypeptide ? '|polypeptide' : q{}; #y parse polypeptides?
    $permitted_types .= '|polypeptide' if $polypeptide; #y parse polypeptides?
    $permitted_types = qr{$permitted_types};

    if ($non_protein_coding_types && $non_protein_coding_types !~ /^(\w+)(,\w+)*$/) { 
        Exception::GffDoc->throw( stage => 'Options', type => 'IlegalNonProteinCodingTypesString',
          error => 'The string passed to the -non_protein_coding_types should be a list of comma-separated types' 
        );
    } else {
        $non_protein_coding_types =~ s/,/|/g;
        $non_protein_coding_types = '|'.$non_protein_coding_types;
    }

    for my $btc (@biotype_conversions) {

        if ($btc !~ /^(\w+:\w+)$/) { 
            Exception::GffDoc->throw( stage => 'Options', type => 'IlegalBiotypeResolutionString',
            error => 'The string passed to the -biotype_conflict should be of form old_biotype:new_biotype' 
            );
        } else {

            my ($old, $new) = split ':', $btc;
            $biotype_conversions{$old} = $new;
            
        }
    }
    for my $bt_r (@biotype_resolutions) {

        if ($bt_r !~ /^(\w+,\w+)(,\w+)*$/) { 
            Exception::GffDoc->throw( stage => 'Options', type => 'IlegalBiotypeResolutionString',
            error => 'The string passed to the -biotype_resolution should be a list of comma-separated biotypes with principal type first' 
            );
        } else {

            my @bt_rs = split ',', $bt_r;
            my $dom_bt = $bt_rs[0];
            @bt_rs = sort { $a cmp $b } @bt_rs;
            my $str = join ':', @bt_rs;
            $biotype_resolutions{$str} = $dom_bt;
            
        }
    }
};
my $e; #y plain lazy
if ( $e = Exception::Class->caught('Exception') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing shoulbe be caught at this level.', $@ );
    exit;
}

#=}}}

#b/ PREPARSE/PARSE ####################################################={{{1

$log4->info('Creating preparser object');

my $preparser = GffDoc::Preparser->new(
    gff3file        => $file,
    leafnonunique   => $leafnonunique,
    ignore_id_duplicates   => $ignore_id_duplicates,
    modules         => scalar @modules,
    polypeptide     => $polypeptide,
    types_regexp    => $permitted_types,
);

$log4->info('Preparsing '.$file);
eval { 
    $preparser->run() || Exception::GffDoc->throw( stage => 'PreParser', type => 'GffFile::Format', error => 'Address PreParser issues!' ); 
    #$parser->_preparseFile() || Exception::GffDoc->throw( stage => 'PreParser', type => 'GffFile::Format', error => 'Address PreParser issues!' ); 
};
if ( $e = Exception::Class->caught('Exception') ) { #y here till can be bothered to fill in details of preparser
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    moan ( 'This is a Bug. Nothing shoulbe be caught at this level.', $@ );
    exit;
}

$log4->info('Creating parser object');

#y late and lazy addition
$GffDoc::Feature::Role::additional_types_string = $non_protein_coding_types if ($non_protein_coding_types);

my $parser = GffDoc::Parser->new(
    gff3file        => $file,
    gtf             => $gtf,
    leafnonunique   => $leafnonunique,
    modules         => scalar @modules,
    polypeptide     => $polypeptide,
    pseudo_string          => $pseudo_string,
    types           => \%types,
    types_regexp    => $permitted_types,
    gff_feature_num => $preparser->gff_feature_num(),
    biotype_key     => $biotype_key,
    biotype_conversions => \%biotype_conversions,
);

$log4->info('Parsing '.$file);

my $gffdoc = $parser->run();

#=}}}

#b/ FETCH/SLICES ######################################################={{{1

#g do this now rather than go through the hassle of reconstruction and having it throw later over this
#g have a disable option if running just for checking file e.g. -noslice

# up till now we catch here - i.e. allow exceptions to propagate. from here, the steps are too slow to not make 
# this really irritating. catch exceptions from now on within the modules themselves.

my $build_log = ${file}.'_'.${tstamp}.'_gff.log.txt';
my $log4_details = get_logger('GffDoc::Build');
$log4_details->level($DEBUG);
my $appender_details = Log::Log4perl::Appender->new(
    "Log::Dispatch::File",
    filename => $build_log,
    #mode     => "append",
);
$appender_details->layout($layout);
$log4_details->add_appender($appender_details);

my $eslices;
if (!$validate) {
    eval {
        $log4->info('Fetching e! slices');
        $eslices = GffDoc::eSlices->new(
            cs_name => $cs_name,
            cs_version => $cs_version,
            dbad => $db_ad, 
            gffdoc => $gffdoc,
            log4 => $log4, 
        );
        $eslices = $eslices->run();
    };
    if ( $e = Exception::Class->caught('Exception::GffDoc') ) { 
        moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
        exit;
    }
} else {
    $log4->info('Skipping e! slice fetch (validation mode)');
}
#=}}}

#y put a regexp for 'single' ':' and use to pass args - thus GTF module uses args for whether or not to clone CDS to exons too
#y keep just -gtf_names option - should allow gff names with GTF with some changes?!?

#b/ CONVERT/MANIPULATE... #############################################={{{1

#my @modules = qw/ GffDoc::Convert::Dummy /;
eval {
    for my $module (@modules) { # -M option

        my $arg_str = q{};
        if ($module =~ /(.*?[^:])(:([^:].*))/) { $module = $1; $arg_str = $3 }
        my @args = split ',', $arg_str if $arg_str;
        
        $log4->info('Running additional module '.$module);

        # require qq{${convMod}}; # need to add '.' path
        my $file = $module;
        $file =~ s/::/\//g;
         
        my $p = $main::INC{'GffDoc/Parser.pm'};
        $p =~ s/GffDoc\/Parser\.pm$//;

        Exception::GffDoc->throw( stage => 'Convert', type => 'ModuleNotFound', error => 'Cannot find module '.$file )
          unless (-e $p.$file.'.pm');

        #Exception::GffDoc->throw( stage => 'Convert', type => 'Module::Error', error => 'Cannot find module '.$file )

        require $file.'.pm'; #b forces path '/'
        my $converter = $module->new(
            gffdoc => $gffdoc,
            log4 => $log4, 
        );

        $converter->run(@args);
    }
};
if ( $e = Exception::Class->caught('Exception') ) { #y here till can be bothered to fill in details of preparser
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
} elsif ( $e = Exception::Class->caught() ) { 
    print  "\n-------------------- EXCEPTION --------------------\n"
      .'External Module Problem'.qq{\n}.$@
      ."---------------------------------------------------\n";
    exit;
}
#=}}}

#r/################# nothing till now is caught/recoverable #############

#b/ REFCHECK (also generates info for e! storage) ########################={{{1

$log4->info('Quick check of referential integrity');

my $refcheck = GffDoc::RefCheck->new(
    gffdoc              => $gffdoc,
    log4                => $log4, 
    non_protein_coding  => $non_protein_coding,
);

eval { $refcheck->run(); };
if ( $e = Exception::Class->caught('Exception::GffDoc::ReferentialIntegrity') ) { 
    moan_ref ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, 
      $e->trace->as_string, scalar @{$e->list}, @{$e->list}[0..3] 
    );
    my $err = $e->error;
    $err =~ s/\s?This is not legal.*$//s;
    $log4_details->error($err.' e.g. feature '.$_) for (@{$e->list});
    exit if (!$ignore_ref_check1);
} elsif ( $e = Exception::Class->caught('Exception::GffDoc::NonUniq') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit if (!$ignore_ref_check1);
} elsif ( $e = Exception::Class->caught('Exception') ) { 
    moan_e ( ref $e, $e->stage, $e->type, $e->error, $e->package, $e->file, $e->line, $e->trace->as_string );
    exit;
}
#=}}}

#b/ RECONSTRUCT #######################################################={{{1

my $prob_genes = [];

#y now we have what should be proper gff we re-build the genes - i.e. total separation of parsing and building
$log4->info('Logging messages to: '.$build_log);
$log4->info('Rebuilding genes and running further consistency checks - this may take some time');
my $builder = GffDoc::Build->new(
    gffdoc      => $gffdoc,
    log4        => $log4_details, 
    prob_genes  => $prob_genes,
    ignore_ref_check => $ignore_ref_check2,
);

$gffdoc = $builder->run();
#=}}}

#b/ CHECK (also generates info for e! storage) ########################={{{1

$log4->info('Checking genes');

my $check = GffDoc::Check->new(
    cs_name     => $cs_name,
    gffdoc      => $gffdoc,
    log4        => $log4_details, 
    logicname   => $logic_name,
    prob_genes  => $prob_genes,
    biotype_resolutions => \%biotype_resolutions,
    mRNA_biotype_dominant => $mRNA_biotype_dominant,
    gene_biotype_dominant => $gene_biotype_dominant,
);

$gffdoc = $check->run();
#=}}}

#b/ STORE #############################################################={{{1

#y at the moment not bothering to skip problem genes here - prolly should (i.e. next on existence of name in prob_genes
my $out;
if (!$validate) {
    $log4->info('Storing genes');
    my $store = GffDoc::eStore->new(
        cs_name => $cs_name,
        dbad => $db_ad, 
        gffdoc => $gffdoc,
        log4 => $log4_details, 
        logicname => $logic_name,
        non_coding_cds => $non_coding_cds,
        # prob_genes  => $prob_genes,
        slices => $eslices,
        version =>  $version, 
    );

    $out = $gffdoc = $store->run();
} else {
    $log4->info('Skipping e! store (validation mode)');
    $out = [0,0,0];
}

#y whenever skipping print $gfh $gene->id().qq{\n};
my $grep_file = ${file}.'_'.${tstamp}.'_grep.txt';
if (!$no_grep_file) { #y make an option
    open my $gfh, '>', $grep_file;
    print $gfh $_.qq{\n} for (@{$prob_genes});
}

$log4->info('Saved '.$out->[0].' genes.');
my $msg = $out->[1] ? ' - See log!' : q{.};
$log4->info('Skipped '.$out->[1].' genes'.$msg);
$log4->info('Recategorised '.$out->[2].' genes and transcripts to non_coding_cds biotype due to presence of stop codons in translation.');
#=}}}
