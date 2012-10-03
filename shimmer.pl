#!/usr/bin/perl -w
#########################################################
# Author:	Nancy F. Hansen
# Program:	"SHiMMer"
# Function:	shimmer
#               Program to use hypothesis testing within
#               the context of an HMM to predict regions
#               of copy number alteration in tumor
#               sequence data (by comparing to a matched
#               normal sample's sequence.)  SHiMMer also
#               predicts single-nucleotide somatic
#               mutations using Fisher's Exact Test with
#               multiple testing correction.
##########################################################

use strict;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw /tempfile tempdir /;

use vars qw($VERSION);

our %Opt;
our $R_EXE = 'R';
our $ANNOVAR_EXE = 'annotate_variation.pl';
our $INSTALL_DIR = '/home/nhansen/projects/gtb';

our @SOM_COLORS = ( {'power' => 0.0, 'color' => '255,0,0'}, # red
                    {'power' => 0.4, 'color' => '255,100,0'}, # orange
                    {'power' => 0.6, 'color' => '255,255,0'}, # yellow
                    {'power' => 0.8, 'color' => '0,255,0'}, # green
                    {'power' => 0.9, 'color' => '0,0,255'} ); # blue

$| = 1; # let's see output right away!

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr(q$Rev:10000$, 4)/10000;

my $program_name = $0;
if ($program_name !~ /^\//) { # add path
    my $path = `pwd`;
    chomp $path;
    $program_name = "$path/$program_name";
}
my $short_program_name = $program_name;
$short_program_name =~ s:.*/::; # without path

my $Usage = "Usage: $short_program_name <--region chr1:1000-2000> <--ref reference fasta> <bam file from normal sample> <bam file from mutated sample>\nFor more information, type \"perldoc $short_program_name\".";

my $printcounts_exe = "$INSTALL_DIR/c/printCompCounts/printCompCounts";
my $baumwelch_exe = "$INSTALL_DIR/c/umdhmm-v1.02-gauss/sh.run_esthmm";
my $viterbi_exe = "$INSTALL_DIR/c/umdhmm-v1.02-gauss/sh.run_viterbi";
my $hmm_model = "$INSTALL_DIR/c/umdhmm-v1.02-gauss/default_ng_model.hmm";

process_commandline();

($#ARGV == 1)
    or die "$Usage";

my $bam1 = $ARGV[0];
my $bam2 = $ARGV[1];
my $ref_fasta = $Opt{'ref'};
my $region = $Opt{'region'};
my $minqual = $Opt{'minqual'}; # disregard any base with quality lower than minqual

if ($Opt{'counts'}) {

    # check files:
    foreach my $file ($ref_fasta, $bam1, $bam2) {
        if (!(-r $file)) {
            die "Can\'t read file $file!\n";
        }
    }
    my $som_file = $Opt{'som_file'};
    my $het_file = $Opt{'het_file'};

    # print counts of different alleles at all interesting positions
    print_counts($ref_fasta, $bam1, $bam2, $region, $minqual, $som_file, $het_file, $printcounts_exe);
    # test counts for significance with Fisher's exact test (without mult testing corr)
    test_counts($som_file, $het_file);
}
elsif ($Opt{'bh'}) {

    # apply Benjamini-Hochberg to determine cutoff for significance
    my $max_q = $Opt{'max_q'};
    my $test_fof = $Opt{'test_fof'};
    my $outfile = $Opt{'outfile'};
    my $vsfile = $Opt{'vs_file'};
    my $vcffile = $Opt{'vcf_file'};
    if (!(-r $test_fof)) {
        die "Cannot read file $test_fof to apply B-H correction!\n";
    }
    bh_correct_tests($test_fof, $max_q, $outfile);
    write_vs_file($outfile, $vsfile) if ($vsfile);
    write_vcf_file($outfile, $vcffile) if ($vcffile);
}
elsif ($Opt{'annotate'}) {
    my $annovar_db = $Opt{'annovardb'};
    my $vs_file = $Opt{'vs_file'};
    my $ann_vs_file = $Opt{'outfile'};
    my $buildver = $Opt{'buildver'};
    annotate_variants($annovar_db, $buildver, $vs_file, $ann_vs_file);
}
elsif ($Opt{'symbols'}) {

    my $input_file = $Opt{'input'}
        or die "Must specify file of counts with --input option when running $0 --symbols\n";

    my $output_file = $Opt{'outfile'}
        or die "Must specify output file with --outfile option when running $0 --symbols\n";

    my $outdir = $Opt{'outdir'}
        or die "Must specify output directory with --outdir option when running $0 --symbols\n";

    my $max_q = $Opt{'max_q'};

    my $ra_symbols = write_symbols($input_file, $output_file, $outdir, $max_q, $baumwelch_exe, $viterbi_exe, $hmm_model);

    my $plot = $Opt{'plots'};
    if ($plot) {
        plot_symbols($ra_symbols, $outdir);
    }
}
elsif ($Opt{'covg'}) {

    my $input_file = $Opt{'input'}
        or die "Must specify input file with --input option when running $0 -covg\n";

    my $outdir = $Opt{'outdir'}
        or die "Must specify output directory with --outdir option when running $0 --covg\n";

    calc_power_file($input_file, $outdir, $bam1, $bam2, $ref_fasta, $printcounts_exe);

}
else {
    run_shimmer($program_name, $ref_fasta, $bam1, $bam2, $region, $minqual);
}

## Subroutine to run all of the steps of SHiMMer

sub run_shimmer {

    my $shimmer = shift;
    my $ref_fasta = shift;
    my $bam1 = shift;
    my $bam2 = shift;
    my $region = shift;
    my $minqual = shift;

    if ($Opt{'outdir'} && !(-e $Opt{'outdir'})) {
        mkdir $Opt{'outdir'}
            or die "Couldn\'t create $Opt{'outdir'}: $!\n";
    }

    my $tmpdir = $Opt{'outdir'} || tempdir( "run_shimmer_XXXXXX", DIR => '.' );
    if ($tmpdir !~ /^\//) {
        my $pwd = `pwd`;
        chomp $pwd;
        $tmpdir = "$pwd/$tmpdir";
    }
   
    if ((!$Opt{'som_file'}) || (!$Opt{'het_file'})) { 
        my $command = "$shimmer --counts --ref $ref_fasta $bam1 $bam2 --min_som_reads $Opt{min_som_reads} --som_file $tmpdir/som_counts.txt --het_file $tmpdir/het_counts.txt";
        $command .= " --region $region" if ($region);
        $command .= " --minqual $minqual" if ($minqual);
        system($command) == 0
            or die "Failed to run $shimmer --counts!\n";
    }

    if (!$Opt{'skip_tests'}) {
        my $max_q = $Opt{'max_q'};
        open SOMFOF, ">$tmpdir/som_counts.fof"
            or die "Couldn\'t open $tmpdir/som_counts.fof for writing: $!\n";
        print SOMFOF "$tmpdir/som_counts.tests.txt\n";
        close SOMFOF;
        
        (system("$shimmer --max_q $max_q --test_fof $tmpdir/som_counts.fof --bh --vs_file $tmpdir/somatic_diffs.vs --vcf_file $tmpdir/somatic_diffs.vcf --outfile $tmpdir/som_counts.bh.txt $bam1 $bam2") == 0)
            or die "Failed to run $shimmer --test_fof $tmpdir/som_counts.fof --bh --vs_file $tmpdir/somatic_diffs.vs --outfile $tmpdir/som_counts.bh.txt $bam1 $bam2!\n";
    
        if ($Opt{'annovardb'}) {
            my $annovardb = $Opt{'annovardb'};
            my $buildver = $Opt{'buildver'};
            (system("$shimmer --annotate --buildver $buildver --annovardb $annovardb --annovar $ANNOVAR_EXE --vs_file $tmpdir/somatic_diffs.vs --outfile $tmpdir/somatic_diffs.ANN.vs $bam1 $bam2")==0)
                or die "Failed to run $shimmer --annotate --buildver $buildver --annovardb $annovardb --annovar $ANNOVAR_EXE --vs_file $tmpdir/somatic_diffs.vs --outfile $tmpdir/somatic_diffs.ANN.vs $bam1 $bam2: $!\n";
        }

        if (!$Opt{'skip_cna'}) {    
            open HETFOF, ">$tmpdir/het_counts.fof"
                or die "Couldn\'t open $tmpdir/het_counts.fof for writing: $!\n";
            print HETFOF "$tmpdir/het_counts.tests.txt\n";
            close HETFOF;
            
            (system("$shimmer --max_q 0 --test_fof $tmpdir/het_counts.fof --bh --outfile $tmpdir/het_counts.bh.txt $bam1 $bam2") == 0)
                or die "Failed to run $shimmer --test_fof $tmpdir/het_counts.fof!\n";
        
            my $run_string = "$shimmer --max_q 0.05 --symbols --input $tmpdir/het_counts.bh.txt --outfile $tmpdir/cnv.symbols.txt --outdir $tmpdir $bam1 $bam2";
            $run_string .= " --plots" if ($Opt{'plots'});

            system("$run_string") == 0
                or die "Couldn\'t run $run_string\n";
        }
    }

    if ($Opt{'power'}) {
        my $run_string = "$shimmer --input $tmpdir/som_counts.bh.txt --outdir $tmpdir --covg --ref $ref_fasta $bam1 $bam2";

        system("$run_string") == 0
            or die "Couldn\'t run $run_string\n";
    }
   
} ## end run_shimmer

sub print_counts {

    my $ref_fasta = shift;
    my $bam1 = shift;
    my $bam2 = shift;
    my $region = shift;
    my $minqual = shift;
    my $som_file = shift;
    my $het_file = shift;
    my $printcounts_exe = shift;

    my $insert = $Opt{'insert'};
    my $min_max_counts_file = $Opt{'min_max'};
    my $min_tumor_reads = $Opt{'min_som_reads'};

    my $ra_het_count_limits = read_min_max_file($min_max_counts_file); # store limits of different genotypes

    # call mpileup (via the c-script "printCompCounts"), and select sites that are independent, but have enough coverage/diversity to be informative
    my $printcounts_call = "$printcounts_exe -bam1 $bam1 -bam2 $bam2 -fasta $ref_fasta";
    $printcounts_call .= " -region $region" if ($region);
    $printcounts_call .= " -minqual $minqual" if ($minqual);
    
    print "Calling $printcounts_call\n";

    open COUNTS, "$printcounts_call | "
        or die "Couldn\'t execute $printcounts_call!\n";

    open SOM, ">$som_file"
        or die "Couldn\'t open $som_file for writing: $!\n";

    open CNV, ">$het_file"
        or die "Couldn\'t open $het_file for writing: $!\n";

    my ($cur_chr, $cur_pos, $cur_win_start, $cur_win_end, $best_geno, $best_cov, $best_string);

    while (<COUNTS>) {
        chomp;
        my ($chr, $pos, $ref_base, $base1, $normal1_count, $tumor1_count, $base2, $normal2_count, $tumor2_count) = split /\t/, $_;
        $ref_base = uc $ref_base;
        my $total_norm = $normal1_count + $normal2_count;
        my $total_tumor = $tumor1_count + $tumor2_count;
        my $total_alt = ($base2 eq $ref_base) ? $normal1_count + $tumor1_count : $normal2_count + $tumor2_count;
        my $normal_ref = ($base2 eq $ref_base) ? $normal2_count : $normal1_count;
        my $tumor_ref = ($base2 eq $ref_base) ? $tumor2_count : $tumor1_count;
        my $normal_alt = ($base2 eq $ref_base) ? $normal1_count : $normal2_count;
        my $tumor_alt = ($base2 eq $ref_base) ? $tumor1_count : $tumor2_count;
        my $alt_base = ($base2 eq $ref_base) ? $base1 : $base2;
        my $first_base = ($base2 eq $ref_base) ? $base2 : $base1;

        my $geno = call_genotype($normal_alt, $total_norm, $ra_het_count_limits);

        # first check for potential somatic alterations:
        if (($total_norm >= $min_tumor_reads) && ($total_tumor >= $min_tumor_reads) && ($total_alt >= $min_tumor_reads)) {
            print SOM "$chr\t$pos\t$ref_base\t$first_base\t$normal_ref\t$tumor_ref\t$alt_base\t$normal_alt\t$tumor_alt\t$geno\n";
        }

        next if ($geno eq 'und');
        # then check if this position is a potential "best position" in its window:
        if (!$cur_win_start) { 
            $cur_chr = $chr;
            $cur_win_start = $pos;
            $cur_win_end = $pos + $insert;
        }

        if (($chr ne $cur_chr) || ($pos > $cur_win_end)) { # process old window, create new
            print CNV "$best_string";
            $best_string = '';
            $best_geno = '';
            $best_cov = '';

            $cur_chr = $chr;
            $cur_win_start = $pos;
            $cur_win_end = $pos + $insert;
        }

        # replace best string, if appropriate:

        if ((!$best_string) || ($geno eq $best_geno && $total_norm > $best_cov) || ($best_geno ne 'het' && $geno eq 'het')) {
            $best_string = "$chr\t$pos\t$ref_base\t$base1\t$normal1_count\t$tumor1_count\t$base2\t$normal2_count\t$tumor2_count\t$geno\n";
            $best_geno = $geno;
            $best_cov = $total_norm;
        }
    }

    # print final window:
    print CNV "$best_string";

    close COUNTS;
    close CNV;
    close SOM;

}

sub read_min_max_file {
    my $file = shift;
    my $ra_min_max = [];

    open MINMAX, $file
        or die "Couldn\'t open $file to read min and max counts for genotypes!\n";

    while (<MINMAX>) {
        next if ((/^#/) || (/^COV/));
        chomp;
        if (/^(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my ($cov, $max_hom, $min_het, $max_het, $min_hnr) = ($1, $2, $3, $4, $5);
            $ra_min_max->[$cov] = { 'max_hom' => $max_hom,
                                    'min_het' => $min_het,
                                    'max_het' => $max_het,
                                    'min_hnr' => $min_hnr }; 
        }
    } 
    close MINMAX;
    return $ra_min_max;
}

sub call_genotype {
    my $alt_count = shift;
    my $total_count = shift;
    my $ra_count_limits = shift;

    my $rh_limits;
    if ($total_count > $#{$ra_count_limits}) {
        $rh_limits = {'max_hom' => 17,
                      'min_het' => 0.20*$total_count,
                      'max_het' => 0.80*$total_count };
    }
    else {
        $rh_limits = $ra_count_limits->[$total_count];
    }
    if ($rh_limits) {
        if ($alt_count <= $rh_limits->{'max_hom'}) {
            return 'hom';
        }
        elsif ($alt_count >= $rh_limits->{'min_het'} && 
                $alt_count <= $rh_limits->{'max_het'}) {
            return 'het';
        }
    }
    return 'und';
}

sub test_counts {

    # test counts for significance with Fisher's exact test (without mult testing corr)
    my $som_file = shift;
    my $het_file = shift;

    foreach my $file ($som_file, $het_file) {
        my $type = ($file eq $som_file) ? 'som' : 'het';

        my $r_command_file = "$file.r";
    
        open COM, ">$r_command_file"
            or die "Couldn\'t open $r_command_file for writing: $!\n";
    
        print COM <<"DOC";

library("statmod");
con <- file("$file", "r");
while (length(input <- readLines(con, n=1000)) > 0) {
    for (i in 1:length(input)) {
        line <- input[i];
        linevec <- strsplit(line, split="\\t");
        allele_counts <- c(linevec[[1]][[5]], linevec[[1]][[6]], linevec[[1]][[8]], linevec[[1]][[9]]);
        allele_counts <- as.numeric(allele_counts);
        geno <- linevec[[1]][[10]];
        dim(allele_counts) <- c(2,2);
        if ((("$type" == "som") && (geno == "hom")) || (("$type" == "het") && (geno == "het"))) {
            exact_result <- fisher.test(allele_counts);
            pvalue <- exact_result["p.value"];
            output <- paste(linevec[[1]][[1]], linevec[[1]][[2]], linevec[[1]][[3]], linevec[[1]][[4]], linevec[[1]][[5]], linevec[[1]][[6]], linevec[[1]][[7]], linevec[[1]][[8]], linevec[[1]][[9]], linevec[[1]][[10]], pvalue, sep=":");
        }
        else {
            output <- paste(linevec[[1]][[1]], linevec[[1]][[2]], linevec[[1]][[3]], linevec[[1]][[4]], linevec[[1]][[5]], linevec[[1]][[6]], linevec[[1]][[7]], linevec[[1]][[8]], linevec[[1]][[9]], linevec[[1]][[10]], "NA", sep=":");
        }
        print(output, quote=FALSE, max.levels=0);
    }
}

DOC
 
        close COM
            or die "Couldn\'t close file $r_command_file: $!\n";
        
        my $r_pipe = "$R_EXE --file=$r_command_file | ";
        open ROUTPUT, "$r_pipe"
            or die "Couldn\'t open pipe to $r_pipe!\n";
       
        my $test_output = $file;
        $test_output =~ s/\.txt$//;
        $test_output .= '.tests.txt'; 

        open TEST, ">$test_output"
            or die "Couldn\'t open $test_output for writing: $!\n";

        my $last = 0;
        my @vals = ();
        while (<ROUTPUT>) {
            if (/^\[1\]\s(\S+)$/) {
                 my $output = $1;
                 $output =~ s/:/\t/g;
                 print TEST "$output\n";
            }
        }
        while (<ROUTPUT>) {
            next;
        }

        close TEST;
        close ROUTPUT;
    }
}

# subroutine to apply Benjamini-Hochberg procedure to p-values to obtain q values

sub bh_correct_tests {
    my $test_fof = shift;
    my $max_q = shift;
    my $outfile = shift;

    my @files = ();
    open FILES, "$test_fof"
        or die "Couldn\'t open $test_fof for reading: $!\n";
    while (<FILES>) {
        chomp;
        push @files, $_;
    }
    close FILES;

    my @lines = ();
    my @pvalues = ();
    my $no_tests = 0;
    my $line_index = 0;
    foreach my $test_file (@files) {
        open TESTS, $test_file
            or die "Couldn\'t open $test_file for reading: $!\n";
    
        while (<TESTS>) {
            my $line = $_;
            push @lines, $line;
            chomp $line;
            my @fields = split /\s/, $line;
            my $geno = $fields[$#fields - 1];
            my $pvalue = $fields[$#fields];
            $no_tests++ if (($pvalue ne 'NA') || ($geno eq 'und'));;
            if (($pvalue ne 'NA') && ((!$max_q) || ($pvalue <= $max_q))) {
                push @pvalues, {'pvalue' => $pvalue, 'line_index' => $line_index};
            }
            $line_index++;
        }
        close TESTS;
    }
   
    print "Applying Benjamini-Hochberg correction to somatic change predictions with $no_tests tests.  Results in $outfile.\n"; 

    # now order and assign q-values:
    
    my $index = 1;
   
    my @lines_to_print = ();
    my @sorted_pvalues = sort bypthenna @pvalues;
    foreach my $rh_line (@sorted_pvalues) {
        my $pvalue = $rh_line->{'pvalue'};
        my $qvalue = ($pvalue eq "NA") ? "NA" : $no_tests*$rh_line->{'pvalue'}/$index;
        last if (($max_q) && ($qvalue ne "NA") && ($qvalue > $max_q));
    
        $index++;
        $lines[$rh_line->{'line_index'}] =~ s/\n/\t$qvalue\n/;
        my $this_line = $lines[$rh_line->{'line_index'}];
        push @lines_to_print, {'line_index' => $rh_line->{'line_index'}, 'line' => $lines[$rh_line->{'line_index'}]};
    }

    sub bypthenna {
        my $avalue = $a->{'pvalue'};
        my $bvalue = $b->{'pvalue'};
        if ($avalue eq "NA") {
            return 1;
        }
        elsif ($bvalue eq "NA") {
            return -1;
        }
        else {
            return $avalue <=> $bvalue;
        }
    }
    
    # and write out the lines to stdout:
   
    open OUT, ">$outfile"
        or die "Couldn\'t open $outfile for writing: $!\n";
 
    foreach my $line (sort {$a->{'line_index'} <=> $b->{'line_index'}} @lines_to_print) {
        print OUT $line->{'line'};
    }
    close OUT;

} # end bh_correct_tests

sub write_vs_file {
    my $outfile = shift;
    my $vsfile = shift;

    open SOM, "$outfile"
        or die "Couldn\'t open $outfile for reading: $!\n";

    open VS, ">$vsfile"
        or die "Couldn\'t open $vsfile for writing: $!\n";

    print VS "Index\tChr\tLeftFlank\tRightFlank\tref_allele\tvar_allele\tmuttype\tnormal_covg\ttumor_covg\tnormal_ratio\ttumor_ratio\tq_value\n";
    my $index = 1;
    while (<SOM>) {
        chomp;
        my ($chr, $pos, $ref, $allele1, $norm1_count, $tumor1_count, $allele2, $norm2_count, $tumor2_count, $gen, $pvalue, $qvalue) = split /\t/, $_;
        $ref = uc $ref;
        my $normal_covg = $norm1_count + $norm2_count;
        my $tumor_covg = $tumor1_count + $tumor2_count;
        my $normal_ratio = $norm2_count/$normal_covg;
        my $tumor_ratio = $tumor2_count/$tumor_covg;
        my $lfe = $pos - 1;
        my $rfs = $pos + 1;
        print VS "$index\t$chr\t$lfe\t$rfs\t$ref\t$allele2\tSNP\t$normal_covg\t$tumor_covg\t$normal_ratio\t$tumor_ratio\t$qvalue\n";
        $index++;
    }

    close SOM;
    close VS;

} ## end write_vs_file
  
sub write_vcf_file {
    my $outfile = shift;
    my $vcffile = shift;

    open SOM, "$outfile"
        or die "Couldn\'t open $outfile for reading: $!\n";

    open VCF, ">$vcffile"
        or die "Couldn\'t open $vcffile for writing: $!\n";

    print VCF "##fileformat=VCFv4.0\n";
    my ($sec, $min, $hour, $mday, $mon, $year ) = localtime();
    $mon++;
    $year += 1900;
    printf VCF "##fileDate=%d%02d%02d\n", $year, $mon, $mday;
    print VCF "##source=SHiMMer\n";

    # include info for RM flag for repeat-masked sequence:
    print VCF "##INFO=<ID=RM,Number=0,Type=Flag,Description=\"Lower-case reference\">\n";

    # included genotype id's:
    print VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    #print VCF "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    print VCF "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";

    # print required eight fields:
    print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n";

    my $index = 1;
    while (<SOM>) {
        chomp;
        my ($chr, $pos, $ref, $allele1, $norm1_count, $tumor1_count, $allele2, $norm2_count, $tumor2_count, $gen, $pvalue, $qvalue) = split /\t/, $_;
        my $normal_covg = $norm1_count + $norm2_count;
        my $tumor_covg = $tumor1_count + $tumor2_count;
        my $normal_ratio = $norm2_count/$normal_covg;
        my $tumor_ratio = $tumor2_count/$normal_covg;
        my $info_flag = ($ref eq uc $ref) ? '.' : 'RM';
        $ref = uc $ref;
        print VCF "$chr\t$pos\t$ref\t$allele2\t$pvalue\t.\t$info_flag\tGT:DP\t0/0:$normal_covg\t0/1:$tumor_covg\n";
        $index++;
    }

    close SOM;
    close VCF;

} ## end write_vcf_file

sub annotate_variants {
    my $annovar_db = shift;
    my $buildver = shift;
    my $vs_file = shift;
    my $ann_vs_file = shift;

    # write Annovar input file:
    open VS, "$vs_file"
        or die "Couldn\'t open $vs_file: $!\n";
  
    my $anv = "$vs_file.annovar.txt";
    open ANV, ">$anv"
        or die "Couldn\'t open $anv for writing: $!\n";

    while (<VS>) {
        chomp;
        my ($index, $chr, $lfe, $rfs, $ref, $var, $rest) = split /\t/, $_;
        next if ($index eq 'Index');
        my $pos = $lfe + 1;
        print ANV "$chr\t$pos\t$pos\t$ref\t$var\t$index\n";
    }

    close ANV; 
    close VS;

    my $annovar_cmd = "$ANNOVAR_EXE --buildver $buildver --geneanno --dbtype knowngene --hgvs $anv $annovar_db";
    (system("$annovar_cmd") == 0)
        or die "Some problem running $ANNOVAR_EXE!\n";

    # parse output, add columns to VS file:
    my $rh_function = {};
    open VAR, "$anv.variant_function"
        or die "Couldn\'t open $anv.variant_function: $!\n";
    while (<VAR>) {
        chomp;
        my @fields = split /\t/;
        $rh_function->{$fields[7]} = {'loc_type' => $fields[0], 'gene' => $fields[1]};
    }
    close VAR;

    my $rh_cons = {};
    open EXON, "$anv.exonic_variant_function"
        or die "Couldn\'t open $anv.exonic_function: $!\n";
    while (<EXON>) {
        chomp;
        my @fields = split /\t/;
        my $type = $fields[1];
        $type =~ s/\s/_/g;
        $rh_cons->{$fields[8]} = {'type' => $type, 'cons' => $fields[2]};
    }
    close EXON;

    print "About to write annotated file $vs_file!\n";
    # now write annotated VS file with extra columns:
    open VS, "$vs_file"
        or die "Couldn\'t open $vs_file: $!\n";
  
    open ANN, ">$ann_vs_file"
        or die "Couldn\'t open $ann_vs_file for writing: $!\n";

    while (<VS>) {
        chomp;
        my @fields = split /\t/, $_;
        my @new_fields = @fields[0..6];
        if ($new_fields[0] eq 'Index') {
            push @new_fields, qw( type Gene_name loc_type consequence );
        }
        else {
            my $index = $new_fields[0];
            my $type = $rh_cons->{$index}->{'type'} || $rh_function->{$index}->{'loc_type'} || 'NA';
            my $gene = $rh_function->{$index}->{'gene'} || 'NA';
            my $loc_type = $rh_function->{$index}->{'loc_type'} || 'NA';
            my $consequence = $rh_cons->{$index}->{'cons'} || 'NA';
            push @new_fields, ($type, $gene, $loc_type, $consequence);
        }
        push @new_fields, @fields[7..$#fields];
        my $annotated_string = join "\t", @new_fields;
        print ANN "$annotated_string\n";
    }

    close ANN;
    close VS;

} ## end annotate_variants 

sub write_symbols {
    my $input_file = shift; 
    my $output_file = shift; 
    my $outdir = shift;
    my $max_q = shift;
    my $bw_exe = shift;
    my $viterbi_exe = shift;
    my $hmm_model = shift;

    my $depth_normalization = calculate_depth_norm($input_file);
    open DIFFS, $input_file
        or die "Couldn\'t open $input_file: $!\n";

    my @symbols = ();
    my $symbol_string = '';
    my $gauss_string = '';
    my $no_symbols = 0;
    my $no_homs = 0;
    while (<DIFFS>) {
        chomp;
        my @fields = split /\t/, $_;
        my $no_fields = @fields;
        if ($no_fields == 12) {
            my $chr = $fields[0];
            my $pos = $fields[1];
            my $normal_a = $fields[4];
            my $tumor_a = $fields[5];
            my $normal_b = $fields[7];
            my $tumor_b = $fields[8];
            my $geno = $fields[9];
            my $pvalue = $fields[10];
            my $qvalue = $fields[11];
            next if (!($normal_a + $normal_b)); # need some normal reads!
            my $trr = ($tumor_a + $tumor_b)/($normal_a + $normal_b) * $depth_normalization;
            next if (($Opt{'gauss'}) && ($trr == 0));
            my $logtrr = log($trr) if ($Opt{'gauss'});
            my $sig = (($geno eq 'hom') || ($qvalue > $max_q)) ? 'I' : 'S';

            my $trr_index = ($trr < 0.7) ? 0 : ($trr < 1.3) ? 1 : ($trr < 1.7) ? 2 : 3;
            my $geno_index = ($geno eq 'hom') ? 1 : ($sig eq 'I') ? 2 : 3;
            my $symbol_number = ($Opt{'gauss'}) ? $geno_index : 3*$trr_index + $geno_index;
            $symbol_string .= "$symbol_number ";
            if ($Opt{'gauss'}) { # generate Gaussian emission string, track hom ratio
                $gauss_string .= "$trr ";
                $no_homs++ if ($geno_index == 1);
            }
            push @symbols, {'chr' => $chr, 'pos' => $pos, 'geno' => $geno, 'sig' => $sig, 'trr' => $trr };
            $no_symbols++;
        }
    }
    open SYMB, ">$output_file"
        or die "Couldn\'t open $output_file for writing: $!\n";
    print SYMB "T=$no_symbols\n$symbol_string\n";
    print SYMB "$gauss_string\n" if ($Opt{'gauss'});
    close SYMB;
    close DIFFS;

    # run baum-welch
    my $inithmm = "$outdir/initialhmm.hmm";
    my $outhmm = "$outdir/trainedhmm.hmm";
    if ($Opt{'gauss'}) {
        my $perc_hom = $no_homs/$no_symbols;
        my $neut_in = (1 - $perc_hom)*0.99;
        my $neut_sig = (1 - $perc_hom)*0.01;
        my $alt_in = (1 - $perc_hom)*0.8;
        my $alt_sig = (1 - $perc_hom)*0.2;
        open INIT, ">$inithmm"
            or die "Couldn\'t open $inithmm for writing: $!\n";
        print INIT "M= 3\nN= 4\nA:\n0.9999997 0.0000001 0.0000001 0.0000001\n0.0000001 0.9999997 0.0000001 0.0000001\n0.0000001 0.0000001 0.9999997 0.0000001\n0.0000001 0.0000001 0.0000001 0.9999997\nB:\n$perc_hom $alt_in $alt_sig\n$perc_hom $neut_in $neut_sig\n$perc_hom $alt_in $alt_sig\n$perc_hom $alt_in $alt_sig\nGB:\n0.6 0.2\n0.9 0.2\n0.9 0.2\n1.3 0.2\npi:\n0.001 0.997 0.001 0.001\n";
        close INIT;
        system("$bw_exe $inithmm $output_file > $outhmm");
    }

    # run viterbi here
    if ($Opt{'gauss'}) {
        $hmm_model = $outhmm;
    }

    my $outstatefile = ($Opt{'gauss'}) ? "$outdir/cnv.gaussstates.txt" : "$outdir/cnv.states.txt";
    my $gaussstring = ($Opt{'gauss'}) ? '-g' : '';
    system("$viterbi_exe $hmm_model $output_file $gaussstring > $outstatefile")
        or die "Couldn\'t run $viterbi_exe on $output_file with model $hmm_model.\n";

    # store predicted states in @symbols:
    open STATES, "$outstatefile"
        or die "Couldn\'t read $outstatefile: $!\n";
    while (<STATES>) {
        next if (!/log probabilities/);
        for (my $i=1; $i<=4; $i++) {
            $_ = <STATES>;
        }
        chomp;
        my @states = split /\s/, $_;
        my $no_states = @states;
        my $no_symbols = @symbols;
        if ($no_states != $no_symbols) {
            die "Failure in Viterbi algorithm--wrong number of states.\n";
        }
        for (my $i=0; $i<=$#symbols; $i++) {
            $symbols[$i]->{'state'} = $states[$i];
        }
    }
    close STATES;
    
    return [@symbols];

} # end write_symbols

sub calc_power_file {
    my $input_file = shift;
    my $outdir = shift;
    my $bam1 = shift;
    my $bam2 = shift;
    my $ref_fasta = shift;
    my $printcounts_exe = shift;

    # call mpileup (via the c-script "printCompCounts"), and record bed regions for somatic and CNA power values.
    my $printcounts_call = "$printcounts_exe -bam1 $bam1 -bam2 $bam2 -fasta $ref_fasta";

    my $perc_sum = 0;
    my $perc_no = 0;
    open SOM, "$input_file"
        or die "Couldn\'t open $input_file: $!\n";
    while (<SOM>) {
        chomp;
        my ($chr, $pos, $ref, $base1, $norm1, $tumor1, $base2, $norm2, $tumor2, $geno) = split /\t/, $_;
        $perc_sum += ($norm1 + $norm2)/($norm1 + $tumor1 + $norm2 + $tumor2);
        $perc_no++;
    }
    close SOM;
    my $tumor_perc_estimate = ($perc_no) ? int(100*$perc_sum/$perc_no*2) : 100;
    $tumor_perc_estimate = 100 if ($tumor_perc_estimate > 100);

    my $purity = int($perc_sum/$perc_no*100);

    my $rh_somatic_power = read_somatic_power($purity);
    open COUNTS, "$printcounts_call | "
        or die "Couldn\'t execute $printcounts_call!\n";

    my ($current_chr, $current_start, $current_pos, $current_color);
    open SOMBED, ">$outdir/somatic_power.bed"
        or die "Couldn\'t open $outdir/somatic_power.bed for writing: $!\n";

    while (<COUNTS>) {
        my ($chr, $pos, $ref_base, $base1, $normal1_count, $tumor1_count, $base2, $normal2_count, $tumor2_count) = split /\t/, $_;
        my $total_norm = $normal1_count + $normal2_count;
        my $total_tumor = $tumor1_count + $tumor2_count;

        my $power_color = power_color($rh_somatic_power, $total_norm, $total_tumor);
        if (!$current_color || !$current_chr || $chr ne $current_chr || $pos != $current_pos + 1 || $power_color ne $current_color) {
            # write out old entry:
            if ($current_color) {
                print SOMBED "$current_chr\t$current_start\t$current_pos\t-\t-\t-\t-\t-\t$current_color\n";
            }
            $current_start = $pos - 1;
            $current_color = $power_color;
        }
        $current_chr = $chr;
        $current_pos = $pos;
    }
    close COUNTS;
    if ($current_color) {
        print SOMBED "$current_chr\t$current_start\t$current_pos\t-\t-\t-\t-\t-\t$current_color\n";
    }

    close SOMBED;
    close COUNTS;

} # end calc_power_file

sub read_somatic_power {
    my $purity = shift;
    my $som_power_file = "/home/nhansen/projects/shimmer/power_graphs/somatic_power_table.txt";
    open SOM, $som_power_file
        or die "Couldn\'t open $som_power_file: $!\n";
    my $rh_som_power = {};
    my $best_purity;
    while (<SOM>) {
        chomp;
        my ($perc_tumor, $normal_reads, $tumor_reads, $power) = split;
        next if ($perc_tumor > $purity);
        if (!$best_purity) {
            $best_purity = $perc_tumor;
        }
        elsif ($perc_tumor != $best_purity) {
            last;
        }
        $rh_som_power->{$normal_reads}->{$tumor_reads} = $power;
        $rh_som_power->{$normal_reads}->{'max'} = $tumor_reads;
        $rh_som_power->{'max'} = $normal_reads;
    }
    close SOM;

    return $rh_som_power;
}

sub power_color {
    my $rh_power = shift;
    my $normal_total = shift;
    my $tumor_total = shift;

    my $normal_tens = 10*int($normal_total/10); 
    my $tumor_tens = 10*int($tumor_total/10);

    $normal_tens = $rh_power->{'max'} if ($normal_tens > $rh_power->{'max'});

    my $power = 0;
    if ($rh_power->{$normal_tens} && $rh_power->{$normal_tens}->{$tumor_tens}) {
        $power = $rh_power->{$normal_tens}->{$tumor_tens};
    }
    elsif ($rh_power->{$normal_tens}) {
        if ($tumor_tens > $rh_power->{$normal_tens}->{'max'}) {
            $power = $rh_power->{$normal_tens}->{'max'};
        }
    }
    my $power_color;
    for (my $i=0; $i<= $#SOM_COLORS; $i++) {
        if ($power >= $SOM_COLORS[$i]->{'power'}) {
            $power_color = $SOM_COLORS[$i]->{'color'};
        }
    }
    return $power_color;
}

sub calculate_depth_norm {
    my $file = shift;

    open DIFFS, $file
        or die "Couldn\'t open $file: $!\n";

    my $sum_ratio = 0;
    my $sumsq_ratio = 0;
    my $total_points = 0;
    my @norm_ratios = ();
    while (<DIFFS>) {
        chomp;
        my @fields = split /\t/, $_;
        my $line = $_;
        my $no_fields = @fields;
        if ($no_fields == 12) { # data line 
            my $pvalue = $fields[10];
    
            if (($pvalue eq 'NA') || ($pvalue < 0.5)) {
                next;
            }

            my $normal_a = $fields[4];
            my $tumor_a = $fields[5];
            my $normal_b = $fields[7];
            my $tumor_b = $fields[8];
            #next if (!($normal_a + $normal_b)); # need some tumor reads to assess anything
            next if (($normal_a + $normal_b) < 100 || ($tumor_a + $tumor_b) < 100); # need some reads to assess anything
            my $this_ratio = ($normal_a + $normal_b) / ($tumor_a + $tumor_b);
            $sum_ratio += $this_ratio;
            $sumsq_ratio += $this_ratio**2;
            push @norm_ratios, $this_ratio;
            $total_points++;
        }
    }
    close DIFFS;
    if (!$total_points) {
        die "Unable to calculate a normalization ratio for depths--too few \"normal\" points!\n";
    }
    my @sorted_norm_ratios = sort {$a <=> $b} @norm_ratios;
    my $median_ratio = $sorted_norm_ratios[int($#sorted_norm_ratios/2)];
    
    my $avg_ratio = $sum_ratio/$total_points;
    print "Calculated average read depth ratio of $avg_ratio, median $median_ratio (normal divided by tumor).\n";

    return $avg_ratio;

} # end calculate_depth_norm

sub plot_symbols {
    my $ra_symbols = shift;
    my $tmpdir = shift;

    my $ylow = 0.0;
    my $yhigh = 3.2;
    my $stateline = 3.0;

    # write files of points to plot

    my @all_chrs = ();
    my $last_chr;
    for (my $i=0; $i<=$#{$ra_symbols}; $i++) {
        my $rh_symbol = $ra_symbols->[$i];
        my $chr = $rh_symbol->{'chr'};
        my $symbol_file = ($Opt{'gauss'}) ? "$tmpdir/$chr.gausssymbols.txt" : "$tmpdir/$chr.symbols.txt";
        if (!$last_chr) {
            open POINTS, ">$symbol_file"
                or die "Unable to open $symbol_file for writing: $!\n";
            push @all_chrs, $chr;
        }
        elsif (($last_chr) && ($chr ne $last_chr)) {
            close POINTS;
            open POINTS, ">$symbol_file"
                or die "Unable to open $symbol_file for writing: $!\n";
            push @all_chrs, $chr;
        }
        my $pos = $rh_symbol->{'pos'};
        my $trr = $rh_symbol->{'trr'};
        my $geno = $rh_symbol->{'geno'};
        my $sig = $rh_symbol->{'sig'};
        my $state = $rh_symbol->{'state'};
        print POINTS "$pos\t$trr\t$geno\t$sig\t$state\n";
        $last_chr = $chr;
    }
    close POINTS;

    # now plot:
    foreach my $chr (@all_chrs) {
        # first colored plots of coverage, genotype, significance:
        my $r_cmd_file = ($Opt{'gauss'}) ? "$tmpdir/$chr.gausssymbols.r" : "$tmpdir/$chr.plot.r";
        my $symbol_file = ($Opt{'gauss'}) ? "$tmpdir/$chr.gausssymbols.txt" : "$tmpdir/$chr.symbols.txt";
        my $png_file = ($Opt{'gauss'}) ? "$tmpdir/$chr.gausssymbols.png" : "$tmpdir/$chr.symbols.png";
        open RCMD, ">$r_cmd_file"
            or die "Couldn\'t open $r_cmd_file for writing: $!\n";
        print RCMD <<"DOC";

datapoints <- read.table("$symbol_file", sep="\\t");
data <- as.matrix(datapoints);
bitmap(file="$png_file", type="png16m");
pos <- data[,1];
ratio <- as.numeric(data[,2]);
geno <- data[,3];
sig <- data[,4];
state <-data[,5];
threes <- rep($stateline, length(pos));
plot(pos, ratio, col=ifelse(sig=="S", "white", ifelse(geno=="het", "green", "lightblue")), ylim=c($ylow, $yhigh), pch=ifelse(sig=="S", 16, 1), cex=ifelse(sig=="S", 0.0, ifelse(geno=="het", 0.0, 0.1)));
points(pos, ratio, col=ifelse(sig=="S", "white", ifelse(geno=="het", "green", "white")), pch=ifelse(sig=="S", 16, 1), cex=ifelse(sig=="S", 0.0, ifelse(geno=="het", 0.1, 0.0)));
points(pos, ratio, col=ifelse(sig=="S", "red", ifelse(geno=="het", "white", "white")), pch=ifelse(sig=="S", 16, 1), cex=ifelse(sig=="S", 0.2, ifelse(geno=="het", 0.0, 0.0)));
points(pos, threes, col=ifelse(state==1, "red", ifelse(state==2, "green", ifelse(state==3, "yellow", "blue"))));
dev.off();

DOC
 
        close RCMD
            or die "Couldn\'t close file $r_cmd_file!\n";
        my $r_cmd = "$R_EXE --file=$r_cmd_file";
        system("$r_cmd")==0
            or warn "Failed to plot symbols for chromosome $chr!\n";

        # now histograms of normalized read depths:
        $r_cmd_file = "$tmpdir/$chr.hist.r";
        open RCMD, ">$r_cmd_file"
            or die "Couldn\'t open $r_cmd_file for writing: $!\n";
        print RCMD <<"DOC";

datapoints <- read.table("$tmpdir/$chr.symbols.txt", sep="\\t");
data <- as.matrix(datapoints);
bitmap(file="$tmpdir/$chr.covghist.png", type="png16m");
ratio <- as.numeric(data[,2]);
hist(ratio, breaks=50);
dev.off();

DOC

        close RCMD
            or die "Couldn\'t close file $r_cmd_file!\n";
        $r_cmd = "$R_EXE --file=$r_cmd_file";
        system("$r_cmd")==0
            or warn "Failed to plot coverage histogram for chromosome $chr!\n";

    } 

} # end plot_symbols

sub process_commandline {
    
    # Set defaults here
    %Opt = ( 
             max_q => 0.05, insert => 300, min_max => '/home/nhansen/projects/shimmer/min_max_counts.txt', min_som_reads => 10, skip_cna => 1
           );
    GetOptions(\%Opt, qw(
                region=s ref=s counts som_file=s
                het_file=s bh vs_file=s vcf_file=s max_q=f test_fof=s
                outfile=s outdir=s symbols input=s plots power covg minqual=i min_som_reads=i
                viterbi insert=i min_max=s annovar=s annovardb=s buildver=s skip_tests skip_cna
                annotate gauss help+ version verbose 
               )) || pod2usage(0);
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}); }
    if ($Opt{version}) { die "$0, ", q$Revision: $, "\n"; }

    # argument checking:
    if ($Opt{'annovardb'}) {
        $ANNOVAR_EXE = $Opt{'annovar'} || `which $ANNOVAR_EXE`;
        chomp $ANNOVAR_EXE;
        if (!$ANNOVAR_EXE || !(-e $ANNOVAR_EXE)) {
            die "Cannot find $ANNOVAR_EXE--ignoring opt --annovardb!\n";
        }
        if (!$Opt{'buildver'}) {
            die "Must specify a build version (e.g., hg18) for annovar with --buildver\n";
        }
    }
}

=pod

=head1 NAME

B<shimmer> - Program to call somatic single base changes and copy number alterations from matched tumor and normal next generation sequences.

=head1 SYNOPSIS

B<shimmer> I<normal sample bam file> I<tumor sample bam file>

=head1 DESCRIPTION

This script uses samtools to process two BAM formatted files (http://samtools.sourceforge.net) and call differences between the two files across a specified region.

=head1 INPUT

The first argument to shimmer is the path of a fasta-formatted file for the reference sequence.  This fasta file must have a corresponding samtools index file with the same name except for an appended ".fai".

The second argument to bam2mpg is the path of a BAM-formatted file of aligned sequencing reads.  This file must be sorted and indexed using samtools prior to running bam2mpg.

=head1 OUTPUT

=head1 OPTIONS

=over 5

=item B<--region> I<chr>

=item B<--region> I<chr:start-end>

This option specifies a region as a reference entry optionally followed by a position range, and causes variants to be called only in that region.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut

