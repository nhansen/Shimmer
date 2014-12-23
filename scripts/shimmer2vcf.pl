#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(:ALL);
use GTB::FASTA;
use GTB::File qw(Open :use_bgzip);
use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr( q$Rev: 0$, 4 ) / 10000;

=head1 NAME

shimmer2vcf.pl - script to convert shimmer calls to VCFv4.0/TCGA1.2 for the TCGA Benchmark4 somatic caller evaluation.

=head1 SYNOPSIS

  shimmer2vcf.pl --mpg <normal mpg.out> --snvs <shimmer snv varsifter> --indels <shimmer indel varsifter> --ref <fasta> --sample1 <sample_name> --sample2 <sample_name> <--bam1 bamfilename> <--bam2 bamfilename> [-notabix]

=head1 DESCRIPTION

This script reads in MPG output from bam2mpg and/or somatic calls from a Shimmer Varsifter file and writes a VCFv4.0 file containing the same calls.  If output files end in '.gz' it will compress with bgzip and create tabix indexes. To prevent index creation pass -notabix option

=cut

#------------
# Begin MAIN
#------------

process_commandline();
my $mpg_file = $Opt{'mpg'};
my $snv_file = $Opt{'snv'};
my $indel_file = $Opt{'indel'};
my $bamfile1 = $Opt{'bam1'};
my $bamfile2 = $Opt{'bam2'};
my $ref_fasta   = $Opt{'ref'};
my $sample1_name = $Opt{'sample1'};
my $sample2_name = $Opt{'sample2'};
my $outfile = $Opt{'outfile'};

my $fasta_db = GTB::FASTA->new($ref_fasta);

my $mode = $Opt{'append'} ? 'a' : 'w';

print_header( $sample1_name, $sample2_name ) unless ( $Opt{'noheader'} );

# first process germline variants:
my $rh_germline_variants = ($mpg_file) ? read_germline_variants($mpg_file) : {};
my $rh_somatic_variants = read_somatic_variants($snv_file, $indel_file);
write_vcf_lines($rh_germline_variants, $rh_somatic_variants);

#if ( $outfile =~ /gz$/ && $Opt{tabix} ) {
#    unlink "$snv_outfile.tbi";  # tabix won't overwrite, so delete first
#    my $cmd = " $Opt{tabix} -p vcf $snv_outfile ";
#    my $w = `$cmd`;
#}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ref => '/scratch/fasta/hg18/hg18.mfa', shorten => 1, sample1 => 'Unknown', sample2 => 'Unknown', min_score => 10, min_ratio => 0.5, hicov => 500  );
    GetOptions(
        \%Opt, qw(
            manual help+ version
            verbose ref=s mpg=s snv=s indel=s bam=s
            sample1=s sample2=s outfile=s chr=s noq20 passonly
            append gzip noheader tabix! min_score=i min_ratio=f min_somatic_score=i
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "shimmer2vcf.pl, ", q$Revision: 4129 $, "\n";
	}

	#if ( !exists( $Opt{mpg} ) ) {
	#	die "Specify input normal mpg filename with option --mpg.\n";
	#}

	if ( !exists( $Opt{snv} ) ) {
		die "Specify input tumor shimmer filename with option --snv.\n";
	}

	if ( !$Opt{outfile} ) {
		$Opt{outfile} = "$Opt{sample2}.mpv.vcf";
		if ($Opt{gzip}) {
			$Opt{outfile} .= ".gz";
		}
	}

	if ( !defined $Opt{tabix} ) {
		$Opt{tabix} = `which tabix`;
		chomp $Opt{tabix};
	}
}

sub print_header {
	my $sample1_name = shift;
	my $sample2_name = shift;

	print "##fileformat=VCFv4.1\n";
	print "##tcgaversion=1.2\n";
	my ( $sec, $min, $hour, $mday, $mon, $year ) = localtime();
	$mon++;
	$year += 1900;
	printf "##fileDate=%d%02d%02d\n", $year, $mon, $mday;
	print "##center=\"NHGRI\"\n";
	print "##reference=GRCh37\n";
        print "##vcfProcessLog=<InputVCFSource=<shimmer.pl,shimmer2vcf.pl>,InputVCFVer=<0.1,0.1>>\n";

        print "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">\n";

        # filtering info:
        print "##FILTER=<ID=n2p,Description=\"Allele frequency greater than 2\% in normal\">\n";
        print "##FILTER=<ID=n1p,Description=\"Allele frequency greater than 1\% in normal\">\n";
        if (!$Opt{'noq20'}) {
            print "##FILTER=<ID=q20,Description=\"False discovery rate greater than 20\%\">\n";
        }
        if ($Opt{'hicov'}) {
            print "##FILTER=<ID=hicov,Description=\"Combined coverage greater than $Opt{'hicov'}\">\n";
        }

	# included genotype id's:
	print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	print "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
	print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position in the sample\">\n";
        print "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Depth of reads supporting alleles 0/1/2/3...\">\n"; 
        print "##FORMAT=<ID=BQ,Number=.,Type=Integer,Description=\"Average base quality for reads supporting alleles\">\n";
        print "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown\">\n";
        print "##FORMAT=<ID=SSC,Number=1,Type=Integer,Description=\"Somatic score between 0 and 255\">\n";

        print "##SAMPLE=<ID=$sample1_name>\n";
        print "##SAMPLE=<ID=$sample2_name>\n";

	# print required eight fields:
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample1_name\t$sample2_name\n";
}

sub read_germline_variants {
    my $mpg_file = shift;

    my $mpg_fh = Open($mpg_file);

    my $rh_germline = {}; 
    while (<$mpg_fh>) {
        if (/^MPG(_SNV){0,1}\s(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\d+)\s([01])(\s(\d+)){0,1}$/) # SNP!
        {
            my ( $chr, $pos, $ref_allele, $genotype, $score, $ref_nonref, $coverage ) = ( $2, $3, $4, $5, $6, $7, $9 );
            if (($Opt{'min_score'} && $score < $Opt{'min_score'}) || ($Opt{'min_ratio'} && $score < $Opt{'min_ratio'}*$coverage)) {
                next;
            }
            my $uc_ref_allele = uc $ref_allele;
            my $info_flag = '.';
            $ref_allele = $uc_ref_allele;
            my @alleles = split //, $genotype;
            my %allele_hash = map { $_ => 1 } @alleles;
            my @nonref_alleles = grep { $_ ne $ref_allele } keys %allele_hash;
            @nonref_alleles = sort @nonref_alleles;
            my $alt_alleles = (@nonref_alleles) ? join ',', @nonref_alleles : '.';
            my @all_alleles = @nonref_alleles;
            unshift @all_alleles, $ref_allele;    # include all alleles
            my %index_hash = map { $all_alleles[$_] => $_ } ( 0 .. $#all_alleles );
            my @allele_indices = map { $index_hash{$_} } @alleles;
            @allele_indices = sort {$a <=> $b} @allele_indices;
            my $index_genotype = join '/', @allele_indices;
            my $filter = 'PASS';
            my $vcf_line = "$chr\t$pos\t.\t$ref_allele\t$alt_alleles\t$score\t$filter\t$info_flag\tGT:GQ:SS:DP:AD:BQ\t$index_genotype:$score:1:$coverage:.:.\t.:.:.:.:.:.\n";
            $rh_germline->{$chr}->{$pos}->{'SNV'} = $vcf_line;
        }
        elsif ( /^MPG(_DIV){0,1}\s(\S+)\s(\d+):(\d+)\s(\S+)\s(\S+)\s(\d+)\s([01])(\s(\d+)){0,1}$/ )
        {
            my ( $chr, $lfe, $rfs, $ref_allele, $genotype, $score, $ref_nonref, $coverage ) = ( $2, $3, $4, $5, $6, $7, $8, $10 );
            if (($Opt{'min_score'} && $score < $Opt{'min_score'}) || ($Opt{'min_ratio'} && $score < $Opt{'min_ratio'}*$coverage)) {
                next;
            }
            if ( ( !$Opt{'keep_ref_div'} ) && ( $ref_nonref == 0 ) ) {
                next;
            }
            next if ( $lfe < 1 );    # something weird going on!
            $ref_allele = uc $ref_allele;
            $ref_allele = '' if ( $ref_allele eq '*' );
            my @alleles = split /:/, $genotype;
            foreach my $allele (@alleles) {    # get rid of asterisks
                $allele =~ s:\*::g;
            }
    
            my $pos = $lfe + 1;
    
            # shorten repetitive indels:
            if ( $Opt{'shorten'} ) {
                while ( ( $ref_allele ne '' ) && ( !grep { $_ eq '' } @alleles ) ) {
                    my $chop = 1;
                    my $last_char = ( $ref_allele =~ /(.)$/ ) ? $1 : '';
                    foreach my $allele (@alleles) {
                        my $allele_last_char = ( $allele =~ /(.)$/ ) ? $1 : '';
                        if ( $allele_last_char ne $last_char ) {
                            $chop = 0;
                            last;
                        }
                    }
                    if ($chop) {
                        chop $ref_allele;
                        foreach my $allele (@alleles) {
                            chop $allele;
                        }
                    }
                    else {
                        last;
                    }
                }
            }
    
            # now move coordinate one to the left as required by specs:
	    $pos--;
            my $ref_base = uc $fasta_db->seq( $chr, $pos, $pos );
            $ref_allele = $ref_base . $ref_allele;
            foreach my $allele (@alleles) {
                $allele = $ref_base . $allele;
            }
            my %allele_hash = map { $_ => 1 } @alleles;
            my @nonref_alleles = grep { $_ ne $ref_allele } keys %allele_hash;
            @nonref_alleles = sort @nonref_alleles;
            my $alt_alleles = (@nonref_alleles) ? join ',', @nonref_alleles : '.';
            my @all_alleles = @nonref_alleles;
            unshift @all_alleles, $ref_allele;    # include all alleles
            my %index_hash = map { $all_alleles[$_] => $_ } ( 0 .. $#all_alleles );
            my @allele_indices = map { $index_hash{$_} } @alleles;
            @allele_indices = sort {$a <=> $b} @allele_indices;
            my $index_genotype = join '/', @allele_indices;
            my $filter = 'PASS';
            my $vcf_line = "$chr\t$pos\t.\t$ref_allele\t$alt_alleles\t$score\t$filter\t.\tGT:GQ:DP:SS:AD:BQ\t$index_genotype:$score:$coverage:1:.:.\t.:.:.:.:.:.\n";
        
            push @{$rh_germline->{$chr}->{$pos}->{'DIV'}}, $vcf_line;
            #print "$vcf_line";
        }
    }
    return $rh_germline;
}

sub read_somatic_variants {
    my $shimmer_snv_file = shift;
    my $shimmer_indel_file = shift;

    my $somatic_fh = Open($shimmer_snv_file);
    my $rh_somatic = {};
    while (<$somatic_fh>) {
        next if /^Index/;
        chomp;
        my @fields = split /\t/, $_;
        my $chr = $fields[1];
        next if ($Opt{'chr'} && $chr ne $Opt{'chr'});
        my $lfe = $fields[2];
        my $rfs = $fields[3];
        my $ref_allele = $fields[4];
        my $pos = $lfe + 1;
        $ref_allele = uc $ref_allele;
        my $alt_allele = $fields[5];
        my $normal_depth = $fields[7];
        my $tumor_depth = $fields[8];
        my $normal_alt_freq = $fields[9];
        my $tumor_alt_freq = $fields[10];
        my @all_alleles = ($ref_allele, $alt_allele);
        my %index_hash = map { $all_alleles[$_] => $_ } ( 0 .. $#all_alleles );
        my $qvalue = $fields[11];
        my @filters = ();
        if (($qvalue > 0.2) && (!$Opt{'noq20'})) {
            push @filters, 'q20';
        }
        if ($normal_alt_freq > 0.02) {
            push @filters, 'n2p';
        }
        if ($normal_alt_freq > 0.01) {
            push @filters, 'n1p';
        }
        if (($Opt{'hicov'}) && ($normal_depth+$tumor_depth > $Opt{'hicov'})) {
            push @filters, 'hicov';
        }
        my $filter_field = @filters ? join ';', @filters : 'PASS';
        my $somatic_score = ($qvalue > 0) ? int(-10.0*log($qvalue)/log(10.0)) : '999';
        my @info_flags = ('SOMATIC');
        my $info_flag = join ';', @info_flags;
        my $tumor_gt = ($tumor_alt_freq > 0.75) ? '1/1' : '0/1';

        if (!$Opt{'passonly'} || $filter_field eq 'PASS') {
            my $vcf_line = "$chr\t$pos\t.\t$ref_allele\t$alt_allele\t$somatic_score\t$filter_field\t$info_flag\tGT:GQ:DP:SS:SSC:AD:BQ\t0/0:.:$normal_depth:.:.:.:.\t$tumor_gt:.:$tumor_depth:2:$somatic_score:.:.\n";
            $rh_somatic->{$chr}->{$pos}->{'SNV'} = $vcf_line;
        }
    }
    if ($shimmer_indel_file) {
        my $indel_fh = Open($shimmer_indel_file);
        while (<$indel_fh>) {
            next if /^Index/;
            chomp;
            my @fields = split /\t/, $_;
            my $chr = $fields[1];
            next if ($Opt{'chr'} && $chr ne $Opt{'chr'});
            my $lfe = $fields[2];
            my $rfs = $fields[3];
            my $ref_allele = $fields[4];
            $ref_allele = '' if ($ref_allele eq '*');
            my $alt_allele = $fields[5];
            $alt_allele = '' if ($alt_allele eq '*');
            # retrieve base one to left for specification:
            my $pos = $lfe;
            my $ref_base = uc $fasta_db->seq( $chr, $pos, $pos );
            $ref_allele = $ref_base . $ref_allele;
            $alt_allele = $ref_base . $alt_allele;
            my $normal_depth = $fields[7];
            my $tumor_depth = $fields[8];
            my $normal_alt_freq = $fields[9];
            my $tumor_alt_freq = $fields[10];
            next if ($normal_alt_freq > 0.05);
            my @all_alleles = ($ref_allele, $alt_allele);
            my %index_hash = map { $all_alleles[$_] => $_ } ( 0 .. $#all_alleles );
            my $qvalue = $fields[11];
            my @filters = ();
            if (($qvalue > 0.2) && (!$Opt{'noq20'})) {
                push @filters, 'q20';
            }
            if ($normal_alt_freq > 0.02) {
                push @filters, 'n2p';
            }
            if ($normal_alt_freq > 0.01) {
                push @filters, 'n1p';
            }
            if (($Opt{'hicov'}) && ($normal_depth+$tumor_depth > $Opt{'hicov'})) {
                push @filters, 'hicov';
            }
            my $filter_field = @filters ? join ';', @filters : 'PASS';
            my $somatic_score = ($qvalue > 0) ? int(-10.0*log($qvalue)/log(10.0)) : '999';
            my @info_flags = ('SOMATIC');
            my $info_flag = join ';', @info_flags;
            my $tumor_gt = ($tumor_alt_freq > 0.75) ? '1/1' : '0/1';
   
            if (!$Opt{'passonly'} || $filter_field eq 'PASS') { 
                my $vcf_line = "$chr\t$pos\t.\t$ref_allele\t$alt_allele\t$somatic_score\t$filter_field\t$info_flag\tGT:GQ:DP:SS:SSC:AD:BQ\t0/0:.:$normal_depth:.:.:.:.\t$tumor_gt:.:$tumor_depth:2:$somatic_score:.:.\n";
                push @{$rh_somatic->{$chr}->{$pos}->{'DIV'}}, $vcf_line;
            }
        }
    }
    return $rh_somatic;
}

sub write_vcf_lines {
    my $rh_germline = shift;
    my $rh_somatic = shift;

    my $rh_all_lines = {};
    my %all_chroms = ();
    foreach my $chr (keys %{$rh_germline}, keys %{$rh_somatic}) {
        $all_chroms{$chr} = 1;
    }

    foreach my $chr (sort keys %all_chroms) {
        my %all_pos = ();
        if ($rh_germline->{$chr}) {
            foreach my $pos (keys %{$rh_germline->{$chr}}) {
                $all_pos{$pos} = 1;
            }
        }
        if ($rh_somatic->{$chr}) {
            foreach my $pos (keys %{$rh_somatic->{$chr}}) {
                $all_pos{$pos} = 1;
            }
        }
        foreach my $pos (sort {$a <=> $b} keys %all_pos) {
            if ($rh_germline->{$chr} && $rh_germline->{$chr}->{$pos}) {
                print $rh_germline->{$chr}->{$pos}->{'SNV'} if ($rh_germline->{$chr}->{$pos}->{'SNV'});
            }
            if ($rh_somatic->{$chr} && $rh_somatic->{$chr}->{$pos}) {
                print $rh_somatic->{$chr}->{$pos}->{'SNV'} if ($rh_somatic->{$chr}->{$pos}->{'SNV'});
            }
            if ($rh_germline->{$chr} && $rh_germline->{$chr}->{$pos}) {
                if ($rh_germline->{$chr}->{$pos}->{'DIV'}) {
                    my $linestring = join '', @{$rh_germline->{$chr}->{$pos}->{'DIV'}};
                    print $linestring;
                }
            }
            if ($rh_somatic->{$chr} && $rh_somatic->{$chr}->{$pos}) {
                if ($rh_somatic->{$chr}->{$pos}->{'DIV'}) {
                    my $linestring = join '', @{$rh_somatic->{$chr}->{$pos}->{'DIV'}};
                    print $linestring;
                }
            }
        }
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--mpg> file.mpg.out

Specify the (optional) MPG output file to be converted to VCFv4.0.

=item B<--ref>

Specify the reference fasta file that was used in the bam2mpg run that resulted
in the MPG output file.

=item B<--sample>

Specify the sample name from which these genotypes come.  This will be used
inside the VCF file as the column header, and as the base name for the
default output file names.

=item B<--snv_outfile>

Specifies the name of the output VCF file for SNV variants.  If the filename ends
in "gz", will zip with bgzip. (default: sample.mpg.snv.vcf)

=item B<--div_outfile>

Specifies the name of the output VCF file for DIV variants.  If the filename ends
in "gz", will zip with bgzip.  If it is the same as the file specified in
--snv_outfile, will write all variants to a single file.
(default: sample.mpg.div.vcf)

=item B<--shorten>

When the --shorten option is specified, repetitive DIV variants will be shortened
to their shortest variant. (Note: currently, vcf2mpg.pl will not re-create the
original MPG file when this option is used).

=item B<--notabix>

By default when a output file ends in 'gz', tabix will create the corresponding index.
To override this behaviour specify --notabix.

=item B<--min_score> N

Minimum MPG score to create VCF record.  Variants with scores below this
threshold will not be included.  Default is no minimum score.

=item B<--gzip>

Generate compressed output files (suffixed with .gz extension).  Has no
effect when C<--snv_outfile> or C<--div_outfile> are explicitly specified.

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
