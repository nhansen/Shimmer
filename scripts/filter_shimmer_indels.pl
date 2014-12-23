#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::FASTA;
use GTB::File qw(Open :use_bgzip);
use GTB::Var::Polymorphism;
use GTB::BAM::Slice;
use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr( q$Rev: 0$, 4 ) / 10000;

=head1 NAME

filter_shimmer_indels.pl - script to re-assess frequency of alleles for called indel polymorphisms by requiring reads to span the entirety of the indel signature.

=head1 SYNOPSIS

  filter_shimmer_indels.pl --indels <shimmer indel varsifter> --bam1 <BAM file from normal sample> --bam2 <BAM file from tumor sample> --ref <fasta>

=head1 DESCRIPTION

This script reads in a Shimmer Varsifter file with indel calls, and examines the BAM files to re-calculate the number of reads for each indel, requiring reads to span the entire indel signature.  It then outputs a VarSifter file containing altered allele frequency fields.

=cut

#------------
# Begin MAIN
#------------

process_commandline();
my $indel_file = $Opt{'indels'};
my $bamfile1 = $Opt{'bam1'};
my $bamfile2 = $Opt{'bam2'};
my $ref_fasta   = $Opt{'ref'};

filter_somatic_indels($indel_file, $bamfile1, $bamfile2, $ref_fasta);

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( ref => '/scratch/fasta/hg18/hg18.mfa', maxnorm => 0.05, maxq => 0.20 );
    GetOptions(
        \%Opt, qw(
            manual help+ version maxnorm=f maxq=f min_som_reads=i
            verbose ref=s indels=s bam1=s bam2=s
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "filter_shimmer_indels.pl, ", q$Revision: 4129 $, "\n";
	}

}

sub filter_somatic_indels {
    my $shimmer_indel_file = shift;
    my $normal_bam = shift;
    my $tumor_bam = shift;
    my $ref_fasta = shift;

    my $fasta_db = GTB::FASTA->new($ref_fasta);
    my $normal_bam_obj = Bio::DB::Sam->new(
                            -bam => $normal_bam,
                            -fasta => $ref_fasta );

    my $normal_slice = GTB::BAM::Slice->new(-bam_obj => $normal_bam_obj);

    my $tumor_bam_obj = Bio::DB::Sam->new(
                            -bam => $tumor_bam,
                            -fasta => $ref_fasta);

    my $tumor_slice = GTB::BAM::Slice->new(-bam_obj => $tumor_bam_obj);

    my $indel_fh = Open($shimmer_indel_file);
    while (<$indel_fh>) {
        if (/^Index/) {
            print;
            next;
        }
        chomp;
        my @fields = split /\t/, $_;
        my $chr = $fields[1];
        my $lfe = $fields[2];
        my $rfs = $fields[3];

        my $ref_allele = $fields[4];
        $ref_allele = '' if ($ref_allele eq '*');
        my $alt_allele = $fields[5];
        $alt_allele = '' if ($alt_allele eq '*');

        my $normal_alt_freq = $fields[9];
        my $ref_alt_freq = $fields[10];
        my $qvalue = $fields[11];

        next if ($Opt{maxnorm} && $normal_alt_freq > $Opt{maxnorm});
        next if ($Opt{maxq} && $qvalue > $Opt{maxq});

        # widen if necessary:

        my $chrseq = $fasta_db->seq($chr);
        my $ref_end = length($chrseq);

        my $flank_bases_to_include = 1000;

        my $lf_offset = ($lfe - $flank_bases_to_include >= 0) ? $lfe - $flank_bases_to_include : 0;
        my $lf_include = $lfe - $lf_offset;
        my $rf_include = ($ref_end - $rfs + 1 < $flank_bases_to_include) ? $ref_end - $rfs + 1 : $flank_bases_to_include;
        my $left_flank_seq = uc substr($chrseq, $lf_offset, $lf_include);
        my $right_flank_seq = uc substr($chrseq, $rfs - 1, $rf_include);

        my $type = (length($ref_allele) > length($alt_allele)) ? 'Deletion' : 'Insertion';
        my $poly = GTB::Var::Polymorphism->new(
                                 -left_flank_end => $lfe,
                                 -right_flank_start => $rfs,
                                 -left_flank_seq => $left_flank_seq,
                                 -right_flank_seq => $right_flank_seq,
                                 -type => $type,
                                 -allele_seqs => [$ref_allele, $alt_allele] );


        my $new_left_flank_end = $poly->left_flank_end();
        my $new_right_flank_start = $poly->right_flank_start();
        my $new_ra_alleles = $poly->allele_seqs();
        my $new_ref_allele = $new_ra_alleles->[0];
        my $new_alt_allele = $new_ra_alleles->[1];

        my $ra_norm_alleles = $normal_slice->fetch_alleles(
                                 -chromosome => $chr,
                                 -start => $new_left_flank_end,
                                 -end => $new_right_flank_start);

        my $rh_norm_allele_length_counts = {}; # will only count alleles of different lengths
        foreach my $rh_align (@{$ra_norm_alleles}) {
            my $span = $rh_align->{'span'};
            next if (!$span);

            my $allele = $rh_align->{'read_seq'};
            $rh_norm_allele_length_counts->{length($allele)}++;
        }
        my $norm_ref_count = $rh_norm_allele_length_counts->{length($new_ref_allele)} || 0;
        my $norm_alt_count = $rh_norm_allele_length_counts->{length($new_alt_allele)} || 0;
        next if ($norm_alt_count + $norm_ref_count == 0); # no spanning reads
        next if ($Opt{'maxnorm'} && $norm_alt_count > $Opt{maxnorm} * ($norm_alt_count + $norm_ref_count));

        my $ra_tumor_alleles = $tumor_slice->fetch_alleles(
                                 -chromosome => $chr,
                                 -start => $new_left_flank_end,
                                 -end => $new_right_flank_start);

        my $rh_tumor_allele_length_counts = {}; # will only count alleles of different lengths
        foreach my $rh_align (@{$ra_tumor_alleles}) {
            my $span = $rh_align->{'span'};
            next if (!$span);

            my $allele = $rh_align->{'read_seq'};
            $rh_tumor_allele_length_counts->{length($allele)}++;
        }
        my $tumor_ref_count = $rh_tumor_allele_length_counts->{length($new_ref_allele)} || 0;
        my $tumor_alt_count = $rh_tumor_allele_length_counts->{length($new_alt_allele)} || 0;
        next if (($Opt{'min_som_reads'} && $tumor_alt_count + $norm_alt_count < $Opt{'min_som_reads'}) || 
                             ($tumor_alt_count + $tumor_ref_count == 0)); # not enough spanning reads

        # recalculate and print:
        $fields[9] = $norm_alt_count/($norm_alt_count + $norm_ref_count);
        $fields[10] = $tumor_alt_count/($tumor_alt_count + $tumor_ref_count);
        my $newline = join "\t", @fields;
        print "$newline\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=item B<--ref>

Specify the reference fasta file that was used in the shimmer run that resulted
in the indel VarSifter file.

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
