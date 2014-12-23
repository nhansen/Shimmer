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

filter_messy_indels.pl - script to examine reads from the normal sample in regions of predicted somatic indels to remove presumed artifacts of alignment

=head1 SYNOPSIS

  filter_messy_indels.pl --indels <shimmer indel varsifter> --normalbam <BAM file from normal sample>

=head1 DESCRIPTION

This script reads in a Shimmer Varsifter file with indel calls, and examines the normal BAM file to look for reads with evidence of poor alignment.

=cut

#------------
# Begin MAIN
#------------

process_commandline();
my $indel_file = $Opt{'indels'};
my $bamfile = $Opt{'normalbam'};
my $maxmessy = $Opt{'maxmessy'};

filter_messy_indels($indel_file, $bamfile);

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( maxmessy => 3 );
    GetOptions(
        \%Opt, qw(
            manual help+ version maxmessy=i
            verbose indels=s normalbam=s
		  )
	) || pod2usage(0);
	if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
	if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
	if ( $Opt{version} ) {
		die "filter_messy_indels.pl, ", q$Revision: 4129 $, "\n";
	}
	if ( !$Opt{'indels'} || !$Opt{'normalbam'} ) {
		die "filter_messy_indels.pl --normalbam <bam file for normal sample> --indels <indel varsifter file>\n";
	}

}

sub filter_messy_indels {
    my $shimmer_indel_file = shift;
    my $normal_bam = shift;

    my $indel_fh = Open($shimmer_indel_file);
    print STDERR "Chr\tLeftFlank\tRightFlank\tPairUnaligned\tInserted\tDeleted\tClippedNoIndel\tUnremarkable\n";
    while (<$indel_fh>) {
        if (/^Index/) {
            print;
            next;
        }
        my $line = $_;
        chomp;
        my @fields = split /\t/, $_;
        my $chr = $fields[1];
        my $lfe = $fields[2];
        my $rfs = $fields[3];

        open SAM, "samtools view $normal_bam $chr:$lfe-$rfs |"
            or die "Couldn\'t run samtools view $normal_bam $chr:$lfe-$rfs!\n";

        my ($unaligned_pair, $inserted, $deleted, $clipped_noindel, $unremarkable) = (0, 0, 0, 0, 0);
        while (<SAM>) {
            chomp;
            my @fields = split /\t/, $_;
            my $cigar = $fields[5];
            if ($cigar eq '*') { # unaligned pair
                 $unaligned_pair++;
            }
            elsif ($cigar =~ /I/) {
                 $inserted++;
            }
            elsif ($cigar =~ /D/) {
                 $deleted++;
            }
            elsif ($cigar =~ /S/) {
                 $clipped_noindel++;
            }
            else {
                 $unremarkable++;
            }
        }
        close SAM;

        next if ($unaligned_pair+$inserted+$deleted >= $maxmessy);
        print STDERR "$chr\t$lfe\t$rfs\t$unaligned_pair\t$inserted\t$deleted\t$clipped_noindel\t$unremarkable\n";

        # survived filtering--print
        print "$line";
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
