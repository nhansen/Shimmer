# $Id$
# t/02_call.t - check that shimmer successfully calls variants

use strict;
use Test::More;
use Module::Build;

my $script = 'shimmer.pl';

# Direct useless output to STDERR, to avoid confusing Test::Harness
my $stdin = select STDERR;
# Restore STDOUT as default filehandle
select $stdin;

plan tests => 4;

my $out;
my $testref = 't/testref.fa';
my $testbam1 = 't/testbam1.bam';
my $testbam2 = 't/testbam2.bam';
my $outdir = 't/testout';
$ENV{PATH} = "./printCompCounts:$ENV{PATH}";
system("perl -w -I lib $script --ref $testref --outdir $outdir $testbam1 $testbam2 > calltest.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/testout/som_counts.txt`;
like $out, qr/43/, "$script count";
$out = `awk '\$1==1 {print \$3}' t/testout/somatic_diffs.vs`;
like $out, qr/11589021/, "$script variant";
system("rm -rf $outdir");
system("perl -w -I lib $script --ref $testref --outdir $outdir $testbam1 $testbam2 --mapqual 40 > calltest.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/testout/som_counts.txt`;
like $out, qr/11/, "$script map qual count";
system("rm -rf $outdir");
system("perl -w -I lib $script --ref $testref --outdir $outdir $testbam1 $testbam2 --minqual 30 > calltest.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/testout/som_counts.txt`;
like $out, qr/41/, "$script base qual count";
system("rm -rf $outdir");
