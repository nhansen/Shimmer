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

plan tests => 7;

my $out;
my $testref = 't/testref.fa';
my $testbam1 = 't/testbam1.bam';
my $testbam2 = 't/testbam2.bam';
my $testbedfile = 't/testbed.bed';
$ENV{PATH} = "./printCompCounts:$ENV{PATH}";
system("perl -w -I lib $script --ref $testref --outdir t/test1out $testbam1 $testbam2 > t/calltest1.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/test1out/som_counts.txt`;
like $out, qr/43/, "$script count";
$out = `awk '\$1==1 {print \$3}' t/test1out/somatic_diffs.vs`;
like $out, qr/11589021/, "$script variant";
system("rm -rf t/test1out");
system("perl -w -I lib $script --ref $testref --outdir t/test2out $testbam1 $testbam2 --mapqual 40 > t/calltest2.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/test2out/som_counts.txt`;
like $out, qr/11/, "$script map qual count";
system("rm -rf t/test2out");
system("perl -w -I lib $script --ref $testref --outdir t/test3out $testbam1 $testbam2 --minqual 30 > t/calltest3.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/test3out/som_counts.txt`;
# altered test output condition because github version of samtools is lowering some quality scores to 0
like $out, qr/^4[01]$/, "$script base qual count";
$out = `awk '\$2==120396876 {print \$4, \$5, \$6}' t/test3out/somatic_diffs.vcf`;
like $out, qr/^C\sT\s4(69|70)$/, "$script vcf file";
system("rm -rf t/test3out");
system("printCompCounts -fasta $testref -bam1 $testbam1 -bam2 $testbam2 -bedfile $testbedfile > t/calltest4.out 2>&1");
$out = `grep -v 'mpileup' t/calltest4.out | grep -v 'bam_header' | wc -l`;
like $out, qr/^20\s*/, "$script bedfile1";
system("perl -w -I lib $script --ref $testref --outdir t/test4out --bedfile $testbedfile $testbam1 $testbam2 > t/calltest5.out 2>&1");
$out = `awk '{print \$12}' t/test4out/som_counts.bh.txt`;
like $out, qr/e\-49/, "$script bedfile2";
system("rm -rf t/test4out");
