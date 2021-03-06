# $Id:$
# To create build script, run, "perl Build.PL"
use Module::Build;

my $class = Module::Build->subclass(
	class => 'My::Build',
	code => q{
		sub ACTION_code {
			my $self = shift;
                        $self->SUPER::ACTION_code();
			$self->do_system(qw(make));
		}
		} );

my $build = $class->new(
        module_name   => 'shimmer',
        dist_version  => 0.10,
        dist_author   => 'GTB Perl Programmers',
        dist_abstract => 'Script for calling somatic mutations from matched tumor and normal BAM files',
        license       => 'perl',
        scripts       => [qw(shimmer.pl)],
        requires      => { },
        );
if (!`which samtools`) {
    warn "ERROR: samtools executable is not found.\nThis script requires "
        . "the samtools binary to be in the user's PATH\n";
}
if (!`which R`) {
    warn "ERROR: R executable is not found.\nThis script requires "
        . "the R executable to be in the user's PATH\n";
}
my $statmod = `R --file=check_statmod.R`;
if ($statmod !~ /TRUE/) {
    warn "ERROR: You have R installed, but it doesn't seem to have the statmod library.\n"
        . "Shimmer requires the R statmod library to be installed.\n";
}
$build->create_build_script;

