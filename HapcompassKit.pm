# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::HapcompassKit

=head1 SYNOPSIS

Running Hapcompass and analyze it's output

=head1 DESCRIPTION



=head1 FEEDBACK

=head2 Support

Please send you questions or bug reports to Email:

I<lufuhao@gmail.com>

=head1 AUTHORS - Fu-Hao Lu

Email: lufuhao@gmail.com (Always)
       Fu-Hao.Lu@jic.ac.uk (2012-2016)

=head1 CONTRIBUTORS

None

=head1 APPENDIX

Fight against Bioinformatics with Perl ^_^

=cut

#Coding starts

package FuhaoPerl5Lib::HapcompassKit;
use strict;
use warnings;
use Exporter;
use FuhaoPerl5Lib::CmdKit qw(exec_cmd_return);
use FuhaoPerl5Lib::VcfKit;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '20150603';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 Both    => [qw()]);



my $HapcompassKit_success=1;
my $HapcompassKit_failure=0;
mu $HapcompassKit_debug=0;



### RunHapCompass
### &RunHapCompass()
### Global: $path_java, $path_hapcompassjar, $path_hc2vcf
### Require: hapcompass.jar hc2vcf.jar
### Dependancy: freebayes
### Note:
sub RunHapCompass {
	my ($RHCref, $RHCbam, $RHCvcf, $RHCploidy, $RHCoutput_vcf_prefix)=@_;
	
	my $RHsubinfo='SUB(HapcompassKit::RunHapCompass)';
##COMMIT: check input files for HapCompass
	unless (-s $RHCref and -s $RHCbam and -s $RHCvcf) {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass inputs error: \n";
		return $HapcompassKit_failure;
	}
	unless ($RHCploidy>=2) {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass ploidy error: \n";
		return $HapcompassKit_failure;
	}
	if (! defined $RHCoutput_vcf_prefix or $RHCoutput_vcf_prefix eq '') {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass output prefix error: \n";
		return $HapcompassKit_failure;
	}
##COMMIT: Run HapCompass, output: $RHCoutput_vcf_prefix
##									_frags.txt
##									_MWER_solution.txt
##									_phasedSolution.txt
##									_reads.sam
##									_reduced_representation.sam
##									_reduced_representation.vcf
	my $RHCcmd_hapcompass='';
#	$RHCcmd_hapcompass="$path_java -jar $path_hapcompassjar --reference $RHCref --bam $RHCbam --vcf $RHCvcf --ploidy $RHCploidy --output $RHCoutput_vcf_prefix";###HapCompass does not need reference
	$RHCcmd_hapcompass="$path_java -jar $path_hapcompassjar --bam $RHCbam --vcf $RHCvcf --ploidy $RHCploidy --output $RHCoutput_vcf_prefix";
	if (! &exec_cmd_return($RHCcmd_hapcompass)) {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass running error\n";
		return $HapcompassKit_failure;
	}
##COMMIT: check HapCompass output
	if (! -s $RHCoutput_vcf_prefix.'_MWER_solution.txt') {
		print STDERR "SUB(RunHapCompass)Reference: HapCompass output *_MWER_solution.txt error\n";
		return $HapcompassKit_failure;
	}
##COMMIT: run hc2vcf to convert hapcompass output to vcf
##			output $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf'
	my $RHCcmd_hc2vcf='';
	$RHCcmd_hc2vcf="$path_java -jar $path_hc2vcf ${RHCoutput_vcf_prefix}_MWER_solution.txt $RHCvcf $RHCploidy";
	if (! &exec_cmd_return($RHCcmd_hc2vcf)) {
		print STDERR "SUB(RunHapCompass)Reference: hc2vcf running error\n";
		return $HapcompassKit_failure;
	}
	if (! -s $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf') {
		print STDERR "SUB(RunHapCompass)Reference: hc2vcf running error\n";
		return $HapcompassKit_failure;
	}
	else {
		return ($HapcompassKit_success, $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf');
	}
}



###
sub ReadCompassResults {
	my ($RCRfixallele_hashindex, $RCRploidy, $RCRcompassout, $RCHgeno_ploidy)=@_;
	
	local *HAPCOMPASSOUT;
	my %fixallele_bef=%{$RCRfixallele_hashindex}; undef $RCRfixallele_hashindex;
	my %fixallele_aft=();
	my $SCHploidy=length($RCHgeno_ploidy);
	if ($SCHploidy !~ m/^[123]{1}$/) {
		return $HapcompassKit_failure;
	}
	if (! -s $RCRcompassout) {
		print STDERR "SUB(ReadCompassResults)Error: hapcompass output not exist\n";
		return $HapcompassKit_failure;
	}
	close HAPCOMPASSOUT if (defined fileno(HAPCOMPASSOUT));
	unless (open (HAPCOMPASSOUT, "<$RCRcompassout") ) {
		return $HapcompassKit_failure;
	}
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	exampleBAM.bam
#chr1	1	rs5766699	A	G	99	PASS	.	GT:GQ	0/1:100
	my ($RCRtest_readvcf, $RCRreadvcf_hashindex)=&ReadVcf($RCRcompassout);
	if ($RCRtest_readvcf) {
		return $HapcompassKit_failure;
	}
	my %RCRreadvcf=%{$RCRreadvcf_hashindex};undef $RCRreadvcf_hashindex;
	foreach my $RCRchrom (keys %RCRreadvcf) {
		###PAUSE
	}
	close HAPCOMPASSOUT;
	return $HapcompassKit_success;
}



1;
