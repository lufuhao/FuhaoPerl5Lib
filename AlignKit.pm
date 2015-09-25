# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::AlignKit

=head1 SYNOPSIS

File operations

=head1 Requirements



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
package FuhaoPerl5Lib::AlignKit;
use strict;
use warnings;
#use FuhaoPerl5Lib::CmdKit;
#use FuhaoPerl5Lib::FileKit;
use Exporter;
use Cwd;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150617';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 ALL    => [qw()]);


my $AlignKit_success=1;
my $AlignKit_failure=0;
my $AlignKit_debug=0; #checge to 1 if debug


### run bowtie2-build to index reference sequences
### RunBowtie2Index(reference.fasta, index_name)
### Global: $pathbowtie2build
### Dependency: &exec_cmd_return
### Note:
sub RunBowtie2Index {
	my ($RBIreference, $RBIindex)=@_;
	unless (defined $RBIreference and $RBIreference ne '' and -s $RBIreference) {
		print STDERR "SUB(RunBowtie2Index)Error: invalid fasta file for bowtie2-build\n";
		return 1;
	}
	unless (defined $RBIindex and $RBIindex ne '') {
		print STDERR "SUB(RunBowtie2Index)Error: invalid index name for bowtie2-build\n";
		return 1;
	}
##clean existing index
	unlink "$RBIindex.1.bt2" if (-e "$RBIindex.1.bt2");
	unlink "$RBIindex.2.bt2" if (-e "$RBIindex.2.bt2");
	unlink "$RBIindex.3.bt2" if (-e "$RBIindex.3.bt2");
	unlink "$RBIindex.4.bt2" if (-e "$RBIindex.4.bt2");
	unlink "$RBIindex.rev.1.bt2" if (-e "$RBIindex.rev.1.bt2");
	unlink "$RBIindex.rev.2.bt2" if (-e "$RBIindex.rev.2.bt2");
## run bowtie2-build
	if (! &exec_cmd_return("$pathbowtie2build -f $RBIreference $RBIindex")) {
		print STDERR "SUB(RunBowtie2Index)Error: bowtie2-build running error\n";
		return 1;
	}
	elsif (! -s "$RBIindex.1.bt2" or ! -s "$RBIindex.2.bt2" or ! -s "$RBIindex.3.bt2" or ! -s "$RBIindex.4.bt2" or ! -s "$RBIindex.rev.1.bt2" or ! -s "$RBIindex.rev.2.bt2") {
		print STDERR "SUB(RunBowtie2Index)Error: bowtie2-build output error\n";
		return 1;
	}
	return 0;
}


### Run bowtie2 and map read to reference, and index
### RunBowtie2p(reference, fq)
### Global: $pathbowtie2p, $fq_qual_score, $numthreads
### Dependency: 
### Note: $fq_qual_score='phred33'; $maxinsert=800
sub RunBowtie2p {
	my ($RBPindex, $RBPreadgroup, $RBPrep, $RBPfq_pair1, $RBPfq_pair2, $RBPbamout)=@_;
	unless (defined $RBPfq_pair1 and -s $RBPfq_pair1 and defined $RBPfq_pair2 and -s $RBPfq_pair2) {
		print STDERR "SUB(RunBowtie2p)Error: invalid bowtie2 R1 and R2\n";
		return 1;
	}
	(my $RBPbamout_base=$RBPbamout)=~s/\.\w+$//;
#Par${tissue}${RBPrep}_rep${RBPrep}.vs.$bt2index.sort
	unlink $RBPbamout if (-e $RBPbamout_base);
	unlink "$RBPbamout.bai" if (-e "$RBPbamout_base.bai");
	if (! &exec_cmd_return("$pathbowtie2p -q --$fq_qual_score --threads $numthreads --maxins $maxinsert --rg-id Par${RBPreadgroup}${RBPrep} --rg \"SM:Par${RBPreadgroup}\" --rg 'PL:ILLUMINA' --rg \"LB:Par${RBPreadgroup}\" -x $RunDir/4.mapping/$RBPindex -1 $RBPfq_pair1 -2 $RBPfq_pair2 | $path_samtools view -b -h -S - | $path_samtools sort - $RBPbamout_base" )) {
		print STDERR "SUB(RunBowtie2p)Error: bowtie2 running error\n";
		return 1;
	}
	elsif (! -s "Par${RBPreadgroup}${RBPrep}_rep${RBPrep}.vs.$RBPindex.sort.bam") {
		print STDERR "SUB(RunBowtie2p)Error: bowtie2 output error\n";
		return 1;
	}
	if (! &IndexBam("Par${RBPreadgroup}${RBPrep}_rep${RBPrep}.vs.$RBPindex.sort.bam")) {
		print STDERR "SUB(RunBowtie2p)Error: BAM index error\n";
		return 1;
	}
	return 0;
}

1;
