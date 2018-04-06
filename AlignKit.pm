# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::AlignKit

=head1 SYNOPSIS

Alignment operations

=head1 Requirements

Perl Modules:

    FuhaoPerl5Lib::BamKit qw/IndexBam/;



=head1 DESCRIPTION

=over 4

=item Bowtie2Index(reference.fasta, index_name, [$pathbowtie2build])

    * Run bowtie2-build to index reference sequences
    * Dependency: &exec_cmd_return
    * Note:
    * Return: 1=Success    0=Failure

=item Bowtie2p($index, R1.fq, R2.fq, $out.bam, $BPoptions_bowtie2p, $BPpath_bowtie2p, $BPpath_samtools)

    * Run bowtie2 and map read to reference, sort and index
    * Dependency: FuhaoPerl5Lib::BamKit=IndexBam
    * $BPoptions_bowtie2p default : 
                    -q --phred33 --threads 1 --maxins 800 \
                    --rg-id MySample1 --rg \"SM:My\" --rg \"PL:ILLUMINA\" --rg \"LB:Sample\"
    * Return: 1=Success    0=Failure

=back

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
use FuhaoPerl5Lib::BamKit qw/IndexBam/;
use Exporter;
use Cwd;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20161103';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(Bowtie2Index Bowtie2p MummerChainNet);
%EXPORT_TAGS = ( DEFAULT => [qw(Bowtie2Index Bowtie2p MummerChainNet)],
                 ALL    => [qw(Bowtie2Index Bowtie2p MummerChainNet)]);


my $AlignKit_success=1;
my $AlignKit_failure=0;
my $AlignKit_debug=0; #checge to 1 if debug


### run bowtie2-build to index reference sequences
### RunBowtie2Index(reference.fasta, index_name, [$pathbowtie2build])
### Global: $AlignKit_success=1; $AlignKit_failure=0;
### Dependency: &exec_cmd_return
### Note:
sub Bowtie2Index {
	my ($BIreference, $BIindex, $BIpath_bowtie2build)=@_;

	my $BIsubinfo='SUB(AlignKit::RunBowtie2Index)';
	$BIpath_bowtie2build='bowtie2-build' unless (defined $BIpath_bowtie2build);


	unless (defined $BIreference and $BIreference ne '' and -s $BIreference) {
		print STDERR $BIsubinfo, "Error: invalid fasta file for bowtie2-build\n";
		return $AlignKit_failure;
	}
	unless (defined $BIindex and $BIindex ne '') {
		print STDERR $BIsubinfo, "Error: invalid index name for bowtie2-build\n";
		return $AlignKit_failure;
	}

##clean existing index
	unlink "$BIindex.1.bt2" if (-e "$BIindex.1.bt2");
	unlink "$BIindex.2.bt2" if (-e "$BIindex.2.bt2");
	unlink "$BIindex.3.bt2" if (-e "$BIindex.3.bt2");
	unlink "$BIindex.4.bt2" if (-e "$BIindex.4.bt2");
	unlink "$BIindex.rev.1.bt2" if (-e "$BIindex.rev.1.bt2");
	unlink "$BIindex.rev.2.bt2" if (-e "$BIindex.rev.2.bt2");

## run bowtie2-build
	if (! &exec_cmd_return("$BIpath_bowtie2build -f $BIreference $BIindex")) {
		print STDERR $BIsubinfo, "Error: bowtie2-build running error\n";
		return $AlignKit_failure;
	}
	elsif (! -s "$BIindex.1.bt2" or ! -s "$BIindex.2.bt2" or ! -s "$BIindex.3.bt2" or ! -s "$BIindex.4.bt2" or ! -s "$BIindex.rev.1.bt2" or ! -s "$BIindex.rev.2.bt2") {
		print STDERR $BIsubinfo, "Error: bowtie2-build output error\n";
		return $AlignKit_failure;
	}

	return $AlignKit_success;
}



### Run bowtie2 and map read to reference, and index
### Bowtie2p($index, R1.fq, R2.fq, $out.bam, $BPoptions_bowtie2p, $BPpath_bowtie2p, $BPpath_samtools)
### Global: $AlignKit_success=1; $AlignKit_failure=0;
### Dependency: FuhaoPerl5Lib::BamKit=IndexBam
### Note:
sub Bowtie2p {
	my ($BPindex, $BPfq_pair1, $BPfq_pair2, $BPbamout, $BPoptions_bowtie2p, $BPpath_bowtie2p, $BPpath_samtools)=@_;

### Default
	my $BPsubinfo='SUB(AlignKit::Bowtie2p)';
	$BPpath_bowtie2p='bowtie2' unless (defined $BPpath_bowtie2p);
	$BPoptions_bowtie2p=' -q --phred33 --threads 1 --maxins 800 --rg-id MySample1 --rg \"SM:My\" --rg \"PL:ILLUMINA\" --rg \"LB:Sample\" ';
	$BPpath_samtools='samtools' unless (defined $BPpath_samtools);

### Check input and output
	unless (defined $BPfq_pair1 and -s $BPfq_pair1) {
		print STDERR $BPsubinfo, "Error: invalid R1 fastq for bowtie2\n";
		return $AlignKit_failure;
	}
	unless (defined $BPfq_pair2 and -s $BPfq_pair2) {
		print STDERR $BPsubinfo, "Error: invalid R2 fastq for bowtie2\n";
		return $AlignKit_failure;
	}
	(my $BPbamout_base=$BPbamout)=~s/\.\w+$//;
	unlink $BPbamout if (-e $BPbamout_base);
	unlink "$BPbamout.bai" if (-e "$BPbamout_base.bai");

### Running bowtie and sort
	if (! &exec_cmd_return("$BPpath_bowtie2p $BPoptions_bowtie2p -x $BPindex -1 $BPfq_pair1 -2 $BPfq_pair2 | $BPpath_samtools view -b -h -S - | $BPpath_samtools sort - $BPbamout_base" )) {
		print STDERR $BPsubinfo, "Error: bowtie2 running error\n";
		return $AlignKit_failure;
	}
	elsif (! -s "$BPbamout_base.bam") {
		print STDERR $BPsubinfo, "Error: bowtie2 output error\n";
		return $AlignKit_failure;
	}

### Index BAM
	unless (&IndexBam("$BPbamout_base.bam")) {
		print STDERR $BPsubinfo, "Error: BAM index error: \n";
		return $AlignKit_failure;
	}
	
	return $AlignKit_success;
}



### Chain/Net 13 columns output of show-coords in mummer
### MummerChainNet($coords_in, $coords_out, [$gap:1000], [$min_identity_percentage:0], [$min_alignlength:0])
### Dependency:
### Global:
### Return: 1=success, 0=failure
### Note: 
### Input is the output of `show-coords -c -r -l -T xxx.delta
###      13 column and the first 4 line (headers) would be ignored
### $min_identity_percentage: say 90, 95, etc
sub MummerChainNet {
	my ($MCNcoords_in, $MCNcoords_out, $MCNmax_allowed_gap, $MCNmin_identity, $MCNmin_alignlength)=@_;

	my $MCNsubinfo='SUB(AlignKit::MummerChainNet)';
	$MCNmax_allowed_gap=1000 unless (defined $MCNmax_allowed_gap);
	$MCNmin_identity=0 unless (defined $MCNmin_identity);
	$MCNmin_alignlength=0 unless (defined $MCNmin_alignlength);
	my $MCNlinenum=0;
	local *MCNCOORDSIN; local *MCNCOORDSOUT;
	
	unless (defined $MCNcoords_in and -s $MCNcoords_in) {
		print STDERR $MCNsubinfo, "Error: invalid coordinate input\n----Please use show-coords -c -r -l -T xxx.delta > xxx.coord to prepare your input\n";
		return $AlignKit_failure;
	}
	unless ($MCNcoords_out) {
		print STDERR $MCNsubinfo, "Error: invalid coordinate output\n";
		return $AlignKit_failure;
	}
	unlink $MCNcoords_out if (-e $MCNcoords_out);
	unless ($MCNmax_allowed_gap=~/^\d+$/) {
		print STDERR $MCNsubinfo, "Error: invalid MCNmax_allowed_gap, should be INT\n";
		return $AlignKit_failure;
	}
	unless ($MCNmin_identity=~/^\d+\.{0,1}\d*$/) {
		print STDERR $MCNsubinfo, "Error: invalid MCNmin_identity, should be INT/FLOAT\n";
		return $AlignKit_failure;
	}
	unless ($MCNmin_alignlength=~/^\d+$/) {
		print STDERR $MCNsubinfo, "Error: invalid MCNmin_alignlength, should be INT\n";
		return $AlignKit_failure;
	}
	
	close MCNCOORDSIN if (defined fileno(MCNCOORDSIN));
	unless (open MCNCOORDSIN, "<", $MCNcoords_in) {
		print STDERR $MCNsubinfo, "Error: can not open coordinate input\n";
		return $AlignKit_failure;
	}
	close MCNCOORDSOUT if (defined fileno(MCNCOORDSOUT));
	unless (open MCNCOORDSOUT, ">", $MCNcoords_out) {
		print STDERR $MCNsubinfo, "Error: can not write coordinate output\n";
		return $AlignKit_failure;
	}
	my $MCNline=<MCNCOORDSIN>; print MCNCOORDSOUT $MCNline;$MCNlinenum++;### ignore first 4 header lines
	$MCNline=<MCNCOORDSIN>; print MCNCOORDSOUT $MCNline;$MCNlinenum++;### ignore first 4 header lines
	$MCNline=<MCNCOORDSIN>; print MCNCOORDSOUT $MCNline;$MCNlinenum++;### ignore first 4 header lines
	$MCNline=<MCNCOORDSIN>; print MCNCOORDSOUT $MCNline;$MCNlinenum++;### ignore first 4 header lines
	while ($MCNline=<MCNCOORDSIN>) {
		chomp $MCNline;
		$MCNlinenum++;
		my @MCNarr=split(/\t/, $MCNline);
		unless (scalar(@MCNarr)==13) {
			print STDERR $MCNsubinfo, "Error: NOT 13 column at line($MCNlinenum): $MCNcoords_in\n";
			return $AlignKit_failure;
		}
		### to be continued
	}
	close MCNCOORDSIN;
	
	
	
	
	close MCNCOORDSOUT;
	return $AlignKit_success;
}


#my $BIsubinfo='SUB(AlignKit::RunBowtie2Index)';
#my $AlignKit_success=1;
#my $AlignKit_failure=0;
1;
__END__
