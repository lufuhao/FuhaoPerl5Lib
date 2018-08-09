# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::BamKit

=head1 SYNOPSIS

  use samtools to index, extract, merge, rehead BAM files

=head1 Requirements

Perl Modules:
        Cwd;
        Bio::DB::Sam;
        FuhaoPerl5Lib::CmdKit;
        FuhaoPerl5Lib::FileKit;
        FuhaoPerl5Lib::MiscKit;

=head1 DESCRIPTION

=over 2

=item Bam2FastQ ($bamin, $fastqout, map_code, MAPQ_code, [path_samtools])

    * Convert bam files into fastq
    * Map_code: 0=all    1=mapped_only    2=unmapped_only
    * MAPQ_code: any alignment less than MAPQ would be ignored
    * Return: 1=Sucesss    0=Failure

=item Bam2FastqProg($BAM, $fastq_prefix, $path_bam2fastq)

    * Convert BAM to FASTQ using bam2fastq program (https://github.com/jts/bam2fastq)
    * Return (1/0, \%hash=('R1' => $fastqR1, 'R2' => $fastqR2, 'M' => $fastqUnpaired)

=item BamAtacShift($in.bam, $out.bam)

    * Modify BAM file to adjust 9 bp for ATAC-seq: forward strand start +4 and reverse strans start -5
    * Require: samtools
    * Output: Default BAMin.ATAC9shift.bam if not specified
    * Return: 1=Success    0= Failure

=item BamExtractReadsUsingBed($inputbam, $file_region, $outreadname[, $path_samtools])

    * Extract reads in a BAM file using specified region
    * $file_region: region file format: 
            seq_ID\tstart\tend
            seqID
            seqID:Num1-Num2
    * Note: bed format is half open, start is less 1 than real start
    * Return: 1=Sucesss    0=Failure

=item BamFilterReadsByNames($inputbam, $file_readnames, $filter_code, $outbam[, $path_samtools])

    * Filter BAM reads by a file of read name list
    * $file_readnames: per readname per line
    * $filter_code: 1=filter    0=filter_out
    * Return: 1=Sucesss    0=Failure

=item BamKeepBothMates ($input.R1.bam, $input.R2.bam, $paired.R1.bam, $paired.R2.bam, [$unpaired.R1.bam], [$unpaired.R2.bam], [path_samtools])

    * Keep alignments that have both mates mapped
    * Return: 1=Success    0=Failure

=item BamKeepNumAlignemnts(In.name-sorted.bam, $k [default:1], out.bam)

    * Remove alignment >k mapping times from BAM
    * Require: samtools
    * Note: need to double times for paired end, so k2 means unique for paired reads
    * Return: 1=Success    0=Failure

=item BamRestoreSplit($bamin, $bed, $bamout, [$path_samtools])

    * Restore BAM coordintes due to splited reference
    * Return: 1=Success    0= Failure
    * BED file 6 columns: 
        chr1A_part1     0       471304005       chr1A   0       471304005
        chr1A_part2     0       122798051       chr1A   471304005       594102056
        chr1B_part1     0       438720154       chr1B   0       438720154
        chr1B_part2     0       251131716       chr1B   438720154       689851870
        chr1D_part1     0       452179604       chr1D   0       452179604
        chr1D_part2     0       43273582        chr1D   452179604       495453186

=item CalCigarRefLength (CIGAR)

    * Calculate reference length based cigar
    * Return the reference length of that BAM alignment

=item CalCigarReadLength (CIGAR)

    * Calculate read length based cigar
    * Return the read length of that BAM alignment

=item CalcFPKM(BAMin, $total_mapped_reads, readgroup aware(0/1), [$CFpath_samtools])

=item ExpressFpkm($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, \@EFsamfiles, \@seq_ids, [path_express])

    * For special

=item ExtactBam($input.bam, $seqid_string, $out.bam, [$samtools_path])

    * Return: 1=Sucesss    0=Failure

=item IndexBam($baminput, $path_samtools)

    * Return: 1=Sucesss    0=Failure

=item BamMarkPairs ($input.R1.bam, $input.R2.bam, $output.R1.bam, $output.R2.bam, [path_samtools])

    * Mark separately mapped FLAG to R1 67 R2 131
    * Return: 1=Success    0=Failure

=item ReadSam($file.sam,$ref,fa, ReturnCode(1/2/3))

    * Dependancy: Bio::DB::Sam
    * ReturnCode: 1=Bio::DB::Sam objective; 2=reference name array, 3=to be defined
    * Return: (0/1, $CodeAssigned)

=item SamCleanHeader($input.bam, $out.bam, [out.noRG.bam], [\@seqids], [path2samtools])

    * Return: 1=Sucesss    0=Failure
    * [\@seqids] could be a reference to a constant or array
    * [\@seqids] could be a string delimited by space or comma

=item SortBam($input.bam, $out.bam, [$samtools_path], [$samtools_sort_options])

    * 1=Sucesss    0=Failure
    * samtools_sort_options: -n, -l Num, -@ Num

=item SplitCigar(SamCigar)

    * Return: double array ((2, M), (5, D), ....)

=item VerifyCigarLength ($CIGAR, $READLENGTH)

    * Verify Cigar length == readlength
    * Return: 1=equal    0=NOTequal

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

package FuhaoPerl5Lib::BamKit;
use strict;
use warnings;
use Cwd;
use Exporter;
use FuhaoPerl5Lib::CmdKit qw(exec_cmd_return);
use Bio::DB::Sam;
use FuhaoPerl5Lib::FileKit;
use FuhaoPerl5Lib::MiscKit qw(IsReference);
use Data::Dumper qw /Dumper/;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '20180809';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(IndexBam ExtactBam SplitCigar SamCleanHeader Bam2FastQ SortBam CalcFPKM ReduceReadNameLength ReadSam Bam2FastqProg VerifyCigarLength CalCigarRefLength BamFilterReadsByNames BamExtractReadsUsingBed BamKeepBothMates BamMarkPairs BamKeepNumAlignemnts BamAtacShift BamRestoreSplit);
%EXPORT_TAGS = ( DEFAULT => [qw(IndexBam ExtactBam SplitCigar SamCleanHeader Bam2FastQ SortBam CalcFPKM ReduceReadNameLength ReadSam Bam2FastqProg VerifyCigarLength CalCigarRefLength BamFilterReadsByNames BamExtractReadsUsingBed BamKeepBothMates BamMarkPairs BamKeepNumAlignemnts BamAtacShift BamRestoreSplit)],
                 Both    => [qw(IndexBam ExtactBam SplitCigar SamCleanHeader Bam2FastQ SortBam CalcFPKM ReduceReadNameLength ReadSam Bam2FastqProg VerifyCigarLength CalCigarRefLength BamFilterReadsByNames BamExtractReadsUsingBed BamKeepBothMates BamMarkPairs BamKeepNumAlignemnts BamAtacShift BamRestoreSplit)]);



my $BamKit_success=1;
my $BamKit_failure=0;
my $BamKit_debug=0;


### samtools index sorted BAM
### IndexBam($baminput, $path_samtools)
### Global: $BamKit_success; $BamKit_failure;
### Dependency: &exec_cmd_return 
### Note:
sub IndexBam {
	my ($IBbamin, $IBpath_samtools)=@_;
	
	my $CFsubinfo='SUB(BamKit::IndexBam)';
	$IBpath_samtools='samtools' unless (defined $IBpath_samtools);
	
	unless (defined $IBbamin and -s $IBbamin) {
		print STDERR $CFsubinfo, "Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unlink "$IBbamin.bai" if (-e "$IBbamin.bai");
	
	if (! &exec_cmd_return("$IBpath_samtools index $IBbamin")) {
		print STDERR $CFsubinfo, "Error: SAMtools index running error\n";
		return $BamKit_failure;
	}
	elsif (! -s "$IBbamin.bai") {
		print STDERR $CFsubinfo, "Error: SAMtools index output error\n";
		return $BamKit_failure;
	}
	return $BamKit_success;
}



### Extract subsets from BAM
### &ExtactBam($input.bam, $seqid_string, $out.bam, [$samtools_path])
### Global:
### Dependency: &exec_cmd_return
### Note:
sub ExtractBam {
	my ($EBbamin, $EBseqids, $EBbamout, $EBpath_samtools)=@_;
	
	my $EBsubinfo='SUB(BamKit::ExtractBam)';
	$EBpath_samtools='samtools' unless (defined $EBpath_samtools);
	my $EBcmd="";
	
	unless (defined $EBbamin and -s $EBbamin) {
		print STDERR $EBsubinfo, "Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $EBseqids and $EBseqids !~/^\s*$/) {
		print STDERR $EBsubinfo, "Error: invalid seqids\n";
		return $BamKit_failure;
	}
	unless (defined $EBbamout) {
		print STDERR $EBsubinfo, "Error: invalid BAM output\n";
		return $BamKit_failure;
	}
	unlink $EBbamout if (-e $EBbamout);
	
	if ($EBbamin=~/\.sam$/i) {
		$EBcmd="$EBpath_samtools view -b -h -S ";
	}
	elsif ($EBbamin=~/\.bam$/i) {
		$EBcmd="$EBpath_samtools view -b -h ";
	}
	else {
		print STDERR $EBsubinfo, "Error: unknown input BAM/SAM extension \n";
		return $BamKit_failure;
	}
	unless (-s "$EBbamin.bai") {
		unless (&IndexBam($EBbamin)) {
			print STDERR $EBsubinfo, "Error: can not Index BAM: $EBbamin\n";
			return $BamKit_failure;
		}
	}
	if (! &exec_cmd_return("$EBcmd $EBbamin $EBseqids > $EBbamout")) {
		print STDERR $EBsubinfo, "Error: BAM extract failed\n";
		return $BamKit_failure;
	}
	else {
		return $BamKit_success;
	}
}



### Sort bams
### &SortBam($input.bam, $out.bam, [$samtools_path], [$samtools_sort_options])
### Global: $BamKit_failure;
### Dependency: &exec_cmd_return
### Note:
sub SortBam {
	my ($SBbamin, $SBbamout, $SBpath_samtools, $SBsamtools_options)=@_;
	
	my $SBsubinfo='SUB(BamKit::SortBam)';
	$SBpath_samtools='samtools' unless (defined $SBpath_samtools);
	$SBsamtools_options=' ' unless (defined $SBsamtools_options);
	(my $SBbamout_prefix=$SBbamout)=~s/\.bam$//i;
	my $SBcmd='';
	
	unless (defined $SBbamin and -s $SBbamin) {
		print STDERR $SBsubinfo, "Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $SBbamout and $SBbamout=~/^\S+$/) {
		print STDERR $SBsubinfo, "Error: invalid BAM output\n";
		return $BamKit_failure;
	}
	else {
		unlink $SBbamout if (-e $SBbamout);
		unlink "$SBbamout.bai" if (-e "$SBbamout.bai");
	}
	
	if ($SBbamin=~/\.sam$/i) {
		$SBcmd="$SBpath_samtools view -b -h -S $SBbamin | $SBpath_samtools sort $SBsamtools_options - $SBbamout_prefix";
	}
	elsif ($SBbamin=~/\.bam$/i) {
		$SBcmd="$SBpath_samtools sort $SBsamtools_options $SBbamin $SBbamout_prefix";
	}
	else {
		print STDERR ${SBsubinfo}, "Error: unknown input BAM/SAM extension \n";
		return $BamKit_failure;
	}
#	if (! -s "$SBbamin.bai") { ### Index BAM
#		if (! &IndexBam($SBbamin)) {
#			print STDERR "${SBsubinfo}Error: can not Index BAM: $SBbamin\n";
#			return $BamKit_failure;
#		}
#	}
	unless (exec_cmd_return("$SBcmd")) {
		print STDERR $SBsubinfo, "Error: BAM sort failed\n";
		return $BamKit_failure;
	}
	return $BamKit_success;
}





### extract batch of sequence alignment from Bams, and index
### &ExtactBam(bam_arr, sequence_arr, out.bam, [out.no_RG.bam], [samtools_path])
### Global:
### Dependency: bam_clean_header.pl &exec_cmd_return
### Note:
sub ExtractBam_old {
	my ($EBfiles_bam_index, $EBseq_obj, $EBoutput, $EBoutnogroup, $EBpath_samtools)=@_;

	my $CFsubinfo='SUB(BamKit::ExtractBam)';
	$EBpath_samtools='samtools' unless (defined $EBpath_samtools);
	my $EBi=0;
	my @EBtemp_bams=();
	my $EBsam_seq='';
	my $test_outnogroup=0;
	
	if (scalar(@{$EBseq_obj})<1) {###return error if no sequences
		print STDERR "${CFsubinfo}Error: fasta ID array empty to be extracted to extract\n";
		return $BamKit_failure;
	}
	else {
		$EBsam_seq=join (' ', @{$EBseq_obj});
		if ($EBsam_seq=~/^\s*$/) {
			print STDERR "${CFsubinfo}Error: empty fasta ID list to extract\n";
			return $BamKit_failure;
		}
	}
	unlink $EBoutput if (-e $EBoutput);
	unlink "$EBoutput.bai" if (-e "$EBoutput.bai");
	if (defined $EBoutnogroup) {
		$test_outnogroup=1;
		unlink $EBoutnogroup if (-e $EBoutnogroup);
		unlink "$EBoutnogroup.bai" if ( -e "$EBoutnogroup.bai");
	}
	
##COMMENT: extract bam from each
#	print STDERR "SUB(ExtractBam)Info: sequence to be extracted: $EBsam_seq\n";###for test###
	foreach my $EBind_bam (@{$EBfiles_bam_index}) {
		$EBi++;
		my $EBcmd="$EBpath_samtools view -bh $EBind_bam $EBsam_seq > $EBoutput.$EBi.bam";
		if (! &exec_cmd_return($EBcmd)) {
			print STDERR "${CFsubinfo}Error: bam extract $EBsam_seq from $EBind_bam failed\n";
			return $BamKit_failure;
		}
		else {
			push (@EBtemp_bams, "$EBoutput.$EBi.bam");
		}
	}
##COMMENT: merge bams 
	my $EBmerge_bams=join(' ', @EBtemp_bams);
	my $EBcmd2="$EBpath_samtools merge $EBoutput.multihead.bam $EBmerge_bams";
	if (! &exec_cmd_return($EBcmd2)) {
		print STDERR "${CFsubinfo}Error: bam merge from $EBmerge_bams failed\n";
		return $BamKit_failure;
	}
	else {
		if (-s "$EBoutput.multihead.bam") {
			map {unlink($_)} @EBtemp_bams;###delete temporary files ### For test ###
		}
		else{
			print STDERR "${CFsubinfo}Error: bam merge from $EBmerge_bams not found\n";
			return $BamKit_failure;
		}
	}
	
	my $EBcmd4="bam_clean_header.pl -i $EBoutput.multihead.bam -o $EBoutput -n $EBoutnogroup";
	if (! &exec_cmd_return($EBcmd4)) {
		print STDERR "${CFsubinfo}Error: bam header clean from $EBmerge_bams failed\n";
		return $BamKit_failure;
	}
	else {
		if ( -s $EBoutput and -s $EBoutnogroup) {
			unlink ("$EBoutput.multihead.bam") if (-s "$EBoutput.multihead.bam");
			if (&IndexBam($EBoutput)) {
				print STDERR "${CFsubinfo}Error: bam index from $EBmerge_bams failed\n";
				return $BamKit_failure;
			}
			if (&IndexBam($EBoutnogroup)) {
				print STDERR "${CFsubinfo}Error: bam index2 from $EBmerge_bams failed\n";
				return $BamKit_failure;
			}
		}
		else{
			print STDERR "${CFsubinfo}Error: bam headerclean output from $EBmerge_bams not found\n";
			return $BamKit_failure;
		}
	}
	return $BamKit_success;
}



### Clean those sequences in header but not exist in alignments
### SamCleanHeader($input.bam, $out.bam, [out.noRG.bam], [\@seqids], [path2samtools])
### Global: 
### Dependency: 
### Note: 
sub SamCleanHeader {
	my ($SCHbam_input, $SCHbam_withrg, $SCHbam_norg, $SCHseq_obj, $SCHpathsamtools)=@_;

	local *SCHBAMINPUT1; local *SCHBAMINPUT2; local *SCHBAMOUTPUT;
	my $SCHsubinfo='SUB(BamKit::SamCleanHeader)';
	my $SCHnew_readgroup=0;
	my %SCHreahgroup=();
	my $BamKit_debug=0;
	$SCHpathsamtools='samtools' unless (defined $SCHpathsamtools);
	my %SCHseqID=();
	my $SCHgetseqids_from_bam=1;

	unless (defined $SCHbam_input and -s $SCHbam_input) {
		print STDERR $SCHsubinfo, "Error: BAM input file not found\n";
		return $BamKit_failure;
	}
	unless (defined $SCHbam_withrg and $SCHbam_withrg ne '' and $SCHbam_withrg !~ m/\s+/) {
		my $SCHbam_prefix=&RetrvNoExt($SCHbam_input);
		$SCHbam_withrg=$SCHbam_prefix.".CleanHeader.bam";
	}
	unlink ($SCHbam_withrg) if (-e $SCHbam_withrg);
	my $SCHtest_norg=0;
	if (defined $SCHbam_norg) {
		$SCHtest_norg=1;
		unlink $SCHbam_norg if (-e $SCHbam_norg);
		unlink "$SCHbam_norg.bai" if (-e "$SCHbam_norg.bai");
	}
	if (defined $SCHseq_obj) {
		my @SCHids=();
		my $SCHtest_id_reference=&IsReference($SCHseq_obj);
		if ($SCHtest_id_reference>0) {
			if ($SCHtest_id_reference==1) {
				push (@SCHids, ${$SCHseq_obj});
			}
			elsif ($SCHtest_id_reference==2) {
				@SCHids=@{$SCHseq_obj;};
			}
			else {
				print STDERR $SCHsubinfo, "Error: unknown seqids reference type\n";
				return $BamKit_failure;
			}
		}
		elsif ($SCHtest_id_reference==0) {
			@SCHids=split(/\s+|,/, $SCHseq_obj);
		}
#		undef $SCHseq_obj;
		foreach (@SCHids){
			if (/^\S+$/) {
				$SCHseqID{$_}++;
				$SCHgetseqids_from_bam=0;
			}
		}
	}

##COMMENT: read bam/sam
	close SCHBAMINPUT1 if (defined fileno(SCHBAMINPUT1));
	if ($SCHbam_input=~/\.sam$/i) {
		unless (open (SCHBAMINPUT1, "$SCHpathsamtools view -S $SCHbam_input|")) {
			print STDERR $SCHsubinfo, "Error: sam open error\n";
			return $BamKit_failure;
		}
	}
	elsif ($SCHbam_input=~/\.bam$/i) {
		unless (open (SCHBAMINPUT1, "$SCHpathsamtools view $SCHbam_input|")) {
			print STDERR $SCHsubinfo, "Error: bam open error\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $SCHsubinfo, "Error: unknown BAM format\n";
		return $BamKit_failure;
	}
	
	while (my $SCHline1=<SCHBAMINPUT1>) {
		my @SCHarr1=split(/\t/, $SCHline1);
		$SCHseqID{$SCHarr1[2]}++;
		if ($SCHline1=~/RG:Z:(\S*)/) {
#			print $1."\n";###For test###
			if (defined $1 and $1 ne '') {
				$SCHnew_readgroup=1;
				$SCHreahgroup{$1}++;
			}
		}
	}
	close SCHBAMINPUT1;
	print STDERR $SCHsubinfo, "Info: total of ".scalar(keys %SCHseqID)." REF sequences detected\n" if ($BamKit_debug);
	print STDERR $SCHsubinfo, "Info: total of ".scalar(keys %SCHreahgroup)." readgroups detected\n" if ($BamKit_debug);
	unless (scalar(keys %SCHseqID)>=1) {### For BAMs without alignments ###
		print STDERR $SCHsubinfo, "Error: no alignments\n";
		return 2;
	}
##COMMENT: write output
	if ($SCHbam_input=~/\.sam$/i) {
		unless (open (SCHBAMINPUT2, "$SCHpathsamtools view -h -S $SCHbam_input|")) {
			print STDERR $SCHsubinfo, "Error: sam open error\n";
			return $BamKit_failure;
		}
	}
	elsif ($SCHbam_input=~/\.bam$/i) {
		unless (open (SCHBAMINPUT2, "$SCHpathsamtools view -h $SCHbam_input|")) {
			print STDERR $SCHsubinfo, "Error: bam open error\n";
			return $BamKit_failure;
		}
	}
	unless (open (SCHBAMOUTPUT, "| $SCHpathsamtools view -bS - > $SCHbam_withrg")) {
		print STDERR $SCHsubinfo, "Error: bam write error\n";
		return $BamKit_failure;
	}
	if ($SCHtest_norg==1) {
		unless (open (SCHBAMOUTNOGROUP, "| $SCHpathsamtools view -bS - > $SCHbam_norg")) {
			print STDERR $SCHsubinfo, "Error: bam2 write error\n";
			return $BamKit_failure;
		}
	}
	
	while (my $SCHline=<SCHBAMINPUT2>) {
		if ($SCHline=~m/^\@/) {
##COMMENT: clean header
			if ($SCHline=~m/^\@SQ\s+SN:(\S+)\s+LN:\d+/) {
				if (! exists $SCHseqID{$1}) {
					print STDERR $SCHsubinfo, "Info: ignore $SCHline\n" if ($BamKit_debug);
					next;
				}
				print SCHBAMOUTPUT $SCHline;
				print SCHBAMOUTNOGROUP $SCHline if ($SCHtest_norg==1);
			}
##Remove RG in header
			elsif ($SCHline=~m/^\@RG/) {
				if (scalar(keys %SCHreahgroup) >0) {
					if ($SCHnew_readgroup==1) {
						foreach (keys %SCHreahgroup) {
							print SCHBAMOUTPUT '@RG'."\tID:$_\tSM:$_\tPL:Illumina\tLB:$_\n";
						}
						$SCHnew_readgroup=0;
					}
				}
				else {
					print SCHBAMOUTPUT $SCHline;
				}
				next;
			}
			else {
				print SCHBAMOUTPUT $SCHline;
				print SCHBAMOUTNOGROUP $SCHline if ($SCHtest_norg==1);
			}
		}
		else {
			print SCHBAMOUTPUT $SCHline;
			if ($SCHtest_norg==1) {
				$SCHline=~s/\tRG:Z:\S+//g;
				print SCHBAMOUTNOGROUP $SCHline;
			}
		}
	}
	close SCHBAMOUTPUT;
	close SCHBAMINPUT2;
	close SCHBAMOUTNOGROUP if (defined $SCHbam_norg);
	
	return $BamKit_success;
}



### split cigar into double array ((2, M), (5, D), ....)
### SplitCigar(SamCigar)
### Global:
### Dependancy:
### Note:
sub SplitCigar {
	my $SCcigar_string = shift;
	my @SCreturn_cigar_arr=();
	while ($SCcigar_string =~ /(\d+)([MIDNSHP=X])/g) {
		my @SCoperation=($1, $2);
		push @SCreturn_cigar_arr, \@SCoperation;
	}
	return \@SCreturn_cigar_arr;
}



### Verify Cigar length == readlength
### VerifyCigarLength ($CIGAR, $READLENGTH);
### Global: None
### Dependancy:
### Note:
sub VerifyCigarLength {
	my ($VCLcigar, $VCLlength)=@_;
	
	my $VCLsubinfo='SUB(BamKit::VerifyCigarLength)';
	
	my $VCLcigarOperations = &SplitCigar($VCLcigar);
	my $VCLcigar_cal_length=0;
	foreach (@{$VCLcigarOperations}){#calcular cigar length
		unless (defined $_->[0] and defined $_->[1]) {
			print STDERR $VCLsubinfo, "Warnings: invalid cigar: $VCLcigar\n";
			return 0;
		}
		if ($_->[1] =~/^[MZIS=X]{1}$/) {
			$VCLcigar_cal_length+=$_->[0];
		}
	}
	if ($VCLcigar_cal_length == $VCLlength) {
		return 1;
	}
	else {
		print STDERR $VCLsubinfo, "Warnings: cigar: $VCLcigar length != read length $VCLcigar\n";
		return 0;
	}
}



### calculate reference length based cigar
### CalCigarRefLength (CIGAR)
### Global: None
### Dependancy:
### Return the reference length of that alignment in BAM
sub CalCigarRefLength {
	my $CCRLcigar=shift;
	
	my $CCRLsubinfo='SUB(BamKit::CalCigarRefLength)';
	
	my $CCRLcigarOperations = &SplitCigar($CCRLcigar);
	my $CCRLcigar_cal_length=0;
	foreach (@{$CCRLcigarOperations}){#calcular cigar length
		unless (defined $_->[0] and defined $_->[1]) {
			print STDERR $CCRLsubinfo, "Warnings: invalid cigar: $CCRLcigar\n";
			return 0;
		}
		if ($_->[1] =~/^[MD=N]{1}$/) {
			$CCRLcigar_cal_length+=$_->[0];
		}
	}
	return $CCRLcigar_cal_length;
}



### calculate read length based cigar
### CalCigarReadLength (CIGAR)
### Global: None
### Dependancy:
### Return the read length of that BAM alignment
sub CalCigarReadLength {
	my $CCRLcigar=shift;
	
	my $CCRLsubinfo='SUB(BamKit::CalCigarReadLength)';
	
	my $CCRLcigarOperations = &SplitCigar($CCRLcigar);
	my $CCRLcigar_cal_length=0;
	foreach (@{$CCRLcigarOperations}){#calcular cigar length
		unless (defined $_->[0] and defined $_->[1]) {
			print STDERR $CCRLsubinfo, "Warnings: invalid cigar: $CCRLcigar\n";
			return 0;
		}
		if ($_->[1] =~/^[MIS=X]{1}$/) {
			$CCRLcigar_cal_length+=$_->[0];
		}
	}
	return $CCRLcigar_cal_length;
}



### convert bam files into fastq
### Bam2FastQ ($bamin, $fastqout, map_code, MAPQ_code, [path_samtools])
### Global:
### Dependency:
### Note: map_code (0=all, 1=mapped, 2=unmapped)
sub Bam2FastQ {
	my ($BFQbamin, $BFQfqout, $BFQmapcode, $BFQmapq, $BFQpath_samtools)=@_;
	
	local *BAMIN; local *FQOUT;
	my $BFQsubinfo='SUB(BamKit::Bam2FastQ)';
	$BFQpath_samtools='samtools' unless (defined $BFQpath_samtools);
	$BFQmapcode=0 unless (defined $BFQmapcode);
	$BFQmapq=0 unless (defined $BFQmapq);
	
	unless (defined $BFQbamin and -s $BFQbamin) {
		print STDERR "${BFQsubinfo}Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	if (! defined $BFQfqout) {
		print STDERR "${BFQsubinfo}Error: undefined Fastq output\n";
		return $BamKit_failure;
	}
	else {
		unlink $BFQfqout if (-e $BFQfqout);
	}
	
	close BAMIN if (defined fileno(BAMIN));
	unless (open (BAMIN, "$BFQpath_samtools view $BFQbamin | ")) {
		print STDERR "${BFQsubinfo}Error: open BAM: $BFQbamin \n";
		return $BamKit_failure;
	}
	close FQOUT if (defined fileno(FQOUT));
	unless (open (FQOUT, ">$BFQfqout")) {
		print STDERR "${BFQsubinfo}Error: write FQ: $BFQfqout \n";
		return $BamKit_failure;
	}
	my $BFQnumline=0;
	while (my $BFQline1=<BAMIN>) {
		$BFQnumline++;
		chomp $BFQline1;
		my @BFQarr=split(/\t/, $BFQline1);
#Check column number
		if (scalar(@BFQarr)<11) {
			print STDERR "${BFQsubinfo}Warnings: col<11 at line $BFQnumline (Readid: $BFQarr[0]) in BAM $BFQbamin\n";
			next;
		}
#check if mapped
		if ($BFQmapcode==1) {
			next if ($BFQarr[1] & 0x0004);
		}
		elsif ($BFQmapcode==2) {
			next unless ($BFQarr[1] & 0x0004);
		}
##check if MAPQ threshold
		next unless (defined $BFQarr[4] and $BFQarr[4]>=$BFQmapq);
##check if mapped to reverse strand
		my $BFQread_id=$BFQarr[0];
		my $BFQreadseq=$BFQarr[9];
		my $BFQreadqual=$BFQarr[10];
		if ($BFQarr[1] & 0x0010) {###0x0010 = reverse strand
			$BFQreadseq=reverse ($BFQreadseq);
			$BFQreadseq=~tr/ATCGatcg/TAGCtagc/;
			$BFQreadqual=reverse ($BFQreadqual);
		}
##check if one of pair and check R1 or R2
		if ($BFQarr[1] & 0x0001) {
			if ($BFQarr[1] & 0x0040) {
				$BFQread_id='@'.$BFQread_id.'/1';
				print FQOUT "$BFQread_id\n$BFQreadseq\n+\n$BFQreadqual\n";
			}
			elsif ($BFQarr[1] & 0x0080) {
				$BFQread_id='@'.$BFQread_id.'/2';
				print FQOUT "$BFQread_id\n$BFQreadseq\n+\n$BFQreadqual\n";
			}
			else {
				print STDERR "${BFQsubinfo}Warnings: unknown R1 or R2 (FLAG: $BFQarr[1]) at line $BFQnumline (Readid: $BFQarr[0]) in BAM $BFQbamin\n";
				next;
			}
		}
		else {
			$BFQread_id='@'.$BFQread_id;
			print FQOUT "$BFQread_id\n$BFQreadseq\n+\n$BFQreadqual\n";
		}
	}
	close BAMIN;
	close FQOUT;
	return $BamKit_success;
}



### convert BAM to FASTQ using bam2fastq program
### Bam2FastqProg($BAM, $fastq_prefix, $path_bam2fastq)
### Global: $BamKit_success $BamKit_failure
### Dependency: 
### Note:
### Return (1/0, \%hash=('R1' => $fastqR1, 'R2' => $fastqR2, 'M' => $fastqUnpaired)
sub Bam2FastqProg {
	my ($BFPbamfile, $BFPprefix, $BFPpath_bam2fastq)=@_;
	
	my $BFPsubinfo='SUB(Bam2FastqProg)';
	$BFPpath_bam2fastq='bam2fastq' unless (defined $BFPpath_bam2fastq);
	$BFPprefix="MyFastq" unless (defined $BFPprefix and $BFPprefix=~/^\S+$/);
	my %BFQfastqfiles=();
	my $BFPtestout=0;
	
	unless (defined $BFPbamfile and -s $BFPbamfile) {
		print STDERR $BFPsubinfo, "Error: invalid BAM file for fastq: $BFPbamfile\n";
		return $BamKit_failure;
	}
	unless (exec_cmd_return("$BFPpath_bam2fastq -o $BFPprefix.R#.fq --quiet --aligned $BFPbamfile")) {
		print STDERR $BFPsubinfo, "Error: bam2fastq running error\n";
		return $BamKit_failure;
	}
	if (-s "$BFPprefix.R_1.fq") {
		$BFPtestout++;
		$BFQfastqfiles{'R1'}="$BFPprefix.R_1.fq";
	}
	elsif (-e "$BFPprefix.R_1.fq") {
		unlink "$BFPprefix.R_1.fq";
	}
	if (-s "$BFPprefix.R_2.fq") {
		$BFPtestout++;
		$BFQfastqfiles{'R2'}="$BFPprefix.R_2.fq";
	}
	elsif (-e "$BFPprefix.R_2.fq") {
		unlink "$BFPprefix.R_2.fq";
	}
	if (-s "$BFPprefix.R_M.fq") {
		$BFPtestout++;
		$BFQfastqfiles{'M'}="$BFPprefix.R_M.fq";
	}
	elsif (-e "$BFPprefix.R_M.fq") {
		unlink "$BFPprefix.R_M.fq";
	}
	
	unless ($BFPtestout >0) {
		print STDERR $BFPsubinfo, "Error: bam2fastq output error\n";
		return $BamKit_failure;
	}
	
	return ($BamKit_success, \%BFQfastqfiles);
}



### Express RPKM
### ExpressFpkm($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, \@EFsamfiles, \@seq_ids, [path_express])
### Global: $BamKit_debug
### Dependancy: exec_cmd_return, FileKit::RetrieveBasename, FileKit::DeletePath
### Note: bamfiles should be the intact to get the total number of reads
sub ExpressFpkm {
	my ($EFref, $EFfrag_len_mean, $EFfrag_len_stddev, $EFmax_read_len, $EFsamfiles_index, $EFcluster_seqids_arrindex, $EFpath_express)=@_;
	
	local *FPKM;
	my $EFsubinfo='SUB(BamKitExpressFpkm)';
	my @EFfpkms=();
	my $EFi=0;
	$EFpath_express='express' unless (defined $EFpath_express);
#Format: @return_arr=(AABBDD(1/0), AABB(1/0), AA(1/0), DD(1/0))
	my @return_arr=();
	my $EFcurdir=getcwd;
	
	if (scalar(@{$EFcluster_seqids_arrindex}) <1) {
		print STDERR "${EFsubinfo}Error: empty cluster IDs\n";
		return $BamKit_failure;
	}
	if (scalar(@{$EFsamfiles_index}) < 1) {
		print STDERR "${EFsubinfo}Error: empty bam files\n";
		return $BamKit_failure;
	}
	
	unless (-d "$EFcurdir/Express") {
		unless (mkdir ("$EFcurdir/Express", 0766)) {
			print STDERR "${EFsubinfo}Error: can not create folder $EFcurdir/Express\n";
			return $BamKit_failure;
		}
	}
	unlink glob "$EFcurdir/Express/*";

#@EFsamfiles=(bam_AABBDD, bam_AABB, bam_AA, bam_DD) in this order
	foreach my $EFbamfile (@{$EFsamfiles_index}) {
		my $EFbamfile_base=RetrieveBasename($EFbamfile);
		DeletePath("$EFcurdir/Express/$EFbamfile_base") if (-d "$EFcurdir/Express/$EFbamfile_base");
		unless (mkdir ("$EFcurdir/Express/$EFbamfile_base", 0766)) {
			print STDERR "${EFsubinfo}Error: can not create folder $EFcurdir/Express/$EFbamfile_base\n";
			return $BamKit_failure;
		}
#		unlink glob "$EFcurdir/Express/$EFbamfile_base/*.xprs";###delete last-run files
		if (-s $EFbamfile) {
			my $EFcmd="$EFpath_express --frag-len-mean $EFfrag_len_mean --frag-len-stddev $EFfrag_len_stddev --max-read-len $EFmax_read_len --output-dir $EFcurdir/Express/$EFbamfile_base $EFref $EFbamfile 2> /dev/stderr";
			if (! &exec_cmd_return($EFcmd)) {
				return $BamKit_failure;
			}
			elsif (! -s "$EFcurdir/Express/$EFbamfile_base/results.xprs") {
				return $BamKit_failure;
			}
			close FPKM if (defined fileno(FPKM));
			unless (open (FPKM, "$EFcurdir/Express/$EFbamfile_base/results.xprs")) {
				print STDERR "${EFsubinfo}Error: can not open file $EFcurdir/Express/$EFbamfile_base/results.xprs\n";
				return $BamKit_failure;
			}
			my $EFline=<FPKM>;###firstline is header
			my $EFtest_expressed=0;
			while ($EFline=<FPKM>) {
				chomp $EFline;
				my @arr=split(/\t/, $EFline);
#$arr[1] is reference sequence ID
#$arr[10] is the estimated relative abundance in units of fragments per kilobase per million mapped
				if ($arr[10]>0) {###express if any reference FPKM >0 in this cluster
					$EFtest_expressed=1;
				}
				${$EFfpkms[$EFi]}{$arr[1]}=$arr[10];
			}
			close FPKM;
			push (@return_arr, $EFtest_expressed);
			DeletePath("$EFcurdir/Express/$EFbamfile_base") if (-d "$EFcurdir/Express/$EFbamfile_base");#delete last-run files
			#unlink glob "$EFbamfile_base/*.xprs";###delete last-run files
			####Delete directory###PAUSE
			$EFi++;
		}
		else {
			print STDERR "${EFsubinfo}Error: bam file $EFbamfile empty or non-exists\n";
			return $BamKit_failure;
		}
	}
#Check if each cluster_seqids have a FPKM value, 0 if not;
#Output FPKM to $fpkm_log file
	foreach my $EFseqid (@{$EFcluster_seqids_arrindex}) {
		my $EFfpkm_logprint=$EFseqid;
		for (my $i=0; $i<scalar(@EFfpkms); $i++) {
			if (! exists ${$EFfpkms[$i]}{$EFseqid}) {
				${$EFfpkms[$i]}{$EFseqid}=0;
			}
			$EFfpkm_logprint.="\t".${$EFfpkms[$i]}{$EFseqid};
			print "${EFsubinfo}Test: Bam ${$EFsamfiles_index}[$i]: Ref $EFseqid: ${$EFfpkms[$i]}{$EFseqid}\n" if ($BamKit_debug); ###Test###
		}
#		print FPKMLOG $EFfpkm_logprint."\n";
	}
#Return format
#@EFfpkms=(
#			(ref1 => fpkm1, ref2 => fpkm2, ...), 
#			(ref1 => fpkm1, ref2 => fpkm2, ...), 
#			...)
#$EFfpkms[$ith]->{reference_id}=fpkm
	unless (scalar(@{$EFsamfiles_index}) == scalar(@return_arr)) {
		print STDERR "${EFsubinfo}Error: unequal list\n";
		return $BamKit_failure;
	}
	return ($BamKit_success, \@return_arr);
}






### Calculate FPKM
### CalcFPKM(BAMin, $total_mapped_reads, readgroup aware(0/1), [$CFpath_samtools])
### Global: 
### Dependency:
### Note:
sub CalcFPKM {
	my ($CFbamin, $CFtotalmapreads, $CFrg_aware, $CFpath_samtools)=@_;
	
	local *FPKMIN;
	my $CFsubinfo='SUB(BamKit::CalcFPKM)';
	$CFpath_samtools='samtools' unless (defined $CFpath_samtools);
	
	unless (defined $CFbamin and -s $CFbamin) {
		print STDERR "${CFsubinfo}Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $CFrg_aware and ($CFrg_aware==1 or $CFrg_aware==0)) {
		print STDERR "${CFsubinfo}Warnings: ignore read group info when calculate FPKM\n";
		$CFrg_aware=0;
	}
	unless (defined $CFtotalmapreads and $CFtotalmapreads=~/^\d+$/) {
		print STDERR "${CFsubinfo}Error: Unknown total number of reads\n";
		return $BamKit_failure;
	}
	
	close FPKMIN if (defined fileno(FPKMIN));
	unless (open (FPKMIN, "$CFpath_samtools view -h $CFbamin | ")) {
		print STDERR "${CFsubinfo}Error: BAM open error: $CFbamin\n";
		return $BamKit_failure;
	}
	my %CFhashcount=();
	my %CFreflen=();
	while (my $CFline=<FPKMIN>) {
		chomp $CFline;
		if ($CFline=~/^\@/) {
			if ($CFline=~/^\@SQ\s+SN:(\S+)\s+LN:(\d+)$/) {
				$CFreflen{$1}=$2;
			}
			next;
		}
		my $readgroup='unknown';
		
		if ($CFline=~/RG:Z:(\S+)/) {
			$readgroup=$1;
		}
		my @CFarr=split(/\t/, $CFline);
		next unless ($CFarr[1] & 0x0001);##paired
		next if ($CFarr[1] & 0x0004 or $CFarr[1] & 0x0008);##both pair mapped
		if ($CFarr[1] & 0x0040) {##which mate: 1/2
			${${${$CFhashcount{$CFarr[2]}}{$readgroup}}{$CFarr[0]}}{1}++;
		}
		elsif ($CFarr[1] & 0x0080) {
			${${${$CFhashcount{$CFarr[2]}}{$readgroup}}{$CFarr[0]}}{2}++;
		}
		else {
			print STDERR "${CFsubinfo}Warnings: unknown mate 1/2: FLAG: $CFarr[1] of read $CFarr[0] at chr:pos : $CFarr[2]:$CFarr[3] in $CFbamin\n";
		}
	}
	my %CFfpkm=();
	foreach my $CFchrom (sort (keys %CFhashcount)) {
		foreach my $CFindrg (sort (keys %{$CFhashcount{$CFchrom}})) {
			my $CFreadpaircount=0;
			foreach my $CFreadid (keys %{${$CFhashcount{$CFchrom}}{$CFindrg}}) {
				if (exists ${${$CFhashcount{$CFchrom}}{$CFindrg}}{$CFreadid}) {
					if (scalar(keys %{${${$CFhashcount{$CFchrom}}{$CFindrg}}{$CFreadid}}) ==2) {
						$CFreadpaircount++;
					}
				}
			}
			my $CFind_fpkm=$CFreadpaircount * 10^9/($CFreflen{$CFchrom} * $CFtotalmapreads);
			${$CFfpkm{$CFchrom}}{$CFindrg}=$CFind_fpkm;
		}
	}
##Format: $CFfpkm{chromosome} -> {readgroup} = FPKM_value
	return ($BamKit_success, \%CFfpkm);
}


### ReduceReadNameLengthBy read group
### &ReduceReadNameLength($BAMin, $BAMout, $path_samtools)
### Global: $BamKit_failure;
### Dependency: 
### Note: 
### Return: 
sub ReduceReadNameLength {
	my ($RRNLbamin, $RRNLbamout, $RRNLpath_samtools)=@_;
	
	local *RRNLBAMIN; local *RRNLBAMOUT;
	my $RRNLsubinfo='BamKit::ReduceReadNameLength';
	$RRNLpath_samtools='samtools' unless (defined $RRNLpath_samtools);
#	print "${RRNLsubinfo}Test:\n\tBAMin: $RRNLbamin\n\tBAMout: $RRNLbamout\n\tSamtools: $RRNLpath_samtools\n";### For test ###
	
	unless (defined $RRNLbamin and -s $RRNLbamin) {
		print STDERR "${RRNLsubinfo}Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $RRNLbamout and $RRNLbamout=~/^\S+$/) {
		print STDERR "${RRNLsubinfo}Error: invalid BAM output\n";
		return $BamKit_failure;
	}
	unlink $RRNLbamout if (-e $RRNLbamout);
	
	close RRNLBAMIN if (defined fileno(RRNLBAMIN));
	unless (open(RRNLBAMIN, "$RRNLpath_samtools view -h $RRNLbamin | ")) {
		print STDERR "${RRNLsubinfo}Error: can not open BAM: $RRNLbamin\n";
		return $BamKit_failure;
	}
	close RRNLBAMOUT if (defined fileno(RRNLBAMOUT));
	unless (open (RRNLBAMOUT, " | $RRNLpath_samtools view -b -h -S - > $RRNLbamout")) {
		print STDERR "${RRNLsubinfo}Error: can not write BAM: $RRNLbamout\n";
		return $BamKit_failure;
	}
	my $RRNLnum_line=0;
	while (my $RRNLline=<RRNLBAMIN>) {
		$RRNLnum_line++;
		if ($RRNLline=~/^\@/) {
#			print "Header: ".$RRNLline;### For test ###
			print RRNLBAMOUT $RRNLline;
			next;
		}
		my $RRNLrg='';
		if ($RRNLline=~/\tRG:Z:(\S+)/) {
			$RRNLrg=$1;
#			print "ReadGroup: $RRNLrg";### For test ###
		}
		else {
			print STDERR "${RRNLsubinfo}Error: NO RG:Z at line $RRNLnum_line of BAM: $RRNLbamin\n";
			close RRNLBAMIN;
			close RRNLBAMOUT;
			unlink $RRNLbamout if (-s $RRNLbamout);
		}
		if ($RRNLline=~s/^[a-zA-Z0-9-_]+:[a-zA-Z0-9-_]+:[a-zA-Z0-9-_]+/$RRNLrg/) {
#			print "BAMbody: ".$RRNLline;### For test ###
			print RRNLBAMOUT $RRNLline;
		}
		else {
			print STDERR "${RRNLsubinfo}Error: can not replace at line $RRNLnum_line of BAM: $RRNLbamin\n$RRNLline";
			close RRNLBAMIN;
			close RRNLBAMOUT;
			unlink $RRNLbamout if (-s $RRNLbamout);
		}
	}
	close RRNLBAMIN;
	close RRNLBAMOUT;
	return $BamKit_success;
}



### ReadSam
###&ReadSam($file.sam,$ref,fa, ReturnCode(1/2/3))
###Global: $BamKit_success; $BamKit_failure;
###Dependancy: Bio::DB::Sam
###Note: 1=Bio::DB::Sam objective; 2=reference name array, 3=to be defined
sub ReadSam {
	my ($RSbam_file, $RSref_file, $RSret_code)=@_;
	
	my $RSsubinfo='SUB(BamKit::ReadSam)';
	
	unless (defined $RSbam_file and -s $RSbam_file) {
		print STDERR "${RSsubinfo}Error: invalid BAM file\n";
		return $BamKit_failure;
	}
	unless (defined $RSref_file and -s $RSref_file) {
		print STDERR "${RSsubinfo}Error: invalid Fasta file\n";
		return $BamKit_failure;
	}
	unless (defined $RSret_code and $RSret_code=~/^[12]{1}$/) {
		print STDERR "${RSsubinfo}Error: invalid return code\n";
		return $BamKit_failure;
	}
	
	my $RSsam_obj=Bio::DB::Sam->new(  -bam => "$RSbam_file", 
							-fasta => "$RSref_file",
							-expand_flags => 1,
							-split_splices => 1,
							-autoindex => 1,
						);
	if ($RSret_code==1) {
		return ($BamKit_success, $RSsam_obj);
	}
	elsif ($RSret_code==2) {
		my @arr=$RSsam_obj->seq_ids;
		return ($BamKit_success, @arr);
	}
	elsif ($RSret_code==3) {
		##do sth
		return $BamKit_failure;
	}
}



### Filter Bam Reads By Names
### &BamFilterReadsByNames($inputbam, $file_readnames, $filter_code, $outbam[, $path_samtools])
### Global: $BamKit_success; $BamKit_failure;
### Dependancy:
### Note: $file_readnems: per readname perl line
### Note: $code: 1=filter / 0=filter_out
sub BamFilterReadsByNames {
	my ($BFRBNbamin, $BFRBNreadnames, $BFRBNfiltercode, $BFRBNbamout, $BFRBNpath_samtools)=@_;
	
	my $BFRBNsubinfo='SUB(BamKit::BamFilterReadsByNames)';
	$BFRBNpath_samtools='samtools' unless (defined $BFRBNpath_samtools);
	my %BFRBNreadnamehash=();
	my $BFRBNlinenum=0;
	my $BFRBNinclude=0;
	my $BFRBNexclude=0;
	my $BFRBNheader=0;
	local *BFRBNINPUTBAM; local *BFRBNREADNAMEFILE; local *BFRBNOUTPUTBAM;
	
	unless (defined $BFRBNbamin and -s $BFRBNbamin) {
		print STDERR $BFRBNsubinfo, "Error: invalid input BAM/SAM file\n";
		return $BamKit_failure;
	}
	unless ($BFRBNbamin=~/\.[bs]am$/i) {
		print STDERR $BFRBNsubinfo, "Error: input BAM/SAM extension error, should be .sam or .bam [case-insensitive]\n";
		return $BamKit_failure;
	}
	unless (defined $BFRBNreadnames and -s $BFRBNreadnames) {
		print STDERR $BFRBNsubinfo, "Error: invalid input read name file [per readname per line]\n";
		return $BamKit_failure;
	}
	unless (defined $BFRBNfiltercode and $BFRBNfiltercode==0 or $BFRBNfiltercode==1) {
		print STDERR $BFRBNsubinfo, "Error: invalid filter_code, should be 0 or 1\n";
		return $BamKit_failure;
	}
	unless (defined $BFRBNbamout and $BFRBNbamout=~/^\S+$/i) {
		print STDERR $BFRBNsubinfo, "Error: invalid output BAM/SAM file\n";
		return $BamKit_failure;
	}
	unless ($BFRBNbamout=~/\.[bs]am$/i) {
		print STDERR $BFRBNsubinfo, "Error: output BAM/SAM extension error, should be .sam or .bam [case-insensitive]\n";
		return $BamKit_failure;
	}
	
	close BFRBNREADNAMEFILE if (defined fileno(BFRBNREADNAMEFILE));
	unless (open (BFRBNREADNAMEFILE, " < $BFRBNreadnames")) {
		print STDERR $BFRBNsubinfo, "Error: can not open readname file: $BFRBNreadnames\n";
		return $BamKit_failure;
	}
	while (my $BFRBNline=<BFRBNREADNAMEFILE>) {
		chomp $BFRBNline;
		$BFRBNreadnamehash{$BFRBNline}++;
	}
	close BFRBNREADNAMEFILE;
	close BFRBNINPUTBAM if (defined fileno(BFRBNINPUTBAM));
	if ($BFRBNbamin=~/\.bam$/i) {
		unless (open (BFRBNINPUTBAM, "$BFRBNpath_samtools view -h $BFRBNbamin | ")) {
			print STDERR $BFRBNsubinfo, "Error: can not open BAM input: $BFRBNbamin\n";
			return $BamKit_failure;
		}
	}
	elsif ($BFRBNbamin=~/\.sam$/i) {
		unless (open (BFRBNINPUTBAM, "$BFRBNpath_samtools view -S -h $BFRBNbamin | ")) {
			print STDERR $BFRBNsubinfo, "Error: can not open SAM input: $BFRBNbamin\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BFRBNsubinfo, "Error: invalid input BAM/SAM extension: $BFRBNbamin\n";
		return $BamKit_failure;
	}
	close BFRBNOUTPUTBAM if (defined fileno(BFRBNOUTPUTBAM));
	if ($BFRBNbamout=~/\.bam$/i) {
		unless (open (BFRBNOUTPUTBAM, " | $BFRBNpath_samtools view -S -b -h - > $BFRBNbamout")) {
			print STDERR $BFRBNsubinfo, "Error: can not write BAM input: $BFRBNbamout\n";
			return $BamKit_failure;
		}
	}
	elsif ($BFRBNbamout=~/\.sam$/i) {
		unless (open (BFRBNOUTPUTBAM, " | $BFRBNpath_samtools view -S -h - > $BFRBNbamout")) {
			print STDERR $BFRBNsubinfo, "Error: can not write SAM input: $BFRBNbamout\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BFRBNsubinfo, "Error: invalid output BAM/SAM extension: $BFRBNbamout\n";
		return $BamKit_failure;
	}
	while (my $BFRBNline=<BFRBNINPUTBAM>) {
		chomp $BFRBNline;
		$BFRBNlinenum++;
		if ($BFRBNline=~/^\@/) {
			print BFRBNOUTPUTBAM $BFRBNline, "\n";
			$BFRBNheader++;
			next;
		}
		my @BFRBNarr=split(/\t/, $BFRBNline);
		if (defined $BFRBNarr[0] and $BFRBNarr[0]=~/^\S+$/) {
			if (exists $BFRBNreadnamehash{$BFRBNarr[0]}) {
				if ($BFRBNfiltercode==1) {
					print BFRBNOUTPUTBAM $BFRBNline, "\n";
					$BFRBNinclude++;
				}
				else {
					$BFRBNexclude++;
				}
			}
			else {
				if ($BFRBNfiltercode==0) {
					print BFRBNOUTPUTBAM $BFRBNline, "\n";
					$BFRBNinclude++;
				}
				else {
					$BFRBNexclude++;
				}
			}
		}
		else {
			print STDERR $BFRBNsubinfo, "Wanings: invalid line($BFRBNlinenum): $BFRBNline\n";
		}
	}
	close BFRBNINPUTBAM;
	close BFRBNOUTPUTBAM;
	
	print $BFRBNsubinfo, "\n\n\n###SUMMARY: \n";
	print $BFRBNsubinfo, "BAM/SAM input:      $BFRBNbamin\n";
	print $BFRBNsubinfo, "Total lines:        $BFRBNlinenum\n";
	print $BFRBNsubinfo, "    header:         $BFRBNheader\n";
	print $BFRBNsubinfo, "    filtered:       $BFRBNinclude\n";
	print $BFRBNsubinfo, "    filtered out:   $BFRBNexclude\n";

	return $BamKit_success;
}



### Extract reads in a BAM file using specified region in BED format
### &BamExtractReadsUsingBed($inputbam, $file_region, $outreadname[, $path_samtools])
### Global: $BamKit_success; $BamKit_failure;
### Dependancy:
### Note: region file format: 
###		seq_ID\tstart\tend
###		seqID
###		seqID:Num1-Num2
### Note:
sub BamExtractReadsUsingBed {
	my ($BERUBbamin, $BERUBbedfiles, $BERUBreadnameout, $BERUBpath_samtools)=@_;
	
	my $BERUBsubinfo='SUB(BamKit::BamExtractReadsUsingBed)';
	$BERUBpath_samtools='samtools' unless (defined $BERUBpath_samtools);
	my $BERUBlinenum1=0;
	my $BERUBlinenum2=0;
	my %BERUBreadnamehash=();
	my $BERUBoutnum=0;
	local *BERUBREGIONFILE; local *BERUBALIGNMENT; local *BERUBREADNAMEOUT;
	
	unless (defined $BERUBbamin and -s $BERUBbamin) {
		print STDERR $BERUBsubinfo, "Error: invalid input BAM\n";
		return $BamKit_failure;
	}
	unless (-s "$BERUBbamin.bai") {
		print $BERUBsubinfo, "Info: BAM not indexed: $BERUBbamin\n";
		print $BERUBsubinfo, "Info: Now trying to index\n";
		unless (&IndexBam($BERUBbamin, $BERUBpath_samtools)) {
			print STDERR $BERUBsubinfo, "Error: unable to index BAM: $BERUBbamin\n";
			return $BamKit_failure;
		}
	}
	unless (defined $BERUBbedfiles and -s $BERUBbedfiles) {
		print STDERR $BERUBsubinfo, "Error: invalid bed files\n";
		return $BamKit_failure;
	}
	unless (defined $BERUBreadnameout) {
		print STDERR $BERUBsubinfo, "Error: invalid readname output file name\n";
		return $BamKit_failure;
	}
	unlink $BERUBreadnameout if (-e $BERUBreadnameout);
	
	
	print $BERUBsubinfo, "############# SUMMARY ####################\n";
	print $BERUBsubinfo, "###   Input    BAM:    $BERUBbamin\n";
	print $BERUBsubinfo, "###   Input    REGION: $BERUBbedfiles\n";
	print $BERUBsubinfo, "###   Output   READs:  $BERUBreadnameout\n";
	print $BERUBsubinfo, "###   SAMTOOLS path:   $BERUBpath_samtools\n";
	
	close BERUBREGIONFILE if (defined fileno(BERUBREGIONFILE));
	unless (open (BERUBREGIONFILE, "< $BERUBbedfiles")) {
		print STDERR $BERUBsubinfo, "Error: unable to open region file\n";
		return $BamKit_failure;
	}
	close BERUBREADNAMEOUT if (defined fileno(BERUBREADNAMEOUT));
	unless (open (BERUBREADNAMEOUT, " > $BERUBreadnameout ")) {
		print STDERR $BERUBsubinfo, "Error: unable to write region name file: $BERUBreadnameout\n";
		return $BamKit_failure;
	}
	while (my $BERUBline=<BERUBREGIONFILE>) {
		chomp $BERUBline;
		$BERUBlinenum1++;
		my $BERUBregion;
		if ($BERUBline=~/^(\S+)\t(\d+)\t(\d+)/) {### Bed format first number is half open
			$BERUBregion=$1.':'.($2+1).'-'.$3;
		}
		elsif ($BERUBline=~/^(\S+):(\d+)-(\d+)$/) {
			$BERUBregion=$1.':'.$2.'-'.$3;
		}
		elsif ($BERUBline=~/^(\S+)$/) {
			$BERUBregion=$1;
		}
		else {
			print STDERR $BERUBsubinfo, "Error: unable to understand region line: BED $BERUBbedfiles line($BERUBlinenum1) $BERUBline\n";
			return $BamKit_failure;
		}
		close BERUBALIGNMENT if (defined fileno(BERUBALIGNMENT));
#		print $BERUBsubinfo, "Test: $BERUBpath_samtools view $BERUBbamin $BERUBregion\n"; ### For test ###
		unless (open (BERUBALIGNMENT, "$BERUBpath_samtools view $BERUBbamin $BERUBregion | ")) {
			print STDERR $BERUBsubinfo, "Error: unable to open BAM $BERUBbamin REGION $BERUBregion using cmd: $BERUBpath_samtools view $BERUBbamin $BERUBregion\n";
			return $BamKit_failure;
		}
		$BERUBlinenum2=0;
		while (my $BERUBline2=<BERUBALIGNMENT>) {
			chomp $BERUBline2;
			if ($BERUBline2=~/^(\S+)\t/) {
				$BERUBlinenum2++;
				unless (exists $BERUBreadnamehash{$1}) {
					print BERUBREADNAMEOUT $1, "\n";
					$BERUBoutnum++;
				}
				$BERUBreadnamehash{$1}++;
			}
			else {
				print STDERR $BERUBsubinfo, "Warnings: unable to know readname: BED $BERUBbedfiles line($BERUBlinenum1) $BERUBline BAM $BERUBbamin REGION $BERUBregion line($BERUBlinenum2) $BERUBline2\n";
			}
		}
		close BERUBALIGNMENT;
#		if ($BERUBlinenum2==0) {print STDERR $BERUBsubinfo, "Warnings: no alignments: BED $BERUBbedfiles line($BERUBlinenum1) $BERUBline BAM $BERUBbamin REGION $BERUBregion\n";}### For debug
	}
	close BERUBREGIONFILE;
	close BERUBREADNAMEOUT;
	
	print $BERUBsubinfo, "###   Output   READs lines:  $BERUBoutnum\n";
	print $BERUBsubinfo, "###   Output   READs hashes: ", scalar(keys %BERUBreadnamehash), "\n";
	
	return $BamKit_success;
}



### Keep alignments that have both mates mapped
### BamKeepBothMates ($input.R1.bam, $input.R2.bam, $paired.R1.bam, $paired.R2.bam, $unpaired.R1.bam, $unpaired.R2.bam, [path_samtools])
### Global: $BamKit_success; $BamKit_failure;
### Dependency:
### Note:
### Return: 1=Success    0=Failure
sub BamKeepBothMates {
	my ($BKBMbamin1, $BKBMbamin2, $BKBMbamout1, $BKBMbamout2, $BKBMbamunpaird1, $BKBMbamunpaird2, $BKBMpath_samtools)=@_;

#Default
	my $BKBMsubinfo='SUB(BamKit::BamKeepBothMates)';
	$BKBMpath_samtools = 'samtools' unless (defined $BKBMpath_samtools);
	my %BKBMreadhash=();
	my $BKBMnumkeep=0;
	my $BKBMnumtotal=0;
	my $BKBMnumsingle=0;
	my $BKBMtest_unpaired1=0;
	my $BKBMtest_unpaired2=0;
	local *BKBMBAMIN1; local *BKBMBAMIN2;
	local *BKBMBAMOUT1; local *BKBMBAMOUT2;
	local *BKBMUNPAIR1; local *BKBMUNPAIR2;

#Input and Output
	unless (defined $BKBMbamin1 and -s $BKBMbamin1) {
		print STDERR $BKBMsubinfo, "Error: invalid input R1 BAM\n";
		return $BamKit_failure;
	}
	unless (defined $BKBMbamin2 and -s $BKBMbamin2) {
		print STDERR $BKBMsubinfo, "Error: invalid input R2 BAM\n";
		return $BamKit_failure;
	}
	unless (defined $BKBMbamout1 and $BKBMbamout1=~/^\S+$/) {
		print STDERR $BKBMsubinfo, "Error: invalid output paired R1 BAM\n";
		return $BamKit_failure;
	}
	unlink $BKBMbamout1 if (-e $BKBMbamout1);
	unless (defined $BKBMbamout2 and $BKBMbamout2=~/^\S+$/) {
		print STDERR $BKBMsubinfo, "Error: invalid output paired R2 BAM\n";
		return $BamKit_failure;
	}
	unlink $BKBMbamout2 if (-e $BKBMbamout2);
	if (defined $BKBMbamunpaird1) {
		unless ($BKBMbamunpaird1=~/^\S+$/) {
			print STDERR $BKBMsubinfo, "Error: invalid output unpaired R2 BAM\n";
			return $BamKit_failure;
		}
		else {
			$BKBMtest_unpaired1=1;
		}
	}
	if (defined $BKBMbamunpaird2) {
		unless ($BKBMbamunpaird2=~/^\S+$/) {
			print STDERR $BKBMsubinfo, "Error: invalid output unpaired R2 BAM\n";
			return $BamKit_failure;
		}
		else {
			$BKBMtest_unpaired2=1;
		}
	}
	
	print "\n", $BKBMsubinfo, "########### SUMMARY 1 ###############\n";
	print $BKBMsubinfo, " InPut:  input R1 BAM:     $BKBMbamin1\n";
	print $BKBMsubinfo, "         input R2 BAM:     $BKBMbamin2\n";
	print $BKBMsubinfo, " Output: output R1 BAM:    $BKBMbamout1\n";
	print $BKBMsubinfo, "         output R2 BAM:    $BKBMbamout2\n\n";
	
	close BKBMBAMIN1 if (defined fileno(BKBMBAMIN1));
	if ($BKBMbamin1=~/\.bam$/i) {
		unless (open (BKBMBAMIN1, "$BKBMpath_samtools view $BKBMbamin1 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R1 BAM: $BKBMbamin1\n";
			return $BamKit_failure;
		}
	}
	elsif ($BKBMbamin1=~/\.sam$/i) {
		unless (open (BKBMBAMIN1, "$BKBMpath_samtools view -S $BKBMbamin1 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R1 SAM: $BKBMbamin1\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BKBMsubinfo, "Error: can NOT open input R1 SAM/BAM: $BKBMbamin1\n";
		return $BamKit_failure;
	}
	while (my $BKBMline=<BKBMBAMIN1>) {
		chomp $BKBMline;
		next if ($BKBMline=~/^\@/);
		my @BKBMarr=split(/\t/, $BKBMline);
		unless ($BKBMarr[1] & 0x0004) {
			$BKBMreadhash{$BKBMarr[0]}{1}++;
		}
	}
	close BKBMBAMIN1;
	
	close BKBMBAMIN2 if (defined fileno(BKBMBAMIN2));
	if ($BKBMbamin2=~/\.bam$/i) {
		unless (open (BKBMBAMIN2, "$BKBMpath_samtools view $BKBMbamin2 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R2 BAM: $BKBMbamin2\n";
			return $BamKit_failure;
		}
	}
	elsif ($BKBMbamin2=~/\.sam$/i) {
		unless (open (BKBMBAMIN2, "$BKBMpath_samtools view -S $BKBMbamin2 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R2 SAM: $BKBMbamin2\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BKBMsubinfo, "Error: can NOT open input R2 SAM/BAM: $BKBMbamin2\n";
		return $BamKit_failure;
	}
	while (my $BKBMline=<BKBMBAMIN2>) {
		chomp $BKBMline;
		next if ($BKBMline=~/^\@/);
		my @BKBMarr=split(/\t/, $BKBMline);
		unless ($BKBMarr[1] & 0x0004) {
			$BKBMreadhash{$BKBMarr[0]}{2}++;
		}
	}
	close BKBMBAMIN2;

#Mark both mapped mates
	foreach my $BKBMread (keys %BKBMreadhash) {
		$BKBMnumtotal++;
		if (exists $BKBMreadhash{$BKBMread}{1} and $BKBMreadhash{$BKBMread}{1}>0 and exists $BKBMreadhash{$BKBMread}{2} and $BKBMreadhash{$BKBMread}{2}>0) {
			$BKBMreadhash{$BKBMread}{'keep'}++;
			$BKBMnumkeep++;
		}
		else {
			$BKBMnumsingle++;
		}
	}
	print "\n", $BKBMsubinfo, "########### SUMMARY 2 ###############\n";
	print $BKBMsubinfo, " Total number of reads:     $BKBMnumtotal\n";
	print $BKBMsubinfo, " Number to keep:            $BKBMnumkeep\n";
	print $BKBMsubinfo, " Number to exclude:         $BKBMnumsingle\n\n";

#Filter R1
	close BKBMBAMIN1 if (defined fileno(BKBMBAMIN1));
	if ($BKBMbamin1=~/\.bam$/i) {
		unless (open (BKBMBAMIN1, "$BKBMpath_samtools view -h $BKBMbamin1 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R1 BAM: $BKBMbamin1\n";
			return $BamKit_failure;
		}
	}
	elsif ($BKBMbamin1=~/\.sam$/i) {
		unless (open (BKBMBAMIN1, "$BKBMpath_samtools view -S -h $BKBMbamin1 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R1 SAM: $BKBMbamin1\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BKBMsubinfo, "Error: can NOT open input R1 SAM/BAM\n";
		return $BamKit_failure;
	}
	close BKBMBAMOUT1 if (defined fileno(BKBMBAMOUT1));
	if ($BKBMbamout1=~/\.bam$/i) {
		unless (open (BKBMBAMOUT1, " | $BKBMpath_samtools view -S -h -b - > $BKBMbamout1")) {
			print STDERR $BKBMsubinfo, "Error: can NOT write output paired R1 BAM: $BKBMbamout1\n";
			return $BamKit_failure;
		}
	}
	elsif ($BKBMbamout1=~/\.sam$/i) {
		unless (open (BKBMBAMOUT1, " | $BKBMpath_samtools view -S -h - > $BKBMbamout1")) {
			print STDERR $BKBMsubinfo, "Error: can NOT write output paired R1 SAM: $BKBMbamout1\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BKBMsubinfo, "Error: can NOT write output paired R1 SAM/BAM\n";
		return $BamKit_failure;
	}
	if ($BKBMtest_unpaired1==1) {
		close BKBMUNPAIR1 if (defined fileno(BKBMUNPAIR1));
		if ($BKBMbamunpaird1=~/\.bam$/i) {
			unless (open (BKBMUNPAIR1, " | $BKBMpath_samtools view -S -h -b - > $BKBMbamunpaird1")) {
				print STDERR $BKBMsubinfo, "Error: can NOT write output unpaired R1 BAM: $BKBMbamunpaird1\n";
				return $BamKit_failure;
			}
		}
		elsif ($BKBMbamunpaird1=~/\.sam$/i) {
			unless (open (BKBMUNPAIR1, " | $BKBMpath_samtools view -S -h - > $BKBMbamunpaird1")) {
				print STDERR $BKBMsubinfo, "Error: can NOT write output unpaired R1 SAM: $BKBMbamunpaird1\n";
				return $BamKit_failure;
			}
		}
		else {
			print STDERR $BKBMsubinfo, "Error: can NOT write output unpaired R1 SAM/BAM\n";
			return $BamKit_failure;
		}
	}
	while (my $BKBMline=<BKBMBAMIN1>) {
		chomp $BKBMline;
		if ($BKBMline=~/^\@/) {
			print BKBMBAMOUT1 $BKBMline, "\n";
			print BKBMUNPAIR1 $BKBMline, "\n" if ($BKBMtest_unpaired1==1);
			next;
		}
		else{
			my @BKBMarr=split(/\t/, $BKBMline);
			if (exists $BKBMreadhash{$BKBMarr[0]} and exists $BKBMreadhash{$BKBMarr[0]}{'keep'} and $BKBMreadhash{$BKBMarr[0]}{'keep'}>0) {
				print BKBMBAMOUT1 $BKBMline, "\n";
			}
			else {
				print BKBMUNPAIR1 $BKBMline, "\n" if ($BKBMtest_unpaired1==1);
			}
		}
	}
	close BKBMBAMIN1;
	close BKBMBAMOUT1;
	close BKBMUNPAIR1 if ($BKBMtest_unpaired1==1);
	unless (-s $BKBMbamout1) {
		print STDERR $BKBMsubinfo, "Error: output R1 SAM/BAMerror: $BKBMbamout1\n";
		return $BamKit_failure;
	}
#Filter R2
	close BKBMBAMIN2 if (defined fileno(BKBMBAMIN2));
	if ($BKBMbamin2=~/\.bam$/i) {
		unless (open (BKBMBAMIN2, "$BKBMpath_samtools view -h $BKBMbamin2 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R2 BAM: $BKBMbamin2\n";
			return $BamKit_failure;
		}
	}
	elsif ($BKBMbamin2=~/\.sam$/i) {
		unless (open (BKBMBAMIN2, "$BKBMpath_samtools view -S -h $BKBMbamin2 | ")) {
			print STDERR $BKBMsubinfo, "Error: can NOT open input R2 SAM: $BKBMbamin2\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BKBMsubinfo, "Error: can NOT open input R2 SAM/BAM: $BKBMbamin2\n";
		return $BamKit_failure;
	}
	close BKBMBAMOUT2 if (defined fileno(BKBMBAMOUT2));
	if ($BKBMbamout2=~/\.bam$/i) {
		unless (open (BKBMBAMOUT2, " | $BKBMpath_samtools view -S -h -b - > $BKBMbamout2")) {
			print STDERR $BKBMsubinfo, "Error: can NOT write output R2 BAM: $BKBMbamout2\n";
			return $BamKit_failure;
		}
	}
	elsif ($BKBMbamout2=~/\.sam$/i) {
		unless (open (BKBMBAMOUT2, " | $BKBMpath_samtools view -S -h - > $BKBMbamout2")) {
			print STDERR $BKBMsubinfo, "Error: can NOT write output R2 SAM: $BKBMbamout2\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BKBMsubinfo, "Error: can NOT write output R2 SAM/BAM\n";
		return $BamKit_failure;
	}
	
	if ($BKBMtest_unpaired2==1) {
		close BKBMUNPAIR2 if (defined fileno(BKBMUNPAIR2));
		if ($BKBMbamunpaird2=~/\.bam$/i) {
			unless (open (BKBMUNPAIR2, " | $BKBMpath_samtools view -S -h -b - > $BKBMbamunpaird2")) {
				print STDERR $BKBMsubinfo, "Error: can NOT write output unpaired R2 BAM: $BKBMbamunpaird2\n";
				return $BamKit_failure;
			}
		}
		elsif ($BKBMbamunpaird2=~/\.sam$/i) {
			unless (open (BKBMUNPAIR2, " | $BKBMpath_samtools view -S -h - > $BKBMbamunpaird2")) {
				print STDERR $BKBMsubinfo, "Error: can NOT write output unpaired R2 SAM: $BKBMbamunpaird2\n";
				return $BamKit_failure;
			}
		}
		else {
			print STDERR $BKBMsubinfo, "Error: can NOT write output unpaired R2 SAM/BAM\n";
			return $BamKit_failure;
		}
	}
	while (my $BKBMline=<BKBMBAMIN2>) {
		chomp $BKBMline;
		if ($BKBMline=~/^\@/) {
			print BKBMBAMOUT2 $BKBMline, "\n";
			print BKBMUNPAIR2 $BKBMline, "\n" if ($BKBMtest_unpaired2==1);
			next;
		}
		else {
			my @BKBMarr=split(/\t/, $BKBMline);
			if (exists $BKBMreadhash{$BKBMarr[0]} and exists $BKBMreadhash{$BKBMarr[0]}{'keep'} and $BKBMreadhash{$BKBMarr[0]}{'keep'}>0) {
				print BKBMBAMOUT2 $BKBMline, "\n";
			}
			else {
				print BKBMUNPAIR2 $BKBMline, "\n" if ($BKBMtest_unpaired2==1);
			}
		}
	}
	close BKBMBAMIN2;
	close BKBMBAMOUT2;
	close BKBMUNPAIR2 if ($BKBMtest_unpaired2==1);
	unless (-s $BKBMbamout2) {
		print STDERR $BKBMsubinfo, "Error: output R2 SAM/BAMerror: $BKBMbamout2\n";
		return $BamKit_failure;
	}
	
	return $BamKit_success;
}




### Mark separately mapped FLAG to 67 131
### BamMarkPairs ($input.R1.bam, $input.R2.bam, $output.R1.bam, $output.R2.bam, [path_samtools])
### Global: $BamKit_success; $BamKit_failure;
### Dependency:
### Note:
### Return: 1=Success    0=Failure
sub BamMarkPairs {
	my ($BMPbamin1, $BMPbamin2, $BMPbamout1, $BMPbamout2, $BMPpath_samtools)=@_;
#Default
	my $BMPsubinfo='SUB(BamKit::BamMarkPairs)';
	$BMPpath_samtools = 'samtools' unless (defined $BMPpath_samtools);
	local *BMPBAMIN1; local *BMPBAMIN2;
	local *BMPBAMOUT1; local *BMPBAMOUT2; 

#Input and Output
	unless (defined $BMPbamin1 and -s $BMPbamin1) {
		print STDERR $BMPsubinfo, "Error: invalid input R1 BAM\n";
		return $BamKit_failure;
	}
	unless (defined $BMPbamin2 and -s $BMPbamin2) {
		print STDERR $BMPsubinfo, "Error: invalid input R2 BAM\n";
		return $BamKit_failure;
	}
	unless (defined $BMPbamout1 and $BMPbamout1=~/^\S+$/) {
		print STDERR $BMPsubinfo, "Error: invalid output R1 BAM\n";
		return $BamKit_failure;
	}
	unlink $BMPbamout1 if (-e $BMPbamout1);
	unless (defined $BMPbamout2 and $BMPbamout2=~/^\S+$/) {
		print STDERR $BMPsubinfo, "Error: invalid output R2 BAM\n";
		return $BamKit_failure;
	}
	unlink $BMPbamout2 if (-e $BMPbamout2);
	
	print "\n", $BMPsubinfo, "########### SUMMARY 1 ###############\n";
	print $BMPsubinfo, " InPut:  input R1 BAM:     $BMPbamin1\n";
	print $BMPsubinfo, "         input R2 BAM:     $BMPbamin2\n";
	print $BMPsubinfo, " Output: output R1 BAM:    $BMPbamout1\n";
	print $BMPsubinfo, "         output R2 BAM:    $BMPbamout2\n\n";
	
	close BMPBAMIN1 if (defined fileno(BMPBAMIN1));
	if ($BMPbamin1=~/\.bam$/i) {
		unless (open (BMPBAMIN1, "$BMPpath_samtools view -h $BMPbamin1 | ")) {
			print STDERR $BMPsubinfo, "Error: can NOT open input R1 BAM: $BMPbamin1\n";
			return $BamKit_failure;
		}
	}
	elsif ($BMPbamin1=~/\.sam$/i) {
		unless (open (BMPBAMIN1, "$BMPpath_samtools view -S -h $BMPbamin1 | ")) {
			print STDERR $BMPsubinfo, "Error: can NOT open input R1 SAM: $BMPbamin1\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BMPsubinfo, "Error: can NOT open input R1 SAM/BAM: $BMPbamin1\n";
		return $BamKit_failure;
	}
	if ($BMPbamout1=~/\.bam$/i) {
		unless (open (BMPBAMOUT1, " | $BMPpath_samtools view -h -b -S - > $BMPbamout1 ")) {
			print STDERR $BMPsubinfo, "Error: can NOT write output R1 BAM: $BMPbamout1\n";
			return $BamKit_failure;
		}
	}
	elsif ($BMPbamout1=~/\.sam$/i) {
		unless (open (BMPBAMOUT1, " | $BMPpath_samtools view -S -h - > $BMPbamout1 ")) {
			print STDERR $BMPsubinfo, "Error: can NOT write output R1 SAM: $BMPbamout1\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BMPsubinfo, "Error: can NOT write output R1 SAM/BAM: $BMPbamout1\n";
		return $BamKit_failure;
	}
	while (my $BMPline=<BMPBAMIN1>) {
		chomp $BMPline;
		if ($BMPline=~/^\@/) {
			print BMPBAMOUT1 $BMPline, "\n";
		}
		else  {
			my @BMParr=split(/\t/, $BMPline);
			unless ($BMParr[1] & 0x0004) {
				$BMParr[1]+=67;
			}
			print BMPBAMOUT1 join ("\t", @BMParr), "\n";
		}
	}
	close BMPBAMIN1;
	close BMPBAMOUT1;
	unless (-s $BMPbamout1) {
		print STDERR $BMPsubinfo, "Error: output R1 SAM/BAM: $BMPbamout1\n";
		return $BamKit_failure;
	}

	close BMPBAMIN2 if (defined fileno(BMPBAMIN2));
	if ($BMPbamin2=~/\.bam$/i) {
		unless (open (BMPBAMIN2, "$BMPpath_samtools view -h $BMPbamin2 | ")) {
			print STDERR $BMPsubinfo, "Error: can NOT open input R2 BAM: $BMPbamin2\n";
			return $BamKit_failure;
		}
	}
	elsif ($BMPbamin2=~/\.sam$/i) {
		unless (open (BMPBAMIN2, "$BMPpath_samtools view -S -h $BMPbamin2 | ")) {
			print STDERR $BMPsubinfo, "Error: can NOT open input R2 SAM: $BMPbamin2\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BMPsubinfo, "Error: can NOT open input R2 SAM/BAM: $BMPbamin2\n";
		return $BamKit_failure;
	}
	if ($BMPbamout2=~/\.bam$/i) {
		unless (open (BMPBAMOUT2, " | $BMPpath_samtools view -h -b -S - > $BMPbamout2 ")) {
			print STDERR $BMPsubinfo, "Error: can NOT write output R2 BAM: $BMPbamout2\n";
			return $BamKit_failure;
		}
	}
	elsif ($BMPbamout2=~/\.sam$/i) {
		unless (open (BMPBAMOUT2, " | $BMPpath_samtools view -S -h - > $BMPbamout2 ")) {
			print STDERR $BMPsubinfo, "Error: can NOT write output R2 SAM: $BMPbamout2\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BMPsubinfo, "Error: can NOT write output R2 SAM/BAM: $BMPbamout2\n";
		return $BamKit_failure;
	}
	while (my $BMPline=<BMPBAMIN2>) {
		chomp $BMPline;
		if ($BMPline=~/^\@/) {
			print BMPBAMOUT2 $BMPline, "\n";
		}
		else  {
			my @BMParr=split(/\t/, $BMPline);
			unless ($BMParr[1] & 0x0004) {
				$BMParr[1]+=131;
			}
			print BMPBAMOUT2 join ("\t", @BMParr), "\n";
		}
	}
	close BMPBAMIN2;
	close BMPBAMOUT2;
	unless (-s $BMPbamout2) {
		print STDERR $BMPsubinfo, "Error: output R2 SAM/BAM: $BMPbamout2\n";
		return $BamKit_failure;
	}
	
	return $BamKit_success;
}



### remove alignment >k mapping times from BAM
### BamKeepNumAlignemnts(In.name-sorted.bam, $times[1], out.bam)
### Global: $BamKit_success; $BamKit_failure;
### Dependency:
### CMD: samtools
### Note: bam_col1 count >k; so k2 means unique for paired reads
### Return: 1=Success    0= Failure
sub BamKeepNumAlignemnts {
	my ($BAMbamin, $BAMmax, $BAMbamout)=@_;
	
	my $BAMsubinfo="SUB(BamKit::BamKeepNumAlignemnts)";
	my @BAMprint_content=();
	my $BAMheaderline=0;
	my $BAMnumalign_before=0;
	my $BAMnumalign_after=0;
	my $BAMcount=0;
	my $BAMlast_read="";
	unless (defined $BAMmax and $BAMmax=~/^\d+$/ and $BAMmax>0) {
		print $BAMsubinfo, "Warnings: none-INT mapping times, use default 1\n";
		print STDERR $BAMsubinfo, "Warnings: none-INT mapping times, use default 1\n";
		$BAMmax=1;
	}
	local *BAMBAMINPUT; local *BAMBAMOUTPUT;
	
	unless (defined $BAMbamin and -s $BAMbamin) {
		print STDERR $BAMsubinfo, "Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $BAMbamout) {
		print STDERR $BAMsubinfo, "Error: invalid BAM output name\n";
		return $BamKit_failure;
	}
	unlink $BAMbamout if (-e $BAMbamout);
	
	close BAMBAMINPUT if (defined fileno(BAMBAMINPUT));
	unless (open BAMBAMINPUT, "samtools view -h $BAMbamin | ") {
		print STDERR $BAMsubinfo, "Error: can not open BAM input: $BAMbamin\n";
		return $BamKit_failure;
	}
	close BAMBAMOUTPUT if (defined fileno(BAMBAMOUTPUT));
	unless (open BAMBAMOUTPUT, " | samtools view -h -b -S - > $BAMbamout") {
		print STDERR $BAMsubinfo, "Error: can not write BAM output: $BAMbamout\n";
		return $BamKit_failure;
	}
	while (my $BAMline=<BAMBAMINPUT>) {
		chomp $BAMline;
		$BAMnumalign_before++;
		if ($BAMline=~/^\@/) {
			$BAMheaderline++;
			print BAMBAMOUTPUT $BAMline, "\n";
			next;
		}
		my @BAMarr=split(/\t/, $BAMline);
		if ($BAMlast_read eq "" and $BAMcount==0) {
			$BAMlast_read=$BAMarr[0];
			$BAMcount=1;
			push (@BAMprint_content, $BAMline);
			next;
		}
		if ($BAMarr[0] eq $BAMlast_read) {
			push (@BAMprint_content, $BAMline);
			$BAMcount++;
		}
		else {
			if ($BAMcount<=$BAMmax) {
				print BAMBAMOUTPUT join("\n", @BAMprint_content), "\n";
				$BAMnumalign_after+=$BAMcount;
			}
			$BAMlast_read=$BAMarr[0];
			$BAMcount=1;
			@BAMprint_content=();
			push (@BAMprint_content, $BAMline);
		}
	}
	if (scalar(@BAMprint_content)>0 and $BAMcount<=$BAMmax) {
		print BAMBAMOUTPUT join("\n", @BAMprint_content), "\n";
		$BAMnumalign_after+=$BAMcount;
		@BAMprint_content=();
		$BAMcount=0;
	}
	close BAMBAMINPUT;
	close BAMBAMOUTPUT;
	
	print $BAMsubinfo, "######## SUMMARY #########\n";
	print $BAMsubinfo, "Max alignment         : $BAMmax\n";
	print $BAMsubinfo, "Num header lines      : $BAMheaderline\n";
	print $BAMsubinfo, "Num alignments before : $BAMnumalign_before\n";
	print $BAMsubinfo, "Num alignments after  : $BAMnumalign_after\n\n";
	
	return $BamKit_success;
}



### submodules for BamAtacShift
sub trimFirstCIGARunit{
	my $TFCUcigar = shift;
	my $TFCUfirstChar = (split/\d+/, $TFCUcigar)[1];
	#print STDERR "$TFCUcigar\n$TFCUfirstChar\n";
	my $TFCUfirstNum = (split/$TFCUfirstChar/, $TFCUcigar)[0];
	my $TFCUlength =length($TFCUfirstNum.$TFCUfirstChar);
	return(substr($TFCUcigar,$TFCUlength), $TFCUfirstNum, $TFCUfirstChar);
}

# this cuts off the last number and letter from a CIGAR string
# and returns the trimmed CIGAR plus the number and letter trimmed off
sub trimLastCIGARunit{
	my $TLCUcigar = shift;
	my $TLCUlastChar = substr($TLCUcigar,-1); # split was acting weird, so this should grab the last character
	my $TLCUlastNum = (split/[MINDS]/, $TLCUcigar)[-1]; # these are the possible letters equalling Match, Insert, N, Deletion, Shoft clip
	my $TLCUlength =length($TLCUlastNum.$TLCUlastChar);
	return(substr($TLCUcigar,0,(-1*$TLCUlength)), $TLCUlastNum, $TLCUlastChar);
}

# this will take a cigar string, the strand the read maps to and a number of reads to cut off
# it returns the properly trimmed CIGAR string
sub CIGARtrim {
	my ($CTcigar, $CTstrand, $CTnum) = @_;
	
	my $CTtoReturn;
	my $CTsubinfo="SUB(BamKit::CIGARtrim)";
	
	if ($CTstrand eq '-') {
		my ($CTcigarToKeep, $CTlastNum, $CTlastChar) = &trimLastCIGARunit($CTcigar);
		my $CTnewLastNum = $CTlastNum - $CTnum;

		if ($CTlastChar eq 'D'){ # skip to the next letter because the deletion is completely skipped over
			($CTcigarToKeep, $CTlastNum, $CTlastChar) = &trimLastCIGARunit($CTcigarToKeep);
			$CTnewLastNum = $CTlastNum - $CTnum;
		}

		if ($CTnewLastNum==0) {
			my $CTnextChar = substr($CTcigarToKeep,-1);
			if ($CTnextChar eq 'D') { # skip to the next letter because a read shouldn't end in a deletion
				($CTcigarToKeep, $CTlastNum, $CTlastChar) = &trimLastCIGARunit($CTcigarToKeep);
			}
			$CTtoReturn = $CTcigarToKeep;
		}
		elsif ($CTnewLastNum<0) { # Not everything has been trimmed off, so iterate over this again
			$CTtoReturn=&CIGARtrim($CTcigarToKeep,$CTstrand,-1*$CTnewLastNum);
		}
		else {
			$CTtoReturn=$CTcigarToKeep.$CTnewLastNum.$CTlastChar;
		}
	}
	elsif($CTstrand eq '+'){
		# print STDERR "$CTcigar\n";
		my ($CTcigarToKeep, $CTfirstNum, $CTfirstChar) = &trimFirstCIGARunit($CTcigar);
		my $CTnewFirstNum = $CTfirstNum - $CTnum;
		if ($CTfirstChar eq 'D') { # skip to the next letter because the deletion is completely skipped over
			($CTcigarToKeep, $CTfirstNum, $CTfirstChar) = &trimFirstCIGARunit($CTcigarToKeep);
			$CTnewFirstNum = $CTfirstNum - $CTnum;
		}
		if ($CTnewFirstNum==0) {
			my $CTnextChar = substr($CTcigarToKeep,-1);
			if ($CTnextChar eq 'D'){ # skip to the next letter because a read shouldn't end in a deletion
				($CTcigarToKeep, $CTfirstNum, $CTfirstChar) = &trimFirstCIGARunit($CTcigarToKeep);
			}
			$CTtoReturn = $CTcigarToKeep;
		}
		elsif ($CTnewFirstNum<0) { # Not everything has been trimmed off, so iterate over this again
			$CTtoReturn=&CIGARtrim($CTcigarToKeep,$CTstrand,-1*$CTnewFirstNum);
		}
		else {
			$CTtoReturn=$CTnewFirstNum.$CTfirstChar.$CTcigarToKeep;
		}
	}

	return ($CTtoReturn);
}
### Modify BAM file to adjust 9 bp for ATAC-seq: forward strand start +4 and reverse strans start -5
### BamAtacShift($in.bam, $out.bam)
### Global:
### CMD: samtools
### Output: ### Default BAMout=BAMin.ATAC9shift.bam if not specified
### Return: 1=Success    0= Failure
sub BamAtacShift {
	my ($BASbam_in, $BASbam_out)=@_;
	
	my $BASsubinfo="SUB(BamKit::BamAtacShift)";
	my @BASbamarr=();
	my $BASnum_lines_header=0;
	my $BASnum_lines_alignment_before=0;
	my $BASnum_lines_alignment_after=0;
	local *BASBAMINPUT; local *BASBAMOUTPUT;
	
	unless (defined $BASbam_in and -s $BASbam_in) {
		print STDERR $BASsubinfo, "Error: invalid BAM input\n";
		return $BamKit_failure;
	}
	unless (defined $BASbam_out) {
		$BASbam_out="$BASbam_in.ATAC9shift.bam";
		print STDERR $BASsubinfo, "Warnings: BAM output not specified, use default: $BASbam_out\n";
	}
	unlink $BASbam_out if (-e $BASbam_out);
	
	close BASBAMINPUT if (defined fileno(BASBAMINPUT));
	if ($BASbam_in=~/\.bam$/i) {
		unless (open BASBAMINPUT, "samtools view -h $BASbam_in | ") {
			print STDERR $BASsubinfo, "Error: can not open BAM input: $BASbam_in\n";
			return $BamKit_failure;
		}
	}
	elsif ($BASbam_in=~/\.sam$/i) {
		unless (open BASBAMINPUT, "samtools view -S -h $BASbam_in | ") {
			print STDERR $BASsubinfo, "Error: can not open SAM input: $BASbam_in\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BASsubinfo, "Error: can not guess BAM/SAM input format: $BASbam_in\n";
		return $BamKit_failure;
	}
	close BASBAMOUTPUT if (defined fileno(BASBAMOUTPUT)); 
	if ($BASbam_out=~/\.bam$/i) {
		unless (open BASBAMOUTPUT, " | samtools view -S -b -h - > $BASbam_out") {
			print STDERR $BASsubinfo, "Error: can not write BAM output: $BASbam_out\n";
			return $BamKit_failure;
		}
	}
	elsif ($BASbam_out=~/\.sam$/i) {
		unless (open BASBAMOUTPUT, " | samtools view -S -h - > $BASbam_out") {
			print STDERR $BASsubinfo, "Error: can not write SAM output: $BASbam_out\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BASsubinfo, "Error: can not guess BAM/SAM output format: $BASbam_out\n";
		return $BamKit_failure;
	}
	
	while (my $BASline=<BASBAMINPUT>) {
		if ($BASline=~/^\@/) {
			print BASBAMOUTPUT $BASline;
			$BASnum_lines_header++;
			next;
		}
		$BASnum_lines_alignment_before++;
		chomp $BASline;
		@BASbamarr=();
		@BASbamarr=split(/\t/, $BASline);
		unless (defined $BASbamarr[1] and $BASbamarr[1]=~/^\d+$/ and $BASbamarr[1]>0) {
			print STDERR $BASsubinfo, "Warnings: Ignored line as bitwise FLAG: $BASline\n";
			next;
		}
		next unless ($BASbamarr[1] & 0x0002);
		next if ($BASbamarr[1] & 0x0004);### unmapped
		next if ($BASbamarr[1] & 0x0008);### mate unmapped
		next if ($BASbamarr[1] & 0x0100);### not primary
		next if ($BASbamarr[1] & 0x0200);### fail QC
		next if ($BASbamarr[1] & 0x0400);### PCR duplicate
		next if ($BASbamarr[1] & 0x0800);### supplementary alignments
		if($BASbamarr[8]==0){ # just in case
			print STDERR $BASsubinfo, "Warnings: ", $BASbamarr[0]." has a length of 0\n$_\n"; # just in case
		}

		if($BASbamarr[1] & 0x0010){ # if the read is mapped to the reverse strand
			$BASbamarr[7]+=4; # move the mate start 4 bp
			$BASbamarr[8]-=9; # adjust the inferred insert size
			$BASbamarr[9]=substr($BASbamarr[9],0,-5); # remove 5 bp from the reads
			$BASbamarr[10]=substr($BASbamarr[10],0,-5);# remove the 5bp from quality as well
			# this is ugly, but this is a way to change the CIGAR to reflect that I'm shortening the read
			# for the negative stran I want to subtract the 5p from the last entry (which should be M, for match)
			$BASbamarr[5] = &CIGARtrim($BASbamarr[5],'-',5); # from the CIGAR string, trim 5, and we're on the negtive strand
		}
		else { # if the read is mapped to the positive strand
			$BASbamarr[3]+=4; # move the start 4 bp
			# these were removed, because they were incorrect, you don't actually move this end of the negative read (which this would be, because it's the mate of a positive strand read).
			# for a negative read you just trim from the other end, the start position stays the same.
			#$BASbamarr[7]-=5; # move the mate start 5 bp 
			#$BASbamarr[7]=1 if ($BASbamarr[7]<=0); # just to make sure we didn't go off the end
			$BASbamarr[8]-=9; # adjust the inferred insert size
			$BASbamarr[9] =substr($BASbamarr[9],4);# remove 4bp from the reads
			$BASbamarr[10] =substr($BASbamarr[10],4);# remove the 4bp from quality as well
			# this is ugly, but this is a way to change the CIGAR to reflect that I'm shortening the read
			# for the positive strand I want to subtract the 4 from the first entry (which should be M, for match)
			$BASbamarr[5] = &CIGARtrim($BASbamarr[5],'+',4); # from the CIGAR string, trim 4, and we're on the positive strand
		}

		print BASBAMOUTPUT join("\t",@BASbamarr)."\n";
		$BASnum_lines_alignment_after++;
	}

	close BASBAMINPUT;
	close BASBAMOUTPUT;

	print $BASsubinfo, "Info: ############ SUMMARY ##############\n";
	print $BASsubinfo, "Info: BAM input  : $BASbam_in\n";
	print $BASsubinfo, "Info: BAM output : $BASbam_out\n";
	print $BASsubinfo, "Info: Num headers: $BASnum_lines_header\n";
	print $BASsubinfo, "Info: Num Alignments before : $BASnum_lines_alignment_before\n";
	print $BASsubinfo, "Info: Num Alignments after  : $BASnum_lines_alignment_after\n\n";

	return $BamKit_success;
}



### Restore BAM coordintes due to splited reference
### BamRestoreSplit($BRSbamin, $BRSbed, $BRSbamout)
### Global: $BamKit_success; $BamKit_failure;
### Dependency:
### Note:
### Return: 1=Success    0= Failure
### chr1A_part1     0       471304005       chr1A   0       471304005
### chr1A_part2     0       122798051       chr1A   471304005       594102056
### chr1B_part1     0       438720154       chr1B   0       438720154
### chr1B_part2     0       251131716       chr1B   438720154       689851870
### chr1D_part1     0       452179604       chr1D   0       452179604
### chr1D_part2     0       43273582        chr1D   452179604       495453186
sub BamRestoreSplit {
	my ($BRSbamin, $BRSbed, $BRSbamout, $BRSpath_samtools)=@_;
	
	my $BRSsubinfo="SUB(BamKit::BamRestoreSplit)";
	my %BRShash=();
	my %BRSfulllength=();
	my $BRSnumline=0;
	my %BRSlength_print=();
	$BRSpath_samtools='samtools' unless (defined $BRSpath_samtools);
	local *BRSBAMINPUT;
	local *BRSBEDCOORDS;
	local *BRSBAMOUT;
	
	unless (defined $BRSbamin and -s $BRSbamin) {
		print STDERR $BRSsubinfo, "Error: invalid input BAM\n";
		return $BamKit_failure;
	}
	unless (defined $BRSbed and -s $BRSbed) {
		print STDERR $BRSsubinfo, "Error: invalid input BED\n";
		return $BamKit_failure;
	}
	
	close BRSBEDCOORDS if (defined fileno(BRSBEDCOORDS));
	if ($BRSbed=~/\.bed\.gz$/i) {
		unless (open BRSBEDCOORDS, "zcat $BRSbed | ") {
			print STDERR $BRSsubinfo, "Error: can not open gzipped BED coordinates\n";
			return $BamKit_failure;
		}
	}
	elsif ($BRSbed=~/\.bed$/i) {
		unless (open BRSBEDCOORDS, '<', $BRSbed) {
			print STDERR $BRSsubinfo, "Error: can not open flat BED coordinates\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BRSsubinfo, "Error: BED file should have a suffix of .bed or .bed.gz\n";
		return $BamKit_failure;
	}
	while (my $BRSline=<BRSBEDCOORDS>) {
		$BRSnumline++;
		next if ($BRSline=~/^#/);
		chomp $BRSline;
		my @BRSarr=split(/\t/, $BRSline);
		unless (scalar(@BRSarr)>=6) {
			print STDERR $BRSsubinfo, "Warnings: invalid BED line ($BRSnumline): ncol<6\n$BRSline\n";
			next;
		}
		if (exists $BRShash{$BRSarr[0]}) {
			print STDERR $BRSsubinfo, "Error: repeated ID (col1 at line$BRSnumline) in BED file: $BRSarr[0]\n";
			return $BamKit_failure;
		}
		$BRShash{$BRSarr[0]}{'ref'}=$BRSarr[3];
		$BRShash{$BRSarr[0]}{'start'}=$BRSarr[4];
		unless (exists $BRSfulllength{$BRSarr[3]} and $BRSfulllength{$BRSarr[3]}<=$BRSarr[5]) {
			$BRSfulllength{$BRSarr[3]}=$BRSarr[5];### seq_full_length
		}
	}
	close BRSBEDCOORDS;
	
	print $BRSsubinfo, "Info: Read lines: $BRSnumline\n";
	print $BRSsubinfo, "Info: Read parts: ". scalar(keys %BRShash) ."\n";
#	print $BRSsubinfo, "Test: \%BRShash\n"; print Dumper \%BRShash; print "\n";### For test ###
#	print $BRSsubinfo, "Test: \%BRSfulllength\n"; print Dumper \%BRSfulllength; print "\n";### For test ###
	
	close BRSBAMINPUT if (defined fileno(BRSBAMINPUT));
	if ($BRSbamin=~/\.bam$/i) {
		unless (open BRSBAMINPUT, "$BRSpath_samtools view -h $BRSbamin | ") {
			print STDERR $BRSsubinfo, "Error: can not open BAM input\n";
			return $BamKit_failure;
		}
	}
	elsif ($BRSbamin=~/\.sam$/i) {
		unless (open BRSBAMINPUT, "$BRSpath_samtools view -h -S $BRSbamin | ") {
			print STDERR $BRSsubinfo, "Error: can not open SAM input\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BRSsubinfo, "Error: BAM input should have a suffix of .bam or .sam\n";
		return $BamKit_failure;
	}
	close BRSBAMOUT if (defined fileno(BRSBAMOUT));
	if ($BRSbamout=~/\.bam$/i) {
		unless (open BRSBAMOUT, " | $BRSpath_samtools view -h -S -b - > $BRSbamout") {
			print STDERR $BRSsubinfo, "Error: can not write BAM input\n";
			return $BamKit_failure;
		}
	}
	elsif ($BRSbamout=~/\.sam$/i) {
		unless (open BRSBAMOUT, " | $BRSpath_samtools view -h -S - > $BRSbamout") {
			print STDERR $BRSsubinfo, "Error: can not write SAM input\n";
			return $BamKit_failure;
		}
	}
	else {
		print STDERR $BRSsubinfo, "Error: BAM output should have a suffix of .bam or .sam\n";
		return $BamKit_failure;
	}
	while (my $BRSline=<BRSBAMINPUT>) {
		chomp $BRSline;
		if ($BRSline=~/^\@/) {
			if ($BRSline=~/^\@SQ/) {
				my $BRSrefname=$BRSline;
				$BRSrefname=~s/^\@SQ\tSN://;$BRSrefname=~s/\s+LN:.*$//;
				unless (exists $BRShash{$BRSrefname} and exists $BRShash{$BRSrefname}{'ref'}) {
					print STDERR "Warnings: no conversion for seqID /$BRSrefname/\n";
#					print Dumper $BRShash{$BRSrefname};### For test ###
					print BRSBAMOUT $BRSline, "\n";
					next;
				}
				my $BRSnewref=$BRShash{$BRSrefname}{'ref'};
				unless (exists $BRSfulllength{$BRSnewref}) {
					print STDERR "Warnings: conversion for seqID $BRSrefname got no full-length\n";
					return $BamKit_failure;
				}
				unless (exists $BRSlength_print{$BRSnewref}) {
					$BRSline="\@SQ\tSN:".$BRSnewref."\tLN:".$BRSfulllength{$BRSnewref};
					$BRSlength_print{$BRSnewref}++;
				}
			}
		}
		else {
			my @BRSarr2=split(/\t/,$BRSline);
			my $BRSprev_ref1=$BRSarr2[2];
			if (exists $BRShash{$BRSprev_ref1} and defined $BRSarr2[3] and $BRSarr2[3]=~/^\d+$/ and $BRSarr2[3]>0) {
				$BRSarr2[2]=$BRShash{$BRSprev_ref1}{'ref'};
				$BRSarr2[3]+=$BRShash{$BRSprev_ref1}{'start'};
			}
			my $BRSprev_ref2=$BRSarr2[6];
			if (exists $BRShash{$BRSprev_ref2} and defined $BRSarr2[7] and $BRSarr2[7]=~/^\d+$/ and $BRSarr2[7]>0) {
				$BRSarr2[6]=$BRShash{$BRSprev_ref2}{'ref'};
				$BRSarr2[7]+=$BRShash{$BRSprev_ref2}{'start'};
			}
			$BRSline=join("\t", @BRSarr2);
		}
		print BRSBAMOUT $BRSline, "\n";
	}
	close BRSBAMINPUT;
	close BRSBAMOUT;
	
	return $BamKit_success;
}






### 
### 
### Global: $BamKit_success; $BamKit_failure;
### Dependency:
### Note:
### Return: 1=Success    0= Failure
#my $BMPsubinfo='SUB(BamKit::BamMarkPairs)';
1;
__END__
