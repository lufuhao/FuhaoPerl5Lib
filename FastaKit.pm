# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::FastaKit

=head1 SYNOPSIS

Fasta -related tools

=head1 Requirements

    Bio::SeqIO
    FuhaoPerl5Lib::FileKit qw/MoveFile RetrieveDir MergeFiles/;
    FuhaoPerl5Lib::CmdKit;
    FuhaoPerl5Lib::MiscKit qw/IsReference FullDigit/;

=head1 DESCRIPTION

=over 4

=item CdbFasta (fasta_file, [path_cdbfasta])

    * cdbfasta index fasta, return fasta index name for cdbyank
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: (1/0, $indexname)
    * Note: would NOT regenerate if index already exists

=item CdbYank (cdbfasta_index, $out.fa, \@seq_ids, [path_cdbyank])

    * Retrieve fasta sequences using input
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/;
    * Return: 1=Success    0=Failure
    * Note: NOT works very well for massive IDs

=item CdbYankFromFile ($cdbfasta_index, $out.fa, $fasta_id_file, [path_cdbyank])

    * Retrieve fasta sequences using a file containing a list of IDs, 1 ID perl line
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 1=Success    0=Failure

=item CdHitEst(input.fasta, output.fasta, add_cmd, [path_cdhitest])

    * Run cd-hit-est to derep the fasta file
    * Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Note: design your $CHEaddition_cmd='-c 1.00 -n 10 -T 0 -r 1 -d 0 -M 30000';
    * Return: 1=Success    0=Failure

=item Code2AA ($codon)

    * Convert 3-base codon to single-letter AA
    * Return: single-letter AA code

=item CreateFastaRegion (xxx.fa/xxx.fai, bin, xxx.out, [path_samtools])

    * Create fasta region file for freebayes-parallel
    * Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 1=Success    0=Failure

=item ExtractFastaSamtools ($input.fa, $out.fa, $fasta_id_file, [path_samtools])

    * Extract fasta sequences using a file containing a list of IDs, 1 ID perl line
    * For a short list only
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 1=Success    0=Failure

=item ExtractFastaSamtoolsID ($input.fa, $out.fa, "id1 id2", [path_samtools])

    * Extract fasta sequences using s string of IDs
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 1=Success    0=Failure

=item ExtractFastaSamtoolsList ($input.fa, $out.fa, $fasta_id_file, [path_samtools])

    * Extract fasta sequences using a file containing a list of IDs, 1 ID perl line
    * Support long list
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 1=Success    0=Failure

=item ExtractFastaSeqtk (my.fa, output.fa, $id_file, [path_seqtk]);

    * Extract fasta sequences using seqtk with a list of ids/bed, better for large list
    * Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 0=fail; 1=success

=item FastaDedup ($input_fasta, $output_fasta)

    * Remove duplicated fasta
    * Dependancy: Bio::SeqIO;
    * Return: 1=Success    0=Failure

=item Fastq2Fasta ($input.fastq, $out.fa)

    * Convert fastq to fasta
    * Return: 1=Success    0=Failure

=item Frame3Translation ($cdsseq)

    * Translation CDS to protein in 6 frames
    * Return: $prot={ 0=> AA, 1 => AA, 2 => AA}
    * Note: 0,1,2 forsword strand: shift 0/1/2 base from 5 end

=item Frame6Translation ($cdsseq)

    * Translation CDS to protein in 6 frames
    * Return: $prot={ 0=> AA, 1 => AA, 2 => AA, 3 => AA, 4 => AA, 5 => AA}
    * Note: 0,1,2 forsword strand: shift 0/1/2 base from 5 end
    * Note: 3,4,5 reverse complement strand: shift 0/1/2 base from 5 end

=item GuessFormat ($input.seq)

    * Guess sequence file format from extension names

=item IndexFasta (input.fa, [path_samtools])

    * Index fasta using samtools index
    * Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
    * Return: 1=Success    0=failure

=item NumSeq ($input.fa)

    * Count number sequences in fasta
    * Return: 0/number_of_sequence

=item ReadFastaLength ($input.fa[.gz])

     Read fasta seq length into hash \$hash={('seqid1' => length1, 'seqid2' => length2,)}
     Return: (1=Success/0=failure, $hash)

=item RenameFasta ($fastain, $fastqout, $prefix, $num_digit, $RFDesc)

    * Rename fasta using regex
    * Dependency: FuhaoPerl5Lib::MiscKit qw/FullDigit/;
    * Return: 1=Success    0=failure

=item RunCap3(reads[.fq|.fasta], output.fasta, [path_cap3])

    * RunCap3 to assemble fasta/fastq
    * Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/; FuhaoPerl5Lib::FileKit qw/MergeFiles/;
    * Note: fastq input would be convert into fasta first
    * Return: 1=Success    0=failure

=item RunMira4 (fastq, /path/manifest_name, output_fasta, seq_prefix, [Min_alternative_count], [$RMpath_mira4], [num_threads], [$RMtmpdir])

    * Run mira4 to assemble fastq
    * Dependency: Cwd; FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/; FuhaoPerl5Lib::FileKit qw/RetrieveDir/
    * Note: Be sure to chdir back in case of any error
    * Default: 
                [Min_alternative_count]    3
                [$RMpath_mira4]            mira
                [num_threads]              1
                [$RMtmpdir]                ./miratmp
    * Return: 1=Success    0=Failure

=item RunFqTrinity ($in.fq, $out.fa, $RQTadd_cmd, [$path_trinity:Trinity])

    * Run Trinity v2.0.6 for assembly
    * Example: $RQTadd_cmd='--max_memory 2G --run_as_paired --CPU 1 --group_pairs_distance 800 --full_cleanup --min_kmer_cov 3 --min_glue 3'
    * Dependency: Cwd; FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/; FuhaoPerl5Lib::FileKit qw/MoveFile/;
    * Return: 1=Success    0=failure

=item SeqRevComp ($seq_string)

    * Seq Reverse Complement
    * Return: reverse complement seq

=item SplitFastaByLength (input.fasta, total_length, output.prefix)

    * Split Fasta file into several files by total_length
    * Return: (1/0, \@arr(filenames))

=item SplitFastaByNumber (input.fasta, number, output.prefix)

    * Split Fasta file into several files by number of sequence
    * Return: (1/0, \@arr(filenames))



Codon2AA CountFasta
CheckFastaIdDup
RunEmbossStretcher
AnalysisEmbossStretcherOutput

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
package FuhaoPerl5Lib::FastaKit;
use strict;
use warnings;
use Exporter;
use Cwd;
use Data::Dumper qw /Dumper/;
use FuhaoPerl5Lib::FileKit qw/MoveFile RetrieveDir MergeFiles/;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::MiscKit qw/IsReference FullDigit/;
use Bio::SeqIO;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20171109';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(CdbFasta CdbYank CdbYankFromFile ExtractFastaSamtools ExtractFastaSamtoolsID IndexFasta CreateFastaRegion RunMira4 CdHitEst RenameFasta RunFqTrinity SplitFastaByNumber RunCap3 Fastq2Fasta SeqRevComp Codon2AA CountFasta CheckFastaIdDup RunEmbossStretcher AnalysisEmbossStretcherOutput NumSeq FastaDedup ExtractFastaSeqtk Frame3Translation Frame6Translation SplitFastaByLength ReadFastaLength ExtractFastaSamtoolsList);
%EXPORT_TAGS = ( DEFAULT => [qw(CdbFasta CdbYank CdbYankFromFile ExtractFastaSamtools ExtractFastaSamtoolsID IndexFasta CreateFastaRegion RunMira4 RenameFasta RunFqTrinity SplitFastaByNumber RunCap3 Fastq2Fasta SeqRevComp Codon2AA CountFasta CheckFastaIdDup RunEmbossStretcher AnalysisEmbossStretcherOutput NumSeq FastaDedup ExtractFastaSeqtk Frame3Translation Frame6Translation SplitFastaByLength ReadFastaLength ExtractFastaSamtoolsList)],
                 ALL    => [qw(CdbFasta CdbYank CdbYankFromFile ExtractFastaSamtools IndexFasta CreateFastaRegion ExtractFastaSamtoolsID RunMira4 RenameFasta RunFqTrinity SplitFastaByNumber RunCap3 Fastq2Fasta SeqRevComp Codon2AA CountFasta CheckFastaIdDup RunEmbossStretcher AnalysisEmbossStretcherOutput NumSeq FastaDedup ExtractFastaSeqtk Frame3Translation Frame6Translation SplitFastaByLength ReadFastaLength ExtractFastaSamtoolsList)]);

my $FastaKit_success=1;
my $FastaKit_failure=0;
my $FastaKit_debug=0;



my %genetic_code = (
	'UCA' => 'S', # Serine
	'UCC' => 'S', # Serine
	'UCG' => 'S', # Serine
	'UCU' => 'S', # Serine
	'UUC' => 'F', # Phenylalanine
	'UUU' => 'F', # Phenylalanine
	'UUA' => 'L', # Leucine
	'UUG' => 'L', # Leucine
	'UAC' => 'Y', # Tyrosine
	'UAU' => 'Y', # Tyrosine
	'UAA' => '*', # Stop
	'UAG' => '*', # Stop
	'UGC' => 'C', # Cysteine
	'UGU' => 'C', # Cysteine
	'UGA' => '*', # Stop
	'UGG' => 'W', # Tryptophan
	'CUA' => 'L', # Leucine
	'CUC' => 'L', # Leucine
	'CUG' => 'L', # Leucine
	'CUU' => 'L', # Leucine
	'CCA' => 'P', # Proline
	'CAU' => 'H', # Histidine
	'CAA' => 'Q', # Glutamine
	'CAG' => 'Q', # Glutamine
	'CGA' => 'R', # Arginine
	'CGC' => 'R', # Arginine
	'CGG' => 'R', # Arginine
	'CGU' => 'R', # Arginine
	'AUA' => 'I', # Isoleucine
	'AUC' => 'I', # Isoleucine
	'AUU' => 'I', # Isoleucine
	'AUG' => 'M', # Methionine
	'ACA' => 'T', # Threonine
	'ACC' => 'T', # Threonine
	'ACG' => 'T', # Threonine
	'ACU' => 'T', # Threonine
	'AAC' => 'N', # Asparagine
	'AAU' => 'N', # Asparagine
	'AAA' => 'K', # Lysine
	'AAG' => 'K', # Lysine
	'AGC' => 'S', # Serine
	'AGU' => 'S', # Serine
	'AGA' => 'R', # Arginine
	'AGG' => 'R', # Arginine
	'CCC' => 'P', # Proline
	'CCG' => 'P', # Proline
	'CCU' => 'P', # Proline
	'CAC' => 'H', # Histidine
	'GUA' => 'V', # Valine
	'GUC' => 'V', # Valine
	'GUG' => 'V', # Valine
	'GUU' => 'V', # Valine
	'GCA' => 'A', # Alanine
	'GCC' => 'A', # Alanine
	'GCG' => 'A', # Alanine
	'GCU' => 'A', # Alanine
	'GAC' => 'D', # Aspartic Acid
	'GAU' => 'D', # Aspartic Acid
	'GAA' => 'E', # Glutamic Acid
	'GAG' => 'E', # Glutamic Acid
	'GGA' => 'G', # Glycine
	'GGC' => 'G', # Glycine
	'GGG' => 'G', # Glycine
	'GGU' => 'G'  # Glycine
);





### cdbfasta index fasta, return fasta index name for cdbyank
### my ($test_cdbfasta, $cdbfasta_index)=&CdbFasta(fasta_file, [path_cdbfasta])
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Return: (1/0, $indexname)
### Note: would NOT regenerate if index already exists
sub CdbFasta {
	my ($CSfasta_file, $CFpath_cdbfasta)=@_;
	
	my $CFsubinfo='SUB(FastaKit::CdbFasta)';
	$CFpath_cdbfasta='cdbfasta' unless (defined $CFpath_cdbfasta);
	
	unless (defined $CSfasta_file and -s $CSfasta_file) {
		print STDERR "{$CFsubinfo}Error: fasta not exist or empty\n";
		return $FastaKit_failure;
	}
	
	if (-s "$CSfasta_file.cdbz" and -s "$CSfasta_file.cdbz.cidx") {
		return ($FastaKit_success, "$CSfasta_file.cdbz.cidx");
	}
	elsif (-s $CSfasta_file) {
		unless (exec_cmd_return("$CFpath_cdbfasta $CSfasta_file -z $CSfasta_file.cdbz")) {
			print STDERR "${CFsubinfo}Error: cdbfasta running error\n";
			return $FastaKit_failure;
		}
		if (-s "$CSfasta_file.cdbz" and -s "$CSfasta_file.cdbz.cidx") {
			return ($FastaKit_success, "$CSfasta_file.cdbz.cidx");
		}
		else {
			print STDERR "${CFsubinfo}Error: cdbfasta output error\n";
			return $FastaKit_failure;
		}
	}
}





### Retrieve fasta sequences using input
### &CdbYank(cdbfasta_index, $CYoutput, $CYseq_ids_arr_index, [path_cdbyank]);
### Return: 0=fail; 1=success
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: NOT works very well for massive IDs
sub CdbYank {
	my ($CYindex, $CYoutput, $CYseq_ids_index, $CYpath_cdbyank)=@_;
	
	$CYpath_cdbyank='cdbyank' unless (defined $CYpath_cdbyank);
	my $CYsubinfo='SUB(FastaKit::CdbYank)';
	
	unless (defined $CYindex and -s $CYindex) {
		print STDERR "${CYsubinfo}Error: invalid cdbfasta index: $CYindex\n";
		return $FastaKit_failure;
	}
	unlink $CYoutput if (-e $CYoutput);
#Returncode (0=NotRef, 1=ScalarRef, 2=ArrayRef, 3=HashRef, 4=Unknown)
	my $test_id_reference=&IsReference($CYseq_ids_index);
	my $CYseqids_join='';
	if ($test_id_reference>0) {
		if ($test_id_reference==1) {
			$CYseqids_join=${$CYseq_ids_index};
			if ($CYseqids_join=~/\s+/) {
				$CYseqids_join=~s/\s+/\\n/g;
			}
			while ($CYseqids_join=~/\\n\\n/) {
				$CYseqids_join=~s/\\n\\n/\\n/g;
			}
		}
		elsif ($test_id_reference==2) {
			$CYseqids_join=join('\n', @{$CYseq_ids_index});
		}
		else {
			print STDERR "${CYsubinfo}Error: unknown seqids reference type\n";
			return $FastaKit_failure;
		}
	}
	
	unless (defined $CYseqids_join and $CYseqids_join=~/\S+/) {
		print STDERR "${CYsubinfo}Error: invalid cdbfasta IDs: $CYseqids_join\n";
		return $FastaKit_failure;
	}
	
	my $CYcmd="perl -e 'print \"$CYseqids_join\\n\";' | $CYpath_cdbyank $CYindex -o $CYoutput -w";
	unless (exec_cmd_return($CYcmd)) {
		print STDERR "${CYsubinfo}Error: cdbyank running error\n";
		return $FastaKit_failure;
	}
	unless (-s $CYoutput) {
		print STDERR "${CYsubinfo}Error: cdbyank output error\n";
		return $FastaKit_failure;
	}

	return $FastaKit_success;
}





### Retrieve fasta sequences using a file containing a list of IDs, 1 ID perl line
### CdbYankFromFile(cdbfasta_index, $CYoutput, $CYfasta_id_file, [path_cdbyank]);
### Return: 0=fail; 1=success
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: 
sub CdbYankFromFile {
	my ($CYFFindex, $CYFFoutput, $CYFFseqidsfile, $CYFFpath_cdbyank)=@_;
	
	local *CYFFID;
	$CYFFpath_cdbyank='cdbyank' unless (defined $CYFFpath_cdbyank);
	my $CYFFsubinfo='SUB(FastaKit::CdbYankFromFile)';
	
	unless (defined $CYFFindex and -s $CYFFindex) {
		print STDERR $CYFFsubinfo, "Error: invalid cdbfasta index: $CYFFindex\n";
		return $FastaKit_failure;
	}
	unlink $CYFFoutput if (-e $CYFFoutput);
	unless (-s $CYFFseqidsfile) {
		print STDERR $CYFFsubinfo, "Error: Invalid Fasta ID list file: $CYFFseqidsfile\n";
		return $FastaKit_failure;
	}
	
	my $CYFFcmd="perl -ne 'print;' $CYFFseqidsfile | $CYFFpath_cdbyank $CYFFindex -o $CYFFoutput -w";
	unless (exec_cmd_return($CYFFcmd)) {
		print STDERR $CYFFsubinfo, "Error: cdbyank running error\n";
		return $FastaKit_failure;
	}
	unless (-s $CYFFoutput) {
		print STDERR $CYFFsubinfo, "Error: cdbyank output error\n";
		return $FastaKit_failure;
	}

	return $FastaKit_success;
}





### Extract fasta sequences using a file containing a list of IDs, 1 ID perl line
### &ExtractFastaSamtools(my.fa, output.fa, fasta_id_file, [path_samtools]);
### Return: 0=fail; 1=success
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: 
sub ExtractFastaSamtools {
	my ($EFSinputfa, $EFSoutputfa, $EFSidlist, $EFSpath_samtools)=@_;
	
	my $EFSsubinfo='SUB(FastaKit::ExtractFastaSamtools)';
	$EFSpath_samtools='samtools' unless (defined $EFSpath_samtools);
	
	unless (defined $EFSinputfa and -s $EFSinputfa) {
		print STDERR $EFSsubinfo, "Error: invalid input fasta file: $EFSinputfa\n";
		return $FastaKit_failure;
	}
	unless (-s "$EFSinputfa.fai") {
		unless (&IndexFasta($EFSinputfa, $EFSpath_samtools)) {
			print STDERR $EFSsubinfo, "Error: unable to index fasta file: $EFSinputfa\n";
			return $FastaKit_failure;
		}
	}
	unless (defined $EFSidlist and -s $EFSidlist) {
		print STDERR $EFSsubinfo, "Error: invalid fasta ID list file: $EFSidlist\n";
		return $FastaKit_failure;
	}
	unless (defined $EFSoutputfa and $EFSoutputfa=~/^\S+$/) {
		print STDERR $EFSsubinfo, "Error: invalid output fasta file: $EFSoutputfa\n";
		return $FastaKit_failure;
	}
	unlink $EFSoutputfa if (-e $EFSoutputfa);
	
	unless (exec_cmd_return("$EFSpath_samtools faidx $EFSinputfa `cat $EFSidlist` > $EFSoutputfa")) {
		print STDERR $EFSsubinfo, "Error: samtools extract fasta running error\n";
		return $FastaKit_failure;
	}
	unless (-s $EFSoutputfa) {
		print STDERR $EFSsubinfo, "Error: samtools extract fasta output error: $EFSoutputfa\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}





### Extract fasta sequences using a file containing a list of IDs, 1 ID perl line
### &ExtractFastaSamtools(my.fa, output.fa, fasta_id_file, [path_samtools]);
### Return: 0=fail; 1=success
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: 
sub ExtractFastaSamtoolsList {
	my ($EFSLinputfa, $EFSLoutputfa, $EFSLidlist, $EFSLpath_samtools)=@_;
	
	my $EFSLsubinfo='SUB(FastaKit::ExtractFastaSamtoolsList)';
	$EFSLpath_samtools='samtools' unless (defined $EFSLpath_samtools);
	local *EFSLLIST; local *EFSLOUTPUT;
	
	unless (defined $EFSLinputfa and -s $EFSLinputfa) {
		print STDERR $EFSLsubinfo, "Error: invalid input fasta file: $EFSLinputfa\n";
		return $FastaKit_failure;
	}
	unless (-s "$EFSLinputfa.fai") {
		unless (&IndexFasta($EFSLinputfa, $EFSLpath_samtools)) {
			print STDERR $EFSLsubinfo, "Error: unable to index fasta file: $EFSLinputfa\n";
			return $FastaKit_failure;
		}
	}
	unless (defined $EFSLidlist and -s $EFSLidlist) {
		print STDERR $EFSLsubinfo, "Error: invalid fasta ID list file: $EFSLidlist\n";
		return $FastaKit_failure;
	}
	unless (defined $EFSLoutputfa and $EFSLoutputfa=~/^\S+$/) {
		print STDERR $EFSLsubinfo, "Error: invalid output fasta file: $EFSLoutputfa\n";
		return $FastaKit_failure;
	}
	unlink $EFSLoutputfa if (-e $EFSLoutputfa);
	
	close EFSLLIST if (defined fileno(EFSLLIST));
	unless (open (EFSLLIST, " < $EFSLidlist ")) {
		print STDERR $EFSLsubinfo, "Error: open list file error\n";
		return $FastaKit_failure;
	}
#	close EFSLOUTPUT if(defined fileno(EFSLOUTPUT));
#	unless (open (EFSLOUTPUT, " > $EFSLoutputfa")) {
#		print STDERR $EFSLsubinfo, "Error: can not write output fasta\n";
#		return $FastaKit_failure;
#	}
	while (my $EFSLline=<EFSLLIST>) {
		chomp $EFSLline;
		unless (exec_cmd_return("$EFSLpath_samtools faidx $EFSLinputfa $EFSLline >> $EFSLoutputfa")) {
			print STDERR $EFSLsubinfo, "Error: samtools extract fasta running error: $EFSLline\n";
			return $FastaKit_failure;
		}
	}
	
#	close EFSLLIST;
	unless (-s $EFSLoutputfa) {
		print STDERR $EFSLsubinfo, "Error: samtools extract fasta output error: $EFSLoutputfa\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}





### Extract fasta sequences using s string of IDs
### &ExtractFastaSamtoolsID(my.fa, output.fa, "id1 id2", [path_samtools]);
### Return: 0=fail; 1=success
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: 
sub ExtractFastaSamtoolsID {
	my ($EFSinputfa, $EFSoutputfa, $EFSidlist, $EFSpath_samtools)=@_;
	
	my $EFSsubinfo='SUB(FastaKit::ExtractFastaSamtoolsID)';
	$EFSpath_samtools='samtools' unless (defined $EFSpath_samtools);
	
	unless (defined $EFSinputfa and -s $EFSinputfa) {
		print STDERR $EFSsubinfo, "Error: invalid input fasta file: $EFSinputfa\n";
		return $FastaKit_failure;
	}
	unless (-s "$EFSinputfa.fai") {
		unless (&IndexFasta($EFSinputfa, $EFSpath_samtools)) {
			print STDERR $EFSsubinfo, "Error: unable to index fasta file: $EFSinputfa\n";
			return $FastaKit_failure;
		}
	}
	unless (defined $EFSidlist and $EFSidlist=~/\S+/) {
		print STDERR $EFSsubinfo, "Error: invalid fasta ID list file: $EFSidlist\n";
		return $FastaKit_failure;
	}
	unless (defined $EFSoutputfa and $EFSoutputfa=~/^\S+$/) {
		print STDERR $EFSsubinfo, "Error: invalid output fasta file: $EFSoutputfa\n";
		return $FastaKit_failure;
	}
	unlink $EFSoutputfa if (-e $EFSoutputfa);
	
	unless (&exec_cmd_return("$EFSpath_samtools faidx $EFSinputfa $EFSidlist > $EFSoutputfa")) {
		print STDERR $EFSsubinfo, "Error: samtools extract fasta running error\n";
		return $FastaKit_failure;
	}
	unless (-s $EFSoutputfa) {
		print STDERR $EFSsubinfo, "Error: samtools extract fasta output error: $EFSoutputfa\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}





### Extract fasta sequences using seqtk with a list of ids/bed, better for large list
### &ExtractFastaSeqtk (my.fa, output.fa, $id_file, [path_seqtk]);
### Return: 0=fail; 1=success
### Global:
### Dependancy: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: 
sub ExtractFastaSeqtk {
	my ($EFSinputfa, $EFSoutputfa, $EFSidfile, $EFSpath_seqtk)=@_;
	
	my $EFSsubinfo='SUB(FastaKit::ExtractFastaSeqtk)';
	$EFSpath_seqtk='seqtk' unless (defined $EFSpath_seqtk);
	
	unless (defined $EFSinputfa and -s $EFSinputfa) {
		print STDERR $EFSsubinfo, "Error: invalid input fasta file: $EFSinputfa\n";
		return $FastaKit_failure;
	}
	unless (-s "$EFSinputfa.fai") {
		unless (&IndexFasta($EFSinputfa, $EFSpath_seqtk)) {
			print STDERR $EFSsubinfo, "Error: unable to index fasta file: $EFSinputfa\n";
			return $FastaKit_failure;
		}
	}
	unless (defined $EFSidfile and $EFSidfile=~/^\S+$/) {
		print STDERR $EFSsubinfo, "Error: invalid fasta ID list file: $EFSidfile\n";
		return $FastaKit_failure;
	}
	unless (defined $EFSoutputfa and $EFSoutputfa=~/^\S+$/) {
		print STDERR $EFSsubinfo, "Error: invalid output fasta file: $EFSoutputfa\n";
		return $FastaKit_failure;
	}
	unlink $EFSoutputfa if (-e $EFSoutputfa);
	
	unless (&exec_cmd_return("$EFSpath_seqtk subseq $EFSinputfa $EFSidfile > $EFSoutputfa")) {
		print STDERR $EFSsubinfo, "Error: seqtk extract fasta running error\n";
		return $FastaKit_failure;
	}
	unless (-s $EFSoutputfa) {
		print STDERR $EFSsubinfo, "Error: seqtk extract fasta output error: $EFSoutputfa\n";
		return $FastaKit_failure;
	}
	
	return $FastaKit_success;
}




### Index fasta using samtools index
### IndexFasta (input.fa, [path_samtools])
### Global:
### Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Return: 1=Success    0=failure
sub IndexFasta {
	my ($IFfasta, $IFpath_samtools)=@_;
	
	my $IFsubinfo='SUB(FastaKit::IndexFasta)';
	$IFpath_samtools='samtools' unless (defined $IFpath_samtools);
	
	unless (defined $IFfasta and -s $IFfasta) {
		print STDERR "${IFsubinfo}Error: invalid fasta file\n";
		return $FastaKit_failure;
	}
	unlink "$IFfasta.fai" if (-e "$IFfasta.fai");
	
	unless (exec_cmd_return("$IFpath_samtools faidx $IFfasta")) {
		print STDERR "${IFsubinfo}Error: samtools index running error\n";
		return $FastaKit_failure;
	}
	unless (-s "$IFfasta.fai") {
		print STDERR "${IFsubinfo}Error: samtools index output error\n";
		return $FastaKit_failure;
	}
	
	return $FastaKit_success;
}





### Count number sequences in fasta
### NumSeq($input.fa)
### Global:
### Dependency:
### Note: 
### Return: number_of_sequence
sub NumSeq {
	my $NSseqfile=shift;

	my $NSsubinfo='SUB(FastaKit::NumSeq)';
	my $NSnumberseq=0;
	local *NSSEQFILE;

	unless (-s $NSseqfile) {
		print STDERR $NSsubinfo, "Error: sequence file not existed: $NSseqfile\n";
		return 0;
	}

	close NSSEQFILE if (defined fileno(NSSEQFILE));
	unless (open (NSSEQFILE, "< $NSseqfile")) {
		print STDERR $NSsubinfo, "Error: can not open sequence file: $NSseqfile\n";
		return 0;
	}
	while (my $NSline=<NSSEQFILE>) {
		chomp $NSline;
		if ($NSline=~/^>\S+/) {
			$NSnumberseq++;
		}
	}
	close NSSEQFILE;
	return $NSnumberseq;
}





### Guess sequence file format from extension names
### GuessFormat ($input.seq)
sub GuessFormat {
	my $GFseqfile = shift @_;
#	print $GFseqfile."\n";     #For TEST
	$GFseqfile =~ s/^.*\.(\w+)$/$1/;
#	print $GFseqfile."\n";   #For TEST
	my $GFfailed = 0;
	my $GFseq_format;
 	SW: {
	if ($GFseqfile =~ /(^fasta)|(^fast)|(^fst)|(^fsa)|(^ft)|(^fs)|(^fa)|(^fas)/i) {$GFseq_format = 'Fasta'; last SW};
	if ($GFseqfile =~ /(^fastq)|(^fq)/i) {$GFseq_format = 'Fastq'; last SW};
	if ($GFseqfile =~ /(lfasta)|(lfast)|(lfst)|(lfsa)|(lft)|(lfs)/i) {$GFseq_format = 'LabeledFasta'; last SW};
	if ($GFseqfile =~ /(embl)|(emb)|(em)|(eml)/i) {$GFseq_format = 'EMBL'; last SW};
	if ($GFseqfile =~ /(genebank)|(genbank)|(genb)|(geneb)|(gbank)|(gb)/i) {$GFseq_format = 'GenBank'; last SW};
	if ($GFseqfile =~ /(swissprot)|(sprt)|(swissp)|(sprot)|(sp)|(spr)/i) {$GFseq_format = 'Swissprot'; last SW};
	if ($GFseqfile =~ /pir/i) {$GFseq_format = 'PIR'; last SW};
	if ($GFseqfile =~ /gcg/i) {$GFseq_format = 'GCG'; last SW};
	if ($GFseqfile =~ /scf/i) {$GFseq_format = 'SCF'; last SW};
	if ($GFseqfile =~ /ace/i) {$GFseq_format = 'Ace'; last SW};
	if ($GFseqfile =~ /phd/i) {$GFseq_format = 'phd'; last SW};
	if ($GFseqfile =~ /phred|phd/i) {$GFseq_format = 'phred'; last SW};
	if ($GFseqfile =~ /raw/i) {$GFseq_format = 'raw'; last SW};
	if ($GFseqfile =~ /bsml/i) {$GFseq_format = 'bsml'; last SW};
	$GFfailed++;
	}
	return eval{$GFfailed ? 0 : $GFseq_format};
}





### Create fasta region file for freebayes-parallel
### CreateFastaRegion (xxx.fa/xxx.fai, bin, xxx.out, [path_samtools])
### Global:
### Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note:
### Return: 1=Success    0=Failure
sub CreateFastaRegion {
	my ($CFRinput, $CFRbinzise, $CFRoutput, $CFpath_samtools)=@_;
	
	local *REGIONIN; local *REGIONOUT;
	my $CFsubinfo='SUB(FastaKit::CreateFastaRegion)';
	$CFpath_samtools='samtools' unless (defined $CFpath_samtools);
	
	unless (defined $CFRinput and -s $CFRinput) {
		print STDERR $CFsubinfo, "Error: invalid fasta index\n";
		return $FastaKit_failure;
	}
	if ($CFRinput !~ /\.fai$/i) {
		if ($CFRinput =~ /\.(fa)|(fasta)|(fas)$/i) {
			print $CFsubinfo, "Info: input format fasta\n" if ($FastaKit_debug); ### For test ###
			unless (-s "$CFRinput.fai") {
				print $CFsubinfo, "Info: $CFRinput.fai not exists; use samtools faidx to generates ..." if ($FastaKit_debug);### For test ###
				unless (&IndexFasta($CFRinput, $CFpath_samtools)) {
					print STDERR $CFsubinfo, "Error: samtools faidx $CFRinput running error\n";
					return $FastaKit_failure;
				}
			}
			$CFRinput.='.fai';
			unless (-s $CFRinput) {
				print STDERR $CFsubinfo, "Error: samtools faidx $CFRinput output error\n";
				return $FastaKit_failure;
			}
		}
		else {
			print STDERR $CFsubinfo, "Error: unknown input format for samtools faidx error\n";
			return $FastaKit_failure;
		}
	}
	
	unless (defined $CFRoutput and $CFRoutput !~ /^\s*$/) {
		print STDERR $CFsubinfo, "Error: invalid region file output";
		return $FastaKit_failure;
	}
	unlink $CFRoutput if (defined $CFRoutput and -e $CFRoutput);
	
	close REGIONIN if (defined fileno(REGIONIN));
	unless (open(REGIONIN, "<$CFRinput")) {
		print STDERR $CFsubinfo, "Error: open $CFRinput error\n";
		return $FastaKit_failure;
	}
	close REGIONOUT if (defined fileno(REGIONOUT));
	unless (open(REGIONOUT, ">$CFRoutput")) {
		close REGIONIN;
		print STDERR $CFsubinfo, "Error: write $CFRoutput error\n";
		return $FastaKit_failure;
	}
	while (my $CFRline=<REGIONIN>) {
		chomp $CFRline;
		my @CFRarr=();
		@CFRarr=split(/\t/, $CFRline);
#		if (scalar(@CFRarr) != 5) {###Force to check fai column number
#			print STDERR $CFsubinfo, "Error: samtools faidx $CFRinput running error\n";
#			print STDERR $CFsubinfo, "Error: fai colnum != 5\n";
#			return $FastaKit_failure;
#		}
		if ($CFRbinzise==0) {
			print REGIONOUT $CFRarr[0]."\n";
		}
		elsif ($CFRbinzise>0) {
			my $CFRstart=0;
			my $CFRend=0;
			while ($CFRstart < $CFRarr[1]) {
				$CFRend=($CFRstart+$CFRbinzise > $CFRarr[1]) ? $CFRarr[1] : $CFRstart+$CFRbinzise;
				print REGIONOUT $CFRarr[0].":$CFRstart-$CFRend\n";
				$CFRstart=$CFRend;
			}
		}
		else {
			print STDERR $CFsubinfo, "Error: Wrong bin size\n";
			return $FastaKit_failure;
		}
	}
	close REGIONIN;
	close REGIONOUT;
	unless (-s $CFRoutput) {
		print STDERR $CFsubinfo, "Error: output $CFRoutput not exists\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}





### Run mira4 to assemble fastq
### RunMira4(fastq, /path/manifest_name, output_fasta, seq_prefix, [Min_alternative_count], [$RMpath_mira4], [num_threads], [$RMtmpdir])
### Global: $FastaKit_success $FastaKit_failure
### Dependency: Cwd; FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/; FuhaoPerl5Lib::FileKit qw/RetrieveDir/;
### Note: Be sure to chdir back in case of any error
### Return: 1=Success    0=Failure
### Default: 
###                [Min_alternative_count]    3
###                [$RMpath_mira4]            mira
###                [num_threads]              1
###                [$RMtmpdir]                ./miratmp
sub RunMira4 {
	my ($RMfastq, $RMmira_manifest, $RMassembly_fasta, $RMseq_prefix, $RMmin_alternative_count, $RMpath_mira4, $RMnum_threads, $RMtmpdir)=@_;

	local *MANIFEST; local *RMOUT; local *RMFA;
	my $RMsubinfo='SUB(FastaKit::RunMira)';
	my $RMreturn_fasta;
	$RMmin_alternative_count=3 unless (defined $RMmin_alternative_count);
	$RMpath_mira4='mira' unless (defined $RMpath_mira4);
	$RMnum_threads=1 unless (defined $RMnum_threads);
	my $RMcurpath=getcwd;
	$RMseq_prefix='MyHapAssem_' unless (defined $RMseq_prefix);
	$RMtmpdir="$RMcurpath/miratmp" unless (defined $RMtmpdir and $RMtmpdir ne '');

	unless (defined $RMmira_manifest and $RMmira_manifest=~/^\S+$/) {
		print STDERR $RMsubinfo, "Error: invalid MIRA4 manifest file\n";
		return $FastaKit_failure;
	}
	my $RMrundir=RetrieveDir($RMmira_manifest);
	unless (-d $RMrundir) {
		print STDERR $RMsubinfo, "Error: can not retrieve manifest path for $RMmira_manifest\n";
		return $FastaKit_failure;
	}
#	print "\n", $RMsubinfo, "Test: \n###Fastq: $RMfastq\n###Manifest: $RMmira_manifest\n###Outfasta: $RMassembly_fasta\n###SeqPrefix: $RMseq_prefix\n###MinCount: $RMmin_alternative_count\n###MIRA4_path: $RMpath_mira4\n###Threads: $RMnum_threads\n\n";### For test ###
	
	unless (-s $RMfastq) {
		print STDERR $RMsubinfo, "Error: can not find fastq input\n";
		return $FastaKit_failure;
	}
	
	close MANIFEST if (defined fileno(MANIFEST));
	unless (open (MANIFEST, ">$RMmira_manifest")) {
		print STDERR $RMsubinfo, "Error: can not open manifest file: $RMmira_manifest\n";
		return $FastaKit_failure;
	}
	my $RMproject='Ta';
	my $RMjob='denovo,est,accurate';
	my $RMparameters="COMMON_SETTINGS -GENERAL:number_of_threads=$RMnum_threads -NAG_AND_WARN:cnfs=warn:check_template_problems=no:check_maxreadnamelength=no -CO:mark_repeats=yes:assume_snp_instead_repeat=yes:name_prefix=$RMseq_prefix -OUT:output_result_caf=no:output_result_tcs=no:output_result_maf=no ";
	$RMparameters.=" -DI:tmp_redirected_to=$RMtmpdir " if ($RMtmpdir ne '');
	$RMparameters.=" SOLEXA_SETTINGS -CO:min_reads_per_group=$RMmin_alternative_count -AS:minimum_reads_per_contig=3";
	print MANIFEST "project = $RMproject\njob = $RMjob\nparameters = $RMparameters\n###Readgroup\n";
	print MANIFEST "readgroup = wheat\nautopairing\ndata = $RMfastq\ntechnology = solexa\ntemplate_size = 50 1000 autorefine\nsegment_placement = ---> <---\n";
	close MANIFEST;
	unless (chdir $RMrundir) {
		print STDERR $RMsubinfo, "Error: can not chdir to manifest path: $RMrundir\n";
		return $FastaKit_failure;
	}
	unless (exec_cmd_return("$RMpath_mira4 $RMmira_manifest >mira.log 2>mira.err")) {
		print STDERR $RMsubinfo, "Error: MIRA return non-zero code\n";
		chdir $RMcurpath;
		return $FastaKit_failure;
	}
#	my @RMfasta_files=glob ("$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}*out.unpadded.fasta");#, "$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}_LargeContigs_out.unpadded.fasta");
	my @RMfasta_files=glob ("$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}*out.padded.fasta");
	if (scalar(@RMfasta_files)<1) {
		print STDERR $RMsubinfo, "Error: MIRA output empty\n";
		return $FastaKit_failure;
	}
	elsif (scalar(@RMfasta_files)==1) {
		$RMreturn_fasta=shift @RMfasta_files;
	}
	else {
		$RMreturn_fasta=$RMrundir."/merge.fa";
		close RMOUT unless (defined fileno(RMOUT));
		unless (open (RMOUT, ">$RMreturn_fasta")) {
			print STDERR $RMsubinfo, "Error: merge fasta1: open $RMreturn_fasta\n";
			chdir $RMcurpath;
			return $FastaKit_failure;
		}
		foreach (@RMfasta_files) {
			unless (open (RMFA, " < $_ ")) {
				print STDERR "${RMsubinfo}Error: merge fasta2: open $_\n";
				chdir $RMcurpath;
				return $FastaKit_failure;
			}
			while (<RMFA>) {
				print RMOUT $_;
			}
			close RMFA;
		}
		close RMOUT;
	}
	unless (rename ($RMreturn_fasta, $RMassembly_fasta) and -s $RMassembly_fasta) {
		print STDERR $RMsubinfo, "Error: move fasta error: $RMreturn_fasta to $RMassembly_fasta\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}





### RunCap3 to assemble fasta
### RunCap3(reads[.fq|.fasta], output.fasta, [path_cap3])
### Global: $FastaKit_failure, $FastaKit_success
### Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/; FuhaoPerl5Lib::FileKit qw/MergeFiles/;
### Note: fastq input would be convert into fasta first
sub RunCap3 {
	my ($RCreads, $RCout, $RCpath_cap3)=@_;
	
	my $RCsubinfo='SUB(FastaKit::RunCap3)';
	$RCpath_cap3='cap3' unless (defined $RCpath_cap3);
	my @RCtempsuf=('.cap.ace', '.cap.contigs', '.cap.contigs.links', '.cap.contigs.qual', '.cap.info', '.cap.singlets');
	
	unless (defined $RCreads and -s $RCreads) {
		print STDERR "Error: invalid CAP3 input\n";
		return $FastaKit_failure;
	}
	if ($RCreads=~/(\.fasta$)|(\.fa$)|(\.fas$)/i) {
		print $RCsubinfo, "Info: Cap3 input in fasta format: $RCreads\n"; ### For test ###
	}
	elsif ($RCreads=~/(\.fastq$)|(\.fq$)/i) {
		print $RCsubinfo, "Info: Cap3 input in fastq format: $RCreads\n"; ### For test ###
		unless (&Fastq2Fasta($RCreads, "$RCreads.fq2fa.fa")) {
			print STDERR $RCsubinfo, "Error: fastq2fasta failed: $RCreads\n";
			return $FastaKit_failure;
		}
		$RCreads="$RCreads.fq2fa.fa";
	}
	else {
		print STDERR $RCsubinfo, "Error: unable to guess CAP3 input format: $RCreads\n";
		return $FastaKit_failure;
	}
	unless (defined $RCout and $RCout=~/^\S+$/) {
		print STDERR $RCsubinfo, "Error: invalid output fasta name: $RCout\n";
		return $FastaKit_failure;
	}
	unlink $RCout if (-e $RCout);
	
	unless (exec_cmd_return("$RCpath_cap3 $RCreads > cap3.log 2>cap3.err")) {
		print STDERR $RCsubinfo, "Error: cap3 running error\n";
		return $FastaKit_failure;
	}
	my @RCoutputs = map{$RCreads.$_}(".cap.contigs", ".cap.singlets");
#	print $RCsubinfo, "Test: @RCoutputs\n"; ### For test ###
	unless (MergeFiles($RCout, @RCoutputs)) {
		print STDERR $RCsubinfo, "Error: merge CAP3 out error\n";
		return $FastaKit_failure;
	}
#	unlink glob "$RCreads.*";
	foreach (@RCtempsuf) {
		my $RCtempfile=$RCreads.$_;
		unlink $RCtempfile if (-e $RCtempfile);
	}
	unless (defined $RCout and -s $RCout) {
		print STDERR $RCsubinfo, "Error: RunCap3 output error\n";
		return $FastaKit_failure;
	}
	
	return $FastaKit_success;
}





### Convert fastq to fasta
### Fastq2Fasta ($input.fastq, $out.fa)
### Global: $FastaKit_failure;$FastaKit_success;
### Dependency: 
### Note: 
sub Fastq2Fasta {
	my ($FFfq, $FFfa)=@_;
	
	local *FFFQIN; local *FFFAOUT;
	my $FFsubinfo='SUB(FastaKit::Fastq2Fasta)';
	
	unless (defined $FFfq and -s $FFfq) {
		print STDERR $FFsubinfo, "Error: invalid input FASTQ file: $FFfq\n";
		return $FastaKit_failure;
	}
	unless (defined $FFfa and $FFfa=~/^\S+$/) {
		print STDERR $FFsubinfo, "Error: invalid output FASTA name: $FFfa\n";
		return $FastaKit_failure;
	}
	unlink $FFfa if (-e $FFfa);

	close FFFQIN if (defined fileno(FFFQIN));
	unless (open (FFFQIN, "< $FFfq")) {
		print STDERR $FFsubinfo, "Error: unable to open FASTQ: $FFfq\n";
		return $FastaKit_failure;
	}
	close FFFAOUT if (defined fileno(FFFAOUT));
	unless (open (FFFAOUT, "> $FFfa")) {
		print STDERR $FFsubinfo, "Error: unable to wirte FASTA: $FFfa\n";
		return $FastaKit_failure;
	}
	while (my $FFline=<FFFQIN>) {
		chomp $FFline;
		$FFline=~s/^\@/>/;
		print FFFAOUT $FFline."\n";;
		$FFline=<FFFQIN>;
		print FFFAOUT $FFline;
		$FFline=<FFFQIN>;
		$FFline=<FFFQIN>;
	}
	close FFFQIN;
	close FFFAOUT;

	unless (defined $FFfa and -s $FFfa) {
		print STDERR $FFsubinfo, "Error: output NOT exists: $FFfa\n";
		return $FastaKit_failure;
	}

	return $FastaKit_success;
}





### run cd-hit-est to derep the fasta file MIRA4 assembled
### &CdHitEst(input.fasta, output.fasta, add_cmd, [path_cdhitest])
### Global: $FastaKit_failure;$FastaKit_success;
### Dependency: FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/
### Note: $CHEaddition_cmd='-c 1.00 -n 10 -T 0 -r 1 -d 0 -M 30000';
### Return: 1=Success    0=Failure
sub CdHitEst {
	my ($CHEfastain, $CHEfastaout, $CHEaddition_cmd, $CHEpath_cdhitest)=@_;
	
	my $CHEsubinfo='SUB(FastaKit::CdHitEst)';
	$CHEpath_cdhitest='cd-hit-est' unless (defined $CHEpath_cdhitest);
	$CHEaddition_cmd='';

	unless (defined $CHEfastain and -s $CHEfastain) {
		print STDERR $CHEsubinfo, "Error: fasta input for CDHIT not defined or exists\n";
		return $FastaKit_failure;
	}
	unless (defined $CHEfastaout and $CHEfastaout=~/^\S+$/) {
		print STDERR $CHEsubinfo, "Error: fasta output for CDHIT not defined\n";
		return $FastaKit_failure;
	}
	if (-e $CHEfastaout) {
		print STDERR $CHEsubinfo, "Warnings: fasta output for CDHIT existed but deleted\n";
		unlink $CHEfastaout;
	}

	unless (exec_cmd_return("$CHEpath_cdhitest -i $CHEfastain -o $CHEfastaout $CHEaddition_cmd")) {
		print STDERR $CHEsubinfo, "Error: CDHIT running error: $CHEfastain\n";
		return $FastaKit_failure;
	}
	unless (defined $CHEfastaout and -s $CHEfastaout) {
		print STDERR $CHEsubinfo, "Error: CDHIT output error: $CHEfastain\n";
		return $FastaKit_failure;
	}

	return $FastaKit_success;
}



### Run Trinity v2.0.6 for assembly
### RunFqTrinity($FastqIn, $fastaout, $RQTadd_cmd, $path_trinity)
### $RQTadd_cmd='--max_memory 2G --run_as_paired --CPU 1 --group_pairs_distance 800 --full_cleanup --min_kmer_cov 3 --min_glue 3'
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: Cwd; FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/; FuhaoPerl5Lib::FileKit qw/MoveFile/;
### Note:
### Return: 1=Success    0=failure
sub RunFqTrinity {
	my ($RQTfastq, $RQTfastaout, $RQTadd_cmd, $path_trinity)=@_;
	
	my $RQTsubinfo='SUB(FastaKit::RunTrinity)';
	$path_trinity='Trinity' unless (defined $path_trinity);
	$RQTadd_cmd=' ' unless (defined $RQTadd_cmd);
	my $RQTseqtype='fq';
	my $RQToutput='Trinity';
	my $RQTcurdir=getcwd;
	
	unless (defined $RQTfastq and -s $RQTfastq) {
		print STDERR $RQTsubinfo, "Error: invalid FastQ/Fasta input\n";
		return $FastaKit_failure;
	}
	
	if ($RQTfastq=~/(\.fq$)|(\.fastq$)/i) {
		$RQTseqtype='fq';
	}
	elsif ($RQTfastq=~/(\.fa$)|(\.fasta$)|(\.fas$)/i) {
		$RQTseqtype='fa';
	}
	else {
		print STDERR $RQTsubinfo, "Error: guess FastQ/Fasta format error\n";
		return $FastaKit_failure;
	}
	
	unless (exec_cmd_return("$path_trinity --seqType $RQTseqtype --single $RQTfastq $RQTadd_cmd --output $RQToutput > trinity.log 2> trinity.err")) {
		print STDERR $RQTsubinfo, "Error: Trinity running error\n";
		return $FastaKit_failure;
	}
	my @RQTfasta_files=glob "$RQTcurdir/$RQToutput.*.f*a";
	unless (scalar(@RQTfasta_files)==1) {
		print STDERR $RQTsubinfo, "Error: Trinity output detecting error: @RQTfasta_files\n";
		return $FastaKit_failure;
	}
	unlink $RQTfastaout if (-e $RQTfastaout);
	unless (MoveFile($RQTfasta_files[0], $RQTfastaout)) {
		print STDERR $RQTsubinfo, "Error: Trinity output rename error: $RQTfasta_files[0] to $RQTfastaout\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}





### Rename fasta using regex
### RenameFasta($fastain, $fastqout, $prefix, $num_digit, $RFDesc)
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: FuhaoPerl5Lib::MiscKit qw/FullDigit/;
### Note:
### Return: 1=Success    0=failure
sub RenameFasta {
	my ($RFfastain, $RFfaout, $RFprefix, $RFnumdigit, $RFseq_desc)=@_;
	
	local *RFFASTAIN; local *RFFASTAOUT;
	my $RFsubinfo='SUB(FastaKit::RenameFasta)';
	my $RFstartnumber=0;
	my $RFnew_seqid='';
	my $RFnew_desc='';
	$RFseq_desc='' unless (defined $RFseq_desc);
	my %RFcheckid_duplicated=();
	
	unless (defined $RFfastain and -s $RFfastain) {
		print STDERR $RFsubinfo, "Error: invalid Fasta input\n";
		return $FastaKit_failure;
	}
	unless (defined $RFfaout and $RFfaout=~/^\S+$/) {
		print STDERR $RFsubinfo, "Error: invalid Fasta output\n";
		return $FastaKit_failure;
	}
	unlink $RFfaout if (-e $RFfaout);

	close RFFASTAIN if (defined fileno(RFFASTAIN));
	unless (open (RFFASTAIN, " < $RFfastain")) {
		print STDERR $RFsubinfo, "Error: can not open fasta: $RFfastain\n";
		return $FastaKit_failure;
	}
	close RFFASTAOUT if (defined fileno(RFFASTAOUT));
	unless (open(RFFASTAOUT, " > $RFfaout")) {
		print STDERR $RFsubinfo, "Error: can not write fasta: $RFfaout\n";
		return $FastaKit_failure;
	}
	my $RFnum_line=0;
	while (my $RFline=<RFFASTAIN>) {
		$RFnum_line++;
		if ($RFline=~/^>/){
			
			$RFstartnumber++;
			chomp $RFline;
			$RFline=~s/^>//;
			$RFnew_desc=$RFseq_desc;
			$RFnew_desc.=' '.$RFline;
			if ($RFline=~/^(\S+)\s*.*$/) {
				if (defined $RFprefix) {
					$RFnew_seqid=$RFprefix;
					if (defined $RFnumdigit) {
						$RFnew_seqid.=FullDigit($RFstartnumber, $RFnumdigit);
					}
					else {
						$RFnew_seqid.=$RFstartnumber;
					}
				}
				else {
					if (defined $RFnumdigit) {
						$RFnew_seqid.=FullDigit($RFstartnumber, $RFnumdigit);
					}
					else {
						$RFnew_seqid.=$RFstartnumber;
					}
				}
			}
			else {
				print STDERR $RFsubinfo, "Error: invalid fastaid: $RFline\n";
				close RFFASTAOUT;
				unlink $RFfaout if (-e $RFfaout);
				close RFFASTAIN;
				return $FastaKit_failure;
			}
			if (exists $RFcheckid_duplicated{$RFnew_seqid}) {
				print STDERR $RFsubinfo, "Error: seqid duplicated: $RFnew_seqid\n";
				close RFFASTAOUT;
				unlink $RFfaout if (-e $RFfaout);
				close RFFASTAIN;
				return $FastaKit_failure;
			}
			else {
				$RFcheckid_duplicated{$RFnew_seqid}++;
			}
			print RFFASTAOUT '>'.$RFnew_seqid.' '.$RFnew_desc."\n";
		}
		else {
			print RFFASTAOUT $RFline;
			next;
		}
	}
	
	close RFFASTAOUT;
	close RFFASTAIN;
	return $FastaKit_success;
}





### split Fasta file into several files by number of sequence
### SplitFastaByNumber (input.fasta[.gz], number, output.prefix)
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: 
### Note: 
### Return: (1/0, \@arr(filenames))
sub SplitFastaByNumber {
	my ($SFBNfasta, $SFBNnum, $SFBNoutprefix)=@_;
	
	local *SFBNFASTAIN; local *SFBNFASTAOUT;
	my $SFBNsubinfo='SUB(FastaKit::SplitFastaByNumber)';
	$SFBNnum=100 unless (defined $SFBNnum and $SFBNnum=~/^\d+$/ and $SFBNnum>0);
	unless (defined $SFBNoutprefix) {
		$SFBNoutprefix=$SFBNfasta;
		$SFBNoutprefix=~s/^.*\///;
	}

	unless (defined $SFBNfasta and -s $SFBNfasta) {
		print STDERR $SFBNsubinfo, "Error: invalid fasta input\n";
		return $FastaKit_failure;
	}
	my $SFBNcount=0;
	my $SFBNcountfileno='00000001';
	my @SFBNfilelist=();

	close SFBNFASTAIN if (defined fileno(SFBNFASTAIN));
	if ($SFBNfasta=~/\.gz$/i) {
		unless (open (SFBNFASTAIN, " zcat $SFBNfasta | ")) {
			print STDERR "${SFBNsubinfo}Error: can not open gzipped fasta: $SFBNfasta\n";
			return $FastaKit_failure;
		}
	}
	else {
		unless (open (SFBNFASTAIN, " < $SFBNfasta")) {
			print STDERR "${SFBNsubinfo}Error: can not open fasta: $SFBNfasta\n";
			return $FastaKit_failure;
		}
	}
	close SFBNFASTAOUT if (defined fileno(SFBNFASTAOUT));
	unless (open (SFBNFASTAOUT, " > $SFBNoutprefix.$SFBNcountfileno.fa")) {
		print STDERR $SFBNsubinfo, "Error: can not write fasta: $SFBNoutprefix.$SFBNcountfileno.fa\n";
		return $FastaKit_failure;
	}
	while (my $SFBNline=<SFBNFASTAIN>) {
		if ($SFBNline=~/^>/) {
			$SFBNcount++;
			if ($SFBNcount>$SFBNnum) {
				close SFBNFASTAOUT if (defined fileno(SFBNFASTAOUT));
				if (-s "$SFBNoutprefix.$SFBNcountfileno.fa") {
					push (@SFBNfilelist, "$SFBNoutprefix.$SFBNcountfileno.fa");
				}
				$SFBNcountfileno++;
				unlink ("$SFBNoutprefix.$SFBNcountfileno.fa") if (-e "$SFBNoutprefix.$SFBNcountfileno.fa");
				unless (open (SFBNFASTAOUT, ">$SFBNoutprefix.$SFBNcountfileno.fa")) {
					print STDERR $SFBNsubinfo, "Error: can not write fasta: $SFBNoutprefix.$SFBNcountfileno.fa\n";
					return $FastaKit_failure;
				}
				$SFBNcount=1;
			}
			print SFBNFASTAOUT $SFBNline;
		}
		else {
			print SFBNFASTAOUT $SFBNline;
		}
	}
	close SFBNFASTAIN;
	close SFBNFASTAOUT;
	if (-s "$SFBNoutprefix.$SFBNcountfileno.fa") {
		push (@SFBNfilelist, "$SFBNoutprefix.$SFBNcountfileno.fa");
	}
	return ($FastaKit_success, \@SFBNfilelist);
}





### split Fasta file into several files by total length
### SplitFastaByNumber (input.fasta, length, output.prefix)
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: 
### Note: 
### Return: (1/0, \@arr(filenames))
sub SplitFastaByLength {
	my ($SFBLfasta, $SFBLnum, $SFBLoutprefix)=@_;
	
	local *SFBLFASTAIN; local *SFBLFASTAOUT;
	my $SFBLsubinfo='SUB(FastaKit::SplitFastaByLength)';
	unless (defined $SFBLnum and $SFBLnum=~/^\d+$/ and $SFBLnum>0) {
		print STDERR $SFBLsubinfo, "Error: invalid total length\n";
		return $FastaKit_failure;
	}
	unless (defined $SFBLoutprefix) {
		$SFBLoutprefix=$SFBLfasta;
		$SFBLoutprefix=~s/^.*\///;
	}
	unless (defined $SFBLfasta and -s $SFBLfasta) {
		print STDERR $SFBLsubinfo, "Error: invalid fasta input\n";
		return $FastaKit_failure;
	}
	my $SFBLcountfileno='00000001';
	my @SFBLfilelist=();
	my $SFBLindfilelength=0;
	my $SFBLindfilenum=0;
	
	my ($SFBLtest, $SFBLseqlenhash)=&ReadFastaLength($SFBLfasta);
	unless ($SFBLtest) {
		print STDERR $SFBLsubinfo, "Error: ReadFastaLength failed: $SFBLfasta\n";
		return $FastaKit_failure;
	}
	
	foreach my $SFBLindseq (sort keys %{$SFBLseqlenhash}) {
		if (${$SFBLseqlenhash}{$SFBLindseq} !~/^\d+/) {
			print STDERR $SFBLsubinfo, "Error: invalid seq : $SFBLindseq\n";
			return $FastaKit_failure;
		}
		elsif (${$SFBLseqlenhash}{$SFBLindseq}<=0) {
			print STDERR $SFBLsubinfo, "Warnings: empty seq : $SFBLindseq\n";
		}
		elsif (${$SFBLseqlenhash}{$SFBLindseq}>$SFBLnum) {
			print STDERR $SFBLsubinfo, "Error: larger seq : $SFBLindseq", ${$SFBLseqlenhash}{$SFBLindseq}, "larger than required total length $SFBLnum\n";
			return $FastaKit_failure;
		}
	}
	
	if ($SFBLfasta=~/\.gz$/i) {
		unless (open (SFBLFASTAIN, " zcat $SFBLfasta | ")) {
			print STDERR $SFBLsubinfo, "Error: can not open gzipped fasta: $SFBLfasta\n";
			return $FastaKit_failure;
		}
	}
	else {
		unless (open (SFBLFASTAIN, " < $SFBLfasta")) {
			print STDERR $SFBLsubinfo, "Error: can not open fasta: $SFBLfasta\n";
			return $FastaKit_failure;
		}
	}
	close SFBLFASTAOUT if (defined fileno(SFBLFASTAOUT));
	unless (open (SFBLFASTAOUT, " > $SFBLoutprefix.$SFBLcountfileno.fa")) {
		print STDERR $SFBLsubinfo, "Error: can not write fasta: $SFBLoutprefix.$SFBLcountfileno.fa\n";
		return $FastaKit_failure;
	}
	while (my $SFBLline=<SFBLFASTAIN>) {
		if ($SFBLline=~/^>/) {
			my $SFBLid=$SFBLline; 
			if ($SFBLline=~/^>(\S+)\s*/) {
				$SFBLid=$1;
			}
			else {
				print STDERR $SFBLsubinfo, "Error: invalid seq header (a line $.): $SFBLline\n";
				return $FastaKit_failure;
			}
			unless (exists ${$SFBLseqlenhash}{$SFBLid}) {
				print STDERR $SFBLsubinfo, "Error: seqid no length: $SFBLid\n";
				return $FastaKit_failure;
			}
			if (($SFBLindfilelength + ${$SFBLseqlenhash}{$SFBLid})>$SFBLnum) {
				close SFBLFASTAOUT if (defined fileno(SFBLFASTAOUT));
				if (-s "$SFBLoutprefix.$SFBLcountfileno.fa") {
					push (@SFBLfilelist, "$SFBLoutprefix.$SFBLcountfileno.fa");
					print $SFBLsubinfo, "Info: num $SFBLindfilenum total $SFBLindfilelength bp $SFBLoutprefix.$SFBLcountfileno.fa\n";
				}
				
				$SFBLcountfileno++;
				unlink ("$SFBLoutprefix.$SFBLcountfileno.fa") if (-e "$SFBLoutprefix.$SFBLcountfileno.fa");
				unless (open (SFBLFASTAOUT, ">$SFBLoutprefix.$SFBLcountfileno.fa")) {
					print STDERR $SFBLsubinfo, "Error: can not write fasta: $SFBLoutprefix.$SFBLcountfileno.fa\n";
					return $FastaKit_failure;
				}
				$SFBLindfilenum=0;
				$SFBLindfilelength=0;
			}
			print SFBLFASTAOUT $SFBLline;
			$SFBLindfilelength = $SFBLindfilelength + ${$SFBLseqlenhash}{$SFBLid};
			$SFBLindfilenum++;
		}
		else {
			print SFBLFASTAOUT $SFBLline;
		}
		
	}
	close SFBLFASTAIN;
	close SFBLFASTAOUT;
	if (-s "$SFBLoutprefix.$SFBLcountfileno.fa") {
		push (@SFBLfilelist, "$SFBLoutprefix.$SFBLcountfileno.fa");
		print $SFBLsubinfo, "Info: num $SFBLindfilenum total $SFBLindfilelength bp $SFBLoutprefix.$SFBLcountfileno.fa\n";
	}
	print $SFBLsubinfo, "Info: total splited files: ", scalar(@SFBLfilelist), "\n";
	return ($FastaKit_success, \@SFBLfilelist);
}


	my $SFBLindfilelength=0;
	my $SFBLindfilenum=0;


### Seq Reverse Complement
### SeqRevComp($seq_string)
### Global: 
### Dependency: 
### Note:
### Return: reverse complement seq
sub SeqRevComp {
	my $SRCoriseq=shift;
	my $SRCseq_rev=reverse $SRCoriseq;
	$SRCseq_rev=~tr/ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv/TACGAtacgaNnRrYySsWwMmKkVvHhDdBb/;
	return $SRCseq_rev;
}





### Codon2AA
### &Code2AA('$codon')
### Global:
### Dependency: 
### Note:
sub Codon2AA {
	my $CAcodon=shift;
	
	$CAcodon=uc $CAcodon;
	$CAcodon=~tr/T/U/;
	
	if (exists $genetic_code{$CAcodon}) {
		return $genetic_code{$CAcodon};
	}
	else {
		return '?';
	}
}
### The Redundancy of the Genetic Code
sub codon2aa2 {
    my($CAcodon) = @_;
	
	$CAcodon=UC $CAcodon;
	$CAcodon=~tr/T/U/;
	
	if ( $CAcodon =~ /GC./i)        { return 'A' }    # Alanine
	elsif ( $CAcodon =~ /TG[TC]/i)     { return 'C' }    # Cysteine
    elsif ( $CAcodon =~ /GA[TC]/i)     { return 'D' }    # Aspartic Acid
    elsif ( $CAcodon =~ /GA[AG]/i)     { return 'E' }    # Glutamic Acid
    elsif ( $CAcodon =~ /TT[TC]/i)     { return 'F' }    # Phenylalanine
    elsif ( $CAcodon =~ /GG./i)        { return 'G' }    # Glycine
    elsif ( $CAcodon =~ /CA[TC]/i)     { return 'H' }    # Histidine
    elsif ( $CAcodon =~ /AT[TCA]/i)    { return 'I' }    # Isoleucine
    elsif ( $CAcodon =~ /AA[AG]/i)     { return 'K' }    # Lysine
    elsif ( $CAcodon =~ /TT[AG]|CT./i) { return 'L' }    # Leucine
    elsif ( $CAcodon =~ /ATG/i)        { return 'M' }    # Methionine
    elsif ( $CAcodon =~ /AA[TC]/i)     { return 'N' }    # Asparagine
    elsif ( $CAcodon =~ /CC./i)        { return 'P' }    # Proline
    elsif ( $CAcodon =~ /CA[AG]/i)     { return 'Q' }    # Glutamine
    elsif ( $CAcodon =~ /CG.|AG[AG]/i) { return 'R' }    # Arginine
    elsif ( $CAcodon =~ /TC.|AG[TC]/i) { return 'S' }    # Serine
    elsif ( $CAcodon =~ /AC./i)        { return 'T' }    # Threonine
    elsif ( $CAcodon =~ /GT./i)        { return 'V' }    # Valine
    elsif ( $CAcodon =~ /TGG/i)        { return 'W' }    # Tryptophan
    elsif ( $CAcodon =~ /TA[TC]/i)     { return 'Y' }    # Tyrosine
    elsif ( $CAcodon =~ /TA[AG]|TGA/i) { return '*' }    # Stop
    else {
        return '?';
    }
}




### Get 3 Frame translation
### $prot=Frame3Translation($cds_seq)
### Global:
### Dependency: &Codon2AA
### Return $prot={ 0=> AA, 1 => AA, 2 => AA}
### Note: 0,1,2 forsword strand: shift 0/1/2 base from 5 end
sub Frame3Translation {
	my $FTcdsseq=shift;
	
	my $FTsubinfo='SUB(FastaKit::Frame3Translation)';
	my $FTaa={};
	
	for (my $FTframe=0; $FTframe<3; $FTframe++) {
		${$FTaa}{$FTframe}='';
		for (my $Ftseqstart=$FTframe; $Ftseqstart<=(length($FTcdsseq)-3); $Ftseqstart+=3) {
			my $FTcodon=substr($FTcdsseq, $Ftseqstart, 3);
			${$FTaa}{$FTframe}.=&Codon2AA($FTcodon);
		}
	}
	
	return $FTaa;
}



### Get 6 Frame translation
### $prot=Frame6Translation($cds_seq)
### Global:
### Dependency: &SeqRevComp, &Codon2AA
### Return $prot={ 0=> AA, 1 => AA, 2 => AA, 3 => AA, 4 => AA, 5 => AA}
### Note: 0,1,2 forsword strand: shift 0/1/2 base from 5 end
### Note: 3,4,5 reverse complement strand: shift 0/1/2 base from 5 end
sub Frame6Translation {
	my $FTcdsseq=shift;
	
	my $FTsubinfo='SUB(FastaKit::Frame3Translation)';
	my $FTaa={};
	
	for (my $FTframe=0; $FTframe<3; $FTframe++) {
		${$FTaa}{$FTframe}='';
		for (my $Ftseqstart=$FTframe; $Ftseqstart<(length($FTcdsseq)-2); $Ftseqstart+=3) {
			my $FTcodon=substr($FTcdsseq, $Ftseqstart, 3);
			${$FTaa}{$FTframe}.=&Codon2AA($FTcodon);
		}
	}
	
	$FTcdsseq=&SeqRevComp($FTcdsseq);
	for (my $FTframe=0; $FTframe<3; $FTframe++) {
		${$FTaa}{$FTframe+3}='';
		for (my $Ftseqstart=$FTframe; $Ftseqstart<(length($FTcdsseq)-2); $Ftseqstart+=3) {
			my $FTcodon=substr($FTcdsseq, $Ftseqstart, 3);
			${$FTaa}{$FTframe+3}.=&Codon2AA($FTcodon);
		}
	}
	
	return $FTaa;
}





### Calculate and get fasta parameters
### CountFasta(my.fasta, code)
### Global: 
### Dependency: 
### Note: comma delimited codes: Code 1 = total_num_seqs
### Return: (1/0, Assigned)
sub CountFasta {
	my ($CFinputfa, $CFcode)=@_;
	
	my $CFsubinfo="SUB(FastaKit::CountFasta)";
	my %CFcodehash=(); ### For code dependency in future
	my $CFnumseqs=0;
	local *CFINPUTFASTA;
	my @CFreturnarr=();
	
	unless (defined $CFcode and $CFcode=~/^\S+$/) {
		print STDERR  $CFsubinfo, "Error: invalid fasta count code\n";
		return $FastaKit_failure;
	}
	my @CTcodearr=split(/,/, $CFcode);
	
	close CFINPUTFASTA if (defined fileno(CFINPUTFASTA));
	unless (open (CFINPUTFASTA, "< $CFinputfa")) {
		print STDERR  $CFsubinfo, "Error: unable to open fasta file: $CFinputfa\n";
		return $FastaKit_failure;
	}
	while (my $CFline=<CFINPUTFASTA>) {
		chomp $CFline;
		if ($CFline=~/^>\S+/) {
			$CFnumseqs++;
		}
	}
	close CFINPUTFASTA;
	
	foreach my $CFindcode (@CTcodearr) {
		push (@CFreturnarr, $CFnumseqs) if ($CFindcode==1);
		### For development
	}
	return ($FastaKit_success, @CFreturnarr);
}



### check if duplicated sequence ID in fasta
### CheckFastaIdDup (my.fasta)
### Global: 
### Dependency: 
### Note:
### Return: 1=NO_duplicate; 0=Error_or_have duplicate
sub CheckFastaIdDup {
	my $CFIDfasta=shift;
	
	my $CFIDsubinfo="SUB(FastaKit::CheckFastaIdDup)";
	my %CFIDfastaID=();
	my $CFIDnum_duplicated_ID=0;
	local *CFIDFASTA;
	
	unless (defined $CFIDfasta and -s $CFIDfasta) {
		print STDERR $CFIDsubinfo, "Error: fasta file not found\n";
		return $FastaKit_failure;
	}
	close CFIDFASTA if (defined fileno(CFIDFASTA));
	unless (open (CFIDFASTA, "<$CFIDfasta")) {
		print STDERR $CFIDsubinfo, "Error: can not open fasta file: $CFIDfasta\n";
		return $FastaKit_failure;
	}
	while (my $CFIDline=<CFIDFASTA>) {
		if ($CFIDline=~/^>(\S+)\s*.*/) {
			if (exists $CFIDfastaID{$1}) {
				print "Duplicate\t".$1."\n";
				$CFIDnum_duplicated_ID++;
			}
			else {
				$CFIDfastaID{$1}++;
			}
		}
	}
	close CFIDFASTA;
	if ($CFIDnum_duplicated_ID>0) {
		print $CFIDsubinfo, "Info: total number of duplicates: $CFIDnum_duplicated_ID\n";
		return $FastaKit_failure;
	}
	elsif ($CFIDnum_duplicated_ID==0) {
		return $FastaKit_success;
	}
}



### Try 6 frame translation and select longest for each
sub Codon6EL2aa {
	###
}
### Try 6 frame translation and select longest for all
sub Codon6AL2aa {
	###
}
### Try 6 frame translation
sub Codon6aa {
	###
}


### Calculate the global identity of two sequences
### RunGgsearch ($seq1, $seq2, $path_ggsearch)
### Global:
### Dependency:
### Note: only two sequences
### Return: (0/1, $identity)


###GGEARCH
sub RunGgsearch {
	my ($RGfasta01, $RGfasta02, $RGpath_ggsearch)=@_;
	
	my $RGsubinfo="SUB(FastaKit::RunGgsearch)";
	$RGpath_ggsearch = 'ggsearch' unless (defined $RGpath_ggsearch);
	my $RGoutput = 'temp_lufuhao.out';
	local *RGGGSEARCHOUT;
	my $RGpc_id=-100;
	
	unlink $RGoutput if (defined $RGoutput and -e $RGoutput);
	
	unless (exec_cmd_return("$RGpath_ggsearch $RGfasta01 $RGfasta02 > $RGoutput")) {
		print STDERR $RGsubinfo, "Error: ggsearch running failed\n";
		return $FastaKit_failure;
	}
	unless (-s $RGoutput) {
		print STDERR $RGsubinfo, "Error: ggsearch output failed\n";
		return $FastaKit_failure;
	}
	close RGGGSEARCHOUT if (defined fileno(RGGGSEARCHOUT));
	unless (open(RGGGSEARCHOUT," < $RGoutput")) {
		print STDERR  $RGsubinfo, "ERROR: can not open ggsearch: $RGoutput\n";
		return $FastaKit_failure;
	}
	while(my $RGline=<RGGGSEARCHOUT>) {
		chomp $RGline;
		if ($RGline =~ /\% identity/) {
# eg. global/global (N-W) score: -158; 18.3% identity (45.9% similar) in 987 aa overlap (1-912:1-930)
			my @RGtemp=split(/\s+/,$RGline);
			if ($RGpc_id == -100) {# THERE MAY BE MULTIPLE ALIGNMENTS IN THE FILE, BUT WE WANT TO JUST TAKE THE TOP (BEST) ALIGNMENT
				$RGpc_id=$RGtemp[4]; # 18.3%
				if (substr($RGpc_id,length($RGpc_id)-1,1) eq '%') {
					chop($RGpc_id);
				}
				else {
					print STDERR $RGsubinfo, "ERROR: pc_id $RGpc_id\n";
					return $FastaKit_failure;
				}
			}
		}
		elsif ($RGline =~ /!! No sequences with E\(\) < \d*/) {# THE SEQUENCES CAN'T BE ALIGNED USING ggsearch
			$RGpc_id=0.0;
		}
	}
	close RGGGSEARCHOUT;
	
	return ($FastaKit_success, $RGpc_id);
}



### run global alignment using EMBOSS stretcher for protein/DNA
### RunEmbossStretcher ($seq1.fa, $seq2.fa, nucl/prot, $output, [$path_stretcher], [$additional_options])
### Global:
### Dependency: stretcher software in EMBOSS
### Note: return 0=failed 1=success
sub RunEmbossStretcher {
	my ($RESseq1, $RESseq2, $REStype, $RESoutput, $RESpath_stretcer, $RESoptions)=@_;

##Defaults
	my $RESsubinfo='SUB(FastaKit::RunEmbossStretcher)';
	$RESpath_stretcer='stretcher' unless (defined $RESpath_stretcer);
	my ($RESdatafile, $RESgapopen, $RESgapextend, $RESs1type, $RESs2type);
	my $REScmd='';
	$RESoptions='' unless (defined $RESoptions);

## Input and output
	unless (defined $RESseq1 and -s $RESseq1) {
		print STDERR $RESsubinfo, "Error: invalid seqfile 1\n";
		return $FastaKit_failure;
	}
	unless (defined $RESseq2 and -s $RESseq2) {
		print STDERR $RESsubinfo, "Error: invalid seqfile 2\n";
		return $FastaKit_failure;
	}
	if (defined $REStype) {
		if ($REStype=~/(^nucl$)/i) {
			$RESdatafile='EDNAFULL';
			$RESgapopen=16;
			$RESgapextend=4;
			$RESs1type=' -snucleotide1 ';
			$RESs2type=' -snucleotide2 '
		}
		elsif ($REStype=~/(^prot$)/i) {
			$RESdatafile='EBLOSUM62';
			$RESgapopen=12;
			$RESgapextend=2;
			$RESs1type=' -sprotein1 ';
			$RESs2type=' -sprotein2 '
		}
		else {
			print STDERR $RESsubinfo, "Error: invalid sequence type, should be either 'nucl' or 'prot'\n";
			return $FastaKit_failure;
		}
	}
	else {
		print STDERR $RESsubinfo, "Error: not defined sequence type, should be either 'nucl' or 'prot'\n";
		return $FastaKit_failure;
	}

##Running 
	$REScmd="$RESpath_stretcer -auto -asequence $RESseq1 -bsequence $RESseq2 -outfile $RESoutput -datafile $RESdatafile -gapopen $RESgapopen -gapextend $RESgapextend -aformat3 pair $RESs1type $RESs2type ";
	unless (exec_cmd($REScmd)) {
		print STDERR $RESsubinfo, "Error: running stretcher failed:\n$REScmd\n###\n";
		return $FastaKit_failure;
	}
	unless (-s $RESoutput) {
		print STDERR $RESsubinfo, "Error: stretcher output failed: $RESoutput\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}



### get similarity from global alignment generated by EMBOSS stretcher for protein/DNA
### AnalysisEmbossStretcherOutput ($seq1.fa, $seq2.fa, nucl/prot, $output, [$path_stretcher], [$additional_options])
### Global:
### Dependency:
### Note: return (0=failed/1=success, $identity, $similarity, $gaps, $score)
###	
sub AnalysisEmbossStretcherOutput {
	my $ARESOstretcher_output=shift;
	
	my $ARESOsubinfo='SUB(FastaKit::AnalysisEmbossStretcherOutput)';
	my $ARESOidentity='NaN';
	my $ARESOsimilarity='NaN';
	my $ARESOgaps='NaN';
	my $ARESOscore='NaN';
	my $ARESOcheck_stretcher_format=0;
	my $ARESOcheck_seq=0;
	my $ARESOseq1='';
	my $ARESOseq2='';
	local *ARESOSTRECHEROUTPUT;
	
	unless (defined $ARESOstretcher_output and -s $ARESOstretcher_output) {
		print STDERR $ARESOsubinfo, "Error: invalid EMBOSS output file\n";
		return $FastaKit_failure;
	}
	close ARESOSTRECHEROUTPUT if (defined fileno(ARESOSTRECHEROUTPUT));
	unless (open (ARESOSTRECHEROUTPUT, " < $ARESOstretcher_output")) {
		print STDERR $ARESOsubinfo, "Error: unable to open EMBOSS output file\n";
		return $FastaKit_failure;
	}
	while (my $ARESOline=<ARESOSTRECHEROUTPUT>) {
		chomp $ARESOline;
		if ($ARESOline=~/^#\s+Program:\s+stretcher$/) {
			$ARESOcheck_stretcher_format=1;
		}
		next unless ($ARESOcheck_stretcher_format==1);
		if ($ARESOline=~/^#\s+Aligned_sequences:\s+2$/) {
			$ARESOcheck_seq=1;
		}
		next unless ($ARESOcheck_seq==1);
		if ($ARESOline=~/^#\s+1:\s+(\S+)$/) {
			if ($ARESOseq1 eq '') {
				$ARESOseq1=$1;
			}
			else {
				print STDERR $ARESOsubinfo, "Warnings: multiple seq1 values in EMBOSS output file, only keep the first one\n";
			}
		}
		if ($ARESOline=~/^#\s+2:\s+(\S+)$/) {
			if ($ARESOseq2 eq '') {
				$ARESOseq2=$1;
			}
			else {
				print STDERR $ARESOsubinfo, "Warnings: multiple seq2 values in EMBOSS output file, only keep the first one\n";
			}
		}
		if ($ARESOline=~/^#\s+Identity:\s+(\d+\/\d+\s+\(\S+\%\))$/) {
			if ($ARESOidentity eq 'NaN') {
				$ARESOidentity=$1;
			}
			else {
				print STDERR $ARESOsubinfo, "Warnings: multiple identity values in EMBOSS output file, only keep the first one\n";
			}
		}
		if ($ARESOline=~/^#\s+Similarity:\s+(\d+\/\d+\s+\(\S+\%\))$/) {
			if ($ARESOsimilarity eq 'NaN') {
				$ARESOsimilarity=$1;
			}
			else {
				print STDERR $ARESOsubinfo, "Warnings: multiple similarity values in EMBOSS output file, only keep the first one\n";
			}
		}
		if ($ARESOline=~/^#\s+Gaps:\s+(\d+\/\d+\s+\(\S+\%\))$/) {
			if ($ARESOgaps eq 'NaN') {
				$ARESOgaps=$1;
			}
			else {
				print STDERR $ARESOsubinfo, "Warnings: multiple gaps values in EMBOSS output file, only keep the first one\n";
			}
		}
		if ($ARESOline=~/^#\s+Score:\s+(\S+)$/) {
			if ($ARESOscore eq 'NaN') {
				$ARESOscore=$1;
			}
			else {
				print STDERR $ARESOsubinfo, "Warnings: multiple score values in EMBOSS output file, only keep the first one\n";
			}
		}
	}
	close ARESOSTRECHEROUTPUT;
	return ($FastaKit_success, 	$ARESOidentity, $ARESOsimilarity, $ARESOgaps, $ARESOscore);
}



### Deduplicate fasta based on 100% sequence similarity, and will warn if different sequence have the same seqID
### FastaDedup ($input_fasta, $output_fasta)
### Global:
### Dependancy: Bio::SeqIO;
### Note:
sub FastaDedup {
	my ($FDfastain, $FDfastaout)=@_;
	
	my %FDidhash=();
	my %FDseqhash=();
	my $FDsubinfo='SUB(FastaKit::FastaDedup)';
	
	unless (defined $FDfastain and -s $FDfastain) {
		print STDERR $FDsubinfo, "Error: invalid input fasta\n";
		return $FastaKit_failure;
	}
	unless (defined $FDfastaout and $FDfastaout=~/^[\-\._\w\/\\]+$/){
		print STDERR $FDsubinfo, "Error: invalid output fasta\n";
		return $FastaKit_failure;
	}
	unlink $FDfastaout if (-e $FDfastaout);
	
	my $FDseqio  = Bio::SeqIO->new(-file => $FDfastain, -format => "fasta");
	my $FDoutseq = Bio::SeqIO->new(-file => " > $FDfastaout", -format => "fasta");
	
	while(my $FDseqs = $FDseqio->next_seq) {
		my $FDid  = $FDseqs->display_id;
		my $FDseq = $FDseqs->seq;
		unless (exists($FDseqhash{$FDseq})) {
			$FDoutseq->write_seq($FDseqs);
			$FDseqhash{$FDseq}++;
			if (exists $FDidhash{$FDid}) {
				print STDERR "Warnings: duplicated ID but different seq: $FDid\n";
			}
			$FDidhash{$FDid}++;
		}
	}
	
	return $FastaKit_success;
}




### Read [gzipped fasta get the length into hash]
### ReadFastaLength($fasta[.gz])
### 
sub ReadFastaLength {
	my $RFLfasta=shift;
	
	my $RFLsubinfo='SUB(FastaKit::ReadFastaLength)';
	my %RFLseqlength=();
	my $RFLcount=0;
	local *RFLFASTAIN;
	my $RFLseqid='';
	my $RFLseqidvlen=0;

	unless (defined $RFLfasta and -s $RFLfasta) {
		print STDERR $RFLsubinfo, "Error: invalid fasta\n";
		return $FastaKit_failure;
	}
	
	close RFLFASTAIN if (defined fileno(RFLFASTAIN));
	if ($RFLfasta=~/\.gz$/i) {
		unless (open (RFLFASTAIN, " gzip -cd $RFLfasta | ")) {
			print STDERR $RFLsubinfo, "Error: can not open gzipped fasta: $RFLfasta\n";
			return $FastaKit_failure;
		}
	}
	else {
		unless (open (RFLFASTAIN, " < $RFLfasta")) {
			print STDERR $RFLsubinfo, "Error: can not open fasta: $RFLfasta\n";
			return $FastaKit_failure;
		}
	}

	while (my $RFLline=<RFLFASTAIN>) {
		chomp $RFLline;
		if ($RFLline=~/^>/ or eof(RFLFASTAIN)) {
			if ($RFLcount>0 or eof(RFLFASTAIN)) {
				if ($RFLseqidvlen <=0) {
					print STDERR $RFLsubinfo, "Warnings: seqid $RFLseqid length $RFLseqidvlen <=0\n";
				}
				if (eof(RFLFASTAIN)) {
					unless ($RFLline=~/^>/) {
						$RFLseqidvlen = $RFLseqidvlen + length($RFLline);
					}
				}
				if (exists $RFLseqlength{$RFLseqid}) {
					print STDERR $RFLsubinfo, "Warnings: duplicated seqid: $RFLseqid\n";
					print Dumper \%RFLseqlength;
					return $FastaKit_failure;
				}
				else {
					$RFLseqlength{$RFLseqid}=$RFLseqidvlen;
				}
			}
			if ($RFLline=~/^>/) {
				if ($RFLline=~/^>(\S+)\s*/) {
					$RFLseqid=$1;
					$RFLcount++;
				}
				else {
					print STDERR $RFLsubinfo, "Error: invalid seq name (at line $.): $RFLline\n";
					return $FastaKit_failure;
				}
			}
			$RFLseqidvlen=0;
		}
		else {
			$RFLseqidvlen = $RFLseqidvlen + length($RFLline);
		}

	}
	close RFLFASTAIN;
	print $RFLsubinfo, "Sum: read total seqs ", scalar(keys %RFLseqlength), "\n";
	return ($FastaKit_success, \%RFLseqlength);
}




#$FastaKit_failure; $FastaKit_success;



1;
__END__
