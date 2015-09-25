# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::FastaKit

=head1 SYNOPSIS

Fasta -related tools

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
package FuhaoPerl5Lib::FastaKit;
use strict;
use warnings;
use Exporter;
use Cwd;
use FuhaoPerl5Lib::FileKit qw/MoveFile RetrieveDir MergeFiles/;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::MiscKit qw/IsReference FullDigit/;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150603';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(CdbFasta CdbYank IndexFasta CreateFastaRegion RunMira4 CdHitEst RenameFasta RunFqTrinity SplitFastaByNumber RunCap3 Fastq2Fasta SeqRevComp);
%EXPORT_TAGS = ( DEFAULT => [qw(&CdbFasta &CdbYank &IndexFasta &CreateFastaRegion &RunMira4 &RenameFasta &RunFqTrinity &SplitFastaByNumber RunCap3 Fastq2Fasta SeqRevComp)],
                 ALL    => [qw(&CdbFasta &CdbYank &IndexFasta &CreateFastaRegion &RunMira4 &RenameFasta &RunFqTrinity &SplitFastaByNumber RunCap3 Fastq2Fasta SeqRevComp)]);

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
### Dependancy: &exec_cmd_return
### Note:
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
		if (! &exec_cmd_return("$CFpath_cdbfasta $CSfasta_file -z $CSfasta_file.cdbz")) {
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
### Dependancy: &exec_cmd_return
### Note: 
sub CdbYank {
	my ($CYindex, $CYoutput, $CYseq_ids_index, $CYpath_cdbyank)=@_;
	
	$CYpath_cdbyank='cdbyank' unless (defined $CYpath_cdbyank);
	my $CYsubinfo='SUB(FastaKit::CdbYank)';
	
	unless (defined $CYindex and -s $CYindex) {
		print STDERR "${CYsubinfo}Error: invalid cdbfasta index: $CYindex\n";
		return $FastaKit_failure;
	}
	unlink $CYoutput if (-e $CYoutput);
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
		undef $CYseq_ids_index;
	}
	
	unless (defined $CYseqids_join and $CYseqids_join=~/\S+/) {
		print STDERR "${CYsubinfo}Error: invalid cdbfasta IDs: $CYseqids_join\n";
		return $FastaKit_failure;
	}
	
	my $CYcmd="perl -e 'print \"$CYseqids_join\\n\";' | $CYpath_cdbyank $CYindex -o $CYoutput -w";
	if (! &exec_cmd_return($CYcmd)) {
		print STDERR "${CYsubinfo}Error: cdbyank running error\n";
		return $FastaKit_failure;
	}
	elsif ( ! -s $CYoutput) {
		print STDERR "${CYsubinfo}Error: cdbyank output error\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}



### Index fasta using samtools index
### IndexFasta(input.fa, [path_samtools])
### Global:
### Dependency:
### Note:
sub IndexFasta {
	my ($IFfasta, $IFpath_samtools)=@_;
	
	my $IFsubinfo='SUB(FastaKit::IndexFasta)';
	$IFpath_samtools='samtools' unless (defined $IFpath_samtools);
	
	unless (defined $IFfasta and -s $IFfasta) {
		print STDERR "${IFsubinfo}Error: invalid fasta file\n";
		return $FastaKit_failure;
	}
	unlink "$IFfasta.fai" if (-e "$IFfasta.fai");
	
	if (! &exec_cmd_return("$IFpath_samtools faidx $IFfasta")) {
		print STDERR "${IFsubinfo}Error: samtools index running error\n";
		return $FastaKit_failure;
	}
	elsif (! -s "$IFfasta.fai") {
		print STDERR "${IFsubinfo}Error: samtools index output error\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}



sub GuessFormat {
	my $s = shift @_;
#	print $s."\n";     #For TEST
	$s =~ s/^.*\.(.+)$/$1/;
#	print $s."\n";   #For TEST
	my $failed = 0;
	my $seq_format;
 	SW: {
	if ($s =~ /(^fasta)|(^fast)|(^fst)|(^fsa)|(^ft)|(^fs)|(^fa)|(^fas)/i) {$seq_format = 'Fasta'; last SW};
	if ($s =~ /(^fastq)|(^fq)/i) {$seq_format = 'Fastq'; last SW};
	if ($s =~ /(lfasta)|(lfast)|(lfst)|(lfsa)|(lft)|(lfs)/i) {$seq_format = 'LabeledFasta'; last SW};
	if ($s =~ /(embl)|(emb)|(em)|(eml)/i) {$seq_format = 'EMBL'; last SW};
	if ($s =~ /(genebank)|(genbank)|(genb)|(geneb)|(gbank)|(gb)/i) {$seq_format = 'GenBank'; last SW};
	if ($s =~ /(swissprot)|(sprt)|(swissp)|(sprot)|(sp)|(spr)/i) {$seq_format = 'Swissprot'; last SW};
	if ($s =~ /pir/i) {$seq_format = 'PIR'; last SW};
	if ($s =~ /gcg/i) {$seq_format = 'GCG'; last SW};
	if ($s =~ /scf/i) {$seq_format = 'SCF'; last SW};
	if ($s =~ /ace/i) {$seq_format = 'Ace'; last SW};
	if ($s =~ /phd/i) {$seq_format = 'phd'; last SW};
	if ($s =~ /phred|phd/i) {$seq_format = 'phred'; last SW};
	if ($s =~ /raw/i) {$seq_format = 'raw'; last SW};
	if ($s =~ /bsml/i) {$seq_format = 'bsml'; last SW};
	$failed++;
	}
	return eval{$failed ? 0 : $seq_format};
}



### Create fasta region file for freebayes-parallel
### &CreateFastaRegion(xxx.fa/xxx.fai, bin, xxx.out, [path_samtools])
### Global: 
### Dependency: &exec_cmd_return
### Note: 
sub CreateFastaRegion {
	my ($CFRinput, $CFRbinzise, $CFRoutput, $CFpath_samtools)=@_;
	
	local *REGIONIN; local *REGIONOUT;
	my $CFsubinfo='SUB(FastaKit::CreateFastaRegion)';
	$CFpath_samtools='samtools' unless (defined $CFpath_samtools);
	
	unless (defined $CFRinput and -s $CFRinput) {
		print STDERR "${CFsubinfo}Info: invalid fasta index\n";
		return $FastaKit_failure;
	}
	if ($CFRinput !~ /\.fai$/i) {
		if ($CFRinput =~ /\.(fa)|(fasta)|(fas)$/i) {
			print "${CFsubinfo}Info: input format fasta\n" if ($FastaKit_debug); ### For test ###
			unless (-s "$CFRinput.fai") {
				print "${CFsubinfo}Info: $CFRinput.fai not exists; use samtools faidx to generates ..." if ($FastaKit_debug);### For test ###
				if (! &exec_cmd_return("$CFpath_samtools faidx $CFRinput")) {
					print STDERR "${CFsubinfo}Error: samtools faidx $CFRinput running error\n";
					return $FastaKit_failure;
				}
			}
			$CFRinput.='.fai';
			if (! -s $CFRinput) {
				print STDERR "${CFsubinfo}Error: samtools faidx $CFRinput output error\n";
				return $FastaKit_failure;
			}
		}
		else {
			print STDERR "${CFsubinfo}Error: unknown input format for samtools faidx error\n";
			return $FastaKit_failure;
		}
	}
	
	unless (defined $CFRoutput and $CFRoutput !~ /^\s*$/) {
		print STDERR "${CFsubinfo}Error: invalid region file output";
		return $FastaKit_failure;
	}
	unlink $CFRoutput if (defined $CFRoutput and -e $CFRoutput);
	
	close REGIONIN if (defined fileno(REGIONIN));
	unless (open(REGIONIN, "<$CFRinput")) {
		print STDERR "${CFsubinfo}Error: open $CFRinput error\n";
		return $FastaKit_failure;
	}
	close REGIONOUT if (defined fileno(REGIONOUT));
	unless (open(REGIONOUT, ">$CFRoutput")) {
		close REGIONIN;
		print STDERR "${CFsubinfo}Error: write $CFRoutput error\n";
		return $FastaKit_failure;
	}
	while (my $CFRline=<REGIONIN>) {
		chomp $CFRline;
		my @CFRarr=();
		@CFRarr=split(/\t/, $CFRline);
#		if (scalar(@CFRarr) != 5) {print STDERR "${CFsubinfo}Error: samtools faidx $CFRinput running error\n";
#			print STDERR "${CFsubinfo}Error: fai colnum != 5\n";
#			return $FastaKit_failure;
#		}###Force to check fai column number
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
			print STDERR "${CFsubinfo}Error: Wrong bin size\n";
			return $FastaKit_failure;
		}
	}
	close REGIONIN;
	close REGIONOUT;
	if (! -s $CFRoutput) {
		print STDERR "${CFsubinfo}Error: CreateFastaRegion output $CFRoutput not exists\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}



### running mira
### &RunMira4(fastq, /path/manifest_name, output_fasta, seq_prefix, seq_desc, [Min_alternative_count], [$RMpath_mira4], [num_threads])
### Global: 
### Dependency: $FastaKit_success $FastaKit_failure
### Note: Be sure to chdir back in case of any error
### Return: 0=failure; 1=success
sub RunMira4 {
	my ($RMfastq, $RMmira_manifest, $RMassembly_fasta, $RMseq_prefix, $RMmin_alternative_count, $RMpath_mira4, $RMnum_threads)=@_;
#	return $FastaKit_failure; ### For test RunCap3 ###
	local *MANIFEST; local *RMOUT; local *RMFA;
	my $RMsubinfo='SUB(FastaKit::RunMira)';
	my $RMreturn_fasta;
	$RMmin_alternative_count=3 unless (defined $RMmin_alternative_count);
	$RMpath_mira4='mira' unless (defined $RMpath_mira4);
	$RMnum_threads=1 unless (defined $RMnum_threads);
	my $RMcurpath=getcwd;
	$RMseq_prefix='MyHapAssem_' unless (defined $RMseq_prefix);
	
	unless (defined $RMmira_manifest and $RMmira_manifest=~/^\S+$/) {
		print STDERR "${RMsubinfo}Error: invalid MIRA4 manifest file\n";
		return $FastaKit_failure;
	}
	my $RMrundir=&RetrieveDir($RMmira_manifest);
	unless (-d $RMrundir) {
		print STDERR "${RMsubinfo}Error: can not retrieve manifest path for $RMmira_manifest\n";
		return $FastaKit_failure;
	}
#	print "\n${RMsubinfo}Test: \n###Fastq: $RMfastq\n###Manifest: $RMmira_manifest\n###Outfasta: $RMassembly_fasta\n###SeqPrefix: $RMseq_prefix\n###MinCount: $RMmin_alternative_count\n###MIRA4_path: $RMpath_mira4\n###Threads: $RMnum_threads\n\n";### For test ###
	
	if (! -s $RMfastq) {
		print STDERR "${RMsubinfo}Error: can not find fastq input\n";
		return $FastaKit_failure;
	}
	
	close MANIFEST if (defined fileno(MANIFEST));
	unless (open (MANIFEST, ">$RMmira_manifest")) {
		print STDERR "${RMsubinfo}Error: can not open manifest file: $RMmira_manifest\n";
		return $FastaKit_failure;
	}
	my $RMproject='Ta';
	my $RMjob='denovo,est,accurate';
	my $RMparameters="COMMON_SETTINGS -GENERAL:number_of_threads=$RMnum_threads -NAG_AND_WARN:cnfs=warn:check_template_problems=no:check_maxreadnamelength=no -CO:mark_repeats=yes:assume_snp_instead_repeat=yes:name_prefix=$RMseq_prefix -OUT:output_result_caf=no:output_result_tcs=no:output_result_maf=no SOLEXA_SETTINGS -CO:min_reads_per_group=$RMmin_alternative_count -AS:minimum_reads_per_contig=3";
	print MANIFEST "project = $RMproject\njob = $RMjob\nparameters = $RMparameters\n###Readgroup\n";
	print MANIFEST "readgroup = wheat\nautopairing\ndata = $RMfastq\ntechnology = solexa\ntemplate_size = 50 1000 autorefine\nsegment_placement = ---> <---\n";
	close MANIFEST;
	unless (chdir $RMrundir) {
		print STDERR "${RMsubinfo}Error: can not chdir to manifest path: $RMrundir\n";
		return $FastaKit_failure;
	}
	if (! &exec_cmd_return("$RMpath_mira4 $RMmira_manifest >mira.log 2>mira.err")) {
		print STDERR "${RMsubinfo}Error: MIRA return non-zero code\n";
		chdir $RMcurpath;
		return $FastaKit_failure;
	}
#	my @RMfasta_files=glob ("$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}*out.unpadded.fasta");#, "$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}_LargeContigs_out.unpadded.fasta");
	my @RMfasta_files=glob ("$RMrundir/${RMproject}_assembly/${RMproject}_*_results/${RMproject}*out.padded.fasta");
	if (scalar(@RMfasta_files)<1) {
		print STDERR "${RMsubinfo}Error: MIRA output empty\n";
		return $FastaKit_failure;
	}
	elsif (scalar(@RMfasta_files)==1) {
		$RMreturn_fasta=shift @RMfasta_files;
	}
	else {
		$RMreturn_fasta=$RMrundir."/merge.fa";
		close RMOUT unless (defined fileno(RMOUT));
		unless (open (RMOUT, ">$RMreturn_fasta")) {
			print STDERR "${RMsubinfo}Error: merge fasta1: open $RMreturn_fasta\n";
			chdir $RMcurpath;
			return $FastaKit_failure;
		}
		foreach (@RMfasta_files) {
			unless (open (RMFA, "<$_")) {
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
		print STDERR "${RMsubinfo}Error: move fasta error: $RMreturn_fasta to $RMassembly_fasta\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}



### RunCap3
### RunCap3(reads[.fq|.fasta], output.fasta, [path_cap3])
### Global: $FastaKit_failure, $FastaKit_success
### Dependency: 
### Note:
sub RunCap3 {
	my ($RCreads, $RCout, $RCpath_cap3)=@_;
	
	my $RCsubinfo='SUB(FastaKit::RunCap3)';
	$RCpath_cap3='cap3' unless (defined $RCpath_cap3);
	
	unless (-s $RCreads) {
		print STDERR "Error: invalid CAP3 input\n";
		return $FastaKit_failure;
	}
	if ($RCreads=~/(\.fasta$)|(\.fa$)/i) {
#		print $RCsubinfo, "Info: Cap3 input in fasta format\n"; ### For test ###
	}
	elsif ($RCreads=~/(\.fastq$)|(\.fq$)/i) {
#		print $RCsubinfo, "Info: Cap3 input in fasta format\n"; ### For test ###
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
	
	unless (exec_cmd_return("$RCpath_cap3 $RCreads > cap3.log 2>cap3.err")) {
		print STDERR $RCsubinfo, "Error: cap3 running error\n";
		return $FastaKit_failure;
	}
	my @RCoutputs = map{$RCreads.$_}(".cap.contigs", ".cap.singlets");
	print $RCsubinfo, "Test: @RCoutputs\n";
	unless (MergeFiles($RCout, @RCoutputs)) {
		print STDERR $RCsubinfo, "Error: merge CAP3 out error\n";
		return $FastaKit_failure;
	}
	if (-s $RCout) {
		return $FastaKit_success;
	}
	else {
		return $FastaKit_failure;
	}
}



### 
### 
### Global: 
### Dependency: $FastaKit_failure;$FastaKit_success;
### Note: 
sub Fastq2Fasta {
	my ($FFfq, $FFfa)=@_;
	
	local *FFFQIN; local *FFFAOUT;
	my $FFsubinfo='SUB(FastaKit::Fastq2Fasta)';
	
	unless (-s $FFfq) {
		print STDERR $FFsubinfo, "Error: invalid FASTQ: $FFfq\n";
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
	
	if (-s $FFfa) {
		return $FastaKit_success;
	}
	else {
		return $FastaKit_failure;
	}
}



### run cd-hit-est to derep the fasta file MIRA4 assembled
### &CdHitEst(input.fasta, output.fasta, add_cmd, [path_cdhitest])
### Global:
### Dependency:$FastaKit_failure;$FastaKit_success;
### Note: $CHEaddition_cmd='-c 1.00 -n 10 -T 0 -r 1 -d 0 -M 30000';
### Return: 0=failure, 1=success
sub CdHitEst {
	my ($CHEfastain, $CHEfastaout, $CHEaddition_cmd, $CHEpath_cdhitest)=@_;
	
	my $CHEsubinfo='SUB(FastaKit::CdHitEst)';
	$CHEpath_cdhitest='cd-hit-est' unless (defined $CHEpath_cdhitest);
	$CHEaddition_cmd='';
	unless (defined $CHEfastain and -s $CHEfastain) {
		print STDERR "${CHEsubinfo}Error: fasta input for CDHIT not defined or exists\n";
		return $FastaKit_failure;
	}
	unless (defined $CHEfastaout) {
		print STDERR "${CHEsubinfo}Error: fasta output for CDHIT not defined\n";
		return $FastaKit_failure;
	}
	if (-e $CHEfastaout) {
		print STDERR "${CHEsubinfo}Warnings: fasta output for CDHIT existed but deleted\n";
		unlink $CHEfastaout;
	}
	
	if (! &exec_cmd_return("$CHEpath_cdhitest -i $CHEfastain -o $CHEfastaout $CHEaddition_cmd")) {
		print STDERR "${CHEsubinfo}Error: CDHIT running error: $CHEfastain\n";
		return $FastaKit_failure;
	}
	elsif (-s $CHEfastaout) {
		return $FastaKit_success;
	}
	else {
		print STDERR "${CHEsubinfo}Error: CDHIT output error: $CHEfastain\n";
		return $FastaKit_failure;
	}
}



### Run Trinity v2.0.6 for assembly
### RunFqTrinity($FastqIn, $fastaout, $RQTadd_cmd, $path_trinity)
### $RQTadd_cmd='--max_memory 2G --run_as_paired --CPU 1 --group_pairs_distance 800 --full_cleanup --min_kmer_cov 3 --min_glue 3'
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: 
### Note:
### Return: 0=failure; 1=success
sub RunFqTrinity {
	my ($RQTfastq, $RQTfastaout, $RQTadd_cmd, $path_trinity)=@_;
	
	my $RQTsubinfo='SUB(FastaKit::RunTrinity)';
	$path_trinity='Trinity' unless (defined $path_trinity);
	$RQTadd_cmd=' ' unless (defined $RQTadd_cmd);
	my $RQTseqtype='fq';
	my $RQToutput='Trinity';
	my $RQTcurdir=getcwd;
	
	unless (defined $RQTfastq and -s $RQTfastq) {
		print STDERR "{$RQTsubinfo}Error: invalid FastQ/Fasta input\n";
		return $FastaKit_failure;
	}
	
	if ($RQTfastq=~/(\.fq$)|(\.fastq$)/i) {
		$RQTseqtype='fq';
	}
	elsif ($RQTfastq=~/(\.fa$)|(\.fasta$)|(\.fas$)/i) {
		$RQTseqtype='fa';
	}
	else {
		print STDERR "{$RQTsubinfo}Error: guess FastQ/Fasta format error\n";
		return $FastaKit_failure;
	}
	
	unless (exec_cmd_return("$path_trinity --seqType $RQTseqtype --single $RQTfastq $RQTadd_cmd --full_cleanup --output $RQToutput > trinity.log 2> trinity.err")) {
		print STDERR $RQTsubinfo, "Error: Trinity running error\n";
		return $FastaKit_failure;
	}
	my @RQTfasta_files=glob "$RQTcurdir/$RQToutput.*.f*a";
	unless (scalar(@RQTfasta_files)==1) {
		print STDERR "{$RQTsubinfo}Error: Trinity output detecting error: @RQTfasta_files\n";
		return $FastaKit_failure;
	}
	unlink $RQTfastaout if (-e $RQTfastaout);
	unless (MoveFile($RQTfasta_files[0], $RQTfastaout)) {
		print STDERR "{$RQTsubinfo}Error: Trinity output rename error: $RQTfasta_files[0] to $RQTfastaout\n";
		return $FastaKit_failure;
	}
	return $FastaKit_success;
}



### Rename fasta using regex
### RenameFasta($fastain, $fastqout, $prefix, $num_digit, $RFDesc)
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: 
### Note:
### Return: 0=failure; 1=success
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
		print STDERR "${RFsubinfo}Error: invalid Fasta input\n";
		return $FastaKit_failure;
	}
	unless (defined $RFfaout and $RFfaout=~/^\S+$/) {
		print STDERR "${RFsubinfo}Error: invalid Fasta output\n";
		return $FastaKit_failure;
	}
	unlink $RFfaout if (-e $RFfaout);
	close RFFASTAIN if (defined fileno(RFFASTAIN));
	unless (open (RFFASTAIN, "< $RFfastain")) {
		print STDERR "${RFsubinfo}Error: can not open fasta: $RFfastain\n";
		return $FastaKit_failure;
	}
	close RFFASTAOUT if (defined fileno(RFFASTAOUT));
	unless (open(RFFASTAOUT, "> $RFfaout")) {
		print STDERR "${RFsubinfo}Error: can not write fasta: $RFfaout\n";
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
				print STDERR "${RFsubinfo}Error: invalid fastaid: $RFline\n";
				close RFFASTAOUT;
				unlink $RFfaout if (-e $RFfaout);
				close RFFASTAIN;
				return $FastaKit_failure;
			}
			if (exists $RFcheckid_duplicated{$RFnew_seqid}) {
				print STDERR "${RFsubinfo}Error: seqid duplicated: $RFnew_seqid\n";
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
### SplitFastaByNumber(input.fasta, number, output.prefix)
### Global: $FastaKit_failure; $FastaKit_success;
### Dependency: 
### Note: 
sub SplitFastaByNumber {
	my ($SFBNfasta, $SFBNnum, $SFBNoutprefix)=@_;
	
	local *SFBNFASTAIN; local *SFBNFASTAOUT;
	my $SFBNsubinfo='SUB(FastaKit::SplitFastaByNumber)';
	$SFBNnum=100 unless (defined $SFBNnum and $SFBNnum=~/^\d+$/ and $SFBNnum>0);
	$SFBNoutprefix=$SFBNfasta unless (defined $SFBNoutprefix);
	
	unless (defined $SFBNfasta and -s $SFBNfasta) {
		print STDERR "${SFBNsubinfo}Error: invalid fasta input\n";
		return $FastaKit_failure;
	}
	my $SFBNcount=0;
	my $SFBNcountfileno=0;
	my @SFBNfilelist=();
	close SFBNFASTAIN if (defined fileno(SFBNFASTAIN));
	unless (open (SFBNFASTAIN, "< $SFBNfasta")) {
		print STDERR "${SFBNsubinfo}Error: can not open fasta: $SFBNfasta\n";
		return $FastaKit_failure;
	}
	unless (open (SFBNFASTAOUT, ">$SFBNoutprefix.$SFBNcountfileno.fa")) {
		print STDERR "${SFBNsubinfo}Error: can not write fasta: $SFBNoutprefix.$SFBNcountfileno.fa\n";
		return $FastaKit_failure;
	}
	while (my $SFBNline=<SFBNFASTAIN>) {
		if ($SFBNline=~/^>/) {
			$SFBNcount++;
			if ($SFBNcount>$SFBNnum) {
				close SFBNFASTAOUT if (defined fileno(SFBNFASTAOUT));
				if (-s "$SFBNoutprefix.$SFBNcountfileno.fa") {
					push (@SFBNfilelist, "$SFBNoutprefix.$SFBNcountfileno.fa")
				}
				$SFBNcountfileno++;
				unlink ("$SFBNoutprefix.$SFBNcountfileno.fa") if (-e "$SFBNoutprefix.$SFBNcountfileno.fa");
				unless (open (SFBNFASTAOUT, ">$SFBNoutprefix.$SFBNcountfileno.fa")) {
					print STDERR "${SFBNsubinfo}Error: can not write fasta: $SFBNoutprefix.$SFBNcountfileno.fa\n";
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
	close SFBNFASTAOUT if (defined fileno(SFBNFASTAOUT));
	return $FastaKit_success;
}



### Seq Reverse Complement
### SeqRevComp($seq_string)
### Global: 
### Dependency: 
### Note:
sub SeqRevComp {
	my $SRCoriseq=shift;
	my $SRCseq_rev=reverse $SRCoriseq;
	$SRCseq_rev=~tr/actgACTG/tgacTGAC/;
	return $SRCseq_rev;
}



### Codon2AA
### &Code2AA('str')
### Global:
### Dependency: 
### Note:
sub Codon2AA {
	my $CAcodon=shift;
	
	$CAcodon=UC $CAcodon;
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
    my($codon) = @_;
	
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



1;
