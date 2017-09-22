# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::GffKit

=head1 SYNOPSIS

GFF -related tools

=head1 Requirements

Perl Modules:

=head1 DESCRIPTION

=over 2

=item GffReverseStrand (input.gff3, reference.fa/fai, output.gff3, samtools_path)

    * Convert GFF file to reverse complement strand
    * Return: 1=Sucesss    0=Failure
    * Note: reference.fa needs samtools while reference.fa.fai donot


=item WriteGff3 (output.gff3, \%(gene))

=item ExonerateGff3Reformat (\$cDNA_gff3_input, \$CDS_gff3_input, \$file_gff3_output)

    * Reformat exonerate GFF3
    * Return: 1=success, 0=failure

=back

=head1 FEEDBACK

=head2 Support

Please send you questions or bug reports to Email:

I<lufuhao@gmail.com>

=head1 AUTHORS - Fu-Hao Lu

Email: lufuhao@gmail.com (Always)
       Fu-Hao.Lu@jic.ac.uk (2012-2018)

=head1 CONTRIBUTORS

None

=head1 APPENDIX

Fight against Bioinformatics with Perl ^_^

=cut

#Coding starts
package FuhaoPerl5Lib::GffKit;
use strict;
use warnings;
use Exporter;
use Cwd;
use FuhaoPerl5Lib::FastaKit qw/IndexFasta Frame3Translation SeqRevComp/;
use Bio::DB::Fasta;
use Data::Dumper qw/Dumper/;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20170914';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(GffReverseStrand ReadGff3 WriteGff3 ExonerateGff3Reformat);
%EXPORT_TAGS = ( DEFAULT => [qw(GffReverseStrand ReadGff3 WriteGff3 ExonerateGff3Reformat)],
                 ALL    => [qw(GffReverseStrand ReadGff3 WriteGff3 ExonerateGff3Reformat)]);

my $GffKit_success=1;
my $GffKit_failure=0;
my $GffKit_debug=0;


### Convert GFF file to reverse complement strand
### $test_gffrc=&GffReverseStrand(GFF_input, reference.fa/fai, GFF_output, samtools_path)
### Global:
### Dependancy: 
### Note: 1. samools path needed if using fasta and index not exists
###

sub GffReverseStrand {
	my ($GRSgffin, $GRSfasta, $GRSgffout, $GRSpath_samtools)=@_;
	
	my $GRSsubinfo='SUB(GffKit::GffReverseStrand)';
	my %GRSseqlength=();
	$GRSpath_samtools='samtools' unless (defined $GRSpath_samtools);
	local *GRSINDEX; local *GRSGFFIN; local *GRSGFFOUT;

### Check input and output
	unless (defined $GRSgffin and -s $GRSgffin) {
		print STDERR $GRSsubinfo, "Error: invalid GFF input\n";
		return $GffKit_failure;
	}
	unless (defined $GRSfasta and -s $GRSfasta) {
		print STDERR $GRSsubinfo, "Error: invalid fasta (index)\n";
		return $GffKit_failure;
	}
	unless (defined $GRSgffout and $GRSgffout=~/^\S+$/) {
		print STDERR $GRSsubinfo, "Error: invalid GFF output name\n";
		return $GffKit_failure;
	}
	unlink $GRSgffout if (-e $GRSgffout);
	close GRSINDEX if (defined fileno(GRSINDEX));
	close GRSGFFIN if (defined fileno(GRSGFFIN));
	close GRSGFFOUT if (defined fileno(GRSGFFOUT));

### check if fasta or fai, only need fasta index file
	if ($GRSfasta=~/\.fai$/i) {
		print $GRSsubinfo, "Info: Fasta index detected, continue...\n";
	}
	elsif ($GRSfasta=~/(\.fa$)|(\.fas$)|(\.fasta$)/i) {
		print $GRSsubinfo, "Info: Fasta detected, need to check existing index\n";
		if (-e "$GRSfasta.fai") {
			print $GRSsubinfo, "Info: Fasta detected, use existing index\n";
			$GRSfasta.='.fai';
		}
		else {###Need to create fasta index
			print $GRSsubinfo, "Info: Fasta detected, and no existing index\n";
			unless (IndexFasta($GRSfasta, $GRSpath_samtools)) {
				print STDERR $GRSsubinfo, "Error: creating fasta index failed\n";
				return $GffKit_failure;
			}
			unless (-s "$GRSfasta.fai") {
				print STDERR $GRSsubinfo, "Error: samtools faidx output not exists\n";
				return $GffKit_failure;
			}
			$GRSfasta.='.fai';
		}
	}

	unless(open (GRSINDEX, "< $GRSfasta")) {
		print STDERR $GRSsubinfo, "Error: can not open fasta index file: $GRSfasta\n";
		return $GffKit_failure;
	}
	while (my $GRSline=<GRSINDEX>) {
		chomp $GRSline;
		my @GRSarr=split(/\t/, $GRSline);
		if (defined $GRSarr[0] and $GRSarr[0]=~/^\S+$/ and defined $GRSarr[1] and $GRSarr[1]=~/^\d+$/) {
			if (exists $GRSseqlength{$GRSarr[0]}) {
				print STDERR $GRSsubinfo, "Warnings: duplicated seqID: $GRSarr[0], checking if the same length\n";
				if ($GRSseqlength{$GRSarr[0]}=~/^\d+$/ and $GRSseqlength{$GRSarr[0]}==$GRSarr[1]) {
					print STDERR $GRSsubinfo, "Warnings: duplicated seqID: $GRSarr[0], has the same length\n";
				}
				else {
					print STDERR $GRSsubinfo, "Error: duplicated seqID: $GRSarr[0], has different length\n";
					return $GffKit_failure;
				}
			}
			else {
				$GRSseqlength{$GRSarr[0]}=$GRSarr[1];
			}
		}
		else {
			print STDERR $GRSsubinfo, "Error: invalid fasta index format: $GRSfasta\n";
			return $GffKit_failure;
		}
	}
	close GRSINDEX;


	unless (open (GRSGFFIN, " < $GRSgffin")) {
		print STDERR $GRSsubinfo, "Error: can not open GFF input: $GRSgffin\n";
		return $GffKit_failure;
	}
	unless (open (GRSGFFOUT, " > $GRSgffout")) {
		print STDERR $GRSsubinfo, "Error: can not write GFF output: $GRSgffin\n";
		return $GffKit_failure;
	}
	while (my $GRSline=<GRSGFFIN>) {
		chomp $GRSline;
		if ($GRSline=~/^#/) {
			print GRSGFFOUT $GRSline, "\n";
			next;
		}
		my @GRSarr=split(/\t/, $GRSline);
		unless (defined $GRSarr[3] and $GRSarr[3]=~/^\d+$/ and defined $GRSarr[4] and $GRSarr[4]=~/^\d+$/) {
			print STDERR $GRSsubinfo, "Error: invalid line ($.): $GRSline\n";
			return $GffKit_failure;
		}
		unless (exists $GRSseqlength{$GRSarr[0]} and $GRSseqlength{$GRSarr[0]}=~/^\d+$/) {
			print STDERR $GRSsubinfo, "Error: seqID $GRSarr[0] in GFF input not found in fasta index at line ($.): $GRSline\n";
			return $GffKit_failure;
		}
		if ($GRSarr[6] eq '+') {
			$GRSarr[6]='-';
		}
		elsif ($GRSarr[6] eq '-') {
			$GRSarr[6]='+';
		}
		else {
			print STDERR $GRSsubinfo, "Error: unknown strand at line ($.): $GRSline\n";
			return $GffKit_failure;
		}
		my $temp=$GRSarr[3];
		$GRSarr[3]=$GRSseqlength{$GRSarr[0]}+1-$GRSarr[4];
		$GRSarr[4]=$GRSseqlength{$GRSarr[0]}+1-$temp;
		print GRSGFFOUT join("\t", @GRSarr), "\n";
	}
	close GRSGFFIN;
	close GRSGFFOUT;
	return $GffKit_success;
}



### Read GFF3 into hash
### ($GffKit_success, $RGreferenceids, $RGgene, $RGgene2mrna, $RGmrnas, $RGexons, $RGcds)=ReadGff3($RGgffin, $RGfasta)
###
### %referenceids  => ( $reference_id => $gene_start_pos => $gene_id => num++)
### %gene2mrna     => ( $gene_id => $mrna_id => num++ )
### %gene=($geneid => ('reference' => $arr[0],
###                    'start'     => $arr[3], 
###                    'end'       => $arr[4],
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5],
###                    'Note'      => $note ###Not necessarily exist
###                   )
###       )
### %mrnas=($geneid => ('reference' => $arr[0],
###                     'start'      => $arr[3], 
###                     'end'        => $arr[4],
###                     'strand'     => $arr[6],
###                     'score'      => $arr[5],
###                     'Note'       => $note ###Not necessarily exist
###                    )
###       )
### %exon=($mrnaid => ('reference' => $arr[0],
###                    'exon' => ({$arr[3]} => ($arr[4] => $exonid)),
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5]
###                    )
###       )
### %cds=($mrnaid => ('reference' => $arr[0],
###                   'cds'       => ({$arr[3]} => ($arr[4] => num++)),
###                   'strand'    => $arr[6],
###                   'score'     => $arr[5],
###                   'phase'     => ({$arr[3]} => ($arr[4] => $arr[7]))
###                   )
###       )
sub ReadGff3 {
	my ($RGgffin, $RGfasta)=@_;
	
	my $RGsubinfo='SUB(GffKit::ReadGff3)';
	my $RGreferenceids={};
	my $RGgene={};
	my $RGcds={};
	my $RGgene2mrna={};
	my $RGmrnas={};
	my $RGexons={};
	my %RGmrna2gene=();
	my %RGmrna_border=();
	my $RGcheck_phase=0;
	my $RGseqdb;
	local *RGGFFIN;
	
	unless (defined $RGgffin and -s $RGgffin) {
		print STDERR $RGsubinfo, "Error: invalid GFF3 input\n";
		return $GffKit_failure;
	}
	if (defined $RGfasta and -s $RGfasta) {
		$RGseqdb=Bio::DB::Fasta->new($RGfasta);
		$RGcheck_phase=1;
	}
	close RGGFFIN if (defined fileno(RGGFFIN));
	if ($RGgffin=~/gff.*.gz/i) {
		unless (open (RGGFFIN, " gzip -dc $RGgffin | ")) {
			print STDERR $RGsubinfo, "Error: can not open gzippped GFF3 input: $RGgffin\n";
			return  $GffKit_failure;
		}
	}
	elsif ($RGgffin=~/(\.gff$)|(\.gff3$)/i) {
		unless (open (RGGFFIN, "< $RGgffin")) {
			print STDERR $RGsubinfo, "Error: can not open GFF3 input: $RGgffin\n";
			return  $GffKit_failure;
		}
	}
	else{
		print STDERR $RGsubinfo, "Error: can not guess GFF format, must be suffexed with .gff/.gff3\n";
		return  $GffKit_failure;
	}
	print $RGsubinfo, "Test: Reading GFF3: $RGgffin\n" if ($GffKit_debug); ### For test ###
	
	while (my $RGline=<RGGFFIN>) {
		chomp $RGline;
		next if ($RGline=~/^#/);
		my @RGarr=split(/\t/, $RGline);
		unless (scalar(@RGarr)==9) {
			print STDERR $RGsubinfo, "Error: invalid GFF line (at $.): $RGline\n";
			return  $GffKit_failure;
		}
		unless ($RGarr[3]=~/^\d+$/ and $RGarr[4]=~/^\d+$/ and $RGarr[3]<=$RGarr[4]) {
			print STDERR $RGsubinfo, "Error: border error at line $.: $RGline\n";
			return  $GffKit_failure;
		}
		
		unless (defined $RGarr[6] and $RGarr[6]=~/(^\+$)|(^\-$)/) {
			print STDERR $RGsubinfo, "Warnings: unknown strand at line($.): ".$RGline."\n" ;
			next;
		}
		if ($RGarr[2] =~ /^gene$/i) {
			my $RGthisgeneid=$RGarr[8];
			$RGthisgeneid=~s/^.*ID=//;$RGthisgeneid=~s/;.*$//;
			if (exists ${$RGgene}{$RGthisgeneid}) {
				print STDERR $RGsubinfo, "Error: duplicated gene ID: $RGthisgeneid\n";
				return  $GffKit_failure;
			}
			${$RGreferenceids}{$RGarr[0]}{$RGarr[3]}{$RGthisgeneid}++;
			${$RGgene}{$RGthisgeneid}{'reference'}=$RGarr[0];
			${$RGgene}{$RGthisgeneid}{'start'}=$RGarr[3];
			${$RGgene}{$RGthisgeneid}{'end'}=$RGarr[4];
			${$RGgene}{$RGthisgeneid}{'strand'}=$RGarr[6];
			${$RGgene}{$RGthisgeneid}{'score'}=$RGarr[5];
			if ($RGarr[8]=~/Note=(\S+)/) {
				my $RGnotes=$1;
				$RGnotes=~s/;.*$//;
				$RGnotes=~s/\s+$//;
				$RGnotes=~s/^\s+//;
				${$RGgene}{$RGthisgeneid}{'note'}=$RGnotes;
			}
#			print $RGsubinfo, "Test: geneid: ",$RGthisgeneid , "\tgenenote: ",${$RGgene}{$RGthisgeneid}{'note'},  "\n" if ($GffKit_debug); ### For test ###
		}
		elsif ($RGarr[2] =~ /^mRNA$/i) {
			my $RGthismrnaid=$RGarr[8];
			$RGthismrnaid=~s/^.*ID=//; $RGthismrnaid=~s/;.*$//;
			if (exists ${$RGmrnas}{$RGthismrnaid}) {
				print STDERR $RGsubinfo, "Error: duplicated mRNA ID: $RGthismrnaid\n";
				return  $GffKit_failure;
			}
			${$RGmrnas}{$RGthismrnaid}{'start'}=$RGarr[3];
			${$RGmrnas}{$RGthismrnaid}{'end'}=$RGarr[4];
			${$RGmrnas}{$RGthismrnaid}{'strand'}=$RGarr[6];
			${$RGmrnas}{$RGthismrnaid}{'score'}=$RGarr[5];
			${$RGmrnas}{$RGthismrnaid}{'reference'}=$RGarr[0];
			my $RGparent=$RGarr[8];
			$RGparent=~s/^.*Parent=//; $RGparent=~s/;.*$//;
			my @RGtemparr=split(/,/, $RGparent);
			foreach my $RGindpar (@RGtemparr) {
				${$RGgene2mrna}{$RGindpar}{$RGthismrnaid}++;
				$RGmrna2gene{$RGthismrnaid}{$RGindpar}++;
			}
			if ($RGarr[8]=~/Note=(\S+)/) {
				my $RGnotes=$1;
				$RGnotes=~s/;.*$//;
				$RGnotes=~s/\s+$//;
				$RGnotes=~s/^\s+//;
				${$RGmrnas}{$RGthismrnaid}{'note'}=$RGnotes;
			}
		}
		elsif ($RGarr[2] =~ /^exon$/i) {
			my $RGparent=$RGarr[8];
			$RGparent=~s/^.*Parent=//; $RGparent=~s/;.*$//;
#			print "Test: exon $RGarr[3]-$RGarr[4] \tParent: $RGparent\n"  if ($GffKit_debug); ### For test ###
			my @RGtemparr=split(/,/, $RGparent);
			my $RGexonid=$RGarr[8];
			$RGexonid=~s/^.*ID=//; $RGexonid=~s/;.*$//;
			foreach my $RGindpar (@RGtemparr) {
				if (exists ${$RGexons}{$RGindpar} and exists ${$RGexons}{$RGindpar}{'exon'} and exists ${$RGexons}{$RGindpar}{'exon'}{$RGarr[3]} and exists ${$RGexons}{$RGindpar}{'exon'}{$RGarr[3]}{$RGarr[4]}) {
					print STDERR $RGsubinfo, "Error: repeated exons ".$RGarr[3].' - '.$RGarr[4]. " for mRNA $RGindpar\n";
					return  $GffKit_failure;
				}
				${$RGexons}{$RGindpar}{'exon'}{$RGarr[3]}{$RGarr[4]}=$RGexonid;
				if (exists ${$RGexons}{$RGindpar}{'strand'}) {
					if (${$RGexons}{$RGindpar}{'strand'} ne $RGarr[6]) {
#						print STDERR $RGsubinfo, "Error: exon strand problem at line($.): \n$RGline\n";### For test ###
						print $RGsubinfo, "Error: inconsistant exon strand at line($.) : ", ${$RGexons}{$RGindpar}{'strand'}, "\t", $RGarr[6], "\n";
						return $GffKit_failure;
					}
				}
				else {
					${$RGexons}{$RGindpar}{'strand'}=$RGarr[6];
				}
				if ($RGarr[5]=~/^\d+\.*\d*$/) {
					if (exists ${$RGexons}{$RGindpar}{'score'}) {
						if (${$RGexons}{$RGindpar}{'score'}=~/^\d+\.*\d*$/){
							if ($RGarr[5] > ${$RGexons}{$RGindpar}{'score'}) {
								${$RGexons}{$RGindpar}{'score'}=$RGarr[5];
							}
						}
						${$RGexons}{$RGindpar}{'score'}=$RGarr[5];
					}
					else {
						${$RGexons}{$RGindpar}{'score'}=$RGarr[5];
					}
				}
				else {
					${$RGexons}{$RGindpar}{'score'}=$RGarr[5];
				}
				if (${$RGexons}{$RGindpar}{'reference'}) {
					if (${$RGexons}{$RGindpar}{'reference'} ne $RGarr[0]) {
						print STDERR $RGsubinfo, "Error: exon ref problem at line($.): \n$RGline\n";
						return $GffKit_failure;
					}
				}
				else {
					${$RGexons}{$RGindpar}{'reference'}=$RGarr[0];
				}
			}
		}
		elsif ($RGarr[2] =~ /^CDS$/i) {
			unless (defined $RGarr[7] and $RGarr[7]=~/^[0-2.]{1,1}$/) {
				print STDERR $RGsubinfo, "Error: CDS phase problem at line($.): \n$RGline\n";
				return $GffKit_failure;
			}
			my $RGparent=$RGarr[8];
			$RGparent=~s/^.*Parent=//; $RGparent=~s/;.*$//;
			print "Test: CDS $RGarr[3]-$RGarr[4] \tParent: $RGparent\n" if ($GffKit_debug);### For test ###
			my @RGtemparr=split(/,/, $RGparent);
			foreach my $RGindpar (@RGtemparr) {
				${$RGcds}{$RGindpar}{'cds'}{$RGarr[3]}{$RGarr[4]}++;
				if (exists ${$RGcds}{$RGindpar}{'strand'}) {
					if (${$RGcds}{$RGindpar}{'strand'} ne $RGarr[6]) {
						print STDERR $RGsubinfo, "Error: CDS strand problem at line($.): \n$RGline\n";
						return $GffKit_failure;
					}
				}
				else {
					${$RGcds}{$RGindpar}{'strand'}=$RGarr[6];
				}
				if ($RGarr[5]=~/^\d+\.*\d*$/) {
					if (exists ${$RGcds}{$RGindpar}{'score'}) {
						if (${$RGcds}{$RGindpar}{'score'}=~/^\d+\.*\d*$/){
							if ($RGarr[5] > ${$RGcds}{$RGindpar}{'score'}) {
								${$RGcds}{$RGindpar}{'score'}=$RGarr[5];
							}
						}
						${$RGcds}{$RGindpar}{'score'}=$RGarr[5];
					}
					else {
						${$RGcds}{$RGindpar}{'score'}=$RGarr[5];
					}
				}
				else {
					${$RGcds}{$RGindpar}{'score'}=$RGarr[5];
				}
				if (${$RGcds}{$RGindpar}{'reference'}) {
					if (${$RGcds}{$RGindpar}{'reference'} ne $RGarr[0]) {
						print STDERR $RGsubinfo, "Error: CDS ref problem at line($.): \n$RGline\n";
						return $GffKit_failure;
					}
				}
				else {
					${$RGcds}{$RGindpar}{'reference'}=$RGarr[0];
				}
				if (exists ${$RGcds}{$RGindpar}{'phase'} and exists ${$RGcds}{$RGindpar}{'phase'}{$RGarr[3]} and exists ${$RGcds}{$RGindpar}{'phase'}{$RGarr[3]}{$RGarr[4]}) {
					print STDERR $RGsubinfo, "Error: Repeated CDS at line($.): \n$RGline\n";
					return $GffKit_failure;
				}
				${$RGcds}{$RGindpar}{'phase'}{$RGarr[3]}{$RGarr[4]}=$RGarr[7];
			}
		}
	}
	close RGGFFIN;



### Check
### 1. check if each mRNA get a unique gene ID
	foreach my $RGind_mrna (sort keys %RGmrna2gene) {
		if (scalar(keys %{$RGmrna2gene{$RGind_mrna}})==0) {
			print STDERR $RGsubinfo, "Error: mRNA get no Parent: $RGind_mrna\n";
			return $GffKit_failure;
		}
		elsif (scalar(keys %{$RGmrna2gene{$RGind_mrna}})>1) {
			print STDERR $RGsubinfo, "Error: mRNA gets invalid gene PARENT: $RGind_mrna\n";
			return $GffKit_failure;
		}
		elsif (scalar(keys %{$RGmrna2gene{$RGind_mrna}})==1) {
			my @RGtemparr=keys %{$RGmrna2gene{$RGind_mrna}};
			unless (exists ${$RGgene}{$RGtemparr[0]}) {
				print STDERR $RGsubinfo, "Error: mRNA get no GENE: $RGind_mrna\n";
				print Dumper $RGgene;
				print Dumper \%RGmrna2gene;
				return $GffKit_failure;
			}
			$RGmrna2gene{$RGind_mrna}=$RGtemparr[0];
		}
	}
### 2. check if each gene gets a mRNA and gene border
	foreach my $RGind_gene (sort keys %{$RGgene}) {
		unless (exists ${$RGgene2mrna}{$RGind_gene} and scalar(keys %{${$RGgene2mrna}{$RGind_gene}})>0) {
			print STDERR $RGsubinfo, "Error: gene gets NO mRNAs: $RGind_gene\n";
			return $GffKit_failure;
		}
		my @RGtemparr=();
		foreach my $RGind_mrna (sort keys %{${$RGgene2mrna}{$RGind_gene}}) {
			unless (exists ${$RGmrnas}{$RGind_mrna}) {
				print STDERR $RGsubinfo, "Error: gene $RGind_gene gets NO mRNAs $RGind_mrna\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGmrnas}{$RGind_mrna}{'start'} and exists ${$RGmrnas}{$RGind_mrna}{'end'}) {
				print STDERR $RGsubinfo, "Error: invalid mRNAs $RGind_mrna start-end\n";
				return $GffKit_failure;
			}
			push (@RGtemparr, ${$RGmrnas}{$RGind_mrna}{'start'});
			push (@RGtemparr, ${$RGmrnas}{$RGind_mrna}{'end'});
		}
		@RGtemparr=sort {$a<=>$b} @RGtemparr;
		unless (exists ${$RGgene}{$RGind_gene}{'start'} and ${$RGgene}{$RGind_gene}{'start'}==$RGtemparr[0]) {
			if (${$RGgene}{$RGind_gene}{'start'}>$RGtemparr[0]) {
				print STDERR $RGsubinfo, "Warnings: invalid gene start: $RGind_gene ", ${$RGgene}{$RGind_gene}{'start'}, ' => ', $RGtemparr[0], "\n";
#				return $GffKit_failure;
			}
			else {
				print STDERR $RGsubinfo, "Warnings: Problematic gene start: $RGind_gene ", ${$RGgene}{$RGind_gene}{'start'}, ' => ', $RGtemparr[0], "\n";
			}
		}
		
		unless (exists ${$RGgene}{$RGind_gene}{'end'} and ${$RGgene}{$RGind_gene}{'end'}==$RGtemparr[-1]) {
			if (${$RGgene}{$RGind_gene}{'end'}<$RGtemparr[-1]) {
				print STDERR $RGsubinfo, "Warnings: invalid gene end: $RGind_gene ",${$RGgene}{$RGind_gene}{'end'}, ' => ', $RGtemparr[-1], "\n";
#				return $GffKit_failure;
			}
			else {
				print STDERR $RGsubinfo, "Warnings: Problematic gene end: $RGind_gene ",${$RGgene}{$RGind_gene}{'end'}, ' => ', $RGtemparr[-1], "\n";
			}
		}
		
		@RGtemparr=();
	}
	### %mrnas=($geneid => ('reference' => $arr[0],
###                     'start'      => $arr[3], 
###                     'end'        => $arr[4],
###                     'strand'     => $arr[6],
###                     'score'      => $arr[5],
###                     'Note'       => $note ###Not necessarily exist
###                    )
### 3. check if all the mRNA gets exons
	foreach my $RGind_mrna (sort keys %{$RGmrnas}) {
		unless (exists ${$RGexons}{$RGind_mrna}) {
			print STDERR $RGsubinfo, "Error: mRNA no exons: $RGind_mrna\n";
			return $GffKit_failure;
		}
	}

### 4.Check if all exons
	foreach my $RGind_mrna (sort keys %{$RGexons}) {
		unless (exists ${$RGmrnas}{$RGind_mrna}) {### if EXON mrna exists
			print STDERR $RGsubinfo, "Error: exon get no mRNA: $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (exists ${$RGexons}{$RGind_mrna}{'reference'}) {### if EXON reference exists
			print STDERR $RGsubinfo, "Error: exon get no reference: $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (exists ${$RGexons}{$RGind_mrna}{'strand'}) {### if EXON reference exists
			print STDERR $RGsubinfo, "Error: exon get no reference: $RGind_mrna\n";
			return $GffKit_failure;
		}
		### check with gene
		if (exists $RGmrna2gene{$RGind_mrna}) {
			my $RGindgene_id=$RGmrna2gene{$RGind_mrna};
			unless (exists ${$RGgene}{$RGindgene_id}) {### if EXON gene exists
				print STDERR $RGsubinfo, "Error: exon get no 'GENE': $RGind_mrna\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGgene}{$RGindgene_id}{'reference'} and (${$RGgene}{$RGindgene_id}{'reference'} eq ${$RGexons}{$RGind_mrna}{'reference'})) {### check ref if consistant with gene
				print STDERR $RGsubinfo, "Error: exon inconsistant ref: $RGind_mrna : GENE ", ${$RGgene}{$RGindgene_id}{'reference'}, " while EXON ", ${$RGexons}{$RGind_mrna}{'reference'}, "\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGgene}{$RGindgene_id}{'strand'} and (${$RGgene}{$RGindgene_id}{'strand'} eq ${$RGexons}{$RGind_mrna}{'strand'})) {### check strand if consistant with gene
				print STDERR $RGsubinfo, "Error: exon inconsistant strand: $RGind_mrna : GENE ", ${$RGgene}{$RGindgene_id}{'strand'}, " while EXON ", ${$RGexons}{$RGind_mrna}{'strand'}, "\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGexons}{$RGind_mrna}{'score'} and ${$RGexons}{$RGind_mrna}{'score'}=~/^\d+\.*\d*$/) {### check score if consistant with gene
				if (exists ${$RGgene}{$RGindgene_id}{'score'} and ${$RGgene}{$RGindgene_id}{'score'}=~/^\d+\.*\d*$/) {
					${$RGexons}{$RGind_mrna}{'score'}=${$RGgene}{$RGindgene_id}{'score'};
				}
				else {
					${$RGexons}{$RGind_mrna}{'score'}='.';
				}
			}
		}
		else {
			print STDERR $RGsubinfo, "Error: exon mRNA no GENE: $RGind_mrna\n";
			return $GffKit_failure;
		}
		### check with mRNA
		unless (exists ${$RGmrnas}{$RGind_mrna}{'reference'} and (${$RGmrnas}{$RGind_mrna}{'reference'} eq ${$RGexons}{$RGind_mrna}{'reference'})) {### check ref if consistant with mRNA
			print STDERR $RGsubinfo, "Error: exon inconsistant ref: mRNA $RGind_mrna ", ${$RGmrnas}{$RGind_mrna}{'reference'}, " while EXON ", ${$RGexons}{$RGind_mrna}{'reference'}, "\n";
				return $GffKit_failure;
		}
		unless (exists ${$RGmrnas}{$RGind_mrna}{'strand'} and (${$RGmrnas}{$RGind_mrna}{'strand'} eq ${$RGexons}{$RGind_mrna}{'strand'})) {### check ref if consistant with mRNA
			print STDERR $RGsubinfo, "Error: exon inconsistant strand: mRNA $RGind_mrna ", ${$RGmrnas}{$RGind_mrna}{'strand'}, " while EXON ", ${$RGexons}{$RGind_mrna}{'strand'}, "\n";
				return $GffKit_failure;
		}
		### check overlap
		my ($RGtest_temp, $RGexonarray)=&GetCoordArray(${$RGexons}{$RGind_mrna}{'exon'});
		unless ($RGtest_temp) {
			print STDERR $RGsubinfo, "Error: GetCoordArray running failed: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (scalar(@{$RGexonarray})==scalar(keys %{${$RGexons}{$RGind_mrna}{'exon'}})) {
			print STDERR $RGsubinfo, "Error: GetCoordArray number failed: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (&CheckFeatureOverlap($RGexonarray)) {
			print STDERR $RGsubinfo, "Error: exon overlap: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		### check if exon covered by mRNA
		unless (${$RGexonarray}[0][0] == ${$RGmrnas}{$RGind_mrna}{'start'}) {###
			print STDERR $RGsubinfo, "Warnings: mRNA left border: mRNA $RGind_mrna\n";
		}
		unless (${$RGexonarray}[-1][-1] == ${$RGmrnas}{$RGind_mrna}{'end'}) {###
			print STDERR $RGsubinfo, "Warnings: mRNA right border: mRNA $RGind_mrna\n";
		}
		$RGmrna_border{$RGind_mrna}=[${$RGexonarray}[0][0], ${$RGexonarray}[-1][-1]];
	}
### 5.Check if all CDS
	foreach my $RGind_mrna (sort keys %{$RGcds}) {
		unless (exists ${$RGmrnas}{$RGind_mrna}) {### if CDS mrna exists
			print STDERR $RGsubinfo, "Error: CDS get no mRNA: $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (exists ${$RGcds}{$RGind_mrna}{'reference'}) {### if CDS reference exists
			print STDERR $RGsubinfo, "Error: CDS get no reference: $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (exists ${$RGcds}{$RGind_mrna}{'strand'}) {### if CDS reference exists
			print STDERR $RGsubinfo, "Error: CDS get no reference: $RGind_mrna\n";
			return $GffKit_failure;
		}
		### check with gene
		if (exists $RGmrna2gene{$RGind_mrna}) {
			my $RGindgene_id=$RGmrna2gene{$RGind_mrna};
			unless (exists ${$RGgene}{$RGindgene_id}) {### if CDS gene exists
				print STDERR $RGsubinfo, "Error: CDS get no 'GENE': $RGind_mrna\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGgene}{$RGindgene_id}{'reference'} and (${$RGgene}{$RGindgene_id}{'reference'} eq ${$RGcds}{$RGind_mrna}{'reference'})) {### check ref if consistant with gene
				print STDERR $RGsubinfo, "Error: CDS inconsistant ref: $RGind_mrna : GENE ", ${$RGgene}{$RGindgene_id}{'reference'}, " while CDS ", ${$RGcds}{$RGind_mrna}{'reference'}, "\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGgene}{$RGindgene_id}{'strand'} and (${$RGgene}{$RGindgene_id}{'strand'} eq ${$RGcds}{$RGind_mrna}{'strand'})) {### check strand if consistant with gene
				print STDERR $RGsubinfo, "Error: CDS inconsistant strand: $RGind_mrna : GENE ", ${$RGgene}{$RGindgene_id}{'strand'}, " while CDS ", ${$RGcds}{$RGind_mrna}{'strand'}, "\n";
				return $GffKit_failure;
			}
			unless (exists ${$RGcds}{$RGind_mrna}{'score'} and ${$RGcds}{$RGind_mrna}{'score'}=~/^\d+\.*\d*$/) {### check score if consistant with gene
				if (exists ${$RGgene}{$RGindgene_id}{'score'} and ${$RGgene}{$RGindgene_id}{'score'}=~/^\d+\.*\d*$/) {
					${$RGcds}{$RGind_mrna}{'score'}=${$RGgene}{$RGindgene_id}{'score'};
				}
				else {
					${$RGcds}{$RGind_mrna}{'score'}='.';
				}
			}
		}
		else {
			print STDERR $RGsubinfo, "Error: CDS mRNA no GENE: $RGind_mrna\n";
			return $GffKit_failure;
		}
		### check overlap
		my ($RGtest_temp1, $RGcdsarray)=&GetCoordArray(${$RGcds}{$RGind_mrna}{'cds'});
		unless ($RGtest_temp1) {
			print STDERR $RGsubinfo, "Error: GetCoordArray CDS running failed: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (scalar(@{$RGcdsarray})==scalar(keys %{${$RGcds}{$RGind_mrna}{'cds'}})) {
			print STDERR $RGsubinfo, "Error: GetCoordArray number failed: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (&CheckFeatureOverlap($RGcdsarray)) {### check if CDS feature overlap
			print STDERR $RGsubinfo, "Error: exon overlap: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		my ($RGtest_temp2, $RGexonarray)=&GetCoordArray(${$RGexons}{$RGind_mrna}{'exon'});
		unless ($RGtest_temp2) {
			print STDERR $RGsubinfo, "Error: GetCoordArray exon2 running failed: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		unless (&CdsCoordConfirm($RGexonarray, $RGcdsarray)) {### Confirm all CDS are in exons
			print STDERR $RGsubinfo, "Error: CDS out of exons: mRNA $RGind_mrna\n";
			return $GffKit_failure;
		}
		my $RGexistsingphase=0;
		my $RGnum_notphased=0;
		my $RGphases=[];
#		print $RGsubinfo, "Test: mRNA $RGind_mrna \$RGcdsarray\n"; print Dumper $RGcdsarray;### For test ###
#		print $RGsubinfo, "Test: mRNA $RGind_mrna \${\$RGcds}{\$RGind_mrna}{'cds'}\n"; print Dumper ${$RGcds}{$RGind_mrna}{'cds'};### For test ###
		foreach my $RGcdsstart (@{$RGcdsarray}) {
			if (exists ${$RGcds}{$RGind_mrna}{'phase'} and exists ${$RGcds}{$RGind_mrna}{'phase'}{$RGcdsstart->[0]} and exists ${$RGcds}{$RGind_mrna}{'phase'}{$RGcdsstart->[0]}{$RGcdsstart->[1]}) {
				if (${$RGcds}{$RGind_mrna}{'phase'}{$RGcdsstart->[0]}{$RGcdsstart->[1]}=~/^[012]{1,1}$/) {
					push (@{$RGphases}, ${$RGcds}{$RGind_mrna}{'phase'}{$RGcdsstart->[0]}{$RGcdsstart->[1]});
					$RGexistsingphase++;
				}
				else {
					push (@{$RGphases}, '.');
					$RGnum_notphased++;
				}
			}
			else {
				print STDERR $RGsubinfo, "Error: mRNA $RGind_mrna invalid CDS phase\n";
				return $GffKit_failure;
			}
		}
		
#		print $RGsubinfo, "Test: mRNA $RGind_mrna \$RGphases 1\n"; print Dumper $RGphases;### For test ###
#		print $RGsubinfo, "Test: mRNA $RGind_mrna phased: ", $RGexistsingphase, " Not phased ", $RGnum_notphased, "\n"; ### For test ###
		if ($RGnum_notphased>0 and $RGexistsingphase>0 and ($RGexistsingphase<scalar(keys %{${$RGcds}{$RGind_mrna}{'phase'}}))) {
#			print $RGsubinfo, "Test: mRNA $RGind_mrna Partially phased\n";### For test ###
#			print Dumper $RGphases;### For test ###
			unless (&GetPartialPhase($RGcdsarray, $RGphases, ${$RGcds}{$RGind_mrna}{'strand'})) {
				print STDERR $RGsubinfo, "Warnings: Error to auto complete phase for mRNA: $RGind_mrna\n";
			}
			else {
#				print $RGsubinfo, "Test: mRNA $RGind_mrna after autocomplete\n";### For test ###
#				print Dumper $RGphases; ### For test ###
			}
		}
#		elsif ($RGnum_notphased==0 and ($RGexistsingphase==scalar(keys %{${$RGcds}{$RGind_mrna}{'phase'}}))) {### For test ###
#			print $RGsubinfo, "Test: mRNA $RGind_mrna Fully phased\n";
#		}
#		elsif (($RGnum_notphased==scalar(keys %{${$RGcds}{$RGind_mrna}{'phase'}})) and $RGexistsingphase==0) {### For test ###
#			print $RGsubinfo, "Test: mRNA $RGind_mrna NOT phased\n";
#		}
#		print $RGsubinfo, "Test: mRNA $RGind_mrna \$RGphases 2\n"; print Dumper $RGphases;### For test ###
		RGLOOP2: {if ($RGcheck_phase) {### use genome sequence to correct phases
#			print $RGsubinfo, "Test: Get phase from fasta\n";### For test ###
			my $RGcdsseq='';
			my $RGlength=0;
			my $RGprot={};
			foreach my $RGindcds (@{$RGcdsarray}) {
				$RGcdsseq.=$RGseqdb->seq(${$RGcds}{$RGind_mrna}{'reference'}, $RGindcds->[0] => $RGindcds->[1]);
				$RGlength=$RGlength+($RGindcds->[1]-$RGindcds->[0]+1);
			}
			unless ($RGlength=~/^\d+$/ and $RGlength>0 and length($RGcdsseq) == $RGlength) {
				print STDERR "Warnings: unable to get CDS sequence: mRNA ", $RGind_mrna, " length ", $RGlength, " REF ", ${$RGcds}{$RGind_mrna}{'reference'}, "\n";
				last RGLOOP2;
			}
#			print $RGsubinfo, "Test: strand", ${$RGcds}{$RGind_mrna}{'strand'}, "\n";### For test ###
			if (${$RGcds}{$RGind_mrna}{'strand'} eq '+') {
				$RGprot=Frame3Translation($RGcdsseq);
				my @RGbestframe=();
				foreach (my $RGframe=0; $RGframe<3; $RGframe++) {
					if (exists ${$RGprot}{$RGframe} and (${$RGprot}{$RGframe} =~/\S+/)) {
						unless (${$RGprot}{$RGframe} =~ /\S+\*\S+/) {
							push (@RGbestframe, $RGframe);
						}
						else {
							if ($GffKit_debug) {
								print $RGsubinfo, "Test: Frame $RGframe internal stop: "; print ${$RGprot}{$RGframe}, "\n";
							}
						}
					}
				}
#				print $RGsubinfo, "Test: 111 @RGbestframe\n";### For test ###
				if (scalar(@RGbestframe)>1) {
					my @RGtempframe=();
					my $RGknownframse=0;
#					print $RGsubinfo, "Test: 444\n";### For test ###
					RGLOOP3: {foreach my $RGframe (@RGbestframe) {
						if (defined ${$RGphases}[0] and ${$RGphases}[0]=~/^[012]{1,1}$/) {
							if ($RGframe==${$RGphases}[0]) {
								@RGbestframe=(); push (@RGbestframe, $RGframe);
								$RGknownframse=1;
								last RGLOOP3;
							}
						}
						else {
							if (${$RGprot}{$RGframe}=~/\S+\*$/) {
								push (@RGtempframe, $RGframe);
							}
							elsif (${$RGprot}{$RGframe}=~/^M/) {
								push (@RGtempframe, $RGframe);
							}
						}
					}}###RGLOOP3
					unless ($RGknownframse) {
						@RGbestframe=@RGtempframe;
						@RGtempframe=();
					}
				}
				elsif (scalar(@RGbestframe)==0) {
					print $RGsubinfo, "Error: strand + internal stop: mRNA $RGind_mrna\n";
					print Dumper $RGphases;
					print Dumper $RGprot;
#					return $GffKit_failure;
				}
				
				if (scalar(@RGbestframe)==1) {
#					print $RGsubinfo, "Test: 333\n";### For test ###
					$RGphases=[];
#					print $RGsubinfo, "Test: mRNA $RGind_mrna Get CDS phase +\n";### For Test ###
					my $z=shift @RGbestframe;
					foreach my $RGy (@{$RGcdsarray}) {
						push (@{$RGphases}, $z);
						$z=3 - ($RGy->[1] - $RGy->[0] +1 - $z)%3;
						$z=0 if ($z==3);
					}
				}
			}
			elsif (${$RGcds}{$RGind_mrna}{'strand'} eq '-') {
				$RGcdsseq=SeqRevComp($RGcdsseq);
				$RGprot=Frame3Translation($RGcdsseq);
				my @RGbestframe=();
				foreach (my $RGframe=0; $RGframe<3; $RGframe++) {
					if (exists ${$RGprot}{$RGframe} and (${$RGprot}{$RGframe} =~/\S+/)) {
						unless (${$RGprot}{$RGframe} =~ /\S+\*\S+/) {
							push (@RGbestframe, $RGframe);
						}
						else {
							if ($GffKit_debug) {
#								print $RGsubinfo, "Test: Frame $RGframe internal stop: "; print ${$RGprot}{$RGframe}, "\n";### For test ###
							}
						}
					}
				}
				if (scalar(@RGbestframe)>1) {
					my @RGtempframe=();
					my $RGknownframse=0;
					RGLOOP4: {foreach my $RGframe (@RGbestframe) {
						if (defined ${$RGphases}[-1] and ${$RGphases}[-1]=~/^[012]{1,1}$/) {
							if ($RGframe==${$RGphases}[-1]) {
								@RGbestframe=(); push (@RGbestframe, $RGframe);
								$RGknownframse=1;
								last RGLOOP4;
							}
						}
						else {
							if (${$RGprot}{$RGframe}=~/\S+\*$/) {
								push (@RGtempframe, $RGframe);
							}
							elsif (${$RGprot}{$RGframe}=~/^M/) {
								push (@RGtempframe, $RGframe);
							}
						}
					}}###RGLOOP4
					unless ($RGknownframse) {
						@RGbestframe=@RGtempframe;
						@RGtempframe=();
					}
				}
				elsif (scalar(@RGbestframe)==0) {
					print $RGsubinfo, "Error: strand - internal stop: mRNA  $RGind_mrna\n";
					print Dumper $RGphases;
					print Dumper $RGprot;
#					return $GffKit_failure;
				}
				if (scalar(@RGbestframe)==1) {
					$RGphases=[];
#					print $RGsubinfo, "Test: mRNA $RGind_mrna Get CDS phase -\n";### For Test ###
					my $z=shift @RGbestframe;
					for (my $RGx=scalar(@{$RGcdsarray})-1; $RGx>=0; $RGx--) {
						unshift (@{$RGphases}, $z);
						my $RGy=${$RGcdsarray}[$RGx];
						$z=3 - ($RGy->[1] - $RGy->[0] +1 - $z)%3;
						$z=0 if ($z==3);
					}
				}
			}
			else {
				print STDERR "Warnings: invalid strand mRNA ", $RGind_mrna, "\n";
				last RGLOOP2;
			}
			$RGcdsseq='';
		}}###RGLOOP2
#		print $RGsubinfo, "Test: mRNA $RGind_mrna \$RGphases 3\n"; print Dumper $RGphases;### For test ###
		unless (scalar(@{$RGphases})==scalar(keys %{${$RGcds}{$RGind_mrna}{'phase'}})) {
			print STDERR "Error: invalid phases ", $RGind_mrna, "\n";
			print Dumper $RGphases;
			return $GffKit_failure;
		}
		else {
			### Check if all phased
			for (my $RGi=0; $RGi<scalar(@{$RGphases}); $RGi++) {
				unless (defined ${$RGphases}[$RGi] and ${$RGphases}[$RGi]=~/^[012]{1,1}$/) {
					print STDERR "Error: can not determine phased value ", $RGind_mrna, "\n";
					print Dumper $RGphases;
					return $GffKit_failure;
				}
			}
			### Transfer nre phases
			for (my $RGi=0; $RGi<scalar(@{$RGcdsarray}); $RGi++) {
				my ($RGx, $RGy)=@{${$RGcdsarray}[$RGi]};
				unless (exists ${$RGcds}{$RGind_mrna}{'phase'}{$RGx} and exists ${$RGcds}{$RGind_mrna}{'phase'}{$RGx}{$RGy}) {
					print STDERR "Error: CDS not exists ", $RGind_mrna, " $RGx - $RGy\n";
					print Dumper ${$RGcds}{$RGind_mrna}{'phase'};
					return $GffKit_failure;
				}
				${$RGcds}{$RGind_mrna}{'phase'}{$RGx}{$RGy}=${$RGphases}[$RGi];
			}
		}
	}
	
	if ($GffKit_debug) {### For test ###
#		print "Test: \$RGRGreferenceids\n"; print Dumper $RGreferenceids; print "\n"; ### For test ###
#		print "Test: \$RGgene\n"; print Dumper $RGgene; print "\n"; ### For test ###
#		print "Test: \$RGgene2mrna\n"; print Dumper $RGgene2mrna; print "\n"; ### For test ###
#		print "Test: \$RGmrnas\n"; print Dumper $RGmrnas; print "\n"; ### For test ###
#		print "Test: \$RGexons\n"; print Dumper $RGexons; print "\n"; ### For test ###
#		print "Test: \$RGcds\n"; print Dumper $RGcds; print "\n"; ### For test ###
	}
	%RGmrna2gene=();
	%RGmrna_border=();
	return ($GffKit_success, $RGreferenceids, $RGgene, $RGgene2mrna, $RGmrnas, $RGexons, $RGcds);
}
### Convert exons coordinates hash into array, and check if exonID unique
### &GetCoordArray ($hash)
### Input:  %{$hash}={$start1 => $end1 => exonid1, $start2 => $end2 => exonid2, ...}
### Output: $arr=[[$start1, $end1], [$start2, $end2]...]
###            array is '$start'-sorted (small -> large: $start1<=$start2) for overlap checking
sub GetCoordArray {
	my $GCAhash=shift;
	
	my $GCAret_array=[];
	my $GCAsubinfo='SUB(GffKit::GetCoordArray)';
	my %GCAexonid=();
	
	foreach my $GCAstart (sort {$a<=>$b} keys %{$GCAhash}) {
		my @GCAtemp_arr=();
		@GCAtemp_arr=keys %{${$GCAhash}{$GCAstart}};
		unless ($GCAstart=~/^\d+$/ and scalar(@GCAtemp_arr)==1 and $GCAtemp_arr[0]=~/^\d+$/) {
			print STDERR $GCAsubinfo, "Error: invalid coordinates\n";
			print Dumper $GCAhash;
			return $GffKit_failure;
		}
		push (@{$GCAret_array}, [$GCAstart, $GCAtemp_arr[0]]);
		$GCAexonid{${$GCAhash}{$GCAstart}{$GCAtemp_arr[0]}}++;
	}
	foreach my $GCAexonname (sort keys %GCAexonid) {### Check unique exonID
		unless ($GCAexonid{$GCAexonname}==1 or $GCAexonid{$GCAexonname}==scalar(keys %{$GCAhash})) {
			print STDERR $GCAsubinfo, "Error: exon/CDS ID not unique: $GCAexonname\n";
			return $GffKit_failure;
		}
	}
	
	return ($GffKit_success, $GCAret_array);
}
### Check if coordinates get overlap
### &CheckFeatureOverlap($arr)
### Input: $arr=[[$start1, $end1], [$start2, $end2]...]
### Output: 0=overlap; 1=NO_overlap
sub CheckFeatureOverlap {
	my $CFOtest_arr=shift;
	
	my $CROsubinfo='SUB(GffKit::CheckFeatureOverlap)';
	
	for (my $CFOi=0; $CFOi<scalar(@{$CFOtest_arr}); $CFOi++) {
		unless (scalar(@{${$CFOtest_arr}[$CFOi]})==2 and defined ${$CFOtest_arr}[$CFOi][0] and ${$CFOtest_arr}[$CFOi][0]=~/^\d+$/ and defined ${$CFOtest_arr}[$CFOi][1] and ${$CFOtest_arr}[$CFOi][1]=~/^\d+$/) {
			print STDERR $CROsubinfo, "Error: invalid feature\n";
			print Dumper $CFOtest_arr;
			return $GffKit_failure;
		}
		if ($CFOi>1) {
			unless (${$CFOtest_arr}[$CFOi][0] > ${$CFOtest_arr}[$CFOi-1][1]) {
				print STDERR $CROsubinfo, "Error: feature overlap: ", ${$CFOtest_arr}[$CFOi-1][0], '-', ${$CFOtest_arr}[$CFOi-1][1], ' and ', ${$CFOtest_arr}[$CFOi][0], '-', ${$CFOtest_arr}[$CFOi][1], "\n";
				return $GffKit_failure;
			}
		}
	}
	return $GffKit_success;
}
### check if CDS totally included in exons
### CdsCoordConfirm ($exon_arr, $cds_arr)
### arr is return by GetCoordArray
### Return: 1=perfect, 0=problematic
sub CdsCoordConfirm {
	my ($CCCexons, $CCCcds)=@_;
	
	my $CCCsubinfo='SUB(GffKit::CdsCoordConfirm)';
	my $CCCcdsmin;
	my $CCCcdsmax;
	my $CCCorfstart=0;
	my $CCCorfend=0;
	my @CCCexpectedcds=();
	
	
	unless (scalar(@{$CCCcds})>0) {
		print STDERR $CCCsubinfo, "Error: no CDS coordinates\n";
		return $GffKit_failure;
	}
	unless (defined ${$CCCcds}[0] and defined ${$CCCcds}[0][0] and ${$CCCcds}[0][0]=~/^\d+$/) {
		print STDERR $CCCsubinfo, "Error: CDS array problem1\n";
		print Dumper $CCCcds;
		return $GffKit_failure;
	}
	$CCCcdsmin=${$CCCcds}[0][0];
	unless (defined ${$CCCcds}[-1] and defined ${$CCCcds}[-1][-1] and ${$CCCcds}[-1][-1]=~/^\d+$/) {
		print STDERR $CCCsubinfo, "Error: CDS array problem2\n";
		print Dumper $CCCcds;
		return $GffKit_failure;
	}
	$CCCcdsmax=${$CCCcds}[-1][-1];
	CCCLOOP1 : {foreach my $CCexon (@{$CCCexons}) {
		if ($CCCorfstart==0) {
			if ($CCCcdsmin>=$CCexon->[0] and $CCCcdsmin<=$CCexon->[1]) {
				if ($CCCcdsmax<=$CCexon->[1]) {
					push (@CCCexpectedcds, [$CCCcdsmin, $CCCcdsmax]);
					$CCCorfstart=1;
					$CCCorfend=1;
				}
				else {
					push (@CCCexpectedcds, [$CCCcdsmin, $CCexon->[1]]);
					$CCCorfstart=1;
				}
			}
			elsif ($CCCcdsmin < $CCexon->[0]) {
				print STDERR $CCCsubinfo, "Error: CDS coordinate problem 1: might missing exons: $CCCcdsmin vs ", $CCexon->[0], "\n";
#				print $CCCsubinfo, "Test: $CCCcdsmin\n"; ### for Test ###
				return $GffKit_failure if ($CCCcdsmin==0);
			}
		}
		elsif ($CCCorfstart==1 and $CCCorfend==0) {
			if ($CCCcdsmax>$CCexon->[1]) {
				push (@CCCexpectedcds, [$CCexon->[0], $CCexon->[1]]);
			}
			elsif ($CCCcdsmax>=$CCexon->[0] and $CCCcdsmax<=$CCexon->[1]) {
				push (@CCCexpectedcds, [$CCexon->[0], $CCCcdsmax]);
				$CCCorfend=1;
			}
		}
		elsif ($CCCorfend==1) {
			last CCCLOOP1;
		}
	}}###CCCLOOP1
	unless ($CCCorfstart==1 and $CCCorfend==1) {
		print STDERR $CCCsubinfo, "Error: CDS coordinate problem 2: start / end\n";
		return $GffKit_failure;
	}
	unless (scalar(@CCCexpectedcds)==scalar(@{$CCCcds})) {
		print STDERR $CCCsubinfo, "Error: CDS coordinate problem 3\n";
		print $CCCsubinfo, "Test: \@CCCexpectedcds\n"; print Dumper \@CCCexpectedcds; ### for Test ###
		print $CCCsubinfo, "Test: \$CCCcds\n"; print Dumper $CCCcds; ### for Test ###
		return $GffKit_failure;
	}
	for (my $CCCi=0; $CCCi<scalar(@CCCexpectedcds); $CCCi++) {
		unless (($CCCexpectedcds[$CCCi][0] == ${$CCCcds}[$CCCi][0]) and ($CCCexpectedcds[$CCCi][1] == ${$CCCcds}[$CCCi][1])) {
			print STDERR $CCCsubinfo, "Error: CDS coordinate problem 4: CDS ", ${$CCCcds}[$CCCi][0], '-', ${$CCCcds}[$CCCi][1],  ' expected: ', $CCCexpectedcds[$CCCi][0], '-', $CCCexpectedcds[$CCCi][1], "\n";
			return $GffKit_failure;
		}
	}
	return $GffKit_success;
}
###	Calculate phase based on partial phase
### Input: $arr=[[$start1, $end1], [$start2, $end2]...]
sub GetPartialPhase {
	my ($GPPcdsarr, $GCPphases, $GPPstrand)=@_;
	
	my $GPPsubinfo='SUB(GffKit::GetPartialPhase)';
	my $GPPret_phase_arr=[];
	my $GPPfit=0;
	my $GPPunfit=0;
	
	unless (scalar(@{$GPPcdsarr})==scalar(@{$GCPphases})) {
		print STDERR $GPPsubinfo, "Error: phase num != CDS num\n";
		return $GffKit_failure;
	}
	
	if ($GPPstrand eq '+') {
		GPPLOOP1: { for (my $GPPi=0; $GPPi<3; $GPPi++) {
			my $GPPphase=$GPPi;
			$GPPret_phase_arr=[];
			$GPPfit=0;
			$GPPunfit=0;
			GPPLOOP2 : {for (my $GPPj=0; $GPPj<scalar(@{$GPPcdsarr}); $GPPj++) {
				if (${$GCPphases}[$GPPj]=~/^[012]{1,1}$/) {
					if (${$GCPphases}[$GPPj]==$GPPphase) {
						$GPPfit=1;
					}
					else {
						$GPPunfit=1;
						next GPPLOOP1;
					}
				}
				push (@{$GPPret_phase_arr}, $GPPphase);
				$GPPphase=3-((${$GPPcdsarr}[$GPPj][1]-${$GPPcdsarr}[$GPPj][0]+1-$GPPphase)%3);
				$GPPphase=0 if ($GPPphase==3);
			}}###GPPLOOP2
			if ($GPPfit==1 and $GPPunfit==0) {
				last GPPLOOP1;
			}
		}}###GPPLOOP1
	}
	elsif ($GPPstrand eq '-') {
		GPPLOOP3: { for (my $GPPi=0; $GPPi<3; $GPPi++) {
			my $GPPphase=$GPPi;
			$GPPret_phase_arr=[];
			$GPPfit=0;
			$GPPunfit=0;
			GPPLOOP4 : { for (my $GPPj=0; $GPPj<scalar(@{$GPPcdsarr}); $GPPj++) {
				if (${$GCPphases}[$GPPj]=~/^[012]{1,1}$/) {
					if (${$GCPphases}[$GPPj]==$GPPphase) {
						$GPPfit=1;
					}
					else {
						$GPPunfit=1;
						next GPPLOOP3;
					}
				}
				$GPPphase=(${$GPPcdsarr}[$GPPj][1]-${$GPPcdsarr}[$GPPj][0]+1-$GPPphase)%3;
				$GPPphase=0 if ($GPPphase==3);
				push (@{$GPPret_phase_arr}, $GPPphase);
			}}###GPPLOOP4
			if ($GPPfit==1 and $GPPunfit==0) {
				last GPPLOOP3;
			}
		}}###GPPLOOP3
	}
	else {
		print STDERR $GPPsubinfo, "Error: invalid strand\n";
		return $GffKit_failure;
	}
	unless (scalar(@{$GCPphases})==scalar(@{$GPPret_phase_arr})) {
		print STDERR $GPPsubinfo, "Error: phase complete error\n";
		return $GffKit_failure;
	}
	
	@{$GCPphases}=@{$GPPret_phase_arr};
	
	return $GffKit_success;
}




### Write GFF3 output
### WriteGff3($EGoutgff3, $WGref2gene, $WGgene2mrna, $WGgene, $WGmrna, $WGexon, $WGcds)
###
### %referenceids  => ( $reference_id => $gene_start_pos => $gene_id => num++)
### %gene2mrna     => ( $gene_id => $mrna_id => num++ )
### %gene=($geneid => ('reference' => $arr[0],
###                    'start'     => $arr[3], 
###                    'end'       => $arr[4],
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5],
###                    'note'      => $note ###Not necessarily exist
###                       )
###       )
### %mrnas=($mrnaid => ('reference' => $arr[0],
###                    'start'     => $arr[3], 
###                    'end'       => $arr[4],
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5],
###                    'note'      => $note ###Not necessarily exist
###       )
### %exons=($mrnaid => ('reference' => $arr[0],
###                    'exon' => ({$arr[3]} => ($arr[4] => $exonid)),
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5],
###                    )
###       )
### %cds=($mrnaid => ('reference' => $arr[0],
###                   'cds' => ({$arr[3]} => ($arr[4] => num++)),
###                   'strand'    => $arr[6],
###                   'score'     => $arr[5],
###                   'phase'     => ({$arr[3]} => ($arr[4] => $arr[7]))
###                   )
###       )

### Return: 1=success, 0=failed
sub WriteGff3 {
	my ($EGoutgff3, $WGref2gene, $WGgene2mrna, $WGgene, $WGmrna, $WGexon, $WGcds) =@_;
	
	my $WGsubinfo='SUB(GffKit::WriteGff3)';
	local *WGGFF3OUT;
	
	close WGGFF3OUT if (defined fileno(WGGFF3OUT));
	unless (open (WGGFF3OUT, " > $EGoutgff3 ")) {
		print STDERR $WGsubinfo, "Error: can not write GFF3 out: $EGoutgff3\n";
		return $GffKit_failure;
	}
	print WGGFF3OUT '##gff-version 3', "\n";
	foreach my $WGref (sort keys %{$WGref2gene}) {
		foreach my $WGgenepos (sort {$a<=>$b} keys %{${$WGref2gene}{$WGref}}) {
			foreach my $WGgeneid (sort keys %{${$WGref2gene}{$WGref}{$WGgenepos}}) {
				unless (exists ${$WGgene}{$WGgeneid}) {
					print STDERR $WGsubinfo, "Error: gene not exists: GENE $WGgeneid\n";
					return $GffKit_failure;
				}
				unless (exists ${$WGgene}{$WGgeneid}{'reference'} and 
						( ${$WGgene}{$WGgeneid}{'reference'} eq $WGref )) {
					print STDERR $WGsubinfo, "Error: inconsistent ref $WGref: GENE $WGgeneid\n";
					return $GffKit_failure;
				}
				unless (exists ${$WGgene}{$WGgeneid}{'start'} and 
						${$WGgene}{$WGgeneid}{'start'} =~/^\d+$/) {
					print STDERR $WGsubinfo, "Error: invalid gene start: GENE $WGgeneid\n";
					return $GffKit_failure;
				}
				unless (exists ${$WGgene}{$WGgeneid}{'end'} and 
						${$WGgene}{$WGgeneid}{'end'} =~/^\d+$/) {
					print STDERR $WGsubinfo, "Error: invalid gene end: GENE $WGgeneid\n";
					return $GffKit_failure;
				}
				unless (exists ${$WGgene}{$WGgeneid}{'strand'} and 
						${$WGgene}{$WGgeneid}{'strand'} =~/^(\+)|(\-)$/) {
					print STDERR $WGsubinfo, "Error: invalid gene strand: GENE $WGgeneid\n";
					return $GffKit_failure;
				}
				unless (exists ${$WGgene}{$WGgeneid}{'score'} and ${$WGgene}{$WGgeneid}{'score'} =~/^\d+\.*\d*$/) {
						${$WGgene}{$WGgeneid}{'score'}='.';
				}
				my @WGarr=();
				@WGarr=($WGref, '.', 'gene', ${$WGgene}{$WGgeneid}{'start'}, ${$WGgene}{$WGgeneid}{'end'}, ${$WGgene}{$WGgeneid}{'score'}, ${$WGgene}{$WGgeneid}{'strand'}, '.', "ID=$WGgeneid;Name=$WGgeneid");
				foreach my $WGkey (keys %{${$WGgene}{$WGgeneid}}) {
#seqid	source	type	start	end	score	strand	phase	attributes
##gff-version 3
#ctg123 . mRNA            1300  9000  .  +  .  ID=mrna0001;Name=sonichedgehog
					next if ($WGkey =~/^(reference)|(start)|(end)|(strand)|(score)$/);
					if ($WGkey eq 'score' and ${$WGgene}{$WGgeneid}{'score'} =~/^\d+\.*\d*$/) {
						$WGarr[5]=${$WGgene}{$WGgeneid}{'score'};
					}
					elsif ($WGkey eq 'note' and ${$WGgene}{$WGgeneid}{'note'} =~/^\S+/) {
						$WGarr[8]=$WGarr[8].';Note='.${$WGgene}{$WGgeneid}{'note'};
					}
					else {
						$WGarr[8]=$WGarr[8].';'.$WGkey.'='.${$WGgene}{$WGgeneid}{$WGkey};
					}
				}
				print WGGFF3OUT "##sequence-region   ", $WGarr[0], " ", $WGarr[3], " ", $WGarr[4], "\n";
				print WGGFF3OUT join ("\t", @WGarr), "\n";
				### write mRNA
				@WGarr=();
				if (exists ${$WGgene2mrna}{$WGgeneid} and scalar(keys %{${$WGgene2mrna}{$WGgeneid}})>0) {
					my %WGexonidhash=();
					foreach my $WGmrnaid (sort keys %{${$WGgene2mrna}{$WGgeneid}}) {
						unless (exists ${$WGmrna}{$WGmrnaid}) {
							print STDERR "Warnings; mRNA no data: $WGmrnaid\n";
							next;
						}
						unless (exists ${$WGmrna}{$WGmrnaid}{'reference'} and 
								( ${$WGmrna}{$WGmrnaid}{'reference'} eq $WGref )) {
							print STDERR $WGsubinfo, "Error: inconsistent ref: mRNA $WGmrnaid\n";
							return $GffKit_failure;
						}
						unless (exists ${$WGmrna}{$WGmrnaid}{'start'} and 
								${$WGmrna}{$WGmrnaid}{'start'} =~/^\d+$/) {
							print STDERR $WGsubinfo, "Error: invalid mRNA start: mRNA $WGmrnaid\n";
							return $GffKit_failure;
						}
						unless (exists ${$WGmrna}{$WGmrnaid}{'end'} and 
								${$WGmrna}{$WGmrnaid}{'end'} =~/^\d+$/) {
							print STDERR $WGsubinfo, "Error: invalid mRNA end: mRNA $WGmrnaid\n";
							return $GffKit_failure;
						}
						unless (exists ${$WGmrna}{$WGmrnaid}{'strand'} and 
								${$WGmrna}{$WGmrnaid}{'strand'} =~/^(\+)|(\-)$/ and 
								(${$WGmrna}{$WGmrnaid}{'strand'} eq ${$WGgene}{$WGgeneid}{'strand'}) ) {
							print STDERR $WGsubinfo, "Error: inconsistent strand: REF $WGref GENE $WGgeneid ".${$WGgene}{$WGgeneid}{'strand'}." mRNA ".${$WGmrna}{$WGmrnaid}{'strand'}."\n";
							return $GffKit_failure;
						}
						unless (exists ${$WGmrna}{$WGmrnaid}{'score'} and ${$WGmrna}{$WGmrnaid}{'score'} =~/^(\+)|(\-)$/) {
							if (exists ${$WGgene}{$WGgeneid}{'score'} and ${$WGgene}{$WGgeneid}{'score'} =~/^\d+\.*\d*$/) {
								${$WGmrna}{$WGmrnaid}{'score'}=${$WGgene}{$WGgeneid}{'score'};
							}
							else {
								${$WGmrna}{$WGmrnaid}{'score'}='.';
							}
						}
						
						@WGarr=();
						@WGarr=($WGref, '.', 'mRNA', ${$WGmrna}{$WGmrnaid}{'start'}, ${$WGmrna}{$WGmrnaid}{'end'}, ${$WGmrna}{$WGmrnaid}{'score'}, ${$WGmrna}{$WGmrnaid}{'strand'}, '.', "ID=$WGmrnaid;Parent=$WGgeneid");
						foreach my $WGkey2 (keys %{${$WGmrna}{$WGmrnaid}}) {
							next if ($WGkey2 =~/^(reference)|(start)|(end)|(strand)|(score)$/);
#							print $WGsubinfo, "Test: gene $WGgeneid \${\$WGgene}{\$WGgeneid}\n"; print Dumper ${$WGgene}{$WGgeneid}; ### For test ###
							if (($WGkey2 eq 'score') and ${$WGmrna}{$WGmrnaid}{$WGkey2} =~/^\d+\.*\d*$/) {
								$WGarr[5]=${$WGmrna}{$WGmrnaid}{$WGkey2};
							}
							elsif (($WGkey2 eq 'note') and ${$WGmrna}{$WGmrnaid}{$WGkey2} =~/^\S+$/) {
								$WGarr[8]=$WGarr[8].';Note='.${$WGmrna}{$WGmrnaid}{$WGkey2};
							}
							else {
								$WGarr[8]=$WGarr[8].';'.$WGkey2.'='.${$WGmrna}{$WGmrnaid}{$WGkey2};
							}
						}
						print WGGFF3OUT join ("\t", @WGarr), "\n";
					
						### write exons
						if (exists ${$WGexon}{$WGmrnaid}) {
							unless (exists ${$WGexon}{$WGmrnaid}{'reference'} and 
									(${$WGexon}{$WGmrnaid}{'reference'} eq $WGref)
							) {
								print STDERR $WGsubinfo, "Error: inconsistent exon ref ($WGref): mRNA $WGmrnaid REF ", ${$WGexon}{$WGmrnaid}{'reference'}, "\n";
								
								return $GffKit_failure;
							}
							unless (exists ${$WGexon}{$WGmrnaid}{'strand'} and 
									(${$WGexon}{$WGmrnaid}{'strand'} eq ${$WGgene}{$WGgeneid}{'strand'})
							) {
								print STDERR $WGsubinfo, "Error: inconsistent exon strand: mRNA $WGmrnaid\n";
								return $GffKit_failure;
							}
							my $WGexonscore='.';
							if (exists ${$WGexon}{$WGmrnaid}{'score'} and ${$WGexon}{$WGmrnaid}{'score'} =~/^\d+$/) {
								$WGexonscore=${$WGexon}{$WGmrnaid}{'score'};
							}
							if (exists ${$WGexon}{$WGmrnaid}{'exon'} and 
								(scalar(keys %{${$WGexon}{$WGmrnaid}{'exon'}})>0)
							){
								foreach my $WGx (keys %{${$WGexon}{$WGmrnaid}{'exon'}}) {
									foreach my $WGy (keys %{${$WGexon}{$WGmrnaid}{'exon'}{$WGx}}) {
										my $WGz=${$WGexon}{$WGmrnaid}{'exon'}{$WGx}{$WGy};
										if (exists $WGexonidhash{$WGz}) {###exonid must be unique
											if (($WGexonidhash{$WGz}[0] != $WGx) and ($WGexonidhash{$WGz}[0] != $WGy)) {
												print STDERR $WGsubinfo, "Error: repeated exonID: $WGz for mRNA $WGmrnaid\n";
												return $GffKit_failure;
											}
											else {
												next;
											}
										}
										unless ($WGy=~/^\d+$/ and $WGx=~/^\d+$/ and $WGy>=$WGx) {
											print STDERR $WGsubinfo, "Error: exon end > start: $WGz for mRNA $WGmrnaid\n";
											return $GffKit_failure;
										}
	#									$WGexonidhash{$WGz}=[$WGx, $WGy];
										@WGarr=();
										@WGarr=($WGref, '.', 'exon', $WGx, $WGy, $WGexonscore, ${$WGmrna}{$WGmrnaid}{'strand'}, '.', "ID=$WGz;Parent=$WGmrnaid");
										print WGGFF3OUT join ("\t", @WGarr), "\n";
										@WGarr=();
									}
								}
							}
							else {
								print STDERR $WGsubinfo, "Error: mRNA no exon regions: mRNA $WGmrnaid\n";
								return $GffKit_failure;
							}
						}
						else {
							print STDERR $WGsubinfo, "Error: mRNA no exons: mRNA $WGmrnaid\n";
							return $GffKit_failure;
						}
						### write CDS
						if (exists ${$WGcds}{$WGmrnaid}) {
							unless (exists ${$WGcds}{$WGmrnaid}{'reference'} and 
									(${$WGcds}{$WGmrnaid}{'reference'} eq $WGref)
							) {
								print STDERR $WGsubinfo, "Error: inconsistent CDS ref: mRNA $WGmrnaid\n";
								return $GffKit_failure;
							}
							unless (exists ${$WGcds}{$WGmrnaid}{'strand'} and 
									(${$WGcds}{$WGmrnaid}{'strand'} eq ${$WGgene}{$WGgeneid}{'strand'})
							) {
								print STDERR $WGsubinfo, "Error: inconsistent CDS strand: mRNA $WGmrnaid\n";
								return $GffKit_failure;
							}
							my $WGcdsscore='.';
							if (exists ${$WGcds}{$WGmrnaid}{'score'} and ${$WGcds}{$WGmrnaid}{'score'} =~/^\d+$/) {
								$WGcdsscore=${$WGcds}{$WGmrnaid}{'score'};
							}
							if (exists ${$WGcds}{$WGmrnaid}{'cds'} and 
								(scalar(keys %{${$WGcds}{$WGmrnaid}{'cds'}})>0)
							){
								foreach my $WGx (keys %{${$WGcds}{$WGmrnaid}{'cds'}}) {
									foreach my $WGy (keys %{${$WGcds}{$WGmrnaid}{'cds'}{$WGx}}) {
										unless ($WGy=~/^\d+$/ and $WGx=~/^\d+$/ and $WGy>=$WGx) {
											print STDERR $WGsubinfo, "Error: CDS end > start: $WGx-$WGy for mRNA $WGmrnaid\n";
											return $GffKit_failure;
										}
										my $WGcdsphase='.';
										if (exists ${$WGcds}{$WGmrnaid}{'phase'} and 
											exists ${$WGcds}{$WGmrnaid}{'phase'}{$WGx} and 
											exists ${$WGcds}{$WGmrnaid}{'phase'}{$WGx}{$WGy} and 
											${$WGcds}{$WGmrnaid}{'phase'}{$WGx}{$WGy}=~/^(0)|(1)|(2)$/
										) {
											$WGcdsphase=${$WGcds}{$WGmrnaid}{'phase'}{$WGx}{$WGy};
										}
										else {
											print STDERR $WGsubinfo, "Warnings: unknown CDS phase ",$WGx, '-', $WGy, " for mRNA $WGmrnaid\n";
										}
										@WGarr=();
										@WGarr=($WGref, '.', 'CDS', $WGx, $WGy, $WGcdsscore, ${$WGmrna}{$WGmrnaid}{'strand'}, $WGcdsphase, "ID=$WGmrnaid.cds;Parent=$WGmrnaid");
										print WGGFF3OUT join ("\t", @WGarr), "\n";
										@WGarr=();
									}
								}
							}
							else {
								print STDERR $WGsubinfo, "Error: mRNA no CDS region: mRNA $WGmrnaid\n";
								return $GffKit_failure;
							}
						}
						else {
							print STDERR $WGsubinfo, "warnings: mRNA no CDS: mRNA $WGmrnaid\n";
						}
					}
				}
				else {
					print STDERR $WGsubinfo, "Warnings: no mRNAs and exon/CDS not written: GENE $WGgeneid\n";
				}
			}
		}
	}

	return $GffKit_success;
}



### Exonerate GFF to final GFF
### &ExonerateGff3Reformat($cDNA_gff3_input, $CDS_gff3_input, $file_gff3_output)
### Global:
### Dependency:
### Note: 
### Return: 1=success, 0=failure
sub ExonerateGff3Reformat {
	my ($EGRgff_in_cDNA, $EGRgff_in_CDS, $EGRgff_out)=@_;
	
	my $EGRsubinfo='SUB(GffKit::ExonerateGff3Reformat)';
	my $EGRis_there_cds=0;
	my $EGRlinenum=0;
	local *EGRCDNAIN1; local *EGRCDNAIN2;
	local *EGRCDSIN1; local *EGRCDSIN2;
	local *EGROUT;
#	%EGRcdnareformat=($ref => ( ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene => (ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene.trans0001 => ++,
#(ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene.trans0002 => ++))
#	
	my %EGRgene2mrna=();
#	%EGRcdnareformat=($ref => ( gene000000007 => (gene => ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene;
#                                             mrna => ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene.trans0001
#                                             minc => 111
#                                             maxc => 222
#                                             )
#                          )
#                );
	my %EGRcdnareformat=();
	my %EGRcdsreformat=();
#	%EGRfinalreformat=($ref => ( gene ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene => (genemin)
#                          )
#                );
	my %EGRfinalreformat=();
#	%EGRcoordinates=($ref => ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene => ('minc' => $start,
#																					'maxc' => $end))
#							ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene.trans0001 => ('minc' => $start,
#																					'maxc' => $end))
	my %EGRcoordinates=();
#	%EGRprinted_genes=($ref => ( ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene => ++)
	my %EGRprinted_genes=();
#	%EGRexonid=($ref => ( $mRNA => $mRNA.exon0001++)
	my %EGRexonid=();
	
	unless (defined $EGRgff_in_cDNA and -s $EGRgff_in_cDNA) {
		print STDERR $EGRsubinfo, "Error: invalid cDNA GFF3 input file\n";
		return $GffKit_failure;
	}
	if (defined $EGRgff_in_CDS and -s $EGRgff_in_CDS) {
		$EGRis_there_cds=1;
	}
	else {
		print STDERR $EGRsubinfo, "Warnings: invalid CDS GFF3 input file\n";
	}
	unless (defined $EGRgff_out) {
		print STDERR $EGRsubinfo, "Error: invalid GFF3 output file\n";
		return $GffKit_failure;
	}
	unlink "$EGRgff_out" if (-e $EGRgff_out);
	
### read cDNA
	close EGRCDNAIN1 if (defined fileno(EGRCDNAIN1));
	unless (open (EGRCDNAIN1, "< $EGRgff_in_cDNA")) {
		print STDERR $EGRsubinfo, "Error: can not open cDNA GFF3 input file: $EGRgff_in_cDNA\n";
		return $GffKit_failure;
	}
	while (my $EGRline=<EGRCDNAIN1>) {
		$EGRlinenum++;
		chomp $EGRline;
		if ($EGRline=~/^#/) {
			next;
		}
		my @EGRarr=split(/\t/, $EGRline);
		unless (scalar(@EGRarr)==9) {
			print STDERR $EGRsubinfo, "Error: col !=9 at cDNA line($EGRlinenum), sure GFF3 format?\n";
			return $GffKit_failure;
		}
		
		if ($EGRarr[2]=~/^gene$/i) {
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*ID=//;$EGRgenename=~s/;.*$//;
			my $EGRmrnaname=$EGRarr[8];$EGRmrnaname=~s/^.*Name=//;$EGRmrnaname=~s/;.*$//;
			my $EGRgenefmt=$EGRmrnaname;$EGRgenefmt=~s/\.\w+$//;
			$EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'gene'}=$EGRgenefmt;
			$EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'mrna'}=$EGRmrnaname;
		}
		elsif ($EGRarr[2]=~/^exon$/i) {
			unless ($EGRarr[3]=~/^\d+$/ and $EGRarr[4]=~/^\d+$/ and $EGRarr[3]<=$EGRarr[4]) {
				print STDERR $EGRsubinfo, "Error: invalid cDNA coordinate at line($EGRlinenum)\n";
				return $GffKit_failure;
			}
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*Parent=//;$EGRgenename=~s/;.*$//;
			unless (exists $EGRcdnareformat{$EGRarr[0]} and exists $EGRcdnareformat{$EGRarr[0]}{$EGRgenename} and exists $EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'minc'} and $EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'minc'}<$EGRarr[3]) {
				$EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'minc'}=$EGRarr[3];
			}
			unless (exists $EGRcdnareformat{$EGRarr[0]} and exists $EGRcdnareformat{$EGRarr[0]}{$EGRgenename} and exists $EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'maxc'} and $EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'maxc'}>$EGRarr[4]) {
				$EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'maxc'}=$EGRarr[4];
			}
		}
		else {
			print STDERR $EGRsubinfo, "Error: unknown cDNA format at line($EGRlinenum)\n";
			return $GffKit_failure;
		}
	}
	close EGRCDNAIN1;
	
### read CDS
	$EGRlinenum=0;
	close EGRCDSIN1 if (defined fileno(EGRCDSIN1));
	unless (open (EGRCDSIN1, "< $EGRgff_in_CDS")) {
		print STDERR $EGRsubinfo, "Error: can not open CDS GFF3 input file: $EGRgff_in_CDS\n";
		return $GffKit_failure;
	}
	while (my $EGRline=<EGRCDSIN1>) {
		$EGRlinenum++;
		chomp $EGRline;
		if ($EGRline=~/^#/) {
			next;
		}
		my @EGRarr=split(/\t/, $EGRline);
		unless (scalar(@EGRarr)==9) {
			print STDERR $EGRsubinfo, "Error: col !=9 at CDS line($EGRlinenum), sure GFF3 format?\n";
			return $GffKit_failure;
		}
		
		if ($EGRarr[2]=~/^gene$/i) {
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*ID=//;$EGRgenename=~s/;.*$//;
			my $EGRmrnaname=$EGRarr[8];$EGRmrnaname=~s/^.*Name=//;$EGRmrnaname=~s/;.*$//;
			my $EGRgenefmt=$EGRmrnaname;$EGRgenefmt=~s/\.\w+$//;
			$EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'gene'}=$EGRgenefmt;
			$EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'mrna'}=$EGRmrnaname;
		}
		elsif ($EGRarr[2]=~/^exon$/i) {
			unless ($EGRarr[3]=~/^\d+$/ and $EGRarr[4]=~/^\d+$/ and $EGRarr[3]<=$EGRarr[4]) {
				print STDERR $EGRsubinfo, "Error: invalid CDS coordinate at line($EGRlinenum)\n";
				return $GffKit_failure;
			}
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*Parent=//;$EGRgenename=~s/;.*$//;
			unless (exists $EGRcdsreformat{$EGRarr[0]} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'minc'} and $EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'minc'}>$EGRarr[3]) {
				$EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'minc'}=$EGRarr[3];
			}
			unless (exists $EGRcdsreformat{$EGRarr[0]} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'maxc'} and $EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'maxc'}<$EGRarr[4]) {
				$EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'maxc'}=$EGRarr[4];
			}
		}
		else {
			print STDERR $EGRsubinfo, "Error: unknown format at line($EGRlinenum)\n";
			return $GffKit_failure;
		}
	}
	close EGRCDSIN1;
### Calculate Gene and mRNA coordinates
	foreach my $EGRchrom (keys %EGRcdnareformat) {
		foreach my $EGRgeneID1 (keys %{$EGRcdnareformat{$EGRchrom}}) {
			unless (exists $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'gene'} and $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'gene'}=~/^\S+$/) {
				print STDERR $EGRsubinfo, "Error: gene name REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'mrna'} and $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'mrna'}=~/^\S+$/) {
				print STDERR $EGRsubinfo, "Error: mRNA name REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'minc'} and $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'minc'}=~/^\d+$/) {
				print STDERR $EGRsubinfo, "Error: minc REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'maxc'} and $EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'maxc'}=~/^\d+$/) {
				print STDERR $EGRsubinfo, "Error: maxc REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			my $EGRcdna2gene=$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'gene'};
			my $EGRcdna2mrna=$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'mrna'};
			$EGRgene2mrna{$EGRchrom}{$EGRcdna2gene}{$EGRcdna2mrna}++;
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2gene} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'minc'} and $EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'minc'}<$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'minc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'minc'}=$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'minc'};
			}
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2gene} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'maxc'} and $EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'maxc'}>$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'maxc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcdna2gene}{'maxc'}=$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'maxc'};
			}
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'minc'} and $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'minc'}<$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'minc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'minc'}=$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'minc'};
			}
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna} and exists $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'maxc'} and $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'maxc'}>$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'maxc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcdna2mrna}{'maxc'}=$EGRcdnareformat{$EGRchrom}{$EGRgeneID1}{'maxc'};
			}
		}
	}
	foreach my $EGRchrom (keys %EGRcdsreformat) {
		foreach my $EGRgeneID1 (keys %{$EGRcdsreformat{$EGRchrom}}) {
			unless (exists $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'gene'} and $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'gene'}=~/^\S+$/) {
				print STDERR $EGRsubinfo, "Error: gene name REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'mrna'} and $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'mrna'}=~/^\S+$/) {
				print STDERR $EGRsubinfo, "Error: mRNA name REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'minc'} and $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'minc'}=~/^\d+$/) {
				print STDERR $EGRsubinfo, "Error: minc REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'maxc'} and $EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'maxc'}=~/^\d+$/) {
				print STDERR $EGRsubinfo, "Error: maxc REF $EGRchrom ID $EGRgeneID1\n";
				return $GffKit_failure;
			}
			my $EGRcds2gene=$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'gene'};
			my $EGRcds2mrna=$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'mrna'};
			$EGRgene2mrna{$EGRchrom}{$EGRcds2gene}{$EGRcds2mrna}++;
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2gene} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'minc'} and $EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'minc'}<$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'minc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'minc'}=$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'minc'};
			}
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2gene} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'maxc'} and $EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'maxc'}>$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'maxc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcds2gene}{'maxc'}=$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'maxc'};
			}
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2mrna} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'minc'} and $EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'minc'}<$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'minc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'minc'}=$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'minc'};
			}
			unless (exists $EGRcoordinates{$EGRchrom} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2mrna} and exists $EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'maxc'} and $EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'minc'}=~/^\d+$/ and $EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'maxc'}>$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'maxc'}) {
				$EGRcoordinates{$EGRchrom}{$EGRcds2mrna}{'maxc'}=$EGRcdsreformat{$EGRchrom}{$EGRgeneID1}{'maxc'};
			}
		}
	}
	
### Debug
	if (0) {print "Test: \%EGRcoordinates\n"; print Dumper \%EGRcoordinates; print "\n";}
	
### Start reformat
	close EGROUT if (defined fileno(EGROUT));
	unless (open(EGROUT, "> $EGRgff_out")) {
		print STDERR $EGRsubinfo, "Error: can not write GFF3 output file: $EGRgff_out\n";
		return $GffKit_failure;
	}
	print EGROUT "##gff-version 3\n";


### Reformat cDNA
	close EGRCDNAIN2 if (defined fileno(EGRCDNAIN2));
	unless (open (EGRCDNAIN2, "< $EGRgff_in_cDNA")) {
		print STDERR $EGRsubinfo, "Error: can not open cDNA GFF3 input file2: $EGRgff_in_cDNA\n";
		return $GffKit_failure;
	}
	while (my $EGRline=<EGRCDNAIN2>) {
		chomp $EGRline;
		if ($EGRline=~/^#/) {
			next;
		}
		my @EGRarr=split(/\t/, $EGRline);
		if ($EGRarr[2]=~/^gene$/i) {
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*ID=//;$EGRgenename=~s/;.*$//;
			my $EGRmrnaname=$EGRarr[8];$EGRmrnaname=~s/^.*Name=//;$EGRmrnaname=~s/;.*$//;
			my $EGRgenefmt=$EGRmrnaname;$EGRgenefmt=~s/\.\w+$//;
			unless (exists $EGRcoordinates{$EGRarr[0]} and exists $EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}) {
				print STDERR $EGRsubinfo, "Error: invalid GENE coordinates: REF $EGRarr[0] GENE $EGRgenename mRNA $EGRmrnaname GENENAME $EGRgenefmt\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'minc'} and $EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'minc'}=~/^\d+$/) {
				print STDERR $EGRsubinfo, "Error: invalid GENE start: REF $EGRarr[0] GENE $EGRgenename mRNA $EGRmrnaname GENENAME $EGRgenefmt\n";
				return $GffKit_failure;
			}
			unless (exists $EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'maxc'} and $EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'maxc'}=~/^\d+$/) {
				print STDERR $EGRsubinfo, "Error: invalid GENE end: REF $EGRarr[0] GENE $EGRgenename mRNA $EGRmrnaname GENENAME $EGRgenefmt\n";
				return $GffKit_failure;
			}
			unless ($EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'minc'}<=$EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'maxc'}) {
				print STDERR $EGRsubinfo, "Error: invalid GENE start>end: REF $EGRarr[0] GENE $EGRgenename mRNA $EGRmrnaname GENENAME $EGRgenefmt\n";
				return $GffKit_failure;
			}
			unless (exists $EGRprinted_genes{$EGRarr[0]} and exists $EGRprinted_genes{$EGRarr[0]}{$EGRgenefmt}) {
				print EGROUT $EGRarr[0], "\t", $EGRarr[1], "\tgene\t",$EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'minc'}, "\t",$EGRcoordinates{$EGRarr[0]}{$EGRgenefmt}{'maxc'}, "\t.\t",$EGRarr[6], "\t.\tID=$EGRgenefmt;Name=$EGRgenefmt\n";
				$EGRprinted_genes{$EGRarr[0]}{$EGRgenefmt}++;
			}

			
			unless (exists $EGRcoordinates{$EGRarr[0]} and exists $EGRcoordinates{$EGRarr[0]}{$EGRmrnaname} and exists $EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'minc'} and $EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'minc'}=~/^\d+$/ and exists $EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'maxc'} and $EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'maxc'}=~/^\d+$/ and $EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'minc'}<=$EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'maxc'}) {
				print STDERR $EGRsubinfo, "Error: invalid mRNA coordinates: REF $EGRarr[0] GENE $EGRgenename mRNA $EGRmrnaname GENENAME $EGRgenefmt\n";
				return $GffKit_failure;
			}
			print EGROUT $EGRarr[0], "\t", $EGRarr[1], "\tmRNA\t",$EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'minc'}, "\t",$EGRcoordinates{$EGRarr[0]}{$EGRmrnaname}{'maxc'}, "\t.\t",$EGRarr[6], "\t.\tID=$EGRmrnaname;Parent=$EGRgenefmt\n";
		}
		elsif ($EGRarr[2]=~/^exon$/i) {
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*Parent=//;$EGRgenename=~s/;.*$//;
			unless (exists $EGRcdnareformat{$EGRarr[0]} and exists $EGRcdnareformat{$EGRarr[0]}{$EGRgenename} and exists $EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'mrna'}) {
				print STDERR $EGRsubinfo, "Error: unknown mRNA name for exon: REF $EGRarr[0] GENE $EGRgenename\n";
				return $GffKit_failure;
			}
			my $EGRthismRNA=$EGRcdnareformat{$EGRarr[0]}{$EGRgenename}{'mrna'};
			unless (exists $EGRexonid{$EGRarr[0]} and exists $EGRexonid{$EGRarr[0]}{$EGRthismRNA}) {
				$EGRexonid{$EGRarr[0]}{$EGRthismRNA}=$EGRthismRNA.'.exon0001';
			}
			print EGROUT $EGRarr[0], "\t", $EGRarr[1], "\texon\t",$EGRarr[3], "\t",$EGRarr[4], "\t.\t",$EGRarr[6], "\t.\tID=",$EGRexonid{$EGRarr[0]}{$EGRthismRNA}, ";Parent=$EGRthismRNA\n";
			my $EGRexonid;
			if ($EGRexonid{$EGRarr[0]}{$EGRthismRNA}=~/\.(exon\d+)$/) {
				$EGRexonid=$1;
			}
			else {
				print STDERR $EGRsubinfo, "Error: invalid exon ID for exon: REF $EGRarr[0] GENE $EGRgenename\n";
				return $GffKit_failure;
			}
			$EGRexonid++;
			$EGRexonid{$EGRarr[0]}{$EGRthismRNA}=$EGRthismRNA.'.'.$EGRexonid;
		}
	}
	close EGRCDNAIN2;
	
### Reformat CDS
	close EGRCDSIN2 if (defined fileno(EGRCDSIN2));
	unless (open (EGRCDSIN2, "< $EGRgff_in_CDS")) {
		print STDERR $EGRsubinfo, "Error: can not open cDNA GFF3 input file2: $EGRgff_in_CDS\n";
		return $GffKit_failure;
	}
	while (my $EGRline=<EGRCDSIN2>) {
		chomp $EGRline;
		if ($EGRline=~/^#/) {
			next;
		}
		my @EGRarr=split(/\t/, $EGRline);
		if ($EGRarr[2]=~/^gene$/i) {
			next;
		}
		elsif ($EGRarr[2]=~/^exon$/i) {
			my $EGRgenename=$EGRarr[8];$EGRgenename=~s/^.*Parent=//;$EGRgenename=~s/;.*$//;
			unless (exists $EGRcdsreformat{$EGRarr[0]} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'mrna'} and exists $EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'gene'}) {
				print STDERR $EGRsubinfo, "Error: unknown mRNA name for CDS: REF $EGRarr[0] GENE $EGRgenename\n";
				return $GffKit_failure;
			}
			my $EGRthismRNA=$EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'mrna'};
			my $EGRgenefmt=$EGRcdsreformat{$EGRarr[0]}{$EGRgenename}{'gene'};
			unless (exists $EGRprinted_genes{$EGRarr[0]} and exists $EGRprinted_genes{$EGRarr[0]}{$EGRgenefmt}) {
				print STDERR $EGRsubinfo, "Warnings: CDS no gene & mRNA: REF ", $EGRarr[0], " GENE $EGRgenefmt\n";
			}
			print EGROUT $EGRarr[0], "\t", $EGRarr[1], "\tCDS\t",$EGRarr[3], "\t",$EGRarr[4], "\t.\t",$EGRarr[6], "\t.\tID=",$EGRthismRNA, ".cds;Parent=$EGRthismRNA\n";
		}
	}
	close EGRCDSIN2;
	close EGROUT;
}
###cDNA
##gff-version 3
#chr3DL_BAC_00000198	exonerate:est2genome	gene	287556	290776	9647	-	.	#ID=gene000000007;Name=ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene.trans0001
#chr3DL_BAC_00000198	exonerate:est2genome	exon	290049	290776	.	-	.	ID=exon000000050;Parent=gene000000007
#chr3DL_BAC_00000198	exonerate:est2genome	exon	289690	289869	.	-	.	ID=exon000000051;Parent=gene000000007
#chr3DL_BAC_00000198	exonerate:est2genome	exon	289301	289402	.	-	.
#CDS
##gff-version 3
#chr3DL_BAC_00000198	exonerate:est2genome	gene	288095	290366	5360	-	.	ID=gene000000005;Name=ctg00006_019_O17.flt500_scaffold_1_6000-9700.gene.trans0001
#chr3DL_BAC_00000198	exonerate:est2genome	exon	290049	290366	.	-	.	ID=exon000000036;Parent=gene000000005
#chr3DL_BAC_00000198	exonerate:est2genome	exon	289690	289869	.	-	.	ID=exon000000037;Parent=gene000000005
#chr3DL_BAC_00000198	exonerate:est2genome	exon	289301	289402	.	-	.
#my $GRSsubinfo='SUB(GffKit::GffReverseStrand)';
#my $GffKit_success=1; $GffKit_failure=0; $GffKit_debug=0;






#my $GRSsubinfo='SUB(GffKit::GffReverseStrand)';
#my $GffKit_success=1; $GffKit_failure=0; $GffKit_debug=0;
1;
