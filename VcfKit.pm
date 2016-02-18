# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::VcfKit

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
package FuhaoPerl5Lib::VcfKit;
use strict;
use warnings;
use Data::Dumper;
use FuhaoPerl5Lib::BamKit qw/SplitCigar/;
use FuhaoPerl5Lib::CmdKit;
use FuhaoPerl5Lib::FileKit;
use FuhaoPerl5Lib::FastaKit qw(CreateFastaRegion);
use FuhaoPerl5Lib::StatsKit qw /CallCorrelation/;
use Bio::DB::Sam;
use Exporter;
use Cwd;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150820';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(IndexVcf BgzipVcf ExtractVcf ExtractVcfFromFile ReadVcf RunHapCompass LoadVcf RunFreebayes HcFragments ReadVariantType HapcompassVcf GroupFragments CorrectAlleles ReadHcOut);
%EXPORT_TAGS = ( DEFAULT => [qw(IndexVcf BgzipVcf ExtractVcf ExtractVcfFromFile ReadVcf RunHapCompass LoadVcf RunFreebayes HcFragments ReadVariantType HapcompassVcf )],
                 ALL    => [qw(IndexVcf BgzipVcf ExtractVcf ExtractVcfFromFile ReadVcf RunHapCompass LoadVcf RunFreebayes HcFragments ReadVariantType HapcompassVcf)]);


my $VcfKit_success=1;
my $VcfKit_failure=0;
my $VcfKit_debug=0; #checge to 1 if debug


### Index Vcf using tabix
### &IndexVcf (input.vcf.gz, [path_tabix])
### Global: path_tabix
### Dependency: &exec_cmd_return
### Note: input.vcf.gz must be bgzip-ped
sub IndexVcf {
	my ($IVvcfin, $IVpath_tabix)=@_;
	
	my $IVsubinfo='SUB(VcfKit::IndexVcf)';
	$IVpath_tabix='tabix' unless (defined $IVpath_tabix);
	
	unless (defined $IVvcfin and -s $IVvcfin) {
		print STDERR "${IVsubinfo}Error: VCF input not found\n";
		return $VcfKit_failure;
	}
	unlink "$IVvcfin.tbi" if (-e "$IVvcfin.tbi");
	
	if (! &exec_cmd_return("$IVpath_tabix -p vcf $IVvcfin")) {
		print STDERR "${IVsubinfo}Error: VCF input not found\n";
		return $VcfKit_failure;
	}
	elsif (! -s "$IVvcfin.tbi") {
		print STDERR "${IVsubinfo}Error: VCF input not found\n";
		return $VcfKit_failure;
	}
	return $VcfKit_success;
}



### Compress VCF into .vcf.gz using bgzip
### &BgzipVcf (input.vcf, output.vcf.gz, [path_bgzip])
### Global: 
### Dependency: &exec_cmd_return
### Note:
sub BgzipVcf {
	my ($BVvcfin, $BVvcfout, $BVpath_bgzip)=@_;
	
	my $BVsubinfo='SUB(VcfKit::BgzipVcf)';
	$BVpath_bgzip='bgzip' unless (defined $BVpath_bgzip);
	
	unless (defined $BVvcfin and -s $BVvcfin) {
		print STDERR "${BVsubinfo}Error: VCF input not found: $BVvcfin\n";
		return $VcfKit_failure;
	}
	unlink $BVvcfout if (-e $BVvcfout);
	
	if (&exec_cmd_return("$BVpath_bgzip -c $BVvcfin > $BVvcfout")) {
		print STDERR "${BVsubinfo}Error: bgzip running error\n";
		return $VcfKit_failure;
	}
	elsif (! -s $BVvcfout) {
		print STDERR "${BVsubinfo}Error: bgzip output error\n";
		return $VcfKit_failure;
	}
	return $VcfKit_success;
}



### Extract VCF subset with tabix
### &ExtractVcf (indexed.vcf.gz, fasta_ids, output.vcf.gz, [path_tabix], [path_bgzip])
### Global: path_tabix, $path_bgzip
### Dependency: &exec_cmd_return
### Note: input.vcf.gz must be bgzip-ped and indexed with 'tabix -p vcf *'
### Note: $EVids are space delimited
sub ExtractVcf {
	my ($EVvcfin, $EVids, $EVvcfout, $EVpath_tabix, $EVpath_bgzip)=@_;
	
	my $EVsubinfo='SUB(VcfKit::ExtractVcf)';
	$EVpath_tabix='tabix' unless (defined $EVpath_tabix);
	$EVpath_bgzip='bgzip' unless (defined $EVpath_bgzip);
	
	unless (defined $EVvcfin and -s $EVvcfin) {
		print STDERR "${EVsubinfo}Error: VCF input not found: $EVvcfin\n";
		return $VcfKit_failure;
	}
	unless ($EVids=~/\S+/) {
		print STDERR "${EVsubinfo}Error: invalid extract ids\n";
		return $VcfKit_failure;
	}
	if ($EVvcfin=~/\.vcf$/i) {
		if (! &BgzipVcf($EVvcfin, "$EVvcfin.gz", $EVpath_bgzip)) {
			print STDERR "${EVsubinfo}Error: VCF bgzip error: $EVvcfin\n";
			return $VcfKit_failure;
		}
		$EVvcfin="$EVvcfin.gz";
	}
	unless (-s "$EVvcfin.tbi") {
		if (! &IndexVcf($EVvcfin, $EVpath_tabix)) {
			print STDERR "${EVsubinfo}Error: VCF bgzip error: $EVvcfin\n";
			return $VcfKit_failure;
		}
	}
	unlink $EVvcfout if (-e $EVvcfout);
	
	if (! &exec_cmd_return("$EVpath_tabix $EVvcfin $EVids | $EVpath_bgzip > $EVvcfout")) {
		print STDERR "${EVsubinfo}Error: tabix extract running error\n";
		return $VcfKit_failure;
	}
	elsif ( ! -e $EVvcfout) {
		print STDERR "${EVsubinfo}Error: tabix extract output error\n";
		return $VcfKit_failure;
	}
	return $VcfKit_success;
}




### Extract VCF subset with tabix from a file containing ID list, 1 ID per line
### &ExtractVcf (indexed.vcf.gz, fasta_id.input.file, output.vcf.gz, [path_tabix], [path_bgzip])
### Global: path_tabix, $path_bgzip
### Dependency: &exec_cmd_return
### Note: input.vcf.gz must be bgzip-ped and indexed with 'tabix -p vcf *'
### Note: fasta_id.input.file: 1 ID per line
sub ExtractVcfFromFile {
	my ($EVFFvcfin, $EVFFidfile, $EVFFvcfout, $EVFFpath_tabix, $EVFFpath_bgzip)=@_;
	
	my $EVFFsubinfo='SUB(VcfKit::ExtractVcfFromFile)';
	$EVFFpath_tabix='tabix' unless (defined $EVFFpath_tabix);
	$EVFFpath_bgzip='bgzip' unless (defined $EVFFpath_bgzip);
	
	unless (defined $EVFFvcfin and -s $EVFFvcfin) {
		print STDERR $EVFFsubinfo, "Error: VCF input not found: $EVFFvcfin\n";
		return $VcfKit_failure;
	}
	if ($EVFFvcfin=~/\.vcf$/i) {
		unless (&BgzipVcf($EVFFvcfin, "$EVFFvcfin.gz", $EVFFpath_bgzip)) {
			print STDERR $EVFFsubinfo, "Error: VCF bgzip error: $EVFFvcfin\n";
			return $VcfKit_failure;
		}
		$EVFFvcfin="$EVFFvcfin.gz";
	}
	unless (-s $EVFFidfile) {
		print STDERR $EVFFsubinfo, "Error: invalid ID list file\n";
		return $VcfKit_failure;
	}
	unless (-s "$EVFFvcfin.tbi") {
		unless (&IndexVcf($EVFFvcfin, $EVFFpath_tabix)) {
			print STDERR $EVFFsubinfo, "Error: VCF bgzip error: $EVFFvcfin\n";
			return $VcfKit_failure;
		}
	}
	unlink $EVFFvcfout if (-e $EVFFvcfout);
	
	unless (exec_cmd_return("$EVFFpath_tabix -h $EVFFvcfin `cat $EVFFidfile` | $EVFFpath_bgzip > $EVFFvcfout")) {
		print STDERR $EVFFsubinfo, "Error: tabix extract running error\n";
		return $VcfKit_failure;
	}
	unless (-e $EVFFvcfout) {
		print STDERR $EVFFsubinfo, "Error: tabix extract output error\n";
		return $VcfKit_failure;
	}
	return $VcfKit_success;
}




###Read vcf into hash, return index
###&ReadVcf(vcf_file, [minmqm])
###Global: 
###Dependancy: 
###VCF
#0		1	2	3	4	5		6		7		8		9
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ParLeaf
#0	ParMerge_MblContig65
#1	1106
#2	.
#3	C
#4	T
#5	310.79
#6	.
#7		AB=0.44;ABP=3.79203;AC=1;AF=0.333333;AN=3;AO=11;CIGAR=1X;DP=25;DPB=25;DPRA=0;EPP=19.0002;EPPR=5.49198;GTI=0;LEN=1;MEANALT=1;MQM=41.4545;MQMR=41.8571;NS=1;NUMALT=1;ODDS=4.15888;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=431;QR=517;RO=14;RPL=8;RPP=7.94546;RPPR=5.49198;RPR=3;RUN=1;SAF=4;SAP=4.78696;SAR=7;SRF=8;SRP=3.63072;SRR=6;TYPE=snp;technology.ILLUMINA=1	
#8	GT:DP:RO:QR:AO:QA
#9	0/0/1:25:14:517:11:431
sub ReadVcf {
	my ($RVvcf, $RVvcfqual)=@_;
	
	local *VCF;
	my $RVsubinfo='SUB(VcfKit::ReadVcf)';
	my %RVret_vcf=();
	$RVvcfqual=20 unless (defined $RVvcfqual);
	my $RVgeno_delimiter=',';
##COMMENT: guess format
	unless (defined $RVvcf and -s $RVvcf) {
		print STDERR "${RVsubinfo}Error: invalid VCF file\n";
		return $VcfKit_failure;
	}
	close VCF if (defined fileno(VCF));
	if ($RVvcf=~/\.gz$|\.gzip/i) {
		unless (open (VCF, "zcat $RVvcf|")) {
			print STDERR "${RVsubinfo}Error: can not open vcf.gz file\n";
			return $VcfKit_failure;
		}
	}
	elsif ($RVvcf=~/\.vcf$/i) {
		unless (open (VCF, "cat $RVvcf|")) {
			print STDERR "${RVsubinfo}Error: can not open vcf file\n";
			return $VcfKit_failure;
		}
	}
	else {
		print STDERR "${RVsubinfo}Error: can not guess vcf file format\n";
		return $VcfKit_failure;
	}
##COMMENT: reading VCF
	my $RVnum_line=0;
	while (my $RVline=<VCF>) {
		$RVnum_line++;
		next if ($RVline =~/^#/);
		chomp $RVline;
		my @RVarr=();
		@RVarr=split(/\t/, $RVline);
##COMMENT: comment next line if column number <10
		if (scalar(@RVarr)<10) {
			print STDERR "${RVsubinfo}Error: Wrong variation at line $RVnum_line in $RVvcf\n";
			next;
		}
##COMMENT: 
		my ($RVchr, $RVpos, $RVref, $RVvar)=('', '', '', '');
		($RVchr, $RVpos, $RVref, $RVvar)=($RVarr[0], $RVarr[1], $RVarr[3], $RVarr[4]);
		if (exists ${$RVret_vcf{$RVchr}}{$RVpos} ) {
			print STDERR "${RVsubinfo}Error: repeated $RVchr - $RVpos at line $RVnum_line in $RVvcf\n";
			return $VcfKit_failure;
		}
		@{${$RVret_vcf{$RVchr}}{$RVpos}}=($RVref, $RVvar);
##COMMENT: Parse vcf
		for (my $i=9; $i<scalar(@RVarr); $i++) {
			#GT:DP:RO:QR:AO:QA=0/0/1:90:47:1816:42:1060cd Clu							
			if ($RVarr[$i] eq '.') {
				push (@{${$RVret_vcf{$RVchr}}{$RVpos}}, '.');
				next;
			}
			my @RVarr2=();
			@RVarr2=split(/:/, $RVarr[$i]);
			#geno: $RVarr2[0]
			#readdepth: $RVarr2[1]
			#ref count: $RVarr2[2]
			#var count: $RVarr2[4]
##COMMENT: following two lines if not ignore those sites with only one var allele
#			unless (&IsZeroIn($RVarr2[2])) {
#				delete ${$RVret_vcf{$RVchr}}{$RVpos};
#				last;
#			}
#			unless (&IsZeroIn($RVarr2[4])){
#				delete ${$RVret_vcf{$RVchr}}{$RVpos};
#				last;
#			}
			if (defined $RVarr[5] and $RVarr[5]<$RVvcfqual) {##COMMENT: check mapping quality; would take the allele with max count
				my @arr3=split(/$RVgeno_delimiter/, $RVvar);
				my @arr4=split(/$RVgeno_delimiter/, $RVarr2[4]);
				unless (scalar(@arr3) == scalar(@arr4)) {
					print STDERR "${RVsubinfo}Error: VCF allele count error $RVchr - $RVpos at line $RVnum_line in $RVvcf\n";
					return $VcfKit_failure;
				}
				my %RVallele_count=();
				for (my $RVj=0; $RVj<scalar(@arr4); $RVj++) {
					$RVallele_count{$RVj+1}=$arr4[$RVj];
				}
				my $RVbest_allele=0;
				my $RVmax_count=$RVarr2[2];
				foreach (keys %RVallele_count) {
					if ($RVallele_count{$_} > $RVmax_count) {
						$RVbest_allele=$_;
						$RVmax_count=$RVallele_count{$_};
					}
				}
				my @RVarr5=split(/\//, $RVarr2[0]);
				my @RVarr6=();
				foreach (@RVarr5) {
					push (@RVarr6, $RVbest_allele);
				}
				$RVarr2[0]=join('/', @RVarr6);
			}
			push (@{${$RVret_vcf{$RVchr}}{$RVpos}}, $RVarr2[0]);
		}
	}
	close VCF;
	
	#	print "vcf input\n";###test###
	if ($VcfKit_debug) {
		foreach my $RVchrom (keys %RVret_vcf) {
			my @RVarr3=keys %{$RVret_vcf{$RVchrom}};
			@RVarr3=sort {$a<=>$b} @RVarr3;
			print "${RVsubinfo}Reference: $RVchrom\n";
			print "${RVsubinfo}Number of variations on $RVchrom: ", scalar(@RVarr3), "\n";
			foreach (@RVarr3) {
				print $RVchrom, "\t", "$_\t${${$RVret_vcf{$RVchrom}}{$_}}[0]\t${${$RVret_vcf{$RVchrom}}{$_}}[1]\t${${$RVret_vcf{$RVchrom}}{$_}}[2]\n";
			}
		}
	}
	return ($VcfKit_success, \%RVret_vcf);
###return
###%vcf=(chr1 => (pos1 => (ref, var, gen/gen2, ...),
###				  pos2 => (ref, var, gen/gen2, ...),
###				...
###				),
###		chr2 => (pos1 => (ref, var, gen/gen2, ...),
###				  pos2 => (ref, var, gen/gen2, ...),
###				...
###				),
###		...
###		);
###$vcf{chr}->{pos}=(ref, var, gen/gen2, ...);
###Be noted: var could be comma-delimited string, in case of two variants
}



### Load vcf into hash, return index, Upgrade version of ReadVcf
### &LoadVcf(vcf_file, minMQM, ReadCol10Only)
### Global: 
### Dependancy:
### Note: 20150617
### VCF
#0		1	2	3	4	5		6		7		8		9
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ParLeaf
#0	ParMerge_MblContig65
#1	1106
#2	.
#3	C
#4	T
#5	310.79
#6	.
#7		AB=0.44;ABP=3.79203;AC=1;AF=0.333333;AN=3;AO=11;CIGAR=1X;DP=25;DPB=25;DPRA=0;EPP=19.0002;EPPR=5.49198;GTI=0;LEN=1;MEANALT=1;MQM=41.4545;MQMR=41.8571;NS=1;NUMALT=1;ODDS=4.15888;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=431;QR=517;RO=14;RPL=8;RPP=7.94546;RPPR=5.49198;RPR=3;RUN=1;SAF=4;SAP=4.78696;SAR=7;SRF=8;SRP=3.63072;SRR=6;TYPE=snp;technology.ILLUMINA=1	
#8	GT:DP:RO:QR:AO:QA
#9	0/0/1:25:14:517:11:431
###OUT
##CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
#20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.
#%{hash}=( chr => (pos => (readgroup => ('ref' => ref_allele,
#									     'alt' => var_allele(s),
#						   		     	 'qual' => MQM/mapping quality,
#						   				 'filter' => Filter,
#						   				 'info' => information,
#						   				 'readgroup/colnum' => ('genotypes' => genotypes
#															'GT'	=> genometype,
#						  				 					'GQ' => Genotype_quality,
#						   				 					'DP'	=> Read_depth,
#						   				 					'RO' => Reference allele observation count,
#						   				 					'GR' => Sum of quality of the reference observations,
#						   				 					'AO' => Alternate allele observation count,
#						   				 					'QA' => Sum of quality of the alternate observations, ...
#														)
#										'num2allele' => (	0 => ref_allele,
#															1 => allele1, ...
#														)
#										'allele2num' => (	re_allele => 0,
#															allele1 => 1, ...
#														)
#						   				)
#						   )
#				)
#		)
sub LoadVcf {
	my ($LVvcf, $LVvcfqual, $LVonly10thcol)=@_;
	
	local *VCF;
	my $LVsubinfo='SUB(VcfKit::LoadVcf)';
	my %LVret_vcf=();
	$LVvcfqual=20 unless (defined $LVvcfqual);
	my $LVgeno_delimiter=',';
	$LVonly10thcol=1 unless (defined $LVonly10thcol);
	
##COMMENT: guess format
	unless (defined $LVvcf and -s $LVvcf) {
		print STDERR "${LVsubinfo}Error: invalid VCF file\n";
		return $VcfKit_failure;
	}
	close VCF if (defined fileno(VCF));
	if ($LVvcf=~/(\.gz$)|(\.gzip$)/i) {
		unless (open (VCF, "zcat $LVvcf |")) {
			print STDERR "${LVsubinfo}Error: can not open vcf.gz file\n";
			return $VcfKit_failure;
		}
	}
	elsif ($LVvcf=~/\.vcf$/i) {
		unless (open (VCF, "cat $LVvcf|")) {
			print STDERR "${LVsubinfo}Error: can not open vcf file\n";
			return $VcfKit_failure;
		}
	}
	else {
		print STDERR "${LVsubinfo}Error: can not guess vcf file format\n";
		return $VcfKit_failure;
	}

##COMMENT: reading VCF
	my $LVnum_line=0;
	my @LVreadgroup=();
	my @LVcolgroup=();
	my $LVreadheader=0;
	while (my $LVline=<VCF>) {
		$LVnum_line++;
		chomp $LVline;
		if ($LVline =~/^#/) {
			if ($LVline =~ /^#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT.*$/ and $LVreadheader ==0) {
#				print "${LVsubinfo}Test: header starts\n";### For test ###
				my @LVarr3=split(/\t/, $LVline);
				LU1: {for (my $LVi=0; $LVi<(scalar(@LVarr3)-9); $LVi++) {
					if (defined $LVarr3[$LVi+9] and $LVarr3[$LVi+9] ne '') {
#						print "${LVsubinfo}Test: header: ".$LVarr3[ $LVi + 9 ]."\n"; ### For test ###
						$LVreadgroup[$LVi]=$LVarr3[$LVi+9];
					}
					else {
						$LVreadgroup[$LVi]="col".$LVi+10;
					}
					last LU1 if ($LVonly10thcol==1);
				}}###LU1
				$LVreadheader=1;
			}
			next;
		}
		
		my @LVarr=();
		@LVarr=split(/\t/, $LVline);
##COMMENT: comment next line if column number <10
		if (scalar(@LVarr)<10) {
			print STDERR "${LVsubinfo}Error: Wrong variation at $LVnum_line in $LVvcf\n";
			next;
		}
##COMMENT: check if existing sites
		if (exists ${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}) {
			print STDERR "${LVsubinfo}Error: repeated $LVarr[0] - $LVarr[1] at $LVnum_line in $LVvcf\n";
			return $VcfKit_failure;
		}
		${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'ref'}=$LVarr[3];
		${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'var'}=$LVarr[4];
		${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'qual'}=$LVarr[5];
#		${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'filter'}=$LVarr[6]; ### Only uncomment when you need this column
#		${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'info'}=$LVarr[7]; ### Only uncomment when you need this column
		${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'num2allele'}}{0}=$LVarr[3];
		${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'allele2num'}}{$LVarr[3]}=0;
		my @LVvar_alleles=split(/,/, $LVarr[4]);
#		print "${LVsubinfo}Test: alleles: $LVarr[3], @LVvar_alleles\n"; ### For test ###
		for (my $LVindvar=1;$LVindvar<=scalar(@LVvar_alleles);$LVindvar++) {
			if (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'num2allele'}}{$LVindvar}) {
				print STDERR "${LVsubinfo}Error: existing num2allele: $LVindvar at Chr:pos $LVarr[0]:$LVarr[1] of $LVnum_line in $LVvcf\n";
				return $VcfKit_failure;
			}
			else {
				${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'num2allele'}}{$LVindvar}=$LVvar_alleles[$LVindvar-1];
			}
			if (${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'allele2num'}}{$LVvar_alleles[$LVindvar-1]}) {
				print STDERR "${LVsubinfo}Error: existing allele2num: ".$LVvar_alleles[$LVindvar-1]." at Chr:pos $LVarr[0]:$LVarr[1] of $LVnum_line in $LVvcf\n";
				return $VcfKit_failure;
			}
			else {
				${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'allele2num'}}{$LVvar_alleles[$LVindvar-1]}=$LVindvar;
			}
		}
		my @LVformat=split(/:/, $LVarr[8]);
##COMMENT: Parse vcf
		LU2: {for (my $LVx=0; $LVx<(scalar(@LVarr)-9); $LVx++) {
			#GT:DP:RO:QR:AO:QA=0/0/1:90:47:1816:42:1060
			my $LVrg = (defined $LVreadgroup[$LVx] and $LVreadgroup[$LVx]=~/^\S+$/) ? $LVreadgroup[$LVx] : 'col'.($LVx + 10);
			if ($LVarr[$LVx+9] eq '.') {
				${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{'GT'}='.';
				${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{'genotypes'}='.';
				unless (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'GT'}) {
					${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'GT'}='.';
				}
				unless (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'genotypes'}) {
					${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'genotypes'}='.';
				}
				if ($LVonly10thcol==1) {
					last LU2;
				}
				else {
					next;
				}
			}
			my @LVarr2=();
			@LVarr2=split(/:/, $LVarr[$LVx+9]);
			#geno: $LVarr2[0]
			#readdepth: $LVarr2[1]
			#ref count: $LVarr2[2]
			#var count: $LVarr2[4]
##COMMENT: following two lines if not ignore those sites with only one var allele
#			unless (&IsZeroIn($LVarr2[2])) {
#				delete ${$LVret_vcf{$LVarr[0]}}{$LVarr[1]};
#				last;
#			}
#			unless (&IsZeroIn($LVarr2[4])){
#				delete ${$LVret_vcf{$LVarr[0]}}{$LVarr[1]};
#				last;
#			}
			unless (scalar(@LVarr2) == scalar(@LVformat)) {
				print STDERR "${LVsubinfo}Error: FORMAT != genometype line at $LVnum_line in $LVvcf\n";
				return $VcfKit_failure;
			}
#			print "${LVsubinfo}Test: readgroup: $LVrg\n"; ### For test ###
			for (my $LVy=0; $LVy<scalar(@LVformat);$LVy++) {
				if (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{$LVformat[$LVy]}) {
					print STDERR "${LVsubinfo}Error: existing hash: Format: $LVformat[$LVy]; readgroup: $LVrg; chr:pos $LVarr[0]: $LVarr[1] at $LVnum_line in $LVvcf\n";
					return $VcfKit_failure;
				}
				else {
					${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{$LVformat[$LVy]}=$LVarr2[$LVy];
					
				}
				unless (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{$LVformat[$LVy]}) {
					${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{$LVformat[$LVy]}=$LVarr2[$LVy];
				}
			}
			if (defined $LVarr[5] and $LVarr[5]<$LVvcfqual) {##COMMENT: check mapping quality; would take the allele with max count
				my @arr3=split(/,/, $LVarr[4]);
				my @arr4=split(/,/, ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{'AO'});
				unless (scalar(@arr3) == scalar(@arr4)) {
					print STDERR "${LVsubinfo}Error: VCF allele count $LVarr[0] - $LVarr[1] at $LVnum_line in $LVvcf\n";
					return $VcfKit_failure;
				}
				
				my %LVallele_count=();
				for (my $LVj=0; $LVj<scalar(@arr4); $LVj++) {
					$LVallele_count{$LVj+1}=$arr4[$LVj];
				}
				my $LVbest_allele=0;
				my $LVmax_count=$LVarr2[2];
				foreach (keys %LVallele_count) {
					if ($LVallele_count{$_} > $LVmax_count) {
						$LVbest_allele=$_;
					}
				}
				my @LVarr5=split(/\//, $LVarr2[0]);
				my @LVarr6=();
				foreach (@LVarr5) {
					push (@LVarr6, $LVbest_allele);
				}
				$LVarr2[0]=join('/', @LVarr6);
			}
			my $LVjoin_arr2=join(':', @LVarr2);
			${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{'GT'}=$LVarr2[0];
			${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{$LVrg}}{'genotypes'}=$LVjoin_arr2;
			unless (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'GT'}) {
				${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'GT'}=$LVarr2[0];
			}
			unless (exists ${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'genotypes'}) {
				${${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col'.($LVx + 10)}}{'genotypes'}=$LVjoin_arr2;
			}
			$LVcolgroup[$LVx]='col'.($LVx + 10);
			last LU2 if ($LVonly10thcol==1);
		}}##LU2
#		print "${LVsubinfo}Test: readgroup list: @LVreadgroup\n"; ### For test ###
#		print "${LVsubinfo}Test: Column list: @LVcolgroup\n"; ### For test ###
		@{${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'read_groups'}}=@LVreadgroup;
		@{${${$LVret_vcf{$LVarr[0]}}{$LVarr[1]}}{'col_groups'}}=@LVcolgroup;
	}
	close VCF;
	
	#	print "vcf input\n";###test###
	if ($VcfKit_debug) {
		foreach my $LVchrom (keys %LVret_vcf) {
			my @LVarr3=keys %{$LVret_vcf{$LVchrom}};
			@LVarr3=sort {$a<=>$b} @LVarr3;
			print "${LVsubinfo}Reference: $LVchrom\n";
			print "${LVsubinfo}Number of variations on $LVchrom: ", scalar(@LVarr3), "\n";
			foreach (@LVarr3) {
				print $LVchrom, "\t", "$_\t${${$LVret_vcf{$LVchrom}}{$_}}[0]\t${${$LVret_vcf{$LVchrom}}{$_}}[1]\t${${$LVret_vcf{$LVchrom}}{$_}}[2]\n";
			}
		}
	}
	return ($VcfKit_success, \%LVret_vcf);
}



### Read variant type based on the vcf genotypes
### &ReadVariantType(Pos, RefAllele, VarAllele, GenoType, Bio::DB::Sam alignment, returncode)
### Global: 
### Dependancy: Bio::DB::Sam
### Note:returncode (1 = return allele code, 2=return allele code and quality)
sub ReadVariantType {
	my ($RVTpos, $RVTref, $RVTvar, $RVTgeno, $RVTsam_align, $RVTret_code)=@_;
	
	$RVTret_code=1 unless (defined $RVTret_code);
	my $RVTgeno_delimiter=',';
	my $RVTdebug=$VcfKit_debug;
	
	my $RVTref_seqid=$RVTsam_align->seq_id;
	my $RVTreadid=$RVTsam_align->name;
	my $RVTquery_seq=$RVTsam_align->query->dna;
	my $RVTreadstart=$RVTsam_align->start;
	my $RVTcigar=$RVTsam_align->cigar_str;
	my $RVTmdstring=$RVTsam_align->get_tag_values('MD');
	my $returnString='';
	my %RVTgenoahsh=();
#	if ($RVTgeno=~/\//) {
#		my @RVTgenoarray=split(/\//, $RVTgeno);
#		foreach (@RVTgenoarray) {
#			if ($_=~/^\d+$/) {
#				$RVTgenoahsh{$_}++;
#			}
#		}
#	}
#	else {
#		return '.';
#	}
	
#	print "\n\n\nSUB(ReadVariantType)Input: \n\tReference ID: ".$RVTref_seqid."\n\tVCF_pos: ".$RVTpos."\n\tGenotype: ".$RVTgeno."\n\tReadID: ".$RVTreadid."\n\tRead Seq: ".$RVTquery_seq."\n\tReadStart: ".$RVTreadstart."\n\tRefAllele: ".$RVTref."\n\tVarAllele: ".$RVTvar."\n\tCIGAR: ".$RVTcigar."\n\tMDstring: ".$RVTmdstring."\n";###test###
##COMMENT: Separate genotype, Assign number to each genotype: ref=0, var1=1, var2=2 ...
	my %RVTallele2geno=();
	$RVTallele2geno{$RVTref}=0;
#	print "SUB(ReadVariantType)Test: ref: $RVTref, $RVTallele2geno{$RVTref}\n";###test###
	my $RVTi=1;
	my @RVTvars=split(/$RVTgeno_delimiter/, $RVTvar);
#	print "VariantAllele: $RVTvar, SplitNumber: ".scalar(@RVTvars)."\n" if ($RVTdebug); ### For test ###
	my $RVTcat_is_deletion=0;
	foreach my $RVind_var (@RVTvars) {
		if ($RVind_var=~ /^[a-zA-Z]+$/) {
			$RVTallele2geno{$RVind_var}=$RVTi;
		}
		elsif ($RVind_var eq '.') {
			$RVTcat_is_deletion=1;
		}
		$RVTi++;
	}
#	map {print "Allele: $_, Geno: $RVTallele2geno{$_}\n" } keys %RVTallele2geno; ### For Test ###
###COMMENT: verify CIGAR length = query length
	my $RVTcigarOperations = &SplitCigar($RVTcigar);
	my $EVTcigar_cal_length=0;
	foreach (@$RVTcigarOperations){#calcular cigar length
		return 'L' unless (defined $_->[0] and defined $_->[1]);
		if ($_->[1] =~/^[MZIS=X]{1}$/) {
			$EVTcigar_cal_length+=$_->[0];
		}
	}
	unless ($EVTcigar_cal_length == length($RVTquery_seq)) {#if legnth(cigar) != length(read_seq), exit
		return 'L';
	}

#	foreach (keys %RVTallele2geno) {print "SUB(ReadVariantType)Test: ref2: $_, $RVTallele2geno{$_}\n";} ###test###
##COMMENT: Determine SNP/InDel vcf position for a read
	my $RVTrefPos = $RVTreadstart;
	my $RVTseqPos = 0;
	my $RVTthis_read_var='';
	my $EVTcapture_start_pos=0;
	my $RVTcapture_end_pos=0;
	my $RVTtest_var_start=0;
	my $RVTtest_var_end=0;
	my $RVTref_var_end=$RVTpos+length($RVTref);
	my $RVT5endloss=0;
	if ($RVTpos<$RVTreadstart) {
		$RVTpos=$RVTreadstart;
		$RVT5endloss=1;
	}
	foreach my $operation (@$RVTcigarOperations) {
		my $cig_length = $operation->[0];
		my $cig_op = $operation->[1];
		print "SUB(ReadVariantType)Test: Cigar: ($cig_length, $cig_op)\n" if ($RVTdebug);
		my ($RVTadd_ref, $RVTadd_seq)=(0, 0);
		my $RVTcigar_type=0;
		if ($cig_op =~ /^D$/) {
			$RVTadd_ref=$cig_length;
			$RVTcigar_type=3;
#			print "Deletion\n";###for test###
		}
		elsif ($cig_op =~ /^I$/) {
			$RVTadd_seq=$cig_length;
			$RVTcigar_type=2;
#			print "Insertion\n";###for test###
		}
		elsif ($cig_op =~ /^[M=]{1}$/) {
			$RVTadd_ref= $cig_length;
			$RVTadd_seq= $cig_length;
			$RVTcigar_type=1;
#			print "match\n";###for test###
		}
		elsif ($cig_op =~ /^S$/) {
			$RVTadd_seq=$cig_length;
			$RVTcigar_type=4;
		}
		elsif ($cig_op =~ /^N$/) {
			$RVTadd_ref= $cig_length;
			$RVTcigar_type=5;
		}
		elsif ($cig_op =~ /^P$/) {
			$RVTcigar_type=6;
		}
		elsif ($cig_op =~ /^H$/) {
			$RVTcigar_type=7;
		}
		else {
			print STDERR "SUB(ReadVariantType)Error: unrecognized Cigar State: $cig_op (Cigar: ". $RVTsam_align->cigar_str. " at position $RVTpos of Reference sequence ".$RVTsam_align->seq_id."\n";
			return '?';
		}
		$EVTcapture_start_pos=$RVTseqPos if ($RVTtest_var_start==0);
		$RVTcapture_end_pos=$RVTseqPos;
		if (($RVTrefPos <= $RVTpos) and ($RVTpos <($RVTrefPos+$RVTadd_ref))) {
			if (($RVTrefPos+$RVTadd_ref)<$RVTref_var_end) {
				if ($RVTtest_var_start==0) {
					if ($RVTtest_var_start==0 and $RVTcigar_type==1) {
						$EVTcapture_start_pos+=($RVTpos-$RVTrefPos);
						$RVTcapture_end_pos+=$RVTadd_seq;
						$RVTtest_var_start=1;
#						print "Starts1: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
					}
					else {
						return 'unknown2';
					}
				}
				else {
					if ($RVTcigar_type=~/^[23]{1}$/) {
						$RVTcapture_end_pos+=$RVTadd_seq;
#						print "Starts2: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
					}
					else {
						return 'unknown3';
					}
				}
			}
			else {
				if ($RVTtest_var_start==0 and $RVTcigar_type==1) {
					$EVTcapture_start_pos+=$RVTpos-$RVTrefPos;
					$RVTcapture_end_pos=$EVTcapture_start_pos+$RVTref_var_end-$RVTpos;
					$RVTtest_var_end=1;
#					print "Starts3: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
				}
				else {
					return 'Unknown1';
				}
			}
		}
		elsif ($RVTrefPos > $RVTpos) {
			if ($RVTcigar_type=~/^[4567]{1}$/) {
				return 'Nreads';
			}
			elsif (($RVTrefPos+$RVTadd_ref)<=$RVTref_var_end) {#insertion&deletion
				 $RVTcapture_end_pos+=$RVTadd_seq;
#				 print "Starts4: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
			}
			elsif ($RVTref_var_end>$RVTrefPos and $RVTref_var_end<($RVTrefPos+$RVTadd_ref)) {
					$RVTcapture_end_pos+=$RVTref_var_end-$RVTrefPos;
					$RVTtest_var_end=1;
#					print "Starts5: $EVTcapture_start_pos - $RVTcapture_end_pos\n";###for test###
			}
		}
		last if ($RVTtest_var_end==1);
		$RVTrefPos+=$RVTadd_ref;
		$RVTseqPos+=$RVTadd_seq;
	}
	print "SUB(ReadVariantType)Test: Read pos: $RVTseqPos\n" if ($RVTdebug);
	print "SUB(ReadVariantType)Test: ReadLength: ".length($RVTquery_seq)."\tSubstr(Start - End): ".$EVTcapture_start_pos." - ".$RVTcapture_end_pos."\n" if ($RVTdebug);###For test###

	if ($RVTcapture_end_pos>length($RVTquery_seq)) {
		print STDERR "SUB(ReadVariantType)Test: Readlength: ".length($RVTquery_seq)."\tSubstrEnd: $RVTcapture_end_pos\n";
		print STDERR "SUB(ReadVariantType)Warnings: use read ($RVTreadid) length at $RVTpos (Ref:Var:Gonotype=$RVTref : $RVTvar : $RVTgeno) of $RVTref_seqid instead\n";###For test###
		$RVTcapture_end_pos=length($RVTquery_seq);
	}
	unless ($EVTcapture_start_pos < length($RVTquery_seq) and $EVTcapture_start_pos < $RVTcapture_end_pos) {
		print STDERR "SUB(ReadVariantType)Error: Unknown2\n";###For test###
		return "unknown";
	}
	my $RVTcapture_string='';
	if (($RVTcapture_end_pos-$EVTcapture_start_pos) >0) {
		$RVTcapture_string=substr($RVTquery_seq, $EVTcapture_start_pos, $RVTcapture_end_pos-$EVTcapture_start_pos);
	}
	elsif (($RVTcapture_end_pos-$EVTcapture_start_pos) ==0 and $RVT5endloss==1) {
		$RVTcapture_string='.';
	}
	if ($RVTcapture_string eq '') {
		print STDERR "SUB(ReadVariantType)Error: substr failed\n";###For test###
		return "unknown2";
	}
#	print "SUB(ReadVariantType)Test: Captured String: $RVTcapture_string\n";# if ($RVTdebug); ### For test ###
	my @returnarr=();
	my $RVTreturn_type;
	EVTBLOCK1: {if ($RVTtest_var_end) {
		if ($RVT5endloss==0) {
			if (exists $RVTallele2geno{$RVTcapture_string}) {
				$RVTreturn_type=$RVTallele2geno{$RVTcapture_string};
				last EVTBLOCK1;
			}
			else {
				return 'NotExist';
			}
		}
		else {
			foreach my $RVTind_allele (keys %RVTallele2geno) {
				if ($RVTind_allele=~/$RVTcapture_string$/) {
					push (@returnarr, $RVTallele2geno{$RVTind_allele});
				}
			}
		}
	}
	else {
		foreach my $RVTind_allele (keys %RVTallele2geno) {
			if ($RVT5endloss==0) {
				if ($RVTind_allele=~/^$RVTcapture_string/) {
					push (@returnarr, $RVTallele2geno{$RVTind_allele});
				}
			}
			else {
				if ($RVTind_allele=~/$RVTcapture_string/) {
					push (@returnarr, $RVTallele2geno{$RVTind_allele});
				}
			}
		}
	}
#	print "@returnarr\n";### For test ###
	
	if (scalar(@returnarr)==1) {
		$RVTreturn_type= $returnarr[0];
	}
	else {
		return 'Multiple';
	}
#	print "ReturnType: $RVTreturn_type\n";### For test ###
	}###EVTBLOCK1
#	unless (exists $RVTgenoahsh{$RVTreturn_type}) {
#		return 'NOTallelic';
#	}
	if ($RVTret_code ==1) {
		return $RVTreturn_type;
	}
	elsif ($RVTret_code ==2) {
#		print "ReturnCode: 2\n"; ### For test ###
		my @RVTqscore=$RVTsam_align->query->qscore;
#		print "QScore: @RVTqscore\n";### For test ###
		my @RVTqsstr=map {chr($_+33)} @RVTqscore;
#		print "ReadID: $RVTreadid\nReadSeq: $RVTquery_seq\nReadQual: ".join('',@RVTqsstr)."\n";### For test ###
		
		my $RVTreturn_score='D';
		my $EVTqs_string='';
#		print "Start-End: $EVTcapture_start_pos-$RVTcapture_end_pos\n";### For test ###
		if (($RVTcapture_end_pos-$EVTcapture_start_pos) >0) {
			my $RVTcount=0;
			my $RVTsum=0;
#			print "Capture quality\n";### For test ###
#			print "Start-end: $EVTcapture_start_pos - $RVTcapture_end_pos\n";### For test ###
			for (my $RVTi=$EVTcapture_start_pos; $RVTi<$RVTcapture_end_pos; $RVTi++) {
				if (0) {### For test ###
					$EVTqs_string.=chr($RVTqscore[$RVTi]+33);
#					print "Capture quality: $EVTqs_string\n";
				}
				$RVTsum+=$RVTqscore[$RVTi];
				$RVTcount++;
#				print "Converting $RVTi\n";### For test ###
			}
			my $RVTreturn_value=int($RVTsum/$RVTcount);
			$RVTreturn_score=chr($RVTreturn_value+33);
			if (0) {### For test ###
				print "QualityString: $EVTqs_string\nAverage_value: $RVTreturn_value\nAverageCode: $RVTreturn_score\n";
			}
		}
		return ($RVTreturn_type, $RVTreturn_score);
	}
	else {
		##to be defined
	}
}



### RunHaptree
### &RunHaptree (AABBDDvcf, phasedVCF_name, HapTreeReadsInput_name, fixgenotype, reads, [path_haptree], [rundir])
### Global: 
### Dependency: FileKit::DeletePath
### Note: chdir back for sure
sub RunHaptree {
	my ($RHvcffilein, $RHvcffileout, $RHreadout, $RHfixgeno_hashindex, $RHreads_hashindex, $RHpath_haptree)=@_;
	
	local *RHVCFIN; local *RHVCFOUT; local *HAPTREEREADS;
	my $RHsubinfo='SUB(VcfKit::RunHaptree)';
	my $RHtestcmd=0;
	my $RHcurdir=getcwd;
	$RHpath_haptree='haptree' unless (defined $RHpath_haptree);
	
	unless (defined $RHvcffilein and -s $RHvcffilein) {
		print STDERR "${RHsubinfo}Error: VCF file not found: $RHvcffilein\n";
		return $VcfKit_failure;
	}
	unless (defined $RHvcffileout) {
		print STDERR "${RHsubinfo}Error: invalid VCF output for HapTree\n";
		return $VcfKit_failure;
	}
#	print "SUB(RunHaptree)Test: vcffile for HapTree: $RHvcffileout\n"; ### For test ###
	unlink $RHvcffileout if (-s $RHvcffileout);
	unless (defined $RHreadout) {
		print STDERR "${RHsubinfo}Error: invalid Read output for HapTree\n";
		return $VcfKit_failure;
	}
	unlink $RHreadout if (-s $RHreadout);
##Format: $RHfixgeno{chr}->{pos}='0/1/0'; ('A/B/D')
	my %RHfixgeno=%{$RHfixgeno_hashindex} if (defined $RHfixgeno_hashindex); undef $RHfixgeno_hashindex;
	if (0) {### For test ###
		foreach my $RHind_chr (keys %RHfixgeno) {
#			print "${RHsubinfo}Test: Chrom: $RHind_chr\n";
			foreach my $RHind_pos (keys %{$RHfixgeno{$RHind_chr}}) {
				print "${RHsubinfo}Test: Chrom: $RHind_chr Pos: $RHind_pos, Allele: ${$RHfixgeno{$RHind_chr}}{$RHind_pos}\n";
			}
		}
	}
##FORMAT: %RHreads=(chr => (pos => (allele1 => @readIDs, allele2 => @readids)))
	my %RHreads=%{$RHreads_hashindex} if (defined $RHreads_hashindex); undef $RHreads_hashindex;
	if (0) {### For test ###
		foreach my $RHchrom (keys %RHreads) {
			foreach my $RHind_pos (keys %{$RHreads{$RHchrom}}) {
				foreach my $RHalleles (keys %{${$RHreads{$RHchrom}}{$RHind_pos}}){
					print "${RHsubinfo}Test: Chrom: $RHchrom Pos: $RHind_pos Allele: $RHalleles NumberofReads: ".scalar(@{${${$RHreads{$RHchrom}}{$RHind_pos}}{$RHalleles}})." \n";
				}
			}		
		
		}
	}
	my $RHoutdir=RetrieveDir($RHvcffileout);
	unless (chdir $RHoutdir) {
		print STDERR "${RHsubinfo}Error: can not chdir: $RHoutdir\n";
		return $VcfKit_failure;
	}
	close RHVCFIN if (defined fileno(RHVCFIN));
	unless (open RHVCFIN, "<", $RHvcffilein) {
		print STDERR "${RHsubinfo}Error: can not open: $RHvcffilein\n";
		return $VcfKit_failure;
	}
	close RHVCFOUT if (defined fileno (RHVCFOUT));
	unless (open RHVCFOUT, ">", $RHvcffileout) {
		print STDERR "${RHsubinfo}Error: can not write: $RHvcffileout\n";
		return $VcfKit_failure;
	}
	close HAPTREEREADS if (defined fileno(HAPTREEREADS));
	unless (open HAPTREEREADS, ">", $RHreadout) {
		print STDERR "${RHsubinfo}Error: can not write: $RHreadout\n";
		return $VcfKit_failure;
	}
	my $RHvcfoutlinenum=0;
##FORMAT: $RHhaptree_reads{readid}=('$RHvcfoutlinenum: allele', '$RHvcfoutlinenum: allele', ...)
	my %RHhaptree_reads=();
	RHBLOCK1: {while (my $RHline=<RHVCFIN>) {
		next if ($RHline=~/^#/);
		chomp $RHline;
		my @RHarr1=split(/\t/, $RHline);
		#chrom: $RHarr1[0]
		#Pos:	$RHarr1[1]
		#geno:	$RHarr1[9]
#		print "${RHsubinfo}Test: Chrom: $RHarr1[0] Pos: $RHarr1[1]\n"; ### For test ###
		unless (exists ${$RHfixgeno{$RHarr1[0]}}{$RHarr1[1]} and ${$RHfixgeno{$RHarr1[0]}}{$RHarr1[1]}=~/\//) {
			print STDERR "${RHsubinfo}Error: fixed genotype not exists at chr:pos $RHarr1[0]: $RHarr1[1]\n";
			$RHtestcmd=1;
			last RHBLOCK1;
		}
		my %RHgenoatpos=();
		my @RHarr2=split(/\//, ${$RHfixgeno{$RHarr1[0]}}{$RHarr1[1]});
		foreach (@RHarr2) {
			$RHgenoatpos{$_}++;
		}
		my @RHarr3=keys %RHgenoatpos; undef %RHgenoatpos;
		if (scalar(@RHarr3)<1) {
			print STDERR "${RHsubinfo}Error: spliting fixed genotype error at chr:pos $RHarr1[0]: $RHarr1[1]\n";
			$RHtestcmd=1;
			last RHBLOCK1;
		}
		my @RHarr4=split(/:/, $RHarr1[9]);
		$RHarr4[0]=${$RHfixgeno{$RHarr1[0]}}{$RHarr1[1]};
		$RHarr1[9]=join(':', @RHarr4);
		$RHarr1[2]=$RHarr1[0].'_'.$RHarr1[1];
		$RHarr1[5]='.';
		$RHarr1[6]='PASS';
		$RHarr1[7]='*';
		print RHVCFOUT join("\t", @RHarr1), "\n";
		foreach my $RHind_allele (@RHarr3) {
			unless (exists ${${$RHreads{$RHarr1[0]}}{$RHarr1[1]}}{$RHind_allele}) {
				print STDERR "${RHsubinfo}Error: allele error at chr:pos:allele $RHarr1[0]: $RHarr1[1] : $RHind_allele\n";
				$RHtestcmd=1;
				last RHBLOCK1;
			}
			foreach my $RHindread (@{${${$RHreads{$RHarr1[0]}}{$RHarr1[1]}}{$RHind_allele}}) {
				push (@{$RHhaptree_reads{$RHindread}}, "$RHvcfoutlinenum: $RHind_allele");
			}
		}
		$RHvcfoutlinenum++;
	}}#RHBLOCK1
	close RHVCFIN;
	close RHVCFOUT;
	return $VcfKit_failure if ($RHtestcmd);
	foreach my $RHind_read2 (keys %RHhaptree_reads) {
		print HAPTREEREADS '{', join(', ', @{$RHhaptree_reads{$RHind_read2}}), '}', "\n";
	}
	close HAPTREEREADS;
	unless (-s $RHvcffileout and -s $RHreadout) {
		print STDERR "${RHsubinfo}Error: HapTree input not exists\n";
		return $VcfKit_failure;
	}
	my $RHhaptree_outdir=$RHoutdir.'/haptreeout';
	DeletePath($RHhaptree_outdir) if (-d $RHhaptree_outdir);###remove haptree outdir as required
	my $RHcmd=$RHpath_haptree." $RHreadout $RHvcffileout $RHhaptree_outdir";
	if (! &exec_cmd_return($RHcmd)) {
		print STDERR "${RHsubinfo}Error: HapTree running error\n";
		return $VcfKit_failure;
	}
	unless (-s "$RHhaptree_outdir/HapTreeSolution") {
		print STDERR "${RHsubinfo}Error: HapTree output error\n";
		return $VcfKit_failure;
	}
	chdir $RHcurdir;
	return ($VcfKit_success, "$RHhaptree_outdir/HapTreeSolution");
}



### Read haptree phased alleles
### ReadHaptreeOut()
### Global: 
### Dependency: 
### Note: 
sub ReadHaptreeOut {
	my ($RHOphasefile, $RHOgeno_ploidy, $RHOfixalleles_hashindex)=@_;
	
	local *PHASED;
	my $RHOsubinfo='SUB(ReadHaptreeOut)';
	my $RHOtest_cmd=0;
	unless (-s $RHOphasefile) {
		print STDERR "${RHOsubinfo}Error: invalid HapTree output: $RHOphasefile\n";
		return 1;
	}
### Initialising ploidy
	my $RHOnum_ploidy=length($RHOgeno_ploidy);
	my ($RHOaa_expressed, $RHObb_expressed, $RHOdd_expressed)=(0, 0, 0);
	my @RHOgenomes=();
	if ($RHOgeno_ploidy=~/A/) {
		$RHOaa_expressed=1;
		push (@RHOgenomes, 'A');
	}
	if ($RHOgeno_ploidy=~/B/) {
		$RHObb_expressed=1;
		push (@RHOgenomes, 'B');
	}
	if ($RHOgeno_ploidy=~/D/) {
		$RHOdd_expressed=1;
		push (@RHOgenomes, 'D');
	}
	if ($RHOnum_ploidy != $RHOaa_expressed + $RHObb_expressed + $RHOdd_expressed or ($RHOnum_ploidy <2 or $RHOnum_ploidy > 3)) {
		print STDERR "${RHOsubinfo}Error: ploidy (=$RHOnum_ploidy) error\n";
		return 1;
	}
##FORMAT: %RHOfixalleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	my %RHOfixalleles=%{$RHOfixalleles_hashindex}; undef $RHOfixalleles_hashindex;
	close PHASED if (defined fileno(PHASED));
	unless (open PHASED, "<", $RHOphasefile) {
		print STDERR "${RHOsubinfo}Error: can not $RHOphasefile\n";
		return 1;
	}
##FORMAT: file format
#BLOCK 	1000	1	1	1
#snp1	1000	0	0	1	
#snp2	1028	0	0	1	
#snp3	1118	1	1	0	
#snp4	1143	0	1	1	
#snp5	1320	0	0	1	
#snp6	1497	1	1	0	
#snp7	1518	1	0	0	
#snp8	1545	0	1	1	
#snp9	1677	1	1	0	
#snp10	1769	1	0	1	
#snp11	1795	0	1	1
	my $RHOblocknum=0;
	my %RHOtemp_block=();
	my %RHOfinal_fixallele=();
	while (my $RHOline=<PHASED>) {
		if ($RHOline=~/^BLOCK/) {
			$RHOblocknum++;
			next;
		}
		chomp $RHOline;
		my @RHOarr1=split(/\t/, $RHOline);
		next if (scalar(@RHOarr1) != ($RHOnum_ploidy+2));
		my ($RHOchrom, $RHOposit)=('', '');
		if ($RHOarr1[0]=~/^([\w+])_(\d+)$/) {
			$RHOchrom=$1;
			$RHOposit=$2;
		}
		else {
			print STDERR "${RHOsubinfo}Error: can not $RHOphasefile\n";
			return 1;
		}
		for (my $RHOi=2; $RHOi<scalar(@RHOarr1); $RHOi++) {
			${${${$RHOtemp_block{$RHOblocknum}}{$RHOchrom}}{$RHOposit}}{$RHOi-1}=$RHOarr1[$RHOi];
		}
	}
	close PHASED;
### count each group to each sungenome in each block
##FORMAT: %RHOtemp_assign=( block => (colnum => (A/B/D => total+count)));
	my %RHOtemp_assign=();
	my %RHOcorrelation=();
	foreach my $RHOblock (sort {$a <=> $b} (keys %RHOtemp_block)) {
		my %RHOfixarr=();
		my %RHOphase=();
		foreach my $RHOchrom2 (keys %{$RHOtemp_block{$RHOblock}}) {
			foreach my $RHOposit2 (keys %{${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}) {
				for (my $RHOi2=1; $RHOi2<=$RHOnum_ploidy; $RHOi2++) {
					push (@{$RHOphase{$RHOi2}}, ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
					if ($RHOaa_expressed) {
						if (exists ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'A'}) {
							${${$RHOtemp_assign{$RHOblock}}{$RHOi2}}{'A'}++ if (${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'A'} eq ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
							push (@{$RHOfixarr{'A'}}, ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'A'});
						}
						else {
							print STDERR "${RHOsubinfo}Error: Final AA assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
							return 1;
						}
					}
					if ($RHObb_expressed) {
						if (exists ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'B'}) {
							${${$RHOtemp_assign{$RHOblock}}{$RHOi2}}{'B'}++ if (${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'B'} eq ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
							push (@{$RHOfixarr{'B'}}, ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'B'});
						}
						else {
							print STDERR "${RHOsubinfo}Error: Final BB assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
							return 1;
						}
					}
					if ($RHOdd_expressed) {
						if (exists ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'D'}) {
							${${$RHOtemp_assign{$RHOblock}}{$RHOi2}}{'D'}++ if (${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'D'} eq ${${$RHOtemp_block{$RHOblock}}{$RHOchrom2}}{$RHOi2});
							push (@{$RHOfixarr{'D'}}, ${${$RHOfixalleles{$RHOchrom2}}{$RHOposit2}}{'D'});
						}
						else {
							print STDERR "${RHOsubinfo}Error: Final DD assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
							return 1;
						}
					}
				}
			}
		}
### Correlation
		foreach my $RHOind_fixallele (keys %RHOfixarr) {
			foreach my $RHOind_phasing (keys %RHOphase) {
				my ($RHOcorr_retcode, $RHOcov, $RHOcor, $RHOnum)=&CallCorrelation(\@{$RHOfixarr{$RHOind_fixallele}}, \@{$RHOphase{$RHOind_phasing}});
				unless ($RHOcorr_retcode) {
					${${$RHOcorrelation{$RHOblock}}{$RHOind_phasing}}{$RHOind_fixallele}=$RHOcor;
				}
			}
		}
	}
### assign each group to each sungenome in each block
	my $RHOtestcmd_count=0;
	RHOBLOCK1: {foreach my $RHOind_block2 (sort {$a <=> $b} (keys %RHOtemp_assign)) {
		my $test_corelation=0;
		my %RHOgroup2genome=();
		my %RHOgenome2group=();
		RHOBLOCK2: {foreach my $RHOind_group (sort {$a <=> $b} (keys %{$RHOtemp_assign{$RHOind_block2}})) {
			my $RHObest_group='unknown';
			my $RHOmax_count=0;
			my $RHOnum_max=0;
			foreach my $RHOgenotype (@RHOgenomes) {
				if (exists ${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype}){
					if (${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype} > $RHOmax_count) {
						$RHOnum_max=1;
						$RHOmax_count=${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype};
						$RHObest_group=$RHOgenotype;
					}
					elsif (${${$RHOtemp_assign{$RHOind_block2}}{$RHOind_group}}{$RHOgenotype} == $RHOmax_count) {
						$RHOnum_max++;
						$RHObest_group.="$RHOgenotype";
					}
				}
			}
			if ($RHOnum_max !=1 or $RHObest_group !~ m/^[ABD]{1}$/) {
				$test_corelation=1;
			}
			$RHOgroup2genome{$RHOind_group}=$RHObest_group;
			if (exists $RHOgenome2group{$RHObest_group}) {
				
				$test_corelation=1;
			}
			else {
				$RHOgenome2group{$RHObest_group}=$RHOind_group;
			}
		}}###RHOBLOCK2
		
	}}###RHOBLOCK1
}



### RunHapCompass
### &RunHapCompass(bam, fragments, vcf, ploidy, out.vcf.prefix, final.phased.vcf, additional_cmds, path_java, path_hapcompass)
### Global:
### Require: hapcompass.jar hc2vcf.jar
### Dependancy: 
### Note:
sub RunHapCompass {
	my ($RHCbam, $RHCfragments, $RHCvcf, $RHCploidy, $RHCoutput_vcf_prefix, $RHCvcfout, $RHCadd_cmds, $RHpath_java, $RHpath_hapcompass)=@_;
	
	my $RHsubinfo='SUB(VcfKit::RunHapCompass)';
	$RHCadd_cmds=' ' unless (defined $RHCadd_cmds);
	$RHpath_java='java' unless (defined $RHpath_java);
	my $RHCcmd_hapcompass='';
	
##COMMIT: check input files for HapCompass
	unless (defined $RHpath_hapcompass and -s $RHpath_hapcompass) {
		print STDERR "${RHsubinfo}Error: HapCompass path error: \n";
		return $VcfKit_failure;
	}
	$RHCcmd_hapcompass="$RHpath_java -jar $RHpath_hapcompass";
	if (defined $RHCfragments and $RHCfragments ne '0' and -s $RHCfragments) {
		$RHCcmd_hapcompass.=" --fragment $RHCfragments ";
	}
	elsif (defined and $RHCbam ne '0' and -s $RHCbam) {
		$RHCcmd_hapcompass.=" --bam $RHCbam ";
	}
	else {
		print STDERR "${RHsubinfo}Error: invalid BAM or fragments file\n";
		return $VcfKit_failure;
	}
	unless (defined $RHCvcf and $RHCvcf ne '0' and -s $RHCvcf) {
		print STDERR "${RHsubinfo}Error: invalid vcf input file\n";
		return $VcfKit_failure;
	}
	$RHCcmd_hapcompass.=" --vcf $RHCvcf ";
	unless ($RHCploidy=~/^\d+$/ and $RHCploidy>=2) {
		print STDERR "${RHsubinfo}Error: HapCompass ploidy error\n";
		return $VcfKit_failure;
	}
	$RHCcmd_hapcompass.=" --ploidy $RHCploidy ";
	if (! defined $RHCoutput_vcf_prefix or $RHCoutput_vcf_prefix eq '') {
		print STDERR "${RHsubinfo}Error: HapCompass output prefix error\n";
		return $VcfKit_failure;
	}
	elsif ($RHCoutput_vcf_prefix !~ /^\//) {
		print STDERR "${RHsubinfo}Error: HapCompass output prefix should be a full path\n";
		return $VcfKit_failure;
	}
	$RHCcmd_hapcompass.=" --output $RHCoutput_vcf_prefix ";
	unless (defined $RHCvcfout) {
		print STDERR "${RHsubinfo}Error: HapCompass output NOT specified\n";
		return $VcfKit_failure;
	}
	unlink $RHCvcfout if (-e $RHCvcfout);
#	$RHCcmd_hapcompass.=" --iterations 10 ";
##COMMIT: Run HapCompass, output: $RHCoutput_vcf_prefix
##									_frags.txt
##									_MWER_solution.txt
##									_phasedSolution.txt
##									_reads.sam
##									_reduced_representation.sam
##									_reduced_representation.vcf
	if (! &exec_cmd_return("$RHCcmd_hapcompass >> hapcompass.log 2>&1")) {
		print STDERR "${RHsubinfo}Error: HapCompass running error\n";
		return $VcfKit_failure;
	}
##COMMIT: check HapCompass output
	if (! -s $RHCoutput_vcf_prefix.'_MWER_solution.txt') {
		print STDERR "${RHsubinfo}Reference: HapCompass output *_MWER_solution.txt error\n";
		return $VcfKit_failure;
	}
	unless (MoveFile("${RHCoutput_vcf_prefix}_MWER_solution.txt", $RHCvcfout)) {
		print STDERR "${RHsubinfo}Error: moving files: ${RHCoutput_vcf_prefix}_MWER_solution.txt to $RHCvcfout\n";
		return $VcfKit_failure;
	}
	return $VcfKit_success;
}



### Comvert Hapcompass output to VCF
### HapcompassOut2Vcf(vcfin, vcfout, hcout, path_hc2vcfjar ('java -jar path_hc2vcf'))
### Global: $VcfKit_failure; $VcfKit_success;
### Dependency: hapcompass.jar hc2vcf.jar
### Note:
### Return: 0=failure; 1= success
sub HapcompassOut2Vcf {
	my ($HOVvcfin, $HOVploidy, $HOVvcfout, $HOVhcout, $HOVhc2vcfjar)=@_;
	
	my $HOVsubinfo='SUB(VcfKit::HapcompassOut2Vcf)';
	
	unless (defined $HOVvcfin, and -s $HOVvcfin) {
		print STDERR "${HOVsubinfo}Error: invalid VCF input\n";
		return $VcfKit_failure;
	}
	unless (defined $HOVvcfout and $HOVvcfout=~/^\S+/) {
		print STDERR "${HOVsubinfo}Error: invalid VCF output\n";
		return $VcfKit_failure;
	}
	unlink $HOVvcfout if (-e $HOVvcfout);
	unless (defined $HOVhcout and -s $HOVhcout) {
		print STDERR "${HOVsubinfo}Error: invalid HapCompass out input\n";
		return $VcfKit_failure;
	}
	unless (defined $HOVploidy and $HOVploidy=~/^\d+$/) {
		print STDERR "${HOVsubinfo}Error: invalid ploidy number\n";
		return $VcfKit_failure;
	}
	##COMMIT: run hc2vcf to convert hapcompass output to vcf
##			output $RHCoutput_vcf_prefix.'_MWER_solution.txt.vcf'
	unless (&exec_cmd_return("$HOVhc2vcfjar $HOVhcout $HOVvcfin $HOVploidy")) {
		print STDERR "${HOVsubinfo}Reference: hc2vcf running error\n";
		return $VcfKit_failure;
	}
	unless (-s "$HOVhcout.vcf") {
		print STDERR "${HOVsubinfo}Error: hc2vcf output error\n";
		return $VcfKit_failure;
	}
	unless (MoveFile("$HOVhcout.vcf", $HOVvcfout)) {
		print STDERR "${HOVsubinfo}Error: hc2vcf output error\n";
		return $VcfKit_failure;
	}
	return $VcfKit_success;
}



### Prepare Hapcompass Fragments
### HcFragments($readid, fragments_out, $HFreadsum)
### Global: 
### Dependency:
### Note:
sub HcFragments {
	my ($HFid_index, $HFfragout, $HFreadsum)=@_;
#$HFid_index is hash reference

	local *FRAGMENTS; local *FRAGMENTPOS;
#Format: %GVreadfragments=(readid => (chrom => (pos => (geno1 => quality, geno2 => quality)));
	my $HFsubinfo='SUB(VcfKit::HcFragments)';
	my %HFallele3sites=();

	unless (defined $HFfragout) {
		print STDERR $HFsubinfo, "Error: undefined gragmentsout\n";
		return $VcfKit_failure;
	}
	
	close FRAGMENTS if (defined fileno(FRAGMENTS));
	close FRAGMENTPOS if (defined fileno(FRAGMENTPOS));
	unlink "$HFfragout" if (-e "$HFfragout");
##Format: %HFreadparts=($read => $chr => $pos =>$strand => $allele => qual)
	my %HFreadparts=();
##Format: %HFcorrectedsites=($chr => $pos => $allele => ++)
	my %HFcorrectedsites=();
##Format: %excludedreads=(read => chr:pos)
	my %excludedreads=();
##Format: %HFreadpart2=($readid => chrom => $partstart => $pos => $allele => $qual);
	my %HFreadpart2=();
##Format: %HFsitereferencer=('pos2ref' => chr => pos => VcfReferenceNumber)
	my %HFsitereferencer=();
	
	unless (open (FRAGMENTS, ">$HFfragout")) {
		print STDERR $HFsubinfo, "Error: can not write to fragments file: $HFfragout\n";
		return $VcfKit_failure;
	}
	unless (open (FRAGMENTPOS, "> $HFfragout.pos")) {
		print STDERR $HFsubinfo, "Error: can not write to POS fragments file: $HFfragout.pos\n";
		return $VcfKit_failure;
	}
	my @HFrealids=keys %{$HFid_index};
	if (scalar(@HFrealids)==0) {
		print $HFsubinfo, "Error: no fragments\n";
		return ($VcfKit_success, \%HFreadpart2, \%HFcorrectedsites);
	}
	foreach my $HFreadid (@HFrealids) {
		my %HFexcluded_region=();
		foreach my $HFchrom (keys %{${$HFid_index}{$HFreadid}}) {
			foreach my $HFpos (sort {$a <=>$b} (keys %{${$HFid_index}{$HFreadid}{$HFchrom}})) {
				my %HFallele=();
				my $HFexcludedatpos=0;
				foreach my $HFstrand (sort keys %{${$HFid_index}{$HFreadid}{$HFchrom}{$HFpos}}) {
					my @HFalles_at_pos=keys %{${$HFid_index}{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}};
#					print STDERR "allele: @HFalles_at_pos\n"; ### For test ###
					foreach my $HFstrandallele (@HFalles_at_pos) {
						if ($HFstrandallele=~/^\d+$/) {
							$HFallele{$HFstrandallele}=${$HFid_index}{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}{$HFstrandallele};
							$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}{$HFstrandallele}=$HFallele{$HFstrandallele};
						}
						else {
							$HFexcludedatpos=1;
							last;
						}
					}
				}
				unless (scalar(keys %HFallele)==1) {
#					print STDERR "${HFsubinfo}Warnings: Read ($HFreadid) has two alleles at Chr:Pos ($HFchrom:$HFpos), ignoring\n";
					$HFexcludedatpos=1;
				}
				if ($HFexcludedatpos==1) {
					foreach my $HFpos3 (sort {$a <=> $b} keys %{${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}}) {
						my $HFendposition=${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}{$HFpos3};
						if ($HFpos>=$HFpos3 and $HFendposition=~/^\d+$/ and $HFpos <= $HFendposition) {
							$HFexcluded_region{$HFchrom}{$HFpos3}=${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}{$HFpos3};
							${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}{$HFpos3}.='N';
						}
					}
				}
			}
		}
		foreach my $HFchrom (keys %{$HFreadparts{$HFreadid}}) {
			my @HFposin=sort keys %{$HFreadparts{$HFreadid}{$HFchrom}};
			foreach my $HFpos (@HFposin) {
				my $HFtest_excluded=0;
				if (exists $HFexcluded_region{$HFchrom}) {
					foreach my $HFpos2 (keys %{$HFexcluded_region{$HFchrom}}) {
						if ($HFpos>=$HFpos2 and $HFpos<=$HFexcluded_region{$HFchrom}{$HFpos2}) {
							$HFtest_excluded=1;
						}
					}
				}
				if ($HFtest_excluded==0){
					foreach my $HFstrand (keys %{$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}}) {
						my @HFallelearr=keys %{$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}};
						my $HFallele=shift @HFallelearr;
						$HFcorrectedsites{$HFchrom}{$HFpos}{$HFallele}++;
					}
				}
				else {
					delete ${$HFreadparts{$HFreadid}{$HFchrom}}{$HFpos};
				}
			}
			if (scalar(keys %{$HFreadparts{$HFreadid}{$HFchrom}})==0) {
				delete ${$HFreadparts{$HFreadid}}{$HFchrom};
			} 
		}
		if (scalar(keys %{$HFreadparts{$HFreadid}}<1)) {
			${$HFreadsum}{$HFreadid}{'excluded'}++;
			delete $HFreadparts{$HFreadid};
		}
	}
	
	if (0) { ### Test %HFreadparts
		print "${HFsubinfo}Test: ".'%HFreadparts'."\n";
		print Dumper \%HFreadparts;
		print "\n";
	}
###Comment: re-defined polymorphic sites
#### Test %HFreadparts
#		print "${HFsubinfo}Test: \%HFcorrectedsites before cleaning\n"; print Dumper \%HFcorrectedsites;print "\n";
	
	foreach my $HFchrom (sort keys %HFcorrectedsites) {
		my $HFreferencer=0;;
		my @HFposits=keys %{$HFcorrectedsites{$HFchrom}};
		foreach my $HFpos (sort {$a <=> $b} @HFposits) {
			my @HFallele4=keys %{$HFcorrectedsites{$HFchrom}{$HFpos}};
			if (scalar(@HFallele4)<2) {
#				print "${HFsubinfo}Test: dropped sites: $HFchrom:$HFpos\n";### for test ###
				delete ${$HFcorrectedsites{$HFchrom}}{$HFpos};
			}
			elsif (scalar(@HFallele4)==2) {
#				print "${HFsubinfo}Test: ==2 alleles for HAPCOMPASS at chr:pos $HFchrom:$HFpos\n";### for test ###
				my $HGtest5=0;
				foreach my $HFallele5 (@HFallele4) {
					if ($HFcorrectedsites{$HFchrom}{$HFpos}{$HFallele5} <3) {
						$HGtest5=1;
						last;
					}
				}
				if ($HGtest5==0) {
					$HFsitereferencer{'pos2ref'}{$HFchrom}{$HFpos}=$HFreferencer;
					$HFreferencer++;
				}
				else {
					delete ${$HFcorrectedsites{$HFchrom}}{$HFpos};
				}
			}
			else {
#				print "${HFsubinfo}Test: more than 2 alleles for HAPCOMPASS at chr:pos $HFchrom:$HFpos\n";### for test ###
#				return $VcfKit_failure;
				$HFallele3sites{$HFchrom}{$HFpos}++;
				delete ${$HFcorrectedsites{$HFchrom}}{$HFpos};
			}
		}
		if (scalar(keys %{$HFcorrectedsites{$HFchrom}})<1) {
			delete $HFcorrectedsites{$HFchrom};
		}
	}


#	print "${HFsubinfo}Test: \%HFcorrectedsites after cleaning\n"; print Dumper \%HFcorrectedsites; print "\n"; ### For test ###
#	print "${HFsubinfo}Test: \%HFvcfreferencer\n"; print Dumper \%HFsitereferencer; print "\n"; ### For test ###

###Comment: cleaning read parts that do not contain corrected sites
	@HFrealids=keys %HFreadparts;
	foreach my $HFreadid (@HFrealids) {
		my $HFnumparts=0;
		my $HFerrorcode=0;
		my %HFincluded_region=();
		foreach my $HFchrom (keys %{${$HFreadsum}{$HFreadid}{'merge'}}) {
			foreach my $HFpos2 (sort {$a <=> $b} keys %{${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}}) {
				if (exists ${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom} and exists ${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}{$HFpos2} and ${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}{$HFpos2}=~/^\d+$/) {
					$HFincluded_region{$HFchrom}{$HFpos2}=${$HFreadsum}{$HFreadid}{'merge'}{$HFchrom}{$HFpos2};
				}
			}
		}
		if (0) {
			print "${HFsubinfo}Test: ";
			print "ReadID: $HFreadid\n";
			print "${HFsubinfo}Test: Merged Range: \n";
			print Dumper \%{${$HFreadsum}{$HFreadid}{'merge'}};
			print "${HFsubinfo}Test: included Range: \n";
			print Dumper \%HFincluded_region;
			print "${HFsubinfo}Test: readparts: \n";
			print Dumper \%{$HFreadparts{$HFreadid}};
		}
		my $HFnumberincludedregion=0;
		map {$HFnumberincludedregion+=scalar(keys %{$HFincluded_region{$_}})} keys %HFincluded_region;
		if (scalar($HFnumberincludedregion)==0) {
			${$HFreadsum}{$HFreadid}{'excluded'}++;
			delete $HFreadparts{$HFreadid};
			delete $HFreadpart2{$HFreadid} if (exists $HFreadpart2{$HFreadid});
#			print "${HFsubinfo}Test: Excluded ", $HFreadid, "\n";### For test ###
			next;
		}
		my $HFpartstart=0;
		foreach my $HFchrom (sort keys %{$HFreadparts{$HFreadid}}) {
			HFBLOCK1: {foreach my $HFpos (sort {$a <=> $b} keys %{$HFreadparts{$HFreadid}{$HFchrom}}) {
				foreach my $HFpos3 (keys %{$HFincluded_region{$HFchrom}}) {
					if ($HFpos>=$HFpos3 and $HFpos<=$HFincluded_region{$HFchrom}{$HFpos3}) {
						$HFpartstart=$HFpos3;
					}
				}
#				print "${HFsubinfo}Test: HFpartstart: ", $HFpartstart."\n"; ###For test ###
				if ($HFpartstart==0) {
					$HFerrorcode=1;
#					print "${HFsubinfo}Test: Errorcode $HFerrorcode at $HFchrom:$HFpos:$HFreadid\n";### For test ###		
					next HFBLOCK1;
				}
				my %HFalleleatpos=();
				foreach my $HFstrand (keys %{$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}}) {
					my @HFallelearr=keys %{$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}};
					my $HFallele=shift @HFallelearr;
#					print "${HFsubinfo}Test: chr:pos:allele:readid $HFchrom:$HFpos:$HFallele:$HFreadid\n"; ### For test ###
					if (exists $HFalleleatpos{$HFallele}) {
						if (ord($HFreadparts{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}{$HFallele})>ord($HFalleleatpos{$HFallele})) {
							$HFalleleatpos{$HFallele}=$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}{$HFallele};
						}
					}
					else {
						$HFalleleatpos{$HFallele}=$HFreadparts{$HFreadid}{$HFchrom}{$HFpos}{$HFstrand}{$HFallele};
					}
				}
#				print Dumper \%HFalleleatpos; ### For test###
				my @HFalleleatposit=keys %HFalleleatpos;
				if (scalar(@HFalleleatposit) ==1) {
					my $HFuniqallele=shift @HFalleleatposit;
#					print "${HFsubinfo}Test: unique allele: $HFuniqallele => $HFalleleatpos{$HFuniqallele}\n"; ### For test ###
					if (exists $HFcorrectedsites{$HFchrom} and $HFcorrectedsites{$HFchrom}{$HFpos} and exists $HFcorrectedsites{$HFchrom}{$HFpos}{$HFuniqallele} and $HFcorrectedsites{$HFchrom}{$HFpos}{$HFuniqallele}>0){
						$HFreadpart2{$HFreadid}{$HFchrom}{$HFpartstart}{$HFpos}{$HFuniqallele}=$HFalleleatpos{$HFuniqallele};
						$HFnumparts++;
					}
				}
			}}###HFBLOCK1
		}
		if  ($HFnumparts==0) {
			if ($HFerrorcode==1) {
#				print "${HFsubinfo}Test: Excluded ", $HFreadid, "\n";### For test ###
				${$HFreadsum}{$HFreadid}{'excluded'}++;
			}
			elsif ($HFerrorcode==0) {
#				print "${HFsubinfo}Test: Shared ", $HFreadid, "\n"; ### For test ###
				${$HFreadsum}{$HFreadid}{'shared'}++;
			}
			delete $HFreadparts{$HFreadid};
			delete $HFreadpart2{$HFreadid} if (exists $HFreadpart2{$HFreadid});
		}
		elsif ($HFnumparts>0) {
#			print "${HFsubinfo}Test: read have parts: $HFnumparts\n"; ### For test ###
#			print Dumer \%{$HFreadpart2{$HFreadid}}; ### For test ###
			my @HFparts=();
			my @HFpartpos=();
			foreach my $HFchrom (sort keys %{$HFreadpart2{$HFreadid}}) {
				foreach my $HFstart (sort {$a <=> $b} keys %{$HFreadpart2{$HFreadid}{$HFchrom}}) {
					my $HFprintline=$HFchrom;
					my $HFprintlinepos=$HFchrom;
					foreach my $HFpos (sort {$a <=> $b} keys %{$HFreadpart2{$HFreadid}{$HFchrom}{$HFstart}}) {
#						print "Pos: $HFpos\n";### For test ###
						if (exists $HFsitereferencer{'pos2ref'} and exists $HFsitereferencer{'pos2ref'}{$HFchrom} and exists $HFsitereferencer{'pos2ref'}{$HFchrom}{$HFpos}) {
							my $HFthissite=' ';
							my $HFthissitepos=' ';
							$HFthissite.=$HFsitereferencer{'pos2ref'}{$HFchrom}{$HFpos};
							$HFthissitepos.=$HFpos;
							my @HFallele=keys %{$HFreadpart2{$HFreadid}{$HFchrom}{$HFstart}{$HFpos}};
							my $HFthisallele=shift @HFallele;
							my $HFthisqual=$HFreadpart2{$HFreadid}{$HFchrom}{$HFstart}{$HFpos}{$HFthisallele};
							$HFthissite.=' '.$HFthisallele.' '.$HFthisqual;
							$HFthissitepos.=' '.$HFthisallele.' '.$HFthisqual;
							if ($HFthissite=~/^\s+\d+\s+\d+\s+\S+$/) {
								$HFprintline.=$HFthissite;
								$HFprintlinepos.=$HFthissitepos;
							}
							else {
								print STDERR "${HFsubinfo}Error: unknown fragment Read:Chr:Pos $HFchrom:$HFpos:$HFthissite\n";
								return $VcfKit_failure;
							}
	#						print "Printline: $HFprintline\n"; ### For test ###
						}
#						else {
#							print STDERR "${HFsubinfo}Error: no existing Read:Chr:Pos $HFreadid:$HFchrom:$HFpos in \$HFsitereferencer{'pos2ref'}\n";
#							return $VcfKit_failure;
#						}
					}
#					print "${HFsubinfo}Test: parts: $HFprintline\n";### For test ###
					if ($HFprintline=~/^$HFchrom\s+\d+\s+\d+\s+\S+/) {
						push (@HFparts, $HFprintline);
						push (@HFpartpos, $HFprintlinepos);
					}
					else {
						print STDERR "${HFsubinfo}Error: unknown ReadID:parts $HFreadid:$HFprintline\n";
					}
				}
			}
#			print "${HFsubinfo}Test: final parts(".scalar(@HFparts)."): @HFparts\n"; ### For test ###
			if (scalar(@HFparts)>0) {
				print FRAGMENTS scalar(@HFparts)."\t".$HFreadid."\t".join("\t",@HFparts)."\n";
				print FRAGMENTPOS scalar(@HFparts)."\t".$HFreadid."\t".join("\t",@HFpartpos)."\n";
#				print "${HFsubinfo}Test: ".scalar(@HFparts)."\t".$HFreadid."\t".join("\t",@HFparts)."\n";### For test###
			}
			else {
#				print "${HFsubinfo}Test: unknown reads: ", $HFreadid, "\n"; ### For test ###
				${$HFreadsum}{$HFreadid}{'excluded'}++;
				delete $HFreadparts{$HFreadid};
				delete $HFreadpart2{$HFreadid} if (exists $HFreadpart2{$HFreadid});
			}
		}
	}

	close FRAGMENTS;
	close FRAGMENTPOS;
	
	if (0) {### Test %HFreadpart2
		print "${HFsubinfo}Test: test \%HFreadpart2\n";
		print Dumper \%HFreadpart2;
		print "\n";
	}
	print "${HFsubinfo}Test: Total Allelic reads: ".scalar(keys %HFreadpart2)."\n";
	
	return ($VcfKit_success, \%HFreadpart2, \%HFcorrectedsites, \%HFsitereferencer);
}



### GroupFragments
### %GroupFragments=();
### Global:
### Dependency:
### Note:
sub GroupFragments {
	my ($GFfragments)=@_;
	
	my $GFsubinfo="SUB(GroupFragments)";
	my %GFcorrmatrix=();
	my %GFallelehash=();
	my %GFblock=();
	my $GFblocknum=0;
	
	#Format: %{$GFfragments}=($readid => chrom => $partstart => $pos => $allele => $qual);
	
	if (scalar(keys %{$GFfragments})==0) {
		return 1;
	}

	foreach my $GFreadid (keys %{$GFfragments}) {
		my @GFreadparts=();
		foreach my $GFchrom (keys %{${$GFfragments}{$GFreadid}}) {
			foreach my $GFblock (keys %{${$GFfragments}{$GFreadid}{$GFchrom}}) {
				foreach my $GFpos (keys %{${$GFfragments}{$GFreadid}{$GFchrom}{$GFblock}}) {
					my @GFalleles=keys %{${$GFfragments}{$GFreadid}{$GFchrom}{$GFblock}{$GFpos}};
					my $GFthisallele=shift @GFalleles;
					my $GFthisqual=${$GFfragments}{$GFreadid}{$GFchrom}{$GFpos}{$GFthisallele};
					$GFallelehash{$GFchrom}{$GFpos}{$GFthisallele}++;
#					print $GFsubinfo, "Test: $GFchrom-$GFpos-$GFthisallele\n";
					push (@GFreadparts, "$GFchrom-$GFpos-$GFthisallele");
				}
			}
		}
		foreach my $GFx (@GFreadparts) {
			foreach my $GFy (@GFreadparts) {
				$GFcorrmatrix{$GFx}{$GFy}++ unless ($GFx eq $GFy);
			}
		}
	}
	if (0) {
		print Dumper \%GFallelehash;
		print Dumper \%GFcorrmatrix;
	}
	my @GFchroms=sort keys %GFallelehash;
	foreach my $GFchrom1 (@GFchroms) {
		my @GFposits=sort {$a <=> $b} keys %{$GFallelehash{$GFchrom1}};
		foreach my $GFpos1 (@GFposits) {
			my @GFcodes=sort {$a <=> $b} keys %{$GFallelehash{$GFchrom1}{$GFpos1}};
			foreach my $GFincode1 (@GFcodes) {
				my @GFarr3=();
				print "$GFchrom1\t$GFpos1\t$GFincode1\t"; ### For test ###
				foreach my $GFchrom2 (@GFchroms) {
					foreach my $GFpos2 (sort {$a <=>$b} keys %{$GFallelehash{$GFchrom2}}) {
						foreach my $GFincode2 (sort {$a <=> $b} keys %{$GFallelehash{$GFchrom2}{$GFpos2}}) {
							if (exists $GFcorrmatrix{"$GFchrom1-$GFpos1-$GFincode1"} and exists $GFcorrmatrix{"$GFchrom1-$GFpos1-$GFincode1"}{"$GFchrom2-$GFpos2-$GFincode2"}) {
								push (@GFarr3, $GFcorrmatrix{"$GFchrom1-$GFpos1-$GFincode1"}{"$GFchrom2-$GFpos2-$GFincode2"});
							}
							else {
								push (@GFarr3, '.');
							}
						}
					}
				}
				print join ("\t", @GFarr3), "\n"; ### For test ###
			}
		}
	}
	return $VcfKit_success;
}


### Prepare VCF for hapcompass
### HapcompassVcf ($vcfin, $vcfout, \%filein_geno, \%polumorphic sites);
### Global: $VcfKit_failure; $VcfKit_success;
### Dependency: 
### Note:
### Return: 
sub HapcompassVcf {
	my ($HVvcfin, $HVvcfout, $HVfillingeno, $HVsites, $HVvcfreferencer)=@_;

	local *HVVCFIN; local *HVVCFOUT;
	my $HVsubinfo='SUB(VcfKit::HapcompassVcf)';
	my $HVnumline=0;
	my %HVcheck_duplicate=();
###Format: $fillingeno{chr}->{pos}='0/1/0'; ('A/B/D')
##Format: %HVref2vcf=($chr => referencer => vcfrecord)
	my %HVref2vcf=();

	unless (defined $HVvcfin and -s $HVvcfin) {
		print STDERR "${HVsubinfo}Error: invalid VCF input\n";
		return $VcfKit_failure;
	}
	unless (defined $HVvcfout and $HVvcfout=~/^\S+$/) {
		print STDERR "${HVsubinfo}Error: invalid VCF input\n";
		return $VcfKit_failure;
	}
	unlink $HVvcfout if (-e $HVvcfout);
	close HVVCFIN if (defined fileno(HVVCFIN));
	if ($HVvcfin=~/\.vcf.gz$/i) {
		unless (open (HVVCFIN, "zcat $HVvcfin | ")) {
			print STDERR "${HVsubinfo}Error: open VCF.gz input: $HVvcfin\n";
			return $VcfKit_failure;
		}
	}
	elsif ($HVvcfin=~/\.vcf$/i) {
		unless (open (HVVCFIN, " < $HVvcfin")) {
			print STDERR "${HVsubinfo}Error: open VCF.gz input: $HVvcfin\n";
			return $VcfKit_failure;
		}
	}
	else {
		print STDERR "${HVsubinfo}Error: unknown VCF format: $HVvcfin\n";
		return $VcfKit_failure;
	}
	close HVVCFOUT if (defined fileno(HVVCFOUT));
	unless (open(HVVCFOUT, " > $HVvcfout")) {
		print STDERR "${HVsubinfo}Error: can not write VCF output: $HVvcfout\n";
		return $VcfKit_failure;
	}
	print HVVCFOUT '##fileformat=VCFv4.1'."\n";
	print HVVCFOUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAABBDD\n";
	while (my $HVline=<HVVCFIN>) {
		$HVnumline++;
		chomp $HVline;
		if ($HVline=~/^#/) {
			if ($HVline=~/^#CHROM\s+POS\s+/) {
				print HVVCFOUT $HVline."\n";
			}
			next;
		}
		my @HVarr=split(/\t/, $HVline);
		if (scalar(@HVarr)<10) {
			print STDERR "${HVsubinfo}Error: ignored VCF line <10 col at line ($HVnumline) of VCF $HVvcfin\n";
			next;
		}
		unless (exists ${$HVsites}{$HVarr[0]} and exists ${$HVsites}{$HVarr[0]}{$HVarr[1]}) {
#			print STDERR "${HVsubinfo}Test: ignored sites at line ($HVnumline) of VCF $HVvcfin\n";### For test ###
			next;
		}
		$HVarr[2]=$HVarr[0].'_'.$HVarr[1];###Rename ID
		if (exists $HVcheck_duplicate{$HVarr[2]}) {
			print STDERR "${HVsubinfo}Error: ignored sites at line ($HVnumline) of VCF $HVvcfin\n";
			return $VcfKit_failure;
		}
		else {
			$HVcheck_duplicate{$HVarr[2]}++;
		}
		$HVarr[6]='PASS';
		$HVarr[7]='.';
		my @HVarr2=split(/:/, $HVarr[9]);
		if (exists ${$HVfillingeno}{$HVarr[0]} and exists ${$HVfillingeno}{$HVarr[0]}{$HVarr[1]} and ${$HVfillingeno}{$HVarr[0]}{$HVarr[1]}=~/\//) {
			$HVarr2[0]=${$HVfillingeno}{$HVarr[0]}{$HVarr[1]};
		}
		else {
			print STDERR "${HVsubinfo}Error: unknown genotypes at line ($HVnumline) of VCF $HVvcfin\n";
			return $VcfKit_failure;
		}
		$HVarr[9]=join(':', @HVarr2);
		if (exists ${$HVvcfreferencer}{'pos2ref'} and exists ${$HVvcfreferencer}{'pos2ref'}{$HVarr[0]} and exists ${$HVvcfreferencer}{'pos2ref'}{$HVarr[0]}{$HVarr[1]} and ${$HVvcfreferencer}{'pos2ref'}{$HVarr[0]}{$HVarr[1]}=~/^\d+$/) {
			my $HVthisreferencer=${$HVvcfreferencer}{'pos2ref'}{$HVarr[0]}{$HVarr[1]};
			if (defined $HVthisreferencer and $HVthisreferencer=~/^\d+$/) {
				$HVref2vcf{$HVarr[0]}{$HVthisreferencer}=join("\t", @HVarr);
			}
			else {
				print "${HVsubinfo}Error: not number referencer for Chr:Pos $HVarr[0]:$HVarr[1]\n";
				return $VcfKit_failure;
			}
		}
		else {
			print "${HVsubinfo}Error: not existed referencer for Chr:Pos $HVarr[0]:$HVarr[1]\n";
			return $VcfKit_failure;
		}
	}
	close HVVCFIN;
	
	foreach my $HVchrom (sort keys %HVref2vcf) {
		my $HVnumrec=0;
		foreach my $HVthisreferencer (sort {$a <=> $b} keys %{$HVref2vcf{$HVchrom}}) {
#			print "${HVsubinfo}Test: $HVchrom:$HVthisreferencer => $HVref2vcf{$HVchrom}{$HVthisreferencer}\n";### For test ###
			if ($HVnumrec eq $HVthisreferencer) {
				print HVVCFOUT $HVref2vcf{$HVchrom}{$HVthisreferencer}, "\n"; ### For test ###
			}
			else {
				print STDERR "${HVsubinfo}Error: missing the $HVnumrec-part\n";
				return $VcfKit_failure;
			}
			$HVnumrec++;
		}
	}
	close HVVCFOUT;
	return $VcfKit_success;
}



### Assign Hapcompass phased alleles to subgenome
### ReadHcOut()
### Global: $VcfKit_failure; $VcfKit_success;
### Dependency: 
### Note: 0=failure; 1=success
sub ReadHcOut {
	my ($RHOphasefile, $RHOgeno_ploidy, $RHOfixalleles, $RHOexcludedchrom)=@_;
	
	local *HAPCOMPASSOUT;
	my $RHOsubinfo='SUB(VcfKit::ReadHcOut)';
##Format: %{$RHOexcludedchrom}=(chrom => ++);

#	print $RHOsubinfo, "Test: \$RHOfixalleles\n"; print Dumper $RHOfixalleles; print "\n"; ### For test ###
	unless (defined $RHOphasefile and -s $RHOphasefile) {
		print STDERR "${RHOsubinfo}Error: invalid HapCompass output\n";
		return $VcfKit_failure;
	}
	unless (defined $RHOgeno_ploidy and $RHOgeno_ploidy=~/^[ABD]{1,3}$/) {
		print STDERR "${RHOsubinfo}Error: invalid ploidy\n";
		return $VcfKit_failure;
	}
### Initialising ploidy
	my $RHOnum_ploidy=length($RHOgeno_ploidy);
	my ($RHOaa_expressed, $RHObb_expressed, $RHOdd_expressed)=(0, 0, 0);
	my @RHOgenomes=();
	if ($RHOgeno_ploidy=~/A/) {
		$RHOaa_expressed=1;
		push (@RHOgenomes, 'A');
	}
	if ($RHOgeno_ploidy=~/B/) {
		$RHObb_expressed=1;
		push (@RHOgenomes, 'B');
	}
	if ($RHOgeno_ploidy=~/D/) {
		$RHOdd_expressed=1;
		push (@RHOgenomes, 'D');
	}

##FORMAT: $RHOfixalleles=(chr => (pos => (A => allele/?, B=> allele/?, D => allele/?)))
	close HAPCOMPASSOUT if (defined fileno(HAPCOMPASSOUT));
	unless (open HAPCOMPASSOUT, "<", $RHOphasefile) {
		print STDERR "${RHOsubinfo}Error: can not open $RHOphasefile\n";
		return $VcfKit_failure;
	}
##FORMAT: file format
#BLOCK   1       62      1       62      612.0001220703125       chr1
#rs5766699       1       1       0       1
#rs4823599       2       2       0       1
#rs7292231       3       3       0       1
#rs5767537       4       4       0       1
#rs738667        5       5       0       1
#rs738668        6       6       0       1
#rs738669        7       7       0       1
#rs5767541       8       8       0       1
#rs5767542       9       9       0       1
#rs5767543       10      10      0       1
#rs738670        11      11      0       1
#rs738671        12      12      0       1
#rs4823601       13      13      0       1
#rs5767545       14      14      0       1
#rs5767546       15      15      0       1
#rs5767547       16      16      0       1
	my $RHOblocknum=0;
##Format: %RHOtemp_block=($blocknum => $chrom => $pos => (1=> 0/1, 2 => 0/1, 3 => 0/1))
	my %RHOtemp_block=();

	my $RHOnumline=0;
	while (my $RHOline=<HAPCOMPASSOUT>) {
		$RHOnumline++;
		if ($RHOline=~/^BLOCK/) {
			$RHOblocknum++;
			next;
		}
		chomp $RHOline;
		next unless ($RHOline=~/\S+/);
		my @RHOarr1=split(/\t/, $RHOline);
		unless (scalar(@RHOarr1) == ($RHOnum_ploidy+3)) {
			close HAPCOMPASSOUT;
			print STDERR "${RHOsubinfo}Error: invalid phaseline\n";
			next;
		}
		my ($RHOchrom, $RHOposit);
		if ($RHOarr1[0]=~/^([\S+]+)_(\d+)$/) {
			$RHOchrom=$1;
			$RHOposit=$2;
		}
		else {
			print STDERR "${RHOsubinfo}Error: can not read chrom and pos from SNP ID: $RHOarr1[0]\n";
			return $VcfKit_failure;
		}
		for (my $RHOi=3; $RHOi<scalar(@RHOarr1); $RHOi++) {
			$RHOtemp_block{$RHOblocknum}{$RHOchrom}{$RHOposit}{$RHOi-2}=$RHOarr1[$RHOi];
		}
	}
	close HAPCOMPASSOUT;
#	print $RHOsubinfo, "Test: \%RHOtemp_block\n"; print Dumper \%RHOtemp_block; print "\n"; ### For test ###
#	print $RHOsubinfo, "Test: \$RHOfixalleles\n"; print Dumper $RHOfixalleles; print "\n"; ### For test ###

### count each group to each subgenome in each block
##FORMAT: %RHOtemp_assign=( block => (colnum => (A/B/D => total+count)));
	my %RHOtemp_assign=();
	my %RHOcorrelation=();
	foreach my $RHOblock (sort {$a <=> $b} (keys %RHOtemp_block)) {
		my %RHOfixarr=();
		my %RHOphase=();
		my $RHOarrindex=0;
		my %RHOnum2pos=();
		foreach my $RHOchrom2 (sort keys %{$RHOtemp_block{$RHOblock}}) {
#			print $RHOsubinfo, "Test: Counting: Chrom: $RHOchrom2\n"; ### for test ###
			foreach my $RHOposit2 (sort {$a <=> $b} keys %{$RHOtemp_block{$RHOblock}{$RHOchrom2}}) {
				my $RHOcount=0;
				$RHOnum2pos{$RHOarrindex}{$RHOchrom2}{$RHOposit2}++;
				for (my $RHOi2=1; $RHOi2<=$RHOnum_ploidy; $RHOi2++) {
					push (@{$RHOphase{$RHOi2}}, $RHOtemp_block{$RHOblock}{$RHOchrom2}{$RHOposit2}{$RHOi2});
					if (exists ${$RHOfixalleles}{$RHOchrom2} and exists ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}) {
						if ($RHOaa_expressed) {
							$RHOtemp_assign{$RHOblock}{$RHOi2}{'A'}=0 unless (exists $RHOtemp_assign{$RHOblock} and exists $RHOtemp_assign{$RHOblock}{$RHOi2} and $RHOtemp_assign{$RHOblock}{$RHOi2}{'A'});
							if (exists ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'A'} and ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'A'}=~/^[?0-9]{1}$/) {
								$RHOtemp_assign{$RHOblock}{$RHOi2}{'A'}++ if (${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'A'} eq $RHOtemp_block{$RHOblock}{$RHOchrom2}{$RHOposit2}{$RHOi2});
								push (@{$RHOfixarr{'A'}}, ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'A'}) if ($RHOcount==0);
							}
							else {
								print STDERR "${RHOsubinfo}Error: Final AA assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
								return $VcfKit_failure;
							}
						}
						if ($RHObb_expressed) {
							$RHOtemp_assign{$RHOblock}{$RHOi2}{'B'}=0 unless (exists $RHOtemp_assign{$RHOblock} and exists $RHOtemp_assign{$RHOblock}{$RHOi2} and $RHOtemp_assign{$RHOblock}{$RHOi2}{'B'});
							if (exists ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'B'}) {
								$RHOtemp_assign{$RHOblock}{$RHOi2}{'B'}++ if (${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'B'} eq $RHOtemp_block{$RHOblock}{$RHOchrom2}{$RHOposit2}{$RHOi2});
								push (@{$RHOfixarr{'B'}}, ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'B'}) if ($RHOcount==0);;
							}
							else {
								print STDERR "${RHOsubinfo}Error: Final BB assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
								return $VcfKit_failure;
							}
						}
						if ($RHOdd_expressed) {
							$RHOtemp_assign{$RHOblock}{$RHOi2}{'D'}=0 unless (exists $RHOtemp_assign{$RHOblock} and exists $RHOtemp_assign{$RHOblock}{$RHOi2} and $RHOtemp_assign{$RHOblock}{$RHOi2}{'D'});
							if (exists ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'D'}) {
								$RHOtemp_assign{$RHOblock}{$RHOi2}{'D'}++ if (${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'D'} eq $RHOtemp_block{$RHOblock}{$RHOchrom2}{$RHOposit2}{$RHOi2});
								push (@{$RHOfixarr{'D'}}, ${$RHOfixalleles}{$RHOchrom2}{$RHOposit2}{'D'}) if ($RHOcount==0);;
							}
							else {
								print STDERR $RHOsubinfo, "Error: Final DD assign error at Blk:Chr:Pos:Col: $RHOblock:$RHOchrom2:$RHOposit2:$RHOi2\n";
								return $VcfKit_failure;
							}
						}
					}
					else {
						print STDERR $RHOsubinfo, "Error: no fixed allele at chr:pos $RHOchrom2:$RHOposit2\n";
						return $VcfKit_failure;
					}
					$RHOcount++;
#					print $RHOsubinfo, "Test: phased $RHOchrom2:$RHOposit2 @{$RHOphase{$RHOi2}}\n";### For test ###
#					print $RHOsubinfo, "Test: A $RHOchrom2:$RHOposit2 @{$RHOfixarr{'A'}}\n";### For test ###
#					print $RHOsubinfo, "Test: B $RHOchrom2:$RHOposit2 @{$RHOfixarr{'B'}}\n";### For test ###
#					print $RHOsubinfo, "Test: D $RHOchrom2:$RHOposit2 @{$RHOfixarr{'D'}}\n";### For test ###
				}
			}
		}
#		print $RHOsubinfo, "Test: \%RHOfixarr\n"; print Dumper \%RHOfixarr; print "\n"; ### For test ###
#		print $RHOsubinfo, "Test: \%RHOphase\n"; print Dumper \%RHOphase; print "\n"; ### For test ###
#		print $RHOsubinfo, "Test: \%RHOnum2pos\n"; print Dumper \%RHOnum2pos; print "\n"; ### For test ###
#		print $RHOsubinfo, "Test: \%RHOtemp_assign\n"; print Dumper \%RHOtemp_assign; print "\n"; ### For test ###
### Correlation
		foreach my $RHOind_fixallele (sort keys %RHOfixarr) {
#			print $RHOsubinfo, "Test: fixedarr1(BLK$RHOblock Geno$RHOind_fixallele): @{$RHOfixarr{$RHOind_fixallele}}\n";
			foreach my $RHOind_phasing (sort {$a <=> $b} keys %RHOphase) {
#				print ${RHOsubinfo}, "Test: phasedarr2(BLK$RHOblock COL$RHOind_phasing): @{$RHOphase{$RHOind_phasing}}\n";
				my ($RHOcorr_retcode, $RHOcor, $RHOcov, $RHOnum)=CallCorrelation(\@{$RHOfixarr{$RHOind_fixallele}}, \@{$RHOphase{$RHOind_phasing}});
				unless ($RHOcorr_retcode) {
#					print STDERR $RHOsubinfo, "Warnings: CallCorrelation failed\n";
				}
				else {
					$RHOcorrelation{$RHOblock}{$RHOind_phasing}{$RHOind_fixallele}=$RHOcov;
#					print $RHOsubinfo, "Test: \nArr1: @{$RHOfixarr{$RHOind_fixallele}}\nArr2: @{$RHOphase{$RHOind_phasing}}\nCov: $RHOcov Cor: $RHOcor\n"; 
				}
			}
		}
	}
#	print $RHOsubinfo, "Test: \%RHOcorrelation\n"; print Dumper \%RHOcorrelation; print "\n"; ### For test ###
#	print $RHOsubinfo, "Test: \%RHOtemp_assign\n"; print Dumper \%RHOtemp_assign; print "\n"; ### For test ###
	
### assign each group to each sungenome in each block
	my $RHOtestcmd_count=0;
	RHOBLOCK1: {foreach my $RHOind_block2 (sort {$a <=> $b} (keys %RHOtemp_assign)) {
		my $test_corelation=0;
		my %RHOgroup2genome=();
		my %RHOgenome2group=();
		RHOBLOCK2: {foreach my $RHOind_group (sort {$a <=> $b} (keys %{$RHOtemp_assign{$RHOind_block2}})) {
			my $RHObest_group='unknown';
			my $RHOmax_count=0;
			my $RHOnum_max=0;
			foreach my $RHOgenotype (@RHOgenomes) {
				if (exists $RHOtemp_assign{$RHOind_block2}{$RHOind_group}{$RHOgenotype}){
					if ($RHOtemp_assign{$RHOind_block2}{$RHOind_group}{$RHOgenotype} > $RHOmax_count) {
						$RHOnum_max=1;
						$RHOmax_count=$RHOtemp_assign{$RHOind_block2}{$RHOind_group}{$RHOgenotype};
						$RHObest_group=$RHOgenotype;
					}
					elsif ($RHOtemp_assign{$RHOind_block2}{$RHOind_group}{$RHOgenotype} == $RHOmax_count) {
						$RHOnum_max++;
						$RHObest_group.="$RHOgenotype";
					}
				}
			}
			if ($RHOnum_max !=1 or $RHObest_group !~ m/^[ABD]{1}$/) {
				$test_corelation=1;
				next;
			}
			$RHOgroup2genome{$RHOind_group}=$RHObest_group;
			if (exists $RHOgenome2group{$RHObest_group}) {
				my $RHOlastgroup=$RHOgenome2group{$RHObest_group};
				if ($RHOtemp_assign{$RHOind_block2}{$RHOind_group}{$RHObest_group} > $RHOtemp_assign{$RHOind_block2}{$RHOlastgroup}{$RHObest_group}) {
					$RHOgenome2group{$RHObest_group}=$RHOind_group;
					delete $RHOgroup2genome{$RHOlastgroup};
					$RHOgroup2genome{$RHOind_group}=$RHObest_group;
				}
				elsif ($RHOtemp_assign{$RHOind_block2}{$RHOind_group}{$RHObest_group} == $RHOtemp_assign{$RHOind_block2}{$RHOlastgroup}{$RHObest_group}) {
					if ($RHOcorrelation{$RHOind_block2}{$RHOind_group}{$RHObest_group} < $RHOcorrelation{$RHOind_block2}{$RHOlastgroup}{$RHObest_group}) {
						$RHOgenome2group{$RHObest_group}=$RHOind_group;
						delete $RHOgroup2genome{$RHOlastgroup};
						$RHOgroup2genome{$RHOind_group}=$RHObest_group;
					}
					elsif ($RHOcorrelation{$RHOind_block2}{$RHOind_group}{$RHObest_group} == $RHOcorrelation{$RHOind_block2}{$RHOlastgroup}{$RHObest_group}) {
						if (exists $RHOtemp_block{$RHOind_block2}) {
							my $RHOtest_same=1;
							foreach my $RHOchrom (keys %{$RHOtemp_block{$RHOind_block2}}) {
								foreach my $RHOpos (keys %{$RHOtemp_block{$RHOind_block2}{$RHOchrom}}) {
									unless (exists $RHOtemp_block{$RHOind_block2}{$RHOchrom}{$RHOpos}{$RHOind_group} and exists $RHOtemp_block{$RHOind_block2}{$RHOchrom}{$RHOpos}{$RHOind_group} and $RHOtemp_block{$RHOind_block2}{$RHOchrom}{$RHOpos}{$RHOind_group} eq $RHOtemp_block{$RHOind_block2}{$RHOchrom}{$RHOpos}{$RHOind_group}) {
										$RHOtest_same=0;
									}
								}
							}
							if ($RHOtest_same==1) {
								$RHOgenome2group{$RHObest_group}=$RHOind_group;
							}
						}
					}
					
				}
				$test_corelation=1;
			}
			else {
				$RHOgenome2group{$RHObest_group}=$RHOind_group;
			}
		}}###RHOBLOCK2
#		print $RHOsubinfo, "Test: BLK($RHOind_block2): \%RHOgenome2group\n"; print Dumper \%RHOgenome2group; print "\n";
#		print $RHOsubinfo, "Test: BLK($RHOind_block2): \%RHOgroup2genome\n"; print Dumper \%RHOgroup2genome; print "\n";
		if (scalar(keys %RHOgenome2group)<($RHOnum_ploidy-1)) {
			foreach my $RHOchrom (keys %{$RHOtemp_block{$RHOind_block2}}) {
#				print $RHOsubinfo, "Test: BLK($RHOind_block2) Exclued(1): $RHOchrom\n"; print Dumper \%RHOgenome2group; print "\n";### For test ###
				${$RHOexcludedchrom}{$RHOchrom}++;
			}
		}
		elsif (scalar(keys %RHOgenome2group)==($RHOnum_ploidy-1)) {
			my %RHOhash1=();
			my %RHOhash2=();
			for (my $RHOi=1; $RHOi<=$RHOnum_ploidy; $RHOi++) {
				$RHOhash1{$RHOi}++;
			}
			foreach (@RHOgenomes) {
				$RHOhash2{$_}++;
			}
			foreach (keys %RHOgenome2group) {
				if (exists $RHOhash2{$_}) {
					delete $RHOhash2{$_};
					my $RHOtarget=$RHOgenome2group{$_};
					delete $RHOhash1{$RHOtarget} if (exists $RHOhash1{$RHOtarget});
				}
			}
			if (scalar(keys %RHOhash1)==1 and scalar(keys %RHOhash2)==1) {
				my @RHOarr1=keys %RHOhash1;
				my @RHOarr2=keys %RHOhash2;
				$RHOgenome2group{$RHOarr2[0]}=$RHOarr1[0];
			}
			else {
				foreach my $RHOchrom (keys %{$RHOtemp_block{$RHOind_block2}}) {
#					print $RHOsubinfo, "Test: BLK($RHOind_block2) Exclued(2): $RHOchrom\n"; print Dumper \%RHOgenome2group; print "\n";### For test ###
					${$RHOexcludedchrom}{$RHOchrom}++;
				}
			}
		}
		if (scalar(keys %RHOgenome2group)==$RHOnum_ploidy) {
			foreach my $RHOchrom (keys %{$RHOtemp_block{$RHOind_block2}}) {
				foreach my $RHOpos (keys %{$RHOtemp_block{$RHOind_block2}{$RHOchrom}}) {
					foreach my $RHOgeno (@RHOgenomes) {
						my $RHOnum=$RHOgenome2group{$RHOgeno};
						if (exists $RHOtemp_block{$RHOind_block2}{$RHOchrom}{$RHOpos}{$RHOnum}){
							${$RHOfixalleles}{$RHOchrom}{$RHOpos}{$RHOgeno}=$RHOtemp_block{$RHOind_block2}{$RHOchrom}{$RHOpos}{$RHOnum};
						}
					}
				}
			}
		}
#		print $RHOsubinfo, "Test: BLOCK $RHOind_block2 \%RHOgroup2genome\n"; print Dumper \%RHOgroup2genome; print "\n"; ### For test ###
#		print $RHOsubinfo, "Test: BLOCK $RHOind_block2 \%RHOgenome2group\n"; print Dumper \%RHOgenome2group; print "\n"; ### For test ###
	}}###RHOBLOCK1
#	print $RHOsubinfo, "Test: \$RHOfixalleles\n"; print Dumper $RHOfixalleles; print "\n"; ### For test ###
#	print $RHOsubinfo, "Test: \$RHOexcludedchrom\n"; print Dumper $RHOexcludedchrom; print "\n"; ### For test ###
	foreach my $RHOchrom (keys %{$RHOfixalleles}) {
		next if (exists ${$RHOexcludedchrom}{$RHOchrom});
		my $RHOerrorcode=0;
		RHOBLOCK3: {foreach my $RHOpos (keys %{${$RHOfixalleles}{$RHOchrom}}) {
			foreach my $RHOcode (@RHOgenomes) {
#				print $RHOsubinfo, "Test: Chr:Pos:Geno: $RHOchrom:$RHOpos:$RHOcode\n"; ### For test ###
				unless (exists ${$RHOfixalleles}{$RHOchrom}{$RHOpos}{$RHOcode} and ${$RHOfixalleles}{$RHOchrom}{$RHOpos}{$RHOcode}=~/^\d+$/) {
					$RHOerrorcode=1;
#					print STDERR $RHOsubinfo, "Test: Excluded $RHOchrom at Chr:Pos:Geno: $RHOchrom:$RHOpos:$RHOcode\n"; ### For test ###
					last RHOBLOCK3;
				}
			}
		}}###RHOBLOCK3;
		${$RHOexcludedchrom}{$RHOchrom}++ if ($RHOerrorcode==1);
	}
	foreach (keys %{$RHOexcludedchrom}) {
		delete ${$RHOfixalleles}{$_} if (exists ${$RHOfixalleles}{$_});
	}
	
#	print $RHOsubinfo, "Test: \$RHOfixalleles\n"; print Dumper $RHOfixalleles; print "\n"; ### For test ###
#	print $RHOsubinfo, "Test: \$RHOexcludedchrom\n"; print Dumper $RHOexcludedchrom; print "\n"; ### For test ###
	if (scalar(keys %{$RHOfixalleles})>0) {
		return $VcfKit_success;
	}
	else {
		print STDERR $RHOsubinfo, "Error: Fixed alleles ==0 \n";
		return $VcfKit_failure;
	}
}



#my $RFadditional_cmd=" --min-coverage $freebayes_min_coverage --min-alternate-count $freebayes_min_alternative_count --min-mapping-quality $min_mapq --genotype-qualities --pooled-discrete ";
### RunFreebayes, return merged vcf and raw AABBDD vcf
### &RunFreebayes($ref.fasta, \@seqids, $bam, $ploidy, $output.vcf.gz, $guided.vcf/0, $additions_cmd, [$run_parallel(1/0)], [Parallel_region file], [$threads], [$fastabin], [path_freebayes], [path_freebayesparallel], [path_vcfsort], [path_concat], [path_tabix], [path_bgzip], [path_samtools])
### 1=success; 0=failure
### Global: 
### Dependancy: $VcfKit_failure;$VcfKit_success;
### Note:
sub RunFreebayes {
	my ($RFfile_reference, $RFseqids_arrindex, $RFfile_bam, $RFploidy, $RFoutput, $RFguide_vcf, $RFadditional_cmd,   $RFrun_parallel, $RFregionfile, $RFthreads, $RFfastabin, $RFpath_freebayes, $RFpath_freebayesparallel, $RFpath_vcfsort, $RFpath_vcfconcat, $RFpath_tabix, $RFpath_bgzip, $RFpath_samtools)=@_;
	
	my $RFsubinfo='SUB(VcfKit::RunFreebayes)';
	my $RFcmd='';
	my $RFtest_createregion=0; ##0 =not creat region file; will set to 0 if defined $RFregionfile
	$RFadditional_cmd=' ' unless (defined $RFadditional_cmd);
	$RFrun_parallel=0 unless (defined $RFrun_parallel);
	$RFthreads=1 unless (defined $RFthreads);
	$RFfastabin=200 unless (defined $RFfastabin);
	$RFpath_freebayes='freebayes' unless (defined $RFpath_freebayes);
	$RFpath_freebayesparallel='freebayes-parallel' unless (defined $RFpath_freebayesparallel);
	$RFpath_vcfsort='vcf-sort' unless (defined $RFpath_vcfsort);
	$RFpath_bgzip='bgzip' unless (defined $RFpath_bgzip);
	$RFpath_samtools='samtools' unless (defined $RFpath_samtools);
	$RFpath_vcfconcat='vcf-concat' unless (defined $RFpath_vcfconcat);
	my $RFrunguided_freebayes=0;
	
	unless (defined $RFfile_reference and -s $RFfile_reference) {
		print STDERR "${RFsubinfo}Error: invalid freebayes reference\n";
		return $VcfKit_failure;
	}
	if (scalar(@{$RFseqids_arrindex})<1) {
		print "${RFsubinfo}Error: no seqs for freebayes";
		return 1;
	}
	unless (defined $RFfile_bam and -s $RFfile_bam) {
		print STDERR "${RFsubinfo}Error: invalid freebayes BAM\n";
		return $VcfKit_failure;
	}
	unless (defined $RFploidy and $RFploidy=~/^\d+$/ and $RFploidy>1) {
		print STDERR "${RFsubinfo}Error: invalid freebayes ploidy\n";
		return $VcfKit_failure;
	}
	unlink $RFoutput if (-e $RFoutput);
	unlink "$RFoutput.tbi" if (-e "$RFoutput.tbi");
	if (defined $RFregionfile and -s $RFregionfile) {
		$RFrun_parallel=1;
		$RFtest_createregion=0;
	}
	elsif ($RFrun_parallel) {
		$RFtest_createregion=1;
	}
	if (defined $RFguide_vcf and $RFguide_vcf ne '0' and -s $RFguide_vcf) {
		$RFrunguided_freebayes=1;
		if ($RFguide_vcf=~/\.vcf$/i) {
			if ( ! &BgzipVcf($RFguide_vcf, "$RFguide_vcf.gz", $RFpath_bgzip)) {
				print STDERR "${RFsubinfo}Error: bgzip VCF for freebayes failed: $RFguide_vcf\n";
				return $VcfKit_failure;
			}
			$RFguide_vcf.='.gz';
		}
		if ($RFguide_vcf=~/\.vcf.gz$/i) {
			if (! -s "$RFguide_vcf.tbi") {
				if (! &IndexVcf ($RFguide_vcf, $RFpath_tabix)) {
					print STDERR "${RFsubinfo}Error: index VCF for freebayes failed: $RFguide_vcf\n";
					return $VcfKit_failure;
				}
			}
		}
		else {
			print STDERR "${RFsubinfo}Error: unknown VCF format: $RFguide_vcf\n";
			return $VcfKit_failure;
		}
	}
	
	my $RFcmd_freebayes=$RFpath_freebayes;
	if ($RFtest_createregion) {
		$RFregionfile="$RFfile_reference.$RFfastabin.region";
		if (! &CreateFastaRegion("$RFfile_reference.fai", $RFfastabin, $RFregionfile, $RFpath_samtools)) {
			print STDERR "${RFsubinfo}Error: CreateFastaRegion for $RFfile_reference\n";
			return $VcfKit_failure;
		}
	}
	if ($RFrun_parallel) {
		$RFcmd_freebayes=$RFpath_freebayesparallel." $RFregionfile $RFthreads ";
	}
	$RFcmd_freebayes.=" --fasta-reference $RFfile_reference --ploidy $RFploidy $RFadditional_cmd ";
	 
	if (! $RFrunguided_freebayes) {
		$RFcmd="$RFcmd_freebayes $RFfile_bam | $RFpath_vcfsort | $RFpath_bgzip > $RFoutput";
		if (! &exec_cmd_return($RFcmd)) {
			print STDERR "${RFsubinfo}Error1:  freebayes running error\n";
			return $VcfKit_failure;
		}
		if (! -s $RFoutput) {
			print STDERR "${RFsubinfo}Error1: freebayes output error\n";
			return $VcfKit_failure;
		}
	}
	elsif ($RFrunguided_freebayes) {
		$RFcmd_freebayes.=" --variant-input $RFguide_vcf --only-use-input-alleles ";
		if ($RFrun_parallel) {
			$RFcmd="$RFcmd_freebayes $RFfile_bam | $RFpath_bgzip > $RFoutput";
			if (! &exec_cmd_return($RFcmd)) {
					print STDERR "${RFsubinfo}Error2:  freebayes running error\n";
					return $VcfKit_failure;
			}
			if (! -s $RFoutput) {
				print STDERR "${RFsubinfo}Error2: freebayes output error\n";
				return $VcfKit_failure;
			}
		}
		else {
			if (scalar(@{$RFseqids_arrindex}) ==1) {
				$RFcmd="$RFcmd_freebayes --region ${$RFseqids_arrindex}[0]  $RFfile_bam | $RFpath_vcfsort | $RFpath_bgzip > $RFoutput";
				if (! &exec_cmd_return($RFcmd)) {
					print STDERR "${RFsubinfo}Error3:  freebayes running error\n";
					return $VcfKit_failure;
				}
				if (! -s $RFoutput) {
					print STDERR "${RFsubinfo}Error3: freebayes output error\n";
					return $VcfKit_failure;
				}
			}
			else {#freebayes can not tolerate more than 1 seqs using --variant-input, so need to separate the seqids
				my @RFtemp_vcf=();
				foreach my $RFind_seq (@{$RFseqids_arrindex}) {
					$RFcmd="$RFcmd_freebayes --region $RFind_seq  $RFfile_bam | $RFpath_vcfsort | $RFpath_bgzip > $RFind_seq.vcf.gz";
					if (! &exec_cmd_return($RFcmd)) {
						print STDERR "${RFsubinfo}Error4:  freebayes running error\n";
						return $VcfKit_failure;
					}
					if (! -s "$RFind_seq.vcf.gz") {
						print STDERR "${RFsubinfo}Error4: freebayes output error\n";
						return $VcfKit_failure;
					}
					push (@RFtemp_vcf, "$RFind_seq.vcf.gz");
				}
				$RFcmd="$RFpath_vcfconcat ".join(' ', @RFtemp_vcf)." | $RFpath_bgzip > $RFoutput";
				if (! &exec_cmd_return($RFcmd)) {
					print STDERR "${RFsubinfo}Error5:  freebayes running error\n";
					return $VcfKit_failure;
				}
				if (! -s $RFoutput) {
					print STDERR "${RFsubinfo}Error5: freebayes output error\n";
					return $VcfKit_failure;
				}
				unlink (@RFtemp_vcf); ###delete temporary files
			}
		}
	}

	if (! &IndexVcf($RFoutput, $RFpath_tabix)) {
		print STDERR "${RFsubinfo}Error: freebayes output index error: $RFoutput\n";
		return $VcfKit_failure;
	}

	return $VcfKit_success;
}



###
sub CorrectAlleles {
	my ($CAgnploidy, $CAfragments, $CAreadsum, $CAalleleref)=@_;
	
	my $CAsubinfo="SUB(VcfKit::CorrectAlleles)";
###Format: %CAalleleassign=(chr => pos => ('A' =>0/1, 'B'=>0/1, 'D'=>0/1))
	my %CAalleleassign=();
##Format: %CAallelecount=(chr => pos => (0 => ('A' => ++, 'B' => ++, 'D' => ++)
##										1 => ...))
	my %CAallelecount=();
##Format: %CAallelecount=(chr => pos => (0 => ('A' => perc, 'B' => perc, 'D' => perc)
##										1 => ...))
	my %CAalleleperc=();
##Format: %genoploidy=('A' => ++, 'B' => ++, 'D' => ++)
	my %CAgenoploidy=();
##Format: %CAgenodepth=(chr => pos => allele => depth)
	my %CAgenodepth=();
	
	if ($CAgnploidy=~/A/) {
		$CAgenoploidy{'A'}++;
	}
	if ($CAgnploidy=~/B/) {
		$CAgenoploidy{'B'}++;
	}
	if ($CAgnploidy=~/D/) {
		$CAgenoploidy{'D'}++;
	}
	
#	$HFreadpart2{$HFreadid}{$HFchrom}{$HFpartstart}{$HFpos}{$HFuniqallele}=$HFalleleatpos{$HFuniqallele};
	foreach my $CAreadid (keys %{$CAfragments}) {
		foreach my $CAchrom (keys %{${$CAfragments}{$CAreadid}}) {
			foreach my $CAstart (keys %{${$CAfragments}{$CAreadid}{$CAchrom}}) {
				foreach my $CApos (keys %{${$CAfragments}{$CAreadid}{$CAchrom}{$CAstart}}) {
					foreach my $CAalelecode (keys %{${$CAfragments}{$CAreadid}{$CAchrom}{$CAstart}{$CApos}}) {
						$CAgenodepth{$CAchrom}{$CApos}{$CAalelecode}++;
						if (exists ${$CAreadsum}{$CAreadid} and exists ${$CAreadsum}{$CAreadid}{'gen'}) {
							if (exists ${$CAreadsum}{$CAreadid}{'gen'}{'A'} and ${$CAreadsum}{$CAreadid}{'gen'}{'A'}>0) {
								$CAallelecount{$CAchrom}{$CApos}{'A'}{$CAalelecode}++;
							}
							if (exists ${$CAreadsum}{$CAreadid}{'gen'}{'B'} and ${$CAreadsum}{$CAreadid}{'gen'}{'B'}>0) {
								$CAallelecount{$CAchrom}{$CApos}{'B'}{$CAalelecode}++;
							}
							if (exists ${$CAreadsum}{$CAreadid}{'gen'}{'D'} and ${$CAreadsum}{$CAreadid}{'gen'}{'D'}>0) {
								$CAallelecount{$CAchrom}{$CApos}{'D'}{$CAalelecode}++;
							}
						}
						else {
							$CAallelecount{$CAchrom}{$CApos}{'A'}{$CAalelecode}=0 unless (exists $CAallelecount{$CAchrom} and exists $CAallelecount{$CAchrom}{$CApos} and exists $CAallelecount{$CAchrom}{$CApos}{'A'} and exists $CAallelecount{$CAchrom}{$CApos}{'A'}{$CAalelecode});
							$CAallelecount{$CAchrom}{$CApos}{'B'}{$CAalelecode}=0 unless (exists $CAallelecount{$CAchrom} and exists $CAallelecount{$CAchrom}{$CApos} and exists $CAallelecount{$CAchrom}{$CApos}{'B'} and exists $CAallelecount{$CAchrom}{$CApos}{'B'}{$CAalelecode});
							$CAallelecount{$CAchrom}{$CApos}{'D'}{$CAalelecode}=0 unless (exists $CAallelecount{$CAchrom} and exists $CAallelecount{$CAchrom}{$CApos} and exists $CAallelecount{$CAchrom}{$CApos}{'D'} and exists $CAallelecount{$CAchrom}{$CApos}{'D'}{$CAalelecode});
						}
					}
				}
			}
		}
	}
#	print $CAsubinfo, "Test: \%CAallelecount\n"; print Dumper \%CAallelecount; print "\n"; ### For test ###
	foreach my $CAchrom (keys %CAallelecount) {
		foreach my $CApos (keys %{$CAallelecount{$CAchrom}}) {
			foreach my $CAgeno (sort keys %{$CAallelecount{$CAchrom}{$CApos}}) {
				 foreach my $CAallele (keys %{$CAallelecount{$CAchrom}{$CApos}{$CAgeno}}) {
#				 	print $CAsubinfo, "Test: chr:pos:geno:allele:count:depth:perc $CAchrom:$CApos:$CAgeno:$CAallele:".$CAallelecount{$CAchrom}{$CApos}{$CAgeno}{$CAallele}.':'.$CAgenodepth{$CAchrom}{$CApos}{$CAallele}.':'.($CAallelecount{$CAchrom}{$CApos}{$CAgeno}{$CAallele}/$CAgenodepth{$CAchrom}{$CApos}{$CAallele})."\n"; ### For test ###
					if (exists $CAgenodepth{$CAchrom} and exists $CAgenodepth{$CAchrom}{$CApos} and $CAgenodepth{$CAchrom}{$CApos}{$CAallele} and $CAgenodepth{$CAchrom}{$CApos}{$CAallele}>0) {
						$CAalleleperc{$CAchrom}{$CApos}{$CAgeno}{$CAallele}=$CAallelecount{$CAchrom}{$CApos}{$CAgeno}{$CAallele}/$CAgenodepth{$CAchrom}{$CApos}{$CAallele};
					}
					else {
						$CAalleleperc{$CAchrom}{$CApos}{$CAgeno}{$CAallele}=0;
					}
				 }
			}
		}
	}
#	print $CAsubinfo, "Test: \%CAalleleperc\n"; print Dumper \%CAalleleperc; print "\n"; ### For test ###
#	print $CAsubinfo, "Test: \%CAgenodepth\n"; print Dumper \%CAgenodepth; print "\n"; ### For test ###

	%CAallelecount=();
	
	foreach my $CAchrom (keys %CAalleleperc) {
		foreach my $CApos (keys %{$CAalleleperc{$CAchrom}}) {
			foreach my $CAgeno (keys %CAgenoploidy) {
				my $CAmaxcount=0;
				my $CAbestallele='?';
				my %CAallelehash=();
				my $CAnumbest=0;
				unless (exists $CAalleleperc{$CAchrom}{$CApos}{$CAgeno} and scalar(keys %{$CAalleleperc{$CAchrom}{$CApos}{$CAgeno}}>0)) {
					if (exists ${$CAalleleref}{$CAchrom} and exists ${$CAalleleref}{$CAchrom}{$CApos} and exists ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno} and ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno}=~/^\d+$/) {
						$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}=${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno};
					}
					else {
						$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}='?';
					}
				}
				else {
					foreach my $CAallelecode (keys %{$CAalleleperc{$CAchrom}{$CApos}{$CAgeno}}) {
						if ($CAalleleperc{$CAchrom}{$CApos}{$CAgeno}{$CAallelecode}>$CAmaxcount) {
							$CAmaxcount=$CAalleleperc{$CAchrom}{$CApos}{$CAgeno}{$CAallelecode};
							%CAallelehash=();
							$CAbestallele=$CAallelecode;
							$CAallelehash{$CAbestallele}++;
							$CAnumbest=1;
						}
						elsif ($CAalleleperc{$CAchrom}{$CApos}{$CAgeno}{$CAallelecode}==$CAmaxcount) {
							$CAnumbest++;
							$CAallelehash{$CAallelecode}++;
						}
					}
					if ($CAnumbest==0) {
						if (exists ${$CAalleleref}{$CAchrom} and exists ${$CAalleleref}{$CAchrom}{$CApos} and exists ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno} and ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno}=~/^\d+$/) {
							$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}=${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno};
						}
						else {
							$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}='?';
						}
					}
					elsif ($CAnumbest==1) {
						if (exists ${$CAalleleref}{$CAchrom} and exists ${$CAalleleref}{$CAchrom}{$CApos} and exists ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno} and ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno}=~/^\d+$/) {
							if ($CAbestallele eq ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno}) {
								$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}=$CAbestallele;
							}
							else {
								$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}='?';
							}
						}
						else {
							$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}=$CAbestallele;
						}
					}
					elsif ($CAnumbest>1) {
						if (exists ${$CAalleleref}{$CAchrom} and exists ${$CAalleleref}{$CAchrom}{$CApos} and exists ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno} and ${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno}=~/^\d+$/) {
							my $CAindref_allele=${$CAalleleref}{$CAchrom}{$CApos}{$CAgeno};
							if (exists $CAallelehash{$CAindref_allele}) {
								$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}=$CAindref_allele;
							}
							else {
								$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}='?';
							}
						}
						else {
							$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}='?';
						}
					}
				}
			}
		}
	}
## check not polymorphic sites
	foreach my $CAchrom (keys %CAalleleassign) {
		foreach my $CApos (keys %{$CAalleleassign{$CAchrom}}) {
			my %CAalleleverify=();
			my $CAtest_unknown=0;
			foreach my $CAgeno (keys %{$CAalleleassign{$CAchrom}{$CApos}}) {
				if ($CAalleleassign{$CAchrom}{$CApos}{$CAgeno} =~ /^\d+$/) {
					$CAalleleverify{$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}}++;
				}
				else {
					$CAtest_unknown=1;
					last;
				}
			}
			if ($CAtest_unknown==0) {
#				print $CAsubinfo, "Test: Chr:Pos $CAchrom:$CApos AlleleNum: ".scalar(keys %CAalleleverify)."\n"; ### For test ###
				unless (scalar(keys %CAalleleverify) == 2) {
					foreach my $CAgeno (keys %{$CAalleleassign{$CAchrom}{$CApos}}) {
						$CAalleleassign{$CAchrom}{$CApos}{$CAgeno}='?';
					}
				}
			}
		}
	}
	if (0) {
		print $CAsubinfo, "Test: \%CAalleleassign\n";
		print Dumper \%CAalleleassign;
		print "\n";
	}
	return ($VcfKit_success, \%CAalleleassign, \%CAgenodepth);
}

1;
