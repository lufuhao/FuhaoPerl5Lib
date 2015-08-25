###Asigning SNP allele 
###&AssignVariationAllele()
###Global: $debug, $verbose
###Dependancy: &ExtractAllele
sub AssignVariationAllele {
	my ($AVAaabbdd_vcfindex, $AVAaabb_vcf_index, $AVAaa_vcf_index, $AVAdd_vcf_index)=@_;

	my $AVAtest_AA_expressed=1;
	$AVAtest_AA_expressed=0 if (scalar(keys %{$AVAaa_vcf_index})==0);
	my $AVAtest_DD_expressed=1;
	$AVAtest_DD_expressed=0 if (scalar(keys %{$AVAdd_vcf_index})==0);
	my $AVAtest_bb_expressed=0;
#Format: %AVAgenome2allele=(chr => (pos => ('A' => allele, 'B' => Allele, 'D' => allele)))
	my %AVAgenome2allele=();
#Format: %AVAallele2genome=(chr => (pos => (allele => 'A', allele => 'B', allele => 'D')))
	my %AVAallele2genome=();
	foreach my $AVAchrom (keys %{$AVAaabbdd_vcfindex}) {
		foreach my $AVApos (sort {$a<=>$b} keys %{${$AVAaabbdd_vcfindex}{$AVAchrom}}) {
			print "SUB(AssignVariationAllele)Test: Chr: $AVAchrom\tPos: $AVApos\t Arr: @{${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}\n" if ($debug or $verbose);
			my ($AVAaabbdd_genohash_index, $AVAaabbdd_genoarr_index)=&ExtractAllele(${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[2], 1);
			if (scalar(@$AVAaabbdd_genoarr_index) ==1) {
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}=${$AVAaabbdd_genoarr_index}[0];
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}=${$AVAaabbdd_genoarr_index}[0];
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}=${$AVAaabbdd_genoarr_index}[0];
				${${$AVAallele2genome{$AVAchrom}}{$AVApos}}{${$AVAaabbdd_genoarr_index}[0]}='ABD';
			}
			else {
				my $AVAaa_allele='';
				my $AVAbb_allele='';
				my $AVAdd_allele='';
###COMMENT: defined dd allele
				if (exists ${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}) {
					my ($AVAdd_genohash_index, $AVAdd_genoarr_index)=&ExtractAllele(${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[2], 1);
					if (scalar(@{$AVAdd_genoarr_index})==1 and ${$AVAdd_genoarr_index}[0]=~/^\d+$/) {
						$AVAdd_allele=${$AVAdd_genoarr_index}[0];
						if (exists ${$AVAaabbdd_genohash_index}{$AVAdd_allele} and ${$AVAdd_genohash_index}{$AVAdd_allele} eq ${$AVAaabbdd_genohash_index}{$AVAdd_allele}) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}=$AVAdd_allele;
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{$AVAdd_allele}.='D';
						}
						else {
							print STDERR "SUB(AssignVariationAllele)Warning: non AABBDD allele $AVAdd_allele in DD at $AVAchrom:$AVApos\n";
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'});
##COMMENT: defined AA allele
				if (exists ${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}) {
					my ($AVAaa_genohash_index, $AVAaa_genoarr_index)=&ExtractAllele(${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[2], 1);
					if (scalar(@$AVAaa_genoarr_index)==1 and $AVAaa_allele=~/^\d+$/) {
						$AVAaa_allele=${$AVAaa_genoarr_index}[0];
						if (exists ${$AVAaabbdd_genohash_index}{$AVAaa_allele} and ${$AVAaa_genohash_index}{$AVAaa_allele} eq ${$AVAaabbdd_genohash_index}{$AVAaa_allele}) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}=$AVAaa_allele;
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{$AVAaa_allele}.='A';
						}
						else {
							print STDERR "SUB(AssignVariationAllele)Warning: non AABBDD allele $AVAaa_allele in AA at $AVAchrom:$AVApos\n";
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'});
##COMMENT: defined B allele---Difficult
				if ($AVAtest_AA_expressed==1 and $AVAaa_allele=~/^\d$/ and exists ${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}) {
					my ($AVAaabb_genohash_index, $AVAaabb_genoarr_index)=&ExtractAllele(${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[0], ${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[1], ${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[2], 1);
					if (scalar(@{$AVAaabb_genoarr_index})==2 and exists ${$AVAaabb_genohash_index}{$AVAaa_allele}) {
						foreach (@{$AVAaabb_genoarr_index}) {
							$AVAbb_allele=$_ unless ($AVAaa_allele eq $_);
						}
						if (exists ${$AVAaabbdd_genohash_index}{$AVAbb_allele} and ${$AVAaabb_genohash_index}{$AVAbb_allele} eq ${$AVAaabbdd_genohash_index}{$AVAbb_allele} and $AVAdd_allele ne $AVAbb_allele) {
							${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}=$AVAbb_allele;
							$AVAtest_bb_expressed=1;
						}
					}
				}
				${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}='?' unless (exists ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'});
##COMMENT: output debug info
				if ($debug) {
					print "chr: $AVAchrom\tPos: $AVApos\tRef: ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[0]\tVar: ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[1]\tAABBDDgeno ${${${$AVAaabbdd_vcfindex}{$AVAchrom}}{$AVApos}}[2]\tAABBgeno: ";
					if (exists ${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}) {
						print ${${${$AVAaabb_vcf_index}{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tAA: ";
					if (exists ${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}) {
						print ${${${$AVAaa_vcf_index}{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tDD: ";
					if (exists ${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}) {
						print ${${${$AVAdd_vcf_index}{$AVAchrom}}{$AVApos}}[2];
					}
					else {
						print '?';
					}
					print "\tRNAphased: ${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'A'}/${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'B'}/${${$AVAgenome2allele{$AVAchrom}}{$AVApos}}{'D'}\n";
				}
			}
		}
	}
	return (\%AVAgenome2allele, $AVAtest_bb_expressed);
}
