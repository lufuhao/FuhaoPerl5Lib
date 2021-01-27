# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::MiscKit

=head1 SYNOPSIS

File operations

=head1 Requirements



=head1 DESCRIPTION

=over 2

=item FullDigit($num, $num_digit[default: 9])

    * Fullfill a number to a certain length: 
    * Add 0 before if < $num_digit: 12345 => 0000012345
    * Or return $num if >= digit: 1234567890 => 1234567890
    * Return: $string


=item GetCascadeList($list1in, $list2in, $list3out)

    * Get Cascade List
    * Note: Ignore lines starting with '#'
    * Return: 1=success, 0-failure
    * List1 input
         ID1	subID1	subID2	...
    * List2 input
         subID1	subsubID1	subsubID2	...
         subID2	subsubID3	subsubID4	...
    * List3 output
         ID1	subsubID1	subsubID2	subsubID3	subsubID4

=item IsReference($var)

    * Detect a reference type
    * Dependancy: Scalar::Util
    * Return: 0=NotRef, 1=ScalarRef, 2=ArrayRef, 3=HashRef, 4=Unknown

=item IsZeroIn($str)

    * Test VCF geno if value or values(delimited by /, example: 1/0/1) contain 0
    * Return: 0=yes, 1=No

=item ListMerger (out.sum, undefined_mark, In1, In2, ...)

    * Merge multiple (>=2) list
    * Return: 1=yes, 0=No
        ### In1
        Name1	Value1
        Name2	value2

        ### In2
        Name1	Value3
        Name4	value4

        ### out.sum
        #     	[In1] 	[IN2]
        Name1	Value1	Value3
        Name2	value2	[undef]
        Name4	[undef]	value4

=item MaxLength(@string_arr)

    * Return Maximum length given a array of strings
    * Return: max_length[INT]

=item MergeRanges($%ranges)

    * Merge overlapped ranges
    * Example: {2=>10, 9=>15, 16=>20} => {2=>20}
    * Return: (1/0, $%ranges)
        1=success, 0=failure

=item MrnaSort (@mRNA_IDs)

    * Sort mRNA by suffix number (.1 .2 .3) or alphabet order (.a .b .c)for one gene
    * Return: (@sorted_mRNA_IDs)

=item ReadConfig($config_file)

    * Read configure files
        # ignore line
        [block1]
        Param1      =        value1
        Param2      =        value2
        [block2]
        Param3      =        value3
    * Return: (1/0, $%ranges{block}{param})
        1/0: 1=success, 0=failure
        ${$ranges}{block1}->{param1}=value1;

=item TestIntersect($start1, $end1, $start2, $end2)

    * Test if two range overlaps
    * Example: GFF start1-end1: 1000-1005 GFF start2-end2: 1003-1004
    * Return: ?= error, 0=NOT overlap, 1=overlap

=item UrlEncode($string)

    * Encode a string to URL code
    * Return: $url_str

=item UrlDecode($url_str)

    * Decode URL code to a string
    * Return: $string

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
package FuhaoPerl5Lib::MiscKit;
use strict;
use warnings;
use Exporter;
use Data::Dumper qw /Dumper/;
use Scalar::Util 'reftype';
use Cwd;
use JSON;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20210127';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(FullDigit GetCascadeList IsReference IsZeroIn KeggJosnLoad ListMerger MaxLength MergeRanges MrnaSort ReadConfig TestIntersect UrlDecode UrlEncode);
%EXPORT_TAGS = ( DEFAULT => [qw(FullDigit GetCascadeList IsReference IsZeroIn KeggJosnLoad ListMerger MaxLength MergeRanges MrnaSort ReadConfig TestIntersect UrlDecode UrlEncode)],
                 ALL    => [qw(FullDigit GetCascadeList IsReference IsZeroIn KeggJosnLoad ListMerger MaxLength MergeRanges MrnaSort ReadConfig TestIntersect UrlDecode UrlEncode)]);


my $MiscKit_success=1;
my $MiscKit_failure=0;
my $MiscKit_debug=0; #checge to 1 if debug


### Detect a reference type
### &IsReference($var)
### Global: $MiscKit_debug
### Dependancy: Scalar::Util
### Note: Returncode (0=NotRef, 1=ScalarRef, 2=ArrayRef, 3=HashRef, 4=Unknown)
sub IsReference {
	my $x=shift;
	
	my $IRsubinfo='SUB(MiscKit::IsReference)';
	
	my $reftype = reftype $x;
	if ( ! defined $reftype ) {
    	print "${IRsubinfo}info: \$x is not a reference.\n" if ($MiscKit_debug);
    	return 0;
	}
	elsif ( $reftype eq 'HASH' ) {
    	print "hash\n" if ($MiscKit_debug);# do something with %$x
    	return 3;
#    	my @arr1=keys %$x;
#    	print scalar(@arr1)."\n";
	}
	elsif ( $reftype eq 'ARRAY' ) {
		print "arr\n" if ($MiscKit_debug);# do something with @$x
		return 2;
	}
	elsif ( $reftype eq 'SCALAR' ) {
		print "Scalar\n" if ($MiscKit_debug);# do something with $$x
		return 1;
	}
	else {
		print "Unknown\n" if ($MiscKit_debug);# do something else
		return 4;
	}
}



###Test if value or values(delimited by /) contain 0
###&IsZeroIn($str)
###Global:none
###Dependancy: $geno_delimiter
###Return: 0=yes, 1=No
sub IsZeroIn {
	my $IZIstr=shift;
	my @IZIarr=split(/\//, $IZIstr);
	my $ISIzero=1;
	foreach (@IZIarr) {
		$ISIzero=0 if ($_ ==0);
	}
	return $ISIzero;
}



### Test if two range overlaps: GFF1 start-end: 1000-1005 GFF2 start-end: 1003-1004
### &TestIntersect(start1, end1, start2, end2)
### ReturnCode: ?= error, 0=NOT overlap, 1=overlap
### Global:
### Dependency:
### Note: 
sub TestIntersect {
	my ($TIgroup1start, $TIgroup1end, $TIgroup2start, $TIgroup2end)=@_;
	
	my $TIsubinfo='SUB(MiscKit::TestIntersect)';
	
	unless (defined $TIgroup1start and $TIgroup1start=~/^\d+$/ and defined $TIgroup1end and $TIgroup1end=~/^\d+$/ and defined $TIgroup2start and $TIgroup2start=~/^\d+$/ and defined $TIgroup2end and $TIgroup2end=~/^\d+$/) {
		print STDERR "${TIsubinfo}Error: invalid numbers\n" if ($MiscKit_debug);
		if (1) {
			if (defined $TIgroup1start) {
				unless ($TIgroup1start=~/^\d+$/) {
					print STDERR "${TIsubinfo}Error: 1-1 ($TIgroup1start) length".length($TIgroup1start)."\n";
					print STDERR "${TIsubinfo}Error: 1-1 ($TIgroup1start) NOT number\n";
				}
			}
			else {
				print STDERR "${TIsubinfo}Error: 1-1 NOT defined\n";
			}
			if (defined $TIgroup1end) {
				unless ($TIgroup1end=~/^\d+$/) {
					print STDERR "${TIsubinfo}Error: 1-2 ($TIgroup1end) length".length($TIgroup1end)."\n";
					print STDERR "${TIsubinfo}Error: 1-2 ($TIgroup1end) NOT number\n";
				}
			}
			else {
				print STDERR "${TIsubinfo}Error: 1-2 NOT defined\n";
			}
			if (defined $TIgroup2start) {
				unless ($TIgroup2start=~/^\d+$/) {
					print STDERR "${TIsubinfo}Error: 2-1 ($TIgroup2start) length".length($TIgroup2start)."\n";
					print STDERR "${TIsubinfo}Error: 2-1 ($TIgroup2start) NOT number\n";
				}
			}
			else {
				print STDERR "${TIsubinfo}Error: 2-1 NOT defined\n";
			}
			if (defined $TIgroup2end) {
				unless ($TIgroup2end=~/^\d+$/) {
					print STDERR "${TIsubinfo}Error: 2-2 ($TIgroup2end) length".length($TIgroup2end)."\n";
					print STDERR "${TIsubinfo}Error: 2-2 ($TIgroup2end) NOT number\n";
				}
			}
			else {
				print STDERR "${TIsubinfo}Error: 2-2 NOT defined\n";
			}
		}
		return '?';
	}
	($TIgroup1start, $TIgroup1end)=sort {$a <=> $b} ($TIgroup1start, $TIgroup1end);
	($TIgroup2start, $TIgroup2end)=sort {$a <=> $b} ($TIgroup2start, $TIgroup2end);
	
	if ($TIgroup1start >=$TIgroup2start and $TIgroup1start <= $TIgroup2end) {
		return 1;
	}
	elsif ($TIgroup2start >=$TIgroup1start and $TIgroup2start <=$TIgroup1end) {
		return 1;
	}
	elsif ($TIgroup1start>=$TIgroup2start and $TIgroup1end <=$TIgroup2end) {
		return 1;
	}
	elsif ($TIgroup1start<=$TIgroup2start and $TIgroup1end >=$TIgroup2end) {
		return 1;
	}
	else {
		return 0;
	}
}



### return Maximum length given a array of strings
### MaxLength(@string)
### Global: 
### Dependency: 
### Note:
### Return: max_length
sub MaxLength {
	my @MLarr=@_;
	my @MLlen=map {length($_)} @MLarr;
	my $MLmax=shift @MLlen;
	foreach (@MLlen) {
		$MLmax=$_ if ($_>$MLmax);
	}
	return $MLmax;
}



### Fullfill a number to a certain length
### &FullDigit(number, digit[9])
### Global: $full_digit
sub FullDigit {
	my ($FDnum, $FDdigit)=@_;

	my $FDsubinfo='SUB(MiscKit::FullDigit)';
	my $FD_8D='';
	$FDdigit=9 unless (defined $FDdigit);

	unless (defined $FDdigit and $FDdigit>0) {
		print STDERR "${$FDsubinfo}Error: Wrong num_digit detected\n" 
	}

	if (length($FDnum)<$FDdigit) {
		$FD_8D='0' x ($FDdigit-length($FDnum)).$FDnum;
	}
	else {
		$FD_8D=$FDnum;
	}
	return $FD_8D;
}



### Merge ramges
### &MergeRanges($%ranges)
### Example: {2=>10, 9=>10} => {2=>15}
### Return (1/0, $%ranges)
sub MergeRanges {
	my $MRrange=shift;
	
	my $MRsubinfo='SUB(MiscKit::MergeRanges)';
	my %MRreturnhash=();
	my $MRtest=0;
	my ($MRlast_start, $MRlast_end);
	
	if (0) {### For test ###
		print $MRsubinfo, "Test: input hash\n";
		print Dumper $MRrange;
		print "\n";
	}
	
	foreach my $MRi (sort {$a<=>$b} keys %{$MRrange}) {
		unless ($MRi=~/^\d+$/ and ${$MRrange}{$MRi}=~/^\d+$/ and $MRi<=${$MRrange}{$MRi}) {
			print STDERR $MRsubinfo, "Warnings: unknown/invalid number pairs: $MRi => ${$MRrange}{$MRi}\n";
			return $MiscKit_failure;
		}
		if ($MRtest==0) {
			$MRlast_start=$MRi;
			$MRlast_end=${$MRrange}{$MRi};
			$MRtest++;
			next;
		}
		if ($MRi>=$MRlast_start and $MRi<=($MRlast_end+1)) {
			if (${$MRrange}{$MRi}<=$MRlast_end) {
#				next;
			}
			elsif (${$MRrange}{$MRi}>$MRlast_end) {
				$MRlast_end=${$MRrange}{$MRi};
#				next;
			}
			else {
				print STDERR $MRsubinfo, "Warnings2: unknown/invalid number pairs: $MRi => ${$MRrange}{$MRi}\n";
				return $MiscKit_failure;
			}
		}
		elsif ($MRi> $MRlast_end) {
			$MRreturnhash{$MRlast_start}=$MRlast_end;
			if (0) {### For test ###
				print $MRsubinfo, "Test: Add into hash: $MRlast_start => $MRlast_end\n";
			}
			$MRlast_start=$MRi;
			$MRlast_end=${$MRrange}{$MRi};
		}
		if (0) {### For test ###
			print $MRsubinfo, "Test: This: $MRi => ${$MRrange}{$MRi}\n";
			print $MRsubinfo, "Last: This: $MRlast_start => $MRlast_end\n";
		}
	}
	$MRreturnhash{$MRlast_start}=$MRlast_end;
	return ($MiscKit_success, \%MRreturnhash);
}



### Encode a string to URL code
### UrlEncode($string)
### Return: $url_str
sub UrlEncode {
	my $UEstr = shift;
	$UEstr =~ s/ /+/g;
	$UEstr =~ s/([^A-Za-z0-9\+-])/sprintf("%%%02X", ord($1))/seg;
	return $UEstr;
}
### Decode URL code to a string
### UrlDecode($url_str)
### Return: $string
sub UrlDecode {
	my $UDstr = shift;
	$UDstr =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
	$UDstr =~ s/\+/ /g;
	return $UDstr;
}
sub url_encode {
	my $rv = shift;
	$rv =~ s/([^a-z\d\Q.-_~ \E])/sprintf("%%%2.2X", ord($1))/geix;
	$rv =~ tr/ /+/;
	return $rv;
}
sub url_decode {
	my $rv = shift;
	$rv =~ tr/+/ /;
	$rv =~ s/\%([a-f\d]{2})/ pack 'C', hex $1 /geix;
	return $rv;
}



### Get Cascade List
### GetCascadeList($list1in, $list2in, $list3out)
### Note: Ignore lines starting with '#'
### Return: 1=success, 0-failure
###list1 input
#ID1	subID1	subID2	...

###List2
#subID1	subsubID1	subsubID2	...
#subID2	subsubID3	subsubID4	...

### List3 output
#ID1	subsubID1	subsubID2	subsubID3	subsubID4
sub GetCascadeList {
	my ($GCLlist1, $GCLlist2, $GCLlist3)=@_;

	my $GCLsubinfo='SUB(MergeRanges)';
	local *GCLLIST1; local *GCLLIST2; local *GCLLIST3;
	my %GCLlist1hash=();
	my %GCLlist2hash=();
	my $GCLlinenum=0;
	my $GCLlinesread=0;

	unless (defined $GCLlist1 and -s $GCLlist1) {
		print STDERR $GCLsubinfo, "Error: invalid list1\n";
		return $MiscKit_failure;
	}
	unless (defined $GCLlist2 and -s $GCLlist2) {
		print STDERR $GCLsubinfo, "Error: invalid list2\n";
		return $MiscKit_failure;
	}
	unless (defined $GCLlist3) {
		print STDERR $GCLsubinfo, "Error: invalid list3\n";
		return $MiscKit_failure;
	}
	unlink $GCLlist3 if (-e $GCLlist3);

	close GCLLIST1 if (defined fileno(GCLLIST1));
	unless (open(GCLLIST1, " < $GCLlist1")) {
		print STDERR $GCLsubinfo, "Error: can not open list1: $GCLlist1\n";
		return $MiscKit_failure;
	}
	while (my $GCLline1=<GCLLIST1>) {
		chomp $GCLline1;
		$GCLlinenum++;
		next if ($GCLline1=~/^#/);
		my @GCLarr1=split(/\t/, $GCLline1);
		unless (scalar(@GCLarr1)>1 and defined $GCLarr1[1] and $GCLarr1[1]=~/\S+/) {
			print STDERR $GCLsubinfo, "Warnings: List1 problematic line($GCLlinenum): $GCLline1\n";
			next;
		}
		$GCLlinesread++;
		my $GCLid1=shift @GCLarr1;
		foreach (@GCLarr1) {
			$GCLlist1hash{$GCLid1}{$_}++;
		}
	}
	close GCLLIST1;
	print $GCLsubinfo, "#### SUMMARY: list1 ######\n";
	print $GCLsubinfo, "#### Total lines: $GCLlinenum\n";
	print $GCLsubinfo, "#### Valid lines: $GCLlinesread\n\n";

	$GCLlinenum=0;
	$GCLlinesread=0;

	close GCLLIST2 if (defined fileno(GCLLIST2));
	unless (open(GCLLIST2, " < $GCLlist2")) {
		print STDERR $GCLsubinfo, "Error: can not open list2: $GCLlist2\n";
		return $MiscKit_failure;
	}
	while (my $GCLline2=<GCLLIST2>) {
		chomp $GCLline2;
		$GCLlinenum++;
		next if ($GCLline2=~/^#/);
		my @GCLarr2=split(/\t/, $GCLline2);
		unless (scalar(@GCLarr2)>1 and defined $GCLarr2[1] and $GCLarr2[1]=~/\S+/) {
			print STDERR $GCLsubinfo, "Warnings: List2 problematic line($GCLlinenum): $GCLline2\n";
			next;
		}
		$GCLlinesread++;
		my $GCLid2=shift @GCLarr2;
		foreach (@GCLarr2) {
			$GCLlist2hash{$GCLid2}{$_}++;
		}
	}
	close GCLLIST2;
	print $GCLsubinfo, "#### SUMMARY: list1 ######\n";
	print $GCLsubinfo, "#### Total lines: $GCLlinenum\n";
	print $GCLsubinfo, "#### Valid lines: $GCLlinesread\n\n";

	$GCLlinenum=0;
	close GCLLIST3 if (defined fileno(GCLLIST3));
	unless (open (GCLLIST3, " > $GCLlist3")) {
		print STDERR $GCLsubinfo, "Error: can not write list3: $GCLlist3\n";
		return $MiscKit_failure;
	}
	foreach my $GCLid3 (sort keys %GCLlist1hash) {
		print GCLLIST3 "$GCLid3";
		$GCLlinenum++;
		my $GCLlinesread1=0;
		my %GCLprintednames=();
		foreach my $GCLid4 (sort keys %{$GCLlist1hash{$GCLid3}}) {
			my $GCLlinesread2=0;
			if (exists $GCLlist2hash{$GCLid4}) {
				foreach my $GCLid5 (sort keys %{$GCLlist2hash{$GCLid4}}) {
					$GCLlinesread1++;
					$GCLlinesread2++;
					unless (exists $GCLprintednames{$GCLid5}) {### control non redundant
						print GCLLIST3 "\t", $GCLid5;
					}
					$GCLprintednames{$GCLid5}++;
				}
			}
			else {
				print STDERR $GCLsubinfo, "Warnings: List1 subID not exists in List2: $GCLid4\n";
			}
			if ($GCLlinesread2==0) {
				print STDERR $GCLsubinfo, "Warnings: List1 subID empty in List2: $GCLid4\n";
			}
		}
		print GCLLIST3 "\n";
		if ($GCLlinesread1==0) {
			print STDERR $GCLsubinfo, "Warnings: List1 ID empty: $GCLid3\n";
		}
	}
	close GCLLIST3;

	return $MiscKit_success;
}




### Merge multiple (>=2) list
### ListMerger (out.sum, undefined_mark, $arr_index={In1, In2, ...})

### In1
#Name1	Value1
#Name2	value2

### In2
#Name1	Value3
#Name4	value4

### out.sum
#     	[In1] 	[IN2]
#Name1	Value1	Value3
#Name2	value2	[undef]
#Name4	[undef]	value4
sub ListMergersingleCol {
	my ($LMout, $LMundef_mark, $LMinputs)=@_;
	
	my $MRsubinfo='SUB(MiscKit::ListMerger)';
	my $LMnumlines=0;
	my %LMhash=();
	$LMundef_mark="undef" unless (defined $LMundef_mark);
	local *LMSUMOUT; local *LMINPUT;
	
	
	foreach my $LMindin (@{$LMinputs}) {
		unless (defined $LMindin and -s $LMindin) {
			print STDERR $MRsubinfo, "Error: invalid file $LMindin\n";
			return $MiscKit_failure;
		}
	}
	
	for (my $LMindex=0; $LMindex<scalar(@{$LMinputs}); $LMindex++) {
		close LMINPUT if (defined fileno(LMINPUT));
		unless (open(LMINPUT, '<', ${$LMinputs}[$LMindex])) {
			print STDERR $MRsubinfo, "Error: can not open ${$LMinputs}[$LMindex]\n";
			return $MiscKit_failure;
		}
		my %LMdup=();
		while (my $LMline=<LMINPUT>) {
			chomp $LMline;
			next if ($LMline=~/^#/);
			my @LMarr=split(/\t/, $LMline);
			
			if (exists $LMdup{$LMarr[0]}) {
				print STDERR $MRsubinfo, "Error: duplicated ID $LMline\n";
				return $MiscKit_failure;
			}
			$LMhash{$LMarr[0]}[$LMindex]=$LMarr[1];
		}
		close LMINPUT;
	}
	
	close LMSUMOUT if (defined fileno(LMSUMOUT));
	unless (open LMSUMOUT, ">", $LMout) {
		print STDERR $MRsubinfo, "Error: can not write output file $LMout\n";
		return $MiscKit_failure;
	}
	print LMSUMOUT "#IDs\t", join("\t", @{$LMinputs}), "\n";
	foreach my $LMindex2 (sort keys %LMhash) {
		my @LMArr2=@{$LMhash{$LMindex2}};
		for (my $LMindex3=0; $LMindex3<scalar(@{$LMinputs});$LMindex3++) {
			if (! defined $LMArr2[$LMindex3]) {
				$LMArr2[$LMindex3]=$LMundef_mark;
			}
#			elsif ($LMArr2[$LMindex3] eq 'NaN') {
#				$LMArr2[$LMindex3]=$LMundef_mark;
#			}
		}
		print LMSUMOUT $LMindex2, "\t", join("\t", @LMArr2), "\n";
	}
	close LMSUMOUT;
	
	return $MiscKit_success;
}
### sipport multiple columns
sub ListMerger {
	my ($LMout, $LMundef_mark, $LMinputs)=@_;
	
	my $LMsubinfo='SUB(MiscKit::ListMergerMultiCols)';
	my $LMnumlines=0;
	my %LMhash=();
	my %LMlength=();
	my @LMheaders=("#IDs");
	$LMundef_mark="undef" unless (defined $LMundef_mark);
	local *LMSUMOUT; local *LMINPUT;
	
	foreach my $LMindin (@{$LMinputs}) {
		unless (defined $LMindin and -s $LMindin) {
			print STDERR $LMsubinfo, "Error: invalid file $LMindin\n";
			return $MiscKit_failure;
		}
	}
	
	for (my $LMindex=0; $LMindex<scalar(@{$LMinputs}); $LMindex++) {
		$LMnumlines=0;
		close LMINPUT if (defined fileno(LMINPUT));
		unless (open(LMINPUT, '<', ${$LMinputs}[$LMindex])) {
			print STDERR $LMsubinfo, "Error: can not open ${$LMinputs}[$LMindex]\n";
			return $MiscKit_failure;
		}
		my %LMdup=();
		while (my $LMline=<LMINPUT>) {
			chomp $LMline; $LMline=~s/[\r\n]$//; 
			$LMnumlines++;
			next if ($LMline=~/^#/);
			my @LMarr=split(/\t/, $LMline);
			
			if (exists $LMdup{$LMarr[0]}) {
				print STDERR $LMsubinfo, "Error: duplicated ID $LMline\n";
				return $MiscKit_failure;
			}
			my $LMidentifier=shift @LMarr;
			if ((! exists $LMlength{$LMindex}) or ($LMlength{$LMindex} < scalar(@LMarr))) {
					$LMlength{$LMindex} = scalar(@LMarr);
			}
			@{$LMhash{$LMidentifier}{$LMindex}}=@LMarr;
		}
		close LMINPUT;
		
		print $LMsubinfo, "Info:  $LMnumlines  lines: ${$LMinputs}[$LMindex]\n";
		
		if ((! exists $LMlength{$LMindex} or ($LMlength{$LMindex}==0))) {
			print STDERR "Error: empty elements for ${$LMinputs}[$LMindex]\n";
			return $MiscKit_failure;
		}
		elsif ($LMlength{$LMindex}==1) {
			push (@LMheaders, ${$LMinputs}[$LMindex]);
		}
		elsif ($LMlength{$LMindex}>1) {
			push (@LMheaders, ${$LMinputs}[$LMindex]);
			@LMheaders=(@LMheaders, ("") x ($LMlength{$LMindex}-1));
		}
	}
	
	print "Info: Total ", scalar(keys %LMhash), "\n";
	
	close LMSUMOUT if (defined fileno(LMSUMOUT));
	unless (open LMSUMOUT, ">", $LMout) {
		print STDERR $LMsubinfo, "Error: can not write output file $LMout\n";
		return $MiscKit_failure;
	}
	print LMSUMOUT join("\t", @LMheaders), "\n";
	foreach my $LMindex2 (sort keys %LMhash) {
		my @LMArr2=();
		my @LMarr_out=($LMindex2);
		for (my $LMindex3=0; $LMindex3<scalar(@{$LMinputs});$LMindex3++) {
			if (! exists $LMhash{$LMindex2}{$LMindex3}) {
				@LMArr2=($LMundef_mark) x $LMlength{$LMindex3};
			}
			else {
				@LMArr2=(@{$LMhash{$LMindex2}{$LMindex3}}, ($LMundef_mark) x ($LMlength{$LMindex3}-scalar(@{$LMhash{$LMindex2}{$LMindex3}})));
			}
			@LMarr_out=(@LMarr_out, @LMArr2);
		}
		print LMSUMOUT join("\t", @LMarr_out), "\n";
	}
	
	close LMSUMOUT;
	
	return $MiscKit_success;
}


### Read configure files
### ReadConfig($config_file)
### Global:
### Dependency:
### Note:
### $config_file
### # ignore line
### [block1]
### Param1      =        value1
### Param2      =        value2
### [block2]
### Param3      =        value3
sub ReadConfig {
	my $RCconf=shift;
	
	my $RCsubinfo='SUB(MiscKit::ReadConfig)';
	my $RCret_config={};
	my $RCblock='';
	local *RCCONFIG;
	
	unless (defined $RCconf and -e $RCconf) {
		print STDERR $RCsubinfo,"Error: invalid configure file\n";
		return $MiscKit_failure;
	}
	close RCCONFIG if (defined fileno(RCCONFIG));
	unless (open RCCONFIG, '<', $RCconf) {
		print STDERR $RCsubinfo,"Error: can not open configure file\n";
		return $MiscKit_failure;
	}
	while (my $RCline=<RCCONFIG>) {
		chomp $RCline;
		$RCline =~ s/^\s+//g;
		$RCline =~ s/\s+$//g;
# skipping comments and empty lines:
		next unless ($RCline =~ /\S+/);
		next if ($RCline =~ /^(\n|\#|;)/);
# parsing the block name:
		if ($RCline =~ /^\s*\[\s*([^\]]+)\s*\]$/) {
			$RCblock=lc($1);
			next;
		}
# parsing key/value pairs
		if ($RCline =~ /^\s*([^=]*\w)\s*=\s*(.*)\s*$/) {
			${$RCret_config}{$RCblock}->{lc($1)}=$2;
			next;
		}
# if we came this far, the syntax couldn't be validated:
		print STDERR $RCsubinfo, "Warnings: syntax error on line($.): \n",$RCline,"\n";
	}
	close RCCONFIG;
	print $RCsubinfo, "Info: parameter set:\n";
	foreach my $RCind_block (sort keys %{$RCret_config}) {
		print "  BLOCK: $RCind_block\n";
		foreach my $RCind_param (sort keys %{${$RCret_config}{$RCind_block}}) {
			print "    $RCind_param  =  ", ${$RCret_config}{$RCind_block}->{$RCind_param},"\n";
		}
	}
	print "\n";
	return ($MiscKit_success, $RCret_config);
}



### Sort mRNA by suffix number or alphabet order for one gene
### MrnaSort (@mRNA_IDs)
### Return: (@sorted_mRNA_IDs)
### Global:
### Dependency:
### Note:
sub MrnaSort {
	my @MSarr_in=@_;
	
	my $MStest_num=1;
	my %MSorder=();
	my @MSarr_out=();
	
	foreach my $MSind_mrna (@MSarr_in) {
		if ($MSind_mrna=~/\.(\d+)$/) {
			my $MSind_num=$1;
			if (exists $MSorder{$MSind_num}) {
				$MStest_num=0;last;
			}
			$MSorder{$MSind_num}=$MSind_mrna;
		}
		else {
			$MStest_num=0;last;
		}
	}
	if ($MStest_num==0) {
		@MSarr_out=sort @MSarr_in;
	}
	else {
		foreach (sort {$a<=>$b} keys %MSorder) {
			push (@MSarr_out, $MSorder{$_})
		}
	}
	
	return @MSarr_out;
}



### return a duplicated N times string
sub DuplicateNString {
	my ($DNSstr, $DNSsep, $DNSnum)=@_;
	
	my $DNSret='';
	my @DNAarr=();
	
	for (my $DNSx=0; $DNSx<$DNSnum; $DNSx++) {
		push (@DNAarr, $DNSstr);
	}
	
	return join($DNSsep, @DNAarr);
}


#$MiscKit_success=1;$MiscKit_failure=0;
#my $MRsubinfo='SUB(MiscKit::MergeRanges)';



sub JsonOpen {
	my $JOjson=shift;
	
	my $JOsubinfo="SUB(MiscKit::JsonOpen)";
	my $JOjsontext="";
	use utf8;
	
	unless (-s $JOjson) {
		print STDERR $JOsubinfo, "Error: invalid json file\n";
		exit 100;
	}
	local $/=undef;
	if (open my $JOfh, " < ", $JOjson) {
		$JOjsontext = <$JOfh>;
		close $JOfh;
	}
	
	return $JOjsontext;
}
sub JsonLoad {
	my $JLjson=shift;
	
	my $JLtext=JsonOpen($JLjson);
	
	my $JLdecode=JSON->new->utf8->decode($JLtext);
#	print Dumper $JLdecode; ### for Test ###
	
	return $JLdecode;
}
sub getJsonName {
	my $JNtext=shift;
	
	my $JNid='NA';
	my $JNannot='NA';
	
	if ($JNtext=~/^(\d+)\s+(.+)$/) {
		$JNid=$1;
		$JNannot=$2;
	}
	elsif ($JNtext=~/^(ko\d+)$/) {
		$JNid=$1;
	}
	elsif ($JNtext=~/^(K\d+)\s*(.*)$/) {
		$JNid=$1;
		$JNannot=$2;
	}
	else {
		print STDERR "Warnings: unknown Name: $JNtext\n";
		return ('NA', 'NA');
	}
	
	return ($JNid, $JNtext);
}
sub KeggJosnLoad {
	my $KJLjson=shift;
	
	my $KJLk2name={};
	my $KJLko2pathway={};
	
	my $KJLsubinfo="SUB(MiscKit::KeggJosnLoad)";
	my $KJLtext=JsonLoad($KJLjson);
	
	unless (exists ${$KJLtext}{'children'}) {
		print STDERR $KJLsubinfo, "Error: invalid KEGG JSON file\n";
		exit 100;
	}
	my ($KJLx1)=getJsonName(${$KJLtext}{'name'});
	for (my $KJLa=0; $KJLa<scalar(@{${$KJLtext}{'children'}}); $KJLa++) {
		my ($KJLx2, $KJLy2)=("NA", "NA");
		if (exists ${$KJLtext}{'children'}[$KJLa]{'name'}) {
			($KJLx2, $KJLy2)=getJsonName(${$KJLtext}{'children'}[$KJLa]{'name'});
			$KJLx2="ko".$KJLx2;
			if (exists ${$KJLk2name}{"$KJLx2"}) {
				print STDERR $KJLsubinfo, "Error: duplicated Level2 pathway ID: $KJLx2\n";
			}
			else {
				${$KJLk2name}{"$KJLx2"}=$KJLy2;
			}
		}
		unless(exists ${$KJLtext}{'children'}[$KJLa]{'children'}) {
			print STDERR $KJLsubinfo, "Warnings: 'children' $KJLa [name: $KJLx2] has no children\n";
			next;
		}
		my $KJLchildren2=${$KJLtext}{'children'}[$KJLa]{'children'};
		for (my $KJLb=0; $KJLb<scalar(@{$KJLchildren2}); $KJLb++) {
			my ($KJLx3, $KJLy3)=("NA", "NA");
			if (exists ${$KJLchildren2}[$KJLb]{'name'}) {
				($KJLx3, $KJLy3)=getJsonName(${$KJLchildren2}[$KJLb]{'name'});
				$KJLx3="ko".$KJLx3;
				if (exists ${$KJLk2name}{"$KJLx3"}) {
					print STDERR $KJLsubinfo, "Error: duplicated Level3 pathway ID: $KJLx3\n";
				}
				else {
					${$KJLk2name}{"$KJLx3"}=$KJLy3;
				}
				
			}
			unless(exists ${$KJLchildren2}[$KJLb]{'children'}) {
				print STDERR $KJLsubinfo, "Warnings: 'children' $KJLa 'children' $KJLb [name: $KJLx3] has no children\n";
				next;
			}
			my $KJLchildren3=${$KJLchildren2}[$KJLb]{'children'};
			for (my $KJLc=0; $KJLc<scalar(@{$KJLchildren3}); $KJLc++) {
				my ($KJLx4, $KJLy4)=("NA", "NA");
				if (exists ${$KJLchildren3}[$KJLc]{'name'}) {
					($KJLx4, $KJLy4)=getJsonName(${$KJLchildren3}[$KJLc]{'name'});
					$KJLx4="ko".$KJLx4;
					if (exists ${$KJLk2name}{"$KJLx4"}) {
						print STDERR $KJLsubinfo, "Error: duplicated Level4 pathway ID: $KJLx4\n";
					}
					else {
						${$KJLk2name}{"$KJLx4"}=$KJLy4;
					}
				}
				unless(exists ${$KJLchildren3}[$KJLc]{'children'}) {
					print STDERR $KJLsubinfo, "Warnings: 'children' $KJLa 'children' $KJLb 'children' $KJLc [name: $KJLx4] has no name/children\n";
					next;
				}
				my $KJLchildren4=${$KJLchildren3}[$KJLc]{'children'};
				for (my $KJLd=0; $KJLd<scalar(@{$KJLchildren4}); $KJLd++) {
					unless(exists ${$KJLchildren4}[$KJLd]{'name'}) {
						print STDERR $KJLsubinfo, "Warnings: 'children' $KJLa 'children' $KJLb 'children' $KJLc 'children' $KJLd has no name\n";
						next;
					}
					my ($KJLx5, $KJLy5)=getJsonName(${$KJLchildren4}[$KJLd]{'name'});
					unless (exists ${$KJLk2name}{"$KJLx5"}) {
						${$KJLk2name}{"$KJLx5"}=$KJLy5;
					}
					${$KJLko2pathway}{$KJLx5}{$KJLx4}++;
				}
			}
		}
	}
#	print Dumper $KJLk2name;
#	print Dumper $KJLko2pathway;
	
	return ($KJLko2pathway)
}



1;
