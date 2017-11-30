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

=item MaxLength(@string_arr)

    * Return Maximum length given a array of strings
    * Return: max_length[INT]

=item MergeRanges($%ranges)

    * Merge overlapped ranges
    * Example: {2=>10, 9=>15, 16=>20} => {2=>20}
    * Return: (1/0, $%ranges)
    * 1=success, 0=failure

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

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20171130';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(IsReference IsZeroIn TestIntersect MaxLength FullDigit MergeRanges UrlEncode UrlDecode GetCascadeList);
%EXPORT_TAGS = ( DEFAULT => [qw(IsReference IsZeroIn TestIntersect MaxLength FullDigit MergeRanges UrlEncode UrlDecode GetCascadeList)],
                 ALL    => [qw(IsReference IsZeroIn TestIntersect MaxLength FullDigit MergeRanges UrlEncode UrlDecode GetCascadeList)]);


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
		foreach my $GCLid4 (sort keys %{$GCLlist1hash{$GCLid3}}) {
			my $GCLlinesread2=0;
			if (exists $GCLlist2hash{$GCLid4}) {
				foreach my $GCLid5 (sort keys %{$GCLlist2hash{$GCLid4}}) {
					$GCLlinesread1++;
					$GCLlinesread2++;
					print GCLLIST3 "\t", $GCLid5;
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


#$MiscKit_success=1;$MiscKit_failure=0;
#my $MRsubinfo='SUB(MiscKit::MergeRanges)';



1;
