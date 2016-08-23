# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::MiscKit

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
package FuhaoPerl5Lib::MiscKit;
use strict;
use warnings;
use Exporter;
use Data::Dumper qw /Dumper/;
use Scalar::Util 'reftype';
use Cwd;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150603';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(IsReference IsZeroIn TestIntersect MaxLength FullDigit MergeRanges);
%EXPORT_TAGS = ( DEFAULT => [qw(IsReference IsZeroIn TestIntersect MaxLength FullDigit MergeRanges)],
                 ALL    => [qw(IsReference IsZeroIn TestIntersect MaxLength FullDigit MergeRanges)]);


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
	
	my $TIsubinfo='SUB(TestIntersect)';
	
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
### Global: $num_digit
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



### Mergeramge {2=>10, 9=>10} => {2=>15}
### 

sub MergeRanges {
	my $MRrange=shift;
	
	my $MRsubinfo='SUB(MergeRanges)';
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
		if ($MRi>=$MRlast_start and $MRi<=$MRlast_end) {
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

#$MiscKit_success=1;$MiscKit_failure=0;



1;
