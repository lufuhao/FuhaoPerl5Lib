# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::StatsKit

=head1 SYNOPSIS

Statistics perl

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
package FuhaoPerl5Lib::StatsKit;
use strict;
use warnings;
use Exporter;
use Data::Dumper qw/Dumper/;
use FuhaoPerl5Lib::CmdKit;
use Statistics::Basic qw/:all/;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150603';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw(CallCorrelation OLD_CallCorrelation);
%EXPORT_TAGS = ( DEFAULT => [qw(&CallCorrelation)],
                 ALL    => [qw(&CallCorrelation)]);

my $StatsKit_success=1;
my $StatsKit_failure=0;
my $debug=0;


### calculate Statistics::Basic correlation
### &CalCorr(\@arr1, \@arr2)
### Global: $verbose, $debug
### Dependency: Statistics::Basic
### Note:
sub OLD_CallCorrelation {
	my ($CCv1_arrindex, $CCv2_arrindex)=@_;
	
	my $CCsubinfo='SUB(StatsKit::CallCorrelation)';
	
	if (scalar(@{$CCv1_arrindex}) <1 or scalar(@{$CCv1_arrindex}) != scalar(@{$CCv2_arrindex})) {
		print STDERR "${CCsubinfo}Error: empty list for correlation\n";
		return $StatsKit_failure;
	}
	my $CCarr1=[];
	my $CCarr2=[];
	for (my $CCi=0; $CCi<scalar(@{$CCv1_arrindex}); $CCi++) {
		if (defined ${$CCv1_arrindex}[$CCi] and ${$CCv1_arrindex}[$CCi] =~/^\d+$/ and defined ${$CCv2_arrindex}[$CCi] and ${$CCv2_arrindex}[$CCi] =~/^\d+$/){
			if (${$CCv1_arrindex}[$CCi]==0) {
				push (@{$CCarr1}, 10);
			}
			else {
				push (@{$CCarr1}, ${$CCv1_arrindex}[$CCi]);
			}
			if (${$CCv2_arrindex}[$CCi]==0) {
				push (@{$CCarr2}, 10);
			}
			else {
				push (@{$CCarr2}, ${$CCv2_arrindex}[$CCi]);
			}
		}
	}
	my $CCVelements=0;
	if (scalar(@{$CCarr1}) == scalar(@{$CCarr2})) {
		$CCVelements=scalar(@{$CCarr1});
	}
	else {
		print STDERR "${CCsubinfo}Error: uneven array list\n";
		return $StatsKit_failure;
	}
	return $StatsKit_failure if ($CCVelements==0);
#	undef $CCv1_arrindex;
#	undef $CCv2_arrindex;
#	my $CCv1=vector(@CCarr1); 
#	my $CCv2=vector(@CCarr2);
#	my ($CCfilv1, $CCfilv2)=handle_missing_values($CCv1, $CCv2);
#	my $CCcov = covariance ($CCfilv1, $CCfilv2);
#	my $CCcor = correlation ($CCfilv1, $CCfilv2);
	print Dumper $CCarr1;
	print Dumper $CCarr2;
	my $CCcov = covariance ($CCarr1, $CCarr2);
	my $CCcor = correlation ($CCarr1, $CCarr2);
	print "SUB(CalCorr)Info: ###Covariance: $CCcov\t### Correlation: $CCcor\n";# if ($debug);
	return ($StatsKit_success, $CCcov, $CCcor, $CCVelements);
}


sub CallMean {
	my $x=shift;
	my $num = scalar(@{$x});
	my $sum_x = '0';
	for (my $i = 0; $i < scalar(@{$x}); $i++){
		$sum_x+= $x->[$i];
	}
	my $mu_x = $sum_x/$num;
	
	return($mu_x);
}
 
### ss = sum of squared deviations to the mean
sub CallSqDev {
	my ($x,$y, $mean_x,$mean_y)=@_;
	my $sum = '0';
	for (my $i=0;$i<scalar(@{$x});$i++){
		$sum += ($x->[$i]-$mean_x)*($y->[$i]-$mean_y);
	}
	return $sum;
}
 
sub CallCorrelation {
	my ($CCv1_arrindex, $CCv2_arrindex)=@_;
	
	my $CCsubinfo='SUB(StatsKit::CallCorrelation)';
	
	if (scalar(@{$CCv1_arrindex}) <1 or scalar(@{$CCv1_arrindex}) != scalar(@{$CCv2_arrindex})) {
		print STDERR "${CCsubinfo}Error: empty list for correlation\n";
		return $StatsKit_failure;
	}
	my $CCarr1=[];
	my $CCarr2=[];
	my %CCarcount1=();
	my %CCarcount2=();
	for (my $CCi=0; $CCi<scalar(@{$CCv1_arrindex}); $CCi++) {
		if (defined ${$CCv1_arrindex}[$CCi] and ${$CCv1_arrindex}[$CCi] =~/^\d+$/ and defined ${$CCv2_arrindex}[$CCi] and ${$CCv2_arrindex}[$CCi] =~/^\d+$/){
			if (${$CCv1_arrindex}[$CCi]==0) {
				push (@{$CCarr1}, 10);
				$CCarcount1{10}++;
			}
			else {
				push (@{$CCarr1}, ${$CCv1_arrindex}[$CCi]);
				$CCarcount1{${$CCv1_arrindex}[$CCi]}++;
			}
			if (${$CCv2_arrindex}[$CCi]==0) {
				push (@{$CCarr2}, 10);
				$CCarcount2{10}++;
			}
			else {
				push (@{$CCarr2}, ${$CCv2_arrindex}[$CCi]);
				$CCarcount2{${$CCv2_arrindex}[$CCi]}++;
			}
		}
	}
	my $CCVelements=0;
	if (scalar(@{$CCarr1}) == scalar(@{$CCarr2})) {
		$CCVelements=scalar(@{$CCarr1});
	}
	else {
		print STDERR "${CCsubinfo}Error: uneven array list\n";
		return $StatsKit_failure;
	}
	return $StatsKit_failure if ($CCVelements==0);
	push (@{$CCarr1}, 20);
	push (@{$CCarr2}, 20);
	
	my $mean_x=&CallMean($CCarr1);
	my $mean_y=&CallMean($CCarr2);
	my $ssxx=&CallSqDev($CCarr1, $CCarr1, $mean_x, $mean_x);
	my $ssyy=&CallSqDev($CCarr2, $CCarr2, $mean_y, $mean_y);
	my $ssxy=&CallSqDev($CCarr1, $CCarr2, $mean_x, $mean_y);
	if ($ssxy==0) {
		return ($StatsKit_success, 1, 0, $CCVelements);
	}
	my $cov=$ssxy/$CCVelements;
	my $correl=&CallCorrel($ssxx,$ssyy,$ssxy);
	my $xcorrel=sprintf("%.4f",$correl);
#	print Dumper $CCarr1;
#	print Dumper $CCarr2;
#	print "Xmean: $mean_x Ymean: $mean_y XX: $ssxx YY: $ssyy XY: $ssxy Cor: $xcorrel Cov: $cov\n";
	return($StatsKit_success, $xcorrel, $cov, $CCVelements);
}

sub CallCorrel {
	my($ssxx,$ssyy,$ssxy)=@_;
#	print "$ssxx,$ssyy,$ssxy\n";
	my $sign=$ssxy/abs($ssxy);
	my $correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));
	return $correl;
}



###Get min/max value from a arr
###&GetValues (min/max, @num_arr)
###Global: 
###Dependency:
###Nnote:
sub GetValues {
	my ($GVindex, @GVnumbers)=@_;
	@GVnumbers=sort {$a<=>$b} @GVnumbers;
	if ($GVindex=~/^min$/i) {
		return defined $GVnumbers[0] ? $GVnumbers[0]:'err';
	}
	elsif ($GVindex=~/^max$/i) {
		return defined $GVnumbers[-1] ? $GVnumbers[-1]:'err';
	}
}

1;
