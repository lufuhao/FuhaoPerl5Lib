# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::CmdKit

=head1 SYNOPSIS

Execuate cmd and return

=head1 DESCRIPTION

&CurrentTime
Require: localtime (Linux)
Output: YYYYMMDD HHMMSS

&exec_cmd ( cmd )
Require: time
Output: 1=success, die if fail

&exec_cmd_return ( cmd )
Require: time
Output: 1=success, die if fail

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
package FuhaoPerl5Lib::CmdKit;
use strict;
use warnings;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150603';
@ISA         = qw(Exporter);
@EXPORT      = qw(CurrentTime exec_cmd exec_cmd_return);
@EXPORT_OK   = qw();
%EXPORT_TAGS = ( DEFAULT => [qw(&CurrentTime &exec_cmd &exec_cmd_return)],
                 ALL    => [qw(&CurrentTime &exec_cmd &exec_cmd_return &CurrentTime)]);

my $CmdKit_success=1;
my $CmdKit_failure=0;
my $CmdKit_debug=0;


### print current date and time
### &CurrentTime()
### Global:
### Dependency:
### Note: localtime
sub CurrentTime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;
}



### Process command
### &exec_cmd(cmd)
### Global:
### Dependency: time (Linux), &CurrentTime
### Note: 
sub exec_cmd {
	my ($cmd) = @_;
	print "#####\n".&CurrentTime()."CMD: $cmd\n";
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
#	if ($return_code == -1) {
#		print “failed to execute: $!\n”;
#	}
#	elsif ($return_code & 127) {
#		printf “child died with signal %d, %s coredump\n”, ($return_code & 127),  ($return_code & 128) ? ‘with’ : ‘without’;
#	}
#	else {
#		printf “child exited with value %d\n”, $return_code >> 8;
#	}
	if ($return_code) {
#		print "Error, cmd: $cmd died with ReturnCode $return_code\n";
		die "SUB(exec_cmd)Error, cmd: $cmd died with ReturnCode $return_code\n";
	}
	else {
		print STDERR "Finished command: $cmd\tat ".&CurrentTime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $CmdKit_success;
	}
}



### Process command
### &exec_cmd(cmd)
### Global:
### Dependency: time (Linux), &CurrentTime
### Note:
sub exec_cmd_return {
	my ($cmd) = @_;
	print "##### CMD Starts\n".&CurrentTime()."CMD: $cmd\n";# if ($CmdKit_debug);
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
	if ($return_code) {
		print STDERR "SUB(exec_cmd_return)Error, cmd: $cmd died with ReturnCode $return_code\n##### CMD Ends\n";# if ($CmdKit_debug);
		return $CmdKit_failure;
	}
	else {
		print STDERR "Finished command: $cmd\tat ".&CurrentTime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n##### CMDEnds\n";# if ($CmdKit_debug);
		return $CmdKit_success;
	}
}





### Test a external cmd exists in PATH
### Global: 
sub CanRun {
	my $CRcmd=shift;
	
	my $CRsubinfo='SUB(MiscKit::CanRun)';
	my $CRpath=which ("$CRcmd");
	
	if (-s $CRpath and -x $CRpath) {
		return $$CmdKit_success;
	}
	else {
		print STDERR "Error: cmd $CRcmd not exists\n";
		return $CmdKit_failure;
	}
}



1;
