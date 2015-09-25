# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::FileKit

=head1 SYNOPSIS

File operations

=head1 Requirements

Cwd

=head1 DESCRIPTION

CloseHandlerIfOpen( filehandlers )

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
package FuhaoPerl5Lib::FileKit;
use strict;
use warnings;
use Exporter;
use File::Copy;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = '20150603';
@ISA         = qw(Exporter);
@EXPORT      = qw(RetrieveDir RetrieveName AddFilePath RetrieveBasename DeletePath MoveFile CopyFile MergeFiles);
@EXPORT_OK   = qw();
%EXPORT_TAGS = ( DEFAULT => [qw(RetrieveDir RetrieveName AddFilePath RetrieveBasename DeletePath &MoveFile CopyFile MergeFiles)],
                 ALL    => [qw(RetrieveDir RetrieveName AddFilePath RetrieveBasename DeletePath &MoveFile CopyFile MergeFiles)]);


my $FileKit_success=1;
my $FileKit_failure=0;


### NOT working
### close file handler if exists 
### CloseHandlerIfOpen (Handlers)
### Global:
### Dependency:
### Note:
sub CloseHandlerIfOpen {
	my @CHIOhanlers=@_;
	
	my $CHIOsubinfo='SUB(FileKit::LuCloseHandlerIfOpen)';
	
	foreach my $CHIOindv_handler (@CHIOhanlers) {
		if (defined fileno($CHIOindv_handler)) {
			unless (close $CHIOindv_handler) {
				print STDERR "${CHIOsubinfo}Error: can not close filehandle for $CHIOindv_handler\n";
				return 0;
			}
		}
	}
	return 1;
}



### Add Path before file
### &AddFilePath ( file1, file2 , ... )
### Global: 
### Dependency: Cwd
### Note: 
sub AddFilePath {
	my @AFPfiles=@_;
	use Cwd 'realpath';
	foreach (@AFPfiles) {
		$_=realpath($_);
	}
	return (@AFPfiles);
}



### Retrieve filename
### &RetrieveName ( file )
### Global: 
### Dependency: &AddFilePath
### Note: 
sub RetrieveName {
	my $RB_ori=shift @_;
	chomp $RB_ori;
	$RB_ori=&AddFilePath($RB_ori);
	my $RB_new='';
	($RB_new=$RB_ori)=~ s/.*\///s;
	return $RB_new;
}



### Retrieve filebasename without extension
### &RetrvNoExt( file )
### Global:
### Dependency:
### Note:
sub RetrieveBasename {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my $RNE_new='';
	my $RNE_base='';
	($RNE_base=$RNE_ori)=~ s/.*\///s;
	($RNE_new=$RNE_base)=~s/^(\S+)\.\w+$/$1/;
	return $RNE_new;
}



### Retrive path for a file
### &RetrvDir(file)
### Global:
### Depedency: &AddFilePath
### Note:
sub RetrieveDir {
	my $RD_ori=shift @_;
	chomp $RD_ori;
	($RD_ori)=&AddFilePath($RD_ori);
	(my $RD_new=$RD_ori)=~ s/(.*)\/.*$/$1/s;
	return ($RD_new);
}



###delete oath and its contents
#&DeletePath(PATH);
#chdir $curDir || die "Error when deleting directory: $curDir\n";
sub DeletePath {
	my $DPpath = shift @_;
	
	my $DPsubinfo="'SUB(FileKit::DeletePath)";
	
	#get all the files in that directory.
	my @DPfiles=glob "$DPpath/*";
	foreach my $DPind_file (@DPfiles){
		if(-d $DPind_file){
		#if the destination file is a directory, go recursion.
			unless (&DeletePath($DPind_file)) {
				print $DPsubinfo, "Error: can not delete PATH: $DPind_file\n";
				return $FileKit_failure;
			}
		}
		else{
			unless (unlink $DPind_file) {
				print $DPsubinfo, "Error: can not delete File: $DPind_file\n";
				return $FileKit_failure;
			}
		}
	}
	#del the destination directory.
	unless (rmdir $DPpath) {
		print $DPsubinfo, "Error: can not delete PATH: $DPpath\n";
		return $FileKit_failure;
	}
	return $FileKit_success;
}



sub BackupFile {
	my $BUori_name=shift @_;
	my $BU_new_name='';
	$BU_new_name=$BUori_name.".bak.".int(rand(10000));
	if (-e "$BU_new_name") {
		&BackupFile($BUori_name);
	}
	rename("$BUori_name" , "$BU_new_name");
	print "The Original file $BUori_name was renamed to $BU_new_name\n";
	return $FileKit_success;
}



### Move file using File::Copy
sub MoveFile {
	my ($MFori, $MFtar)=@_;
	
	my $MFsubinfo='SUB(FileKit::MoveFile)';
	
	unless (defined $MFori and $MFori=~/^\S+$/ and -s $MFori) {
		print STDERR "${MFsubinfo}Error: invalid original name to move\n";
		return $FileKit_failure;
	}
	unless (defined $MFtar and $MFtar=~/^\S+$/) {
		print STDERR "${MFsubinfo}Error: invalid target name to move\n";
		return $FileKit_failure;
	}
	unlink ($MFtar) if (-e $MFtar);
	unless (move($MFori, $MFtar)) {
		return $FileKit_failure;
	}
	return $FileKit_success;
}


### Copy file using File::Copy
sub CopyFile {
	my ($CFori, $CFtar)=@_;
	
	my $CFsubinfo='SUB(FileKit::CopyFile)';
	
	unless (defined $CFori and $CFori=~/^\S+$/ and -s $CFori) {
		print STDERR "${CFsubinfo}Error: invalid original name to move\n";
		return $FileKit_failure;
	}
	unless (defined $CFtar and $CFtar=~/^\S+$/) {
		print STDERR "${CFsubinfo}Error: invalid target name to copy\n";
		return $FileKit_failure;
	}
	unlink ($CFtar) if (-e $CFtar);
	unless (copy($CFori, $CFtar)) {
		return $FileKit_failure;
	}
	return $FileKit_success;
}



###
### MergeFiles
### Global: 
### Dependency: $FileKit_success;$FileKit_failure;

### Note:
sub MergeFiles {
	my $MFout=shift;
	my @MFfiles2merge=@_;
	
	local *MFIN; local *MFOUT;
	my $MFsubinfo='SUB(FileKit::MergeFiles)';
	unlink $MFout if (-e $MFout);
	
	close MFOUT if (defined fileno(MFOUT));
	unless (open(MFOUT, "> $MFout")) {
		print STDERR "Error: unable to write $MFout\n";
		return $FileKit_failure;
	}
	foreach my $MFindfile (@MFfiles2merge) {
		print $MFsubinfo, "Test: Mergeing file: $MFindfile\n"; ### For test ###
		unless (-s $MFindfile) {
			print STDERR $MFsubinfo, "Warnings: not existed file to be merged: $MFindfile\n" unless (-e $MFindfile); ### For test ###
			next;
		}
		close MFIN if (defined fileno(MFIN));
		unless (open (MFIN, "< $MFindfile")) {
			print STDERR $MFsubinfo, "Error: can not open $MFindfile\n";
			return $FileKit_failure;
		}
		while (my $MFline=<MFIN>) {
			print MFOUT $MFline;
		}
		close MFIN;
	}
	print $MFsubinfo, "Test: merging file finished: $MFout\n";
	close MFOUT;
	if (-s $MFout) {
		print $MFsubinfo, "Test: $MFout\n";
		return $FileKit_success;
	}
	else {
		return $FileKit_failure;
	}
}

1;
