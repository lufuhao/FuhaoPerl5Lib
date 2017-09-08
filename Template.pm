# POD documentation - main docs before the code

=head1 NAME

FuhaoPerl5Lib::GffKit

=head1 SYNOPSIS

use samtools to index, extract, merge, rehead BAM files

=head1 Requirements

Perl Modules:


=head1 DESCRIPTION

=over 2

=item Bam2FastQ ($bamin, $fastqout, map_code, MAPQ_code, [path_samtools])

    * Convert bam files into fastq
    * Map_code: 0=all    1=mapped_only    2=unmapped_only
    * MAPQ_code: any alignment less than MAPQ would be ignored
    * Return: 1=Sucesss    0=Failure

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

package FuhaoPerl5Lib::GffKit;
use strict;
use warnings;
use Cwd;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '20170322';
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 Both    => [qw()]);



my $GffKit_success=1;
my $GffKit_failure=0;
my $GffKit_debug=0;






1;
