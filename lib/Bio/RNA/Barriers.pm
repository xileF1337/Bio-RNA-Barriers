package Bio::RNA::Barriers;
our $VERSION = '0.01';

use 5.012;
use strict;
use warnings;

use Bio::RNA::Barriers::Minimum;
use Bio::RNA::Barriers::RateMatrix;
use Bio::RNA::Barriers::Results;

1; # End of Bio::RNA::Barriers


__END__

=pod

=encoding UTF-8

=head1 NAME

Bio::RNA::Barriers - Parse, query and manipulate output of I<Barriers>

=cut

=head1 SYNOPSIS

    use Bio::RNA::Barriers;

    ##### Working with the Barriers results file (*.bar) #####
    $bardat = Bio::RNA::Barriers::Results->new('bar_file.bar');

    print "There are ", $bardat->min_count, " minima.";

    my $min3 = $bardat->get_min(3);
    print $min3->grad_struct_count,
          " structures lead to basin 3 via a gradient walk.\n"
        if $min3->has_bsize;

    my $mfe_min = $bardat->get_global_min();
    print "$mfe_min\n";             # prints minimum as in the results file


    ##### Working with the rate matrix files (rates.{out,bin}) #####
    my $rate_matrix = Bio::RNA::Barriers::RateMatrix->new(
        file_name => '/path/to/rates.bin',
        file_type => 'BIN',
    );

    print "$rate_matrix";               # prints entire matrix in text format
    print 'Dimension of rate matrix is ', $rate_matrix->dim, "\n";
    print 'Rate from state 1 to state 3 is ',
          $rate_matrix->rate_from_to(1, 3),
          "\n";

    $rate_matrix->remove_states(1, 5, 5, 3);    # state 2 becomes state 1 etc.
    $rate_matrix->keep_states(1..3);            # keep only states {1, 2, 3}.
    $rate_matrix->keep_connected();             # remove disconnected states

    open my $out_fh_bin, '>', '/path/to/output/rates.bin';
    $rate_matrix->print_as_bin($out_fh_bin);    # write binary output


=head1 DESCRIPTION

This module provides auxiliary classes to parse, query, manipulate and print
the files generated by I<Barriers>, a tool to compute RNA energy landscapes
and folding kinetics developed at Theoretical Biochemistry Group (TBI) at the
University of Vienna. Note that this module is B<not> written and maintained
by the authors of I<Barriers>.

Supports the result file (written to STDOUT by I<Barriers>) as well as rate
matrices in binary and text format. Properties like the number of minima or,
for each minima, their fathers and children, connectedness, and basin sizes
can be queried.

Rate matrices can be manipulated easily, e. g. to remove or keep certain
states or to convert between binary and text representation.


=head1 CLASSES

This module provides two major classes to handle the results:
L<Bio::RNA::Barriers::Results> and L<Bio::RNA::Barriers::Minimum>, where the
first is an aggregate results object containing objects from the second.
Usually you want to pass a file name or handle to a I<Barriers> file to the
constructor of the results class (i. e.,
C<< Bio::RNA::Barriers::Results->new() >>), and the rest is taken care of
automatically.

For rate matrices, L<Bio::RNA::Barriers::RateMatrix> is provided. It can parse
both text and binary matrices. Make sure to correctly set the file type
argument when passing the path or handle to the constructor.

For a description of the available methods, please refer to the documentation
of each individual class, i. e. run C<perldoc Bio::RNA::Barriers::Results>
etc.


=head1 AUTHOR

Felix Kuehnl, C<< <felix at bioinf.uni-leipzig.de> >>

=head1 BUGS

Please report any bugs or feature requests by raising an issue at
L<https://github.com/xileF1337/Bio-RNA-Barriers/issues>.

You can also do so by mailing to C<bug-bio-rna-barmap at rt.cpan.org>,
or through the web interface at
L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-RNA-BarMap>.  I will be
notified, and then you'll automatically be notified of progress on your bug as
I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::RNA::Barriers


You can also look for information at the official Barriers website:

L<https://www.tbi.univie.ac.at/RNA/Barriers/>


=over 4

=item * Github: the official repository

L<https://github.com/xileF1337/Bio-RNA-Barriers>

=item * RT: CPAN's request tracker (report bugs here)

L<https://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-RNA-Barriers>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-RNA-Barriers>

=item * CPAN Ratings

L<https://cpanratings.perl.org/d/Bio-RNA-Barriers>

=item * Search CPAN

L<https://metacpan.org/release/Bio-RNA-Barriers>

=back


=head1 LICENSE AND COPYRIGHT

Copyright 2019-2021 Felix Kuehnl.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see L<http://www.gnu.org/licenses/>.


=cut

# End of Bio/RNA/Barriers.pm
