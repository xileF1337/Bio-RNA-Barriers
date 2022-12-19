#!perl -T
use 5.012;
use strict;
use warnings FATAL => 'all';

use Test::More;
use Test::NoWarnings;            # produces one additional test!
use Test::Exception;
use File::Spec::Functions qw(catfile);
use File::Slurp qw(read_file);

my $test_count = 7;
plan tests => $test_count + 1;            # +1 for Test::NoWarnings

use Bio::RNA::Barriers;

# Data files
my $rate_matrix_file_txt = catfile qw(t data 12.rates.txt);

my @file_12_txt_rates = qw(
   0.06188    0.01468   0.002897
    0.2387          0          0
    0.2387          0          0
);


##############################################################################
##                              Test functions                              ##
##############################################################################

# runs 3 tests
sub test_rates_from_to {
    my ($rate_matrix, $rates_ref) = @_;

    can_ok $rate_matrix, qw(rate_from_to);

    my $dim = $rate_matrix->dim;
    my @rates = @$rates_ref;                # deep copy

    subtest 'rates from i to j' => sub {
        plan tests => 3**2;
        foreach my $i (1..$dim) {
            foreach my $j (1..$dim) {
                cmp_ok $rate_matrix->rate_from_to($i, $j),
                       '==',
                       shift @rates,
                       "rate from $i to $j"
                       ;
            }
        }
    };

    dies_ok { $rate_matrix->rate_from_to(1, $dim+1) },
            'Rate for out-of-bound states dies';
}

# runs 3 tests
sub test_set_rate_from_to {
    my ($rate_matrix, $rates_ref) = @_;

    can_ok $rate_matrix, qw(set_rate_from_to);

    my ($from, $to, $rate) = (1, 2, 0.12345);
    $rate_matrix->set_rate_from_to(1, 2, $rate);

    cmp_ok $rate_matrix->rate_from_to($from, $to), '==', $rate,
           "Setting rate from $from to $to to $rate";

    dies_ok { $rate_matrix->set_rate_from_to($from, $to, -1*$rate) },
            'Setting negative rate dies';
}


##############################################################################
##                                Call tests                                ##
##############################################################################

-s $rate_matrix_file_txt
    or BAIL_OUT "empty or non-existent data file '$rate_matrix_file_txt'";

# Open input data file.
my $rate_matrix = Bio::RNA::Barriers::RateMatrix->new(
        file_name => $rate_matrix_file_txt,
        file_type => 'TXT',
);

# Check we've got the right amount of rates.
my $dim = $rate_matrix->dim;
cmp_ok $dim*$dim, '==', scalar(@file_12_txt_rates), 'number of rates in matrix';

test_rates_from_to $rate_matrix, \@file_12_txt_rates;       # 3 tests

test_set_rate_from_to $rate_matrix, \@file_12_txt_rates;    # 3 tests

exit 0;

