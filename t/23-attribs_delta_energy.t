#!perl -T
use 5.012;
use strict;
use warnings;

use Test::More;
use Test::NoWarnings;            # produces one additional test!
use File::Spec::Functions;

plan tests => 2 + 1;

use Bio::RNA::Barriers;

# Test whether delta_E of given file matches the passed energy.
sub test_delta_energy {
    my ($barfile, $delta_energy) = @_;

    open my $barfh, '<', $barfile
        or BAIL_OUT "failed to open test data file '$barfile'";

    my $bar_results = Bio::RNA::Barriers::Results->new($barfh);

    cmp_ok $bar_results->delta_energy, '==', $delta_energy,
           "delta energy ($barfile)";
}


# Run tests.
test_delta_energy catfile(qw(t data 12.bar)),           10.0;
test_delta_energy catfile(qw(t data disconnected.bar)),  8  ;


# EOF
