#!perl -T
use 5.012;
use strict;
use warnings;

use Test::More;
use Test::NoWarnings;            # produces one additional test!
use File::Spec::Functions qw(catfile);

my $test_count = 22;
plan tests => $test_count + 1;            # +1 for Test::NoWarnings

use Bio::RNA::Barriers;

# Input files
my $barfile_starred = catfile qw(t data starred.bar);
my %min_has_star = (            # the ground truth
    (map {$_ => 0} 1..18, 21),  # no stars
    (map {$_ => 1} 19..20, 22), # have stars
);

# Read input, construct results.
open my $barfile_fh, '<', $barfile_starred
    or BAIL_OUT "failed to open test data file '$barfile_starred'";

my $bar_res = Bio::RNA::Barriers::Results->new($barfile_fh);

# Return true iff two values are logically equal, i.e. they both evaluate to
# either true or false.
sub logic_eq {
    my ($a, $b) = @_;

    return ($a && $b) || (!$a && !$b);
}

##### Run tests #####

# Check that all mins have a star or not, as intended.
for my $min ($bar_res->mins) {
    my $descr = 'min ' . $min->index . ' has '
                . ($min_has_star{$min->index} ? 'a' : 'no') . ' star';
    ok logic_eq($min->has_star, $min_has_star{$min->index}), $descr;
}


# EOF
