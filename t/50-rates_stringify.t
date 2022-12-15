#!perl -T
use 5.012;
use strict;
use warnings;

use Test::More;
use Test::NoWarnings;            # produces one additional test!
use Test::Exception;
use IO::File;                    # we need autoflush
use File::Spec::Functions qw(catfile);
use File::Slurp qw(read_file);
use File::Temp qw(tempfile);

my $test_count = 10;
plan tests => $test_count + 1;            # +1 for Test::NoWarnings

use Bio::RNA::Barriers;

my ($rate_matrix_data_txt, $rate_matrix_data_bin);
eval {
    $rate_matrix_data_txt = read_file catfile qw(t data 12.rates.txt);
    $rate_matrix_data_bin
        = read_file catfile(qw(t data 12.rates.bin)), binmode => ':raw';
};
BAIL_OUT "failed to open test data file: $@" if $@;

# Construct rate matrix from exact binary data.
open my $rate_matrix_fh_bin, '<', \$rate_matrix_data_bin;
my $rate_matrix = Bio::RNA::Barriers::RateMatrix->new(
    file_handle => $rate_matrix_fh_bin,
    file_type   => 'BIN',
);

##### Functions #####

# Truncate the file and reset file handle position.
sub reset_fh {
    my ($fh) = @_;

    seek $fh, 0, 0;
    truncate $fh, 0;
}

##### Run tests #####

# Stringify & serialize
can_ok $rate_matrix, qw(stringify serialize);
is $rate_matrix->stringify, $rate_matrix_data_txt, 'stringifies correctly';
is "$rate_matrix",          $rate_matrix_data_txt, 'stringify overloading';
is $rate_matrix->serialize, $rate_matrix_data_bin, 'serializes correctly';

# Write to file: text.
can_ok $rate_matrix, qw(print_as_bin print_as_text print_as);

# Create temp file we can write to. Auto-flush to be able to read it again.
my ($temp_fh, $temp_file) = tempfile('50-rates_stringify.t.XXXXX', UNLINK => 1);
$temp_fh->autoflush;

$rate_matrix->print_as_text($temp_fh);
is scalar(read_file($temp_file)), $rate_matrix_data_txt, "print_as_text()";

reset_fh($temp_fh);
$rate_matrix->print_as('TXT', $temp_fh);
is scalar(read_file($temp_file)), $rate_matrix_data_txt, "print_as('TXT')";

# Write to file: bin. Do this after text mode as the print_to_bin() methods
# sets the mode of the given file handle to binary.
reset_fh($temp_fh);
$rate_matrix->print_as_bin($temp_fh);
is scalar(read_file($temp_file, binmode => ':raw')),
   $rate_matrix_data_bin,
   "print_as_bin()";

reset_fh($temp_fh);
$rate_matrix->print_as('BIN', $temp_fh);
is scalar(read_file($temp_file, binmode => ':raw')),
   $rate_matrix_data_bin,
   "print_as('BIN')";

# Catch wrong mode.
dies_ok { $rate_matrix->print_as('FOOBAR', $temp_fh) };


exit 0;

