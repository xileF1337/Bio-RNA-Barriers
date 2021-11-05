package Bio::RNA::Barriers::Results;
our $VERSION = '0.01';

use 5.012;
use strict;
use warnings;

use Moose;
use MooseX::StrictConstructor;
use namespace::autoclean;

use autodie qw(:all);
use Scalar::Util qw( reftype );
use File::Spec::Functions qw(catpath splitpath);

has [qw( seq file_name )] => (
    is       => 'ro',
    required => 1,
);

has _mins => (
    is       => 'ro',
    required => 1,
    init_arg => 'mins'
);

has $_ => (
    is => 'ro',
    predicate => "has_$_",
) foreach qw(_volume _directory);


use overload q{""} => 'stringify';

# Allow various calling styles of the constructor:
# new(barriers_handle): pass only file handle to read data from, file_name
#       will be undef.
# new(barriers_handle, barriers_file_path): read data from handle and set
#       file name from a given path
# new( barriers_file_path): open handle and read data from the passed
#       path, set file name accordingly
# new(hash_ref): manual construction from required attributes
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;

    return $class->$orig(@_) if @_ == 0 or @_ > 2;  # no special handling

    # Read Barriers file from handle and set constructor arguments.
    my $barriers_fh;
    my @args;
    if (@_ == 1 and not reftype $_[0]) {            # file path given
        my $file_path = shift;
        my ($volume, $directory, $file_name) = splitpath $file_path;
        push @args, (
            _volume    => $volume,
            _directory => $directory,
            file_name  => $file_name,
        );

        open $barriers_fh, '<', $file_path;
    }
    elsif (reftype $_[0] eq reftype \*STDIN) {      # file handle given
        if (@_ == 2) {
            # Check if second arg may be a file name
            return $class->$orig(@_) if reftype $_[1];  # it's no file name
            push @args, (file_name => $_[1]);
        }
        else {
            push @args, (file_name => undef);       # no file name given
        }
        $barriers_fh = shift;
    }
    else {
        return $class->$orig(@_);
    }

    # Parse data from handle and create arguments.
    push @args, $class->_read_barriers_file($barriers_fh);
    return $class->$orig(@args);
};

sub _read_barriers_file {
    my ($self, $barriers_fh) = @_;
    my %args;

    # Read sequence from first line
    my $sequence = do {
        my $line = <$barriers_fh>;
        confess 'Could not read sequence from Barriers file handle'
            unless defined $line;
        chomp $line;
        $line =~ s/^\s*//r;             # trim leading space
    };
    $args{seq} = $sequence;

    # Continue reading minima
    my ($father, @mins);
    while (defined (my $line = <$barriers_fh>)) {
        my $min = Bio::RNA::Barriers::Minimum->new($line);
        push @mins, $min;
    }

    confess 'No minima found in Barriers file'
        unless @mins;

    # Set references to father mins (if any). This needs to be done after
    # minima construction because, if several minima share the same
    # energy, a father may have a higher index than the doughter min (and
    # thus not have been parsed when parsing the doughter).
    foreach my $min (@mins) {
        if ($min->has_father) {
            my $father = $mins[$min->father_index - 1];
            confess 'Inconsistent father indexing in Barriers file'
                unless ref $father
                       and $father->index == $min->father_index;
            $min->father($father);
        }
    }

    $args{mins} = \@mins;

    return %args;
}

# Dereferenced list of all minima. Select individual minima using
# get_min[s]().
sub mins {
    my $self = shift;
    return @{ $self->_mins };
}

sub get_mins {
    my ($self, @min_indices) = @_;

    foreach my $min_index (@min_indices) {
        confess "There is no minimum $min_index in this file"
            if $min_index < 1 or $min_index > $self->mins;
        $min_index--;           # Perl is 0-based, barriers 1-based
    }

    # my @mins = ($self->mins)[@min_indices];     # array slice
    my @mins = @{$self->_mins}[@min_indices];     # array slice
    return @mins;
}

sub get_min {
    my ($self, $min_index) = @_;
    my ($min) = $self->get_mins($min_index);
    return $min;
}

sub get_global_min {
    my ($self) = @_;
    my $mfe_min = $self->get_min(1);            # basin 1 == mfe basin
    return $mfe_min;
}

sub min_count {
    my $self = shift;
    my $min_count = $self->mins;
    return $min_count;
}

# List of all minima connected to the mfe minimum (min 1).
sub connected_mins {
    my $self = shift;
    my @connected_mins = grep {$_->is_connected} $self->mins;
    return @connected_mins;
}

# List of indices of all connected minima (cf. connected_mins()).
sub connected_indices {
    my $self = shift;
    my @connected_indices = map {$_->index} $self->connected_mins;
    return @connected_indices;
}

# Re-index minima after deleting some of them.
sub _reindex {
    my $self = shift;

    my $i = 1;
    my @mins = $self->mins;

    # Update min indices.
    $_->index($i++) foreach @mins;

    # Update father indices.
    shift @mins;                        # min 1 is orphan
    $_->father_index($_->father->index) foreach @mins;

    return;
}

# Keep only connected minima and remove all others. The indices are
# remapped to 1..k for k kept minima.
# Returns (old) indices of all kept minima.
sub keep_connected {
    my $self = shift;

    my @connected_indices = $self->connected_indices;
    my @connected_mins    = $self->connected_mins;

    # Update minima list
    @{ $self->_mins } = @connected_mins;
    $self->_reindex;

    return @connected_indices;
}

# Given an ordered list of (all???) indices of connected minima, delete
# all other minima and update their ancesters' basin size information
# accordingly.
# Arguments:
#   ordered_connected_indices: ordered list of indices of (all???)
#       connected minima.
sub update_connected {
    my ($self, @ordered_connected_indices) = @_;

    # Go through all mins and check whether they're next in the connected
    # (==kept) index list. If not, add to removal list.
    my @connected_mins = $self->get_mins(@ordered_connected_indices);
    my @removed_indices;
    for my $min_index (1..$self->min_count) {
        unless (@ordered_connected_indices) {
            # No exclusions left, add rest and stop
            push @removed_indices, $min_index..$self->min_count;
            last;
        }
        if ($min_index == $ordered_connected_indices[0]) {
            shift @ordered_connected_indices;
            next;
        }
        push @removed_indices, $min_index;      # min is deleted
    }

    my @removed_mins  = $self->get_mins(@removed_indices);
    $self->_update_ancestors(@removed_mins);
    @{ $self->_mins } = @connected_mins;
    $self->_reindex;

    return;
}

# Pass a list of ORDERED deleted minima and update their ancestors' bsize
# attributes.
sub _update_ancestors {
    my ($self, @removed_mins) = @_;

    # If the bsize attributes are present, update the basin energy of all
    # (grand)* father basins (substract energy of this one).
    # The minima need to be processed in reversed order because, if an
    # ancester of a removed min is also removed, its merged basin energy
    # includes the basin energy of its child, and thus this contribution
    # would be substracted multiple times from older ancesters. In
    # reversed order, the child contribution is first substracted from the
    # ancestors, and then the contribution of the removed ancestors does
    # not include the child anymore.

    return unless $self->has_bsize;

    foreach my $removed_min (reverse @removed_mins) {
        foreach my $ancestor_min ($removed_min->ancestors) {
            $ancestor_min->_merged_basin_energy(        # private writer
                $ancestor_min->merged_basin_energy
                - $removed_min->grad_basin_energy
            );
        }
    }

    return;
}

# Keep only the first k mins. If there are only k or less mins, do
# nothing.
# WARNING: THIS CAN DISCONNECT THE LANDSCAPE! The bar file will still look
# connected, however, modifying the rate matrix accordingly can lead to
# non-ergodicity (e.g. when basin 3 merged to 2, 2/3 merged to 1 because
# of a possible transition from 3 to 1, and basin 3 is then removed).
# To cope with that, call RateMatrix::keep_connected() and, with its
# return value, Results::update_connected() again.
# Arguments:
#   max_min_count: number of minima to keep.
# Returns a list of all removed mins (may be empty).
sub keep_first_mins {
    my ($self, $max_min_count) = @_;

    my @removed_mins = splice @{ $self->_mins }, $max_min_count;
    $self->_update_ancestors(@removed_mins);

    # # If the bsize attributes are present, update the basin energy of all
    # # (grand)* father basins (substract energy of this one).
    # # The minima need to be processed in reversed order because, if an
    # # ancester of a removed min is also removed, its merged basin energy
    # # includes the basin energy of its child, and thus this contribution
    # # would be substracted multiple times from older ancesters. In
    # # reversed order, the child contribution is first substracted from the
    # # ancestors, and then the contribution of the removed ancestors does
    # # not include the child anymore.
    # foreach my $removed_min (reverse @removed_mins) {
    #     foreach my $ancestor_min ($removed_min->ancestors) {
    #         $ancestor_min->_merged_basin_energy(        # private writer
    #             $ancestor_min->merged_basin_energy
    #             - $removed_min->grad_basin_energy
    #         );
    #     }
    # }

    return @removed_mins;
}


# Check whether the stored minima contain information about the basin
# sizes as computed by Barriers' --bsize switch. Checks only the mfe
# basin (since all or neither min should have this information).
sub has_bsize {
    my ($self) = @_;
    my $has_bsize = $self->get_global_min->has_bsize;
    return $has_bsize;
}


# Construct the file path to this Barriers file. Works only if it was
# actually parsed from a file (of course...).
# Returns the path as a string.
sub path {
    my $self = shift;
    confess 'Cannot construct path, missing volume, directory or file name'
        unless defined $self->file_name
               and $self->has_volume
               and $self->has_directory
               ;
    my ($volume, $dir, $file_name)
        = ($self->_volume, $self->_directory, $self->file_name);
    my $path = catpath($volume, $dir, $file_name);
    return $path;
}


# Convert back to Barriers file.
sub stringify {
    my $self = shift;

    my $header = (q{ } x 5) . $self->seq;
    my @lines  = ($header, map { $_->stringify } $self->mins);

    return join "\n", @lines;
}

__PACKAGE__->meta->make_immutable;

1;                                  # End of Bio::RNA::Barriers::Results


__END__

=pod

=encoding UTF-8

=head1 NAME

Bio::RNA::Barriers::Results - Parse, query and manipulate results
of a I<Barriers> run

=head1 SYNOPSIS

    use Bio::RNA::Barriers;

    # Read in a Barriers output file.
    open my $barriers_handle, '<', $barriers_file;
    my $bardat = Bio::RNA::Barriers::Results->new($barriers_handle);
    # Or, even simpler, pass file name directly:
    $bardat    = Bio::RNA::Barriers::Results->new($barriers_file  );

    # Print some info
    print "There are ", $bardat->min_count, " minima.";
    my $min3 = $bardat->get_min(3);
    print $min3->grad_struct_count,
          " structures lead to basin 3 via a gradient walk.\n";
        if $min3->has_bsize;
    print "Min ", $min3->index, " is ", ($min3->is_connected ? "" : "NOT"),
          " connected to the mfe structure.\n";

    # Print the mfe basin line as in the results file
    my $mfe_min = $bardat->get_min(1);
    print "$mfe_min\n";


=head1 DESCRIPTION

This is what you usually want to use. Pass a file name or a handle to the
constructor and you're golden. When querying a specific minimum by index, a
L<Bio::RNA::Barriers::Minimum> object is returned. For more details on its
methods, please refer to its documentation.

=head1 METHODS

=head3 $res->new()

TODO constructor.

=head3 $res->seq()

Returns the RNA sequence for which the minima have been computed.

=head3 $res->file_name()

Name of the file from which the minima have been read. May be C<undef>.

=head3 $res->mins()

Returns a list of all minima. Useful for iteration in C<for> loops.

=head3 $res->get_min($index)

Return the single minimum given by a 1-based C<$index>.

=head3 $res->get_mins(@indices)

Return a list of minima given by a 1-based list of C<@indices>.

=head3 $res->get_global_min()

Return the basin represented by the (global) mfe structure (i.e. basin 1).

=head3 $res->min_count()

Return the total number of basins.

=head3 $res->connected_mins()

List of all minima connected to the mfe minimum (min 1).

=head3 $res->connected_indices()

List of indices of all connected minima (cf. connected_mins()).

=head3 $res->keep_connected()

Keep only connected minima and remove all others. The indices are remapped to
1..k for k kept minima.  Returns (old) indices of all kept minima.

=head3 $res->update_connected()

Given an ordered list of (all???) indices of connected minima, delete
all other minima and update their ancesters' basin size information
accordingly.

=over 4

=item Arguments:

=over 4

=item ordered_connected_indices:

Ordered list of indices of (all???) connected minima to keep.

=back

=back

=head3 $res->keep_first_mins()

Keep only the first k mins. If there are only k or less mins, do
nothing.

WARNING: THIS CAN DISCONNECT THE LANDSCAPE! The bar file will still look
connected, however, modifying the rate matrix accordingly can lead to
non-ergodicity (e.g. when basin 3 merged to 2, 2/3 merged to 1 because
of a possible transition from 3 to 1, and basin 3 is then removed).

=over 4

=item Arguments:

=over 4

=item max_min_count:

number of minima to keep.

=back

=back

Returns a list of all removed mins (may be empty).

=head3 $res->has_bsize()

Check whether the stored minima contain information about the basin sizes as
computed by Barriers' --bsize switch. Checks only the mfe basin (since all or
neither min should have this information).

=head3 $res->path()

Construct the file path to this Barriers file. Works only if it was actually
parsed from a file (of course...).  Returns the path as a string.

=head3 $res->stringify()

Convert back to Barriers file. Also supports stringification overloading.


=cut


# End of lib/Bio/RNA/Barriers/Results.pm