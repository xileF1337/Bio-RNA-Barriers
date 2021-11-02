package Bio::RNA::Barriers;

use 5.012;
use strict;
use warnings;

our $VERSION = '0.01';

package Bio::RNA::Barriers::Results {
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
}

package Bio::RNA::Barriers::Minimum {
    use Moose;
    use MooseX::StrictConstructor;
    use Moose::Util::TypeConstraints;
    use namespace::autoclean;

    use autodie qw(:all);
    use overload q{""} => 'stringify';

    use Scalar::Util    qw(blessed);
    use List::Util      qw(max);
    use List::MoreUtils qw(zip);

    #### Special types for attribute checking.
    subtype 'RNAStruct' => (
        as 'Str',
        where { m{ ^ [(.)]+ $ }x },
        message {
            "Only '(', ')', and '.' allowed in structure string, found '$_'"
        },
    );

    subtype 'DisconSaddle' => (
        as 'Str',
        where { m{ ^ ~+ $ }x },
        message {
            "Only '~' allowed in disconnected saddle string, found '$_'"
        },
    );

    # index  - index of basins ordered by energy; 1 is lowest
    # struct - struct of lowest energy in minimums
    # mfe    - free energy of the basin's local minimum
    # father_index  - index of father basin (the basin this one is merged to)
    # barrier_height - height of energy barrier (in kcal/mol) to minimum this
    #                 one is merged to (relative to this minimum)
    my @default_attribs = qw( index struct mfe father_index barrier_height);
    my @default_attrib_args = (is => 'rw', required => 1);
    my %default_attrib_isa
        = &zip(\@default_attribs, [qw(Int RNAStruct Num Int Num)]);
    has $_ => (@default_attrib_args, isa => $default_attrib_isa{$_})
        foreach @default_attribs;

    # Return true iff this is the mfe basin 1.
    sub is_global_min {
        my $self = shift;
        my $is_global_min = $self->index == 1;
        return $is_global_min;
    }

    # Optional attributes generated by Barriers options --bsize and --saddle.
    # Descriptions in quotes are from Barriers tutorial at
    # https://www.tbi.univie.ac.at/RNA/tutorial/#sec4_2
    # merged_struct_count - 'numbers of structures in the basin we merge with'
    #       Given is the number of structures in the *current* basin
    #       (including merged ones) *at the time of merging*. For minimum 1,
    #       this is close to the total number of input structures (except for
    #       disconnected structures and other missing ones (???).
    # father_struct_count - 'number of basin which we merge to'
    #       Actually, it's the number of *structures* in the basin that we
    #       merge to (father basin) *at the time of merging*.
    # merged_basin_energy - 'free energy of the basin'
    #       This seems to be the free energy of (the partition function of)
    #       the basin---including all merged basins---at the time this basin
    #       is merged. For minimum 1, this corresponds to the ensemble free
    #       energy as far as it was enumerated by RNAsubopt (excluding
    #       disconnected structures).
    # grad_struct_count - 'number of structures in this basin using gradient walk'
    #       This seems to be the actual number of structures only in this
    #       basin, excluding merged basins. What about --minh merging? Why
    #       doesn't this column sum up to exactly to the total number of
    #       structs if --max==Inf (some are missing)? Issues due to degenerate
    #       energies?
    # grad_basin_energy - 'gradient basin (consisting of all structures where
    #                      gradientwalk ends in the minimum)'
    #       This seems to be free energy of the basin without any merged
    #       basins. Summing up the partition functions corresponding to these
    #       energies, one obtains a free energy almost equal to the ensemble
    #       energy (up to rounding errors due to 6 digit precision).
    my @bsize_attributes = qw(
        merged_struct_count father_struct_count merged_basin_energy
        grad_struct_count grad_basin_energy
    );
    my @opt_attributes = (@bsize_attributes, qw(saddle_struct));
    # Define a bsize predicate only for the first bsize attribute. Ensure in
    # BUILD that either all or none of the attributes are set.
    my @common_attrib_args = (
        is => 'ro',
        lazy => 1,
        default => sub { confess 'attribute undefined, did you use --bsize/--saddle?' },
    );
    # has $bsize_attributes[0] => (is => 'ro', predicate => 'has_bsize' );
    # has $_ => (is => 'ro') foreach @bsize_attributes[1..$#bsize_attributes];
    # has 'saddle_struct' => (is => 'ro', predicate => 'has_saddle_struct');
    has $bsize_attributes[0] => (
        @common_attrib_args,
        isa       => 'Num',
        predicate => 'has_bsize',
        writer    => "_$bsize_attributes[0]",           # private writer
    );
    has $_ => (@common_attrib_args, writer => "_$_")    # private writer
        foreach @bsize_attributes[1..$#bsize_attributes];
    has 'saddle_struct' => (
        @common_attrib_args,
        isa       => 'RNAStruct | DisconSaddle',
        predicate => 'has_saddle_struct',
    );

    # Optional reference to the father minimum.
    has 'father' => (
        is        => 'rw',
        trigger   => \&_check_father,
        # Use name 'has_father_REF' because has_father==false reads as if the
        # min does not have a father at all.
        # We don't need this, we always have it if we have a father.
        # predicate => 'has_father_ref',
    );

    sub _check_father {
        my ($self, $father) = @_;
        confess 'Need a reference to another minimum to set father attribute'
            unless blessed $father and $father->isa( __PACKAGE__ );
        confess "Father's index does not match the index used during construction"
            unless $self->father_index == $father->index;
    }

    # Returns true iff the minimum has a father minimum it has been merged to.
    sub has_father {
        my $self = shift;
        my $has_father = $self->father_index > 0;
        return $has_father;
    }


    # Minimum is connected to basin 1 (mfe).
    has 'is_connected' => (
        is => 'ro',
        lazy => 1,
        init_arg => undef,          # cannot be set manually
        builder => '_build_is_connected',
    );

    sub _build_is_connected {
        my $self = shift;

        return 1 if $self->index == 1;              # this is the mfe basin
        return 0 unless $self->has_father;          # basin has no father

        # confess 'Reference to father minimum has not been set, cannot proceed.'
        #     unless $self->has_father_ref;

        my $is_connected = $self->father->is_connected;
        return $is_connected;
    }

    # Parse passed line read from barriers file.
    around BUILDARGS => sub {
        my $orig  = shift;
        my $class = shift;

        my @args;                               # as passed to the constructor
        if ( @_ == 1 && !ref $_[0] ) {          # process line from bar file
            my $input_line = shift;
            my @fields = split /\s+/, $input_line;
            shift @fields if $fields[0] eq q{}; # drop empty first field

            if (@fields < @default_attribs) {
                confess "Input line has not enough fields: $input_line";
            }

            # Add default args
            push @args, $_ => shift @fields foreach @default_attribs;
            # @args = map { $_ => shift @fields } @attributes;

            # Add saddle struct if present
            if (@fields == 1 or @fields == @opt_attributes) {
                push @args, saddle_struct => shift @fields;
            }

            # Add bsize attributes if present
            if (@fields == @bsize_attributes) {
                push @args, $_ => shift @fields foreach @bsize_attributes;
            }

            confess "Unrecognized number of fields on input line:\n$input_line"
                unless @fields == 0;            # all fields used up?
        }
        else {
            @args = @_;
        }
        return $class->$orig(@args);
    };

    sub BUILD {
        my $self = shift;

        # Ensure presence or absence of all bsize attributes
        my $defined_count = grep {defined $self->{$_}} @bsize_attributes;
        confess "Need to define all or none of the --bsize attributes ",
              join q{, }, @bsize_attributes
            unless $defined_count == 0 or $defined_count == @bsize_attributes;
    }

    # Determine all ancestor minima of this minimum, i.e. this minimum's
    # father, grand father, grand grand father etc. in this order.
    # Returns list of all ancestors (may be empty if min is disconnected).
    sub ancestors {
        my ($self) = @_;

        my $ancestor = $self;
        my @ancestors;
        while ($ancestor->has_father) {
            # confess 'Need father reference to determine ancestors'
            #     unless $ancestor->has_father_ref;
            push @ancestors, $ancestor->father;
            $ancestor = $ancestor->father;
        }

        return @ancestors;
    }

    # Stringify minimum to equal an entry line in the barriers output file.
    # Format strings are taken from Barriers' C source code.
    sub stringify {
        my $self = shift;

        # Default attributes
        my $min_string = $self->brief;

        # Add saddle struct if defined.
        if ($self->has_saddle_struct) {
            $min_string .= q{ } . $self->saddle_struct;
        }

        # Add bsize attributes if defined.
        if ($self->has_bsize) {
            $min_string .= sprintf " %12ld %8ld %10.6f %8ld %10.6f",
                                   map {$self->{$_}} @bsize_attributes;
        }

        return $min_string;
    }

    # Stringification method returning a more brief representation of the
    # minimum containing only the index, min struct, its energy, the father's
    # index, and the barrier height. This equals the output of Barriers if
    # neither --bsize nor --saddle is given.
    sub brief {
        my $self = shift;

        # Default attributes
        my $brief_string = sprintf "%4d %s %6.2f %4d %6.2f",
                                    map {$self->{$_}} @default_attribs;

        return $brief_string;
    }

    # ABSOLUTE energy of the lowest structure connecting this basin to another
    # one. The barrier height, in contrast, is the RELATIVE energy of the same
    # structure w.r.t. to the basin's (local) mfe.
    # BEWARE: if the basin does not have a father (father == 0), then the
    # barrier height is (as reported by Barriers) given with respect to the global exploration
    # threshold.
    # Since this gives unexpected results, the saddle height is set to the
    # global mfe for basin 1, and to Inf for disconnected basins.
    sub saddle_height {
        my $self = shift;

        # Mfe basin is connected to itself with a barrier of 0. Other
        # fatherless basins are disconnected and thus have an unknown saddle
        # height -- set to Inf.
        my $barrier_height =   $self->has_father    ? $self->barrier_height
                             : $self->is_global_min ? 0
                                                    : 'Inf'
                             ;

        # Energy values from Bar file have only 2 digits precision.
        my $saddle_height = sprintf "%.2f", $self->mfe + $barrier_height;
        return $saddle_height;
    }

    # Saddle height as described for saddle_height(), but with respect to the
    # global mfe structure (basin 1).
    sub global_saddle_height {
        my $self = shift;

        # Move up in barrier tree until reachin basin 1 or realizing we are
        # disconnected. Global saddle height is maximal encountered height.
        my $ancestor         = $self;
        my $glob_sadd_height = $self->saddle_height;
        while ($ancestor->father_index > 1) {
            $ancestor = $ancestor->father;
            $glob_sadd_height
                = max $glob_sadd_height, $ancestor->saddle_height;
        }

        return $glob_sadd_height;
    }

    __PACKAGE__->meta->make_immutable;
}   # Bio::RNA::Barriers::Minimum

package Bio::RNA::Barriers::RateMatrix {
    use Moose;
    use MooseX::StrictConstructor;
    use namespace::autoclean;
    use Moose::Util::TypeConstraints qw(enum subtype as where message);

    use autodie qw(:all);
    use overload '""' => \&stringify;
    use Scalar::Util qw( reftype looks_like_number );
    use List::Util qw( all uniqnum );

    enum __PACKAGE__ . 'RateMatrixType', [qw(TXT BIN)];

    # Natural number type directly from the Moose docs.
    subtype 'PosInt',
        as 'Int',
        where { $_ > 0 },
        message { "The number you provided, $_, was not a positive number" };

    has 'file_name' => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_file_name',
    );
    has 'file_type' => (
        is       => 'rw',
        isa      => __PACKAGE__ . 'RateMatrixType',
        required => 1,
    );
    has '_file_handle' => (
        is       => 'ro',
        isa      => 'FileHandle',
        init_arg => 'file_handle',
        lazy     => 1,
        builder  => '_build_file_handle',
    );
    # Splice rate matrix directly when reading the data. This can read big
    # matrices when only keeping a few entries.
    has 'splice_on_parsing' => (
        is        => 'ro',
        isa       => 'ArrayRef[PosInt]',
        predicate => 'was_spliced_on_parsing',
    );
    has '_data' => (is => 'ro', lazy => 1, builder => '_build_data');

    sub BUILD {
        my $self = shift;

        # Enforce data is read from handle immediately despite laziness.
        $self->dim;
    }

    # Read the actual rate data from the input file and construct the
    # matrix from it.
    sub _build_data {
        my $self = shift;

        my $rate_matrix;
        if ($self->file_type eq 'TXT') {
            $rate_matrix = __PACKAGE__->read_text_rate_matrix(
                $self->_file_handle,
                $self->splice_on_parsing,
            );
        }
        elsif ($self->file_type eq 'BIN') {
            $rate_matrix = __PACKAGE__->read_bin_rate_matrix(
                $self->_file_handle,
                $self->splice_on_parsing,
            );
        }
        else {
            confess "Unknown file type, that's a bug...";
        }

        return $rate_matrix;
    }

    # Class method. Reads a rate matrix in text format from the passed file
    # handle and constructs a matrix (2-dim array) from it. Returns a
    # reference to the constructed rate matrix.
    # Arguments:
    #   input_matrix_fh: file handle to text file containing rate matrix
    #   splice_to: ORDERED set of states which are to be kept. The other
    #       states are pruned from the matrix on-the-fly while parsing.
    #       This saves time and memory.
    sub read_text_rate_matrix {
        my ($class, $input_matrix_fh, $splice_to_ref) = @_;

        # During parsing, splice the selected rows / columns. Make 0-based.
        my @splice_to_rows = @{ $splice_to_ref // [] };     # 1-based, modified
        my @splice_to_cols = map {$_ - 1} @splice_to_rows;  # 0-based indices

        my (@rate_matrix, $matrix_dim);
        ROW: while (defined (my $line = <$input_matrix_fh>)) {
            if (defined $splice_to_ref) {
                last unless @splice_to_rows;            # we're done!
                next ROW if $. != $splice_to_rows[0];   # this row is not kept
                shift @splice_to_rows;              # remove the leading index
            }

            my @row = split q{ }, $line;                # awk-style splitting
            # Since the diagonal element may be more or less anything, we need
            # to check it separately (e.g. to not choke on BHGbuilder output).
            my @row_no_diag = @row[0..($.-2), ($.)..$#row];    # $. is 1-based
            confess 'Input file contains non-numeric or negative input on ',
                    "line $.:\n$line"
                unless looks_like_number $row[$.-1]     # diag elem can be <0
                       and all {looks_like_number $_ and $_ >= 0} @row_no_diag;

            # Check that element count is equal in all rows.
            $matrix_dim //= @row;         # first-time init
            confess 'Lines of input file have varying number of elements'
                unless $matrix_dim == @row;

            @row = @row[@splice_to_cols] if defined $splice_to_ref;
            push @rate_matrix, \@row;
            confess 'Input file contains more lines than there are columns'
                if @rate_matrix > $matrix_dim;
        }
        confess 'End of file reached before finding all states requested by ',
                'splicing operation'
            if defined $splice_to_ref and @splice_to_rows > 0;
        confess 'Requested splicing of non-contained state'
            unless all {$_ < $matrix_dim} @splice_to_cols;
        confess 'Input file is empty'
            unless @rate_matrix;
        # Adjust dimension if splicing was applied.
        confess 'Input file contains less lines than there are columns'
            if @rate_matrix < (defined $splice_to_ref ? @splice_to_cols
                                                      : $matrix_dim     );

        return \@rate_matrix;
    }

    sub _transpose_matrix {
        my ($class, $matrix_ref) = @_;

        # Determine dimnensions
        my $max_row = @$matrix_ref - 1;
        return unless $max_row >= 0;
        my $max_col = @{ $matrix_ref->[0] } - 1;    # check elems of first row

        # Swap values
        for my $row (0..$max_row) {
            for my $col (($row+1)..$max_col) {
                my $temp = $matrix_ref->[$row][$col];
                $matrix_ref->[$row][$col] = $matrix_ref->[$col][$row];
                $matrix_ref->[$col][$row] = $temp;
            }
        }
    }

    # Class method. Reads a rate matrix in binary format from the passed file
    # handle and constructs a matrix (2-dim array) from it. Returns a
    # reference to the constructed rate matrix.
    sub read_bin_rate_matrix {
        my ($class, $input_matrix_fh, $splice_to_ref) = @_;

        # During parsing, splice the selected rows / columns. Make 0-based.
        my @splice_to_cols = @{ $splice_to_ref // [] };     # 1-based, modified
        my @splice_to_rows = map {$_ - 1} @splice_to_cols;  # 0-based indices

        # Set read mode to binary
        binmode $input_matrix_fh;

        ##### Read out matrix dimension
        my $size_of_int = do {use Config; $Config{intsize}};
        my $read_count
            = read($input_matrix_fh, my $raw_matrix_dim, $size_of_int);
        confess "Could not read dimension from file, ",
              "expected $size_of_int bytes, got $read_count"
            if $read_count != $size_of_int;

        my $matrix_dim = unpack 'i', $raw_matrix_dim;       # unpack integer

        confess 'Requested splicing of non-contained state'
            unless all {$_ < $matrix_dim} @splice_to_rows;

        ##### Read rate matrix
        my @rate_matrix;
        my $size_of_double = do {use Config; $Config{doublesize}};
        my $bytes_per_column = $size_of_double * $matrix_dim;
        COL: for my $i (1..$matrix_dim) {
            # Each column consists of n=matrix_dim doubles.
            $read_count
                = read($input_matrix_fh, my $raw_column, $bytes_per_column);
            confess "Could not read column $i of file, ",
                  "expected $bytes_per_column bytes, got $read_count"
                if $read_count != $bytes_per_column;

            # Skip column if splicing and column not requested.
            if (defined $splice_to_ref) {
                last unless @splice_to_cols;            # we're done!
                next COL if $i != $splice_to_cols[0];   # this col is not kept
                shift @splice_to_cols;              # remove the leading index
            }

            # Decode raw doubles.
            my @matrix_column = unpack "d$matrix_dim", $raw_column;

            # Splice parsed column if requested.
            @matrix_column = @matrix_column[@splice_to_rows]
                if defined $splice_to_ref;

            push @rate_matrix, \@matrix_column;
        }
        confess 'End of file reached before finding all states requested by ',
                'splicing operation'
            if defined $splice_to_ref and @splice_to_cols > 0;
        confess 'Read data as suggested by dimension, but end of file ',
                'not reached'
            unless defined $splice_to_ref or eof $input_matrix_fh;

        # For whatever reasons, binary rates are stored column-wise instead of
        # row-wise. Transpose to fix that.
        __PACKAGE__->_transpose_matrix(\@rate_matrix);

        return \@rate_matrix;
    }

    sub _build_file_handle {
        my $self = shift;

        confess 'File required if no file handle is passed'
            unless $self->has_file_name;

        open my $handle, '<', $self->file_name;
        return $handle;
    }

    # Get the dimension (= number of rows = number of columns) of the matrix.
    sub dim {
        my $self = shift;

        my $dimension = @{ $self->_data };
        return $dimension;
    }

    # Get the rate from state i to state j. States are 1-based (first state =
    # state 1) just as in the results file.
    sub rate_from_to {
        my ($self, $from_state, $to_state) = @_;

        # Check states are within bounds
        confess "from_state $from_state is out of bounds"
            unless $self->_state_is_in_bounds($from_state);
        confess "to_state $to_state is out of bounds"
            unless $self->_state_is_in_bounds($to_state);

        # Retrieve rate.
        my $rate = $self->_data->[$from_state-1][$to_state-1];
        return $rate;
    }

    # Check whether given state is contained in the rate matrix.
    sub _state_is_in_bounds {
        my ($self, $state) = @_;

        my $is_in_bounds = ($state >= 1 && $state <= $self->dim);
        return $is_in_bounds;
    }

    # Returns a sorted list of all states connected to the (mfe) state 1.
    # Assumes a symmetric transition matrix (only checks path *from* state 1
    # *to* the other states). Quadratic runtime.
    sub connected_states {
        my ($self) = @_;

        # Starting at state 1, perform a traversal of the transition graph and
        # remember all nodes seen.
        my $dim = $self->dim;
        my @cue = (1);
        my %connected = (1 => 1);       # state 1 is connected
        while (my $i = shift @cue) {
            foreach my $j (1..$dim) {
                next if $connected{$j} or $self->rate_from_to($i, $j) <= 0;
                $connected{$j} = 1;         # j is connected to 1 via i
                push @cue, $j;
            }
        }

        # Sort in linear time.
        my @sorted_connected = grep {$connected{$_}} 1..$dim;
        return @sorted_connected;
    }

    # Only keep the states connected to the mfe (as determined by
    # connected_states()).
    sub keep_connected {
        my ($self) = @_;
        my @connected_indices = map {$_ - 1} $self->connected_states;
        return map {$_ + 1} @connected_indices      # none removed.
            if $self->dim == @connected_indices;

        $self->_splice_indices(\@connected_indices);

        return map {$_ + 1} @connected_indices;       # turn into states again
    }

    # Remove all but the passed states from this rate matrix. States are
    # 1-based (first state = state 1) just as in the results file.
    sub keep_states {
        my ($self, @states_to_keep) = @_;

        # We need a sorted, unique list.
        @states_to_keep = uniqnum sort {$a <=> $b} @states_to_keep;

        # Check whether states are within bounds.
        foreach my $state (@states_to_keep) {
            confess "State $state is out of bounds"
                unless $self->_state_is_in_bounds($state);
        }

        return if @states_to_keep == $self->dim;    # keep all == no op

        $_-- foreach @states_to_keep;               # states are now 0-based
        $self->_splice_indices(\@states_to_keep);

        return $self;
    }

    # Only keep the passed states and reorder them as in the passed list. In
    # particular, the same state can be passed multiple times and will then be
    # deep-copied.
    # Arguments:
    #   states: Ordered list of states defining the resulting matrix. May
    #       contain duplicates.
    sub splice {
        my ($self, @states) = @_;

        # Check whether states are within bounds.
        foreach my $state (@states) {
            confess "State $state is out of bounds"
                unless $self->_state_is_in_bounds($state);
        }

        $_-- foreach @states;                       # states are now 0-based
        $self->_splice_indices(\@states);

        return $self;
    }

    # Internal version which performs no boundary checks and assumes REFERENCE
    # to state list.
    sub _splice_indices {
        my ($self, $kept_indices_ref) = @_;

        my $matrix_ref = $self->_data;

        # If no entries are kept, make matrix empty.
        if (@$kept_indices_ref == 0) {
            @$matrix_ref = ();
            return;
        }

        # Splice the matrix.
        # WARNING: This makes a shallow copy of the rows if the same index is
        # passed more than once (e.g. from splice()).
        @$matrix_ref = @{$matrix_ref}[@$kept_indices_ref];  # rows

        # Deep-copy duplicated rows (if any).
        my %row_seen;
        foreach my $row (@$matrix_ref) {
            $row = [@$row] if $row_seen{$row};              # deep-copy array
            $row_seen{$row} = 1;
        }
        @$_ = @{$_}[@$kept_indices_ref]                     # columns
            foreach @$matrix_ref;

        return $self;
    }

    # Remove the passed states from this rate matrix. States are 1-based
    # (first state = state 1) just as in the results file.
    sub remove_states {
        my ($self, @states_to_remove) = @_;

        return unless @states_to_remove;        # removing no states at all

        # Check states are within bounds.
        foreach my $state (@states_to_remove) {
            confess "State $state is out of bounds"
                unless $self->_state_is_in_bounds($state);
        }

        # Invert state list via look-up hash.
        my %states_to_remove = map {$_ => 1} @states_to_remove;
        my @states_to_keep
            = grep {not $states_to_remove{$_}} 1..$self->dim;

        # Let _keep_indices() do the work.
        $_-- foreach @states_to_keep;               # states are now 0-based
        $self->_splice_indices(\@states_to_keep);

        return $self;
    }

    # Print this matrix as text, either to the passed handle, or to STDOUT.
    sub print_as_text {
        my ($self, $text_matrix_out_fh) = @_;
        $text_matrix_out_fh //= \*STDOUT;       # write to STDOUT by default

        my $rate_format = '%10.4g ';            # as in Barriers code

        foreach my $row (@{ $self->_data }) {
            printf {$text_matrix_out_fh} $rate_format, $_ foreach @$row;
            print  {$text_matrix_out_fh} "\n";
        }
    }

    # Print this matrix as binary data, either to the passed handle or to
    # STDOUT.  Data format: matrix dimension as integer, then column by column
    # as double.
    sub print_as_bin {
        my ($self, $rate_matrix_out_fh ) = @_;

        my $rate_matrix_ref = $self->_data;

        # Set write mode to binary
        binmode $rate_matrix_out_fh;

        ##### Print out matrix dimension
        my $matrix_dim = @$rate_matrix_ref;
        my $packed_dim = pack 'i', $matrix_dim;     # machine representation, int
        print {$rate_matrix_out_fh} $packed_dim;

        ##### Print columns of rate matrix
        # For whatever reasons, binary rates are stored column-wise instead of
        # row-wise (Treekin works with the transposed matrix and this way it's
        # easier to slurp the entire file. Treekin transposes the text rates
        # during reading).
        #_transpose_matrix $rate_matrix_ref;
        foreach my $col (0..($matrix_dim-1)) {
            foreach my $row (0..($matrix_dim-1)) {
                # Pack rate as double
                my $packed_rate = pack 'd', $rate_matrix_ref->[$row][$col];
                print {$rate_matrix_out_fh} $packed_rate;
            }
            # my $column = map {$_->[$i]} @$rate_matrix_ref;
            # my $packed_column = pack "d$matrix_dim", @column;
        }
    }

    # Return string containing binary the representation of the matrix (cf.
    # print_as_bin).
    sub serialize {
        my $self = shift;

        # Use print function and capture matrix in a string.
        my $matrix_string;
        open my $matrix_string_fh, '>', \$matrix_string;
        $self->print_as_bin($matrix_string_fh);

        return $matrix_string;
    }

    # Returns a string containing the text representation of the matrix. The
    # overloaded double-quote operator calls this method.
    sub stringify {
        my $self = shift;

        # Use print function and capture matrix in a string. Empty matrices
        # give an empty string (not undef).
        my $matrix_string = q{};
        open my $matrix_string_fh, '>', \$matrix_string;
        $self->print_as_text($matrix_string_fh);

        return $matrix_string;
    }

}           # package Bio::RNA::Barriers::RateMatrix


1; # End of Bio::RNA::Barriers
__END__

=head1 NAME

Bio::RNA::Barriers - Classes for working with I<Barriers> output.

=head1 VERSION

Version 0.01

=cut

=head1 SYNOPSIS

This module provides auxiliary classes to parse, query and print the
files generated by the RNA landscape coarse graining tool I<Barriers>.

Supports the results file (written to STDOUT by I<Barriers> as well as rate
matrices in binary and text format.

    use Bio::RNA::Barriers;

    ##### Working with the results file (bar file) #####
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

    # Parse a single minimum line
    my $min_string = ...;
    Bio::RNA::Barriers::Minimum->new($min_string);


    ##### Working with the rate matrix files (rates.{out,bin}) #####
    # Read a binary rate matrix directly from file. Binary matrices are more
    # precise and smaller than text matrices.
    my $rate_matrix = Bio::RNA::Barriers::RateMatrix->new(
        file_name => '/path/to/rates.bin',
        file_type => 'BIN',
    );

    # Read a text rate matrix from an opened handle.
    open my $rate_matrix_fh_txt, '<', '/path/to/rates.out';
    my $rate_matrix = Bio::RNA::Barriers::RateMatrix->new(
        file_handle => $rate_matrix_fh_txt,
        file_type   => 'TXT',
    );

    # Print out matrix, dimension, and a single rate.
    print "$rate_matrix";
    print 'Dimension of rate matrix is ', $rate_matrix->dim, "\n";
    print 'Rate from state 1 to state 3 is ',
          $rate_matrix->rate_from_to(1, 3),
          "\n";

    # Remove entries for a list of states {1, 3, 5} (1-based as in bar file).
    $rate_matrix->remove_states(1, 5, 5, 3);    # de-dupes automatically
    # Note: former state 2 is now state 1 etc.

    # Keep only states {1, 2, 3}, remove all others. Can also de-dupe.
    $rate_matrix->keep_states(1..3);

    # Write binary matrix to file.
    open my $out_fh_bin, '>', '/path/to/output/rates.bin';
    $rate_matrix->print_as_bin($out_fh_bin);

=head1 CLASSES AND METHODS

This module provides two major classes: C<Bio::RNA::Barriers::Results> and
C<Bio::RNA::Barriers::Minimum>, where the first is an aggregate results object
containing objects from the second. Usually you want to pass a file name or
handle to a I<Barriers> file to the constructor of the results class, and the
rest is taken care of automatically.

There are some other miscellaneous classes included, too. TODO

=head2 C<Bio::RNA::Barriers::Results>

This is what you usually want to use. Pass a file name or a handle to the
constructor and you're golden.

=head3 C<seq()>

Returns the RNA sequence for which the minima have been computed.

=head3 C<file_name()>

Name of the file from which the minima have been read. May be C<undef>.

=head3 C<get_min($index)>

Return the single minimum given by a 1-based C<$index>.

=head3 C<get_mins(@indices)>

Return a list of minima given by a 1-based list of C<@indices>.

=head3 C<get_global_min()>

Return the basin represented by the (global) mfe structure (i.e. basin 1).

=head3 C<min_count()>

Return the total number of basins.

=head3 C<connected_mins()>

List of all minima connected to the mfe minimum (min 1).

=head3 C<connected_indices()>

List of indices of all connected minima (cf. connected_mins()).

=head3 C<keep_connected()>

Keep only connected minima and remove all others. The indices are remapped to
1..k for k kept minima.  Returns (old) indices of all kept minima.

=head3 C<update_connected()>

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

=head3 C<keep_first_mins()>

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

=head3 C<has_bsize()>

Check whether the stored minima contain information about the basin sizes as
computed by Barriers' --bsize switch. Checks only the mfe basin (since all or
neither min should have this information).

=head3 C<path()>

Construct the file path to this Barriers file. Works only if it was actually
parsed from a file (of course...).  Returns the path as a string.

=head3 C<stringify()>

Convert back to Barriers file. Also supports stringification overloading.


=head2 C<Bio::RNA::Barriers::Minimum>

The minima are usually generated and managed by a
C<Bio::RNA::Barriers::Results> object. They do have nice methods, though.

Returns a list of all ancestors (may be empty if min is disconnected).

=head3 C<has_bsize()>

Boolean. True iff the minimum provides information about the basin size as
computed by the I<Barriers> option C<--bsize>. TODO which attributes


=head3 C<is_global_min()>

Returns true iff this is the global minimum (i.e. basin 1).

=head3 C<index()>

1-based index of the minimum (as is the Barriers file).

=head3 C<struct()>

Returns the dot-bracket structure string of the minimum.

=head3 C<mfe()>

B<Local> minimum free energy of the basin (i.e. the minimums energy).

=head3 C<father_index()>

Returns the index of the father minimum (i.e. the one this minimum has been
merged to).

=head3 C<barrier_height()>

Returns the barrier height (B<relative> energy difference of the saddle point
to the local minimum). For the B<absolute> energy of the saddle point, see
C<saddle_height()>.

=head3 C<saddle_height()>

B<Absolute> energy of the lowest structure connecting this basin to another
one. The barrier height, in contrast, is the B<relative> energy of the same
structure w.r.t. to the basin's (local) mfe.

B<Beware>: if the basin does not have a father (father == 0), then the
reported saddle height is given with respect to the global exploration
threshold. This is strange but consistent with original Barriers files.

=head3 C<global_saddle_height()>

Saddle height as described for saddle_height(), but not with respect to
any neighbor minimum, but to the global mfe structure (basin 1).


=head3 C<saddle_struct()>

Returns the saddle structure via which it was merged to its father minimum.
If this attribute was not set (i.e. I<Barriers> was run without the
C<--saddle> option), it croaks when accessed.

Use the C<has_saddle()> predicate to query the status of this attribute.

=head3 C<has_saddle_struct()>

Predicate for the C<saddle> attribute. True iff the minimum provides the
saddle structure via which it was merged to its father minimum, as computed by
the I<Barriers> option C<--saddle>. It can be queried via the C<saddle_struct>
method.

=head3 C<father()>

Returns (a reference to) the father minimum.

=head3 C<has_father()>

Returns true iff the minimum has a father minimum it has been merged to.

=head3 C<is_connected()>

Boolean. True iff the minimum is connected to basin 1 (the mfe basin).

=head3 C<ancestors()>

Determine all ancestor minima of this minimum, i.e. this minimum's father,
grand father, grand grand father etc. in this order.

=head3 C<stringify()>

Stringify minimum to equal an entry line in the barriers output file.
Format strings are taken from Barriers' C source code.

=head3 C<brief()>

Stringification method returning a more brief representation of the
minimum containing only the index, min struct, its energy, the father's
index, and the barrier height. This equals the output of Barriers if
neither C<--bsize> nor C<--saddle> is given.

=head2 C<Bio::RNA::Barriers::RateMatrix>

Parse, modify and print/write rate matrix files, both in text and binary
format.

=head3 C<new()>

Constructor. Reads a rate matrix from a file / handle and creates a new rate
matrix object.

=over 4

=item Arguments:

=over 4

=item file_name | file_handle

Source of the data to read. Pass either or both.

=item file_type

Specifies whether the input data is in binary or text format. Must be either
C<'TXT'> or C<'BIN'>.

=back

=back

=head3 file_name

File from which the data was read. May be undef.

=head3 file_type

Specifies whether the input data is in binary or text format. Must be either
C<'TXT'> or C<'BIN'>.

=head3 C<read_text_rate_matrix($input_matrix_filehandle)>

Class method. Reads a rate matrix in text format from the passed file
handle and constructs a matrix (2-dim array) from it. Returns a
reference to the constructed rate matrix.

=head3 C<read_bin_rate_matrix($input_matrix_filehandle)>

Class method. Reads a rate matrix in binary format from the passed file
handle and constructs a matrix (2-dim array) from it. Returns a
reference to the constructed rate matrix.

=head3 C<dim()>

Get the dimension (= number of rows = number of columns) of the matrix.

=head3 C<rate_from_to($i, $j)>

Get the rate from state i to state j. States are 1-based (first state = state
1) just as in the results file.

=head3 C<remove_states(@indices)>

Remove the passed states from this rate matrix. States are 1-based (first
state = state 1) just as in the results file.

=head3 C<connected_states()>

Returns a sorted list of all states connected to the (mfe) state 1.
Assumes a symmetric transition matrix (only checks path B<from> state 1
B<to> the other states). Quadratic runtime.

=head3 C<keep_connected()>

Only keep the states connected to the mfe (as determined by
C<connected_states()>).

=head3 C<keep_states(@indices)>

Remove all but the passed states from this rate matrix. States are 1-based
(first state = state 1) just as in the results file.

=over 4

=item Arguments:

=over 4

=item indices

List of states to be kept. It can be unordered and may contain duplicates.

=back

=back

=head3 C<splice()>

Only keep the passed states and reorder them as in the passed list. In
particular, the same state can be passed multiple times and will then be
deep-copied.

=over 4

=item Arguments:

=over 4

=item indices

Ordered list of states defining the resulting matrix. May contain duplicates.

=back

=back

=head3 C<print_as_text()>

Print this matrix as text, either to the passed handle, or to STDOUT.

=head3 C<print_as_bin()>

Print this matrix as binary data, either to the passed handle or to
STDOUT.  Data format: matrix dimension as integer, then column by column
as double.

=head3 C<serialize()>

Return string containing binary representation of the matrix (cf.
print_as_bin).

=head3 C<stringify()>

Returns a string containing the text representation of the matrix. The
overloaded double-quote operator calls this method.


=head1 AUTHOR

Felix Kuehnl, C<< <felix at bioinf.uni-leipzig.de> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-rna-barriers at
rt.cpan.org>, or through the web interface at
L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-RNA-Barriers>.  I will
be notified, and then you'll automatically be notified of progress on your bug
as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::RNA::Barriers


You can also look for information at the official Barriers website:

L<https://www.tbi.univie.ac.at/RNA/Barriers/>


=over 4

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

Copyright 2019 Felix Kuehnl.

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


