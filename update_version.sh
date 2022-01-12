#!/bin/bash
# Update $VERSION variable in all *.pm files in lib/.

set -eEuo pipefail

# Abort on errors, displaying error message + code, kill all running jobs.
clean_and_die() {
    error_code="$1"; error_message="$2"
    echo -e "\nERROR: $error_message ($error_code) in script" \
        "'$(basename $0)'" 1>&2
    jobs=$(jobs -pr); [ -z "$jobs" ] || kill $(jobs -pr)
    exit $error_code
}

trap 'clean_and_die $? "terminated unexpectedly at line $LINENO"' ERR
trap 'clean_and_die  1 "interrupted"'           INT
trap 'clean_and_die  1 "caught TERM signal"'    TERM

# Display passed error message on STDERR. warn() continues while die() exits
# with error code 1.
die()  { echo   "ERROR: $*" 1>&2; exit 1; }
warn() { echo "WARNING: $*" 1>&2;         }


##############################################################################
##                              Script options                              ##
##############################################################################

# The line to replace. Must include everything up to quoted version number.
# Spaces are replaced by a \s* pattern.
version_line="our \$VERSION = "


# Search at beginning of line, and replace spaces by \s*.
version_pattern="^$(sed 's/ /\\s*/g' <<< "$version_line")"


##############################################################################
##                                Functions                                 ##
##############################################################################

get_current_version() {
    local module="$1"
    grep "$version_pattern" "$module" |
            head --lines=1 |
            sed "s/$version_pattern//" |
            tr --delete "\"';" ||
        die "No version info found in file ${modules[0]}"
}

##############################################################################
##                                   Main                                   ##
##############################################################################

# Search modules in lib tree.
modules=( $(find lib -iname '*.pm') )
[ ${#modules[@]} -gt 0 ] || die 'No Perl modules (*.pm) found in lib/'

# Get current version from first module.
current_version="$(get_current_version "${modules[0]}")"
echo "Current module version: $current_version"

# Check positional argument for new version number.
[ $# -eq 1 ] || die 'Pass new version number. Aborting.'
new_version="$1"

# Prompt user to allow update.
read -ep  "Update to version $new_version? (y/n) " answer
[ "$answer" == 'y' -o "$answer" == 'Y' ] || die 'Aborting.'

# Update each module.
subst_regex="s/$version_pattern[\"'].*[\"']\s*;/$version_line'$new_version';/"
errors=
for module in "${modules[@]}"; do
    echo -n "    updating $module .. "
    sed --in-place "$subst_regex" "$module"     # replace old version

    # Verify that it worked.
    subst_version="$(get_current_version "$module")"
    [ "$subst_version" == "$new_version" ] &&
        echo 'ok' ||
        { errors=1; echo 'NOT ok'; }
done

echo 'Done.'
[ -z "$errors" ] || warn 'Could not update some module(s).'
exit 0;                                # EOF
