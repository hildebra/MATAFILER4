#!/usr/bin/env bash
# Install and update the environments used by MATAFILER4.

set -Eeuo pipefail
ulimit -c 0 || true

REMOVE_LEGACY_ENVS=0
REFRESH_DATABASES=0
INSTALL_MARKER=""

usage() {
	cat <<'EOF'
Usage: installer.sh [options]

Options:
  --remove-legacy-envs  Remove obsolete MGTK-named micromamba environments.
  --refresh-databases   Redownload the CheckM2 and MetaPhlAn databases.
  -h, --help            Show this help message.
EOF
}

die() {
	echo "ERROR: $*" >&2
	exit 1
}

retry_command() {
	local description=$1
	local max_attempts=$2
	local delay_seconds=$3
	shift 3

	local attempt status
	for ((attempt = 1; attempt <= max_attempts; attempt++)); do
		if "$@"; then
			return 0
		else
			status=$?
		fi

		if ((attempt == max_attempts)); then
			echo "$description failed after $max_attempts attempts." >&2
			return "$status"
		fi

		echo "$description failed (attempt $attempt/$max_attempts); retrying in $delay_seconds seconds." >&2
		sleep "$delay_seconds"
	done
}

while (($#)); do
	case "$1" in
		--remove-legacy-envs) REMOVE_LEGACY_ENVS=1 ;;
		--refresh-databases) REFRESH_DATABASES=1 ;;
		-h|--help) usage; exit 0 ;;
		*) usage >&2; die "Unknown option: $1" ;;
	esac
	shift
done

if [[ -n "${MAMBA_EXE:-}" ]]; then
	if [[ -x "$MAMBA_EXE" ]]; then
		MAMBA_E="$MAMBA_EXE"
	else
		MAMBA_E="$(command -v -- "$MAMBA_EXE" || true)"
	fi
else
	MAMBA_E="$(command -v micromamba || true)"
fi

[[ -n "$MAMBA_E" ]] || die "micromamba could not be found. Set MAMBA_EXE or add micromamba to PATH."
command -v git >/dev/null 2>&1 || die "git could not be found in PATH."

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
MFdir="$(realpath -s -- "$SCRIPT_DIR/../..")"
INSTdir="$MFdir/helpers/install"
DBdir="$MFdir/data/DBs"
INSTALL_MARKER="$INSTdir/runningInstall.sto"

# shellcheck disable=SC2329  # Invoked by the EXIT trap.
cleanup() {
	local status=$?
	if [[ -n "$INSTALL_MARKER" && -e "$INSTALL_MARKER" ]]; then
		if ((status == 0)); then
			rm -f -- "$INSTALL_MARKER"
		else
			echo "Installation failed; leaving $INSTALL_MARKER in place." >&2
		fi
	fi
}

# shellcheck disable=SC2329  # Invoked by the ERR trap.
report_error() {
	local status=$?
	echo "Installer command failed at line $1 with status $status." >&2
	return "$status"
}

trap cleanup EXIT
trap 'report_error "$LINENO"' ERR

env_exists() {
	local target=$1
	local environment_list
	environment_list="$("$MAMBA_E" env list)" || die "Unable to list micromamba environments."
	awk -v target="$target" 'NR > 1 && $1 == target { found=1 } END { exit !found }' <<<"$environment_list"
}

find_in_bashrc() {
	local marker=$1
	[[ -f "$HOME/.bashrc" ]] && grep -Fq -- "$marker" "$HOME/.bashrc"
}

ensure_environment() {
	local name=$1
	local definition=$2
	if env_exists "$name"; then
		echo "Updating $name environment"
		"$MAMBA_E" install --name "$name" --channel-priority flexible -q -y -f "$definition"
	else
		echo "Creating $name environment"
		"$MAMBA_E" create --name "$name" --channel-priority flexible -q -y -f "$definition"
	fi
}

database_current() {
	local marker=$1
	local expected=$2
	((REFRESH_DATABASES == 0)) && [[ -f "$marker" ]] && [[ "$(<"$marker")" == "$expected" ]]
}

echo "MATAFILER4 installer script"
echo "Using micromamba version: $("$MAMBA_E" --version)"

rm -f -- "$INSTdir/progsChecked.sto"
touch -- "$INSTALL_MARKER"
mkdir -p -- "$MFdir/gits" "$DBdir"

if [[ -L "$MFdir/config.txt" && ! -e "$MFdir/config.txt" ]]; then
	die "$MFdir/config.txt is a broken symbolic link; repair or remove it before rerunning the installer."
fi

if [[ ! -e "$MFdir/config.txt" ]]; then
	if [[ ! -e "$MFdir/Mods/MATAFILERcfg.txt" ]]; then
		cp -- "$MFdir/Mods/config.old" "$MFdir/Mods/MATAFILERcfg.txt"
	fi
	ln -s -- "$MFdir/Mods/MATAFILERcfg.txt" "$MFdir/config.txt"
	echo "Created config.txt. Please review its local paths."
fi

if find_in_bashrc "##------------> MG-TK ADDED"; then
	die "$HOME/.bashrc contains an obsolete MG-TK block. Remove that block and rerun the installer."
fi

if ! find_in_bashrc "##------------> MF4 ADDED"; then
	touch -- "$HOME/.bashrc"
	# shellcheck disable=SC2016  # These variables must remain literal in .bashrc.
	printf '\n\n##------------> MF4 ADDED <----------##\nexport MF4DIR=%q\nexport PERL5LIB="${PERL5LIB:+${PERL5LIB}:}${MF4DIR}"\n##------------> MF4 ADDED <----------##\n\n' \
		"$MFdir/" >> "$HOME/.bashrc"
	echo "Added MATAFILER4 paths to ~/.bashrc"
fi

legacy_envs=(MGTK MGTK_R MGTKbinners MGTKcheckm2 MGTKgtdbtk MGTKphylo MGTKsemibin MGTKwhokar)
if ((REMOVE_LEGACY_ENVS)); then
	for legacy_env in "${legacy_envs[@]}"; do
		if env_exists "$legacy_env"; then
			echo "Removing obsolete environment $legacy_env"
			"$MAMBA_E" env remove -y -n "$legacy_env"
		fi
	done
else
	for legacy_env in "${legacy_envs[@]}"; do
		if env_exists "$legacy_env"; then
			echo "Obsolete environment $legacy_env was left untouched; use --remove-legacy-envs to remove it."
		fi
	done
fi

export PIP_USER=false

ensure_environment MF4 "$INSTdir/MF4.yml"

if "$MAMBA_E" run -n MF4 hostile --help >/dev/null 2>&1; then
	echo "Installing/updating the Hostile human reference database"
	retry_command "Hostile database download" 5 15 \
		env HOSTILE_CACHE_DIR="$DBdir/hostile" \
		"$MAMBA_E" run -n MF4 hostile index fetch --name human-t2t-hla
fi

XGTDB_DIR="$MFdir/gits/XGTDB"
XGTDB_REVISION="0e944a700b96bb7b31f12f03ece969d16bcd1a10"
if [[ ! -e "$XGTDB_DIR" ]]; then
	echo "Installing extract_gtdb_mg at revision $XGTDB_REVISION"
	git clone --no-checkout https://github.com/4less/extract_gtdb_mg.git "$XGTDB_DIR"
	git -C "$XGTDB_DIR" checkout --detach "$XGTDB_REVISION"
elif [[ -d "$XGTDB_DIR/.git" ]]; then
	installed_revision="$(git -C "$XGTDB_DIR" rev-parse HEAD)"
	[[ "$installed_revision" == "$XGTDB_REVISION" ]] || \
		die "$XGTDB_DIR is at $installed_revision, expected $XGTDB_REVISION. Preserve any local work, then replace or update this checkout."
	git -C "$XGTDB_DIR" diff --quiet --ignore-submodules -- || \
		die "$XGTDB_DIR contains local changes; restore or preserve them before rerunning the installer."
	git -C "$XGTDB_DIR" diff --cached --quiet --ignore-submodules -- || \
		die "$XGTDB_DIR contains staged local changes; restore or preserve them before rerunning the installer."
else
	die "$XGTDB_DIR exists but is not a Git checkout. Move it aside and rerun the installer."
fi

ensure_environment MF4gtdbtk "$INSTdir/GTDBTK.yml"
ensure_environment MF4semibin "$INSTdir/SemiBin.yml"
ensure_environment MF4binners "$INSTdir/Binners.yml"
ensure_environment MF4genomeface "$INSTdir/MF4genomeface.yml"
ensure_environment MF4scgbinner "$INSTdir/SCGBinner.yml"
ensure_environment MF4checkm2 "$INSTdir/checkm2.yml"

CHECKM2_VERSION="1.0.2"
METAPHLAN_VERSION="4.1"
CM2DB="$DBdir/CM2"
MP4DB="$DBdir/MP4"
CM2_DIAMOND_DB="$CM2DB/CheckM2_database/uniref100.KO.1.dmnd"
CM2_MARKER="$CM2DB/.mf4-checkm2-version"
MP4_MARKER="$MP4DB/.mf4-metaphlan-version"

if ! database_current "$CM2_MARKER" "$CHECKM2_VERSION"; then
	if ((REFRESH_DATABASES == 0)) && [[ -f "$CM2_DIAMOND_DB" ]] && \
		"$MAMBA_E" run -n MF4checkm2 checkm2 database --setdblocation "$CM2_DIAMOND_DB"; then
		echo "Using existing CheckM2 database at $CM2_DIAMOND_DB"
	else
		echo "Installing CheckM2 $CHECKM2_VERSION database"
		mkdir -p -- "$CM2DB"
		retry_command "CheckM2 database download" 5 15 \
			"$MAMBA_E" run -n MF4checkm2 checkm2 database --download --path "$CM2DB"
	fi
	printf '%s\n' "$CHECKM2_VERSION" > "$CM2_MARKER"
fi

if ! database_current "$MP4_MARKER" "$METAPHLAN_VERSION"; then
	echo "Installing MetaPhlAn $METAPHLAN_VERSION database"
	mkdir -p -- "$MP4DB"
	retry_command "MetaPhlAn database download" 5 15 \
		"$MAMBA_E" run -n MF4checkm2 metaphlan --install --bowtie2db "$MP4DB"
	printf '%s\n' "$METAPHLAN_VERSION" > "$MP4_MARKER"
fi

ensure_environment MF4phylo "$INSTdir/phylo.yml"
ensure_environment MF4_R "$INSTdir/MGTK_R.yml"

echo
echo "How to download GTDB and GTDB-Tk databases"
echo
echo "These databases are required for MAG classification. For example:"
echo "    helpers/install/get_gtdb.py all -v 226 -t /path/to/download -d /path/to/extract/to --tk split"
echo "Run 'helpers/install/get_gtdb.py -h' for all options."
echo
echo "Finished MATAFILER4 installation."
echo "Activate the main environment with: micromamba activate MF4"
echo "Then verify the installation with: ./MATAF4.pl -checkInstall"

exit 0
