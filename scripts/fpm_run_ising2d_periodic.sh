set -x -u -e

FCFLAGS="-O3"
# FCFLAGS="${FCFLAGS} -g -fbacktrace -Wall -Wextra"
nx=1000
ny=${nx}
mcs=1000
sample=25
kbt=2.269185314213022
iseed=42
n_skip=0

script_dir="scripts"
source "${script_dir}/fpm_run_ising2d_periodic_core.sh"
