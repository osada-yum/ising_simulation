set -x -u -e

root_dir="./bin/root_$$"

# srcfile="./src/ising2d_periodic_simple_m.f90"
# progfile="./app/mpi_ising2d_periodic_simple_relaxation.f90"
# execname="mpi_ising2d_periodic_simple_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_table_m.f90"
# progfile="./app/mpi_ising2d_periodic_table_relaxation.f90"
# execname="mpi_ising2d_periodic_table_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_probtable_m.f90"
# progfile="./app/mpi_ising2d_periodic_probtable_relaxation.f90"
# execname="mpi_ising2d_periodic_probtable_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_dual_lattice_probtable_m.f90"
# progfile="./app/mpi_ising2d_periodic_dual_lattice_probtable_relaxation.f90"
# execname="mpi_ising2d_periodic_dual_lattice_probtable_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_dual_lattice_probtable_int8_m.f90"
# progfile="./app/mpi_ising2d_periodic_dual_lattice_probtable_int8_relaxation.f90"
# execname="mpi_ising2d_periodic_dual_lattice_probtable_int8_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_dual_lattice_locality_probtable_m.f90"
# progfile="./app/mpi_ising2d_periodic_dual_lattice_locality_probtable_relaxation.f90"
# execname="mpi_ising2d_periodic_dual_lattice_locality_probtable_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_dual_lattice_locality_probtable_globalparam_m.f90"
# progfile="./app/mpi_ising2d_periodic_dual_lattice_locality_probtable_globalparam_relaxation.f90"
# execname="mpi_ising2d_periodic_dual_lattice_locality_probtable_globalparam_relaxation"
# execfile="${root_dir}/bin/${execname}"

# srcfile="./src/ising2d_periodic_dual_lattice_locality_probtable_globalparam_int8_m.f90"
# progfile="./app/mpi_ising2d_periodic_dual_lattice_locality_probtable_globalparam_int8_relaxation.f90"
# execname="mpi_ising2d_periodic_dual_lattice_locality_probtable_globalparam_int8_relaxation"
# execfile="${root_dir}/bin/${execname}"

srcfile="./src/ising2d_periodic_optim_m.f90"
progfile="./app/mpi_ising2d_periodic_optim_relaxation.f90"
execname="mpi_ising2d_periodic_optim_relaxation"
execfile="${root_dir}/bin/${execname}"

output_dir="data/ising2d"
outputfile="${output_dir}/ising2d_periodic_x${nx}_y${ny}_mcs${mcs}_sample${sample}_kbt${kbt}_iseed${iseed}_skip${n_skip}_${execname}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p "${output_dir}"
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

sed -i -e "s/nx = [0-9]*_int64/nx = ${nx}_int64/" \
    -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
    "${srcfile}"

sed -i -e "s/mcs = [0-9]*_int32/mcs = ${mcs}_int32/" \
    -e "s/nsample = [0-9]*_int32/nsample = ${sample}_int32/" \
    -e "s/n_skip = [0-9]*_int64/n_skip = ${n_skip}_int64/" \
    -e "s/iseed = [0-9]*_int64/iseed = ${iseed}_int64/" \
    "${progfile}"

fpm install "${execname}" --prefix="${root_dir}" --verbose --compiler='mpif90' --flag="${FCFLAGS}"
start=$(date +%s)
mpirun -np 2 ${execfile} > "${tmpfile}"
end=$(date +%s)
cp -v "${tmpfile}" "${outputfile}"
echo "output >>> '${outputfile}'"
chmod 400 "${outputfile}"

hour=$( echo "($end - $start) / 3600" | bc)
minute=$( echo "(($end - $start) % 3600) / 60" | bc)
second=$( echo "($end - $start) % 60" | bc)
minute=$(printf "%02d" "${minute}")
second=$(printf "%02d" "${second}")
elapsed_time="${hour}h ${minute}m ${second}s"
#     model,»  size,»     sample,»     mcs,»   kbt, time
echo "2D-Ising_MET(MT19937, ${execname}),  ${nx}x${ny} == $((nx*ny)), ${sample}, ${mcs}, ${kbt}, iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a ising2d.log

rm -rf "${root_dir:?}"
