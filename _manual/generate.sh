#!/bin/bash

BITPIT_SOURCE_DIR="$1"
LOCAL_BUILD="$2"

MANUAL_ROOT_DIR="${PWD}"
JEKYLL_ROOT_DIR="${MANUAL_ROOT_DIR}/.."

MANUAL_BUILD_DIR="${MANUAL_ROOT_DIR}/build"

# Generate build directory
rm -rf "${MANUAL_BUILD_DIR}"
mkdir "${MANUAL_BUILD_DIR}"
mkdir "${MANUAL_BUILD_DIR}/_doxygen_templates"
cd "${MANUAL_BUILD_DIR}"

# Link files needed by jekyll
cp -r "${MANUAL_ROOT_DIR}/_includes" "${MANUAL_BUILD_DIR}"
cp -r "${MANUAL_ROOT_DIR}/_layouts" "${MANUAL_BUILD_DIR}"
cp -r "${MANUAL_ROOT_DIR}/stylesheets" "${MANUAL_BUILD_DIR}"
cp "${JEKYLL_ROOT_DIR}/_includes/"* "${MANUAL_BUILD_DIR}/_includes/"
cp "${JEKYLL_ROOT_DIR}/_layouts/"* "${MANUAL_BUILD_DIR}/_layouts/"
cp "${MANUAL_ROOT_DIR}/doxygen_footer.html" "${MANUAL_BUILD_DIR}"
cp "${MANUAL_ROOT_DIR}/doxygen_header.html" "${MANUAL_BUILD_DIR}"

# Configration files needed by jekyll
JEKYLL_CONFIG_GLOBAL="${JEKYLL_ROOT_DIR}/_config.yml"
JEKYLL_CONFIG_LOCAL="${JEKYLL_ROOT_DIR}/_config.local.yml"
JEKYLL_CONFIG_MANUAL="${MANUAL_BUILD_DIR}/_config.manual.yml"

JEKYLL_CONFIG_MANUAL_CONTENTS="""
layouts_dir:  ${MANUAL_BUILD_DIR}/_layouts

exclude: [generate.sh]

"""
echo "${JEKYLL_CONFIG_MANUAL_CONTENTS}" > "${JEKYLL_CONFIG_MANUAL}"

JEKYLL_CONFIG_LIST="${JEKYLL_CONFIG_GLOBAL}"
if [ -f "${JEKYLL_CONFIG_LOCAL}" -a ${LOCAL_BUILD} -eq 1 ]; then
    JEKYLL_CONFIG_LIST="${JEKYLL_CONFIG_LIST},${JEKYLL_CONFIG_LOCAL}"
fi
JEKYLL_CONFIG_LIST="${JEKYLL_CONFIG_LIST},${JEKYLL_CONFIG_MANUAL}"

# Generate templates needed by doxygen
cd ${JEKYLL_ROOT_DIR}

bundle exec jekyll build \
    --config "${JEKYLL_CONFIG_LIST}" \
    --source ${MANUAL_BUILD_DIR} \
    --destination ${MANUAL_BUILD_DIR}/doxygen_templates

# Generate doxygen files
cd ${MANUAL_BUILD_DIR}
mkdir bitpit
cd bitpit

cmake ${BITPIT_SOURCE_DIR} \
    -DCMAKE_INSTALL_PREFIX:PATH="" \
    -DBITPIT_ENABLE_DOC="ON" \
    -DBITPIT_ENABLE_MPI="ON" \
    -DBLAS_VENDOR="Generic" \
    -DLAPACKE_STATIC="ON" \
    -DPETSC_ARCH="" \
    -DPETSC_DIR="/usr/lib64/mpi/gcc/openmpi4/lib64/petsc/3.14.5/linux-gnu-c-opt" \
    -DMPI_INCLUDE_PATH="/usr/lib64/mpi/gcc/openmpi4/include/" \
    -DMPIEXEC_EXECUTABLE="/usr/lib64/mpi/gcc/openmpi4/bin/mpiexec" \
    -DMPI_CXX_COMPILER="/usr/lib64/mpi/gcc/openmpi4/bin/mpicxx" \
    -DMPI_CXX_HEADER_DIR="/usr/lib64/mpi/gcc/openmpi4/include" \
    -DMPI_CXX_LIB_NAMES="mpi_mpifh;mpi" \
    -DMPI_CXX_LINK_FLAGS="-pthread" \
    -DMPI_mpi_LIBRARY="/usr/lib64/mpi/gcc/openmpi4/lib64/libmpi_mpifh.so;/usr/lib64/mpi/gcc/openmpi4/lib64/libmpi.so" \
    -DMPI_mpi_mpifh_LIBRARY="/usr/lib64/mpi/gcc/openmpi4/lib64/libmpi_mpifh.so" \
    -DDOXY_HTML_HEADER="${MANUAL_BUILD_DIR}/doxygen_templates/doxygen_header.html" \
    -DDOXY_HTML_FOOTER="${MANUAL_BUILD_DIR}/doxygen_templates/doxygen_footer.html" \
    -DDOXY_HTML_EXTRA_STYLESHEET="${MANUAL_BUILD_DIR}/stylesheets/manual.css" \
    -DDOXY_DISABLE_INDEX="NO" \
    -DDOXY_GENERATE_TREEVIEW="NO"

cd doc

make

make DESTDIR=${MANUAL_BUILD_DIR} install/local

# Copy the generated manual into the website
cd ${MANUAL_BUILD_DIR}

BITPIT_VERSION=`ls doc`
BITPIT_VERSION=${BITPIT_VERSION/bitpit-/}
BITPIT_VERSION=${BITPIT_VERSION/\//}

MANUAL_INSTALL_DIR="${JEKYLL_ROOT_DIR}/documentation/manual/${BITPIT_VERSION}"

rm -rf "${MANUAL_INSTALL_DIR}"
mkdir -p "${MANUAL_INSTALL_DIR}"
cp -r "doc/bitpit-${BITPIT_VERSION}/html/"* "${MANUAL_INSTALL_DIR}"
