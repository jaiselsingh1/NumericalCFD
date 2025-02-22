#!/bin/sh
current_dir="$(cd "$(dirname "$0")" && pwd)"
cd ..
project_root_dir="$(cd "$(dirname "$0")" && pwd)"
cd ${current_dir}

# Installing packages
echo "Verifying packages..."
export PATH=/mingw64/bin:$PATH
TOOLCHAIN_PACKAGES="mingw-w64-ucrt-x86_64-gcc"
REQUIRED_PACKAGES="mingw-w64-ucrt-x86_64-blas mingw-w64-ucrt-x86_64-lapack mingw-w64-ucrt-x86_64-vtk
mingw-w64-x86_64-libjpeg-turbo mingw-w64-x86_64-libtiff mingw-w64-x86_64-zlib mingw-w64-x86_64-gnuplot"

is_installed() {
  pacman -Qi "$1" &>/dev/null
}

if ! is_installed "$TOOLCHAIN_PACKAGES"; then
  echo "mingw-w64-ucrt-x86_64-toolchain not found, installing..."
  pacman -S --needed "mingw-w64-ucrt-x86_64-toolchain"
else
  echo "mingw-w64-ucrt-x86_64-toolchain is installed."
fi

for PACKAGE in $REQUIRED_PACKAGES; do
  if ! is_installed "$PACKAGE"; then
    echo "$PACKAGE not found, installing..."
    pacman -S --needed "$PACKAGE"
  else
    echo "$PACKAGE is installed."
  fi
done

# Installing External Libraries
echo -e "\nVerifying external libraries..."
mkdir -p ${project_root_dir}/third_party
cd ${project_root_dir}/third_party

if [ ! -d "matplotplusplus" ]; then
  echo "matplotplusplus not found, installing..."
  git clone https://github.com/alandefreitas/matplotplusplus
  cd ${project_root_dir}/third_party/matplotplusplus
  rm -rf build install
  mkdir -p build install
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/matplotplusplus/install -DCMAKE_BUILD_TYPE=Release ..
  cmake --build . -j 4
  cmake --install .
fi
echo "matplotplusplus is installed!"

cd ${project_root_dir}/third_party

if [ ! -d "eigen" ]; then
  echo "eigen not found, installing..."
  git clone https://gitlab.com/libeigen/eigen.git
  cd ${project_root_dir}/third_party/eigen
  rm -rf build install
  mkdir -p build install
  cd build
  cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/eigen/install -DCMAKE_BUILD_TYPE=Release ..
  cmake --build . -j 4
  cmake --install .
fi
echo "eigen is installed!"


# Compiling
echo -e "\nCompiling..."
cd ${current_dir}
rm -rf ${current_dir}/build
mkdir -p ${current_dir}/build
cd ${current_dir}/build
cmake ..
cmake --build . -j 4
cmake --install .
cd ${current_dir}
mkdir -p ${current_dir}/results
./build/hw -c config.inp > results/results.txt 2>&1
echo "Build complete!"
echo "Outputs are shown in results/results.txt"