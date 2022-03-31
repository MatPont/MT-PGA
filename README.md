# PrincipalGeodesicAnalysisMergeTree

Tested on Ubuntu 20.04.4 LTS.

## Install the dependencies

```bash
sudo apt-get install cmake-qt-gui libboost-system-dev libpython3.8-dev libxt-dev libxcursor-dev libopengl-dev
sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
sudo apt-get install python3-sklearn 
```

## Install Paraview

First, go in the root of the extracted archive then run the following commands:
(replace the 4 in "make -j4" by the number of available cores on your system)

```bash
git clone https://github.com/topology-tool-kit/ttk-paraview.git
cd ttk-paraview
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make install
```

Some warnings are expected when using the `make` command, they should not cause any problems.

## Install cmake 3.21 (Optional)

To install TTK (next step), you will need cmake 3.21 or higher.
If your OS does not allow to install automatically the latest version of cmake (like Ubuntu 20) you can install it with the following commands:
(replace the 4 in "make -j4" by the number of available cores on your system)
(if you do not want to install libssl-dev you can add `-DCMAKE_USE_OPENSSL=OFF` at the end of the command `configure`)

```bash
sudo apt-get install libssl-dev 
wget http://www.cmake.org/files/v3.21/cmake-3.21.6.tar.gz
tar xf cmake-3.21.6.tar.gz
cd cmake-3.21.6
./configure --prefix="./install/"
make -j4
make install
```

## Install TTK

Go in the `ttk-dev2` directory then run the following commands:
(replace the 4 in "make -j4" by the number of available cores on your system)
(if you have needed to install cmake with the step above then replace `cmake` by `../../cmake-3.21.6/install/bin/cmake`)

```bash
mkdir build && cd build
paraviewPath=`pwd`/../../ttk-paraview/install/lib/cmake/paraview-5.10
pluginDir=`pwd`/../../ttk-paraview/install/lib/paraview-5.10
cmake -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath -DTTK_INSTALL_PLUGIN_DIR=$pluginDir ..
make -j4
make install
```

## Get the results

Go in the root directory of the unziped archive and extract the data:

```bash
tar xvJf data.tar.xz
```

### table 1

To reproduce the results of the time table in the paper, please go in the `scripts` directory and enter the following commands:

```bash
for f in *.sh; do chmod u+x $f; done
```

Run the experiments (it will take a LONG time) and print table:
(replace N with the number of available cores on your system)

```bash
./automata2.sh N
./timeTable.sh
```
