# PrincipalGeodesicAnalysisMergeTree

Tested on Ubuntu 20.04.4 LTS.

## Install the dependencies

```bash
sudo apt-get install cmake-qt-gui libboost-system-dev libpython3.8-dev libxt-dev libxcursor-dev libopengl-dev
sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
sudo apt-get install python3-sklearn 
```

## Install Paraview

First, go in the `ttk-paraview` directory then run the following commands:
(replace the 4 in "make -j4" by the number of available cores on your system)

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DPARAVIEW_PYTHON_VERSION=3 -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make install
```

## Install TTK

Go in the `ttk-dev2` directory then run the following commands:
(replace the 4 in "make -j4" by the number of available cores on your system)

```bash
mkdir build && cd build
paraviewPath=`pwd`/../../ttk-paraview/install/lib/cmake/paraview-5.9
cmake -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath ..
make -j4
make install
```

## Get the results

Extract the data:

```bash
tar xvJf data.tar.xz
```

### table 1

To reproduce the results of Table 1 in the paper, please go in the `scripts` directory and enter the following commands:

```bash
for f in *.sh; do chmod u+x $f; done
```

Run the experiments (it will take a LONG time) and print table:
(replace N with the number of available cores on your system)

```bash
./automata2.sh N
./timeTable.sh
```
