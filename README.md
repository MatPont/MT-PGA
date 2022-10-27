# MT-PGA

This github repository contains the exact code used for the reference below.

## Reference

If you plan to use this code to generate results for a scientific document, thanks for referencing the following publication:

"Principal Geodesic Analysis of Merge Trees (and Persistence Diagrams)"  
Mathieu Pont, Jules Vidal, Julien Tierny  
IEEE Transactions on Visualization and Computer Graphics, 2022.  

[Paper](https://arxiv.org/pdf/2207.10960.pdf)

## Installation Note

Tested on Ubuntu 20.04.4 LTS.

### Install the dependencies

```bash
sudo apt-get install cmake-qt-gui libboost-system-dev libpython3.8-dev libxt-dev libxcursor-dev libopengl-dev
sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev libqt5svg5-dev qtxmlpatterns5-dev-tools 
sudo apt-get install python3-sklearn 
sudo apt-get install libsqlite3-dev 
```

### Install cmake 3.21 (Optional)

To install TTK (one of the next step), you will need cmake 3.21 or higher.
If your OS does not allow to install automatically the latest version of cmake (like Ubuntu 20) you can install it with the following commands:  
(while being in the parent directory of where this repository (`MT-PGA`) is)  
(replace the `4` in `make -j4` by the number of available cores on your system)  
(if you do not want to install libssl-dev you can add `-DCMAKE_USE_OPENSSL=OFF` at the end of the command `configure`)

```bash
sudo apt-get install libssl-dev 
wget http://www.cmake.org/files/v3.21/cmake-3.21.6.tar.gz
tar xf cmake-3.21.6.tar.gz
cd cmake-3.21.6
./configure --prefix="./install/"
make -j4
make -j4 install
```

In all the following steps, please replace the command `cmake` by the binary `PATHTO/cmake-3.21.6/install/bin/cmake` (replace `PATHTO` by the absolute path to the `cmake-3.21.6` folder).

### Install Paraview

First, go in the parent directory of where this repository (`MT-PGA`) is and run the following commands:
(replace the `4` in `make -j4` by the number of available cores on your system)

```bash
git clone https://github.com/topology-tool-kit/ttk-paraview.git
cd ttk-paraview
git checkout 5.10.0
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_USE_PYTHON=ON -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make -j4 install
```

Some warnings are expected when using the `make` command, they should not cause any problems.

### Install TTK

Go in the `ttk-dev2` directory then run the following commands:
(replace the `4` in `make -j4` by the number of available cores on your system)

```bash
mkdir build && cd build
paraviewPath=`pwd`/../../../ttk-paraview/install/lib/cmake/paraview-5.9
cmake -DCMAKE_INSTALL_PREFIX=../install -DParaView_DIR=$paraviewPath ..
make -j4
make -j4 install
```

Stay in the build directory and set the environment variables:
(replace `3.8` in `python3.8` by your version of python)

```bash
TTK_PREFIX=`pwd`/../install
export PV_PLUGIN_PATH=$TTK_PREFIX/bin/plugins/TopologyToolKit
export LD_LIBRARY_PATH=$TTK_PREFIX/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$TTK_PREFIX/lib/python3.8/site-packages
```

### Get the results

Go in the root directory of the unziped archive and extract the data:

```bash
tar xvJf data.tar.xz
```

#### table 1

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
