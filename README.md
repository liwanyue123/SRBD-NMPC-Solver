# Manual-Construction-of-an-SRBD-NMPC-Solver

## Requirements
 
### Install Eigen
```bash
sudo apt install libeigen3-dev
```

### Install BLASFEO
For more information, visit the [HPIPM GitHub repository](https://github.com/giaf/hpipm).
```C++
git clone https://github.com/giaf/blasfeo.git 
cd blasfeo
mkdir build 
cd build
cmake ..
make -j8
sudo make install
```

### Install HPIPM
```C++
git clone https://github.com/giaf/hpipm.git
cd hpipm
mkdir build 
cd build
cmake ..
make -j8
sudo make install
```
## Usage
### CMake
```C++
mkdir build
cd build
cmake ..
make -j8
./SRBD_NMPC
```
### VScode
```bash
Compile: Ctrl+shift+B
Run: F5
```
