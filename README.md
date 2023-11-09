# OpenLPTM

A simple *lumped parameter thermal model* in pure **C**.

## Usage

To compile the library and the demo:

```bash
mkdir build
cd build
meson .. --prefix=/usr
ninja 
```

To install:
```bash
sudo ninja install
```

### Demo

Enter into the `build` folder then to execute the demo simply: `demo`  
To generate the `CSV` file with the evolution of the system simply: `demo output.csv` 
You can also generate a binary file using the follwing syntax: `demo -b output.bin`

### Library

All the relevant interface is stored in `ltm.h`.
