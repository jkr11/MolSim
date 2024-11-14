MolSim - Group A
===

# Dependencies
- Cmake 3.10

- Doygen 1.9.8

    ```bash
    sudo apt install doxygen
    ```

- Libxerces
  ```bash
  sudo apt install libxerces-c-dev
  ```

# Build
### Configuration
- Install
  ```bash
  git clone https://github.com/jkr11/MolSim.git
  ```
- Build the project using the provided script by using source
  ```bash
  cd MolSim/scripts
  source build <CMAKE_BUILD_TYPE= Release | Debug | asan | asan-quiet>
  ```
- Creating documentation when doxygen is installed (has to be executed in the specific buildDir <CMAKE_BUILD_TYPE>)
- ```bash
  (cd ../buildDir/<CMAKE_BUILD_TYPE> when starting from /scripts)
  make doc_doxygen 
  ```
- Running the program
  ```bash
  $BUILD <options>
  ``` 
Please note that $BUILD is only available if the script is executed via source and contains the path to the last compiled executable

### Options

```console
  Options:
  --help | -h                     Show this help message
  --file | -f <filename>          Specify the input file
  [--t_end | -t <double>]         Specify the simulation end time (t_end), default=100
  [--delta_t | -d <double>]       Specify the simulation delta time (t_delta), default=0.014
  [--step_size | -s <double>]     Specify how often the output will be written wrt. time(step_size), default=1
                                  Note that this is independent of the time resolution (t_delta) and dependent on the simulation time
  [--loglevel | -l <level>]       Specify the log level, default=info, valid=[off, error, warn, info, debug, trace]
  [--force | -F <forceType>]      Specify what force to use, default=lennardjones, forceType=[lennardjones,gravity]
  [--reader | -R <readerType>]    Specify reader type, default=cuboidreader, readerType=[cuboidreader,defaultreader]Example:
  
  Example Usage Week02:
  $BUILD -f ../input/test.cuboid -t 5 -d 0.0002 -s 0.01 -l info
  
  Example Usage Week01
  $BUILD -f ../input/eingabe-sonne.txt -t 10000 -s 10 --reader defaultreader --force gravity
```

- Output is located in ./output/<current_time>
- --step_size is relative to the passed simulation time and not the number of iterations
- --loglevel debug is only available if compiled with CMAKE_BUILD_TYPE=Debug
- --force gravity and --reader defaultReader are needed to simulate the behaviour of week01 