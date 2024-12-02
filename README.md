MolSim - Group A
===

## Dependencies
- Cmake 3.10
- Doygen 1.9.8 (`sudo apt install doxygen`)
- Libxerces (`sudo apt install libxerces-c-dev`)

## Build
### Configuration
- Install
  ```bash
  git clone https://github.com/jkr11/MolSim.git
  ```
- Build the project using the provided build script by using source, add `-t` to also build and run tests, add `-b` to run the benchmark
  ```bash
  cd MolSim/scripts
  source build <CMAKE_BUILD_TYPE= Release (default) | Debug | asan | asan-quiet>  [-t|--test] [-b|--benchmark]
  ```
- Set the Input file by selecting the corresponding number during the script execution
  ```bash
  source set-input
  ```

- Creating documentation when doxygen is installed (has to be executed in the specific `buildDir/<CMAKE_BUILD_TYPE>`)
  ```bash
  (cd ../buildDir/<CMAKE_BUILD_TYPE> when starting from /scripts)
  make doc_doxygen 
  ```
- Running the program
  ```bash
  $BUILD -f $INPUT <options>
  ``` 
- `$BUILD` contains the location of the last compiled executable
- `$INPUT` contains the location of the selected input file
- Please note that `$BUILD` and `$INPUT` are only available if the scripts are executed via source.

### Options

  ```console
  Options:
  --help | -h                     Show this help message
  --file | -f <filename>          Specify the input file
  [--step_size | -s <double>]     Specify how often the output will be written wrt. time(step_size), default=1
                                    Note that this is independent of the time resolution (t_delta) and dependent on the simulation time
  [--loglevel | -l <level>]       Specify the log level, default=info, valid=[off, error, warn, info, debug, trace]
  Example usage:
  $BUILD -f $INPUT -l <loglevel> -s <number>
  ```
- Output is located in `./output/<current_time>`
- `--step_size` is relative to the passed simulation time and not the number of iterations
- `--loglevel debug` is only available if compiled with CMAKE_BUILD_TYPE=Debug
- all other options are specified in the .xml input file
- old inputs have been migrated to xml and support this pipeline

## LinkedCells vs DirectSum performance
The LinkedCell implementation is more performant than the old DirectSum implementation.

![Benchmark Graph](benchmark/graph.png)