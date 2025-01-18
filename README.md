MolSim - Group A
===

## Dependencies

- Cmake 3.24
- Doygen 1.9.8 (`sudo apt install doxygen`)
- Libxerces (`sudo apt install libxerces-c-dev`)

## Build

### Configuration

- Install
  ```bash
  git clone https://github.com/jkr11/MolSim.git
  ```
- Build the project using the provided build script by using source, add `-t` to also build and run tests, add `-b` to
  enable the BENCHMARK cmake macro
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
  [--checkpoint | -c ]            Specifies if the final particle state will be saved to a checkpoint file.
  Example usage:
  $BUILD -f $INPUT -l <loglevel> -s <number>
  ```

- Output is located in `./output/<current_time>`
- Checkpoint is currently fixed to `checkpoint.xml` also in `./output/`
- `--step_size` is relative to the passed simulation time and not the number of iterations
- `--loglevel debug` is only available if compiled with CMAKE_BUILD_TYPE=Debug
- all other options are specified in the .xml input file
- old inputs have been migrated to xml and support this pipeline

## LinkedCells vs DirectSum performance

### Running benchmark.py

For optimal performance run scripts/build with -b for benchmarking
The python script uses the generated executable in buildDir/Release/src for execution.
Ensure you have python 3.6 or later installed

```bash
cd benchmark
python -m vevn <name>
source <name>/bin/activate
pip3 install argparse matplotlib

python3 benchmark.py <options> 
options:
  -h, --help            show this help message and exit
  -a A_VALUES [A_VALUES ...], --a-values A_VALUES [A_VALUES ...]
                        List of 'a' values (multipliers of 1000 particles) to
                        simulate. Default: [1, 2, 4, 8].
  -o OUTPUT, --output OUTPUT
                        Path to save the output plot. Name it graph.png to
                        push it to git.
  -d, --cached-ds       Use default cached execution times for DirectSum
                        instead of running it.
  -s SAMPLES, --samples SAMPLES
                        Number of runs for each a value. Default = 1
  -c, --cubes           Compare three d files to two d files
  -n, --no-ds           Removes direct sum benchmarks from the plot
```

The LinkedCell implementation is more performant than the old DirectSum implementation.

![Benchmark Graph](benchmark/graph.png)
