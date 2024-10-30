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
- Build remote
  ```bash
  git clone https://github.com/jkr11/MolSim.git
  cd MolSim
  cmake CMakeLists.txt
  ```
- Default
  ```bash
  make
  ```
- With creating documentation when doxygen is installed
- ```bash
  make docs_doxygen 
  ```
- Running the program
  ```bash
  ./MolSim -f <inputfile> -t <t_end> - d <delta_t> -s <output_frequency> -h <help>
  ```