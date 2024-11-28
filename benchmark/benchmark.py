import subprocess
import time
import matplotlib.pyplot as plt
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
executable = os.path.join(script_dir, "../build/src/MolSim")
input_dir = os.path.join(script_dir, "../input")
input_file_template_lc = "task32_{}k_lc.xml"
input_file_template_ds = "task32_{}k_ds.xml"
default_a_values = [1, 2, 4, 8]

if not os.path.exists(executable):
    raise FileNotFoundError(f"Executable not found at {executable}")


def run_ds(a_values, execution_times_ds=None):
    if execution_times_ds is None:
        execution_times_ds = []
    for a in a_values:
        print(f"Running for a = {a}...")
        input_file = os.path.join(input_dir, input_file_template_ds.format(a))
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found at {input_file}")

        start_time = time.perf_counter()
        try:
            subprocess.run([executable, "-f", input_file, str(a), "--loglevel", "off"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error while running the executable for a = {a}: {e}")
            execution_times_ds.append(None)
            continue
        end_time = time.perf_counter()

        elapsed_time = end_time - start_time
        execution_times_ds.append(elapsed_time)
        print(f"Execution time for a = {a}: {elapsed_time:.2f} seconds")
    return execution_times_ds


def run_lc(a_values, execution_times_lc=None):
    if execution_times_lc is None:
        execution_times_lc = []
    for a in a_values:
        print(f"Running for a = {a}...")
        input_file = os.path.join(input_dir, input_file_template_lc.format(a))
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found at {input_file}")

        start_time = time.perf_counter()
        try:
            subprocess.run([executable, "-f", input_file, str(a), "--loglevel", "off"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error while running the executable for a = {a}: {e}")
            execution_times_lc.append(None)
            continue
        end_time = time.perf_counter()

        elapsed_time = end_time - start_time
        execution_times_lc.append(elapsed_time)
        print(f"Execution time for a = {a}: {elapsed_time:.2f} seconds")
    return execution_times_lc


# execution_times_ds = [14.52, 50.81, 192.51, 771.44]  # this is for build with dt = 0.0005 and t = 1


def run_and_plot(a_values):
    #ed = run_ds(a_values)
    ed = [14.52,50.81,192.51,771.44]
    el = run_lc(a_values)
    plt.figure(figsize=(8, 6))
    plt.plot([a * 1000 for a in a_values],
             ed,
             marker='o',
             label="Execution Time DirectSum")
    plt.plot([a * 1000 for a in a_values],
             el,
             marker='o',
             label="Execution Time LinkedCells")
    plt.title(
        "Execution Time vs number of particles running on Intel i7 13700H, 32Gb")
    plt.xlabel("Particles")
    plt.ylabel("Execution Time (seconds)")
    #plt.xscale('log',base=2)
    #plt.yscale('log',base=2)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    output_plot = os.path.join(script_dir, "execution_times.png")
    plt.savefig(output_plot)
    print(f"Plot saved to {output_plot}")
    plt.show()


def run_and_plot_lc(a_values):
    plt.figure(figsize=(8, 6))
    el = run_lc(a_values)
    plt.plot([a * 1000 for a in a_values],
             el,
             marker='o',
             label="Execution Time LinkedCells")
    plt.xlabel("Particles")
    plt.ylabel("Execution Time (seconds)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    output_plot = os.path.join(script_dir, "lctimes.png")
    plt.savefig(output_plot)
    plt.show()


# run_and_plot_lc([1, 2, 4, 8, 16])
run_and_plot([1,2,4,8])