import argparse
import subprocess
import time
import matplotlib.pyplot as plt
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
executable = os.path.join(script_dir, "../buildDir/Release/src/MolSim")  # maybe we should add Benchmark type to script
input_dir = os.path.join(script_dir, "../input")
input_file_template_lc = "task32_{}k_lc.xml"
input_file_template_ds = "task32_{}k_ds.xml"
default_a_values = [1, 2, 4, 8]
default_ds_values = [14.52, 50.81, 192.51, 771.44]

if not os.path.exists(executable):
    raise FileNotFoundError(f"Executable not found at {executable}")


def compute_weighted_average(values: list[float]):
    if not values:
        return None
    valid_values = [v for v in values if v is not None]
    if not valid_values:
        return None
    return sum(valid_values) / len(valid_values)


def run_ds(a_values: list[int], s: int = 1, execution_times_ds: list[float] = None):
    if execution_times_ds is None:
        execution_times_ds = []
    for a in a_values:
        print(f"Running LinkedCells for {a}k particles ...")
        input_file = os.path.join(input_dir, input_file_template_ds.format(a))
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found at {input_file}")

        run_times = []
        for run in range(s):
            print(f"  Run {run + 1}/{s}...")
            start_time = time.perf_counter()
            try:
                subprocess.run([executable, "-f", input_file, str(a), "--loglevel", "off"], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running the executable for a = {a}, run {run + 1}: {e}")
                run_times.append(None)
                continue
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            run_times.append(elapsed_time)
            print(f"    Execution time: {elapsed_time:.2f} seconds")

        # Compute weighted average for this `a` value
        avg_time = compute_weighted_average(run_times)
        execution_times_ds.append(avg_time)
        print(f"Weighted average execution time for a = {a}: {avg_time:.2f} seconds" if avg_time else "No valid runs.")
    return execution_times_ds


def run_lc(a_values: list[int], s: int = 1, execution_times_lc: list[float] = None):
    if execution_times_lc is None:
        execution_times_lc = []
    for a in a_values:
        print(f"Running DirectSum for {a}k particles ...")
        input_file = os.path.join(input_dir, input_file_template_lc.format(a))
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found at {input_file}")

        run_times = []
        for run in range(s):
            print(f"  Run {run + 1}/{s}...")
            start_time = time.perf_counter()
            try:
                subprocess.run([executable, "-f", input_file, str(a), "--loglevel", "off"], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error while running the executable for a = {a}, run {run + 1}: {e}")
                run_times.append(None)
                continue
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            run_times.append(elapsed_time)
            print(f"    Execution time: {elapsed_time:.2f} seconds")

        # Compute weighted average for this `a` value
        avg_time = compute_weighted_average(run_times)
        execution_times_lc.append(avg_time)
        print(f"Weighted average execution time for a = {a}: {avg_time:.2f} seconds" if avg_time else "No valid runs.")
    return execution_times_lc


# execution_times_ds = [14.52, 50.81, 192.51, 771.44]  # this is for build with dt = 0.0005 and t = 1


def plot_results(a_values: list[int], execution_times_ds: list[float], execution_times_lc: list[float], samples: int,
                 output_path: str = "execution_times"):
    assert default_a_values == a_values
    plt.figure(figsize=(8, 6))
    if execution_times_ds is not None:
        plt.plot([a * 1000 for a in a_values],
                 execution_times_ds,
                 marker='o',
                 label="Execution Time DirectSum")
    if execution_times_lc is not None:
        plt.plot([a * 1000 for a in a_values],
                 execution_times_lc,
                 marker='o',
                 label="Execution Time LinkedCells")
    plt.title(
        "Execution Time vs Number of Particles on Intel i7-13700H, 32GB RAM\nDomain = (300,300,1), $r_{cutoff} = 3.0$")
    plt.xlabel("Particles")
    plt.ylabel(f"average Execution Time (seconds) over {samples} runs")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")
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


def main():
    parser = argparse.ArgumentParser(
        description="MolSim benchmark."
    )
    parser.add_argument(
        "-a", "--a-values", type=int, nargs="+", default=default_a_values,
        help="List of 'a' values (multipliers of 1000 particles) to simulate. Default: [1, 2, 4, 8]."
    )
    parser.add_argument(
        "-o", "--output", type=str, default=os.path.join(script_dir, "execution_times.png"),
        help="Path to save the output plot. Name it graph.png to push it to git."
    )
    parser.add_argument(
        "-d", "--cached-ds", action="store_true",
        help="Use default cached execution times for DirectSum instead of running it."
    )
    parser.add_argument(
        "-s", "--samples", type=int, default=1,
        help="Number of runs for each a value. Default = 1"
    )

    args = parser.parse_args()

    a_values = args.a_values
    output_path = args.output
    samples = args.samples

    if args.cached_ds:
        execution_times_ds = default_ds_values
    else:
        execution_times_ds = run_ds(a_values, s=samples)

    execution_times_lc = run_lc(a_values, s=samples)

    plot_results(a_values, execution_times_ds, execution_times_lc, samples, output_path)


if __name__ == "__main__":
    main()
