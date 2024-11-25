import subprocess
import time
import matplotlib.pyplot as plt
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
executable = os.path.join(script_dir, "../build/src/MolSim")
input_dir = os.path.join(script_dir, "../input")
input_file_template = "task32_{}k.xml"
input_file_template_ds = "task32_{}k_ds.xml"
a_values = [1, 2]
execution_times = []

if not os.path.exists(executable):
    raise FileNotFoundError(f"Executable not found at {executable}")

for a in a_values:
    print(f"Running for a = {a}...")
    input_file = os.path.join(input_dir, input_file_template_ds.format(a))
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found at {input_file}")

    start_time = time.perf_counter()  # Start timing
    try:
        subprocess.run([executable, "-f", input_file, str(a)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running the executable for a = {a}: {e}")
        execution_times.append(None)
        continue
    end_time = time.perf_counter()

    elapsed_time = end_time - start_time
    execution_times.append(elapsed_time)
    print(f"Execution time for a = {a}: {elapsed_time:.2f} seconds")

plt.figure(figsize=(8, 6))
plt.plot(a_values*1000, execution_times, marker='o', label="Execution Time")
plt.title("Execution Time vs number of particles")
plt.xlabel("Particles")
plt.ylabel("Execution Time (seconds)")
plt.grid(True)
plt.legend()
plt.tight_layout()

output_plot = os.path.join(script_dir, "execution_times.png")
plt.savefig(output_plot)
print(f"Plot saved to {output_plot}")
plt.show()

