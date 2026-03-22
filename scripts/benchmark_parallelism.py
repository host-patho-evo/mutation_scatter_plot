import time
import subprocess
import os
import multiprocessing

def run_benchmark(threads, input_file, ref_file):
    cmd = [
        "python3", "-m", "mutation_scatter_plot.calculate_codon_frequencies.cli",
        "--alignment-file", input_file,
        "--outfile-prefix", f"/tmp/benchmark_t{threads}",
        "--padded-reference",
        "--reference-infile", ref_file,
        "--threads", str(threads),
        "--overwrite"
    ]
    # Set PYTHONPATH
    env = os.environ.copy()
    env["PYTHONPATH"] = "src"
    
    start = time.time()
    res = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
    end = time.time()
    print(res.stdout)
    return end - start

def main():
    input_base = "tests/inputs/test2_full.fasta"
    ref_file = "tests/inputs/MN908947.3_S_full.fasta"
    bench_file = "/tmp/large_bench.fasta"
    
    # Create a large file with many unique sequences
    print("Generating synthetic dataset (50k unique sequences)...")
    if os.path.exists(bench_file):
        os.remove(bench_file)
    
    with open(bench_file, "w") as out:
        with open(input_base, "r") as f:
            lines = f.readlines()
            header = lines[0].strip()
            seq = lines[1].strip()
            for i in range(50000):
                # Change characters to make it unique
                muted_seq = seq[:i % 1000] + ('G' if seq[i % 1000] == 'A' else 'A') + seq[i % 1000 + 1:]
                out.write(f">{i}_{header[1:]}\n{muted_seq}\n")
    
    print(f"File size: {os.path.getsize(bench_file) / 1024 / 1024:.2f} MB")
    
    cores = [1, 2, 4]
    results = {}
    
    for c in cores:
        print(f"Running benchmark with {c} threads...")
        t = run_benchmark(c, bench_file, ref_file)
        results[c] = t
        print(f"  Result: {t:.2f} seconds")
        
    print("\nParallelization Efficiency Summary:")
    print("-----------------------------------")
    print(f"Threads | Time (s) | Speedup | Efficiency")
    t1 = results[1]
    for c in cores:
        tc = results[c]
        speedup = t1 / tc
        efficiency = speedup / c
        print(f"{c:7d} | {tc:8.2f} | {speedup:7.2f} | {efficiency:10.2f}%")

if __name__ == "__main__":
    main()
