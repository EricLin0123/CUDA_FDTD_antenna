import numpy as np
import cupy as cp
import time

# CPU operation
x_cpu = np.random.rand(100000)
start_cpu = time.time()
l2_cpu = np.linalg.norm(x_cpu)
end_cpu = time.time()
cpu_time = end_cpu - start_cpu

# GPU operation
x_gpu = cp.random.rand(100000)
start_gpu = time.time()
l2_gpu = cp.linalg.norm(x_gpu)
cp.cuda.Stream.null.synchronize()  # Ensure all GPU operations are finished
end_gpu = time.time()
gpu_time = end_gpu - start_gpu

print(f"CPU time: {cpu_time} seconds, L2 norm: {l2_cpu}")
print(f"GPU time: {gpu_time} seconds, L2 norm: {l2_gpu}")
