import matplotlib.pyplot as plt

# Procesadores y tamaño del problema
procesadores = [1, 2, 4, 8, 16, 32]
N = 10000000

# Speedup para bloqueante y no bloqueante desde los datos proporcionados
speedup_block = [1.0, 1.037, 1.078, 1.050, 1.044, 1.005]
speedup_nonblock = [1.0, 1.154, 1.261, 1.279, 1.304, 1.201]

# Graficar Speedup
plt.figure(figsize=(8, 6))
plt.plot(procesadores, speedup_block, marker='o', label='MPI Bloqueante')
plt.plot(procesadores, speedup_nonblock, marker='o', label='MPI No Bloqueante')
plt.plot(procesadores, procesadores, '--', color='gray', label='Speedup Ideal')

plt.title(f'Speedup vs Número de Procesadores (N={N:,})')
plt.xlabel('Número de Procesadores')
plt.ylabel('Speedup')
plt.xticks(procesadores)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('speedup.png')
plt.show()
