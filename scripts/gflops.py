import matplotlib.pyplot as plt

# ================================
# Datos base
# ================================
procesadores = [1, 2, 4, 8, 16, 32]

# Tiempos paralelos no bloqueantes para N = 10 millones (en milisegundos)
tiempos_par_nonblock = [47.444, 43.122, 39.636, 39.986, 41.770, 51.354]

# Suponemos 1 GFLOP equivalente para todos los casos (valor relativo)
# GFLOPS = FLOP / (tiempo en segundos)
gflops = [1 / (t / 1000) for t in tiempos_par_nonblock]

# ================================
# Gráfico de GFLOPS corregido
# ================================
plt.figure(figsize=(10, 6))
plt.plot(procesadores, gflops, '-o', label='GFLOPS', color='tab:blue', linewidth=2, markersize=8)

plt.title('GFLOPS vs Número de Procesadores (N = 10 millones)', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('GFLOPS', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('gflops.png')
plt.show()
