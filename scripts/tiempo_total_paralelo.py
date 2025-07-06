import matplotlib.pyplot as plt
# ================================
# Datos de tiempo de cómputo puro
# ================================
procesadores = [1, 2, 4, 8, 16, 32]
tamanios = [100000, 1000000, 10000000, 100000000, 1000000000]

# =============================
# Estilos visuales del gráfico
# =============================
estilos = ['-o', '-s', '-^', '-D', '-v', '-x']
colores = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

# Datos de tiempo total paralelo con comunicación bloqueante
tiempo_par_block = {
    100000:     [0.518, 0.307, 0.516, 0.480, 1.420, 7.112],
    1000000:    [4.864, 4.795, 4.595, 5.408, 9.272, 9.501],
    10000000:   [47.444, 47.963, 46.378, 48.701, 52.144, 61.342],
    100000000:  [471.237, 483.723, 477.490, 491.209, 511.091, 623.381],
    1000000000: [4808.947, 4810.780, 4804.336, 5049.901, 5212.049, 6443.332],
}

# Datos de tiempo total paralelo con comunicación no bloqueante
tiempo_par_nonblock = {
    100000:     [0.518, 0.443, 0.415, 0.503, 0.730, 0.789],
    1000000:    [4.864, 4.394, 4.080, 4.201, 7.474, 7.810],
    10000000:   [47.444, 43.122, 39.636, 39.986, 41.770, 51.354],
    100000000:  [471.237, 429.345, 401.607, 408.363, 433.673, 508.978],
    1000000000: [4808.947, 4216.784, 4002.156, 4138.011, 4218.533, 5092.162],
}

# Crear gráfico comparando MPI bloqueante y no bloqueante
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_par_block[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'Block - N={N}', linewidth=2, markersize=6)
    plt.plot(procesadores, tiempo_par_nonblock[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], linestyle='--', label=f'NonBlock - N={N}', linewidth=2, markersize=6)

plt.title('Tiempo Total Paralelo: MPI Bloqueante vs No Bloqueante', fontsize=14)
plt.xlabel('Número de Procesadores (P)', fontsize=12)
plt.ylabel('Tiempo Total (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tipo y Tamaño del problema", fontsize=9, title_fontsize=10, loc='upper right')
plt.tight_layout()
plt.savefig('tiempo_computo_mpis.png')
plt.show()
