import matplotlib.pyplot as plt

# ================================
# Datos y configuraciones
# ================================
procesadores = [1, 2, 4, 8, 16, 32]
tamanios = [100000, 1000000, 10000000, 100000000, 1000000000]

estilos = ['-o', '-s', '-^', '-D', '-v', '-x']
colores = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

tiempo_computo = {
    100000: [0.518, 0.250, 0.123, 0.106, 0.040, 0.020],
    1000000: [4.864, 2.398, 1.208, 0.614, 0.370, 0.196],
    10000000: [47.444, 23.683, 11.950, 6.113, 3.246, 1.844],
    100000000: [471.237, 236.688, 124.013, 61.147, 32.380, 19.763],
    1000000000: [4808.947, 2314.814, 1191.064, 615.024, 324.560, 196.675]
}

tiempo_comunicacion_block = {
    100000: [0, 0.057, 0.393, 0.373, 1.380, 7.091],
    1000000: [0, 2.397, 3.386, 4.794, 8.903, 9.304],
    10000000: [0, 24.279, 34.429, 42.588, 48.897, 59.498],
    100000000: [0, 247.034, 353.477, 430.062, 478.711, 603.618],
    1000000000: [0, 2495.966, 3613.272, 4434.876, 4887.489, 6246.657],
}

tiempo_comunicacion_nonblock = {
    100000: [0, 0.192, 0.292, 0.396, 0.690, 0.769],
    1000000: [0, 1.996, 2.872, 3.586, 7.105, 7.614],
    10000000: [0, 19.439, 27.686, 33.874, 38.524, 49.510],
    100000000: [0, 192.657, 277.593, 347.216, 401.292, 489.215],
    1000000000: [0, 1901.970, 2811.093, 3522.987, 3893.973, 4895.487],
}

speedup_block = {10000000: [1.0, 1.037, 1.078, 1.050, 1.044, 1.005]}
speedup_nonblock = {10000000: [1.0, 1.154, 1.261, 1.279, 1.304, 1.201]}
gflops = {10000000: [47.444/1000, 43.122/1000, 39.636/1000, 39.986/1000, 41.770/1000, 51.354/1000]}


# ================================
# Gráficos
# ================================

# Tiempo de Cómputo
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_computo[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'N = {N}', linewidth=2, markersize=8)
plt.title('Tiempo de Cómputo vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Cómputo (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('opti_tiempo_computo.png')
plt.show()

# Tiempo de Comunicación Bloqueante
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_comunicacion_block[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'N = {N}', linewidth=2, markersize=8)
plt.title('Tiempo de Comunicación Bloqueante vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Comunicación (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('opti_tiempo_comunicacion_block.png')
plt.show()

# Tiempo de Comunicación No Bloqueante
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_comunicacion_nonblock[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'N = {N}', linewidth=2, markersize=8)
plt.title('Tiempo de Comunicación No Bloqueante vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Comunicación (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11)
plt.tight_layout()
plt.savefig('opti_tiempo_comunicacion_nonblock.png')
plt.show()

# Comparación Comunicación Block vs NonBlock
plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_comunicacion_block[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'Block - N={N}', linewidth=2, markersize=6)
    plt.plot(procesadores, tiempo_comunicacion_nonblock[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], linestyle='--', label=f'NonBlock - N={N}', linewidth=2, markersize=6)
plt.title('Comparación de Tiempo de Comunicación Bloqueante vs No Bloqueante', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Tiempo de Comunicación (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tipo y Tamaño del problema", fontsize=9, title_fontsize=10, loc='upper left')
plt.tight_layout()
plt.savefig('opti_comparacion_tiempo_comunicacion.png')
plt.show()

# Speedup
plt.figure(figsize=(12, 7))
plt.plot(procesadores, speedup_block[10000000], '-o', color='tab:red', label='MPI Bloqueante', linewidth=2, markersize=8)
plt.plot(procesadores, speedup_nonblock[10000000], '-s', color='tab:blue', label='MPI No Bloqueante', linewidth=2, markersize=8)
plt.plot(procesadores, procesadores, '--', color='gray', label='Speedup Teórico', linewidth=1.5)
plt.title('Speedup vs Número de Procesadores (N=10M)', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Speedup (sin unidades)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('opti_speedup.png')
plt.show()

# GFLOPS
plt.figure(figsize=(12, 7))
plt.plot(procesadores, gflops[10000000], '-o', color='tab:green', label='GFLOPS', linewidth=2, markersize=8)
plt.title('GFLOPS vs Número de Procesadores (N=10M)', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('GFLOPS (GigaFLOPS)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('opti_gflops.png')
plt.show()


# ================================
# Gráfico de Eficiencia (Datos Reales)
# ================================

# Procesadores disponibles
procesadores = [1, 2, 4, 8, 16, 32]

# Eficiencia para N = 10,000,000
eficiencia_block_10M = [1.000, 0.942, 0.910, 0.908, 0.832, 0.392]
eficiencia_nonblock_10M = [1.000, 0.983, 0.982, 0.985, 0.958, 0.462]

# Mostrar en consola
print("\n=== Eficiencia para N = 10,000,000 ===")
for p, e_b, e_nb in zip(procesadores, eficiencia_block_10M, eficiencia_nonblock_10M):
    print(f"P={p:2} | Eficiencia Block: {e_b:.3f} | Eficiencia NonBlock: {e_nb:.3f}")

# Gráfico
plt.figure(figsize=(12, 7))
plt.plot(procesadores, eficiencia_block_10M, '-o', label='MPI Bloqueante', color='tab:red', linewidth=2, markersize=8)
plt.plot(procesadores, eficiencia_nonblock_10M, '-s', label='MPI No Bloqueante', color='tab:blue', linewidth=2, markersize=8)

plt.title('Eficiencia Paralela vs Número de Procesadores (N = 10⁷)', fontsize=14)
plt.xlabel('Número de Procesadores', fontsize=12)
plt.ylabel('Eficiencia', fontsize=12)
plt.ylim(0, 1.1)
plt.xticks(procesadores)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('opti_eficiencia.png')
plt.show()
