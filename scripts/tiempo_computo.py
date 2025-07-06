import matplotlib.pyplot as plt

# ================================
# Datos de tiempo de cómputo puro
# ================================
procesadores = [1, 2, 4, 8, 16, 32]
tamanios = [100000, 1000000, 10000000, 100000000, 1000000000]

tiempo_computo = {
    100000:     [0.518, 0.250, 0.123, 0.106, 0.040, 0.020],
    1000000:    [4.864, 2.398, 1.208, 0.614, 0.370, 0.196],
    10000000:   [47.444, 23.683, 11.950, 6.113, 3.246, 1.844],
    100000000:  [471.237, 236.688, 124.013, 61.147, 32.380, 19.763],
    1000000000: [4808.947, 2314.814, 1191.064, 615.024, 324.560, 196.675],
}

# =============================
# Estilos visuales del gráfico
# =============================
estilos = ['-o', '-s', '-^', '-D', '-v', '-x']
colores = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

plt.figure(figsize=(12, 7))
for i, N in enumerate(tamanios):
    plt.plot(procesadores, tiempo_computo[N], estilos[i % len(estilos)],
             color=colores[i % len(colores)], label=f'N = {N}', linewidth=2, markersize=8)

plt.title('Tiempo de Cómputo Puro (sin comunicación) vs Número de Procesadores', fontsize=14)
plt.xlabel('Número de Procesadores (P)', fontsize=12)
plt.ylabel('Tiempo de Cómputo (ms)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(procesadores)
plt.legend(title="Tamaño del problema (N)", fontsize=10, title_fontsize=11, loc='upper right')
plt.tight_layout()
plt.savefig('tiempo_computo.png')
plt.show()

