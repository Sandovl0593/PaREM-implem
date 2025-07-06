# Makefile para DFA Paralelo MPI - Versión Simplificada
CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Wno-cast-function-type
TARGET_ORIGINAL = dfa_parallel
TARGET_OPTIMIZED = dfa_optimized
TARGET_METRICS = scalability_metrics

# Número de procesos por defecto
NP ?= 4

# Archivos compartidos (sin main)
SHARED_SOURCES = parallel_dfa.cpp communication_tests.cpp dfa_tests.cpp utils.cpp optimized_parallel_dfa.cpp

# Archivos específicos por versión
MAIN_ORIGINAL = main.cpp
MAIN_OPTIMIZED = optimized_main.cpp
MAIN_METRICS = scalability_metrics.cpp

# Fuentes completas por versión
SOURCES_ORIGINAL = $(MAIN_ORIGINAL) $(SHARED_SOURCES)
SOURCES_OPTIMIZED = $(MAIN_OPTIMIZED) $(SHARED_SOURCES)
SOURCES_METRICS = $(MAIN_METRICS) optimized_parallel_dfa.cpp utils.cpp

# Objetivos principales
all: $(TARGET_ORIGINAL) $(TARGET_OPTIMIZED) $(TARGET_METRICS)

# Versión original (usa main.cpp como punto de entrada)
$(TARGET_ORIGINAL): $(SOURCES_ORIGINAL)
	@echo "🔨 Compilando versión original..."
	$(CXX) $(CXXFLAGS) $(SOURCES_ORIGINAL) -o $(TARGET_ORIGINAL)
	@echo "✅ Versión original compilada: $(TARGET_ORIGINAL)"

# Versión optimizada (usa optimized_main.cpp como punto de entrada)
$(TARGET_OPTIMIZED): $(SOURCES_OPTIMIZED)
	@echo "🔨 Compilando versión optimizada..."
	$(CXX) $(CXXFLAGS) $(SOURCES_OPTIMIZED) -o $(TARGET_OPTIMIZED)
	@echo "✅ Versión optimizada compilada: $(TARGET_OPTIMIZED)"

# Compilar solo la versión original
original: $(TARGET_ORIGINAL)

# Compilar solo la versión optimizada  
optimized: $(TARGET_OPTIMIZED)

# Compilar métricas de escalabilidad
$(TARGET_METRICS): $(SOURCES_METRICS)
	@echo "📊 Compilando generador de métricas de escalabilidad..."
	$(CXX) $(CXXFLAGS) $(SOURCES_METRICS) -o $(TARGET_METRICS)
	@echo "✅ Generador de métricas compilado: $(TARGET_METRICS)"

# Compilar solo las métricas
metrics: $(TARGET_METRICS)

# Verificar archivos fuente
check-sources:
	@echo "📂 Archivos para versión original:"
	@echo "$(SOURCES_ORIGINAL)"
	@echo ""
	@echo "📂 Archivos para versión optimizada:"
	@echo "$(SOURCES_OPTIMIZED)"
	@echo ""
	@echo "🔍 Verificando existencia de archivos..."
	@for file in $(SOURCES_ORIGINAL); do \
		if [ -f $$file ]; then \
			echo "✅ $$file"; \
		else \
			echo "❌ $$file (NO ENCONTRADO)"; \
		fi \
	done

# Limpiar
clean:
	@echo "🧹 Limpiando archivos compilados y resultados..."
	rm -f $(TARGET_ORIGINAL) $(TARGET_OPTIMIZED) $(TARGET_METRICS)
	rm -rf resultados_mpi/
	rm -f MetricsGraphics.txt
	@echo "✅ Limpieza completada"

# Ejecutar versión original
run-original: $(TARGET_ORIGINAL)
	@echo "🔄 Ejecutando versión original con $(NP) procesos..."
	mpirun -np $(NP) ./$(TARGET_ORIGINAL)

# Ejecutar versión optimizada
run-optimized: $(TARGET_OPTIMIZED)
	@echo "⚡ Ejecutando versión optimizada con $(NP) procesos..."
	mpirun -np $(NP) ./$(TARGET_OPTIMIZED)

# Comparar ambas versiones
compare: $(TARGET_ORIGINAL) $(TARGET_OPTIMIZED)
	@echo "⚖️ Comparando versiones con $(NP) procesos..."
	@echo "🔄 Versión original:"
	mpirun -np $(NP) ./$(TARGET_ORIGINAL)
	@echo ""
	@echo "⚡ Versión optimizada:"
	mpirun -np $(NP) ./$(TARGET_OPTIMIZED)

# Generar métricas individuales por número de procesos
generate-metrics-p1: scalability_metrics
	@echo "📈 Generando métricas para P = 1..."
	mpirun -np 1 ./scalability_metrics 1
	@echo "✅ Métricas P=1 completadas"

generate-metrics-p2: scalability_metrics
	@echo "📈 Generando métricas para P = 2..."
	mpirun -np 2 ./scalability_metrics 2
	@echo "✅ Métricas P=2 completadas"

generate-metrics-p4: scalability_metrics
	@echo "📈 Generando métricas para P = 4..."
	mpirun -np 4 ./scalability_metrics 4
	@echo "✅ Métricas P=4 completadas"

generate-metrics-p8: scalability_metrics
	@echo "📈 Generando métricas para P = 8..."
	mpirun -np 8 ./scalability_metrics 8
	@echo "✅ Métricas P=8 completadas"

generate-metrics-p16: scalability_metrics
	@echo "📈 Generando métricas para P = 16..."
	mpirun -np 16 ./scalability_metrics 16
	@echo "✅ Métricas P=16 completadas"

generate-metrics-p32: scalability_metrics
	@echo "📈 Generando métricas para P = 32..."
	mpirun -np 32 ./scalability_metrics 32
	@echo "✅ Métricas P=32 completadas"

# Generar todas las métricas secuencialmente
generate-all-metrics: scalability_metrics
	@echo "🚀 Iniciando generación completa de métricas de escalabilidad..."
	@echo "📋 Se generarán métricas para P = 1, 2, 4, 8, 16, 32"
	@echo "📄 Todos los resultados se guardarán en MetricsGraphics.txt"
	@echo ""
	@# Limpiar archivo anterior
	@rm -f MetricsGraphics.txt
	@echo "=== ANÁLISIS DE ESCALABILIDAD FUERTE - DFA PARALELO ===" > MetricsGraphics.txt
	@echo "Generado el: $$(date)" >> MetricsGraphics.txt
	@echo "Sistema: MPI con comunicación bloqueante y no bloqueante" >> MetricsGraphics.txt
	@echo "DFA: |Q| = 3, |Σ| = 2, acepta strings que contienen \"00\"" >> MetricsGraphics.txt
	@echo "" >> MetricsGraphics.txt
	@# Ejecutar cada configuración
	@$(MAKE) generate-metrics-p1
	@$(MAKE) generate-metrics-p2  
	@$(MAKE) generate-metrics-p4
	@$(MAKE) generate-metrics-p8
	@$(MAKE) generate-metrics-p16
	@$(MAKE) generate-metrics-p32
	@echo ""
	@echo "🎉 ¡Generación completa de métricas finalizada!"
	@echo "📊 Revisa el archivo MetricsGraphics.txt para ver todos los resultados"
	@echo "📈 Las métricas incluyen: Speedup, Eficiencia, Tiempos de cómputo y comunicación"

.PHONY: all clean original optimized metrics check-sources run-original run-optimized compare generate-metrics generate-metrics-np generate-metrics-p1 generate-metrics-p2 generate-metrics-p4 generate-metrics-p8 generate-all-metrics 