#!/usr/bin/env python3
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

# Rutas a las carpetas de entrada
test_set_dir = './test_set'
epitopes_dir = './epitopes'

# Comando base para ejecutar
base_command = "python main_epitopes.py -t {} -e {} -pp ../model/retrained/TCRen_TCR_p_mydata_train.csv -g ../structures_annotation/general.txt -p ../pdb_files -s {}_mydata -metric TCRdist"

def run_command(tcr_file, epitope_path, output_file):
    """
    Ejecuta el comando para procesar el archivo TCR y el archivo de epítopos.
    """
    command = base_command.format(tcr_file, epitope_path, output_file)
    print(f"Ejecutando: {command}")
    subprocess.run(command, shell=True)

def main():
    # Lista para almacenar las tareas
    tasks = []

    # Iterar sobre los archivos en la carpeta test_set
    for tcr_file in os.listdir(test_set_dir):
        if tcr_file.endswith('_df.csv'):
            # Obtener el nombre base del archivo TCR
            base_name = tcr_file.replace('_df.csv', '')
            epitope_file = f"{base_name}_epitopes.txt"  # Nombre del archivo de epítopos esperado

            # Construir la ruta completa del archivo de epítopos
            epitope_path = os.path.join(epitopes_dir, epitope_file)

            # Verificar si el archivo de epítopos existe
            if os.path.isfile(epitope_path):
                output_file = f"output_{base_name}.csv"  # Nombre del archivo de salida
                tasks.append((os.path.join(test_set_dir, tcr_file), epitope_path, output_file))
            else:
                print(f"Advertencia: No se encontró el archivo de epítopos para {tcr_file}. Se esperaba {epitope_path}.")

    # Ejecutar las tareas en paralelo
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(run_command, tcr_file, epitope_path, output_file): base_name for tcr_file, epitope_path, output_file in tasks}
        
        for future in as_completed(futures):
            base_name = futures[future]
            try:
                future.result()  # Obtener el resultado, esto lanzará cualquier excepción
                print(f"Tarea completada para: {base_name}")
            except Exception as e:
                print(f"Error al procesar {base_name}: {e}")

    print("Finalizado.")

if __name__ == "__main__":
    main()
