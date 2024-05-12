# 📚 Tarea 1: CC4102 - Diseño y Análisis de Algoritmos

Construcción y búsqueda en M-Tree usando métodos de Ciaccia-Patella (CP) y Sexton-Swinbank (SS).

## 👤 Integrantes (Sección 2): 

- Valentina Alarcón Yañez (valentina.alarcon.y at ug.uchile.cl).
- Naomi Cautivo Bahamóndez (naomi.cautivo at ug.uchile.cl).
- Máximo Flores Valenzuela (mflores at dcc.uchile.cl).

## 👤 Profesores:

- Benjamín Bustos (bebustos at dcc.uchile.cl).
- Gonzalo Navarro Badino (gnavarro at dcc.uchile.cl).

## ✍️ Entregables

En el siguiente repositorio se entregan los siguientes archivos. Todas las implementaciones fueron realizadas en C++:

- 📁 $\texttt{/input/}$ - Contiene las entradas compartidas en formato $\texttt{.txt}$ desde $2^{10}$ a $2^{25}$ datos para la ejecución de CP y SS, además del archivo $\texttt{queries.txt}$ de los puntos de búsqueda. Estos sets fueron generados con la función 'crearSet' que se encuentra en ambos métodos, la cual utiliza una distribución uniforme que selecciona puntos aleatoriamente entre 0.0 y 1.0.
- 📄 $\texttt{sexton-swinbank.cpp}$ - Método de construcción Sexton-Swinbank con su función 'main'.
- 📄 $\texttt{ciaccia-patella.cpp}$ -  Método de construcción Ciaccia-Patella con su función 'main'.
- 📄 $\texttt{informe.pdf}$ - Informe de la tarea en formato PDF.

## 💻 Ejecución

Para realizar la etapa de experimentación, se acordaron ciertos parámetros:

 * El tamaño de un bloque de disco será 4096 Bytes.
 * Luego, el 'B' que se utiliza en el código será 4096/sizeof(Entry)
 * y el 'b' será 2048/sizeof(Entry).
 * El radio r de las consultas será de 0.02.
 * Luego, se utilizarán los archivos de la carpeta input para obtener el set de puntos P y el set de queries Q.

Luego, ambos 'main' ejecutan el mismo proceso:
- Iteran sobre una potencia de 2, encuentran el archivo asociado a esa potencia en la carpeta 'input'
- Detonan la creación del M-Tree con el método que corresponde, midiendo su tiempo de ejecución
- Luego, detona las 100 búsquedas, contando la cantidad de accesos y su tiempo asociado.
- También se calcula la media de los accesos, la desviación estándar y la varianza asociada.
- Finalmente, se libera la memoria asociada al árbol creado.

Para ejecutar la tarea, se debe ejecutar el siguiente comando en 'bash':

g++ ciaccia-patella.cpp -o ciaccia-patella && ./ciaccia-patella && g++ sexton-swinbank.cpp -o sexton-swinbank && ./sexton-swinbank

El cual detona el main de ambos métodos. Para CP, los resultados (la salida estándar de C++) se visualizarán en el archivo 'resultados-cp.txt', y para SS se verán en 'resultados-ss.txt'.

El archivo debería ser considerado como seguro si existe antivirus. Si existe algún problema con el acceso, puede que desactivar temporalmente el activirus sea una posible solución.

Notar además que los resultados pueden variar de acuerdo a lo visto en el informe de acuerdo al sistema que se utilice para correr la tarea.