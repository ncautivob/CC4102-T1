/** TAREA 1, DISEÑO Y ANÁLISIS DE ALGORITMOS
 * Valentina Alarcón
 * Naomi Cautivo
 * Máximo Flores
 * 
 * Main
*/

#include <iostream>
#include <fstream>
#include "sexton-swinbank.h"
#include "ciaccia-patella.h"
#include "search.h"

using namespace std;

/** En este archivo se generará la experimentación.
 * Se testearán ambos métodos de construcción creados con sets de puntos de tamaño creciente,
 * realizand una serie de consultas de búsqueda sobre los árboles que construyen.
 * 
 * Se va a tener un set de puntos 'P' y un set de consultas 'Q'
 * Ambos tienen como coordenadas valores aleatorios de tipo 'double'
 * Los cuales están uniformemente distribuidos en el rango [0,1]
 * 
 * Algunas cosas a definir son:
 * 
 * El tamaño de un bloque de disco será 4096 Bytes.
 * Luego, el 'B' que se utiliza en el código será 4096/sizeof(Entry))
 * y el 'b' será 2048/sizeof(Entry).
 * 
 * Por otro lado, el radio r de las consultas será de 0.02.
 * 
 * La idea será evaluar el resultado para cada n ∈ {2^10,2^11,...,2^25}, siendo este valor la cantidad de puntos en P.
 * Ocupe el mismo conjunto P y Q para evaluar los dos métodos de construcción del M-tree.
 * Se realizarán 100 consultas por cada valor de n, y reportará el intervalo de confianza de sus resultados.
 * 
*/

Points parsePoints(int n) {
    string line;

    string fileName =
        n == N_ITER ? "queries.txt" : "pow2" + to_string(n) + ".txt";
    ifstream File("./input/" + fileName);

    Points result;
    if (File.is_open()) {
        while (getline(File, line)) {
        // Resets the index. The format is (x, y).
            int i = 0;
            string x, y;

            while (line[i] != '(') i++;
            i++;
            while (line[i] != ',') {
                x += line[i];
                i++;
            }

            i += 2;
            while (line[i] != ')') {
                y += line[i];
                i++;
            }

            result.insert({stod(x), stod(y)});
        }
        File.close();
    }

    return result;
}

int main() {
    ofstream outputFile("resultados.txt");

    /** Queremos redireccionar las salidas de stdout al archivo txt */
    streambuf* originalStdout = cout.rdbuf(outputFile.rdbuf());

    /** Primero se testará el método de construcción Ciacia Patella */

    cout << "Construyendo M-Tree con método Ciacia Patella" << endl;

    for(int i = 10; i < 26; i++){ // Se testea desde n=2^10 hasta n=2^25
        cout << "Probando con i = " << i << endl;
        Points testSet = parsePoints(i);
        Points testQueries = parsePoints(N_ITER);

        vector<Point> queryPoints;
        vector<pair<int, double>> response;
        vector<int> accesses(N_ITER);

        for (Point point : testQueries) {
            queryPoints.push_back(point);
        }

        auto start = high_resolution_clock::now();
        /** Nota: cada método tiene su propia definición de la estructura 'Nodo'
         * Pues Ciaccia-Patella requiere de una altura
         * Y Sexton-Swinbank no.
         * Luego, el método 'search' es el mismo, solo que 'search_CCP' recibe 'Nodo'
         * y 'search' recibe 'Node'
        */
        Nodo* root = crear_MTree_CCP(testSet);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);
        cout << "Construcción M-Tree (CCP): " << duration.count() << " ns\n";

        for (int i = 0; i < N_ITER; i++) {
            Point q = queryPoints[i];

            auto start = high_resolution_clock::now();
            search_CCP(root, q, query_radius, accesses, i);
            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<nanoseconds>(stop - start);
            response.push_back({accesses[i], duration.count()});
        }
        double suma = 0.0;
        for (int i = 0; i < N_ITER; i++) {
            cout << i + 1 << " (accesos, tiempo): (" << response[i].first << ", "
                << response[i].second << " ns)\n";
            suma += response[i].first;
        }
        cout << "La media de accesos es " <<  suma/100 << endl;
        double sum2 = 0.0;

        for (int i = 0; i < N_ITER; i++) {
            sum2 += (response[i].first - (suma/100)) * (response[i].first - (suma/100));
        }
        double varianza = (double)sum2 / 100;
        cout << "La varianza es" << varianza << endl;
        double dev_estandar = sqrt(varianza);
        cout << "La desviacion estándar es " << dev_estandar << endl;
        
        deleteTree(root); // Se elimina el árbol y se libera la memoria alocada.
    }

    /** Ahora se testeará Sexton-Swinbank */

    for(int i = 10; i < 16; i++){ // Se testea desde n=2^10 hasta n=2^15
        cout << "Probando con i = " << i << endl;
        Points testSet = parsePoints(i);
        Points testQueries = parsePoints(N_ITER);

        vector<Point> queryPoints;
        vector<pair<int, double>> response;
        vector<int> accesses(N_ITER);

        for (Point point : testQueries) {
            queryPoints.push_back(point);
        }

        auto start = high_resolution_clock::now();
        Node* root = SSAlgorithm(testSet);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);
        cout << "Construcción M-Tree (CCP): " << duration.count() << " ns\n";

        for (int i = 0; i < N_ITER; i++) {
            Point q = queryPoints[i];

            auto start = high_resolution_clock::now();
            search(root, q, query_radius, accesses, i);
            auto stop = high_resolution_clock::now();

            auto duration = duration_cast<nanoseconds>(stop - start);
            response.push_back({accesses[i], duration.count()});
        }
        double suma = 0.0;
        for (int i = 0; i < N_ITER; i++) {
            cout << i + 1 << " (accesos, tiempo): (" << response[i].first << ", "
                << response[i].second << " ns)\n";
            suma += response[i].first;
        }
        cout << "La media de accesos es " <<  suma/100 << endl;
        double sum2 = 0.0;

        for (int i = 0; i < N_ITER; i++) {
            sum2 += (response[i].first - (suma/100)) * (response[i].first - (suma/100));
        }
        double varianza = (double)sum2 / 100;
        cout << "La varianza es" << varianza << endl;
        double dev_estandar = sqrt(varianza);
        cout << "La desviacion estándar es " << dev_estandar << endl;
        
        deleteTree(root); // Se elimina el árbol y se libera la memoria alocada.
    }



    /** Volvemos la salida estándar al original */
    cout.rdbuf(originalStdout);

    outputFile.close();

    return 0;
}