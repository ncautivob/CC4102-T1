/** TAREA 1, DISEÑO Y ANÁLISIS DE ALGORITMOS
 * Valentina Alarcón
 * Naomi Cautivo
 * Máximo Flores
 * 
 * Método: Ciaccia Patella
*/

/** En primer lugar, se incluyen todas las librerías necesarias para la implementación */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <utility>
#include <vector>
#include <fstream>

using namespace std;
using namespace std::chrono;
using Point = pair<double, double>;
using Points = set<Point>;

int B_bytes = 4096;
struct Nodo;
#define N_ITER 100
#define query_radius 0.02

typedef struct Entry {  // una entrada de un nodo
  pair<double, double> p;
  double cr;
  Nodo* a;
} Entry;

int B = B_bytes / sizeof(Entry);  // (128)
int b_min = 0.5 * B;             // capacidad mínima (64)

/** T1-CCP */
/** Establecer semilla para selección aleatoria posterior */
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine e(seed);

/** A continuación, se describe la estructura Nodo */
typedef struct Nodo {
    vector<Entry> entries;
    int altura;
    
    /** Este método es utilizado para insertar una entrada en el vector 'entries' */
    void insertarEntry(const Entry entrada) {
        if (entries.size() < B) { // Si todavía no se alcanza la capacidad máxima, todo bien
            entries.push_back(entrada);
        } else { // Si se excede, error
            cerr << "Nodo lleno" << endl;
        }
    }
} Nodo;

/** Esta función elimina la memoria alocada para los nodos del árbol creado */
void deleteTree(Nodo *nodo) {
    if (nodo == nullptr) // Si el nodo es un puntero nulo, lo ignoramos
        return;
    vector<Entry> entradas = nodo->entries;
    if(entradas[0].a==nullptr){ // Si las entradas del nodo son punteros nulos, lo ignoramos
        delete nodo;
        return;
    }
    else{
        for(Entry entrada : entradas){ // Si no, detonamos la función recursivamente para el hijo
            deleteTree(entrada.a);
            entrada.a = nullptr;
        }
    }

    // Finalmente se elimina el nodo raíz
    delete nodo;
}

/** Esta función calcula la distancia euclideana al cuadrado entre dos puntos (p,q) */
double distancia_cuadrado(const Point& punto1, const Point& punto2){
    double dx = punto1.first - punto2.first;
    double dy = punto1.second - punto2.second;
    double distanciaCuadrada = dx * dx + dy * dy;

    return distanciaCuadrada;
}

/** La siguiente función encuentra, para cada punto de un set de puntos 'points',
 * su punto más cercano dentro de las claves de un mapa,
 * y luego se agrega al conjunto (el valor de su clave en el mapa) asociado a dicho punto.
*/
void punto_mas_cercano(const set<Point>& points, map<Point, set<Point>>& mapa){
    for(const Point& punto: points){
        double mas_cercano = numeric_limits<double>::max(); // El más grande de los doubles
        Point punto_mas_cercano = make_pair(0.0, 0.0);

        // Aplicamos algoritmo
        for(const auto& par : mapa){
            const Point clave = par.first;
            double distancia = distancia_cuadrado(clave, punto); // Se calcula la distancia entre la clave y el punto
            if(distancia < mas_cercano){
                mas_cercano = distancia;
                punto_mas_cercano = clave;
            }
        }
        mapa[punto_mas_cercano].insert(punto); // Se asigna en el mapa
    };
}

/** Función que, dada una raíz de un árbol, encuentra los subárboles de altura h que estén en él.
 * Tras encontrar cada uno, se ingresa al mapa 'arboles' el punto asociado a su raíz (como clave)
 * Y el nodo raíz del subárbol de altura h (como su valor)
 * 
 * Nota: dado que h es la altura mínima de los árboles, entonces nunca ocurrirá que no se encuentre ninguno
*/
void busqueda_h(Nodo *nodo, const int h, map<Point, Nodo*>& arboles) {
    if(nodo==nullptr) return;
    vector<Entry> entradas = nodo->entries;
    if(entradas[0].a == nullptr){ // Nodo es hoja
        return;
    }
    for(Entry entrada : entradas){
        Point clave = entrada.p;
        Nodo* hijo = entrada.a; // No es nullptr
        if(hijo->altura == h){ // Encontramos un subárbol de altura h
            arboles[clave] = hijo;
        }
        else if(hijo->altura < h) continue;
        else{ // Altura mayor a h
            busqueda_h(hijo, h, arboles);
        }
    }
    for(auto& entrada : entradas) {
        entrada.a = nullptr;
    }
    delete nodo;
}

/** Función que seta el radio cobertor de una entrada */
double setear_radio_cobertor(Entry& entry){ // No seteará las hojas porque originalmente ya son 0.0

    double max_radio = 0;
    Point punto = entry.p;
    Nodo *hijo = entry.a;
    if(hijo==nullptr){
        return 0.0;
    }
    // Se maximizará (radio_cobertor_hijo + distancia_al_hijo)
    for (Entry& entrada : hijo->entries) { // Revisamos las entradas del hijo
        double radio_cobertor_hijo = entrada.cr;
        double distancia = sqrt(distancia_cuadrado(punto, entrada.p));
        if ((radio_cobertor_hijo + distancia)>max_radio){
            max_radio = radio_cobertor_hijo + distancia;
        }
    }
    entry.cr = max_radio;


    return entry.cr;
}

/** Función que conecta los subarboles Tj a las hojas de Tsup,
 * Arreglando las alturas en el proceso,
 * Y actualizando los radios cobertores de Tsup.
 */
void conectar_arboles(Nodo* nodo, map<Point,Nodo*>& subarboles){
    if(nodo==nullptr){ // Si el nodo es un puntero nulo, se retorna.
        return;
    }
    vector<Entry>& entries = nodo->entries;
    int altura_max = 0;
    int entries_sz = (int)entries.size();
    if(entries_sz == 0) return;
    if (entries[0].a==nullptr) { // Encontramos una hoja
        nodo->altura = 2;
        for(int j = 0; j<entries_sz; j++){
            if(j>=entries.size()) break;
            Point punto_hoja = entries[j].p;
            entries[j].a = subarboles[punto_hoja];
            setear_radio_cobertor(entries[j]);
        }
    }
    else { // Si el nodo no es una hoja, recorre recursivamente sus entradas
        for (auto& entry : entries) {
            conectar_arboles(entry.a, subarboles);
            int altura = ((entry.a)->altura)+1;
            if (altura > altura_max){
                altura_max = altura;
            }
            setear_radio_cobertor(entry);
        }
        nodo->altura = altura_max; // Se setea la altura para el nodo
    }
}


// Algoritmo CP:
// Input: Un set de puntos P

Nodo *crear_MTree_CCP(const set<Point> points){
    int n = points.size();
    //cout << "n es" << n << endl;
    if(n <= B){ // Entran todos los puntos en un nodo!
        Nodo *T = new Nodo; // Se crea un árbol T (un nodo raíz con vector entries vacío)
        for (Point punto: points){
            // Inicialmente con radio cobertor y a nulos.
            Entry entrada;
            entrada.p = punto;
            entrada.cr = 0.0;
            entrada.a = nullptr;
            T->insertarEntry(entrada);
        }
        T->altura = 1;
        return T; // se retorna T
    }
    else{
        int k = int(ceil(min(double(B), double(n)/B))); // Debe ser el techo, tal que no llegue a 1!
        //cout << "k es" << k << endl;
        map<Point, set<Point>> samples;
        do{
            samples = map<Point, set<Point>>();
            vector<int> indices(n);
            for(int i = 0; i < n; ++i){
                indices[i] = i;
            }
            // Hacemos shuffle
            shuffle(indices.begin(), indices.end(), e);

            // Los k primeros son los seleccionados!
            set<int> iteradores;
            for(int i = 0; i < k; ++i){
                iteradores.insert(indices[i]);
            }

            auto punto_escogido = points.begin();
            int actual = 0;
            for(int iterador : iteradores){
                for(int i=actual; i<iterador; i++){
                    actual++;
                    punto_escogido++;
                }
                samples[*punto_escogido] = set<Point>();
            }

            // Armar los k conjuntos:
            punto_mas_cercano(points, samples);

            // Etapa de redistribución!
            auto it = samples.begin();
            while (it != samples.end()) {
                set<Point>& conjunto = it->second;
                //cout << "el conjunto mide " << conjunto.size() << endl;
                if (conjunto.size() < b_min) {
                    set<pair<double, double>> conjunto_copia = conjunto; // Hacer una copia
                    it = samples.erase(it); // Quitamos pfj de F
                    // Se redistribuye
                    punto_mas_cercano(conjunto_copia, samples);
                } else {
                    ++it; // Se avanza al siguiente elemento
                }
            }
        } while (samples.size() <= 1);
        // Lo anterior se vuelve a realizar si 'samples' finaliza la iteración con tamaño 1

        int h = numeric_limits<int>::max(); // La altura mínima de los árboles T_j.

        map<Point, Nodo*>  subarboles; // Mapa que linkea los puntos raíces con los nodos de sus subárboles asociados

        // Se realizarán las llamadas recursivas con los conjuntos encontrados previamente
        // Y de paso, se determinará la altura mínima de los subárboles (h).
        for(const auto& par : samples){
            Nodo *sub_arbol = crear_MTree_CCP(par.second);
            vector<Entry>& entradas = sub_arbol->entries;

            if(entradas.size() >= b_min){ // Si el nodo tiene más de b_min entradas, todo ok
                subarboles[par.first] = sub_arbol;
                if(sub_arbol->altura < h){
                    h = sub_arbol->altura;
                }
            }
            else{
                /** Tomaremos los subárboles del nodo */
                for(auto entrada = entradas.begin(); entrada != entradas.end(); entrada++){
                    Nodo *sub_subarbol = entrada->a;
                    subarboles[entrada->p] = sub_subarbol;
                    //cout << "la altura de este subarbol es" << sub_subarbol->altura << endl;
                    if(sub_subarbol->altura < h){
                        h = sub_subarbol->altura;
                    }
                    entrada->a = nullptr;
                }
                // for(auto& entrada : entradas) {
                //     entrada.a = nullptr;
                // }
                // // Eliminar el nodo raíz
                delete sub_arbol;
            }
        }
        //cout << "la altura minima encontrada fue " << h << endl;
        
        set<Point> claves_a_eliminar; // Set de las claves que se eliminarán del mapa 'subarboles' posteriormente
        set<Point> keys; // Las llaves del mapa 'subarboles' (el conjunto F del enunciado)
        map<Point, Nodo*> arboles_alt_h; // El mapa que entregará las claves y subárboles a agregar a 'subarboles' posteriormente

        /** A continuación, se balanceará el árbol de acuerdo a la altura mínima hallada */
        for(const auto &par_T_j : subarboles){
            Point clave = par_T_j.first;
            Nodo *T_j = par_T_j.second;
            if(T_j->altura== h){ // Si la altura del subárbol es h, todo sigue ok.
                keys.insert(clave); // Se inserta la clave en el set 'keys'.
            }
            else{
                // 9.1 Se borra el punto pertinente en F.
                claves_a_eliminar.insert(clave);

                // 9.2 Se hace una búsqueda exhaustiva en Tj de todos los subárboles de altura igual a h.
                // Se insertan estos árboles a T′
                busqueda_h(T_j, h, arboles_alt_h); // Se insertan los puntos raíces linkeados con sus subárboles en el mapa 'arboles_alt_h'
                //cout << "termina busqueda" << endl;
            }
        }

        /** Se eliminan las claves acordes a lo anterior */
        for(Point clave : claves_a_eliminar){
            subarboles.erase(clave);
        }

        // 9.3 Se insertan los puntos raíz de T′1, . . . , T′ p: p′ f1, . . . , p′ fp en F
        for (const auto &par_T_j : arboles_alt_h) {
            Point clave = par_T_j.first;
            Nodo *T_j = par_T_j.second;
            keys.insert(clave); // Y la llave acorde a 'keys'
            subarboles[clave] = T_j;
        }

        /** Se crea el árbol Tsup */
        Nodo *Tsup = crear_MTree_CCP(keys);

        /** Se conectan las hojas de Tsup con los subárboles Tj*/
        conectar_arboles(Tsup, subarboles);

        return Tsup;
    }
}

set<Point> crear_set(int n){
    set<Point> ccp_set;
    for(int j=0; j<n; j++){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dis(0.0, 1.0);

        // Generate random double values
        double first = dis(gen);
        double second = dis(gen);

        // Create and return the pair
        ccp_set.insert(make_pair(first, second));
    }
    return ccp_set;
}


/** search method */
void search(Nodo* MTree, Point q, double r, vector<int>& v, int index) {
  if (MTree == nullptr) return; // Si no hay MTree, rip

  // El objetivo es encontrar todos los puntos de T que residen dentro de la bola (q,r)

    vector<Entry> entradas = MTree->entries;
  double r_square = pow(r,2);  // dado que la función de distancia entrega la distancia al cuadrado

  if (entradas[0].a == nullptr) {  // nodo es una hoja
    v[index]++; // Se aumenta en uno los accesos
    // se verifica para cada entrada si p cumple con dist(p, q) ≤ r
        for (Entry entrada : entradas) {
        Point punto = entrada.p;
        if (distancia_cuadrado(punto, q) <= r_square) {
            //resp.insert(punto);
        }
    }
  } else {  // nodo es interno
    v[index]++;
    for (Entry entrada : entradas) {
        Point punto = entrada.p;
        double cr = entrada.cr;
        double cr_square = pow(cr, 2);
        if (distancia_cuadrado(punto, q) <= r_square + cr_square) {
            // se busca en su hijo a posibles respuestas
            search(entrada.a, q, r, v, index);
        }
        // si no, se descarta
        }
    }
}

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

    for(int i = 10; i < 26; i++){
        cout << "probando con i=" << i << endl;
        Points testSet = parsePoints(i);
        Points testQueries = parsePoints(N_ITER);

        vector<Point> queryPoints;
        vector<pair<int, double>> response;
        vector<int> accesses(N_ITER);

        for (Point point : testQueries) {
            queryPoints.push_back(point);
        }

        auto start = high_resolution_clock::now();
        Nodo* root = crear_MTree_CCP(testSet);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<nanoseconds>(stop - start);

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

        double standardDeviation = sqrt(varianza);
        cout << "la desviacion estandar es " << standardDeviation << endl;
        cout << "y la varianza es" << varianza << endl;
        deleteTree(root);
    }

    /** Volvemos la salida estándar al original */
    cout.rdbuf(originalStdout);

    outputFile.close();
    return 0;
}