/** TAREA 1, DISEÑO Y ANÁLISIS DE ALGORITMOS
 * Valentina Alarcón
 * Naomi Cautivo
 * Máximo Flores
 * 
 * Método: Ciaccia Patella
*/

// De acuerdo a esto, diremos que cada nodo tendrá como capacidad B entradas en disco, es decir, que
// el tamaño de un nodo es a lo más B · sizeof(entry).

#include <cstddef>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm> //shuffle
#include <cstdlib>
#include <limits>
#include <cmath>
#include <utility>
#include <map>

#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

using namespace std;

int B; // capacidad máxima de un nodo
typedef struct entry entry;
typedef struct nodo nodo;

typedef struct { // una entrada de un nodo
    pair<double, double> p;
    double cr;
    nodo *a;
} entry;

size_t tamaño_max = B * sizeof(entry); // capacidad: B entradas en disco
int b_min = 0.5 * B; // capacidad mínima
int b_max = B; // capacidad máxima

// establecer semilla
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine e(seed);

typedef struct {
    vector<entry> entries;
    int b_min = 0.5 * B; // capacidad mínima
    int b_max = B; // capacidad máxima

    /** métodos a definir:
    1. insertar entrada nueva:
    - esta función permite agregar una nueva entrada al vector entries ssi no se ha logrado la capacidad máxima
    - no se preocupa de la capacidad mínima, pues entiende que cuando se utiliza este método:
    o estamos hablando de la raíz
    o se accedió a un caso donde el nodo ya ha sido inicializado correctamente (tiene más o exactamente B/2 nodos)
    2. eliminar una entrada del vector de entradas:
    - no necesita revisar la capacidad máxima, pero SÍ la mínima!
    3. bool: esValido
    - revisa si el nodo cumple con las condiciones de capacidad mínima y máxima
    */
    
    // insertar nueva entrada en entries
    void insertarEntry(const entry entrada) {
        if (entries.size() < b_max) { // si todavía no se alcanza la capacidad máxima, todo bien
            entries.push_back(entrada);
        } else {
            cerr << "Nodo lleno" << endl; // error
        }
    }

    // eliminar una entrada del punto IMPLEMENTACIÓN PROVISORIA
    void eliminarEntry() {
        if (entries.size() > b_min) { // no se puede si queda con menos de b/2
            entries.pop_back();
        } else {
            cerr << "Capacidad mínima excedida" << endl;
        }
    }

    bool esValido(){
        if (b_min < entries.size() && entries.size() < b_max){
            return true;
        }
        else{
            return false;
        }
    }
} nodo;

double distancia_cuadrado(const pair<double,double> punto1, const pair<double,double> punto2){
    // diferencia de coordenadas
    double dx = punto1.first - punto2.first;
    double dy = punto1.second - punto2.second;

    // Calcula el cuadrado de la distancia euclidiana
    double distanciaCuadrada = dx * dx + dy * dy;

    // Retorna la raíz cuadrada de la distancia cuadrada para obtener la distancia euclidiana
    return distanciaCuadrada;
}

/** la siguiente función encuentra el punto más cercano dentro de las claves de un mapa,
 * y luego se agrega al conjunto asociado a dicho punto
*/
map<pair<double, double>, set<pair<double, double>>> punto_mas_cercano(set<pair<double, double>> points, map<pair<double, double>, set<pair<double, double>>> mapa){
    for(pair<double,double> punto: points){
        // le asignamos su sample más cercano!

        // para esto podemos utilizar la ventaja del orden
        // dado que en el mapa las llaves están ordenadas, podemos usar búsqueda binaria(?)

        // o simplemente hacer una pasada lineal y ser felices xd (rip cuenteo)

        double mas_cercano = numeric_limits<double>::max(); // el más grande de los doubles
        pair<double,double> punto_mas_cercano = make_pair(0.0, 0.0);
        // aplicamos algoritmo
        for(const auto& par : mapa){
            pair<double,double> clave = par.first;
            double distancia = distancia_cuadrado(clave, punto);
            if(distancia < mas_cercano){
                mas_cercano = distancia;
                punto_mas_cercano = clave;
            }
        }
        // termina algoritmo
        mapa[punto_mas_cercano].insert(punto);
    };
    return mapa;
}


// Algoritmo CP:
// Input: Un set de puntos P
// 1. Si |P| ≤ B, se crea un árbol T , se insertan todos los puntos a T y se retorna T .
// 2. De manera aleatoria se eligen k = min(B, n/B) puntos de P, que los llamaremos samples pf1, . . . , pfk.
// Se insertan en un conjunto F de samples.
// 3. Se le asigna a cada punto en P su sample más cercano. Con eso se puede construir k conjuntos F1, . . . , Fk.
// 4. Etapa de redistribución: Si algún Fj es tal que |Fj| < b:
// 4.1 Quitamos pfj de F
// 4.2 Por cada p ∈ Fj, le buscamos el sample pfl más cercano de F y lo añadimos a su conjunto Fl.
// 5. Si |F| = 1, volver al paso 2.


nodo crear_MTree_CCP(const set<pair<double,double>> points){
    int n = points.size();
    if(n <= B){ // entran todos los puntos en un nodo
        nodo T; // se crea un árbol T (un nodo raíz con vector entries vacío)
        // se insertan todos los puntos a T
        for(pair<double,double> punto: points){
            // crear entrada, inicialmente con radio cobertor y a nulos.
            entry entrada;
            entrada.p = punto;
            entrada.cr = 0.0; // pues double
            entrada.a = nullptr; // pues puntero
            T.insertarEntry(entrada);
        }
        return T; // se retorna T
    }
    else{
        map<pair<double,double>, set<pair<double,double>>> samples;
        do{
            int k = min(B, n/B);

            // IDEA: los primeros k después del shuffle son los seleccionados.
            
            vector<int> indices(n); // un vector de tamaño n
            for(int i = 0; i < n; ++i){
                indices[i] = i; // agregamos los números del 0 al n-1
            }
            // hacemos shuffle
            shuffle(indices.begin(), indices.end(), e);

            // los k primeros son los seleccionados!
            set<int> iteradores;
            for(int i = 0; i < k; ++i){
                iteradores.insert(indices[i]); // recordar que los sets se ordenan de menor a mayor
            }
            // iterador
            auto it = iteradores.begin();
            auto punto_escogido = points.begin();
            int actual = 0;
            while(it!=iteradores.end()){ // iteradores es el que detiene el loop, pues se entiende que nunca tendrá un valor mayor a n-1.
                for(int i=actual; i<*it; i++){ // avanza hasta el índice al que apunta it. comienza en actual pues va secuencialmente
                    actual++;
                    punto_escogido++;
                }
                // se llegó al punto del iterador! se agrega a samples con un un conjunto vacío
                // samples.insert(*punto_escogido)
                samples[*punto_escogido] = set<pair<double,double>>();
                it++;
            }

            // 3. Se le asigna a cada punto en P su sample más cercano. Con eso se puede construir k conjuntos F1, . . . , Fk.
            samples = punto_mas_cercano(points, samples);
            // conjuntos armados!

            // 4. Etapa de redistribución: Si algún Fj es tal que |Fj| < b:
            // 4.1 Quitamos pfj de F
            // 4.2 Por cada p ∈ Fj, le buscamos el sample pfl más cercano de F y lo añadimos a su conjunto Fl.

            for(const auto& par : samples){
                pair<double,double> clave = par.first;
                set<pair<double,double>> conjunto = par.second;
                if(conjunto.size() < b_min){
                    samples.erase(clave); // quitamos pfj de F
                    // para cada elemento que estaba en su conjunto, buscaremos su otro sample más cercano
                    samples = punto_mas_cercano(conjunto, samples);
                }
            }

        } while (samples.size() == 1);

        //5. Si |F| = 1, volver al paso 2.

        // 6. Se realiza recursivamente el algoritmo CP en cada Fj, obteniendo el árbol Tj
        // si llegamos a esta parte, samples tiene más de un conjunto c:
        //nodo raiz; // el nodo al que asociaremos hijos y cosas
        map<pair<double,double>, nodo>  subarboles;
        //nodo raiz;
        for(const auto& par : samples){
            nodo sub_arbol = crear_MTree_CCP(par.second); // se llama recursivamente a crear_MTree_CCP con los conjuntos de cada uno
            // cada llamada devuelve un nodo (raíz de árbol). luego, tendremos varios árboles Tj.
            vector<entry> entradas = sub_arbol.entries;
            if(entradas.size() >= b_min){
                subarboles[par.first] = sub_arbol;
                //raiz.insertarEntry({par.first, 0.0, &sub_arbol});
            }
            else{
                // se quita la raiz, se elimina pfj de F y se trabaja con sus subárboles
                samples.erase(par.first);
                // para acceder a los subárboles, encontraremos los nodos a los que referencia cada entry
                // luego, estos nodos los agregagremos al set de subarboles
                // y finalmente, agregamos las entradas que tenía la raíz del subárbol a samples.
                for(auto entrada = entradas.begin(); entrada != entradas.end(); entrada++){
                    nodo sub_subarbol = *entrada->a;
                    subarboles[entrada.p] = sub_subarbol;
                    //samples[entrada.p] = se;
                    //raiz.insertarEntry({entrada.p, 0.0, &sub_subarbol});
                }
            }
        }

        // 8. Etapa de balanceamiento: Se define h como la altura mínima de los árboles Tj.
        // Se define T"'" inicialmente como un conjunto vacío.

        int h; //la altura mínima de los árboles Tj. IDEA: Cada nodo tiene su altura como propiedad??
        set<nodo> T2;

        // 9. Por cada Tj , si su altura es igual a h, se añade a T"′". Si no se cumple:
        // ** esto podría ser una fn recursiva, tal que la pueda llamar de nuevo!!
        for(const auto& par_T_j : subarboles){
            pair<double,double> clave = par_T_j.first;
            nodo T_j = par_T_j.second;
            if(T_j.altura== h){
                T2.insert(T_j);
            }
            else{
                // 9.1 Se borra el punto pertinente en F.
                subarboles.erase(clave);

                // 9.2 Se hace una búsqueda exhaustiva en Tj de todos los subárboles T"′" 1, . . . , T"′" p de altura igual a h.
                // Se insertan estos árboles a T"′"
                set<nodo> subtrees_h = busqueda_h(T_j); // funcion que busca y retorna un set con los subárboles de altura h que tiene T_j
                for (nodo subtree_h : subtrees_h) {
                    T2.insert(subtree_h);

                    // 9.3 Se insertan los puntos raíz de T"′"1, . . . , T"′" p: p"′" f1, . . . , p"′" fp en F
                    for (entry entrada_h : subtree_h.entries){
                        subarboles[entrada_h.p] = subtree_h;
                    }
                }
            }
        }
        
        set<pair<double,double>> keys; // las llaves de F (finalmente, F)

        // 10. Se define Tsup como el resultado de la llamada al algoritmo CP aplicado a F.
        for (const auto& pair : subarboles) {
            keys.insert(pair.first);
        }
        nodo Tsup = crear_MTree_CCP(keys);
        
        // 11. Se une cada Tj ∈ T"′" a su hoja en Tsup correspondiente al punto pfj ∈ F, obteniendo un nuevo árbol T .

        // 12. Se setean los radios cobertores resultantes para cada entrada en este árbol.

        // 13. Se retorna T .
        nodo T;
        return T;
    }
}

