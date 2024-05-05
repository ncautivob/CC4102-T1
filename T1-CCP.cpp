/** TAREA 1, DISEÑO Y ANÁLISIS DE ALGORITMOS
 * Valentina Alarcón
 * Naomi Cautivo
 * Máximo Flores
 * 
 * Método: Ciaccia Patella
*/

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

int B_solo = 4096; // esto dividido en sizeof(entry) equivale al B real
struct Nodo;

typedef struct Entry { // una entrada de un nodo
    pair<double, double> p;
    double cr;
    Nodo *a;
} Entry;

int B = B_solo / sizeof(Entry); // duda: INT O SIZE_T?
int b_min = 0.5 * B; // capacidad mínima

// establecer semilla
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine e(seed);

typedef struct Nodo {
    vector<Entry> entries;
    int altura;
    int b_min = 0.5 * B; // capacidad mínima
    int b_max = B; // capacidad máxima
    
    // insertar nueva entrada en entries
    void insertarEntry(const Entry entrada) {
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
} Nodo;

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
map<pair<double, double>, set<pair<double, double>>> punto_mas_cercano(const set<pair<double, double>> points, map<pair<double, double>, set<pair<double, double>>> mapa){
    for(pair<double,double> punto: points){
        // le asignamos su sample más cercano!
        // problema con binsearch: 2d

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

/** función que, dada una raíz de un árbol, encuentra los subárboles de altura h que estén en él 
 * nota: dado que h es la altura mínima de los árboles, entonces nunca ocurrirá que no se encuentre ninguno
*/
set<Nodo*> busqueda_h(Nodo *nodo, const int h, set<Nodo*>& arboles) {
    if(nodo->altura == h){
        arboles.insert(nodo);
    }
    else if (nodo->altura > h) { // si su altura es menor a h, sabemos que no tiene hijos de altura h
        for (const auto& hijo : nodo->entries) {
            // cada entry apunta a un nodo
            busqueda_h(hijo.a, h, arboles);
        }
    }
    return arboles;
}

int setear_radio_cobertor(Entry& entry){ // no seteará las hojas porque originalmente ya son 0.0

    int max_radio = 0;
    Nodo *hijo = entry.a;
    if(hijo == nullptr){
        return 0;
    }
    for (Entry entrada : hijo->entries){
        max_radio = max(max_radio, setear_radio_cobertor(entrada) + (int)distancia_cuadrado(entry.p, entrada.p));
    }

    entry.cr = max_radio;

    return max_radio;
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


Nodo *crear_MTree_CCP(const set<pair<double,double>> points){
    int n = points.size();
    if(n <= B){ // entran todos los puntos en un nodo
        Nodo *T = new Nodo; // se crea un árbol T (un nodo raíz con vector entries vacío)
        // se insertan todos los puntos a T
        for (pair<double,double> punto: points){
            // crear entrada, inicialmente con radio cobertor y a nulos.
            Entry entrada;
            entrada.p = punto;
            entrada.cr = 0.0; // pues double
            entrada.a = nullptr; // pues puntero
            (*T).insertarEntry(entrada);
        }
        T->altura = 0;
        return T; // se retorna T
    }
    else{
        map<pair<double,double>, set<pair<double,double>>> samples;
        int k = min(B, n/B);
        do{
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

        } while (samples.size() == 1); //5. Si |F| = 1, volver al paso 2.

        int h = numeric_limits<int>::max(); // la altura mínima de los árboles Tj.

        // 6. Se realiza recursivamente el algoritmo CP en cada Fj, obteniendo el árbol Tj
        // si llegamos a esta parte, samples tiene más de un conjunto c:

        map<pair<double,double>, Nodo*>  subarboles;
        for(const auto& par : samples){
            Nodo *sub_arbol = crear_MTree_CCP(par.second); // se llama recursivamente a crear_MTree_CCP con los conjuntos de cada uno
            // cada llamada devuelve un nodo (raíz de árbol). luego, tendremos varios árboles Tj.
            vector<Entry> entradas = sub_arbol->entries;
            if(entradas.size() >= b_min){
                subarboles[par.first] = sub_arbol;
                if(sub_arbol->altura < h){
                    h = sub_arbol->altura;
                }
                //raiz.insertarEntry({par.first, 0.0, &sub_arbol});
            }
            else{
                // se quita la raiz, se elimina pfj de F y se trabaja con sus subárboles
                //samples.erase(par.first);
                // para acceder a los subárboles, encontraremos los nodos a los que referencia cada entry
                // luego, estos nodos los agregagremos al set de subarboles
                // y finalmente, agregamos las entradas que tenía la raíz del subárbol a samples.
                for(auto entrada = entradas.begin(); entrada != entradas.end(); entrada++){
                    Nodo *sub_subarbol = entrada->a;
                    subarboles[entrada->p] = sub_subarbol;
                    if(sub_subarbol->altura < h){
                        h = sub_arbol->altura;
                    }
                    //samples[entrada.p] = se;
                    //raiz.insertarEntry({entrada.p, 0.0, &sub_subarbol});
                }
            }
        }

        // 8. Etapa de balanceamiento: Se define h como la altura mínima de los árboles Tj.
        // Se define T"'" inicialmente como un conjunto vacío.

        set<Nodo*> T2 = set<Nodo*>();

        // 9. Por cada Tj , si su altura es igual a h, se añade a T"′". Si no se cumple:
        // ** esto podría ser una fn recursiva, tal que la pueda llamar de nuevo!!
        for(const auto& par_T_j : subarboles){
            pair<double,double> clave = par_T_j.first;
            Nodo *T_j = par_T_j.second;
            if(T_j->altura== h){
                T2.insert(T_j);
            }
            else{
                // 9.1 Se borra el punto pertinente en F.
                subarboles.erase(clave);

                // 9.2 Se hace una búsqueda exhaustiva en Tj de todos los subárboles T"′" 1, . . . , T"′" p de altura igual a h.
                // Se insertan estos árboles a T"′"
                set<Nodo*> empty_set;
                set<Nodo *> subtrees_h = busqueda_h(T_j, h, empty_set); // funcion que busca y retorna un set con los subárboles de altura h que tiene T_j
                for (Nodo *subtree_h : subtrees_h) {
                    T2.insert(subtree_h);

                    // 9.3 Se insertan los puntos raíz de T"′"1, . . . , T"′" p: p"′" f1, . . . , p"′" fp en F
                    for (Entry entrada_h : subtree_h->entries){
                        subarboles[entrada_h.p] = entrada_h.a;
                    }
                }
            }
        }
        
        //set<pair<double,double>> keys; // las llaves de F (finalmente, F)

        // 10. Se define Tsup como el resultado de la llamada al algoritmo CP aplicado a F.

        //Nodo *Tsup = crear_MTree_CCP(keys);
        
        // 11. Se une cada Tj ∈ T"′" a su hoja en Tsup correspondiente al punto pfj ∈ F, obteniendo un nuevo árbol T.
        //for(Nodo *T_j : T2) {
            // unirlo a su hoja en Tsup!!
            // cómo ...
        //}

        // 12. Se setean los radios cobertores resultantes para cada entrada en este árbol.

        // 13. Se retorna T .
        Nodo *Tsup;
        int altura_max = 0;
        for(const auto& pair : subarboles){ // uno con los hijos
            (*Tsup).insertarEntry({pair.first, 0.0, pair.second});
            int altura = pair.second->altura;
            if(altura > altura_max){
                altura_max = altura;
            }
        }
        // setear altura de acuerdo a la altura mayor de sus nodos hijos
        Tsup->altura = altura_max + 1;
        for(Entry entrada: Tsup->entries){
            setear_radio_cobertor(entrada);
        }
        return Tsup;
    }
}

set<pair<double,double>> crear_set(int n){
    set<pair<double,double>> ccp_set;
    for(int j=0; j<n; j++){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dis(-100.0, 100.0);

        // Generate random double values
        double first = dis(gen);
        double second = dis(gen);

        // Create and return the pair
        ccp_set.insert(make_pair(first, second));
    }
    return ccp_set;
}

void imprimirArbol(Nodo *arbol)
{
    if (arbol != NULL)
    {

        cout << "Nodo: " << arbol->entries.size() << "\n";
        for (int i = 0; i < arbol->entries.size(); i++) {
            cout << "Entrada: " << arbol->entries[i].p.first << "" << arbol->entries[i].p.second << " con radio igual a " << arbol->entries[i].cr << " \n";
        }

        for (int i = 0; i < arbol->entries.size(); i++) {
            imprimirArbol(arbol->entries[i].a);
        }
    }
    else {
        cout << "Llego a externo \n";
    }
}

int main() {
    set<pair<double, double>> random_pairs = crear_set(5);
    for(pair<double,double> random_pair : random_pairs){
        //cout << "Random pairs: (" << random_pair.first << ", " << random_pair.second << ")" << endl;
    }
    cout << "hola" << endl;
    cout << "B es " << B << endl;
    int n = random_pairs.size();
    if(n <= B){ // entran todos los puntos en un nodo
        Nodo *T = new Nodo; // se crea un árbol T (un nodo raíz con vector entries vacío)
        T->entries = vector<Entry>();
        T->altura = 0;
        cout << T->entries.size() << endl;
        // se insertan todos los puntos a T
        for (pair<double,double> punto: random_pairs){
            // crear entrada, inicialmente con radio cobertor y a nulos.
            Entry entrada;
            entrada.p = punto;
            entrada.cr = 0.0; // pues double
            entrada.a = nullptr; // pues puntero
            (*T).insertarEntry(entrada);
        }
        for (Entry entry: T->entries){
            // crear entrada, inicialmente con radio cobertor y a nulos.
            cout << entry.p.first << " " << entry.p.second << endl;
            cout << entry.cr << endl;
            cout << entry.a << endl;
        }
        delete T;
    }
    Nodo *arbol = crear_MTree_CCP(random_pairs);
    for(Entry entry : arbol->entries){
        cout << "hola" << endl;
        cout << "entry: " << entry.p.first << " " << entry.p.second << endl;
    }
    delete arbol;
    

    return 0;
}