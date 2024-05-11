/** TAREA 1, DISEÑO Y ANÁLISIS DE ALGORITMOS
 * Valentina Alarcón
 * Naomi Cautivo
 * Máximo Flores
 * 
 * Método: Ciaccia Patella
*/

/** En primer lugar, se incluyen todas las librerías necesarias para la implementación */
#include <cstddef>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <utility>
#include <map>
#include <random>
#include <chrono>

using namespace std; // Se setea esto para no tener que escribir 'std::' repetitivamente.
using Point = pair<double, double>; // Esto es para simplificar la notación.
struct Nodo; // Y esto es para poder utilizarlo en las estructuras a definir.


int B_bytes = 4096; // Esto dividido en sizeof(entry) equivale al máximo número de entradas de un nodo.

typedef struct Entry {
    Point p;
    double cr;
    Nodo *a;
} Entry;

int B = B_bytes / sizeof(Entry); // (128)
int b_min = 0.5 * B; // capacidad mínima (64)

// establecer semilla
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine e(seed);

typedef struct Nodo {
    vector<Entry> entries;
    int altura;
    
    // insertar nueva entrada en entries
    void insertarEntry(const Entry entrada) {
        if (entries.size() < B) { // si todavía no se alcanza la capacidad máxima, todo bien
            entries.push_back(entrada);
        } else {
            cerr << "Nodo lleno" << endl; // error
        }
    }
} Nodo;

void deleteTree(Nodo *nodo) {
    if (nodo == nullptr)
        return;
    vector<Entry> entradas = nodo->entries;
    if(entradas[0].a==nullptr){
        return;
    }
    else{
        for(Entry entrada : entradas){
            deleteTree(entrada.a);
        }
    }

    // Delete the current node
    delete nodo;
}

double distancia_cuadrado(const Point& punto1, const Point& punto2){
    double dx = punto1.first - punto2.first;
    double dy = punto1.second - punto2.second;
    double distanciaCuadrada = dx * dx + dy * dy;
    return distanciaCuadrada;
}

/** la siguiente función encuentra el punto más cercano dentro de las claves de un mapa,
 * y luego se agrega al conjunto asociado a dicho punto
*/
void punto_mas_cercano(const set<Point>& points, map<Point, set<Point>>& mapa){
    for(const Point& punto: points){
        double mas_cercano = numeric_limits<double>::max(); // el más grande de los doubles
        Point punto_mas_cercano = make_pair(0.0, 0.0);
        for(const auto& par : mapa){
            const Point clave = par.first;
            double distancia = distancia_cuadrado(clave, punto);
            if(distancia < mas_cercano){
                mas_cercano = distancia;
                punto_mas_cercano = clave;
            }
        }
        mapa[punto_mas_cercano].insert(punto);
    };
}

/** función que, dada una raíz de un árbol, encuentra los subárboles de altura h que estén en él 
 * nota: dado que h es la altura mínima de los árboles, entonces nunca ocurrirá que no se encuentre ninguno
*/
void busqueda_h(Nodo *nodo, const int h, map<Point, Nodo*>& arboles) {
    if(nodo==nullptr) return;
    vector<Entry> entradas = nodo->entries;
    if(entradas[0].a == nullptr){ // nodo es hoja
        return;
    }
    for(Entry entrada : entradas){
        Point clave = entrada.p;
        Nodo* hijo = entrada.a; // no es nullptr
        if(hijo->altura == h){ // encontramos un subárbol de altura h
            arboles[clave] = hijo;
        }
        else if(hijo->altura < h) continue;
        else{ // altura mayor a h
            busqueda_h(hijo, h, arboles);
        }
    }
}

/** Función que seta el radio cobertor de una entrada */
double setear_radio_cobertor(Entry& entry){ // no seteará las hojas porque originalmente ya son 0.0

    double max_radio = 0;
    Point punto = entry.p;
    Nodo *hijo = entry.a;
    if(hijo==nullptr){
        cout <<"no se como ponerle a esto" << endl;
        return 0.0;
    }
    for (Entry& entrada : hijo->entries) { // revisamos las entradas del hijo
        double radio_cobertor_hijo = entrada.cr;
        double distancia = sqrt(distancia_cuadrado(punto, entrada.p));
        if ((radio_cobertor_hijo + distancia)>max_radio){
            max_radio = radio_cobertor_hijo + distancia;
        }
    }

    entry.cr = max_radio;

    return entry.cr;
}

/** Función que conecta los subarboles Tj a las hojas de Tsup, arreglando las alturas en el proceso*/
void conectar_arboles(Nodo* nodo, map<Point,Nodo*>& subarboles){
    if(nodo==nullptr){
        cout << "coso captado" << endl; // pasó una vez
        return;
    }
    vector<Entry>& entries = nodo->entries;
    int altura_max = 0;
    int entries_sz = (int)entries.size();
    if(entries_sz == 0) return;
    if (entries[0].a==nullptr) { // encontramos una hoja
        nodo->altura = 2;
        for(int j = 0; j<entries_sz; j++){
            if(j>=entries.size()) break;
            Point punto_hoja = entries[j].p;
            if(subarboles[punto_hoja]==nullptr){
                //cout << "nompuedeser" << endl;
                continue;
            }
            entries[j].a = subarboles[punto_hoja];
            setear_radio_cobertor(entries[j]);
        }
    }
    else { // seguimos buscando hojas
        // Si el nodo no es una hoja, recorre recursivamente sus entradas
        for (auto& entry : entries) {
            conectar_arboles(entry.a, subarboles);
            int altura = ((entry.a)->altura)+1;
            if (altura > altura_max){
                altura_max = altura;
            }
            setear_radio_cobertor(entry);
        }
        nodo->altura = altura_max;
    }
    //cout << "la altura seteada fue " << nodo->altura << endl;
}


// Algoritmo CP:
// Input: Un set de puntos P

Nodo *crear_MTree_CCP(const set<Point> points){
    int n = points.size();
    //cout << "n es" << n << endl;
    if(n <= B){ // entran todos los puntos en un nodo
        Nodo *T = new Nodo; // se crea un árbol T (un nodo raíz con vector entries vacío)
        // se insertan todos los puntos a T
        for (Point punto: points){
            // crear entrada, inicialmente con radio cobertor y a nulos.
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
        int k = int(ceil(min(double(B), double(n)/B))); // debe ser el techo, tal que no llegue a 1!
        //cout << "k es" << k << endl;
        map<Point, set<Point>> samples;
        do{
            samples = map<Point, set<Point>>();
            //cout << "verificar que samples está vacío " << samples.size() << endl;
            vector<int> indices(n); // un vector de tamaño n
            for(int i = 0; i < n; ++i){
                indices[i] = i; // agregamos los números del 0 al n-1
            }
            // hacemos shuffle
            shuffle(indices.begin(), indices.end(), e);

            // los k primeros son los seleccionados!
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

            // armar los k conjuntos:
            punto_mas_cercano(points, samples);

            // etapa de redistribución

            auto it = samples.begin();
            while (it != samples.end()) {
                //Point clave = it->first;
                set<Point>& conjunto = it->second;
                //cout << "el conjunto mide " << conjunto.size() << endl;
                if (conjunto.size() < b_min) {
                    set<pair<double, double>> conjunto_copia = conjunto; // hacer una copia, espero esto funcione :c
                    it = samples.erase(it); // quitamos pfj de F
                    // redistribuimos
                    punto_mas_cercano(conjunto_copia, samples);
                    //cout << "deberia haberse redistribuido " << endl;
                } else {
                    ++it; // avanzamos al sgte elemento
                }
            }
        } while (samples.size() <= 1);

        int h = numeric_limits<int>::max(); // la altura mínima de los árboles Tj.

        map<Point, Nodo*>  subarboles;
        for(const auto& par : samples){
            Nodo *sub_arbol = crear_MTree_CCP(par.second);
            if(sub_arbol==nullptr){
                    cout << "algo malo" << endl;
            }
            vector<Entry> entradas = sub_arbol->entries;
            if(entradas.size() >= b_min){
                //cout << "todo ok con este arbol "<< endl;
                subarboles[par.first] = sub_arbol;

                if(sub_arbol->altura < h){
                    h = sub_arbol->altura;
                }
            }
            else{
                //cout << "rip, proceso sgte" << endl;
                for(auto entrada = entradas.begin(); entrada != entradas.end(); entrada++){
                    Nodo *sub_subarbol = entrada->a;
                    if(sub_subarbol==nullptr){
                        cout << "otro error pipipi" << endl;
                    }
                    subarboles[entrada->p] = sub_subarbol;
                    //cout << "la altura de este subarbol es" << sub_subarbol->altura << endl;
                    if(sub_subarbol->altura < h){
                        h = sub_subarbol->altura;
                    }
                }
                
            }
        }
        //cout << "la altura minima encontrada fue " << h << endl;

        //set<Nodo*> T2 = set<Nodo*>();

        // for(const auto& par_T_j : subarboles){
        //     if(par_T_j.second == nullptr){
        //         cout <<"ñaoñao" << endl;
        //     }
        //     //cout << "y los subarboles tienen un conj de tamaño" << (subarboles[par_T_j.first]->entries).size() << endl;
        // }
        
        set<Point> claves_a_eliminar;
        set<Point> keys; // las llaves de F (finalmente, F)
        map<Point, Nodo*> arboles_alt_h;

        for(const auto &par_T_j : subarboles){
            Point clave = par_T_j.first;
            Nodo *T_j = par_T_j.second;
            if(T_j == nullptr){
                cout << "nulo" << endl;
                //cout << clave.first << clave.second << endl;
                //continue;
            }
            if(T_j->altura== h){
                //T2.insert(T_j);
                keys.insert(clave);
                //cout << "T_j insertado" << endl;
            }
            else{
                // 9.1 Se borra el punto pertinente en F.
                claves_a_eliminar.insert(clave);
                //cout << "clave a eliminar" << endl;

                // 9.2 Se hace una búsqueda exhaustiva en Tj de todos los subárboles de altura igual a h.
                // Se insertan estos árboles a T′
                //cout << "se entró aquí" << endl;
                busqueda_h(T_j, h, arboles_alt_h); // funcion que busca y retorna un set con los subárboles de altura h que tiene T_j
                cout << "termina busqueda" << endl;
            }
        }

        for(Point clave : claves_a_eliminar){ // eliminamos las claves
            subarboles.erase(clave);
        }

        // 9.3 Se insertan los puntos raíz de T′1, . . . , T′ p: p′ f1, . . . , p′ fp en F
        for (const auto &par_T_j : arboles_alt_h) {
            Point clave = par_T_j.first;
            Nodo *T_j = par_T_j.second;
            keys.insert(clave);
            subarboles[clave] = T_j;
        }

        Nodo *Tsup = crear_MTree_CCP(keys);
        if(Tsup==nullptr){
            cout << "muy muy mal" << endl;
        }

        conectar_arboles(Tsup, subarboles);
        return Tsup;
    }
}

set<Point> crear_set(int n){
    set<Point> ccp_set;
    for(int j=0; j<n; j++){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dis(-100.0, 100.0);

        double first = dis(gen);
        double second = dis(gen);

        ccp_set.insert(make_pair(first, second));
    }
    return ccp_set;
}

/** Función que corrobora que todos los puntos del set original estén en las hojas del MTree */
void corroborar_MTree(Nodo *nodo, set<Point> puntos, set<Point>& puntos_encontrados){
    if(nodo==nullptr){
        cout << "nullptr" << endl;
        return;
    }
    vector<Entry> entradas = nodo->entries;
    for(const Point punto: puntos){ // Para cada punto, se busca en las hojas
        if(entradas[0].a == nullptr){ // nodo es una hoja
        // se verifica para cada entrada si p cumple con dist(p, q) ≤ r
            for(Entry entrada : entradas){
                Point otro_punto = entrada.p;
                if(punto.first==otro_punto.first && punto.second==otro_punto.second){
                    puntos_encontrados.insert(punto);
                }
            }
        }
        else { // no es una hoja... seguimos buscando por hojas!
            for(Entry entrada : entradas){
                corroborar_MTree(entrada.a, puntos, puntos_encontrados);
            }
        }
    }
}

int main() {
    set<pair<double, double>> random_pairs = crear_set(pow(2,19));
    cout << "el tamaño de random_pairs es " << random_pairs.size() << endl;
    //cout << "B es " << B << endl;
    Nodo *arbol = crear_MTree_CCP(random_pairs);
    for(Entry entry : arbol->entries){
        cout << "entry: " << entry.p.first << " " << entry.p.second << endl;
    }
    //set<Point> puntos_encontrados = set<Point>();
    //corroborar_MTree(arbol, random_pairs, puntos_encontrados);
    //cout << "se encontraron " << puntos_encontrados.size() << "puntos" << endl;
    deleteTree(arbol);
    return 0;
}