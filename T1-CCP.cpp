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

/** Por enunciado se mencionan algunas características de los nodos.
 * En primer lugar, cada nodo tendrá como capacidad B entradas en disco
 * Luego, un nodo acepta hasta B / sizeof(entry) entradas.
 * Para la etapa de experimentación se definirá B = 4096 Bytes.
 * 
*/

int B_bytes = 4096; // Esto dividido en sizeof(entry) equivale al máximo número de entradas de un nodo.

/** Definimos la estructura 'Entry' como una de las entradas de un nodo
 * p: un punto
 * cr: radio cobertor (covering radius) de este subárbol.
 * Esto es la máxima distancia que hay entre p y cualquier punto del subárbol relacionado a su entrada.
 * 
*/
typedef struct Entry {
    Point p;
    double cr;
    //double cr2 = 0.0;
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
void punto_mas_cercano(const set<Point>& points, map<Point, set<Point>>& mapa){
    //cout << "se entró a esta parte wii " << endl;
    //cout << "el mapa que se recibió tiene " << mapa.size() << "entradas " << endl;
    for(const Point& punto: points){
        // le asignamos su sample más cercano!
        double mas_cercano = numeric_limits<double>::max(); // el más grande de los doubles
        Point punto_mas_cercano = make_pair(0.0, 0.0);
        // aplicamos algoritmo
        for(const auto& par : mapa){
            const Point clave = par.first;
            double distancia = distancia_cuadrado(clave, punto);
            if(distancia < mas_cercano){
                mas_cercano = distancia;
                punto_mas_cercano = clave;
            }
        }
        // termina algoritmo
        mapa[punto_mas_cercano].insert(punto);
    };
    //for(const auto& par : mapa){
        //cout << "verifiquemos redistribucion" << par.second.size() << endl;
    //}
}

/** función que, dada una raíz de un árbol, encuentra los subárboles de altura h que estén en él 
 * nota: dado que h es la altura mínima de los árboles, entonces nunca ocurrirá que no se encuentre ninguno
*/
void busqueda_h(Nodo *nodo, const int h, set<Nodo*>& arboles) {
    if(nodo==nullptr) return;
    if(nodo->altura == h){
        arboles.insert(nodo);
    }
    else if (nodo->altura > h) { // si su altura es menor a h, sabemos que no tiene hijos de altura h
        for (const auto& hijo : nodo->entries) {
            // cada entry apunta a un nodo
            busqueda_h(hijo.a, h, arboles);
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
    // vamos a maximizar (radio_cobertor_hijo + distancia_al_hijo)
    for (Entry& entrada : hijo->entries) { // revisamos las entradas del hijo
        double radio_cobertor_hijo = entrada.cr;
        double distancia = sqrt(distancia_cuadrado(punto, entrada.p));
        if ((radio_cobertor_hijo + distancia)>max_radio){
            max_radio = radio_cobertor_hijo + distancia;
        }
    }
    
    entry.cr = max_radio;
    //cout << "radio seteado " << entry.cr << endl;

    return entry.cr;
}

double radio_2(Entry& entry, Point p){ // no seteará las hojas porque originalmente ya son 0.0

    double max_radio = 0;
    Nodo *hijo = entry.a;
    vector<Entry> entradas = hijo->entries;
    if(entradas[0].a == nullptr){ // hojas!
        for (Entry entrada : hijo->entries) { // revisamos las hojas del subarbol
            double distancia = sqrt(distancia_cuadrado(p, entrada.p));
            if (distancia>max_radio){
                max_radio = distancia;
            }
        }
    }
    else {
        // buscamos hasta encontrar las hojas de los subarboles
        for (Entry entrada : hijo->entries) { // revisamos las hojas del subarbol
            Nodo *subarbol = entrada.a;
            for (Entry entrada_sub : subarbol->entries){
                double max_local = radio_2(entrada_sub, p);
                if (max_local>max_radio){
                    max_radio = max_local;
                }
            }
        }

    }
    
    //entry.cr2 = max_radio;

    return max_radio;
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
            //cout << "el radio seteado fue" << entries[j].cr << endl;
            //radio_2(entries[j], entries[j].p);
            //cout << "el 2do radio seteado fue" << entries[j].cr2 << endl;
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
            //cout << "el radio seteado fue" << entry.cr << endl;
            //radio_2(entry, entry.p);
            //cout << "el 2do radio seteado fue" << entry.cr2 << endl;
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
                iteradores.insert(indices[i]); // recordar que los sets se ordenan de menor a mayor
            }

            auto punto_escogido = points.begin();
            int actual = 0;
            for(int iterador : iteradores){
                for(int i=actual; i<iterador; i++){ // avanza hasta el índice al que apunta it. comienza en actual pues va secuencialmente
                    actual++;
                    punto_escogido++;
                }
                // se llegó al punto del iterador! se agrega a samples con un un conjunto vacío
                //cout << (*punto_escogido).first << " " << (*punto_escogido).second << endl;
                samples[*punto_escogido] = set<Point>();
                //cout << "verificar size de samples" << samples.size() << endl;
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

        // 6. Se realiza recursivamente el algoritmo CP en cada Fj, obteniendo el árbol Tj

        //cout << " en efecto samples tiene tamaño " << samples.size() << endl;
        map<Point, Nodo*>  subarboles;
        for(const auto& par : samples){
            Nodo *sub_arbol = crear_MTree_CCP(par.second);
            if(sub_arbol==nullptr){
                    cout << "algo malo" << endl;
            }
            // cada llamada devuelve un nodo (raíz de árbol). luego, tendremos varios árboles Tj.
            vector<Entry> entradas = sub_arbol->entries;
            //cout << "este arbol tiene " << entradas.size() << " entradas "<< endl;
            if(entradas.size() >= b_min){
                //cout << "todo ok con este arbol "<< endl;
                subarboles[par.first] = sub_arbol;

                if(sub_arbol->altura < h){
                    h = sub_arbol->altura;
                }
            }
            else{
                //cout << "rip, proceso sgte" << endl;
                // se quita la raiz (en este caso, no se agrega a 'subarboles') y se trabaja con sus subárboles
                // para acceder a los subárboles, encontraremos los nodos a los que referencia cada entry
                // luego, estos nodos los agregaremos al set de subarboles
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
                //cout << "verifiquemos ingreso " << subarboles.size() << endl;
                //for(const auto& par_T_j : subarboles){
                    //cout << "donde los puntos del subarbol son" << par_T_j.second->entries.size() << endl;
                //}
                
            }
        }
        //cout << "la altura minima encontrada fue " << h << endl;

        // 8. Etapa de balanceamiento: Se define h como la altura mínima de los árboles Tj.
        // Se define T' (T2) inicialmente como un conjunto vacío.

        set<Nodo*> T2 = set<Nodo*>();

        // 9. Por cada Tj , si su altura es igual a h, se añade a T′.
        //cout << "subarboles tiene " << subarboles.size() << "subarboles" << endl;
        for(const auto& par_T_j : subarboles){
            if(par_T_j.second == nullptr){
                cout <<"ñaoñao" << endl;
            }
            //cout << "y los subarboles tienen un conj de tamaño" << (subarboles[par_T_j.first]->entries).size() << endl;
        }
        
        set<Point> claves_a_eliminar;
        set<Point> keys; // las llaves de F (finalmente, F)

        for(const auto &par_T_j : subarboles){
            Point clave = par_T_j.first;
            Nodo *T_j = par_T_j.second;
            if(T_j == nullptr){
                cout << "nulo" << endl;
                //cout << clave.first << clave.second << endl;
                //continue;
            }
            if(T_j->altura== h){
                T2.insert(T_j);
                keys.insert(clave);
                //cout << "T_j insertado" << endl;
            }
            else{
                // 9.2 Se hace una búsqueda exhaustiva en Tj de todos los subárboles de altura igual a h.
                // Se insertan estos árboles a T′
                set<Nodo*> empty_set;
                cout << "empieza busqueda" << endl;
                busqueda_h(T_j, h, empty_set); // funcion que busca y retorna un set con los subárboles de altura h que tiene T_j
                for(Nodo* nodo : empty_set){
                    if(nodo == nullptr){
                        cout << "era nullptr" << endl;
                    }
                    if(nodo->altura != h){
                        cout << "no tenia altura h"<< endl;
                    }
                }
                cout << "termina busqueda" << endl;
                cout << "empty_set logró encontrar " << empty_set.size() << "arboles de altura h" << endl;
                for (Nodo *subtree_h : empty_set) {
                    T2.insert(subtree_h);
                    if(subtree_h == nullptr){
                        cout << "era nulo aaa" << endl;
                    }
                    subarboles[clave] = subtree_h;
                    keys.insert(clave);
                    for (Entry entrada_h : subtree_h->entries){
                        Point clave_h = entrada_h.p;
                        // if(entrada_h.a == nullptr){
                        //     //cout << "ñaoooooo" << endl;
                        //     if(subtree_h->altura != 1){
                        //         cout << "algo va mal" << endl;
                        //     }
                        //     cout << "y la altura es" << h << endl;
                        // }

                        subarboles[clave_h] = entrada_h.a;
                        keys.insert(clave_h);

                    }
                    int tamaño = claves_a_eliminar.size();
                    if( tamaño > 0){
                        cout << "ola" << endl;
                        cout << "el tamaño de claves a eliminar es " << claves_a_eliminar.size() << endl;
                    }
                    // 9.3 Se insertan los puntos raíz de T′1, . . . , T′ p: p′ f1, . . . , p′ fp en F
                }

                // 9.1 Se borra el punto pertinente en F.
                claves_a_eliminar.insert(clave);
                cout << "clave a eliminar" << endl;
            }
        }


        
        for(Point clave : claves_a_eliminar){
            //cout << "eliminando clave" << endl;
            subarboles.erase(clave);
        }
        // if (tamaño > 0){
        //     cout << "se terminó de eliminar claves" << endl;
        // }
        

        // 10. Se define Tsup como el resultado de la llamada al algoritmo CP aplicado a F.
        //cout << "el numero de keys es " << keys.size() << endl;
        Nodo *Tsup = crear_MTree_CCP(keys);
        if(Tsup==nullptr){
            cout << "muy muy mal" << endl;
        }

        // 11. Se une cada Tj ∈ T′ a su hoja en Tsup correspondiente al punto pfj ∈ F, obteniendo un nuevo árbol T.

        // Para lograr esto, debo encontrar las hojas de Tsup.
        // Hipótesis útil: si la primera entrada apunta a un puntero nulo, entonces llegamos a una hoja.
        // Tan pronto encuentro un punto, debo asociarlo a su conjunto.
        // Si el nodo es una hoja (no tiene entradas), imprime sus puntos

        conectar_arboles(Tsup, subarboles);
        // 12. Se setean los radios cobertores resultantes para cada entrada en este árbol.
        // 13. Se retorna T .
        return Tsup;
    }
}

set<Point> crear_set(int n){
    set<Point> ccp_set;
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

int main() {
    set<pair<double, double>> random_pairs = crear_set(pow(2,15));
    //cout << "B es " << B << endl;
    Nodo *arbol = crear_MTree_CCP(random_pairs);
    for(Entry entry : arbol->entries){
        cout << "entry: " << entry.p.first << " " << entry.p.second << endl;
    }
    deleteTree(arbol);
    return 0;
}