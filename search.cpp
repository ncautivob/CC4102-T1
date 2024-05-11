/**
Search Method
 * Valentina Alarcón
 * Naomi Cautivo
 * Máximo Flores
*/

#include <algorithm>  //shuffle
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

using namespace std;
using Point = pair<double, double>;

int B_solo = 4096;
struct Nodo;

typedef struct Entry {  // una entrada de un nodo
  pair<double, double> p;
  double cr;
  Nodo* a;
} Entry;

int B = B_solo / sizeof(Entry);  // (128)
int b_min = 0.5 * B;             // capacidad mínima (64)

typedef struct Nodo {
  vector<Entry> entries;
  int altura;

  // insertar nueva entrada en entries
  void insertarEntry(const Entry entrada) {
    if (entries.size() <
        B) {  // si todavía no se alcanza la capacidad máxima, todo bien
      entries.push_back(entrada);
    } else {
      cerr << "Nodo lleno" << endl;  // error
    }
  }
} Nodo;

double distancia_cuadrado(const Point& punto1, const Point& punto2) {
  // diferencia de coordenadas
  double dx = punto1.first - punto2.first;
  double dy = punto1.second - punto2.second;

  // Calcula el cuadrado de la distancia euclidiana
  double distanciaCuadrada = dx * dx + dy * dy;

  // Retorna la raíz cuadrada de la distancia cuadrada para obtener la distancia
  // euclidiana
  return distanciaCuadrada;
}

// La búsqueda tiene como input el M-Tree (T ) donde se buscará una query Q =
// (q, r), donde q es un punto y r es el radio de búsqueda. Es decir, (q, r)
// definen una bola. El objetivo es encontrar todos los puntos de T que residen
// dentro de esta (ver Figura 1). Para realizar la búsqueda, se verifica desde
// la raíz cada nodo hijo de esta: • Si el nodo es una hoja, se verifica para
// cada entrada si p cumple con dist(p, q) ≤ r. Si es así, se agrega p a la
// respuesta. • Si el nodo es interno (ver Figura 2), se verifica para cada
// entrada (p, cr, a) si dist(p, q) ≤ cr + r. Si es así, se busca en su hijo a
// posibles respuestas. Si no se cumple esa desigualdad, se descarta.

set<Point> search(Nodo* MTree, Point q, double r, set<Point>& resp) {
  // El objetivo es encontrar todos los puntos de T que residen dentro de la
  // bola (q,r)
  vector<Entry> entradas = MTree->entries;
  double r_square = pow(
      r,
      2);  // dado que la función de distancia entrega la distancia al cuadrado
  if (entradas[0].a == nullptr) {  // nodo es una hoja
    // se verifica para cada entrada si p cumple con dist(p, q) ≤ r
    for (Entry entrada : entradas) {
      Point punto = entrada.p;
      if (distancia_cuadrado(punto, q) <= r_square) {
        resp.insert(punto);
      }
    }
  } else {  // nodo es interno
    for (Entry entrada : entradas) {
      Point punto = entrada.p;
      double cr = entrada.cr;
      double cr_square = pow(cr, 2);
      if (distancia_cuadrado(punto, q) <= r_square + cr_square) {
        // se busca en su hijo a posibles respuestas
        search(entrada.a, q, r, resp);
      }
      // si no, se descarta
    }
  }
}