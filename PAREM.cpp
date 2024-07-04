#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <omp.h>
#include "Transition.hh"
#include "State.hh"
#include "AFN.hh"


// PAREM algorithm implementation
// - Input: Transition unordered_map, set of initial states, set of final states
// - Output: Set of reachable states