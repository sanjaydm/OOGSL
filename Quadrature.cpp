//Rafael Orozco 
//integrate a function  
//Uses quadrature method 
//Works for Dim 1- and order 2-6
#include "Quadrature.h"

Quadrature::Quadrature(int order, int dim) {
  _order = order;
  _dim   = dim;
  createNodesWeights();
}

Quadrature::~Quadrature(void) {
  
}

Vector& Quadrature :: weights() {
    static Vector weights(_order);
    weights = _weights;
    return weights;
}

vector<Vector>& Quadrature :: nodes()  {
    static vector<Vector> nodes (_nodes);
    return nodes;
} 


void Quadrature :: createNodesWeights() {
    //Grab node and weight values from gausslegendre.h and put them in structures
    _weights.setDim(_order);
    Vector curr_nodes(_order);
    vector<Vector> nodes;
  
    int i;
    int j;
    int num_points = _order / 2;

    if(_order == 2) {
        for(i = 0; i < num_points; i++) {
            _weights(i*2) = w2[i];
            _weights((i*2)+1) = w2[i];
 
            curr_nodes(i*2)     = -x2[i];  
            curr_nodes((i*2)+1) = x2[i]; 
        }
        for(j = 0; j < _dim; j++) {
            _nodes.push_back(curr_nodes);
        }    
    }

    if(_order == 3) {
       
        _weights(0)   =  w3[0];
        curr_nodes(0) =  x3[0];  
        for(i = 1; i < num_points+1; i++) {
            _weights((i*2)-1) = w3[i];
            _weights((i*2))   = w3[i];

            curr_nodes((i*2)-1) = -x3[i]; 
            curr_nodes((i*2))   =  x3[i];
       }     
       
       for(j = 0; j < _dim; j++) {
          _nodes.push_back(curr_nodes);
       }
 
    }

    if(_order == 4) {
        for(i = 0; i < num_points; i++) {
            _weights(i*2) = w4[i];
            _weights((i*2)+1) = w4[i];
    
            curr_nodes(i*2)     = -x4[i];  
            curr_nodes((i*2)+1) = x4[i]; 

       }
       for(j = 0; j < _dim; j++) {
           _nodes.push_back(curr_nodes);
       }
    
    }

    if(_order == 5) {
       
        _weights(0)   =  w5[0];
        curr_nodes(0) =  x5[0];  
        for(i = 1; i < num_points+1; i++) {
            _weights((i*2)-1) = w5[i];
            _weights((i*2)) = w5[i];

            curr_nodes((i*2)-1) = -x5[i]; 
            curr_nodes((i*2)) =  x5[i];
       }     
       
       for(j = 0; j < _dim; j++) {
          _nodes.push_back(curr_nodes);
       }
 
    }


    if(_order == 6) {
        for(i = 0; i < num_points; i++) {
            _weights(i*2) = w6[i];
            _weights((i*2)+1) = w6[i];

            curr_nodes(i*2)     = -x6[i];  
            curr_nodes((i*2)+1) = x6[i]; 
        }
        for(j = 0; j < _dim; j++) {
           _nodes.push_back(curr_nodes);
        }

    }
}

void Quadrature::setOrder(int order) {
  _order = order;

}

// int main(int argc, char** argv) {

//     int order = 6;
//     int dim = 1;
//     Quadrature quad_a(order ,dim);
//     quad_a.createNodesWeights();

//     printf("Test with order %d, and dim %d\n", order, dim);
//     vector<Vector>& myNodes = quad_a.nodes(); // S: should return nodes
   
//     int j;
//     for(j = 0; j < dim; j++) {
//         myNodes[j].print();
//     }

//     //Grab and print weights
//     printf("Grab and print weights  \n");
//     Vector & newVec = quad_a.weights(); 
//     newVec.print();
// }
