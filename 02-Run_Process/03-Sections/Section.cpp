#include <cmath>
#include "Section.hpp"
#include "Definitions.hpp"

//Default constructor.
Section::Section(std::string name) : Name(name) {
    //Does nothing.
}

//Virtual destructor.
Section::~Section(){
    //Erase strings.
    Name.clear();
}

//Gets the Section's Name. 
std::string 
Section::GetName(){
    return Name;
}

//Gets Local Axis Rotation for Line Section according to provided angle.
Eigen::MatrixXd 
Section::GetLineRotationMatrix(double Theta){
    //Transformation Matrix for local axes.
    Eigen::MatrixXd T(3,3);
    T << 1.0,     0.0    ,    0.0    ,
         0.0,  cos(Theta), sin(Theta),
         0.0, -sin(Theta), cos(Theta);

    return T;
}

//Gets Local Axis Rotation for Area Section according to provided angle.
Eigen::MatrixXd 
Section::GetAreaRotationMatrix(double Theta){
    //Transformation Matrix for local axes.
    Eigen::MatrixXd T(3,3);
    T <<  cos(Theta), sin(Theta), 0.0,
         -sin(Theta), cos(Theta), 0.0, 
          0.0       , 0.0       , 1.0;

    return T;
}

//Gets centroid translation axis for Line Section according to insertion point.
Eigen::MatrixXd 
Section::GetLineTranslationMatrix(double h, double b, double zc, double yc, unsigned int ip){
    //Coordinate according to insertion point.
    double x2, x3;
    switch(ip) {
        case 1 : x3 = -zc         ; x2 =  yc        ; break;
        case 2 : x3 =  zc - b/2.0 ; x2 =  yc        ; break;
        case 3 : x3 =  b - zc     ; x2 =  yc        ; break;
        case 4 : x3 = -zc         ; x2 =  yc - h/2.0; break;
        case 5 : x3 =  zc - b/2.0 ; x2 =  yc - h/2.0; break;
        case 6 : x3 =  b - zc     ; x2 =  yc - h/2.0; break;
        case 7 : x3 = -zc         ; x2 =  h - yc    ; break;
        case 8 : x3 =  zc - b/2.0 ; x2 =  h - yc    ; break;
        case 9 : x3 =  b - zc     ; x2 =  h - yc    ; break;
        case 10: x3 =  0.0        ; x2 =  0.0       ; break;
        default : x3 = 0.0  ; x2 =  0.0; 
    }

    //Transformation Matrix for local axes.
    Eigen::MatrixXd A(3,3);
    A << 1.0,  x3,  x2,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

    return A;
}
