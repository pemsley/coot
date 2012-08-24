/*
     util/geomutil.cc: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/


#include <math.h>
#include <stdio.h>
#include <iostream>
#include "cartesian.h"
#include "geomutil.h"

Cartesian Mean(const std::vector<Cartesian> &v){
  if(v.size()==0)
    return Cartesian();

  Cartesian mean(0,0,0,0);
  for(unsigned i=0;i<v.size();i++)
    mean += v[i];
  mean /= v.size();
  return mean;
}

Cartesian PointAtWhichTangentToOneLineIntersectsAnotherLine(const Cartesian &p, const Cartesian &ls, const Cartesian &le,  const Cartesian &ols,  const Cartesian &ole){
  double l = p.DotProduct(ole-ols,le-ls)/p.DotProduct(p-ols,le-ls);
  Cartesian pp = ole * (1.0 /l) + ols * (1.0-1.0/l);
  return pp;
}

std::vector<double> DistanceBetweenPointAndLine(const Cartesian &ls, const Cartesian &le, const Cartesian &p){

  std::vector<double> ret(3);
  ret[0] = -1.0;
  ret[1] = -1.0;
  ret[2] = -1.0;

  double linesize = (le-ls).length();
  if(fabs(linesize)<1e-6){
    printf("Zero length line in DistanceBetweenPointAndLine\n");
    return ret;
  }
  /* t is value in line equation p = p1 + t(p2-p1) */
  double t = p.DotProduct(p-ls,le-ls) / (linesize*linesize);

  Cartesian pt = ls + t*(le-ls);

  ret[0] = (pt-p).length();
  ret[1] = t;
  return ret;

}

std::vector<double> DistanceBetweenTwoLines(const Cartesian &p1, const Cartesian &p2, const Cartesian &p3, const Cartesian &p4){

  Cartesian a1 = p2-p1;
  Cartesian a2 = p4-p3;
  
  std::vector<double> ret(3);
  ret[0] = -1.0;
  ret[1] = -1.0;
  ret[2] = -1.0;

  if(a1.length()==0)
    return ret;
  
  Cartesian dum;

  Cartesian n = dum.CrossProduct(a1,a2);
  if(fabs(n.length())<1e-6){
    return DistanceBetweenPointAndLine(p1,p2,p3);
  }

  if(a2.length()==0)
    return ret;
  
  n.normalize();

  double dist = fabs(dum.DotProduct(p4-p1,n));

  /*
\documentclass{article}
\usepackage[dvips]{color}
\usepackage{epsfig}
\usepackage{amsmath}

\begin{document}
Two lines beginning and ending with points ${\bf p_1}$, ${\bf p_2}$ and
${\bf p_1}$, ${\bf p_2}$, can be expressed:
\begin{align}
 {\bf r_1}  &= {\bf p_1} + t {\bf a_1}\\
 {\bf r_2}  &= {\bf p_3} + u {\bf a_2}
\end{align}
where ${\bf r_n}$ is a point along the line and $t$, $u$ are paramters such that
points with $<=0<t,u<=1$ fall between ${\bf p_1}$, ${\bf p_2}$ and
${\bf p_1}$, ${\bf p_2}$.
\begin{align}
 {\bf a_1}  &= {\bf p_2} - {\bf p_1}\\
 {\bf a_2}  &= {\bf p_4} - {\bf p_3} 
\end{align}
Assuming the two lines are not parallel and non intersecting, there exists 
an infinite number of planes that contain the line ${\bf a_1}$, though only one 
which does not intersect ${\bf a_2}$. ${\bf a_1}$ lies in this plane, and 
${\bf a_2}$ is parallel to the plane, so a normal to the plane $n$ is given by:
\begin{equation}
{\bf n} = {\bf a_1} \times {\bf a_2}
\end{equation}
The distance, d, between this plane and the plane parallel to this one which 
contains ${\bf p_3}$ and ${\bf p_4}$ 
is:
\begin{equation}
d = | ({\bf p_4} - {\bf p_1}) \cdot ({\bf n}/n) | 
\end{equation}

To find the points ${\bf Q_1}$, ${\bf Q_2}$ which are the points on
the two lines where they closest meet:
\begin{align}
 {\bf Q_1}  &= {\bf p_1} + t_Q {\bf a_1}\\
 {\bf Q_2}  &= {\bf p_3} + u_Q {\bf a_2}
\end{align}
we not that ${\bf Q_1} - {\bf Q_2}$ is perpendicular to ${\bf a_1}$
and ${\bf a_2}$
\begin{align}
{\bf a_1} \cdot {\bf Q_1} &= {\bf a_1} \cdot {\bf Q_2}\\
{\bf a_2} \cdot {\bf Q_1} &= {\bf a_2} \cdot {\bf Q_2}
\end{align}
So:
\begin{align}
{\bf a_1} \cdot {\bf p_1} + t_Q a_1^2 &= {\bf a_1} \cdot {\bf p_3} + u_Q {\bf a_1} \cdot {\bf a_2} \label{q1exp}\\
{\bf a_2} \cdot {\bf p_1} + t_Q {\bf a_2} \cdot {\bf a_1} &= {\bf a_2} \cdot {\bf p_3} + u_Q a_2^2 \label{q2exp}
\end{align}
Rearring \ref{q1exp} gives:
\begin{equation}
 t_Q = \frac{{\bf a_1} \cdot {\bf p_3} + u_Q {\bf a_1} \cdot {\bf a_2} - {\bf a_1} \cdot {\bf p_1}}{a_1^2}
\end{equation}
and substituting this into \ref{q2exp} and rerarranging we obtain:
\begin{equation}
 u_Q = \frac{a_1^2({\bf a_2} \cdot {\bf p_3} - {\bf a_2} \cdot {\bf p_1}) + ({\bf a_1} \cdot {\bf a_2})[{\bf a_1} \cdot {\bf p_1} - {\bf a_1} \cdot {\bf p_3}]}{({\bf a_1} \cdot {\bf a_2})^2 - a_1^2 a_2^2}
\end{equation}
Hence ${\bf Q_1}$ and ${\bf Q_1}$ can now be fully expressed in terms of ${\bf p_n}$.

Now we can tell if the intersection of the two lines lies in the intervals, ${\bf p_1}$,
${\bf p_2}$ and ${\bf p_3}$, ${\bf p_4}$:\newline
If $t<0$ or $t>1$ then intersection is outside of ${\bf p_1}$, ${\bf p_2}$.\newline
If $u<0$ or $u>1$ then intersection is outside of ${\bf p_3}$, ${\bf p_4}$.\newline
\end{document}  
*/
  
  double a1sq = a1.length()*a1.length();
  double a2sq = a2.length()*a2.length();

  double u = (a1sq*(dum.DotProduct(a2,p3) - dum.DotProduct(a2,p1)) + 
	     dum.DotProduct(a1,a2)*(dum.DotProduct(a1,p1) - dum.DotProduct(a1,p3)))/
	     (dum.DotProduct(a1,a2)*dum.DotProduct(a1,a2) - a1sq*a2sq);

  double t = (dum.DotProduct(a1,p3) + u * dum.DotProduct(a1,a2) - dum.DotProduct(a1,p1)) / a1sq;

  ret[0] = dist;
  ret[1] = t;
  ret[2] = u;

  /*
    printf("t:%f, u:%f\n",t,u);
    
    cout << "Q1: " << p1 + t*a1;
    cout << "Q2: " << p3 + u*a2;
    
    if(insegment==INLINE1&&(t<0||t>1)){
    return -1.0;
    }
    if(insegment==INLINE2&&(u<0||u>1)){
    return -1.0;
    }
    if(insegment==INLINE1AND2&&(u<0||u>1||t<0||t>1)){
    return -1.0;
  }
  */

  return ret;

}

Cartesian GetCartFrom3Carts(const Cartesian &Atom1, double blength, const Cartesian &Atom2, double angle1, const Cartesian &Atom3, double angle2, int chiral) 
{

  
  double n1dn2,xi,eta,zeta,costheta1,costheta2,temp,sintheta1,sintheta2,ang1,ang2;
  Cartesian r, n1, n2, n3, n1scaled;
  
  ang1=angle1;
  ang2=angle2;
  costheta1=cos(ang1);
  costheta2=cos(ang2);
  sintheta1=sin(ang1);
  sintheta2=sin(ang2);

  if (chiral==0){
    n1 = Atom2 - Atom1;
    n1.normalize();

    n2 = Atom3 - Atom2;
    //n2.normalize();

    double n1n2 = n1.DotProduct(n1,n2);
    n1scaled = n1n2*n1;

    n2 -= n1scaled;
    
    temp=1.0/sqrt(n2.DotProduct(n2,n2));
    
    if(fabs(temp)<1e-6)
      {
	printf("Int2Cart(): Collinearity detected. Coordinates not computed.\n");
	return r;
      };
    
    n2 *= temp;


    n3 = n1.CrossProduct(n1,n2);
    
    temp=1.0/sqrt(n3.DotProduct(n3,n3));
    n3 *= temp;
    
    xi=costheta1;
    eta=sintheta1*costheta2;
    zeta=-sintheta1*sintheta2;
  }else{
    n1 = Atom2 - Atom1;
    n1.normalize();
    
    n2 = Atom3 - Atom1;
    n2.normalize();

    n3 = n1.CrossProduct(n1,n2);
    n3.normalize();

    n1dn2=n1.DotProduct(n1,n2);
    
    if(fabs(n1dn2-1)<1e-6)
      {
	printf("Int2Cart(): Collinearity detected. Coordinates not computed.\n");
	return r;
      };
    
    eta=(n1dn2*costheta1-costheta2)/(n1dn2*n1dn2-1);
    xi=costheta1-eta*n1dn2;
    zeta=chiral*sqrt( fabs( 1-xi*xi-eta*eta-2*eta*xi*n1dn2 ) );
   }


  n1scaled = n1 * xi;
  Cartesian n2scaled = n2 * eta;
  Cartesian n3scaled = n3 * zeta;
  
  r  = n1scaled + n2scaled + n3scaled;
  r *= blength;
  
  r += Atom1;
  
  return r;
}

double LineLength( const Cartesian &A,  const Cartesian &B){

  Cartesian AB=A-B;
  return AB.length();
}

double DihedralAngle( const Cartesian &A, const Cartesian &B,  const Cartesian &C, const Cartesian &D){

  Cartesian AB = B - A;
  Cartesian BC = C - B;
  Cartesian CD = D - C;
  Cartesian Q = A.CrossProduct(AB,BC);
  Cartesian T = A.CrossProduct(BC,CD);
  Cartesian S = A.CrossProduct(Q,T);

  double ss = A.DotProduct(S,S);
  double qt = A.DotProduct(Q,T);

  double angle = atan2(sqrt(ss),qt);

  double sbc = S.get_x()*BC.get_x() + S.get_y()*BC.get_y() + S.get_z()*BC.get_z();

  if(sbc<0.0) {
    return -angle;
  } else {
    return angle;
  }
}

std::vector<double> LeastSquares2D(const std::vector<Cartesian> &p){

  matrix lhs(2,2);
  matrix rhs(2,1);

  double sumx = 0.0;
  double sumxsq = 0.0;
  double sumy = 0.0;
  double sumxy = 0.0;

  for(unsigned int i=0;i<p.size();i++){
    sumx += p[i].get_x();
    sumy += p[i].get_y();
    sumxsq += p[i].get_x()*p[i].get_x();
    sumxy += p[i].get_x()*p[i].get_y();
  }

  lhs(0,0) = sumxsq;
  lhs(0,1) = sumx;
  lhs(1,0) = sumx;
  lhs(1,1) = p.size();

  rhs(0,0) = sumxy;
  rhs(1,0) = sumy;

  std::vector<int> perm;
  int parity;

  matrix lu = matrix::LUDecomposition(lhs,perm,parity);

  matrix abmat = matrix::LUSubstitution(lu,rhs,perm);
  std::vector<double> ab(2);
  ab[0] = abmat(0,0);
  ab[1] = abmat(1,0);

  return ab;
  
}

std::vector<Cartesian> LeastSquaresOrtho3D(const std::vector<Cartesian> &p){

  double sumx=0.0, sumy=0.0, sumz=0.0;

  for(unsigned int i=0;i<p.size();i++){
    sumx += p[i].get_x();
    sumy += p[i].get_y();
    sumz += p[i].get_z();
  }

  Cartesian A(sumx/p.size(),sumy/p.size(),sumz/p.size());


  double a = A.get_x(), b = A.get_y(), c = A.get_z();

  double sum_xma_sq  = 0.0, sum_ymb_sq  = 0.0, sum_zmc_sq = 0.0;
  double sum_xma_ymb = 0.0, sum_xma_zmc = 0.0;
  double sum_ymb_zmc = 0.0;

  for(unsigned int i=0;i<p.size();i++){
     double xi = p[i].get_x(), yi = p[i].get_y(), zi = p[i].get_z();
     sum_xma_sq  += (xi-a)*(xi-a);
     sum_ymb_sq  += (yi-b)*(yi-b);
     sum_zmc_sq  += (zi-c)*(zi-c);
     sum_xma_ymb += (xi-a)*(yi-b);
     sum_xma_zmc += (xi-a)*(zi-c);
     sum_ymb_zmc += (yi-b)*(zi-c);
  }

  double delta = sum_xma_sq + sum_ymb_sq + sum_zmc_sq;

  matrix M(3,kdelta);
  M(0,0) *= delta;
  M(1,1) *= delta;
  M(2,2) *= delta;


  M(0,0) -= sum_xma_sq;
  M(1,1) -= sum_ymb_sq;
  M(2,2) -= sum_zmc_sq;

  M(0,1) -= sum_xma_ymb;
  M(0,2) -= sum_xma_zmc;
  M(1,2) -= sum_ymb_zmc;

  M(1,0) = M(0,1);
  M(2,0) = M(0,2);
  M(2,1) = M(1,2);


  std::vector<matrix> eigen = M.Eigen();
  eigen = M.SortEigenvalues(eigen);

  std::vector<Cartesian> AandD(2);
  AandD[0] = A;
  AandD[1] = Cartesian(eigen[0](0,0),eigen[0](1,0),eigen[0](2,0));

  return AandD;

}

Quat GetStandardRotation(const std::string &plane){
  Quat q;
  Cartesian x_axis(1,0,0);
  Cartesian y_axis(0,1,0);
  double theta = 90;

  if(plane==std::string("XY")){
    return q;
  } else if (plane==std::string("XZ")){
    Cartesian rot_axis(0,1,0);
    Cartesian rot_axis2(1,0,0);
    q = Quat(y_axis,1,theta);
    Quat q2(x_axis,1,theta);
    q.postMult(q2);
    return q;
  } else if (plane==std::string("YZ")){
    Cartesian rot_axis(0,1,0);
    Cartesian rot_axis2(1,0,0);
    q = Quat(x_axis,1,-theta);
    Quat q2(y_axis,1,-theta);
    q.postMult(q2);
    return q;
  }

  std::cout << "Unknown plane specifier: " << plane << "\n";
  return q;
}

/*
 * Algorithm for LeastSquaresPolyFit from
 * http://groups.google.com/group/comp.sys.hp48/browse_frm/thread/1d861a341cb1ec7/31436c9493504d08?lnk=st&q=&rnum=10&hl=en#31436c9493504d08
 *
Polynomial Regression==================== Theory 
======================================== 
Note: symbols between parentheses attached to a variable represent 
subindices, i.e, a(n) = a, sub-n.  Also, the symbol ^ is used to 
represent powers, i.e., x-squared is x^2. 
Polynomial regression consists of finding the coefficients 
b(0), b(1), ...,b(m), that will fit an m-th order polynomial of the form: 
y = b(0) + b(1)·x + b(2)·x^2 +... + b(m-1)·x^(m-1) + b(m)·x^m + e 
[1] 
Here, e represents a random error, which is not typically known. 
Typically you will use n data points given as a table: 
================ 
x(1)   x(2)     ...   x(n) 
y(1)   y(2)    ...   y(n) 
================ 
If we write the vectors y and b as: 
 _        _             _       _ 
|  y(1)    |           |  b(0)   |: 
|  y(2)    |           |  b(1)   | 
|   .        |   and   |   .       |, respectively, 
|   .        |           |   .        | 
|_ y(n) _|           |_ b(m) _| 
and a matrix A given by: 
 _                                        _ 
|  1  x(1)   x(1)^2 . . . x(1)^m | 
|  1  x(2)   x(2)^2 . . . x(2)^m | 
|  .     .      .            .                | 
|  .     .      .            .                | 
|_ 1  x(n)   x(n)^2 . . . x(n)^m_|, 
you can write the multiple-linear regression equation as: 
y = A·b + e,                           [2] 
where e is the random error vector, 
       _       _ 
      |  e(1)   | 
      |  e(2)   | 
e = |   .       | . 
      |   .       | 
      |_ e(n)_| 
Because the actual errors are not known, we can write only an 
approximation to the actual relationship as (the ^ here is used to 
distinguish the original data vector, y, from the approximation 
based on the polynomial in [1], y^. I.e., in [3] ^ does not mean a power): 
y^ = A·b.                          [3] 
The error vector is therefore: 
e = y - y^ = y-A·b                       [4] 
Or, in terms of the vector components: 
e(i) = y(i) - y^(i), i = 1,2,...,n.              [5] 
The sum of squared errors is: 
SSE= Sum(i=1 to n of e(i)^2) = eT·e = (y-A·b)T·(y-A·b)       [6] 
(eT = transpose of e) 
The least-square method requires us to find the (m+1) values of b(i), 
i = 0,2,..,m, that minimize the value of SSE in [6]. This minimization is 
obtained by finding the values of b so that: 
d(SSE)/db = 0.                          [7] 
The procedure to calculate the derivatives in [7] are complicated, 
however, [7] produces the following matrix equation: 
AT·A·b = AT·y,                          [8] (AT = transpose of A) 
which can be solved for b as: 
b = (AT·A)^(-1)·AT·y                           [9] 
 *
 */

std::vector<double> LeastSquaresPolyFit(const std::vector<double> &xs, std::vector<double> &ys, const int order){
  std::vector<double> b;
  if(order==1){
    double mean_x=0.0;
    double mean_y=0.0;
    for(unsigned i=0;i<xs.size();i++){
      mean_x += xs[i];
      mean_y += ys[i];
    }
    mean_x /= xs.size();
    mean_y /= xs.size();
    double xibarsq=0.0,xyibar=0.0;
    for(unsigned i=0;i<xs.size();i++){
      xibarsq += (xs[i] - mean_x)*(xs[i] - mean_x);
      xyibar += (xs[i] - mean_x)*(ys[i] - mean_y);
    }
    double m;
    if(fabs(xibarsq)>1e-8)
      m = xyibar/xibarsq;
    else m = 1e+8;
    double c = mean_y - m * mean_x;
    std::cout << "y = " << m << " * x + " << c << "\n";
    b.push_back(c);
    b.push_back(m);
  }

  if(order>1){
    matrix A(xs.size(),1+order);
    matrix y(xs.size(),1);
    for(unsigned i=0;i<xs.size();i++){
       A(i,0) = 1.0;
       A(i,1) = xs[i];
       A(i,2) = xs[i]*xs[i];
       for(int j=3;j<=order;j++){
         A(i,j) = pow(xs[i],j);
       }
       y(i,0) = ys[i];
    }
    matrix AT = A.Transpose();
    matrix bmat = (AT*A).Inverse()*AT*y;
    //std::cout << (AT*A).Inverse()*(AT*A) << "\n";
    //std::cout << bmat << "\n";
    for(unsigned i=0;i<bmat.get_rows();i++)
      b.push_back(bmat(i,0));
  }

  return b;
}

/* And now matrix A contains columns ABCDEF for polynomial of form:
 * z = Ax^2 + By^2 + Cxy + Dx + Ey + F
 * and amazingly it seems to work at least for test data.
 */

std::vector<double> LeastSquaresQuadraticFit3D(const std::vector<Cartesian> &carts){
  std::vector<double> b;
  matrix A(carts.size(),6);
  matrix z(carts.size(),1);
  for(unsigned i=0;i<carts.size();i++){
     A(i,0) = carts[i].get_x()*carts[i].get_x();
     A(i,1) = carts[i].get_y()*carts[i].get_y();
     A(i,2) = carts[i].get_x()*carts[i].get_y();
     A(i,3) = carts[i].get_x();
     A(i,4) = carts[i].get_y();
     A(i,5) = 1.0;
     z(i,0) = carts[i].get_z();
  }
  matrix AT = A.Transpose();
  matrix bmat = (AT*A).Inverse()*AT*z;
  //std::cout << (AT*A).Inverse()*(AT*A) << "\n";
  //std::cout << bmat << "\n";
  for(unsigned i=0;i<bmat.get_rows();i++)
    b.push_back(bmat(i,0));

  return b;

}

std::vector<double> LeastSquaresCubicFit3D(const std::vector<Cartesian> &carts){
  std::vector<double> b;
  matrix A(carts.size(),10);
  matrix z(carts.size(),1);
  for(unsigned i=0;i<carts.size();i++){
     A(i,0) = carts[i].get_x()*carts[i].get_x();
     A(i,1) = carts[i].get_y()*carts[i].get_y();
     A(i,2) = carts[i].get_x()*carts[i].get_y();
     A(i,3) = carts[i].get_x();
     A(i,4) = carts[i].get_y();
     A(i,5) = 1.0;
     A(i,6) = carts[i].get_x()*carts[i].get_x()*carts[i].get_x();
     A(i,7) = carts[i].get_y()*carts[i].get_y()*carts[i].get_y();
     A(i,8) = carts[i].get_x()*carts[i].get_x()*carts[i].get_y();
     A(i,9) = carts[i].get_x()*carts[i].get_y()*carts[i].get_y();
     z(i,0) = carts[i].get_z();
  }
  matrix AT = A.Transpose();
  matrix bmat = (AT*A).Inverse()*AT*z;
  //std::cout << (AT*A).Inverse()*(AT*A) << "\n";
  //std::cout << bmat << "\n";
  for(unsigned i=0;i<bmat.get_rows();i++)
    b.push_back(bmat(i,0));

  return b;

}
